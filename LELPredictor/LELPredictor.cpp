/* =========================================================================================================== *
 * EPA ToxCast Challenge / TopCoder 
 * PredictLEL v1 [Competition Sensitive]
 * RigelFive - Hudson Ohio USA
 * 12 May 2014 - 14 May 2014
 * =========================================================================================================== *
 * version history:
 * v0: Utilize Mw, atoms, bonds, chains, groups and rings as inputs.
 * v1: Utilize Mw as the only input.
 * 
 * =========================================================================================================== *
 * Atomized Software Specification:
 * Customer Requirements:
 *  1. Predict the LEL value for untested compounds using previous tests of other chemicals.
 * 
 * Functional Requirements:
 *  1. Estimate the LEL value based on information provided from untested compounds.
 * 
 * Interface Requirements:
 *  1. Output a double[] with the negative log value of the LEL estimated.
 *   
 * Performance Requirements:
 *  1. No time or memory requirements as no computations are required in the submitted code
 * 
 * Design Requirements:
 *  1. Utilize C++.
 *  2. Model a NN to calculate the LEL value using a pseudo Joback-like method where 
 *     the main input parameters are: Mw, atoms, bonds, chains, groups and rings.
 *  3. Determine generalized extrinsic parameters can be estimated from specific intrinsic parameters.
 *  4. Do not use external databases related to the EPA or other agencies to determine the LEL values.
 * 
 * Quality Requirements:
 *  1. Maximize overall score of the algorithm (objective: 1,000,000 points).
 * 
 * Verification Requirements:
 *  1. Output 1854 values in a double[] array.
 * 
 * =========================================================================================================== *
 */

#define _LELPredictor_ 0; // output the version number

#include <cstdlib>
#include <cassert>              // remove
#include <cmath>                // remove
#include <fstream>              // remove
#include <sstream>              // remove
#include <algorithm>            // remove
#include <cstdio>
#include <ctime>                // remove
#include <vector>
#include <iostream>
#include <boost/tokenizer.hpp>  // remove

using namespace std;

struct Chemical {
    string chemical_gs_id;
    string chemical_casrn;
    string chemical_name;
    bool lel_value_known;
    double systemic_adjusted_negative_log_lel;
    string chemical_short_gs_id;   // truncate the "DSSTox_GSID_" and keep the id number in a string 
    string chemical_cid;     // matches the id in the DSSTox_Chemical and ToxCast_Chemical data, found in the DSSTox_Chemical_List
    string chemical_SMILES;
    float chemical_Mw;
};
struct DSSTox_Chemical {
     string DSSTox_GSID;
     string DSSTox_CID;
     string TS_ChemName;
     string TS_ChemName_Synonyms;
     string TS_CASRN;
     string ChemNote;
     string STRUCTURE_Shown;
     string STRUCTURE_Formula;
     string STRUCTURE_MW;
     string STRUCTURE_ChemType;
     string STRUCTURE_IUPAC;
     string STRUCTURE_SMILES;
     string STRUCTURE_SMILES_Desalt;
};
struct ToxCast_Chemical {
    string dsstox_cid; //DSSTox_CID 
};

struct Connection {
    double weight;
    double deltaWeight;
};

class Neuron;
typedef vector<Neuron> Layer;
class Neuron {
public:
    Neuron(unsigned numOutputs, unsigned myIndex);
    void setOutputVal(double val) { m_outputVal = val; }
    double getOutputVal(void) const { return m_outputVal; }
    void feedForward(const Layer &prevLayer);
    void calcOutputGradients(double targetVal);
    void calcHiddenGradients(const Layer &nextLayer);
    void updateInputWeights(Layer &prevLayer); 

private:
    static double eta;
    static double alpha;
    static double transferFunction(double x);
    static double transferFunctionDerivative(double x);
    static double randomWeight(void) { return rand() / double(RAND_MAX); }
    double sumDOW(const Layer &nextLayer) const;
    double m_outputVal;
    vector<Connection> m_outputWeights;
    unsigned m_myIndex;
    double m_gradient;
};
double Neuron::eta = 0.15;  
double Neuron::alpha = 0.50;
void Neuron::updateInputWeights(Layer &prevLayer) {
    for (unsigned n = 0; n < prevLayer.size(); ++n) {
        Neuron &neuron = prevLayer[n];
        double oldDeltaWeight = neuron.m_outputWeights[m_myIndex].deltaWeight;
        double newDeltaWeight =
                eta
                * neuron.getOutputVal()
                * m_gradient
                + alpha
                * oldDeltaWeight;
        neuron.m_outputWeights[m_myIndex].deltaWeight = newDeltaWeight;
        neuron.m_outputWeights[m_myIndex].weight += newDeltaWeight;
    }
}
double Neuron::sumDOW(const Layer &nextLayer) const {
    double sum = 0.0;

    for (unsigned n = 0; n < nextLayer.size() - 1; ++n) {
        sum += m_outputWeights[n].weight * nextLayer[n].m_gradient;
    }
    return sum;
}
void Neuron::calcHiddenGradients(const Layer &nextLayer) {
    double dow = sumDOW(nextLayer);
    m_gradient = dow * Neuron::transferFunctionDerivative(m_outputVal);
}
void Neuron::calcOutputGradients(double targetVal) {
    double delta = targetVal - m_outputVal;
    m_gradient = delta * Neuron::transferFunctionDerivative(m_outputVal);
}
double Neuron::transferFunction(double x) {
    return tanh(x);
}
double Neuron::transferFunctionDerivative(double x) {
    return 1.0 - x * x;
}
void Neuron::feedForward(const Layer &prevLayer) {
    double sum = 0.0;

    for (unsigned n = 0; n < prevLayer.size(); ++n) {
        sum += prevLayer[n].getOutputVal() *
                prevLayer[n].m_outputWeights[m_myIndex].weight;
    }
    m_outputVal = Neuron::transferFunction(sum);
}
Neuron::Neuron(unsigned numOutputs, unsigned myIndex) {
    for (unsigned c = 0; c < numOutputs; ++c) {
        m_outputWeights.push_back(Connection());
        m_outputWeights.back().weight = randomWeight();
    }
    m_myIndex = myIndex;
}

// ****************** class Net ******************
class Net
{
public:
    Net(const vector<unsigned> &topology);
    void feedForward(const vector<double> &inputVals);
    void backProp(const vector<double> &targetVals);
    void getResults(vector<double> &resultVals) const;
    double getRecentAverageError(void) const { return m_recentAverageError; }

private:
    vector<Layer> m_layers; // m_layers[layerNum][neuronNum]
    double m_error;
    double m_recentAverageError;
    static double m_recentAverageSmoothingFactor;
};

double Net::m_recentAverageSmoothingFactor = 100.0; // Number of training samples to average over

void Net::getResults(vector<double> &resultVals) const
{
    resultVals.clear();

    for (unsigned n = 0; n < m_layers.back().size() - 1; ++n) {
        resultVals.push_back(m_layers.back()[n].getOutputVal());
    }
}

void Net::backProp(const vector<double> &targetVals)
{
    // Calculate overall net error (RMS of output neuron errors)

    Layer &outputLayer = m_layers.back();
    m_error = 0.0;

    for (unsigned n = 0; n < outputLayer.size() - 1; ++n) {
        double delta = targetVals[n] - outputLayer[n].getOutputVal();
        m_error += delta * delta;
    }
    m_error /= outputLayer.size() - 1; // get average error squared
    m_error = sqrt(m_error); // RMS

    // Implement a recent average measurement

    m_recentAverageError =
            (m_recentAverageError * m_recentAverageSmoothingFactor + m_error)
            / (m_recentAverageSmoothingFactor + 1.0);

    // Calculate output layer gradients

    for (unsigned n = 0; n < outputLayer.size() - 1; ++n) {
        outputLayer[n].calcOutputGradients(targetVals[n]);
    }

    // Calculate hidden layer gradients

    for (unsigned layerNum = m_layers.size() - 2; layerNum > 0; --layerNum) {
        Layer &hiddenLayer = m_layers[layerNum];
        Layer &nextLayer = m_layers[layerNum + 1];

        for (unsigned n = 0; n < hiddenLayer.size(); ++n) {
            hiddenLayer[n].calcHiddenGradients(nextLayer);
        }
    }

    // For all layers from outputs to first hidden layer,
    // update connection weights

    for (unsigned layerNum = m_layers.size() - 1; layerNum > 0; --layerNum) {
        Layer &layer = m_layers[layerNum];
        Layer &prevLayer = m_layers[layerNum - 1];

        for (unsigned n = 0; n < layer.size() - 1; ++n) {
            layer[n].updateInputWeights(prevLayer);
        }
    }
}

void Net::feedForward(const vector<double> &inputVals)
{
    assert(inputVals.size() == m_layers[0].size() - 1);

    // Assign (latch) the input values into the input neurons
    for (unsigned i = 0; i < inputVals.size(); ++i) {
        m_layers[0][i].setOutputVal(inputVals[i]);
    }

    // forward propagate
    for (unsigned layerNum = 1; layerNum < m_layers.size(); ++layerNum) {
        Layer &prevLayer = m_layers[layerNum - 1];
        for (unsigned n = 0; n < m_layers[layerNum].size() - 1; ++n) {
            m_layers[layerNum][n].feedForward(prevLayer);
        }
    }
}

Net::Net(const vector<unsigned> &topology)
{
    unsigned numLayers = topology.size();
    for (unsigned layerNum = 0; layerNum < numLayers; ++layerNum) {
        m_layers.push_back(Layer());
        unsigned numOutputs = layerNum == topology.size() - 1 ? 0 : topology[layerNum + 1];

        // We have a new layer, now fill it with neurons, and
        // add a bias neuron in each layer.
        for (unsigned neuronNum = 0; neuronNum <= topology[layerNum]; ++neuronNum) {
            m_layers.back().push_back(Neuron(numOutputs, neuronNum));
        }

        // Force the bias node's output to 1.0 (it was the last neuron pushed in this layer):
        m_layers.back().back().setOutputVal(1.0);
    }
}

class ExampleChemicalData
{
public:
    ExampleChemicalData(const string filename);    
    void getChemicalList(vector<string> &ChemicalList);
    bool isEof(void) { return m_exampleChemicalDataFile.eof(); }    
private:
    ifstream m_exampleChemicalDataFile;  
};
ExampleChemicalData::ExampleChemicalData(const string filename)
{
    m_exampleChemicalDataFile.open(filename.c_str());
} 
void ExampleChemicalData::getChemicalList(vector<string> &ChemicalList) {
    string line;

    getline(m_exampleChemicalDataFile, line);
    ChemicalList.push_back(line);    
    return;    
}

class DSSToxData
{
public:
    DSSToxData(const string filename);    
    void getDSSToxList(vector<string> &DSSToxList);
    bool isEof(void) { return m_exampleDSSToxDataFile.eof(); }    
private:
    ifstream m_exampleDSSToxDataFile;  
};
DSSToxData::DSSToxData(const string filename)
{
    m_exampleDSSToxDataFile.open(filename.c_str());
} 
void DSSToxData::getDSSToxList(vector<string> &DSSToxList) {
    string line;

    getline(m_exampleDSSToxDataFile, line);
    DSSToxList.push_back(line);    
    return;
}

class ToxCastData
{
public:
    ToxCastData(const string filename);    
    void getToxCastList(vector<string> &ToxCastList);
    bool isEof(void) { return m_exampleToxCastDataFile.eof(); }    
private:
    ifstream m_exampleToxCastDataFile;  
};
ToxCastData::ToxCastData(const string filename)
{
    m_exampleToxCastDataFile.open(filename.c_str());
} 
void ToxCastData::getToxCastList(vector<string> &ToxCastList) {
    string line;

    getline(m_exampleToxCastDataFile, line);
    ToxCastList.push_back(line);    
    return;    
}

string chop_dsstox (string id_to_chop) {
    // truncate the "DSSTox_GSID_" and keep the id number in a string  
    
    stringstream ss(id_to_chop);    
    string blank1;
    string blank2;
    string id;

    getline(ss,blank1,'_') &&       
    getline(ss,blank2,'_') &&
    getline(ss,id); 
    
    return id;    
}

vector<string> tokenize(const string& line )
{
    boost::escaped_list_separator<char> sep( '\\', ',', '"' ) ;
    boost::tokenizer< boost::escaped_list_separator<char> > tokenizer( line, sep ) ;
    return std::vector<std::string>( tokenizer.begin(), tokenizer.end() ) ;
}

void output_c(Chemical c) {
    cout << c.chemical_gs_id << " " << 
            c.chemical_cid << " " <<            
            c.chemical_casrn << " " << 
            c.chemical_name << " " << 
            c.chemical_Mw << " " <<
            c.chemical_SMILES  << " " <<
            c.systemic_adjusted_negative_log_lel << endl;   
}

void output_dc(DSSTox_Chemical dc) {
    cout << dc.DSSTox_GSID << " " <<
            dc.DSSTox_CID << " " <<                           
            dc.TS_ChemName << " " <<
            dc.TS_ChemName_Synonyms << " " <<
            dc.TS_CASRN << " " <<
            dc.ChemNote << " " <<
            dc.STRUCTURE_Shown << " " <<
            dc.STRUCTURE_Formula << " " <<
            dc.STRUCTURE_MW << " " <<
            dc.STRUCTURE_ChemType << " " <<
            dc.STRUCTURE_IUPAC << " " <<                         
            dc.STRUCTURE_SMILES << " " <<
            dc.STRUCTURE_SMILES_Desalt << endl;    
}

void output_tc(ToxCast_Chemical tc) {
    cout << tc.dsstox_cid << " " <<
            endl;
}

void showVectorVals(string label, vector<double> &v)
{
    cout << label << " ";
    for (unsigned i = 0; i < v.size(); ++i) {
        cout << v[i] << " ";
    }

    cout << endl;
}

double sum(vector<double> v) {
    double sum = 0.0;
    
    int vSIZE = v.size();     
    for (int i = 0; i < vSIZE; i++) {
        sum = sum + double(v.at(i));
    }
    return double(sum);
}
double average(vector<double> v) {  
    double sum = 0.0;
    
    int vSIZE = v.size();     
    for (int i = 0; i < vSIZE; i++) {
        sum = sum + abs(v.at(i));
    }
    return double(sum/vSIZE);
}
double stddev2(vector<double> v) {
    double s2;
    double sum = 0.0;
    double sum2 = 0.0;  
    double avg;
    
    int vSIZE = v.size();
    
    for (int i = 0; i < vSIZE; i++) {
        sum = sum + abs(v.at(i));
    }
    avg = double(sum / vSIZE);
    
    for (int i = 0; i < vSIZE; i++) {
        sum2 = sum2 + pow(((double) v.at(i) - avg),2.0);
    }
    
    s2 = (1.0/((double) vSIZE-1.0)) * abs(sum2);
    
    return s2;    
}
double min(vector<double> v) {
    double min1 = 1e+09;
    int vSIZE = v.size();
    
    for (int i = 0; i < vSIZE; i++) {
        if (v.at(i) < min1) { min1 = v.at(i); };
    }

    return double(min1);    
}
double max(vector<double> v) {
    double max1 = 0;

    int vSIZE = v.size();
    
    for (int i = 0; i < vSIZE; i++) {
        if (v.at(i) > max1) { max1 = v.at(i); };
    }

    return double(max1);
}

double parabola(double a, double b, double c, double x) {
    double y;
    
    y = (a*pow((x-b),2.0)) + c;
    
    return y;
}
double parabolic_filter(double x, double x_mean, double x_min, double x_max) {
    double y, norm_x;
    double a1, a2, b, c1, c2;
    
    norm_x = x / x_max;
    
    a1 = 1.0 / pow((x_mean-x_min),2.0);
    a2 = 1.0 / pow((x_max-x_mean),2.0);
    b = x_mean;
    c1 = 0.0;
    c2 = 0.0;

    if (x < 0.0) {
        y = -1.0;
    } else if (x > x_max) {
        y = 1.0;
    } else {
        if (x <= x_mean) {
            y = -parabola(a1, b, c1, x);
        }
        else if (x > x_mean) {
            y = parabola(a2, b, c2, x);
        }
    }
    return y;
}
double parabolic_amplifier(double y, double x_mean, double x_min, double x_max) {
    
    double x;
    double a1, a2, b, c1, c2;
    
    a1 = 1.0 / pow((x_mean-x_min),2.0);
    a2 = 1.0 / pow((x_max-x_mean),2.0);
    b = x_mean;
    c1 = 0.0;
    c2 = 0.0;    

    if (y > 0) {
        x = (1.0/(2.0*a2))*((2.0*a2*b)+(2.0*sqrt(a2*y-a2*c2)));
    } else {
        x = (1.0/(2.0*a1))*((2.0*a1*b)-(2.0*sqrt(-a1*y-a1*c1)));
    }
    return x;
}

int main(int argc, char** argv) {
    
    string path = "/home/scott/Arctria/ARC-32/LELPredictor/";
    ExampleChemicalData exChemicalData("/home/scott/Arctria/ARC-32/LELPredictor/data/ToxCast_MM_Data/ToxRefDB_Challenge_Training.csv");
    DSSToxData exDSSToxData("/home/scott/Arctria/ARC-32/LELPredictor/data/DSSTox/csv/TOX21S_v4a_8599_11Dec2013.csv");
    ToxCastData exToxCastData("/home/scott/Arctria/ARC-32/LELPredictor/data/DSSTox/csv/toxprint_v2_vs_TOX21S_v4a_8599_03Dec2013.csv"); 
    
    // Initialize parameters
    vector<double> the_big_answer;
    
    vector<string> tokenized_string_vector;    
    
    vector<string> ChemicalList;
    vector<string> DSS_ToxData_ChemicalList;
    vector<string> ToxCast_ChemicalList;
    vector<double> LEL_values;
    vector<double> Mw_values;
    vector<double> Predicted_LEL_values;
    
    vector<Chemical> Initial_Chemicals;
    vector<Chemical> Training_Chemicals;
    vector<Chemical> Predict_Chemicals;
    vector<DSSTox_Chemical> DSS_Chemicals;
    vector<ToxCast_Chemical> ToxCast_Chemicals;

    Chemical c;
    DSSTox_Chemical dc;
    ToxCast_Chemical tc;
    Chemical pc;
    
    string in; 
    string blank_line;
    
    int val;
    double d_val;
    
    double Mw_min;
    double Mw_mean;
    double Mw_max;
    
    double LEL_min;
    double LEL_mean;
    double LEL_max;

    // LELPredict Program Timing Parameters for Training NN
    bool startTraining = false;
    bool finishTraining = false;
    bool startCalculations = false;
    bool finishCalculations = false;   
        
    // Program Timing Parameters
    double training_start_time = 0.0;
    double training_finish_time = 0.0;
    double training_elapsed_time = 0.0;
    double calculation_start_time = 0.0;
    double calculation_finish_time = 0.0;
    double calculation_elapsed_time = 0.0;   
    
    // Text analysis
    bool quote_string = false;
    
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if ((arg == "-t") || (arg == "--training")) {
            cout << "LELPredictor v00\nRigelFive - 12 May 2014" << endl;  
            cout << "================================================================================" << endl;

            // Start the training
            if (startTraining == false) {
                training_start_time = clock();               
                startTraining = true;
            } 

            // Read the .csv training file to build the Chemical List
            while (!exChemicalData.isEof()) { 
                exChemicalData.getChemicalList(ChemicalList);
            } 

            // Parse the data using comma separated values from the Chemical List
            for (int i = 1; i < ChemicalList.size()-1; i++) {
                in = ChemicalList.at(i);
                
                tokenized_string_vector.clear();                
                boost::tokenizer<boost::escaped_list_separator<char> > tok(in);
                
                for(boost::tokenizer<boost::escaped_list_separator<char> >::iterator beg=tok.begin(); beg!=tok.end();++beg){
                    tokenized_string_vector.push_back(*beg);
                }                
                
                c.chemical_gs_id = tokenized_string_vector.at(0);
                c.chemical_casrn = tokenized_string_vector.at(1);
                c.chemical_name = tokenized_string_vector.at(2);
                c.chemical_short_gs_id = chop_dsstox(c.chemical_gs_id);
                c.chemical_cid = "-";
                c.chemical_SMILES = "-"; // temporary... changed to real value looked up from DSS_Chemicals
                c.chemical_Mw = 0.0; // temporary... changed to real value looked up from DSS_Chemicals             

                blank_line = tokenized_string_vector.at(3);
                
                if (blank_line.find('?') != string::npos) {
                    c.lel_value_known = false;
                    c.systemic_adjusted_negative_log_lel = -1.0;
                    Predict_Chemicals.push_back(c); 
                    Initial_Chemicals.push_back(c);
                    
                } else {
                    c.lel_value_known = true;               
                    istringstream(blank_line) >> c.systemic_adjusted_negative_log_lel; 
                    LEL_values.push_back(c.systemic_adjusted_negative_log_lel);
                    Training_Chemicals.push_back(c); 
                    Initial_Chemicals.push_back(c);                    
                }
            }

            LEL_min = min(LEL_values);
            LEL_mean = average(LEL_values);
            LEL_max = max(LEL_values);
            
            cout << "Training Set Statistics:" << endl;
            cout << "   Min LEL Value: " << min(LEL_values) << endl;
            cout << "   Max LEL Value: " << LEL_max << endl;
            cout << "   Average LEL Value: " << LEL_mean << endl;
            cout << "   Sigma^2 of LEL Value: " << stddev2(LEL_values) << endl;
            cout << "--------------------------------------------------------------------------------" << endl;              

            // Build DSSTox_Chemical List
            while (!exDSSToxData.isEof()) { 
                exDSSToxData.getDSSToxList(DSS_ToxData_ChemicalList);
            }

            for (int i = 1; i < DSS_ToxData_ChemicalList.size()-1; i++) {
                in = DSS_ToxData_ChemicalList.at(i);
                tokenized_string_vector.clear();                
                
                boost::escaped_list_separator<char> sep( '`', ',', '\"' ) ;         // apparently the character ` is never used!!!       
                boost::tokenizer<boost::escaped_list_separator<char> > tok(in, sep);
                
                for(boost::tokenizer<boost::escaped_list_separator<char> >::iterator beg=tok.begin(); beg!=tok.end();++beg) {
                    tokenized_string_vector.push_back(*beg);
                }  

                dc.DSSTox_GSID = tokenized_string_vector.at(0);
                dc.DSSTox_CID = tokenized_string_vector.at(1);                           
                dc.TS_ChemName = tokenized_string_vector.at(2);
                dc.TS_ChemName_Synonyms = tokenized_string_vector.at(3);
                dc.TS_CASRN = tokenized_string_vector.at(4); 
                dc.ChemNote = tokenized_string_vector.at(5);
                dc.STRUCTURE_Shown = tokenized_string_vector.at(6);
                dc.STRUCTURE_Formula = tokenized_string_vector.at(7);
                dc.STRUCTURE_MW = tokenized_string_vector.at(8);
                dc.STRUCTURE_ChemType = tokenized_string_vector.at(9);
                dc.STRUCTURE_IUPAC = tokenized_string_vector.at(10);                           
                dc.STRUCTURE_SMILES = tokenized_string_vector.at(11);
                dc.STRUCTURE_SMILES_Desalt = tokenized_string_vector.at(12);     
                
                DSS_Chemicals.push_back(dc); 
                
                istringstream(tokenized_string_vector.at(8)) >> d_val;
                Mw_values.push_back(d_val);
            }     
            
            Mw_min = min(Mw_values);
            Mw_mean = average(Mw_values);
            Mw_max = max(Mw_values);
            
            cout << "   Min Mw Value: " << min(Mw_values) << endl;
            cout << "   Max Mw Value: " << Mw_max << endl;
            cout << "   Average Mw Value: " << Mw_mean << endl;
            cout << "   Sigma^2 of Mw Value: " << stddev2(Mw_values) << endl;
            cout << "--------------------------------------------------------------------------------" << endl;                
                       
            // Build ToxCast_Chemical List
            while (!exToxCastData.isEof()) { 
                exToxCastData.getToxCastList(ToxCast_ChemicalList);
            }  
            
            for (int i = 1; i < ToxCast_ChemicalList.size()-1; i++) {
                in = ToxCast_ChemicalList.at(i);
                tokenized_string_vector.clear();                
                
                boost::escaped_list_separator<char> sep( '\\', ',', '\"' ) ;      
                boost::tokenizer<boost::escaped_list_separator<char> > tok(in, sep);
                
                for(boost::tokenizer<boost::escaped_list_separator<char> >::iterator beg=tok.begin(); beg!=tok.end();++beg) {
                    tokenized_string_vector.push_back(*beg);
                }  
                
                tc.dsstox_cid = tokenized_string_vector.at(0);
                
                ToxCast_Chemicals.push_back(tc);

            }            

            // Output basic information about data in memory       
            cout << "   Total Chemicals Listed: " << ChemicalList.size()-2 << endl;     // reduced by 2 for the header line and footer line       
            cout << "   Size of Training Chemicals List: " << Training_Chemicals.size() << endl;
            cout << "   Size of Prediction Chemicals List: " << Predict_Chemicals.size() << endl;
            cout << "--------------------------------------------------------------------------------" << endl;                  
            cout << "   Total DSS Tox Chemicals: " << DSS_Chemicals.size() << endl;
            cout << "   Total ToxCast Chemicals: " << ToxCast_Chemicals.size() << endl;
            cout << "--------------------------------------------------------------------------------" << endl;                 
            
            // Verify that the training set is using a CAS Registry Number (pretty sure in the data guide)
            int nocas = 0;            
            for (int i = 0; i < Training_Chemicals.size(); i++) {
                c = Training_Chemicals.at(i);
                in = c.chemical_casrn;
                if (in.find("NOCAS") != string::npos) {
                    nocas++;
                }
            }
            
            cout << "   Number of non-CAS chemicals in the training set: " << nocas << endl;
            cout << "--------------------------------------------------------------------------------" << endl;   

            cout << "Training Chemical SMILES Structures:" << endl;
            cout << "--------------------------------------------------------------------------------" << endl; 

            int count_training_match_dsstox = 0;
            
            for (int i = 0; i < Training_Chemicals.size(); i++) {
                c = Training_Chemicals.at(i);

                for (int k = 0; k < DSS_Chemicals.size(); k++) {
                    dc = DSS_Chemicals.at(k);
                    if (dc.DSSTox_GSID == c.chemical_short_gs_id) {
                        c.chemical_cid = dc.DSSTox_CID;
                        istringstream(dc.STRUCTURE_MW) >> c.chemical_Mw;
                        c.chemical_SMILES = dc.STRUCTURE_SMILES;
                        Training_Chemicals.at(i) = c;
                        count_training_match_dsstox++;
                    }
                }
                cout << c.chemical_SMILES << endl;
            }

            int count_predict_match_dsstox = 0;
            
            cout << "Predict Chemical SMILES Structures:" << endl;
            cout << "--------------------------------------------------------------------------------" << endl;              
            
            for (int i = 0; i < Predict_Chemicals.size(); i++) {
                c = Predict_Chemicals.at(i);

                for (int k = 0; k < DSS_Chemicals.size(); k++) {
                    dc = DSS_Chemicals.at(k);
                    if (dc.DSSTox_GSID == c.chemical_short_gs_id) {
                        c.chemical_cid = dc.DSSTox_CID;
                        istringstream(dc.STRUCTURE_MW) >> c.chemical_Mw;
                        c.chemical_SMILES = dc.STRUCTURE_SMILES;                        
                        Predict_Chemicals.at(i) = c;
                        count_predict_match_dsstox++;
                    }
                }
                
                cout << c.chemical_SMILES << endl;
            }
            
            vector<unsigned> topology;
            topology.clear();
            topology.push_back(1);  // 1 input:  Mw
            topology.push_back(3); // 3 nodes in the hidden layer
            topology.push_back(1);  // LEL value
            
            cout << "Creating a <" << 
                    topology.at(0) << " " << 
                    topology.at(1) << " " << 
                    topology.at(2) <<  
                    "> neural network" << endl;

            Net myNet(topology);

            vector<double> inputVals, targetVals, resultVals;
            int trainingPass = 0;

            double NN_cycle_start = clock();
            double NN_cycle_stop = NN_cycle_start + (30*1000000.0);
            while (((myNet.getRecentAverageError() > 0.001) || (trainingPass < 5)) && (clock() < NN_cycle_stop)) {
                trainingPass++;
                for (int i = 0; i < Training_Chemicals.size(); i++) {
                    inputVals.clear();
                    targetVals.clear();
                    
                    c = Training_Chemicals.at(i);
                    
                    inputVals.push_back(parabolic_filter((double) c.chemical_Mw, Mw_mean, Mw_min, Mw_max));             
                    targetVals.push_back(parabolic_filter((double) c.systemic_adjusted_negative_log_lel, LEL_mean, LEL_min, LEL_max));

                    myNet.feedForward(inputVals);
                    myNet.getResults(resultVals);

                    assert(targetVals.size() == topology.back()); 
                    myNet.backProp(targetVals); 
                    
                }
            }
            
            cout << "Neural Net Recent Average Error: " << myNet.getRecentAverageError() << endl; 
            cout << "Total Training Passes:  " << trainingPass << endl;
            cout << "--------------------------------------------------------------------------------" << endl;                    

            // finish training calculations
            finishTraining = true;            

            if (finishTraining == true) {
                training_finish_time = clock();
                training_elapsed_time = training_finish_time - training_start_time; 
                cout << "training elapsed time: " << (training_elapsed_time / 1000000.0) << " sec" << endl;       
                finishTraining = true;
                startCalculations = true;
            }
            
            // start calculations
            startCalculations = true;
            if (startCalculations == true) {
                calculation_start_time = clock();        
                startCalculations = false;
            }            
                
            // Calculate the LEL value using the trained NN
            for (int i = 0; i < Predict_Chemicals.size(); i++) {
                inputVals.clear();
                c = Predict_Chemicals.at(i);
                
                inputVals.push_back(parabolic_filter((double)c.chemical_Mw, Mw_mean, Mw_min, Mw_max));

                myNet.feedForward(inputVals);
                myNet.getResults(resultVals);
                
                c.systemic_adjusted_negative_log_lel = parabolic_amplifier(resultVals.at(0), LEL_mean, LEL_min, LEL_max);
                
                Predicted_LEL_values.push_back(c.systemic_adjusted_negative_log_lel);
                Predict_Chemicals.at(i) = c;            
            }     
            
            cout << "--------------------------------------------------------------------------------" << endl;              
     
            LEL_min = min(Predicted_LEL_values);
            LEL_mean = average(Predicted_LEL_values);
            LEL_max = max(Predicted_LEL_values);
            
            cout << "Calculation / Prediction Set Statistics:" << endl;
            cout << "   Min LEL Value: " << LEL_min << endl;
            cout << "   Max LEL Value: " << LEL_max << endl;
            cout << "   Average LEL Value: " << LEL_mean << endl;
            cout << "   Sigma^2 of LEL Value: " << stddev2(Predicted_LEL_values) << endl;
            cout << "--------------------------------------------------------------------------------" << endl;             

            // calculate the answers based on training
            finishCalculations = true;           
            
            cout << "class LELPredictor { " << endl;
            cout << "public:" << endl << "  vector<double> predictLEL(void);" << endl;
            cout << "};" << endl;
            
            cout << "vector<double> LELPredictor::predictLEL(void) {" << endl;
            cout << "   // Output the LEL values" << endl;
            cout << "   vector<double> lel;" << endl << endl;

            for (int i = 0; i < Initial_Chemicals.size(); i++) {
                c = Initial_Chemicals.at(i);
                
                for (int j = 0; j < Predict_Chemicals.size(); j++) {
                    pc = Predict_Chemicals.at(j);
                    
                    if(c.chemical_gs_id == pc.chemical_gs_id) {
                        c.systemic_adjusted_negative_log_lel = pc.systemic_adjusted_negative_log_lel;
                        Initial_Chemicals.at(i) = c;
                    }
                }
                
                cout << "   lel.push_back(" << c.systemic_adjusted_negative_log_lel << ");" << endl;                
            }
            cout << "   return lel;" << endl;
            cout << "}" << endl << endl;
            
            cout << "--------------------------------------------------------------------------------" << endl;               

            if (finishCalculations == true) {
                calculation_finish_time = clock();
                calculation_elapsed_time = calculation_finish_time - training_start_time;                 
                cout << "total elapsed time: " << (calculation_elapsed_time / 1000000.0) << " sec" << endl; 
                cout << "================================================================================" << endl;                
                finishCalculations = false;
            }
        }
    }
    return 0;
}
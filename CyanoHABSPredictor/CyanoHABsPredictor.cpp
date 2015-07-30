/* =========================================================================================================== *
 * EPA Cyano Modeling Challenge | TopCoder 
 * CyanoHABsPredictor [Competition Sensitive]
 * RigelFive - Hudson Ohio USA
 * 10 June 2014 - 28 June 2014
 * =========================================================================================================== *
 * version history:
 * v1-v6:  Establishing of a baseline score using the average of the detected cyanobacteria.
 * v7:  Approach using a <3 4 1> NN with inputs of time, lat, lng and the output as the level of cyanobacteria.
 * v8:  Estimate levels of cyanobacteria using time as an input and the NN topology structured to the map.
 * =========================================================================================================== */

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <string>

using namespace std;

struct AreaInfo_box { 
    string region_name;
    int area_id;
    double lat_min;
    double lat_max;
    double long_min;
    double long_max;
    bool active; 
    double lat; 
    double lng;  
    int latID;
    int lngID;
};
struct areaMerisData_box {
    int year;
    int day;
    int areaID;
    double calculated_cyanobacteria;  
    double lat;
    double lng;
};
struct tessarect {  
    int dateID; 
    double lat;
    double lng;   
    int latID;
    int lngID;   
    double cyanobacteria;     // units:  cells/mL
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
class Analysis {
public:
    double sum(vector<double> v);
    double average(vector<double> v);
    double stddev2(vector<double> v);
    double min(vector<double> v);
    double max(vector<double> v);
    
    double line(double m, double b, double x);
    double parabola(double a, double b, double c, double x);
    double filter(double x, double x_mean, double x_min, double x_max);
    double amplifier(double y, double x_mean, double x_min, double x_max);  

    double distance(double x1, double y1, double x2, double y2); 
    
private:
    
};
double Analysis::sum(vector<double> v) {
    double sum = 0.0;
    
    int vSIZE = v.size();     
    for (int i = 0; i < vSIZE; i++) {
        sum = sum + double(v.at(i));
    }
    return double(sum);
}
double Analysis::average(vector<double> v) {  
    double sum = 0.0;
    
    int vSIZE = v.size();     
    for (int i = 0; i < vSIZE; i++) {
        sum = sum + abs(v.at(i));
    }
    return double(sum/vSIZE);
}
double Analysis::stddev2(vector<double> v) {
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
double Analysis::min(vector<double> v) {
    double min1 = 1e+09;
    int vSIZE = v.size();
    
    for (int i = 0; i < vSIZE; i++) {
        if (v.at(i) < min1) { min1 = v.at(i); };
    }

    return double(min1);    
}
double Analysis::max(vector<double> v) {
    double max1 = 0;

    int vSIZE = v.size();
    
    for (int i = 0; i < vSIZE; i++) {
        if (v.at(i) > max1) { max1 = v.at(i); };
    }

    return double(max1);
}
double Analysis::line(double m, double b, double x) {
    double y;

    y = (m*x) + b;
    
    return y;
}
double Analysis::parabola(double a, double b, double c, double x) {
    double y;
    
    y = (a*pow((x-b),2.0)) + c;
    
    return y;
}
double Analysis::filter(double x, double x_mean, double x_min, double x_max) {
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
double Analysis::amplifier(double y, double x_mean, double x_min, double x_max) { 
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
double Analysis::distance(double x1, double y1, double x2, double y2) {
    double radius;
    radius = 6371.0;
    double pi = 3.14159265;
    x1=(x1/180.0)*pi;
    y1=(y1/180.0)*pi;
    x2=(x2/180.0)*pi;
    y2=(y2/180.0)*pi;
    if (x1==x2 && y1==y2)
        return 0;
    else
    {
        if ((sin(x2)*sin(x1)+cos(x2)*cos(x1)*cos(y2-y1))>1) 
            return radius*acos(1);
        else
            return radius*acos(sin(x2)*sin(x1)+cos(x2)*cos(x1)*cos(y2-y1));
    }
}
class CyanoHABsPredictor {
public:
    int init(string region, double minVInterested, double maxVInterested, int timeSpan, vector<string> &areaInfo);
    int receiveData(int dayOfYear, int year, vector<string> &merisData, vector<string> &areaMerisData, vector<string> &globalClimate, vector<string> &localClimate, vector<string> &waterQuality, vector<string> &cropScape);
    vector<double> predict(int dayOfYear, int year, int timeSpan, vector<int> &requests);    
    void build_B1(string region, vector<string> &areaInfo);
    void build_B3(vector<string> &areaMerisData);
    void build_T1(int dayOfYear, int year);
    void train_T1(void);
    vector<double> feed_T1(int year, int day, int timeSpan, int areaID);
    
    int getDateID(int year, int day);
    vector<double> getLatLong (int areaID);
    double constrainCyano(double value);
    int findClosestPoint(double lat, double lng); 
    string modify_csv(string in);
    double average(vector<double> v);
    void ten_vector_output(vector<double> &vector_set);
    void ten_vector_int_output(vector<int> &vector_set);
    
    // NN methods
    void buildNNtopology(int input, int middle, int output);    
    void feedForward(const vector<double> &inputVals);
    void backProp(const vector<double> &targetVals);
    void getResults(vector<double> &resultVals) const;
    double getRecentAverageError(void) const { return m_recentAverageError; }
    void testingNNalgorithm(void);    

    AreaInfo_box b1;
    areaMerisData_box b3;
    tessarect t1;
    
    vector<AreaInfo_box> B1;
    vector<areaMerisData_box> B3;
    vector<tessarect> T1;
    vector<unsigned> topology;
    vector<double> inputVals, targetVals, resultVals;       
    
    double min_pop;
    double max_pop;
    double min_pop_val = 10232.9;
    double max_pop_val = 3090295.4;
        
    int T1_match = 0;
    int T1_miss = 0;
    int trainingPass = 0;
    
    int GridSize = 16; // total number of inputs for the NN
    int GridWidth = 4;  // square root of the GridSize (i.e. 4x4 grid = 16 grids)
    
private:
    vector<Layer> m_layers;
    double m_error;
    double m_recentAverageError;
    static double m_recentAverageSmoothingFactor;
    
};
void CyanoHABsPredictor::buildNNtopology(int input, int middle, int output) {
    topology.clear();
    topology.push_back(input);
    topology.push_back(middle);
    topology.push_back(output);  
    
    unsigned numLayers = topology.size();
    for (unsigned layerNum = 0; layerNum < numLayers; ++layerNum) {
        m_layers.push_back(Layer());
        unsigned numOutputs = layerNum == topology.size() - 1 ? 0 : topology[layerNum + 1];

        for (unsigned neuronNum = 0; neuronNum <= topology[layerNum]; ++neuronNum) {
            m_layers.back().push_back(Neuron(numOutputs, neuronNum));
        }
        m_layers.back().back().setOutputVal(1.0);
    }  
}
double CyanoHABsPredictor::m_recentAverageSmoothingFactor = 100.0; // Number of training samples to average over
void CyanoHABsPredictor::getResults(vector<double> &resultVals) const {
    resultVals.clear();

    for (unsigned n = 0; n < m_layers.back().size() - 1; ++n) {
        resultVals.push_back(m_layers.back()[n].getOutputVal());
    }
}
void CyanoHABsPredictor::backProp(const vector<double> &targetVals) {
    Layer &outputLayer = m_layers.back();
    m_error = 0.0;

    for (unsigned n = 0; n < outputLayer.size() - 1; ++n) {
        double delta = targetVals[n] - outputLayer[n].getOutputVal();
        m_error += delta * delta;
    }
    m_error /= outputLayer.size() - 1;
    m_error = sqrt(m_error); 


    m_recentAverageError =
            (m_recentAverageError * m_recentAverageSmoothingFactor + m_error)
            / (m_recentAverageSmoothingFactor + 1.0);

    for (unsigned n = 0; n < outputLayer.size() - 1; ++n) {
        outputLayer[n].calcOutputGradients(targetVals[n]);
    }

    for (unsigned layerNum = m_layers.size() - 2; layerNum > 0; --layerNum) {
        Layer &hiddenLayer = m_layers[layerNum];
        Layer &nextLayer = m_layers[layerNum + 1];

        for (unsigned n = 0; n < hiddenLayer.size(); ++n) {
            hiddenLayer[n].calcHiddenGradients(nextLayer);
        }
    }

    for (unsigned layerNum = m_layers.size() - 1; layerNum > 0; --layerNum) {
        Layer &layer = m_layers[layerNum];
        Layer &prevLayer = m_layers[layerNum - 1];

        for (unsigned n = 0; n < layer.size() - 1; ++n) {
            layer[n].updateInputWeights(prevLayer);
        }
    }
}
void CyanoHABsPredictor::feedForward(const vector<double> &inputVals) {
    assert(inputVals.size() == m_layers[0].size() - 1);

    for (unsigned i = 0; i < inputVals.size(); ++i) {
        m_layers[0][i].setOutputVal(inputVals[i]);
    }

    for (unsigned layerNum = 1; layerNum < m_layers.size(); ++layerNum) {
        Layer &prevLayer = m_layers[layerNum - 1];
        for (unsigned n = 0; n < m_layers[layerNum].size() - 1; ++n) {
            m_layers[layerNum][n].feedForward(prevLayer);
        }
    }
}
void CyanoHABsPredictor::build_B1(string region, vector<string> &areaInfo) {
    int i_val, div;
    vector<double> lat;
    vector<double> lng;
    double max_lat, min_lat, max_lng, min_lng;
    
    Analysis a1;
        
    for (int i = 1; i < areaInfo.size(); i++) {
        b1.region_name = region;   
        stringstream ss(modify_csv(areaInfo.at(i)));
        ss >> b1.area_id >> b1.lat_min >> b1.lat_max >> b1.long_min >> b1.long_max >> i_val;
        
        if (i_val == 1) {
            b1.active = true;
        }
        else {
            b1.active = false;
        }
        
        b1.lat = 0.5*(b1.lat_max + b1.lat_min);
        b1.lng = 0.5*(b1.long_max + b1.long_min);
        b1.latID = 0;
        b1.lngID = 0;

        B1.push_back(b1);       
    }
    
    for (int i = 0; i < B1.size(); i++) {
        b1 = B1.at(i);
        lat.push_back(b1.lat);
        lng.push_back(b1.lng);
    }
    
    min_lat = a1.min(lat);
    max_lat = a1.max(lat);
    min_lng = a1.min(lng);
    max_lng = a1.max(lng);
    
    div = GridWidth;
    
    for (int i = 0; i < B1.size(); i++) {
        b1 = B1.at(i);
        b1.latID = (int) round((b1.lat / (max_lat-min_lat)) * (double) div);
        b1.lngID = (int) round((b1.lng / (max_lng-min_lng)) * (double) div);
        B1.at(i) = b1;
    }
}
void CyanoHABsPredictor::build_B3(vector<string> &areaMerisData) {
    for (int i = 0; i < areaMerisData.size(); i++) {
        stringstream ss(modify_csv(areaMerisData.at(i)));
        ss >> b3.year >> b3.day >> b3.areaID >> b3.calculated_cyanobacteria;
        
        for (int j = 0; j < B1.size(); j++) {
            b1 = B1.at(j);
            if (b1.area_id == b3.areaID) {
                b3.lat = b1.lat;
                b3.lng = b1.lng;
            }
        }
        
        B3.push_back(b3);       
    }      
}
int CyanoHABsPredictor::getDateID(int year, int day) {
    int dateID;
    
    dateID = (year-2009)*365 + day;
    return dateID;
}
vector<double> CyanoHABsPredictor::getLatLong (int areaID) {
    double lat;
    double lng;
    vector<double> latlong;
    
    for (int i = 0; i < B1.size(); i++) {
        b1 = B1.at(i);
        if (b1.area_id == areaID) {
            lat = 0.5 * (b1.lat_max + b1.lat_min);
            lng = 0.5 * (b1.long_max + b1.long_min);
        }
    }
    
    latlong.push_back(lat);
    latlong.push_back(lng);
    return latlong;
}
int CyanoHABsPredictor::findClosestPoint(double lat, double lng) {
    Analysis a1;
    
    double dist, closest_dist; 
    int index = 0;
    
    closest_dist = 9999999999999.0;
    
    for (int i = 0; i < B3.size(); i++) {
        b3 = B3.at(i);    
        dist = a1.distance(lat, lng, b3.lat, b3.lng);
        if (dist <= closest_dist) {
            closest_dist = dist;
            index = i;
        }
    }
             
    return index;
}
void CyanoHABsPredictor::build_T1(int dayOfYear, int year) {
    if (B3.size() > 0) {
        t1.dateID = getDateID(year, dayOfYear);
        
        for (int i = 0; i < B3.size(); i++) {
            t1.lat = b3.lat;
            t1.lng = b3.lng;
            b1 = B1.at(findClosestPoint(b3.lat, b3.lng));
            t1.latID = b1.latID;
            t1.lngID = b1.lngID;
            t1.cyanobacteria = b3.calculated_cyanobacteria;
            T1.push_back(t1);
        }
    }
}
void CyanoHABsPredictor::train_T1(void) {

    Analysis a1;
    
    int NN_training_entry_time = clock();   
    int training_time = 1000000; 
    int NN_training_exit_time = NN_training_entry_time + training_time;    
    int sector;
    
    while (clock() < NN_training_exit_time) {
        for (int i = 0; i < T1.size(); i++) {
            t1 = T1.at(i);
            
            inputVals.clear();
            targetVals.clear();
            resultVals.clear();   
            
            for (int i = 0; i < GridSize; i++) {
                sector = (t1.latID*GridWidth) + t1.lngID;
                if (sector == i) {
                    inputVals.push_back(a1.parabola((double) t1.dateID, 638, 0, 1275));
                } else {
                    inputVals.push_back(0.0);
                }
            }            
            
            if (t1.cyanobacteria < 200000.0) { 
                targetVals.push_back(a1.parabola((double) t1.cyanobacteria, 105116.45, 10232.9, 200000.0));
                targetVals.push_back(0);  
                targetVals.push_back(0);  
                targetVals.push_back(0);  
            } else if ((t1.cyanobacteria >= 200000.0) && (t1.cyanobacteria < 300000.0)) { 
                targetVals.push_back(0);
                targetVals.push_back(a1.parabola((double) t1.cyanobacteria, 200000.0, 100000.0, 300000.0));  
                targetVals.push_back(0);  
                targetVals.push_back(0);  
            } else if ((t1.cyanobacteria >= 300000.0) && (t1.cyanobacteria < 1000000.0)) { 
                targetVals.push_back(0);
                targetVals.push_back(0);  
                targetVals.push_back(a1.parabola((double) t1.cyanobacteria, 650000.0, 300000.0, 1000000.0));  
                targetVals.push_back(0);         
            } else if (t1.cyanobacteria >= 1000000.0) {
                targetVals.push_back(0);
                targetVals.push_back(0);  
                targetVals.push_back(0);  
                targetVals.push_back(a1.parabola((double) t1.cyanobacteria, 2045147.7, 1000000.0, 3090295.4));     
            }
            
            feedForward(inputVals);
            getResults(resultVals);
            assert(targetVals.size() == topology.back());
            backProp(targetVals); 
        }
        trainingPass++;        
    }
}
vector<double> CyanoHABsPredictor::feed_T1(int year, int day, int timeSpan, int areaID) {
    inputVals.clear();
    resultVals.clear(); 

    int dateID, sector;
    double lat, lng;
    int latID, lngID;
    
    Analysis a1;
    dateID = getDateID(year, day) + timeSpan;      
    for (int i = 0; i < B1.size(); i++) {
        b1 = B1.at(i);
        if (b1.area_id == areaID) {
            lat = 0.5 * (b1.lat_max + b1.lat_min);
            lng = 0.5 * (b1.long_max + b1.long_min);
            latID = b1.latID;
            lngID = b1.lngID;            
        }
    } 
    for (int i = 0; i < GridSize; i++) {
        sector = (latID*GridWidth) + lngID;
        if (sector == i) {
            inputVals.push_back(a1.parabola((double) dateID, 638, 0, 1275));
        } else {
            inputVals.push_back(0.0);
        }
    }
        
    feedForward(inputVals);
    getResults(resultVals);  

    return resultVals;
}
string CyanoHABsPredictor::modify_csv(string in) {
    string out = in;
    replace(out.begin(), out.end(), ',', ' ');
    return out;
}
double CyanoHABsPredictor::average(vector<double> v) {  
    double sum = 0.0;
  
    int vSIZE = v.size();
    for (int i = 0; i < vSIZE; i++) {
        sum = sum + abs(v.at(i));
    }
    return double(sum/vSIZE);
}
void CyanoHABsPredictor::ten_vector_output(vector<double> &vector_set) {
    for (int i = 0; i < vector_set.size(); i++) {
        cout << " " << vector_set.at(i);
    }
    cout << endl;
}
void CyanoHABsPredictor::ten_vector_int_output(vector<int> &vector_set) {
    for (int i = 0; i < vector_set.size(); i++) {
        cout << " " << vector_set.at(i);
    }
    cout << endl;
}
double CyanoHABsPredictor::constrainCyano(double value) {
    double constrained_value;
    
    if (constrained_value < min_pop_val) {
        constrained_value = min_pop_val;
    }
    else if (constrained_value > max_pop_val) {
        constrained_value = max_pop_val;
    }    
    
    return constrained_value;
}
int CyanoHABsPredictor::init(string region, double minVInterested, double maxVInterested, int timeSpan, vector<string> &areaInfo) {
    
    int inputs, hidden, outputs;
    
    min_pop = minVInterested;
    max_pop = maxVInterested;
    
    B1.clear();
    B3.clear();
    T1.clear();
    
    build_B1(region, areaInfo); 
    
    outputs = 4;    
    inputs = GridSize;
    hidden = inputs+outputs;

    if (m_layers.size() != 3) {
        cout << "Building a <" << inputs << " " << hidden << " " << outputs << "> Neural Network..." << endl;         
        buildNNtopology(inputs,hidden,outputs);
    }   
    
    return 0;
}
int CyanoHABsPredictor::receiveData(int dayOfYear, int year, vector<string> &merisData, 
    vector<string> &areaMerisData, vector<string> &globalClimate, vector<string> &localClimate, 
    vector<string> &waterQuality, vector<string> &cropScape) {

    int max_recycle = 10;
    int recycle = 0;
    
    B3.clear();
    
    if (areaMerisData.size() > 0) {
        build_B3(areaMerisData);
        build_T1(dayOfYear, year); 
        train_T1();
        
        while ((m_recentAverageError > 0.10 ) && (recycle < max_recycle)) {
            train_T1(); 
            recycle++;
        }
    }

    return 0;
}
vector<double> CyanoHABsPredictor::predict(int dayOfYear, int year, int timeSpan, vector<int> &requests) {
    vector<double> the_big_answer;
    Analysis a1;
  
    int areaID;
    vector<double> cyano;
    double c0, c1, c2, c3;
   
    the_big_answer.clear();  
  
    for (int i = 0; i < requests.size(); i++) {
        areaID = requests.at(i);
      	if (T1.size() > 0) {    
            cyano.clear();
            cyano = feed_T1(year, dayOfYear, timeSpan, areaID);
            
            c0 = cyano.at(0);
            c1 = cyano.at(1);
            c2 = cyano.at(2);
            c3 = cyano.at(3);
            
            if ((c0 > c1) && (c0 > c2) && (c0 > c3)) {
                the_big_answer.push_back(constrainCyano(a1.amplifier(c0, 105116.45, 10239.2, 200000.0)));                
            } else if ((c1 > c0) && (c1 > c2) && (c1 > c3)) {
                the_big_answer.push_back(constrainCyano(a1.amplifier(c1, 250000.0, 200000.0, 300000.0)));                                
            } else if ((c2 > c0) && (c2 > c1) && (c2 > c3)) {
                the_big_answer.push_back(constrainCyano(a1.amplifier(c2, 650000.0, 300000.0, 1000000.0)));                                
            } else if ((c3 > c0) && (c3 > c1) && (c3 > c2)) {
                the_big_answer.push_back(constrainCyano(a1.amplifier(c3, 2045147.7, 1000000.0, 3090295.4)));                                
            }

            T1_match++;
      	}
      	else {
            the_big_answer.push_back(min_pop_val);            
            T1_miss++;
        }       
    }

    cout << "Requests:  " << requests.size()
            << "\tB1 size: " << B1.size()
            << "\tB3 size: " << B3.size()
            << "\tT1 size: " << T1.size()            
            << "\tT1 match/miss:  " << T1_match
            << "/" << T1_miss
            << "\tNN error:  " << m_recentAverageError
            << "\tNN-train passes: " << trainingPass
            << endl;

    return the_big_answer;
}
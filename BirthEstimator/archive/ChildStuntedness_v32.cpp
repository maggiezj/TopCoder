/* ChildStuntedness v31 [Competition Sensitive]
 * TopCoder
 * RigelFive - Hudson Ohio USA
 * 11 August 2014 - 21 August 2014
 */

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <tuple>
#include <cassert>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <string>
#include <bitset>
#include <random>

using namespace std;

typedef tuple<int,int> object_tuple;
typedef tuple<int,double> GA_tuple;

bool sort_object_score (const object_tuple &lhs, const object_tuple &rhs){
  return get<1>(lhs) < get<1>(rhs);  // sorts for smallest numbers at beginning of list (output of result to main)
}

bool sort_GA (const object_tuple &lhs, const object_tuple &rhs){
  return get<1>(lhs) > get<1>(rhs);  // sorts for largest numbers at beginning of list (output of result to main)
}

struct _odv {
    bool isNA;
    double odv_value;
    string str;
};
struct _box { // a box is a collection of the data points (cases) for a single parameter (age, ODV5, etc)
    double min;
    double max;
    double mean;
    vector<double> values;
    vector<double> values_no_zero;
    vector<int> category;
    string name;
};
struct _case { // a case is a single data point with a collection of parameters
    int ID;             // 1. Unique Fetus ID
    double age;         // 2. Estimated fetus gestational age from last menstrual recall date
    double sex;         // 3. 0 = Male, 1 = Female
    double nutrition;   // 4. Maternal nutritional status (1 or 2)
    vector<_odv> odv;   // 5-12. Dependent variables: Ultrasound observed measurements
    double weight;      // 13. Birth Weight (w)
    double duration;    // 14. Pregnancy Duration, or Birthday (b)
};
struct _cases { // cases are a collection of data points sorted by the type based on gender and nutrition
    vector<_case> typeI; // sex = male, nutrition = 1
    vector<_case> typeII; // sex = male, nutrition = 2
    vector<_case> typeIII; // sex = female, nutrition = 1
    vector<_case> typeIV; // sex = female, nutrition = 2
    vector<_box> T1; 
    vector<_box> T2; 
    vector<_box> T3; 
    vector<_box> T4;
};
struct _train {
    int ID;
    int category;
    double weight;
    double duration;
};
struct _result {
    int ID;
    double weight;
    double duration;
};
struct Connection {
    double weight;
    double deltaWeight;
};
struct _gn { // structure to build a variable length genome
    int gene_id; // id number for the gene in the genome
    bitset<64> gene;
    double min_value;
    double max_value;
};
struct _gn2 { // structure to build a genome pair
    vector<_gn> genome_1;
    vector<_gn> genome_2;   
};
struct _pop { // structure for a population of genomes
    int id_num;
    int gen_num;
    vector<_gn> genome;
    double fitness;
};
struct _coeff { // structure for wN and dN coefficients    
    vector<double> N; // vector of coefficients. 
    double N_max; // maximum limit for coefficients in GAsearch
    double N_min; // minimum limit for coefficients in GAsearch
};
struct _param {
    double age;
    double ODV5;
    double ODV6;
    double ODV7;
    double ODV8;
    double ODV9;
    double ODV10;
    double ODV11;
    double ODV12;
};
struct _matrix {
    int row;
    int col;
    double val;
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
/*
 * Analysis class
 */
class Analysis {
public:
    double sum(vector<double> v);
    double average(vector<double> v);
    double stddev2(vector<double> v);
    double min(vector<double> v);
    double max(vector<double> v);
    int max_int(vector<int> v);
    
    double line(double m, double b, double x);
    double parabola(double a, double b, double c, double x);
    double filter(double x, double x_mean, double x_min, double x_max);
    double amplifier(double y, double x_mean, double x_min, double x_max);  

    double distance(double x1, double y1, double x2, double y2); 
    vector<int> category(vector<double> v, vector<double> val);
    vector<int> unique(vector<int> list);
    double interpolate(vector<double> x, vector<double> N, double N_0);
    vector<_matrix> covariance(vector<double> x, vector<double> y);
    vector<_matrix> inverse_2x2(vector<_matrix> A);
    
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
    double max1 = -1e20;

    int vSIZE = v.size();
    
    for (int i = 0; i < vSIZE; i++) {
        if (v.at(i) > max1) { max1 = v.at(i); };
    }

    return double(max1);
}
int Analysis::max_int(vector<int> v) {
    int max1 = -1e9;

    int vSIZE = v.size();
    
    for (int i = 0; i < vSIZE; i++) {
        if (v.at(i) > max1) { max1 = v.at(i); };
    }

    return max1;
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
vector<int> Analysis::category(vector<double> v, vector<double> vals) {
    // determines if a value can be categorized as ultra-high, high, low, ultra-low, or average
    // uses 0.674 of the std deviation to partition the data into thirds
    // the average group has 50% of the data within its partition.
    // the ultra group is split for the top 10% and bottom 10%
    // ultra-low = -2, low = -1, average = 0, high = 1, ultra-high = 2

    double stddev = sqrt(stddev2(v));
    double mean = average(v);
    double interval_50pct = 0.950 * stddev; // 0.674490 // 0.950
    double interval_20pct = 1.309 * stddev; // 1.281552 // 1.309
    vector<int> cats;
    int cat = 0;
    double val;
    
    for (int i = 0; i < vals.size(); i++) {
        val = vals.at(i);
        if (val < (mean - interval_20pct)) {
            cat = -2;
        } else if ((val >= (mean - interval_20pct)) && (val < (mean - interval_50pct))) {
            cat = -1;
        } else if ((val >= (mean - interval_50pct)) && (val <= (mean + interval_50pct))) {
            cat = 0;
        } else if ((val > (mean + interval_50pct)) && (val <= (mean + interval_20pct))) {
            cat = 1;
        } else if ((val > (mean + interval_20pct))) {
            cat = 2;
        }
        cats.push_back(cat);
    }

    return cats;
}
vector<int> Analysis::unique(vector<int> list) {
    vector<int> unique_list;
    bool inList;
    for (int i = 0; i < list.size(); i++) {
        inList = false;
        for (int j = 0; j < unique_list.size(); j++) {
            if (unique_list.at(j) == list.at(i)) {
                inList = true;

            }
        }
        if (!inList) {
            unique_list.push_back(list.at(i));
        }
    }
    
    return unique_list;

}
double Analysis::interpolate(vector<double> x, vector<double> N, double N_0) {
    double y = 0;
    
    for (int i = 0; i < x.size(); i++) {
        y += (N.at(i) * x.at(i));
    }

    return y + N_0;
}
vector<_matrix> Analysis::covariance(vector<double> x, vector<double> y) {
    _matrix m;
    vector<_matrix> S;
    
    double x_mean = average(x);
    double y_mean = average(y);
    double x1 = 0;
    double y1 = 0;
    double x1y1 = 0;
    double y1x1 = 0;
    
    //cout << "x: " << x.at(0) << "\ty: " << y.at(0) << endl;
    
    for (int i = 0; i < x.size(); i++) {
        x1 += pow((x.at(i) - x_mean),2.0);
    }

    for (int i = 0; i < y.size(); i++) {
        y1 += pow((y.at(i) - y_mean),2.0);
    }
    
    for (int i = 0; i < x.size(); i++) {
        x1y1 += (x.at(i) - x_mean)*(y.at(i) - y_mean);
    }

    y1x1 = x1y1;

    m.row = 0;
    m.col = 0;
    m.val = x1/x.size();
    S.push_back(m);
    
    //cout << "S [" << m.row << "][" << m.col << "]: " <<  m.val << endl;
    
    m.row = 0;
    m.col = 1;
    m.val = x1y1/x.size();
    S.push_back(m);
    //cout << "S [" << m.row << "][" << m.col << "]: " <<  m.val << endl;
    
    m.row = 1;
    m.col = 0;
    m.val = y1x1/x.size();
    S.push_back(m);
    //cout << "S [" << m.row << "][" << m.col << "]: " <<  m.val << endl;
    
    m.row = 1;
    m.col = 1;
    m.val = y1/x.size();
    S.push_back(m);
    //cout << "S [" << m.row << "][" << m.col << "]: " <<  m.val << endl;

    return S;
}
vector<_matrix> Analysis::inverse_2x2(vector<_matrix> A) {
    vector<_matrix> A_1;
    
    _matrix m;
    
    double a, b, c, d;
    
    a = A.at(0).val;
    b = A.at(1).val;
    c = A.at(2).val;
    d = A.at(3).val;
    
    double adbc = ((a*d) - (b*c));
    
    m.row = 0;
    m.col = 0;
    m.val = d/adbc;
    
    A_1.push_back(m);
    
    cout << "--------------------------------------------------------------------------------" << endl;   
    cout << "S-1 [" << m.row << "][" << m.col << "]: " <<  m.val << endl;
    
    m.row = 0;
    m.col = 1;
    m.val = -b/adbc;
    
    A_1.push_back(m);
    cout << "S-1 [" << m.row << "][" << m.col << "]: " <<  m.val << endl;
    
    m.row = 1;
    m.col = 0;
    m.val = -c/adbc;
    
    A_1.push_back(m);   
    cout << "S-1 [" << m.row << "][" << m.col << "]: " <<  m.val << endl;
    
    m.row = 1;
    m.col = 1;
    m.val = a/adbc;
    
    A_1.push_back(m); 
    cout << "S-1 [" << m.row << "][" << m.col << "]: " <<  m.val << endl;
    
    return A_1;
}
/*
 * ChildStuntedness class
 */
class ChildStuntedness {
public:
    vector<double> predict(vector<string> training, vector<string> testing);
    
    _odv doesStringContainNA(string str); 
    void parseData (vector<string> &dataString, bool isTraining);
    vector<_box> quantifyParameters(vector<_case> Cs);
    int generateCategoryCode(vector<int> codes);
    void loadParameters(vector<_case> type);
    void clearParameters(void);
    vector<_box> buildCases(void);
    vector<int> buildCategories(vector<_box> T);
    vector<int> compareCategories(vector<int> trainingCat, vector<int> testingCat);
    _box loadBox(string name, vector<double> element, bool allowZero);
    vector<_train> buildCategorySet(vector<_case> Cs, bool isTraining);
    vector<_result> lookupOptimizationTable(void);
    vector<_train> generateOptimizationTable(void);
    vector<double> generateResult(vector<_result> resultVector);
    void calcResultStats(void);
    double calcScore(vector<double> w, vector<double> w_real, vector<double> d, vector<double> d_real);
    vector<double> vectorWithoutZeros(vector<double> v);   
    vector<double> makeParameterList(vector<_box> L, int index);

    // NN methods
    void trainNN(double trainingDuration);
    vector<_result> execNN(void);
    void buildNNtopology(int input, int middle, int output);    
    void feedForward(const vector<double> &inputVals);
    void backProp(const vector<double> &targetVals);
    void getResults(vector<double> &resultVals) const;
    double getRecentAverageError(void) const { return m_recentAverageError; }
    void testingNNalgorithm(void);    

    // Genetic Algorithm Methods
    void generate_optimization(void);
    double calculateFitness(int cat, double w, double d);    
    double optimizationError(double targetVal, double currentVal);     
    vector<_pop> create_population(int numW, int numD, int population_size);
    vector<_pop> compete_population(vector<_pop> population, int cat);
    vector<_pop> select_population(vector<_pop> population);
    vector<_pop> evolve_population(vector<_pop> population, int mutate);
    vector<double> create_values(double min, double max, int values);    
    bitset<64> create_64b_gene(double v, double min, double max);
    double create_64b_value(bitset<64> g, double min, double max);
    vector<_gn> mutate_genome(vector<_gn> old_genome, int mutate);
    _gn2 cross_genomes(_gn2 &old_genome_pair);
    void output_genome(vector<_gn> genome); 
    vector<_result> lookupGAOptimizationTable(void);
    
    // Genetic Algorithm #2 Methods
    void generate_optimization_2(_cases C, int mutate, int popSize, int maxCycle);
    vector<double> cycle_population(vector<_case> type, vector<_box> X, int mutate, int popSize, int maxCycle, double min, double max, bool isW);
    vector<double> calculateParameter(vector<double> N, vector<_box> Tx);
    vector<_pop> create_population_2(int population_size, double min, double max);
    vector<_pop> compete_population_2(vector<_pop> population, vector<_train> tX, bool isW);
    double calculateFitness_2(vector<_case> type, int ID, double w, double d);
    void outputCoefficients(vector<_gn> genome, bool isW);
    vector<_box> convertCaseToSmallBox(vector<_box> Tx);
    vector<_result> lookupCoefficientTable(void);
    
    // Multivariate linear interpolation optimization <genie>
    _cases splitCaseData(vector<_case> Cs);
    vector<_case> vectorType(vector<_case> xC, double sex, double nutrition);  
    
    vector<unsigned> topology;
    vector<double> inputVals, targetVals, resultVals;  
    
    _case hC, nC; // historical (training) case, new (test) case
    vector<_case> hCs, nCs;
    vector<int> trainingCategories, testingCategories;
    
    string str_5, str_6, str_7, str_8, str_9, str_10, str_11, str_12;
    vector<_odv> odv_case;
    string readLine;
    
    _box B;
    vector<_box> T; // for training data
    vector<_box> X; // for testing (experimental) data
    _cases tS; // for training cases
    _cases xS; // for test (experimental) cases    
    
    vector<_train> xM; // testing set with data for all cases/data points
    vector<_train> tM; // training set with data for all cases/data points
    vector<_train> tC; // training set with data for optimized categories
    
    //_result result;
    //vector<_result> results;
    
    vector<string> names;
    
    vector<double> age;
    vector<double> sex;
    vector<double> nutrition;
    vector<double> ODV5;
    vector<double> ODV6;
    vector<double> ODV7;
    vector<double> ODV8;
    vector<double> ODV9;
    vector<double> ODV10;
    vector<double> ODV11;
    vector<double> ODV12;
    vector<double> weight;
    vector<double> duration;
    
    vector<int> unique_trainingCategories, unique_testingCategories;

    vector<_coeff> wNs, dNs; // weight and duration coefficients for linear interpolation    
    vector<_matrix> S, S_1;
    
    vector<double> wN_typeI, wN_typeII, wN_typeIII, wN_typeIV;
    vector<double> dN_typeI, dN_typeII, dN_typeIII, dN_typeIV;
    
private:
    vector<Layer> m_layers;
    double m_error;
    double m_recentAverageError;
    static double m_recentAverageSmoothingFactor;
};

/*
 * Neural Network methods
 */
void ChildStuntedness::trainNN (double trainingDuration) {
    Analysis A;
    _box B;
    
    int trainingPass = 0;
    double val, normVal, min, max, mean;
    double NN_cycle_start = clock();
    double NN_cycle_stop = NN_cycle_start + (trainingDuration*1000000.0);

    //cout << "Sizing hC: " << hCs.size() << endl;
    //cout << "Sizing T: " << T.size() << endl;
    while (clock() < NN_cycle_stop) {
        trainingPass++;
        for (int j = 0; j < hCs.size(); j++) {
            
            // generate inputVals for NN
            //  input 1: gestational age (0..1)
            //  input 2: sex (0,1)
            //  input 3: nutrition (1,2)
            //  input 4: ODV6
            //  input 5: ODV7
            //  input 6: ODV8
            //  input 7: ODV9
            //  input 8: ODV10
            //  input 9: ODV11
            //  input 10: ODV12

            inputVals.clear();
            targetVals.clear();
            
            if (T.at(3).values.at(j) == 0) {
                // the data that occurs with ODV5 = 0 may allow ODVD6-12 to predict w, d
                // otherwise, do not perform a training session on ODV5 > 0 and ODV6-12 = 0
                for (int i = 0; i < 3; i++) { // age sex nutn
                    B = T.at(i);
                    val = B.values.at(j);
                    min = B.min;
                    max = B.max;
                    mean = B.mean;
                    normVal = A.filter(val, mean, min, max);
                    //cout << "input i: " << i << "\tj: " << j << "\tval: " << val << "\tmean: " << mean << "\tmin-max: "<< min << "-" << max << "\tnormval: " << normVal << endl;
                    inputVals.push_back(normVal);
                }
                for (int i = 4; i < 11; i++) { // ODV6-12
                    B = T.at(i);
                    val = B.values.at(j);
                    min = B.min;
                    max = B.max;
                    mean = B.mean;
                    normVal = A.filter(val, mean, min, max);
                    //cout << "input i: " << i << "\tj: " << j << "\tval: " << val << "\tmean: " << mean << "\tmin-max: "<< min << "-" << max << "\tnormval: " << normVal << endl;
                    inputVals.push_back(normVal);
                }        
                
                for (int i = 11; i < 13; i++) { // weight, duration
                    B = T.at(i);
                    val = B.values.at(j);
                    min = B.min;
                    max = B.max;
                    mean = B.mean;
                    normVal = A.filter(val, mean, min, max);
                    //cout << "target i: " << i << "\tj: " << j << "\tval: " << val << "\tmean: " << mean << "\tmin-max: "<< min << "-" << max << "\tnormval: " << normVal << endl;
                    targetVals.push_back(normVal);
                }
                
                //cout << "InputVals: " << inputVals.size() << endl;
                //cout << "TargetVals: " << targetVals.size() << endl;
                feedForward(inputVals);
                getResults(resultVals); 
                backProp(targetVals); 
            }
        }
    }
    
    cout << "--------------------------------------------------------------------------------" << endl;   
    cout << "Neural Net Recent Average Error: " << getRecentAverageError() << endl; 
    cout << "Total NN Training Passes:  " << trainingPass << endl;
    
}
vector<_result> ChildStuntedness::execNN (void) {
    Analysis A;
    double NN_cycle_start = clock();
    double min1, mean1, max1, min2, mean2, max2;
    bool output_nC = false;
    vector<double> weights, durations;
    
    _result rNN;
    vector<_result> rNNs;

    for (int j = 0; j < nCs.size(); j++) {
        inputVals.clear();
        resultVals.clear();
        targetVals.clear();
        nC = nCs.at(j);
        
        if (nC.odv.at(0).odv_value == 0) {

            inputVals.push_back(nC.age);
            inputVals.push_back((double) nC.sex);
            inputVals.push_back((double) nC.nutrition);
            //inputVals.push_back(nC.odv.at(0).odv_value);
            inputVals.push_back(nC.odv.at(1).odv_value);
            inputVals.push_back(nC.odv.at(2).odv_value);
            inputVals.push_back(nC.odv.at(3).odv_value);  
            inputVals.push_back(nC.odv.at(4).odv_value);
            inputVals.push_back(nC.odv.at(5).odv_value);
            inputVals.push_back(nC.odv.at(6).odv_value);
            inputVals.push_back(nC.odv.at(7).odv_value);

            feedForward(inputVals);
            getResults(resultVals); 

            if (output_nC) {
                cout << "nC: " << nC.ID
                        << "\tage: " << nC.age
                        << "\tsex: " << nC.sex
                        << "\tnutrition: " << nC.nutrition
                        //<< "\todv5: " << nC.odv.at(0).odv_value
                        << "\todv6: " << nC.odv.at(1).odv_value
                        << "\todv7: " << nC.odv.at(2).odv_value
                        << "\todv8: " << nC.odv.at(3).odv_value
                        << "\todv9: " << nC.odv.at(4).odv_value
                        << "\todv10: " << nC.odv.at(5).odv_value
                        << "\todv11: " << nC.odv.at(6).odv_value
                        << "\todv12: " << nC.odv.at(7).odv_value
                        << "\tweight: " << resultVals.at(0)
                        << "\tduration: " << resultVals.at(1)
                        << endl;
            }

            _box B1 = T.at(T.size()-2);
            _box B2 = T.at(T.size()-1);

            min1 = B1.min;
            max1 = B1.max;
            mean1 = B1.mean;

            double weight = A.amplifier(resultVals.at(0), mean1, min1, max1);

            min2 = B2.min;
            max2 = B2.max;
            mean2 = B2.mean;
            double duration = A.amplifier(resultVals.at(1), mean2, min2, max2);

            rNN.ID = nC.ID;
            rNN.weight = weight;
            rNN.duration = duration;
            rNNs.push_back(rNN);
        }
    }

    for (int i = 0; i < rNNs.size(); i++) {
        weights.push_back(rNNs.at(i).weight);
        durations.push_back(rNNs.at(i).duration);
    }
    
    min1 = A.min(weights);
    mean1 = A.average(weights);
    max1 = A.max(weights);
    min2 = A.min(durations);
    mean2 = A.average(durations);
    max2 = A.max(durations);
    
    double NN_total_time = (clock() - NN_cycle_start) / 1000000.0;

    cout << "--------------------------------------------------------------------------------" << endl;    
    cout << "NN Testing Parameter Statistics:" << endl;
    cout << "Weight:  \tmin1: " << min1 << "\tmean1: " << mean1 << "\tmax1: " << max1 << endl;
    cout << "Duration:\tmin2: " << min2 << "\tmean2: " << mean2 << "\tmax2: " << max2 << endl;
    cout << "Total Number of Results: " << rNNs.size() << endl;
    cout << "Total Time for Testing:  " << NN_total_time << " sec" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;    
    
    return rNNs;

}
void ChildStuntedness::buildNNtopology(int input, int middle, int output) {
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
double ChildStuntedness::m_recentAverageSmoothingFactor = 100.0; // Number of training samples to average over
void ChildStuntedness::getResults(vector<double> &resultVals) const {
    resultVals.clear();

    for (unsigned n = 0; n < m_layers.back().size() - 1; ++n) {
        resultVals.push_back(m_layers.back()[n].getOutputVal());
    }
}
void ChildStuntedness::backProp(const vector<double> &targetVals) {
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
void ChildStuntedness::feedForward(const vector<double> &inputVals) {
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
/*
 * Generate Optimization #1 - optimization of w,d for each category in data <GA1>
 */
void ChildStuntedness::generate_optimization(void) {
    _pop pop;
    vector<_gn> genome;
    vector<_pop> population, population_competed, population_ranked;

    // Genetic Algorithm Parameters
    int mutate = 10; // number of bits to randomly change in the genome
    int popSize = 100; // number of genomes in the population
    int maxCycle = 100; // number of cycles allowed with no improvement to population occurs
    int m = 1;
    
    for (int cat = 0; cat < unique_trainingCategories.size(); cat++) {
        population.clear();
        population = create_population(1, 1, popSize); 
        double score = -1e9;
        int gen = 0;
        int j = 0;   
        
        while (j < maxCycle) {
        
            gen++;
            population_competed.clear();
            population_ranked.clear();
            
            m = round((double) mutate 
                    * (1.0 - cos(3.141592654 
                    * ((double) gen
                    / (double) popSize))) 
                    + 1.0);
            population_competed = compete_population(population, cat);
            population_ranked = select_population(population_competed);
            population = evolve_population(population_ranked, m);
            j++;
            double newScore = population.at(0).fitness;
            if (score < newScore) {
                score = newScore;
                //cout << "Generation #: " << gen << " Mutate: " << mutate << " Fitness: " << score << endl;
                j=0;
            }
        }
        
        cout << "Category: " << cat 
                << "\tMutate: " << m
                << "\tGen: " << gen 
                << "\tBest Fitness: " << score
                //<< "----------------------------------------------------------------------------------" 
                << endl; 
        pop = population.at(0);
        genome = pop.genome;

        _gn gene_1, gene_2;
        _train tc;

        gene_1 = genome.at(0);
        gene_2 = genome.at(1);

        double w = create_64b_value(gene_1.gene, gene_1.min_value, gene_1.max_value);
        double d = create_64b_value(gene_2.gene, gene_2.min_value, gene_2.max_value);

        tc.ID = cat;
        tc.category = unique_trainingCategories.at(cat);
        tc.weight = w;
        tc.duration = d;
        tC.push_back(tc);
        
        /*
        cout << "Category #: " << cat << endl;
        cout << "w: " << w << endl
            << "d: " << d << endl << endl;*/
    }

}
double ChildStuntedness::calculateFitness(int cat, double w, double d) {
    double error = 0;
    double inverseS_00 = S_1.at(0).val; // 3554.42 for the complete set of data
    double inverseS_01 = S_1.at(1).val; // -328.119 for the complete set of data
    double inverseS_10 = S_1.at(2).val; // -328.119 for the complete set of data
    double inverseS_11 = S_1.at(3).val; // 133.511 for the complete set of data   
    
    // take the average w & d value for the category and calc the error for all data points in that category
    int run_category = unique_trainingCategories.at(cat); // category to calculate score
    for (int i = 0; i < tM.size(); i++) {
        if (tM.at(i).category == run_category) {
            double deltaW = w - tM.at(i).weight;
            double deltaD = d - tM.at(i).duration;
            double val1 = deltaD * inverseS_00 + deltaW * inverseS_10;
            double val2 = deltaD * inverseS_01 + deltaW * inverseS_11;
            double ei = (val1 * deltaD) + (val2 * deltaW);
            error += ei;
        }
    }

    // compare the values for w and d with the ChildStuntednewss score function
    
    //double fitness = 100.0 - abs(error);
    double fitness = -log10(error);

    //cout << "fitness: " << fitness << endl;

    return fitness;
}
double ChildStuntedness::optimizationError(double targetVal, double currentVal) {
    double targError;
    
    targError = abs((targetVal-currentVal)/targetVal) * 100.0;
    return targError;

}
vector<_pop> ChildStuntedness::create_population(int numW, int numD, int population_size) {
// generate a population using with parameters defined to the objective function / competition
    
    // for each case:
    //     1) weight from 0-1
    //     2) duration from 0-1

    _gn gene1;
    _gn gene2;

    vector<_gn> genome;
    _pop pop;
    vector<_pop> population;    
    
    vector<double> W;
    double W_min = 0.0; //0.504636;
    double W_max = 1.0; //0.946909;
   
    vector<double> D;
    double D_min = 0.0; //0.380281;
    double D_max = 1.0; //0.467912;

    for (int i = 0; i < population_size; i++) {
        genome.clear();
        
        W = create_values(W_min, W_max, numW);
        D = create_values(D_min, D_max, numD);
        
        for (int j = 0; j < numD; j++) {
            gene1.gene_id = 0;
            gene1.gene = create_64b_gene(W.at(j), W_min, W_max);
            gene1.min_value = W_min;     
            gene1.max_value = W_max;
            genome.push_back(gene1);
            
            gene2.gene_id = 1;
            gene2.gene = create_64b_gene(D.at(j), D_min, D_max);
            gene2.min_value = D_min;       
            gene2.max_value = D_max;
            genome.push_back(gene2); 

        }
        
        pop.id_num = i;
        pop.gen_num = 0;
        pop.genome = genome;
        pop.fitness = 0.0;        
        
        population.push_back(pop);
    }    

    return population;

}
vector<_pop> ChildStuntedness::compete_population(vector<_pop> population, int cat) {
// determine the fitness of each member of the population

    _gn gn;
    _pop pop;
    vector<_pop> new_population;
    double w, d;
    int num_gene_groups = population.at(0).genome.size() / 2;

    for (int i = 0; i < population.size(); i++) {
        pop = population.at(i);
        pop.fitness = 0.0;

        w = create_64b_value(
                pop.genome.at(0).gene,
                pop.genome.at(0).min_value,
                pop.genome.at(0).max_value
            );
        d = create_64b_value(
                pop.genome.at(1).gene,
                pop.genome.at(1).min_value,
                pop.genome.at(1).max_value
            );

        pop.fitness = pop.fitness + calculateFitness(cat, w, d);
        new_population.push_back(pop);
    }    
    
    return new_population;

}
vector<_pop> ChildStuntedness::select_population(vector<_pop> population) {
// select the primary candidates for evolution
    
    _pop pop;
    vector<GA_tuple> test_tuple;    
    int popID;
    double fitness;
    vector<_pop> new_population;
    vector<int> ranked_popID;
    vector<double> ranked_fitness;
    
    // create a tuple with <uniqueID, test_score>
    for (int i = 0; i < population.size(); i++) {
        popID = population.at(i).id_num;
        fitness = population.at(i).fitness;
        test_tuple.push_back(make_tuple(popID, fitness));        
    }
    
    // sort and output the tuple
    sort(test_tuple.begin(),test_tuple.end(),sort_GA); 
    
    // transfer the ranked population information
    for(vector<GA_tuple>::iterator iter = test_tuple.begin(); iter != test_tuple.end(); iter++){
        ranked_popID.push_back(get<0>(*iter)); 
        ranked_fitness.push_back(get<1>(*iter));
    }     
    
    for (int j = 0; j < ranked_popID.size(); j++) {
        new_population.push_back(population.at(ranked_popID.at(j)));
    }
    
    return new_population;
}
vector<_pop> ChildStuntedness::evolve_population(vector<_pop> population, int mutate) {
// evolve the population using the top 10 primary candidates
    _pop pop, pop_1, pop_2;
    _gn2 gene_pair, new_gene_pair;
    
    _gn genome_1;
    vector<_pop> new_population;
    default_random_engine rand_engine;

    int pop_size = population.size();
    int top_ten = (int) pop_size/10;
    
    uniform_real_distribution<double> distribution(0, (int) pop_size/4);   

    // use the top 10% of the population
    for (int j = 0; j < top_ten; j++) {
        pop = population.at(j);
        new_population.push_back(pop);
    }

    int i = top_ten;

    // generate 3x mutant clones of the top 10% 
    for (int j = 0; j < 3*top_ten; j++) {
        int rand_pop1 = (int) (distribution(rand_engine));
        pop = population.at(rand_pop1);        
        pop.genome = mutate_genome(pop.genome, mutate);
        pop.id_num = i;
        pop.gen_num = pop.gen_num  + 1;
        pop.fitness = 0.0;

        new_population.push_back(pop);        
        i++;
    }
    
    // generate the offspring using offspring of the top 25%
    while (i < pop_size) {
 
        int rand_pop1 = (int) (distribution(rand_engine));
        int rand_pop2 = (int) (distribution(rand_engine));
   
        pop_1 = population.at(rand_pop1);
        pop_2 = population.at(rand_pop2);
        
        int gen = (int) max(pop_1.gen_num, pop_2.gen_num);

        gene_pair.genome_1 = mutate_genome(pop_1.genome, mutate);
        gene_pair.genome_2 = mutate_genome(pop_2.genome, mutate);
        new_gene_pair = cross_genomes(gene_pair);
        if (i < pop_size) {
            pop.genome = new_gene_pair.genome_1;  // take the first mutated clone
            pop.id_num = i;
            pop.gen_num = gen + 1;
            pop.fitness = 0.0;

            new_population.push_back(pop);
        }

        i++;
  
        if (i < pop_size) {
            pop.genome = new_gene_pair.genome_2;  // take the second mutated clone
            pop.id_num = i;
            pop.gen_num = gen + 1;
            pop.fitness = 0.0;

            new_population.push_back(pop);
        }
        
        i++;        
        
    }

    return new_population;

}
vector<double> ChildStuntedness::create_values(double min, double max, int values) {
    vector<double> random_values;
    
    default_random_engine rand_engine;
    uniform_real_distribution<double> distribution(min, max);    
 
    for (int i = 0; i < values; i++) {
        random_values.push_back(distribution(rand_engine));
    }
    
    return random_values;
}
bitset<64> ChildStuntedness::create_64b_gene(double v, double min, double max) {    
    bitset<64> b_64;
    bitset<64> b_64_max;
    double b_max;    
    double pct_v;
    uint64_t i_v;    
    
    b_64_max.set();
    b_max = (double) b_64_max.to_ullong();
    pct_v = ((v-min)/(max-min)) * b_max;    
    i_v = (uint64_t) pct_v;
    b_64 = bitset<64> (i_v);
    return b_64;
}
double ChildStuntedness::create_64b_value(bitset<64> g, double min, double max) {
    double value;
    double v, v_max;
    bitset<64> b_64_max;
    
    b_64_max.set();
    v_max = (double) b_64_max.to_ullong();
    v = (double) g.to_ullong();
    value = ((v/v_max) * (max-min)) + min;
    return value;
}
vector<_gn> ChildStuntedness::mutate_genome(vector<_gn> old_genome, int mutate) {
    vector<_gn> new_genome;
    _gn gn;
    bitset<64> gene;
    bitset<1> bit;

    for (int i = 0; i < old_genome.size(); i++) {
        gn = old_genome.at(i);
        new_genome.push_back(gn);
    }
    
    int total_bits = old_genome.size() * 64;
 
    for (int i = 0; i < mutate; i++) { 
        int rand_position = rand() % total_bits;   
        int gene_num = (int) rand_position / 64;
        int position = (int) rand_position - (gene_num * 64);

        gn = new_genome.at(gene_num);
        gene = gn.gene;
        bit[0] = gene[position];
        bit.flip();
        gene[position] = bit[0];    
        gn.gene = gene;

        new_genome.at(gene_num) = gn;
    }

    return new_genome;
}
_gn2 ChildStuntedness::cross_genomes(_gn2 &old_genome_pair) {

    _gn2 new_genome_pair;
    
    vector<_gn> genome_1;
    vector<_gn> genome_2;
    
    bitset<1> bit_1;
    bitset<1> bit_2;
    
    bitset<64> gene_1;
    bitset<64> gene_2;

    genome_1 = old_genome_pair.genome_1;
    genome_2 = old_genome_pair.genome_2;

    int total_bits = genome_1.size() * 64;
    int rand_position = rand() % total_bits;
    int total_genes = (int) total_bits / 64;
    int cross_gene = (int) rand_position / 64;
    int cross_pos = rand_position - (cross_gene * 64);

    // exchange the bits one-by-one in the crossover gene
    gene_1 = genome_1.at(cross_gene).gene;
    gene_2 = genome_2.at(cross_gene).gene;
    
    for (int i = cross_pos; i < 64; i++) {
            bit_1[0] = gene_1[i];
            bit_2[0] = gene_2[i];
            gene_2[i] = bit_1[0];
            gene_1[i] = bit_2[0];
    }
    
    genome_1.at(cross_gene).gene = gene_1;
    genome_2.at(cross_gene).gene = gene_2;
    
    // exchange the genes past the crossover gene
    if (cross_gene < total_genes) {
        for (int i = cross_gene+1; i < total_genes; i++) {
            gene_1 = genome_1.at(i).gene;
            gene_2 = genome_2.at(i).gene;
            
            genome_1.at(i).gene = gene_2;
            genome_2.at(i).gene = gene_1;
        }
    }

    new_genome_pair.genome_1 = genome_1;
    new_genome_pair.genome_2 = genome_2;
    
    return new_genome_pair;
}
void ChildStuntedness::output_genome(vector<_gn> genome) {
    // output the genome with the largest bit to the left, smallest to the right
    bitset<64> gene;
    _gn gn;
    
    int size = genome.size()-1;
    
    for (int i = size; i >= 0; i--) {
        gn = genome.at(i);
        gene = gn.gene;
        cout << gene;
    }
}
vector<_result> ChildStuntedness::lookupGAOptimizationTable(void) {
    _result r;
    vector<_result> results;
    
    Analysis A;
    double NN_cycle_start = clock();
    double min1, mean1, max1, min2, mean2, max2;
    vector<double> weights, durations;
    _train tc;

    for (int j = 0; j < xM.size(); j++) { // need to generate a similar structure to _train for _test     
        for (int i = 0; i < tC.size(); i++) {
            tc = tC.at(i);
            if (xM.at(j).category == tC.at(i).category) {
                r.ID = tC.at(i).category;
                r.weight = tC.at(i).weight;
                r.duration = tC.at(i).duration;
                results.push_back(r);
                break;
            }
        }
    }
    
    for (int i = 0; i < results.size(); i++) {
        weights.push_back(results.at(i).weight);
        durations.push_back(results.at(i).duration);
    }
    
    min1 = A.min(weights);
    mean1 = A.average(weights);
    max1 = A.max(weights);
    min2 = A.min(durations);
    mean2 = A.average(durations);
    max2 = A.max(durations);
    
    double NN_total_time = (clock() - NN_cycle_start) / 1000000.0;
    
    cout << "--------------------------------------------------------------------------------" << endl;    
    cout << "Testing Parameter Statistics:" << endl;
    cout << "Weight:  \tmin1: " << min1 << "\tmean1: " << mean1 << "\tmax1: " << max1 << endl;
    cout << "Duration:\tmin2: " << min2 << "\tmean2: " << mean2 << "\tmax2: " << max2 << endl;
    cout << "Total Number of Results: " << results.size() << endl;
    cout << "Total Time for Testing:  " << NN_total_time << " sec" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;  
    
    return results;
    
}

/*
 * Generate Optimization #2 - optimization of a multivariate linear interpolation equation <GA2>
 */
void ChildStuntedness::generate_optimization_2(_cases C, int mutate, int popSize, int maxCycle) {
    _pop pop;
    
    double w_min = -1.0e0;
    double w_max = 1.0e0;
    double d_min = -1.0e0;
    double d_max = 1.0e0;
    
    cout << "weight TypeI: ";
    wN_typeI = cycle_population(C.typeI, C.T1, mutate, popSize, maxCycle, w_min, w_max, true);
    cout << "weight TypeII: ";
    wN_typeII = cycle_population(C.typeII, C.T2, mutate, popSize, maxCycle, w_min, w_max, true);
    cout << "weight TypeIII: ";
    wN_typeIII = cycle_population(C.typeIII, C.T3, mutate, popSize, maxCycle, w_min, w_max, true);
    cout << "weight TypeIV: ";
    wN_typeIV = cycle_population(C.typeIV, C.T4, mutate, popSize, maxCycle, w_min, w_max, true);

    cout << "duration TypeI: ";
    dN_typeI = cycle_population(C.typeI, C.T1, mutate, popSize, maxCycle, d_min, d_max, false);
    cout << "duration TypeII: "; 
    dN_typeII = cycle_population(C.typeII, C.T2, mutate, popSize, maxCycle, d_min, d_max, false);
    cout << "duration TypeIII: "; 
    dN_typeIII = cycle_population(C.typeIII, C.T3, mutate, popSize, maxCycle, d_min, d_max, false);
    cout << "duration TypeIV: "; 
    dN_typeIV = cycle_population(C.typeIV, C.T4, mutate, popSize, maxCycle, d_min, d_max, false);
}
vector<double> ChildStuntedness::cycle_population(vector<_case> type, vector<_box> Tx, int mutate, int popSize, int maxCycle, double min, double max, bool isW) {
    vector<_pop> pop, pop_c, pop_r;
    vector<_gn> genome;
    double x, y;
    _gn g;
    vector<double> N;
    double n;
    vector<double> Y;
    Analysis A;
    double targetVal;
    vector<_box> X;
    int m;

    pop = create_population_2(popSize, min, max);
     
    vector<int> IDs, uniqueIDs;
    vector<_train> tX;
    
    tX = buildCategorySet(type, true); // type -> {ID, cat, weight, duration}

    double score = -1e99;
    int gen = 0;
    int j = 0;

    while (score < 5.0) {
        gen++;
        pop_c.clear();
        pop_r.clear();
        
        m = round((double) mutate 
                * (1.0 - cos(3.141592654 
                * ((double) gen
                / (double) popSize))) 
                + 1.0);
        //cout << "m: " << m << "\tgen: " << gen << endl;
        pop_c = compete_population_2(pop, tX, isW);
        pop_r = select_population(pop_c);
        pop = evolve_population(pop_r, m);
        j++;

        double newScore = pop.at(0).fitness;

        if (score < newScore) {
            score = newScore;
            //cout << "Generation #: " << gen << " Fitness: " << score << endl;
            j=0;
        }
    }
    
    genome = pop.at(0).genome;

    // build the vector of coefficients N
    for (int i = 0; i < genome.size(); i++) {
        g = genome.at(i);
        n = create_64b_value(g.gene, g.min_value, g.max_value);
        N.push_back(n);
    }   
    
    Y = calculateParameter(N, Tx);
    
    cout << "\tMutate: " << m 
        << "\tGen: " << gen 
        << "\tBest Fitness: " << score
        << endl;

    outputCoefficients(genome, isW);
    cout << "----------------------------------------------------------------------------------" << endl;   
    
    return N;
}
void ChildStuntedness::outputCoefficients(vector<_gn> genome, bool isW) {
    // build the vector of coefficients N
    
    _gn g;
    double n;
    vector<double> N;
    for (int i = 0; i < genome.size(); i++) {
        g = genome.at(i);
        n = create_64b_value(g.gene, g.min_value, g.max_value);
        N.push_back(n);
    } 
    if (isW) {
        cout << "w = "; 
    } else {
        cout << "d = ";
    }
    cout << N.at(0)
            << " + (" << N.at(1) << " * age)"
            << " + (" << N.at(2) << " * ODV5)"
            << " + (" << N.at(3) << " * ODV6)"
            << " + (" << N.at(4) << " * ODV7)"
            << " + (" << N.at(5) << " * ODV8)"
            << " + (" << N.at(6) << " * ODV9)"
            << " + (" << N.at(7) << " * ODV10)"
            << " + (" << N.at(8) << " * ODV11)"   
            << " + (" << N.at(9) << " * ODV12)"
            << endl;
    
}
vector<_box> ChildStuntedness::convertCaseToSmallBox(vector<_box> Tx) {
    vector<_box> X;
    
    X.push_back(Tx.at(0)); // age
    X.push_back(Tx.at(3)); // ODV5
    X.push_back(Tx.at(4)); // ODV6
    X.push_back(Tx.at(5)); // ODV7
    X.push_back(Tx.at(6)); // ODV8
    X.push_back(Tx.at(7)); // ODV9
    X.push_back(Tx.at(8)); // ODV10
    X.push_back(Tx.at(9)); // ODV11
    X.push_back(Tx.at(10)); // ODV12 
    return X;
}
vector<double> ChildStuntedness::calculateParameter(vector<double> N, vector<_box> Tx) {
    vector<double> Y;
    double x, y;
    
    vector<_box> X = convertCaseToSmallBox(Tx);

    int xVal = X.at(0).values.size(); // size of the age vector should be explicit to define of data points to calculate
    //cout << "size N: " << N.size() << "size X: " << X.size() << endl;
    
    // calculate the value using multivariable linear interpolation
    for (int i = 0; i < xVal; i++) { // loop thru each value in the box
        y = N.at(0);
        for (int j = 0; j < X.size(); j++) { // loop thru each box in the tessaract
        //cout << j << "\ty:" << y << endl;
            x = X.at(j).values.at(i); // value in the box
            y += N.at(j+1) * x; // N is offset by one to account for the N0 parameter which is not multiplied by X.
        }

        Y.push_back(y);
    }
    return Y;
}
vector<_pop> ChildStuntedness::create_population_2(int population_size, double min, double max) {
// generate a population for weight coefficients 
// using with parameters defined to the objective function / competition
    
    // for each case:
    //     1) weight from 0-1
    //     2) duration from 0-1

    _gn gene;
    vector<_gn> genome;
    _pop pop;
    vector<_pop> population;  
    _coeff xN;
    vector<_coeff> xNs;
    
    int numGenes = 1; // 4 genes for each type: type I (male/nut=1), type II (male/nut=2), type III (female/nut=1), type IV (female/nut=2)
    int geneID = 0;
    
    for (int i = 0; i < population_size; i++) {
        genome.clear();
        xNs.clear();
        
        for (int j = 0; j < 10; j++) { // N0 (zero intercept), N1 (age), N2 (ODV5), N3 (ODV6), N4 (ODV7), N5 (ODV8), N6 (ODV9), N7 (ODV10), N8 (ODV11), N9 (ODV12)
            xN.N_min = min;
            xN.N_max = max;
            xN.N = create_values(min, max, numGenes); 
            xNs.push_back(xN); // coefficients for weight vs. 4 genes
        }
        
        // load the random coefficients for weight into the genome as a 64 bit gene
        for (int j = 0; j < xNs.size(); j++) {
            for (int k = 0; k < numGenes; k++) {
                gene.gene_id = geneID;
                gene.gene = create_64b_gene(xNs.at(j).N.at(k), xNs.at(j).N_min, xNs.at(j).N_max);
                gene.min_value = xNs.at(j).N_min;
                gene.max_value = xNs.at(j).N_max;
                geneID++;
                genome.push_back(gene);
            }
        }
        
        pop.id_num = i;
        pop.gen_num = 0;
        pop.genome = genome;
        pop.fitness = 0.0;        
        
        population.push_back(pop);
    }    

    return population;
}
vector<_pop> ChildStuntedness::compete_population_2(vector<_pop> population, vector<_train> tX, bool isW) {
// determine the fitness of each member of the population

    _gn gn;
    _pop pop;
    vector<_pop> new_population;
    double w, d, targetVal, err;
    int num_gene_groups = population.at(0).genome.size() / 2;
    
    int size_tX = tX.size(); 
    
    for (int i = 0; i < population.size(); i++) {
        pop = population.at(i);
        err = 0.0;
        if (isW) {
            targetVal = tX.at(i).weight;
        } else {
            targetVal = tX.at(i).duration;
        }
        
        for (int j = 0; j < size_tX; j++) {
            if (isW) {
                w = create_64b_value(
                    pop.genome.at(0).gene,
                    pop.genome.at(0).min_value,
                    pop.genome.at(0).max_value
                    );
                err += optimizationError(targetVal,w);

            } else {
                d = create_64b_value(
                    pop.genome.at(0).gene,
                    pop.genome.at(0).min_value,
                    pop.genome.at(0).max_value
                    );
                err += optimizationError(targetVal,d); 

            }
        }

        pop.fitness = -log10(err/size_tX); // 
        new_population.push_back(pop);
    }  
    
    return new_population;

}
/*
 * <GA2>
 */

_odv ChildStuntedness::doesStringContainNA(string str) {
    bool containsNA = false;
    double value = 0.0;
    _odv odvData;
    size_t found_NA;

    found_NA = str.find("NA");
    
    if (found_NA!=std::string::npos) {
        containsNA = true;
        value = 0.0;
    } else {
        containsNA = false;
        istringstream ss(str);
        ss >> value;
        if (value == 0.0) {
            containsNA = true;
        }
    }
    
    odvData.isNA = containsNA;
    odvData.odv_value = value;
    odvData.str = str;
    
    return odvData;
}
vector<double> ChildStuntedness::vectorWithoutZeros(vector<double> v) {
    vector<double> noZeroVector;
    for (int i = 0; i < v.size(); i++) {
        if (v.at(i) != 0.0) {
            noZeroVector.push_back(v.at(i));
        }
    }
    return noZeroVector;
}
vector<double> makeParameterList(vector<_box> L, int index) {
    
    vector<double> x;
    
    x.push_back(L.at(0).values.at(index)); // age
    x.push_back(L.at(3).values.at(index)); // odv5
    x.push_back(L.at(4).values.at(index)); // odv6
    x.push_back(L.at(5).values.at(index)); // odv7
    x.push_back(L.at(6).values.at(index)); // odv8
    x.push_back(L.at(7).values.at(index)); // odv9
    x.push_back(L.at(8).values.at(index)); // odv10
    x.push_back(L.at(9).values.at(index)); // odv11
    x.push_back(L.at(10).values.at(index)); // odv12

    return x;
}
/* <genie>
 * Optimization of multivariate linear interpolation equation:
 √ 1. configure the data for four sets of data including:  1) men/nutn=1, 2) men/nutn=2, 3) women/nutn=1, 4) women/nutn=2
 √ 2. initialize coefficient sets for weight (wN) and duration (dN) with random values
 √ 3. utilize the create population method to fill the wN and dN coefficients with random initial values.
 * 4. generate score to determine the fitness of interpolation function
 * 5. modify the coefficients to reduce error of estimated weight (eW) and estimated due date (eD) for four sets of of data
 * 6. output the estimated weight and duration values.
 */
_cases ChildStuntedness::splitCaseData(vector<_case> Cs) {
    _cases S;
    
    bool outputStats = false;
    S.typeI = vectorType(Cs, 0, 1);
    S.typeII = vectorType(Cs, 0, 2);
    S.typeIII = vectorType(Cs, 1, 1);
    S.typeIV = vectorType(Cs, 1, 2);
    
    loadParameters(S.typeI); // Cs -> S.typeI -> {age,sex,nutrition,ODV5-12,wt,duration}
    S.T1 = buildCases(); // {age,sex,nutrition,ODV5-12,wt,duration} -> S.T1
    clearParameters(); // clear -> {age,sex,nutrition,ODV5-12,wt,duration}
    
    loadParameters(S.typeII); // Cs -> S.typeII -> {age,sex,nutrition,ODV5-12,wt,duration}
    S.T2 = buildCases(); // {age,sex,nutrition,ODV5-12,wt,duration} -> S.T2
    clearParameters(); // clear -> {age,sex,nutrition,ODV5-12,wt,duration}
    
    loadParameters(S.typeIII); // Cs -> S.typeIII -> {age,sex,nutrition,ODV5-12,wt,duration}
    S.T3 = buildCases(); // {age,sex,nutrition,ODV5-12,wt,duration} -> S.T3
    clearParameters(); // clear -> {age,sex,nutrition,ODV5-12,wt,duration}
    
    loadParameters(S.typeIV); // Cs -> S.typeIII -> {age,sex,nutrition,ODV5-12,wt,duration}
    S.T4 = buildCases(); // {age,sex,nutrition,ODV5-12,wt,duration} -> S.T4
    clearParameters(); // clear -> {age,sex,nutrition,ODV5-12,wt,duration}
    
    if (outputStats) {
        cout << "--------------------------------------------------------------------------------" << endl;   
        cout << "Mean Values from Training Data: " << endl;
        cout << "Type 1: " << S.typeI.size() << "\t";
        for (int i = 0; i < S.T1.size(); i++) {
            cout << S.T1.at(i).mean << " ";
        }
        cout << endl;
        cout << "Type 2: " << S.typeII.size() << "\t";
        for (int i = 0; i < S.T2.size(); i++) {
            cout << S.T2.at(i).mean << " ";
        }
        cout << endl;
        cout << "Type 3: " << S.typeIII.size() << "\t";
        for (int i = 0; i < S.T3.size(); i++) {
            cout << S.T3.at(i).mean << " ";
        }
        cout << endl;
        cout << "Type 4: " << S.typeIV.size() << "\t";
        for (int i = 0; i < S.T4.size(); i++) {
            cout << S.T4.at(i).mean << " ";
        }
        cout << endl;

        cout << "--------------------------------------------------------------------------------" << endl;   
        cout << "Min Values from Training Data: " << endl;
        cout << "Type 1: " << S.typeI.size() << "\t";
        for (int i = 0; i < S.T1.size(); i++) {
            cout << S.T1.at(i).min << " ";

        }
        cout << endl;
        cout << "Type 2: " << S.typeII.size() << "\t";
        for (int i = 0; i < S.T2.size(); i++) {
            cout << S.T2.at(i).min << " ";

        }
        cout << endl;
        cout << "Type 3: " << S.typeIII.size() << "\t";
        for (int i = 0; i < S.T3.size(); i++) {
            cout << S.T3.at(i).min << " ";

        }
        cout << endl;
        cout << "Type 4: " << S.typeIV.size() << "\t";
        for (int i = 0; i < S.T4.size(); i++) {
            cout << S.T4.at(i).min << " ";

        }
        cout << endl;

        cout << "--------------------------------------------------------------------------------" << endl;   
        cout << "Max Values from Training Data: " << endl;
        cout << "Type 1: " << S.typeI.size() << "\t";
        for (int i = 0; i < S.T1.size(); i++) {
            cout << S.T1.at(i).max << " ";

        }
        cout << endl;
        cout << "Type 2: " << S.typeII.size() << "\t";
        for (int i = 0; i < S.T2.size(); i++) {
            cout << S.T2.at(i).max << " ";

        }
        cout << endl;
        cout << "Type 3: " << S.typeIII.size() << "\t";
        for (int i = 0; i < S.T3.size(); i++) {
            cout << S.T3.at(i).max << " ";

        }
        cout << endl;
        cout << "Type 4: " << S.typeIV.size() << "\t";
        for (int i = 0; i < S.T4.size(); i++) {
            cout << S.T4.at(i).max << " ";

        }
        cout << endl;  
    }
    
    return S;
}
vector<_case> ChildStuntedness::vectorType(vector<_case> xCs, double sex, double nutrition) {
    vector<_case> yCs;
    _case c;
    
    for (int i = 0; i < xCs.size(); i++) {
        c = xCs.at(i);
        if ((c.sex == sex) && (c.nutrition == nutrition)) {
            yCs.push_back(xCs.at(i));
        }
    }
    return yCs;
}

/*
 * end of <genie>
 */

void ChildStuntedness::parseData(vector<string> &dataString, bool isTraining) {
            
    for (int i = 0; i < dataString.size(); i++) {
        odv_case.clear();
        readLine = dataString.at(i);
        replace(readLine.begin(), readLine.end(), ',', ' ');
        istringstream ss(readLine);

        if (isTraining) {
            // cout << "Training readline: " << i << " >>> " << readLine << endl;
            if (!(ss >> hC.ID >> hC.age >> hC.sex >> hC.nutrition >> str_5 >> str_6 
                    >> str_7 >> str_8 >> str_9 >> str_10 >> str_11 >> str_12 
                    >> hC.weight >> hC.duration)) {
                break;
            };
        } else {
            // cout << "Testing readline: " << i << " >>> " << readLine << endl;
            if (!(ss >> nC.ID >> nC.age >> nC.sex >> nC.nutrition >> str_5 >> str_6 
                    >> str_7 >> str_8 >> str_9 >> str_10 >> str_11 >> str_12)) {
                break;
            };
            nC.weight = 0.0; // to be determined;
            nC.duration = 0.0; // to be determined;
        }
        ss.str("");
        
        odv_case.push_back(doesStringContainNA(str_5));
        odv_case.push_back(doesStringContainNA(str_6));        
        odv_case.push_back(doesStringContainNA(str_7));        
        odv_case.push_back(doesStringContainNA(str_8));
        odv_case.push_back(doesStringContainNA(str_9));
        odv_case.push_back(doesStringContainNA(str_10));
        odv_case.push_back(doesStringContainNA(str_11));
        odv_case.push_back(doesStringContainNA(str_12));       
          
        if (isTraining) {
            hC.odv = odv_case;
            hCs.push_back(hC);
        } else {
            nC.odv = odv_case;
            nCs.push_back(nC);            
        }
    }
}
_box ChildStuntedness::loadBox(string name, vector<double> element, bool allowZero) {
    Analysis A;
    _box B;
    vector<double> nonZeroElement;
    
    if (!allowZero) {
        nonZeroElement = vectorWithoutZeros(element);
        B.name = name;
        B.min = A.min(nonZeroElement);
        B.mean = A.average(nonZeroElement);
        B.max = A.max(nonZeroElement);
        B.values = element;
        B.values_no_zero = nonZeroElement;
        B.category = A.category(nonZeroElement, element);
    } else {
        B.name = name;
        B.min = A.min(element);
        B.mean = A.average(element);
        B.max = A.max(element);
        B.values = element;
        B.values_no_zero = element;
        B.category = A.category(element, element);
    }
     
    return B;
}
int ChildStuntedness::generateCategoryCode(vector<int> codes) {
    int code = codes.at(0)+2;
    
    for (int i = 1; i < codes.size(); i++) {
        code += (codes.at(i)+2)*pow(5.0,i);
    }
    
    return code;
}
vector<_box> ChildStuntedness::quantifyParameters(vector<_case> Cs) {
    vector<_box> M;    
    loadParameters(Cs); // Cs -> {age,sex,nutrition,ODV5-12,wt,duration}
    M = buildCases(); // {age,sex,nutrition,ODV5-12,wt,duration} -> M
    return M;
}
void ChildStuntedness::loadParameters(vector<_case> type) {
    // type 
    for (int i = 0; i < type.size(); i++) {
        age.push_back(hCs.at(i).age);
        sex.push_back((double) hCs.at(i).sex);
        nutrition.push_back((double) hCs.at(i).nutrition);
        ODV5.push_back(hCs.at(i).odv.at(0).odv_value);
        ODV6.push_back(hCs.at(i).odv.at(1).odv_value);
        ODV7.push_back(hCs.at(i).odv.at(2).odv_value);
        ODV8.push_back(hCs.at(i).odv.at(3).odv_value);
        ODV9.push_back(hCs.at(i).odv.at(4).odv_value);
        ODV10.push_back(hCs.at(i).odv.at(5).odv_value);
        ODV11.push_back(hCs.at(i).odv.at(6).odv_value);
        ODV12.push_back(hCs.at(i).odv.at(7).odv_value);
        weight.push_back(hCs.at(i).weight);
        duration.push_back(hCs.at(i).duration);
    }
}
void ChildStuntedness::clearParameters(void) {
    age.clear();
    sex.clear();
    nutrition.clear();
    ODV5.clear();
    ODV6.clear();
    ODV7.clear();
    ODV8.clear();
    ODV9.clear();
    ODV10.clear();
    ODV11.clear();
    ODV12.clear();
    weight.clear();
    duration.clear();
}
vector<_box> ChildStuntedness::buildCases(void) {
    vector<_box> Bx;
    Bx.push_back(loadBox("age------", age, false));
    Bx.push_back(loadBox("sex------", sex, true));
    Bx.push_back(loadBox("nutrition", nutrition, true));
    Bx.push_back(loadBox("ODV5-----", ODV5, false)); // 0
    Bx.push_back(loadBox("ODV6-----", ODV6, false)); // 1
    Bx.push_back(loadBox("ODV7-----", ODV7, false)); // 2
    Bx.push_back(loadBox("ODV8-----", ODV8, false)); // 3
    Bx.push_back(loadBox("ODV9-----", ODV9, false)); // 4
    Bx.push_back(loadBox("ODV10----", ODV10, false)); // 5
    Bx.push_back(loadBox("ODV11----", ODV11, false)); // 6 
    Bx.push_back(loadBox("ODV12----", ODV12, false)); // 7
    Bx.push_back(loadBox("weight---", weight, false)); 
    Bx.push_back(loadBox("duration-", duration, false));
    
    return Bx;
}
vector<int> ChildStuntedness::buildCategories(vector<_box> M) {
    vector<int> codes;
    vector<int> c_age, c_sex, c_nutn, c_ODV5, c_ODV6, c_ODV7, c_ODV8, c_ODV9, c_ODV10, c_ODV11, c_ODV12;
    int category;
    vector<int> categories;

    c_age = M.at(0).category;
    c_sex = M.at(1).category;
    c_nutn = M.at(2).category;
    c_ODV5 = M.at(3).category;
    c_ODV6 = M.at(4).category;
    c_ODV7 = M.at(5).category; 
    c_ODV8 = M.at(6).category;
    c_ODV9 = M.at(7).category;
    c_ODV10 = M.at(8).category;
    c_ODV11 = M.at(9).category;
    c_ODV12 = M.at(10).category;

    for (int i = 0; i < c_age.size(); i++) {
        codes.clear();
        codes.push_back(c_age.at(i));
        codes.push_back(c_sex.at(i)); 
        codes.push_back(c_nutn.at(i)); 
        codes.push_back(c_ODV5.at(i));
        codes.push_back(c_ODV6.at(i)); 
        codes.push_back(c_ODV7.at(i));
        codes.push_back(c_ODV8.at(i));
        codes.push_back(c_ODV9.at(i));
        codes.push_back(c_ODV10.at(i));
        codes.push_back(c_ODV11.at(i));
        codes.push_back(c_ODV12.at(i));
        category = generateCategoryCode(codes);
        categories.push_back(category);
    }
    
    return categories;
}
vector<int> ChildStuntedness::compareCategories(vector<int> trainingCat, vector<int> testingCat) {
    vector<int> uniqueCat;
    bool isInTraining;
    bool isInUniqueCat;
    
    for (int i = 0; i < testingCat.size(); i++) {
        isInTraining = false;
        for (int j = 0; j < trainingCat.size(); j++) {
            if (testingCat.at(i) == trainingCat.at(j)) {
                isInTraining = true;
            }

        }
        if (!isInTraining) {
            for (int k = 0; k < uniqueCat.size(); k++) {
                isInUniqueCat = false;
                if (uniqueCat.at(k) == testingCat.at(i)) {
                    isInTraining = true;
                }
                if (!isInUniqueCat) {
                    uniqueCat.push_back(testingCat.at(i));
                }
            }
        }
    }

    return uniqueCat;
}
vector<_train> ChildStuntedness::buildCategorySet(vector<_case> Cs, bool isTraining) {
    _train t;
    vector<_train> M;
    _case c;
    
    for (int i = 0; i < Cs.size(); i++) {
        c = Cs.at(i);
        t.ID = c.ID;
        if (isTraining) {
            t.category = trainingCategories.at(i);
        } else {
            t.category = testingCategories.at(i);
        }
        t.weight = c.weight;
        t.duration = c.duration;
        M.push_back(t);
    }

    return M;
}
vector<_train> ChildStuntedness::generateOptimizationTable(void) {
    _train t;
    double weight, duration;
    int dataPts;
    
    vector<_train> trainSet;

    for (int i = 0; i < unique_trainingCategories.size(); i++) {
        // go thru the entire list of data points to determine the average weight/duration for each category
        int cat = unique_trainingCategories.at(i);
        weight = 0.0;
        duration = 0.0;
        dataPts = 0.0;
        
        for (int j = 0; j < hCs.size(); j++) {
            hC = hCs.at(j);
            if (trainingCategories.at(j) == cat) {
                weight += hC.weight;
                duration += hC.duration;
                dataPts++;
            }
        }
        t.ID = i;
        t.category = cat;
        t.weight = weight/dataPts;
        t.duration = duration/(double) dataPts;
        trainSet.push_back(t);
    }
    return trainSet;
}
vector<_result> ChildStuntedness::lookupCoefficientTable(void) {
    _result r;
    vector<_result> rS;
    double NN_cycle_start = clock();
    double dur, wt;
    _box b;
    vector<_box> x;
    
    vector<double> d_typeI, d_typeII, d_typeIII, d_typeIV;
    vector<double> w_typeI, w_typeII, w_typeIII, w_typeIV;
    
    dN_typeI = calculateParameter(dN_typeI, xS.T1);
    dN_typeII = calculateParameter(dN_typeII, xS.T2);
    dN_typeIII = calculateParameter(dN_typeIII, xS.T3);
    dN_typeIV = calculateParameter(dN_typeIV, xS.T4);
   
    wN_typeI = calculateParameter(wN_typeI, xS.T1);
    wN_typeII = calculateParameter(wN_typeII, xS.T2);
    wN_typeIII = calculateParameter(wN_typeIII, xS.T3);
    wN_typeIV = calculateParameter(wN_typeIV, xS.T4);
    
    for (int i = 0; i < xS.typeI.size(); i++) {
        r.ID = xS.typeI.at(i).ID;
        r.duration = dN_typeI.at(i); // need to confirm that the length from T1 is equal to typeI
        r.weight = wN_typeI.at(i);
        rS.push_back(r);
    }
    for (int i = 0; i < xS.typeII.size(); i++) {
        r.ID = xS.typeII.at(i).ID;
        r.duration = dN_typeII.at(i); // need to confirm that the length from T1 is equal to typeI
        r.weight = wN_typeII.at(i);
        rS.push_back(r);
    }
    for (int i = 0; i < xS.typeIII.size(); i++) {
        r.ID = xS.typeIII.at(i).ID;
        r.duration = dN_typeIII.at(i); // need to confirm that the length from T1 is equal to typeI
        r.weight = wN_typeIII.at(i);
        rS.push_back(r);
    }
    for (int i = 0; i < xS.typeIV.size(); i++) {
        r.ID = xS.typeIV.at(i).ID;
        r.duration = dN_typeIV.at(i); // need to confirm that the length from T1 is equal to typeI
        r.weight = wN_typeIV.at(i);
        rS.push_back(r);
    }
    
    return rS;
}
vector<_result> ChildStuntedness::lookupOptimizationTable(void) {
    _result r;
    vector<_result> rS;
  
    double NN_cycle_start = clock();

    for (int j = 0; j < xM.size(); j++) { // need to generate a similar structure to _train for _test     
        for (int i = 0; i < tC.size(); i++) {
            if (xM.at(j).category == tC.at(i).category) {
                r.ID = xM.at(j).ID;
                r.weight = tC.at(i).weight;
                r.duration = tC.at(i).duration;
                rS.push_back(r);
                break;
            }
        }
    }

    double NN_total_time = (clock() - NN_cycle_start) / 1000000.0;
        
    cout << "Total Number of Results: " << rS.size() << endl;
    cout << "Total Time for Testing:  " << NN_total_time << " sec" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;  
    
    return rS;
    
}
vector<double> ChildStuntedness::generateResult(vector<_result> resultVector) {

    _result oR;
    vector<double> outputVector;
    bool ID_inList;
    
    vector<int> IDList;  
    vector<_result> outputResult;
    
    vector<object_tuple> index_tuple;    
    vector<int> sortIDs;
    _result sort_result;
    vector<_result> sortedData;
    
    _result result;

    // make a list of IDs
    for (int i = 0; i < resultVector.size(); i++) {
        ID_inList = false;
        for (int j = 0; j < IDList.size(); j++) {
            result = resultVector.at(i);
            if (result.ID == IDList.at(j)) { 
                ID_inList = true;
            }
        }
        if (!ID_inList) {
            IDList.push_back(resultVector.at(i).ID);
        }
    }

    for (int i = 0; i < IDList.size(); i++) {
        double wt = 0;
        double dur = 0;
        int dataPts = 0;
        for (int j = 0; j < resultVector.size(); j++) {
            if (IDList.at(i) == resultVector.at(j).ID) { 
                // this will only take the last registered values in the calculated/testing data
                //oR.ID = resultVector.at(j).ID;
                //oR.weight = resultVector.at(j).weight;
                //oR.duration = resultVector.at(j).duration;
                
                // averaging
                wt += resultVector.at(j).weight;
                dur += resultVector.at(j).duration;
                dataPts++;
            }
        }
        // this will only take the last registered values in the calculated/testing data
        oR.ID = IDList.at(i);
        oR.weight = wt/dataPts;
        oR.duration = dur / dataPts;
        outputResult.push_back(oR);
        
        //cout << "ID: " << oR.ID << "\tw: " << oR.weight << "\td: " << oR.duration << endl;
    }
    
    for (int i = 0; i < outputResult.size(); i++) {
        index_tuple.push_back(make_tuple(i, outputResult.at(i).ID));    
    }
    
    // sort the index with the lowest IDs first
    sort(index_tuple.begin(),index_tuple.end(),sort_object_score);  

    for(vector<object_tuple>::iterator iter = index_tuple.begin(); iter != index_tuple.end(); iter++){
        sortIDs.push_back(get<0>(*iter));
    }  

    for (int i = 0; i < sortIDs.size(); i++) {
        int indx = sortIDs.at(i);
        sort_result.ID = outputResult.at(indx).ID;
        sort_result.weight = outputResult.at(indx).weight;
        sort_result.duration = outputResult.at(indx).duration;
        sortedData.push_back(sort_result);
        
        if (i < 10) {
            cout << "#: " << sort_result.ID
                 << "\td: " << sort_result.duration
                 << "\tw: " << sort_result.weight
                 << endl;
        }
        
        outputVector.push_back(sort_result.duration); // dur
        outputVector.push_back(sort_result.weight); // wt
    }
    
    return outputVector;

}
void ChildStuntedness::calcResultStats(void) {
    Analysis A;
    vector<double> wt, dur;
    double wt_min, wt_max, wt_mean, dur_min, dur_max, dur_mean;
    wt.clear();
    dur.clear(); 
    
    for (int i = 0; i < hCs.size(); i++) {
        wt.push_back(hCs.at(i).weight);
        dur.push_back(hCs.at(i).duration);
    }
    
    wt_min = A.min(wt);
    wt_mean = A.average(wt);
    wt_max = A.max(wt);
    
    dur_min = A.min(dur);
    dur_mean = A.average(dur);
    dur_max = A.max(dur);
    
    cout << "Weight:  \tmin1: " << wt_min << "\tmean1: " << wt_mean << "\tmax1: " << wt_max << endl;
    cout << "Duration:\tmin2: " << dur_min << "\tmean2: " << dur_mean << "\tmax2: " << dur_max << endl;
    cout << "--------------------------------------------------------------------------------" << endl;  
}
double ChildStuntedness::calcScore(vector<double> w, vector<double> w_real, vector<double> d, vector<double> d_real) {
    Analysis A;
    double error = 0;
    double error_0 = 0;
    
    double inverseS_00 = 3554.42; // for the complete set of data
    double inverseS_01 = -328.119; // for the complete set of data
    double inverseS_10 = -328.119; // for the complete set of data
    double inverseS_11 = 133.511; // for the complete set of data 
    
    double meanD = A.average(d_real);
    double meanW = A.average(w_real);

    for (int i = 0; i < w.size(); i++) {
        double deltaW = w.at(i) - w_real.at(i);
        double deltaD = d.at(i) - d_real.at(i);
        double val1 = deltaD * inverseS_00 + deltaW * inverseS_10;
        double val2 = deltaD * inverseS_01 + deltaW * inverseS_11;
        double ei = (val1 * deltaD) + (val2 * deltaW);
        error += ei;
    }
    
    for (int i = 0; i < w.size(); i++) {
        double deltaW = meanW - w_real.at(i);
        double deltaD = meanD - d_real.at(i);
        double val1 = (deltaD * inverseS_00) + (deltaW * inverseS_10);
        double val2 = (deltaD * inverseS_01) + (deltaW * inverseS_11);
        double ei = (val1 * deltaD) + (val2 * deltaW);
        error_0 += ei;
    }

    double score = 1000000.0 * (1 - (error/error_0));
    
    return score;

}
vector<double> ChildStuntedness::predict(vector<string> training, vector<string> testing) {
    /* Column   Variable   Type    Label/Description
      1     Id         int     Unique Fetus ID
      2     t.ultsnd   float   Estimated fetus gestational age from last menstrual recall date
      3     Sex        int     0 = Male, 1 = Female
      4     Status     int     Maternal nutritional status (1 or 2)
    5-12    Odv        float   Dependent variables: Ultrasound observed measurements
     13     Birth Sz   float   Birth Weight (w)
     14     Duration   float   Pregnancy Duration, or Birthday (b)
    */
    
/*
 * 
 */ 

    Analysis A;

    // initialize parameters and vectors
    vector<double> the_big_answer;
    vector<_result> predictedData;
    vector<int> uniqueCategories;
    vector<double> wt, dur;
    double wt_min, wt_max, wt_mean, dur_min, dur_max, dur_mean;

    // parse the training data
    parseData(training, true); // training -> hCs       
    
    /*
    // construct a NN to estimate the ODV values that are missing
    int inputs = 10;
    int hidden = 12;
    int outputs = 2;
    buildNNtopology(inputs,hidden,outputs);
    
    cout << "--------------------------------------------------------------------------------" << endl;   
    cout << "Building a <" << inputs << " " << hidden << " " << outputs << "> Neural Network..." << endl;*/
    
    // perform training on data in hCs
    T = quantifyParameters(hCs); // hCs => T
    //trainNN(119.0); // perform 10sec NN training.  T -> NN {age, sex, nutrition}

    S = A.covariance(T.at(12).values, T.at(11).values); // (d, w))
    S_1 = A.inverse_2x2(S);
    
    trainingCategories = buildCategories(T);
    unique_trainingCategories = A.unique(trainingCategories);
    tS = splitCaseData(hCs);
    
    cout << "--------------------------------------------------------------------------------" << endl;   
    cout << "Training Data Statistics:" << endl;
    cout << "Total Number of Training Data Points:  " << hCs.size() << endl;
    cout << "Training Categories: " << unique_trainingCategories.size() << endl;
    calcResultStats();
    
    // load the testing data
    parseData(testing, false); // testing -> nCs
    X = quantifyParameters(nCs);
    testingCategories = buildCategories(X);
    unique_testingCategories = A.unique(testingCategories);
    uniqueCategories = compareCategories(trainingCategories, testingCategories);
    xS = splitCaseData(nCs);
    
    // build the training data based on categories
    tM = buildCategorySet(hCs, true);
    xM = buildCategorySet(nCs, false);
    //generate_optimization();
    //generate_optimization_2(tS, 10, 250, 5000);

    // estimate the w,d values using GA generated categorizations
    tC = generateOptimizationTable();
    
    //predictedData = lookupCoefficientTable(); // for generate_optimization_2 method
    predictedData = lookupOptimizationTable(); // for category methods
    //predictedData = lookupGAOptimizationTable(); // for generate_optimization method
    //predictedData = execNN(); // execNN NN {age, sex, nutrition} -> ODV5-12 values that are N/A
    the_big_answer = generateResult(predictedData);


    wt.clear();
    dur.clear();
    
    for (int i = 0; i < the_big_answer.size(); i=i+2) {
        dur.push_back(the_big_answer.at(i));
        wt.push_back(the_big_answer.at(i+1));
    }
    
    wt_min = A.min(wt);
    wt_mean = A.average(wt);
    wt_max = A.max(wt);
    
    dur_min = A.min(dur);
    dur_mean = A.average(dur);
    dur_max = A.max(dur);

    S = A.covariance(dur, wt); // (d, w))
    S_1 = A.inverse_2x2(S);
     
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Testing Statistics: " << endl;
    cout << "Total Number of Testing Data Points:  " << nCs.size() << endl;
    cout << "Weight:  \tmin1: " << wt_min << "\tmean1: " << wt_mean << "\tmax1: " << wt_max << endl;
    cout << "Duration:\tmin2: " << dur_min << "\tmean2: " << dur_mean << "\tmax2: " << dur_max << endl;
    cout << "Testing Categories: " << unique_testingCategories.size() << endl;
    cout << "Total Unique Categories in Test Data: " << uniqueCategories.size() << endl;
    cout << "--------------------------------------------------------------------------------" << endl;  
    cout << "Answer Size: " << the_big_answer.size() << endl;
    cout << "--------------------------------------------------------------------------------" << endl;    
    
    // output the solution
    return the_big_answer;
}
int main(int argc, char** argv) {
    ChildStuntedness CS;
    
    string readLine;
    _case hC, nC;
    vector<_case> hCs, nCs;

    string str_5, str_6, str_7, str_8, str_9, str_10, str_11, str_12;
    vector<_odv> odv_case;
    vector<string> trainingData, testingData;
    vector<double> returnData;
    
    bool output_nC = false; // boolean to output the nC vector during initialization

    // read data from test case file
    string testData = "/Users/scott/Arctria/ARC-32/ChildStuntedness/data/exampleData.csv";
    
    cout << "======================================================================================" << endl;         
    cout << "Opening test data file: " << testData << endl;
    ifstream testDataFile (testData);
    
    while (getline(testDataFile, readLine)) {
        // accumulate test data
        odv_case.clear();
        
        replace(readLine.begin(), readLine.end(), ',', ' ');

        istringstream ss(readLine);
        trainingData.push_back(readLine);
        if (!(ss >> hC.ID >> hC.age >> hC.sex >> hC.nutrition >> str_5 
                >> str_6 >> str_7 >> str_8 >> str_9 >> str_10 >> str_11 >> str_12 
                >> hC.weight >> hC.duration)) {
            break;
        };
        ss.str("");
        
        odv_case.push_back(CS.doesStringContainNA(str_5));
        odv_case.push_back(CS.doesStringContainNA(str_6));        
        odv_case.push_back(CS.doesStringContainNA(str_7));        
        odv_case.push_back(CS.doesStringContainNA(str_8));
        odv_case.push_back(CS.doesStringContainNA(str_9));
        odv_case.push_back(CS.doesStringContainNA(str_10));
        odv_case.push_back(CS.doesStringContainNA(str_11));
        odv_case.push_back(CS.doesStringContainNA(str_12));    
        
        hC.odv = odv_case;
        hCs.push_back(hC);
        
        nC.ID = hC.ID;
        nC.age = hC.age;
        nC.sex = hC.sex;
        nC.nutrition = hC.nutrition;
        nC.odv = hC.odv;
        nC.weight = 0.0; // to be determined
        nC.duration = 0.0; // to be determined
        nCs.push_back(nC);
        
        if (output_nC) {
            cout << "nC: " << nC.ID
                    << "\tage: " << nC.age
                    << "\tsex: " << nC.sex
                    << "\tnutrition: " << nC.nutrition
                    << "\todv5: " << nC.odv.at(0).odv_value
                    << "\todv6: " << nC.odv.at(1).odv_value
                    << "\todv7: " << nC.odv.at(2).odv_value
                    << "\todv8: " << nC.odv.at(3).odv_value
                    << "\todv9: " << nC.odv.at(4).odv_value
                    << "\todv10: " << nC.odv.at(5).odv_value
                    << "\todv11: " << nC.odv.at(6).odv_value
                    << "\todv12: " << nC.odv.at(7).odv_value
                    << "\tweight: " << nC.weight
                    << "\tduration: " << nC.duration
                    << endl;
        }
    }
    
    cout << "Number of data points: " << hCs.size() << endl;
 
    // generate test data string vector
    for (int i = 0; i < nCs.size(); i++) {
        nC = nCs.at(i);
        ostringstream ss;
        ss      << nC.ID << ","
                << nC.age << ","
                << nC.sex << ","
                << nC.nutrition << ","
                << nC.odv.at(0).str << "," // 5
                << nC.odv.at(1).str << "," // 6
                << nC.odv.at(2).str << "," // 7
                << nC.odv.at(3).str << "," // 8
                << nC.odv.at(4).str << "," // 9
                << nC.odv.at(5).str << "," // 10
                << nC.odv.at(6).str << "," // 11
                << nC.odv.at(7).str << ","; // 12
        testingData.push_back(ss.str());
        ss.str("");
    }

    cout << "Sending training and test data to ChildStuntedness class" << endl;
    returnData = CS.predict(trainingData, testingData);
    
    cout << "Number of data pairs returned: " << returnData.size() / 2 << endl;
    
    vector<double> w;
    vector<double> d;
    vector<double> w_real;
    vector<double> d_real;    
    vector<int> IDList;

    for (int k = 0; k < returnData.size(); k=k+2) {
        d.push_back(returnData.at(k));
        w.push_back(returnData.at(k+1));
    } 
    
    for (int j = 0; j < hCs.size(); j++) {
        bool isInList = false;

        for (int i = 0; i < IDList.size(); i++) {
            if (hCs.at(j).ID == IDList.at(i)) {
                isInList = true;
            }
        }
        if (!isInList) {
            IDList.push_back(hCs.at(j).ID);
            w_real.push_back(hCs.at(j).weight);
            d_real.push_back(hCs.at(j).duration);
        }
    }

    double score = CS.calcScore(w, w_real, d, d_real);
    
    cout << "======================================================================================" << endl;         
    cout << "Estimated Score: " << score << endl;
    cout << "======================================================================================" << endl;

    return 0;
}
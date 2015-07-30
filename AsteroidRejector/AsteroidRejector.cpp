/* =========================================================================================================== *
 * NASA Tournament Labs / TopCoder / Planetary Resources
 * Phase I Marathon Match Competition
 * AsteroidRejector v11 [Competition Sensitive]
 * RigelFive - Hudson Ohio USA
 * 18 April 2014 - 3 May 2014
 * =========================================================================================================== *
 * version history:
 * v0:  Tested on TopCoder website successfully 0945EDT 21 April.
 * v1:  Failed TopCoder website test with error due to an out of range condition on a vector.
 * v2:  Segmentation errors are occurring with no specific indications... Java test program is failing.
 * v3:  Segmentation errors are continuing to occur with no specific indications.
 * v4:  Cannot determine the data interfaces with the TopCoder test program.  Need to redevelop local test code.
 * v5:  Successful in generating a first score and completing the algorithm.
 * v6:  Increased score by 120% from simplified algorithm in v5. =D. Submitting for v/t of official score system.
 * v7:  Added NN, OFS, and other input filters.  Revised as the <20 7 1> topology run time exceeds 40 minutes.
 * v8:  Reduced program size to analyze basic image and center image.  Revised NN to <10 5 1> topology.
 * v9:  Modified the code to have a similar structure to v6. 
 * v10: Modified the code to simplify the inputs even further using four inputs.
 * v11: Reconfigured the NN to a <16 20 4>.
 * 
 * =========================================================================================================== *
 * Atomized Software Specification:
 * Customer Requirements:
 *  1. Identify presence of Near Earth Objects [NEO] as small or smaller than 20m with CCD imagery/PA radio telescopes.
 *  2. Improve the process flow to validate imagery taken from ground CCD telescope observatories.
 *  3. Develop advanced open source software to provide a solution before 3 May 2014.
 *  4. Output a sorted list of UniqueIDs "least probable to be rejected objects at the front of the array, 
 *     and those the most probable (to be rejected) at the back."
 * 
 * Functional Requirements:
 *  1. Initialize AsteroidRejector class.
 *  2. Activate [trainingData] function to train a processor to identify images with objects.
 *  3. Activate [testingData] function to identify images with rejected objects.
 *  4. Activate [getAnswer] function to return a sorted list of identified objects.
 *  5. Perform data input operations to provide information for processing.
 *  6. Perform image processing on imagery.
 *  7. Perform processing to determine if objects are existent/non-existent in images.
 *  8. Perform data output operations to generate information required.
 *  9. Terminate AsteroidRejector class.
 * 
 * Interface Requirements:
 *  1. Example image files (for offline testing method only)
 *  2. Example detection object (.dets) data files.
 *   
 * Performance Requirements:
 *  1. Code must complete within 40 minutes with an objective of 30 minutes or less.
 *  2. Utilize 5 GB of memory or less.
 *  3. Compile time of 60 seconds or less.
 *  4. Source code cannot exceed 1MB in length.
 * 
 * Design Requirements:target
 *  1. Utilize Amazon EC2 cloud computer environment.
 *  2. Develop code using GNU g++ 4.8.1 w/-std=c++11 and -O2 compiler declarations.
 *  3. Complete all code development in sufficient time to be finished by 2 May 2014.
 *  4. 16-bit extracted images supplied at 64w x 64h px with the detected object centered in the image.
 *  5. Any data bit in the raw image that falls outside of the range is assigned a value of 65535 in the extracted image.
 *  6. Code must be written as a single file.   5
 * 
 * Quality Requirements:
 *  1. Maximize overall score of the algorithm (objective: 1,000,000 points).
 * 
 * Verification Requirements:
 *  1. The verification of the program may be performed using a Java testing program.
 *  2. Verification may be performed using the TopCoder "Test Program" feature to ensure code compiling.
 * 
 * =========================================================================================================== *
 */

#define _WHAT_DO_YOU_WANT_TO_DOOOO_WITH_YOUR_LIFE_ 2014;
#define _I_WANT_TO_RAWK_ 1;

#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <tuple>
#include <algorithm>
#include <cstdio>
#include <ctime>

using namespace std;
typedef tuple<int,double> object_tuple;

bool sort_object_score (const object_tuple &lhs, const object_tuple &rhs){
  return get<1>(lhs) < get<1>(rhs); 
}
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
double Neuron::eta = 0.150;  
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

class AsteroidRejector {
public:
    ////////////////////////////
    // Methods:
    //////////////////////////// 
    
    // AsteroidRejector class methods
    AsteroidRejector();
    int trainingData(vector<int> imageData, vector<string> detections);
    int testingData(vector<int> imageData, vector<string> detections);
    vector<int> getAnswer();
 
    // Statistic methods
    double sum(vector<int> &image);
    double average(vector<int> &image);
    double stddev2(vector<int> &image);
    double min(vector<int> &image);
    double max(vector<int> &image);
    double max_vector(vector<double> &vector_set);
    double average_vector(vector<double> &vector_data);    
    
    // Geometric methods
    double line(double m, double b, double x);
    double parabola(double a, double b, double c, double x);
    
    // Detections and image data methodstrainingData initiating at
    void acquire_uniqueID(vector<string> &detections); 
    void acquire_target_data(vector<string> &detections);
    void calculate_training_set(vector<int> &imageData);    

    // Image methods
    int image_index(int x, int y, int sqrt_pixels);
    vector<int> create_image(vector<int> &imageData);
    vector<int> create_central_image(vector<int> &image);

    // Filter methods
    double parabola_filter(double x, double x_mean, double x_max);
       
    // Process methods
    void calculate_scores(vector<int> &imageData);
    void calculate_NN_testing_scores(vector<int> &imageData);    
    void pick_best_NN_scores(vector<int> &ID_vector, vector<double> &score_vector);
    
    // NN methods
    void buildNNtopology(int input, int middle, int output);    
    void feedForward(const vector<double> &inputVals);
    void backProp(const vector<double> &targetVals);
    void getResults(vector<double> &resultVals) const;
    double getRecentAverageError(void) const { return m_recentAverageError; }
    void testingNNalgorithm(void);
    
    // Output methods    
    void ten_vector_output(vector<double> &vector_set);
    void ten_vector_int_output(vector<int> &vector_set);  
    void output_test_parameters(void);
    void output_training_parameters(void);   

    ////////////////////////////
    // Variables:
    ////////////////////////////

    // Asteroid Rejector Parameters:
    vector<int> imageData;
    vector<string> detections;    
    
    bool startTraining = false;
    bool startTesting = false;
    bool startAnswer = false;   
    bool topology_built = false;
        
    // Program Timing Parameters
    double training_start_time = 0.0;
    double testing_start_time = 0.0;
    double get_answer_start_time = 0.0;
    double training_elapsed_time = 0.0;
    double testing_elapsed_time = 0.0;
        
    // Image Parameters:
    vector<int> image;
    vector<int> center_image;
    
    // NN Parameters:
    vector<unsigned> topology;
    vector<double> inputVals, targetVals, resultVals;   
    
    // Training Parameters:
    vector<double> tr_sum, tr_avg, tr_stddev2, tr_reject;
    vector<double> tr_sum_center, tr_avg_center, tr_stddev2_center;
    
    // Training Statistics
    vector<double> st_sum, st_avg, st_stddev2;
    vector<double> st_sum_center, st_avg_center, st_stddev2_center;
    
     // Testing Parameters:   
    vector<double> test_sum_center, test_avg_center, test_stddev2_center;
    vector<double> test_sum, test_avg, test_stddev2;
    vector<int> test_uniqueID;       
    vector<double> NN_score;  
    
    // Output Parameters:
    int images = 0;    
    int detection_items = 0;   
    
    // Scoring Parameters:
    vector<int> best_uniqueID;
    vector<double> best_score;  
       
    // Initial Training Statistics - Test 56
    // Note:  Values are for reference only and are recalculated before use.
    double sum_image_mean = 1.79488e+07;
    double sum_image_max = 2.68366e+08;
    double avg_image_mean = 4381.61;
    double avg_image_max = 65535;
    double stddev2_image_mean = 4.66526e+07;
    double stddev2_image_max = 1.02903e+09;

    double sum_center_mean = 1.04848e+06;
    double sum_center_max = 1.6777e+07;
    double avg_center_mean = 4095.13;
    double avg_center_max = 65535;
    double stddev2_center_mean = 1.74199e+07;
    double stddev2_center_max = 1.05302e+09;
   
private:
    vector<Layer> m_layers;
    double m_error;
    double m_recentAverageError;
    static double m_recentAverageSmoothingFactor;
};

AsteroidRejector::AsteroidRejector() { 
}
double AsteroidRejector::sum(vector<int> &image) {
    uint64_t sum = 0.0;
    
    int imageSIZE = image.size();     
    for (int i = 0; i < imageSIZE; i++) {
        sum = sum + uint64_t(abs(image.at(i)));
    }
    return double(sum);
}
double AsteroidRejector::average(vector<int> &image) {  
    uint64_t sum = 0.0;
    
    int imageSIZE = image.size();     
    for (int i = 0; i < imageSIZE; i++) {
        sum = sum + abs(image.at(i));
    }
    return double(sum/imageSIZE);
}
double AsteroidRejector::stddev2(vector<int> &image) {
    double s2;
    uint64_t sum = 0.0;
    double sum2 = 0.0;  
    double avg;
    
    int imageSIZE = image.size();
    
    for (int i = 0; i < imageSIZE; i++) {
        sum = sum + abs(image.at(i));
    }
    avg = double(sum / imageSIZE);
    
    for (int i = 0; i < imageSIZE; i++) {
        sum2 = sum2 + pow(((double) image.at(i) - avg),2.0);
    }
    
    s2 = (1.0/((double) imageSIZE-1.0)) * abs(sum2);
    
    return s2;    
}
double AsteroidRejector::min(vector<int> &image) {
    int min1 = 65536;
    int imageSIZE = image.size();
    
    for (int i = 0; i < imageSIZE; i++) {
        if (image.at(i) < min1) { min1 = image.at(i); };
    }

    return double(min1);    
}
double AsteroidRejector::max(vector<int> &image) {
    int max1 = 0;

    int imageSIZE = image.size();
    
    for (int i = 0; i < imageSIZE; i++) {
        if (image.at(i) > max1) { max1 = image.at(i); };
    }

    return double(max1);
}
double AsteroidRejector::max_vector(vector<double> &vector_set) {
    double max1 = 0;
    
    for (int i = 0; i < vector_set.size(); i++) {
        if (vector_set.at(i) > max1) { max1 = vector_set.at(i); };
    }

    return double(max1);
}
double AsteroidRejector::average_vector(vector<double> &vector_set) {
    double sum;
    double average;
    sum = 0;
    for (int i = 0; i < vector_set.size(); i++) {
        sum = sum + vector_set.at(i);
    }
    average = sum / ((double) vector_set.size());
    return average;
}
double AsteroidRejector::line(double m, double b, double x) {
    double y;

    y = (m*x) + b;
    
    return y;
}
double AsteroidRejector::parabola(double a, double b, double c, double x) {
    double y;
    
    y = (a*pow((x-b),2.0)) + c;
    
    return y;
}
double AsteroidRejector::parabola_filter(double x, double x_mean, double x_max) {
    double y, norm_x;
    double a1, a2, b, c1, c2;
    
    norm_x = x / x_max;
    
    a1 = 1.0 / pow(x_mean,2.0);
    a2 = 1.0 / pow((x_max-x_mean),2.0);
    b = x_mean;
    c1 = 0.0;
    c2 = 0.0;

    if (x < 0.0) {
        y = -1.0;
    } else if (norm_x > 1.0) {
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
int AsteroidRejector::image_index(int x, int y, int sqrt_pixels) {
    int index;
    index = x+y*sqrt_pixels;
    return index;
}
vector<int> AsteroidRejector::create_central_image(vector<int> &image) {
    vector<int> central_image;
    int pos, pixel;
    
    for (int x = 24; x <= 39; x++) {
        for (int y = 24; y <= 39; y++) {
            pos = image_index(x,y,64);
            pixel = image.at(pos);
            central_image.push_back(pixel);
        }
    }
    return central_image;
}
void AsteroidRejector::ten_vector_output(vector<double> &vector_set) {
    for (int i = 0; i < 10; i++) {
        cout << " " << vector_set.at(i);
    }
    cout << endl;
}
void AsteroidRejector::ten_vector_int_output(vector<int> &vector_set) {
    for (int i = 0; i < 10; i++) {
        cout << " " << vector_set.at(i);
    }
    cout << endl;
}
void AsteroidRejector::buildNNtopology(int input, int middle, int output) {
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
double AsteroidRejector::m_recentAverageSmoothingFactor = 100.0;
void AsteroidRejector::getResults(vector<double> &resultVals) const {
    resultVals.clear();
    for (unsigned n = 0; n < m_layers.back().size() - 1; ++n) {
        resultVals.push_back(m_layers.back()[n].getOutputVal());
    }
}
void AsteroidRejector::backProp(const vector<double> &targetVals) {
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
void AsteroidRejector::feedForward(const vector<double> &inputVals) {
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
void AsteroidRejector::acquire_target_data(vector<string> &detections) {
    string line;
    const char *char_line;

    for (int i = 0; i < detections.size(); i++) {
        line = detections.at(i);
        
        string ID;
        string DTN;
        string FN;
        string SN;
        string JD;
        string RA;
        string DEC;
        string X;
        string Y;
        string MAG;
        string FWHM;
        string EL;
        string TH;
        string RMSE;
        string DMU;
        int RJCT;
        stringstream ss;   
        if (line != "") {
            ss.str(line);
            ss>>ID>>DTN>>FN>>SN>>JD>>RA>>DEC>>X>>Y>>MAG>>FWHM>>EL>>TH>>RMSE>>DMU>>RJCT;  
            tr_reject.push_back(RJCT); 
            ss.str("");             
        }
        else {
            cout << "WHOOPS!  Bad detection data line entry->" << line << endl;
        }
    }
}
void AsteroidRejector::calculate_training_set(vector<int> &imageData) {
    int pixel;
    int pixel_count = 0;
    int images = 0;

    for (int i = 0; i < imageData.size(); i++) {
        
        if ((imageData.at(i) >= 0) && (imageData.at(i) <= 65536)) {
            pixel = imageData.at(i);
        } else if (imageData.at(i) < 0) {
            pixel = 0;
        }
        else if (imageData.at(i) > 65536) {
            pixel = 65536;
        }        
        image.push_back(pixel);
        
        if (pixel_count == 4095) {    
            images++;         
            center_image = create_central_image(image);
            
            double sum_image = sum(image);
            double avg_image = average(image);
            double stddev2_image = stddev2(image);

            double sum_center_image = sum(center_image);
            double avg_center_image = average(center_image);
            double stddev2_center_image = stddev2(center_image);

            st_sum.push_back(sum_image);
            st_avg.push_back(avg_image);
            st_stddev2.push_back(stddev2_image);             
            st_sum_center.push_back(sum_center_image);
            st_avg_center.push_back(avg_center_image);
            st_stddev2_center.push_back(stddev2_center_image);

            tr_sum.push_back(abs(parabola_filter(sum_image, sum_image_mean, sum_image_max)));
            tr_avg.push_back(parabola_filter(avg_image, avg_image_mean, avg_image_max));
            tr_stddev2.push_back(-parabola_filter(stddev2_image, stddev2_image_mean, stddev2_image_max));
  
            tr_sum_center.push_back(abs(parabola_filter(sum_center_image, sum_center_mean, sum_center_max)));
            tr_avg_center.push_back(parabola_filter(avg_center_image, avg_center_mean, avg_center_max));
            tr_stddev2_center.push_back(-parabola_filter(stddev2_center_image, stddev2_center_mean, stddev2_center_max));          

            pixel_count = 0;
            image.clear();
            center_image.clear();            
        } 
        pixel_count++;        
    }
    
    // Calculate statistics on images
    sum_image_mean = average_vector(st_sum);
    sum_image_max = max_vector(st_sum);
    avg_image_mean = average_vector(st_avg);
    avg_image_max = max_vector(st_avg);
    stddev2_image_mean = average_vector(st_stddev2);
    stddev2_image_max = max_vector(st_stddev2);
    
    sum_center_mean = average_vector(st_sum_center);
    sum_center_max = max_vector(st_sum_center);
    avg_center_mean = average_vector(st_avg_center);
    avg_center_max = max_vector(st_avg_center);
    stddev2_center_mean = average_vector(st_stddev2_center);
    stddev2_center_max = max_vector(st_stddev2_center);
      
}
void AsteroidRejector::output_training_parameters(void) {
    cout << "==========================================================================================================" << endl; 
    cout << "Calculated Statistical Image Parameters: " << endl << endl;
    cout << "double sum_image_mean = " << sum_image_mean << ";" << endl;
    cout << "double sum_image_max = " << sum_image_max << ";" << endl;
    cout << "double avg_image_mean = " << avg_image_mean << ";" << endl;
    cout << "double avg_image_max = " << avg_image_max << ";" << endl;
    cout << "double stddev2_image_mean = " << stddev2_image_mean << ";" << endl;
    cout << "double stddev2_image_max = " << stddev2_image_max << ";" << endl << endl;;
    
    cout << "double sum_center_mean = " << sum_center_mean << ";" << endl;
    cout << "double sum_center_max = " << sum_center_max << ";" << endl;
    cout << "double avg_center_mean = " << avg_center_mean << ";" << endl;
    cout << "double avg_center_max = " << avg_center_max << ";" << endl;
    cout << "double stddev2_center_mean = " << stddev2_center_mean << ";" << endl;
    cout << "double stddev2_center_max = " << stddev2_center_max << ";" << endl << endl;

    cout << "==========================================================================================================" << endl;
    cout << "Size and content of primary vectors: " << endl;
    cout << "tr_sum: " << tr_sum.size() << " ---> ";
    ten_vector_output(tr_sum);
    cout << "tr_avg: " << tr_avg.size() << " ---> ";
    ten_vector_output(tr_avg);    
    cout << "tr_stddev2: " << tr_stddev2.size() << " ---> ";
    ten_vector_output(tr_stddev2);       
    
    cout << endl;
    cout << "tr_sum_center: " << tr_sum_center.size() << " ---> ";
    ten_vector_output(tr_sum_center);
    cout << "tr_avg_center: " << tr_avg_center.size() << " ---> ";
    ten_vector_output(tr_avg_center);
    cout << "tr_stddev2_center: " << tr_stddev2_center.size() << " ---> ";
    ten_vector_output(tr_stddev2_center);

    cout << endl;
    cout << "tr_reject: " << tr_reject.size() << " ---> ";
    ten_vector_output(tr_reject);
    cout << endl;
    cout << "st_sum:            " << st_sum.size() << " ---> ";
    ten_vector_output(st_sum);
    cout << "st_avg:            " << st_avg.size() << " ---> ";
    ten_vector_output(st_avg);
    cout << "st_stddev2:        " << st_stddev2.size() << " ---> ";
    ten_vector_output(st_stddev2);

    cout << endl;
    cout << "st_sum_center:     " << st_sum_center.size() << " ---> ";
    ten_vector_output(st_sum_center);
    cout << "st_avg_center:     " << st_avg_center.size() << " ---> ";
    ten_vector_output(st_avg_center);
    cout << "st_stddev2_center: " << st_stddev2_center.size() << " ---> ";
    ten_vector_output(st_stddev2_center);

    cout << endl;
    
    cout << "images: " << images << endl;
    cout << "detections: " << detection_items << endl;
}
void AsteroidRejector::calculate_NN_testing_scores(vector<int> &imageData) { 
    int image_ctr = 0;
    int pixel;
    int pixel_count = 0;
    int trainingPass = 0;
    int training_time = 10000000;  // 35 min x 60 sec x 1000000 microsec x 1/500 training cases
    
    // create a neural network
    if (m_layers.size() != 3) {
        cout << "Building a <16 20 4> Neural Network..." << endl;         
        buildNNtopology(16,20,4);
        topology_built = true;
    }    

    int NN_training_entry_time = clock();
    int NN_training_exit_time = NN_training_entry_time + training_time;
    
    // generate optimized NN based on training data
    while (clock() < NN_training_exit_time) {    
        for (int nn = tr_sum.front(); nn < tr_sum.size(); nn=nn+4) {          
            inputVals.clear();
            targetVals.clear();
            resultVals.clear();            
            inputVals.push_back(tr_avg.at(nn));
            inputVals.push_back(tr_stddev2.at(nn));   
            inputVals.push_back(tr_avg_center.at(nn));
            inputVals.push_back(tr_stddev2_center.at(nn));
            
            inputVals.push_back(tr_avg.at(nn+1));
            inputVals.push_back(tr_stddev2.at(nn+1));   
            inputVals.push_back(tr_avg_center.at(nn+1));
            inputVals.push_back(tr_stddev2_center.at(nn+1)); 
            
            inputVals.push_back(tr_avg.at(nn+2));
            inputVals.push_back(tr_stddev2.at(nn+2));   
            inputVals.push_back(tr_avg_center.at(nn+2));
            inputVals.push_back(tr_stddev2_center.at(nn+2)); 
            
            inputVals.push_back(tr_avg.at(nn+3));
            inputVals.push_back(tr_stddev2.at(nn+3));   
            inputVals.push_back(tr_avg_center.at(nn+3));
            inputVals.push_back(tr_stddev2_center.at(nn+3));             
            
            targetVals.push_back(tr_reject.at(nn));
            targetVals.push_back(tr_reject.at(nn+1));    
            targetVals.push_back(tr_reject.at(nn+2));    
            targetVals.push_back(tr_reject.at(nn+3));                

            feedForward(inputVals);
            getResults(resultVals);
            assert(targetVals.size() == topology.back());
            backProp(targetVals);                    
        }
        trainingPass++;
    }  
    
    // calculate scores based on test data
    for (int i = 0; i < imageData.size(); i++) {
        
        if ((imageData.at(i) >=0) && (imageData.at(i) <= 65536)) {
            pixel = imageData.at(i);
        } else if (imageData.at(i) < 0) {
            pixel = 0;
        }
        else if (imageData.at(i) > 65536) {
            pixel = 65536;
        }           
        
        pixel = imageData.at(i);
        image.push_back(pixel);
         
        if (pixel_count == 4095) {    
            images++;
            
            if (image_ctr < 3) {
                center_image = create_central_image(image);

                inputVals.push_back(parabola_filter(average(image), avg_image_mean, avg_image_max));
                inputVals.push_back(-parabola_filter(stddev2(image), stddev2_image_mean, stddev2_image_max));
                inputVals.push_back(parabola_filter(average(center_image), avg_center_mean, avg_center_max));
                inputVals.push_back(-parabola_filter(stddev2(center_image), stddev2_center_mean, stddev2_center_max));
                image_ctr++; 
                pixel_count = 0;
                image.clear();
                center_image.clear();                 
            } else if (image_ctr == 3) {
                center_image = create_central_image(image);

                inputVals.push_back(parabola_filter(average(image), avg_image_mean, avg_image_max));
                inputVals.push_back(-parabola_filter(stddev2(image), stddev2_image_mean, stddev2_image_max));
                inputVals.push_back(parabola_filter(average(center_image), avg_center_mean, avg_center_max));
                inputVals.push_back(-parabola_filter(stddev2(center_image), stddev2_center_mean, stddev2_center_max));

                feedForward(inputVals);
                getResults(resultVals);
                
                NN_score.push_back(resultVals.at(0));
                NN_score.push_back(resultVals.at(1));
                NN_score.push_back(resultVals.at(2));
                NN_score.push_back(resultVals.at(3));  
                
                image_ctr=0; 
                pixel_count = 0;
                image.clear();
                center_image.clear();                 
                inputVals.clear();
                resultVals.clear();                
            }              
        }
        pixel_count++;         
    }      
}
void AsteroidRejector::acquire_uniqueID(vector<string> &detections) {
    string line;
    const char *char_line;
    int ID;
    
    for (int i = 0; i < detections.size(); i++) {
        detection_items++;
        line = detections.at(i);
        char_line = new char(line.length()+1);
        char_line = line.c_str();        
        sscanf(char_line,"%i %*s", &ID);
        test_uniqueID.push_back(ID);
    }
}
void AsteroidRejector::pick_best_NN_scores(vector<int> &ID_vector, vector<double> &score_vector) {
    double score[] = {0.0, 0.0, 0.0, 0.0};
    double avg_top_four_scores, sum_score;
    int uniqueID_val;
    int ctr = 0;
    for (int i = 0; i < score_vector.size(); i++) {
        score[ctr] = score_vector.at(i);
        ctr++;
        if (ctr == 4) {            
            sum_score = score[0] + score[1] + score[2] + score[3];  
            avg_top_four_scores = sum_score/4.0;
            uniqueID_val = ID_vector.at(i);
            best_score.push_back(avg_top_four_scores);
            best_uniqueID.push_back(uniqueID_val);
            ctr = 0;
        }
    }
}
void AsteroidRejector::output_test_parameters(void) {
    
    cout << "==========================================================================================================" << endl; 
    cout << "Testing score vector: " << NN_score.size() << " ---> ";
    ten_vector_output(NN_score);    
    cout << "Testing UniqueID vector: " << test_uniqueID.size() << " ---> ";
    ten_vector_int_output(test_uniqueID);
    cout << endl;
    cout << "Best Scores UniqueID vector: " << best_uniqueID.size() << " ---> ";
    ten_vector_int_output(best_uniqueID); 
    cout << "Best Scores vector: " << best_score.size() << " ---> ";   
    ten_vector_output(best_score);
    cout << endl;
}
int AsteroidRejector::trainingData(vector<int> imageData, vector<string> detections) {
    
    if (startTraining == false) {
        training_start_time = clock(); 
        cout << "trainingData start time: " << (training_start_time/ 1000000.0) << " sec" << endl;        
        cout << "==========================================================================================================" << endl;                
        startTraining = true;
    }   
    
    acquire_target_data(detections);  
    calculate_training_set(imageData);    
    return _I_WANT_TO_RAWK_;
}
int AsteroidRejector::testingData(vector<int> imageData, vector<string> detections) {
    
    if (startTesting == false) {
        testing_start_time = clock();
        training_elapsed_time = (testing_start_time - training_start_time) / 1000000.0; 
        cout << "trainingData elapsed time: " << (training_elapsed_time / 1000000.0) << " sec" << endl;
        cout << "testingData starting time: " << (testing_start_time / 1000000.0) << " sec" << endl;         
        cout << "==========================================================================================================" << endl;        
        startTesting = true;
    }

    calculate_NN_testing_scores(imageData);
    acquire_uniqueID(detections);
    return _I_WANT_TO_RAWK_;
}
vector<int> AsteroidRejector::getAnswer() {
    int ID_val;
    double SCORE_val;
    
    vector<int> the_big_answer;  
    vector<double> the_best_scores; 
    vector<object_tuple> test_tuple;   

    if (startAnswer == false) {
        output_training_parameters();        
        get_answer_start_time = clock();
        testing_elapsed_time = (get_answer_start_time - testing_start_time) / 1000000.0;
        cout << "getAnswer start time: " << (get_answer_start_time / 1000000.0) << " sec" << endl;        
        startAnswer = true;
    }

    pick_best_NN_scores(test_uniqueID, NN_score);

    // create a tuple with <uniqueID, test_score>
    for (int i = 0; i < best_uniqueID.size(); i++) {
        ID_val = best_uniqueID.at(i);
        SCORE_val = best_score.at(i);
        test_tuple.push_back(make_tuple(ID_val, SCORE_val));        
    }
    
    // sort and output the tuple
    sort(test_tuple.begin(),test_tuple.end(),sort_object_score);
    
    // creating a local vector to transfer sorted data to object in main;
    for(vector<object_tuple>::iterator iter = test_tuple.begin(); iter != test_tuple.end(); iter++){
        the_big_answer.push_back(get<0>(*iter)); 
        the_best_scores.push_back(get<1>(*iter));
    } 
    
    output_test_parameters();  
   
    cout << "AsteroidRejector Algorithm Complete: " << _I_WANT_TO_RAWK_ ;
    cout << endl;    
    
    return the_big_answer;    
};

/* ChildStuntedness x4 [Competition Sensitive]
 * TopCoder
 * RigelFive - Hudson Ohio USA
 * 11 August 2014 - 21 August 2014
 * ================================================================
 * Note:  results from this version were able to generate 
 * a score of 474348 using the provided test data offline.
 * This score was calculated with a new covariance matrix, however
 * a higher score of 479880 was generated using the covariance
 * matrix in the ChildStuntedness problem statement.
 * ================================================================
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
typedef tuple<int,double> spearman_tuple;

bool sort_object_score (const object_tuple &lhs, const object_tuple &rhs){
  return get<1>(lhs) < get<1>(rhs);  // sorts for smallest numbers at beginning of list (output of result to main)
}
bool sort_spearman (const spearman_tuple &lhs, const spearman_tuple &rhs){
  return get<1>(lhs) < get<1>(rhs);  // sorts for smallest numbers at beginning of list (output of result to main)
}

struct _odv {
    bool isNA;
    double odv_value;
    string str;
};
struct _case { // a case is a single data point with a collection of parameters
    int ID;             // 1. Unique Fetus ID
    double age;         // 2. Estimated fetus gestational age from last menstrual recall date
    double sex;         // 3. 0 = Male, 1 = Female
    double nutrition;   // 4. Maternal nutritional status (1 or 2)
    vector<_odv> odv;   // 5-12. Dependent variables: Ultrasound observed measurements
    double weight;      // 13. Birth Weight (w)
    double duration;    // 14. Pregnancy Duration, or Birthday (b)
    
    int category;       // defines the category for the data point
};
struct _result {
    int ID;
    double weight;
    double duration;
};
struct _coeff { // structure for wN and dN coefficients    
    vector<double> N; // vector of coefficients. 
    double N_max; // maximum limit for coefficients in GAsearch
    double N_min; // minimum limit for coefficients in GAsearch
};
struct _matrix {
    int row;
    int col;
    double val;
};
struct _w_box {
    double weight;
    vector<int> categories;
};
struct _d_box {
    double duration;
    vector<int> categories;
};
struct _cat_w {
    int category;
    double weight;
};
struct _cat_d {
    int category;
    double duration;
};
/*
 * Analysis class
 */
class Analysis {
public:
    double sum(vector<double> v);
    int sum_integer(vector<int> v_i);
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
    vector<double> spearman(vector<double> x, vector<double> y);
    
private:
    
};
int Analysis::sum_integer(vector<int> v_i) {
    int sum = 0;
    
    int vSIZE = v_i.size();     
    for (int i = 0; i < vSIZE; i++) {
        sum = sum + double(v_i.at(i));
    }
    return sum;
}
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
    double interval_50pct = 0.59 * stddev; // 0.674490 /// 0.59
    double interval_20pct = 0.9069 * stddev; // 1.281552 /// 0.9069
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
vector<double> Analysis::spearman(vector<double> x, vector<double> y) {
    // Calculates Spearman's coefficient of correlation
    // Note:  length of x and y vectors needs to be equal
    double r;
    vector<spearman_tuple> x_tuple, y_tuple; 
    vector<int> rank_x, rank_y;
    vector<double> ranked_x, ranked_y, d, d2, output;
    
    // load the x values into a tuple
    for (int i = 0; i < x.size(); i++) {
        x_tuple.push_back(make_tuple(i, x.at(i)));
        y_tuple.push_back(make_tuple(i, y.at(i)));
    }
    
    // sort the x & y tuples
    sort(x_tuple.begin(),x_tuple.end(),sort_spearman); 
    sort(y_tuple.begin(),y_tuple.end(),sort_spearman); 

    // load the x rank values into a vector
    for(vector<spearman_tuple>::iterator iter = x_tuple.begin(); iter != x_tuple.end(); iter++){
        rank_x.push_back(get<0>(*iter)); 
        ranked_x.push_back(get<1>(*iter));
    } 
    
    // load the y rank values into a vector
    for(vector<spearman_tuple>::iterator iter = y_tuple.begin(); iter != y_tuple.end(); iter++){
        rank_y.push_back(get<0>(*iter)); 
        ranked_y.push_back(get<1>(*iter));
    } 
    
    // calculate the parameters needed to generate Spearman's coefficient
    for (int i = 0; i < rank_x.size(); i++) {
        d.push_back(rank_x.at(i)-rank_y.at(i));
        d2.push_back(d.back()*d.back());
    }
    
    double n = rank_x.size();
    double n2 = n*n;
    double sum_d2 = sum(d2);

    r = 1.0 - ((6.0 * sum_d2)/(n * (n2 - 1)));

    double t_distribution = r * sqrt((n-2.0)/(1.0-(r * r)));
    
    output.push_back(r); // Spearman's coefficient
    output.push_back(t_distribution); // t-distribution (P-value)
    
    
    return output;
}

/*
 * ChildStuntedness class
 */
class ChildStuntedness {
public:
    // methods
    vector<double> predict(vector<string> training, vector<string> testing);
    
    _odv doesStringContainNA(string str); 
    double calcScore(vector<double> w, vector<double> w_real, vector<double> d, vector<double> d_real);
    void parseData (vector<string> &dataString, bool isTraining);
    void clearParameters(void);
    void loadParameters(vector<_case> type);
    vector<double> vectorWithoutZeros(vector<double> v);
    int generateCategoryCode(vector<int> codes);
    vector<int> categorizeVector(vector<double> v, bool allowZero);
    void categorizeCase(vector<_case> Cs, bool isTraining);
    vector<double> generateResult(vector<_result> results);
    void outputStatistics(int numCats);
    void quantifyTrainingData(int numCats);
    void defineMLEValues(void);
    void generateLookupTable(void);
    
    // parameters
    _case hC, nC; // historical (training) case, new (test) case
    vector<_case> hCs, nCs;
    
    string str_5, str_6, str_7, str_8, str_9, str_10, str_11, str_12;
    vector<_odv> odv_case;
    string readLine;
    
    vector<double> age, sex, nutrition, ODV5, ODV6, ODV7, ODV8, ODV9, ODV10, ODV11, ODV12;
    vector<double> weight, duration;
    vector<_w_box> wB;
    vector<_d_box> dB;
    _w_box wb;
    _d_box db;
    
    vector<int> birthEvent, totBirths, weightEntry;
    vector<double> p_Birth;
    vector<double> t_Birth;
    vector<double> p_Weight;
    vector<double> random_values_1, random_values_2;

    vector<_cat_w> w_cat;
    vector<_cat_d> d_cat;
    
    double min_duration, max_duration, mean_duration, min_weight, max_weight, mean_weight;
    double MLE_wt, MLE_dur;
    
    double delta_d, delta_w;
    
    vector<int> training_categories, testing_categories, unique_training_categories, unique_testing_categories;
    double MLE_dur_default, MLE_wt_default;
    
private:

};
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
double ChildStuntedness::calcScore(vector<double> w, vector<double> w_real, vector<double> d, vector<double> d_real) {
    Analysis A;
    double error = 0;
    double error_0 = 0;
    
    vector<_matrix> S, S_1;

    //recalculating the covariance matrix for the specific data set v. the entire data set in the competition
    S = A.covariance(d_real, w_real); // (d, w))
    S_1 = A.inverse_2x2(S);

    double inverseS_00 = 3554.42;
    double inverseS_01 = -328.119;
    double inverseS_10 = -328.119;
    double inverseS_11 = 133.511;
    
    //double inverseS_00 = S_1.at(0).val; // 3554.42 for the complete set of data
    //double inverseS_01 = S_1.at(1).val; // -328.119 for the complete set of data
    //double inverseS_10 = S_1.at(2).val; // -328.119 for the complete set of data
    //double inverseS_11 = S_1.at(3).val; // 133.511 for the complete set of data 
    
    double meanD = A.average(d_real);
    double meanW = A.average(w_real);

    for (int i = 0; i < w.size(); i++) {
        double deltaW = w.at(i) - w_real.at(i);
        double deltaD = d.at(i) - d_real.at(i);
        double val1 = deltaD * inverseS_00 + deltaW * inverseS_10;
        double val2 = deltaD * inverseS_01 + deltaW * inverseS_11;
        double ei = (val1 * deltaD) + (val2 * deltaW);
        error += ei;
        //cout << "i: " << i << "\terror: " << error << endl;
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
    cout << "error: " << error << "\terror_0: " << error_0 << "\tscore: " << score << endl;
    return score;

}
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
void ChildStuntedness::loadParameters(vector<_case> Cs) {
    
    clearParameters();
    
    // type
    for (int i = 0; i < Cs.size(); i++) {
        age.push_back(Cs.at(i).age);
        sex.push_back((double) Cs.at(i).sex);
        nutrition.push_back((double) Cs.at(i).nutrition);
        ODV5.push_back(Cs.at(i).odv.at(0).odv_value);
        ODV6.push_back(Cs.at(i).odv.at(1).odv_value);
        ODV7.push_back(Cs.at(i).odv.at(2).odv_value);
        ODV8.push_back(Cs.at(i).odv.at(3).odv_value);
        ODV9.push_back(Cs.at(i).odv.at(4).odv_value);
        ODV10.push_back(Cs.at(i).odv.at(5).odv_value);
        ODV11.push_back(Cs.at(i).odv.at(6).odv_value);
        ODV12.push_back(Cs.at(i).odv.at(7).odv_value);
        weight.push_back(Cs.at(i).weight);
        duration.push_back(Cs.at(i).duration);
    }
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
int ChildStuntedness::generateCategoryCode(vector<int> codes) {
    int code = codes.at(0)+2;
    
    for (int i = 1; i < codes.size(); i++) {
        code += (codes.at(i)+2)*pow(5.0,i);
    }
    
    return code;
}
vector<int> ChildStuntedness::categorizeVector(vector<double> v, bool allowZero) {
    Analysis A;
    vector<double> v_no_zeros;
    vector<int> category;
    v_no_zeros = vectorWithoutZeros(v);
    
    if (!allowZero) {
        for (int i = 0; i < v.size(); i++) {
            category = A.category(v_no_zeros, v);
        }
    } else {
        for (int i = 0; i < v.size(); i++) {
            category = A.category(v, v);
        }
    }
    
    return category;

}
void ChildStuntedness::categorizeCase(vector<_case> Cs, bool isTraining) {
    Analysis A;
    vector<int> codes, categories;
    int category;
    _case c;
   
    vector<int> c_age, c_sex, c_nutn, c_odv5, c_odv6, c_odv7, c_odv8, c_odv9, c_odv10, c_odv11, c_odv12;
    c_age.clear();
    c_sex.clear();
    c_nutn.clear();
    c_odv5.clear();
    c_odv6.clear();
    c_odv7.clear();
    c_odv8.clear();
    c_odv9.clear();
    c_odv10.clear();
    c_odv11.clear();
    c_odv12.clear();
    
    c_age = categorizeVector(age, false);
    c_sex = categorizeVector(sex, true);
    c_nutn = categorizeVector(nutrition, true);
    c_odv5 = categorizeVector(ODV5, false);
    c_odv6 = categorizeVector(ODV6, false);
    c_odv7 = categorizeVector(ODV7, false);
    c_odv8 = categorizeVector(ODV8, false);
    c_odv9 = categorizeVector(ODV9, false);
    c_odv10 = categorizeVector(ODV10, false);
    c_odv11 = categorizeVector(ODV11, false);
    c_odv12 = categorizeVector(ODV12, false);
    
    unique_training_categories.clear();
    unique_testing_categories.clear();
    categories.clear();
    
    for (int i = 0; i < Cs.size(); i++) {
        codes.clear();
        codes.push_back(c_age.at(i));
        codes.push_back(c_sex.at(i)); 
        codes.push_back(c_nutn.at(i)); 
        codes.push_back(c_odv5.at(i));
        codes.push_back(c_odv6.at(i)); 
        codes.push_back(c_odv7.at(i));
        codes.push_back(c_odv8.at(i));
        codes.push_back(c_odv9.at(i));
        codes.push_back(c_odv10.at(i));
        codes.push_back(c_odv11.at(i));
        codes.push_back(c_odv12.at(i));
        category = generateCategoryCode(codes);
        categories.push_back(category);
    }    
    
    if (isTraining) {
        training_categories.clear();
        for (int i = 0; i < Cs.size(); i++) {
            hCs.at(i).category = categories.at(i);
            training_categories.push_back(categories.at(i));
        }

        unique_training_categories = A.unique(training_categories);
        
        cout << "Training Categories: " << unique_training_categories.size() << endl;
    } else {
        testing_categories.clear();
        for (int i = 0; i < Cs.size(); i++) {
            nCs.at(i).category = categories.at(i);
            testing_categories.push_back(categories.at(i));
        }
        unique_testing_categories = A.unique(testing_categories);
        
        /*cout << "Test Categories: " << unique_testing_categories.size() << endl;
        
        for (int k = 0; k < unique_testing_categories.size(); k++) {
            cout << nCs.at(k).category << " " ;
        }
        cout << endl;*/
    }
}
vector<double> ChildStuntedness::generateResult(vector<_result> results) {
    
    // initialize method
    Analysis A;
    
    // initialize parameters and vectors
    vector<double> answer, avg_wt, avg_dur;
    vector<object_tuple> index_tuple;   
    _result r;
    vector<_result> uniqueID_results, sortID_results;
    vector<int> longIDList, IDList, sortIDIndex;
    bool isInLongIDList;
    
    // identify the list of IDs to be returned
    for (int i = 0; i < results.size(); i++) {
        longIDList.push_back(results.at(i).ID);
    }
    
    // define the unique IDs in the testing data set
    IDList = A.unique(longIDList);

    
    // go thru the list of unique IDs to find the values for duration and weight
    for (int i = 0; i < IDList.size(); i++) {
        isInLongIDList = false;
        avg_wt.clear();
        avg_dur.clear();
    
        double max_wt_val=1e6;
        double max_dur_val=1e6;
        double min_wt_val=0.0;
        double min_dur_val=0.0;
        
        for (int j = 0; j < results.size(); j++) {
            if (IDList.at(i) == results.at(j).ID) {
                if (results.at(j).duration < max_dur_val) {
                    max_wt_val = results.at(j).weight; // no earlier than estimate
                    max_dur_val = results.at(j).duration; // no earlier than estimate
                }
                if (results.at(j).duration > min_dur_val) {
                    min_wt_val = results.at(j).weight; // no later than estimate
                    min_dur_val = results.at(j).duration; // no later than estimate
                }
                avg_dur.push_back(results.at(j).duration);
                avg_wt.push_back(results.at(j).weight);
            }
        }

        r.ID = IDList.at(i);
        
        // List of different approaches to estimate w and d.  Best found was 1/2 of max and min in case category.
        //double est_dur = (max_dur_val+min_dur_val+A.average(avg_dur))/3.0;
        //double est_wt = (max_wt_val+min_wt_val+A.average(avg_wt))/3.0;
        //double est_dur = (max_dur_val+min_dur_val+(2.0*A.average(avg_dur)))/4.0;
        //double est_wt = (max_wt_val+min_wt_val+(2.0*A.average(avg_wt)))/4.0;
        double est_dur = (max_dur_val+min_dur_val)/2.0;
        double est_wt = (max_wt_val+min_wt_val)/2.0;
        //double est_dur = min_dur_val;
        //double est_wt = min_wt_val;
        //double est_dur = max_dur_val;
        //double est_wt = max_wt_val;
        //double est_dur = A.average(avg_dur);
        //double est_wt = A.average(avg_wt);
        
        if ((est_dur < 1.0) && (est_dur > 0.0)) {
            r.duration = est_dur; // A.average(avg_dur);
        } else {
            r.duration = MLE_dur_default;
        }
        
        if ((est_wt < 1.0) && (est_wt > 0)) {
            r.weight = est_wt; //A.average(avg_wt);
        } else {
            r.weight = MLE_wt_default;
        }

        uniqueID_results.push_back(r);
        
    }

    // sort the answers with the lowest IDs first with a sequential output vector {d, w}
    for (int i = 0; i < uniqueID_results.size(); i++) {
        index_tuple.push_back(make_tuple(i, uniqueID_results.at(i).ID));    
    }
    
    sort(index_tuple.begin(),index_tuple.end(),sort_object_score);  

    for(vector<object_tuple>::iterator iter = index_tuple.begin(); iter != index_tuple.end(); iter++){
        sortIDIndex.push_back(get<0>(*iter));
    }  
    
    // load a sorted vector with the ordered data
    for (int i = 0; i < sortIDIndex.size(); i++) {
        int index = sortIDIndex.at(i);
        r.ID = uniqueID_results.at(index).ID;
        r.duration = uniqueID_results.at(index).duration;
        r.weight = uniqueID_results.at(index).weight;
        sortID_results.push_back(r);
    }
    
    /*cout << "ID out: ";
    for (int i = 0; i < 20; i++) {
        cout << sortID_results.at(i).ID << " ";
    }
    cout << endl;*/
    
    // generate the output vector
    for (int i = 0; i < sortID_results.size(); i++) {
        answer.push_back(sortID_results.at(i).duration);
        answer.push_back(sortID_results.at(i).weight);
    }
    
    return answer;
}
void ChildStuntedness::outputStatistics(int numCats) {
    Analysis A;
    
    /*cout << "Test ID in: ";
    for (int i = 0; i < 50; i++) {
        cout << nCs.at(i).ID << " ";
    }
    cout << endl;*/

     
    int births = 0;
    for (int i = 0; i < birthEvent.size(); i++) {
        births += birthEvent.at(i);
        totBirths.push_back(births);
    }
    
    double total_births = (double) totBirths.back();
    double prob;
    for (int i = 0; i < birthEvent.size(); i++) {
        prob = totBirths.at(i)/total_births;
        p_Birth.push_back(prob);
    }
    
    p_Birth.push_back(1.0);
    
    int total_weight = A.sum_integer(weightEntry);
    for (int i = 0; i < weightEntry.size(); i++) {
        prob = (double) weightEntry.at(i)/(double) total_weight;
        p_Weight.push_back(prob);
    }
    
    cout << "--------------------------------------------------------------------------------" << endl;   
    cout << "Number of categories for p-weight and p-duration: " << numCats << endl;
    cout << "--------------------------------------------------------------------------------" << endl;   
    cout << "Most likely weight: " << MLE_wt << endl;
    cout << "Most likely duration: " << MLE_dur << endl;
    cout << "--------------------------------------------------------------------------------" << endl;   
    cout << "Duration Data Points: " << duration.size() 
            << "\tmin: " << min_duration
            << "\tmean: " << mean_duration
            << "\tmax: " << max_duration
            << endl;
   
    cout << "Weight Data Points: " << weight.size() 
            << "\tmin: " << min_weight
            << "\tmean: " << mean_weight
            << "\tmax: " << max_weight
            << endl;
    
    cout << "--------------------------------------------------------------------------------" << endl;   
    /*cout << "Birth events: ";
    for (int i = 0; i < birthEvent.size(); i++) {
        cout << birthEvent.at(i) << " ";
    }
    cout << endl;
    
    cout << "Total births: ";
    for (int i = 0; i < totBirths.size(); i++) {
        cout << totBirths.at(i) << " ";
    }
    cout << endl;
    
    cout << "Probability of births: ";
    for (int i = 0; i < p_Birth.size(); i++) {
        cout << p_Birth.at(i) << " ";
    }
    cout << endl;
    
    cout << "Birth durations: ";
    for (int i = 0; i < t_Birth.size(); i++) {
        cout << t_Birth.at(i) << " ";
    }
    cout << endl;
    cout << "--------------------------------------------------------------------------------" << endl;   
    cout << "Number in weight: ";
    for (int i = 0; i < weightEntry.size(); i++) {
        cout << weightEntry.at(i) << " ";
    }
    cout << endl;  
   
    cout << "Probability of weight class: ";
    for (int i = 0; i < p_Weight.size(); i++) {
        cout << p_Weight.at(i) << " ";
    }
    cout << endl; 
    
    cout << "Category: ";
    for (int i = 0; i < hCs.size(); i++) {
        cout << hCs.at(i).category << " ";
    }
    cout << endl;*/
    
    cout << "r {duration v. weight}: " << A.spearman(duration, weight).at(0) << "\tt-dist: " << A.spearman(duration, weight).at(1) << endl; 
    cout << "--------------------------------------------------------------------------------" << endl;   
}
void ChildStuntedness::quantifyTrainingData(int numCats) {
    Analysis A;
// calculate basic statistical parameters for the weight and duration data
    min_duration = A.min(duration);
    max_duration = A.max(duration);
    mean_duration = A.average(duration);
    min_weight = A.min(weight);
    max_weight = A.max(weight);
    mean_weight = A.average(weight);
   
    delta_d = (max_duration-min_duration)/numCats;
    delta_w = (max_weight-min_weight)/numCats;
    
    // divide the weight and duration data into intervals for the number of categories
    double x0 = min_duration;
    double x1 = x0 + delta_d;
    double y0 = min_weight;
    double y1 = y0 + delta_w;
    
    // quantify the number of births across the span of duration data 
    for (int i = 0; i < numCats; i++) {
        int event = 0;
        int entry = 0;
        
        db.categories.clear();
        for (int j = 0; j < duration.size(); j++) {
            if ((duration.at(j) >= x0) && (duration.at(j) < x1)) {
                event++;
                db.categories.push_back(hCs.at(j).category);
            }
        }
        
        wb.categories.clear();
        for (int j = 0; j < weight.size(); j++) {
            if ((weight.at(j) >= y0) && (weight.at(j) < y1)) {
                entry++;

                wb.categories.push_back(hCs.at(j).category);
            }
        }
        
        t_Birth.push_back(x0);
        
        weightEntry.push_back(entry);
        birthEvent.push_back(event);
        
        //cout << "wB category size: " << wb.categories.size() << "\tdB category size: " << db.categories.size() << endl;
        
        wb.weight = y0 + (0.5*delta_w);
        db.duration = x0 + (0.5*delta_d);

        wB.push_back(wb);
        dB.push_back(db);
        
        //cout << "wB size: " << wB.size() << "\tdB size: " << dB.size() << endl;
        
        x0 += delta_d;
        x1 += delta_d; 
        
        y0 += delta_w;
        y1 += delta_w;
    }
    
    t_Birth.push_back(x0);
}
void ChildStuntedness::defineMLEValues(void) {
    Analysis A;
    vector<double> MLE_values;
    
    // find the most likely birth weight and duration for a default zero training data case
    double max_w_occurance = A.max_int(weightEntry);
    double max_d_occurance = A.max_int(birthEvent);
    
    for (int i = 0; i < weightEntry.size(); i++) {
        if (weightEntry.at(i) == max_w_occurance) {
            MLE_wt = weight.at(i) + (0.5 * delta_w);           
        }
    }

    for (int i = 0; i < birthEvent.size(); i++) {
        if (birthEvent.at(i) == max_d_occurance) {
            MLE_dur = duration.at(i) + (0.5 * delta_d);
        }
    }
}
void ChildStuntedness::generateLookupTable(void) {
    Analysis A;
    _cat_w wc;
    _cat_d dc;
    vector<double> w_table, d_table;
    
    // search thru each box to find the number of occurrences for each category
    for (int j = 0; j < unique_training_categories.size(); j++) {
        w_table.clear();
        for (int i = 0; i < wB.size(); i++) {
            wb = wB.at(i);
            for (int k = 0; k < wb.categories.size(); k++) {
                if (unique_training_categories.at(j) == wb.categories.at(k)) {
                    w_table.push_back(wb.weight);
                }
            }
        }
        
        wc.category = unique_training_categories.at(j);
        
        if (w_table.size() > 0) {
            wc.weight = A.average(w_table);
        } else {
            wc.weight = MLE_wt_default;
        }
        w_cat.push_back(wc);
        
        //cout << " " << wc.weight;
    }
    
    for (int j = 0; j < unique_training_categories.size(); j++) {
        d_table.clear();
        for (int i = 0; i < dB.size(); i++) {
            db = dB.at(i);
            for (int k = 0; k < db.categories.size(); k++) {
                if (unique_training_categories.at(j) == db.categories.at(k)) {
                    d_table.push_back(db.duration);
                }
            }
        }
        
        dc.category = unique_training_categories.at(j);
        
        if (d_table.size() > 0) {
            dc.duration = A.average(d_table);
        } else {
            dc.duration = MLE_dur_default;
        }
        d_cat.push_back(dc);
        
        //cout << " " << dc.duration;
    }

    cout << "w categories: " << w_cat.size() << "\td categories: " << d_cat.size() << endl;
}
vector<double> ChildStuntedness::predict(vector<string> training, vector<string> testing) {
    /* 
    Column  Variable   Type    Label/Description
      1     Id         int     Unique Fetus ID
      2     t.ultsnd   float   Estimated fetus gestational age from last menstrual recall date
      3     Sex        int     0 = Male, 1 = Female
      4     Status     int     Maternal nutritional status (1 or 2)
    5-12    Odv        float   Dependent variables: Ultrasound observed measurements
     13     Birth Sz   float   Birth Weight (w)
     14     Duration   float   Pregnancy Duration, or Birthday (b)
    */

    // initialize classes
    Analysis A;
    
    // initialize parameters
    _result r;
    vector<_result> results;
    vector<double> the_big_answer;
    
    // parse the training and testing data into data structures
    parseData(training, true); // training -> hCs  
    parseData(testing, false); // testing -> nCs

    // create vectors of the training data
    loadParameters(hCs);
    categorizeCase(hCs, true);

    // quantify Training Data
    int numCats = 56; // arbitrarily chosen (values of 16, 24, 56, and 80 seem to work with data provided offline)
    quantifyTrainingData(numCats);    
    defineMLEValues(); // Determines the w and d values to use if there is no data in class for test cases
    generateLookupTable();
    outputStatistics(numCats);
    
    MLE_wt_default = MLE_wt;
    MLE_dur_default = MLE_dur;
    
    // create vectors of the testing data
    loadParameters(nCs);
    categorizeCase(nCs, false);
    
    int w_hit = 0;
    int d_hit = 0;
    
    for (int i = 0; i < nCs.size(); i++) {
        nC = nCs.at(i);
        
        MLE_wt = MLE_wt_default;
        MLE_dur = MLE_dur_default;
        
        // note:  if the category cannot be found, the default value is from above
        for (int j = 0; j < w_cat.size(); j++) {
            if (nC.category == w_cat.at(j).category) {
                w_hit++;
                MLE_wt = w_cat.at(j).weight;
            }
        }  
        
        // note:  if the category cannot be found, the default value is from above
        for (int j = 0; j < d_cat.size(); j++) {
            if (nC.category == d_cat.at(j).category) {
                d_hit++;
                MLE_dur = d_cat.at(j).duration;
            }
        }  
        
        nCs.at(i).weight = MLE_wt;
        nCs.at(i).duration = MLE_dur;
        
        //cout << "ID: " << nC.ID << "\tdur: " << MLE_dur << "\twt: " << MLE_wt << endl;
    }
    
    //cout << "w hit: " << w_hit<< "\td hit: " << d_hit << endl;

    vector<double> w_estimated, d_estimated;
    // generate estimations for duration and weight    
    for (int i = 0; i < nCs.size(); i++) {
        r.ID = nCs.at(i).ID;
        r.duration = nCs.at(i).duration; // MLE_dur_default;
        r.weight = nCs.at(i).weight; // MLE_wt_default;
        results.push_back(r);
        w_estimated.push_back(r.weight);
        d_estimated.push_back(r.duration);
    }
    
    min_weight = A.min(w_estimated);
    max_weight = A.max(w_estimated);
    mean_weight = A.average(w_estimated);
    
    min_duration = A.min(d_estimated);
    max_duration = A.max(d_estimated);
    mean_duration = A.average(d_estimated);
    
    cout << "Test Duration Data Points: " << d_estimated.size() 
            << "\tmin: " << min_duration
            << "\tmean: " << mean_duration
            << "\tmax: " << max_duration
            << endl;
   
    cout << "Test Weight Data Points: " << w_estimated.size() 
            << "\tmin: " << min_weight
            << "\tmean: " << mean_weight
            << "\tmax: " << max_weight
            << endl;
    
    // generate result
    the_big_answer = generateResult(results);
    return the_big_answer;
}

/*
 * main
 */
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
    string testData = "/Users/scott/Arctria/ARC-32/BirthEstimator/data/exampleData.csv";
    
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
    
    cout << "======================================================================================" << endl;  
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
    
    for (int j = 0; j < w.size(); j++) {
        cout << "ID: " << IDList.at(j) 
                << "\tdur: " << d.at(j) << "\tdur(r): " << d_real.at(j) 
                << "\twt: " << w.at(j) << "\twt(r): " << w_real.at(j)
                << endl;
    }
    
    double score = CS.calcScore(w, w_real, d, d_real);       
    cout << "Estimated Score: " << score << endl;
    cout << "======================================================================================" << endl;

    return 0;
}

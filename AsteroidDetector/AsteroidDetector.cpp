/* AsteroidDetector [Competition Sensitive]
 * NASA Asteroid Grand Challenge | NASA | Planetary Resources | TopCoder
 * RigelFive - Hudson Ohio USA
 * 11 August 2014 - 8 September 2014
 */

#define _AsteroidDetector_training_ 0;
#define _AsteroidDetector_testing_ 1;

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

typedef tuple<int,double> object_tuple;
typedef tuple<int,double> spearman_tuple;

bool sort_spearman (const spearman_tuple &lhs, const spearman_tuple &rhs){
  return get<1>(lhs) < get<1>(rhs);  // sorts for smallest numbers at beginning of list (output of result to main)
}
bool sort_object_score (const object_tuple &lhs, const object_tuple &rhs){
  return get<1>(lhs) > get<1>(rhs);  // sorts for largest numbers at beginning of list
}
bool sort_GA (const object_tuple &lhs, const object_tuple &rhs){
  return get<1>(lhs) > get<1>(rhs);  // sorts for largest numbers at beginning of list (output of result to main)
}
bool BothAreSpaces(char lhs, char rhs) { return (lhs == rhs) && (lhs == ' '); }
/*
 * Structures
 */
struct Connection {
    double weight;
    double deltaWeight;
};
struct _detection {
    int uniqueID; // Unique ID - An identifier for what detected object a row belongs to.
    int frameNumber; // Frame Number - which observation is this row relevant to (1, 2, 3 or 4)
    double RA; // RA - right ascension of object in decimal hours
    double DEC; // DEC - declination in decimal degrees
    int X; // X - location in pixels of the object in the FITS image.
    int Y; // Y - location in pixels of the object in the FITS image.
    double magnitude; // Magnitude - brightness of the object in magnitudes.
    int NEO; // Neo - this value will be 1 (true) if the object is a Near Earth Object (NEO), 0 (false) otherwise.
};
struct _image {
    int width;
    int height;
    vector<int> imageData; // 16 bit image values in a vector of size width x height
    vector<string> header;
    vector<double> wcs; // 8 double values to convert (x,y) coordinates to wcs coordinates (RA,DEC)
};
struct _imageSet {
    int imageID;
    vector<_image> images;
    vector<_detection> detections;
};
struct _tile {
    int tileID;
    vector<int> imageData;
};
struct _matrix {
    int row;
    int col;
    double val;
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
struct _top { // table of parameters

    double SMA; // semi-major axis
    double ECC; // eccentricity
    double INC; // inclination
    double APG; // argument of perigee
    double RAAN; // Right ascension of the ascending node
    double TMA; // true-mean anomaly

    string YODA; // year of discovery
    string NAME;
    string Ln;
    double q;
    double Q;
    string EMOID;
    double AVM; // absolute visual magnitude = f(diameter)
    string EPOCH;
};
struct _pos {
    double r;
    double v;
};

/*
 * Astrodynamics class
 */
class Astrodynamics {
public:
    _top Earth;
    
    void calcEarthLocation(int DATE, double UT);
    double trueAnomaly(double E, double e);
    double eccentricAnomaly(double M, double e);
    double meanAnomaly(double a, double t, double M0, double T);
    double Epoch_to_JD(int backwards_date, double UT);
    
private:
    
};
void Astrodynamics::calcEarthLocation(int DATE, double UT) {
    /*
     * http://www.met.rdg.ac.uk/~ross/Astronomy/Planets.html
     */
    // calculates the location of earth in heliocentric coordinates at a specified time
    double JD = Epoch_to_JD(DATE, UT);
    double CY = (JD - 2451545.0)/100.0;
    Earth.SMA = 1.00000011 + (-0.00000005 * CY);
    Earth.ECC = 0.01671022 + (-0.00003804 * CY);
    Earth.INC = 0.00005 + ((-46.94/3600.0) * CY);
    Earth.APG = -11.26064 + ((-18228.25/3600.0) * CY);
    Earth.RAAN = 102.94719 + ((1198.28/3600.0) * CY);
    Earth.TMA = 100.46435 + ((129597740.63/3600.0) * CY); // January 1.5 2000
}
double Astrodynamics::trueAnomaly(double M, double e) {
    double nu;
    double E = eccentricAnomaly(M, e);
    nu = acos((cos(E) - e)/(1.0 - (e * cos(E))));
    return nu;
}
double Astrodynamics::eccentricAnomaly(double M, double e) {
    double E;
    int i = 0;
    double PI = acos(-1.0);
    double old_E = PI; /// just assume E = M to start
    double delta_E = 1.0;
    double new_E;
    double cos_E;
    while (delta_E > 1.0e-12) {
        cos_E = cos(old_E);
        new_E = (M - (e * ((old_E * cos_E) - sin(old_E)))) / (1.0 - (e * cos_E));
        delta_E = abs(new_E - old_E);
        old_E = new_E;
        cout << "E: " << new_E << "\tdelta E: " << delta_E << endl;
    }

    return new_E;
}
double Astrodynamics::meanAnomaly(double a, double DAYS_since_EPOCH, double M0, double EPOCH) {
    double P = pow(a, 1.5);
    double PI = acos(-1.0);
    double JD = Epoch_to_JD(EPOCH, 0.0);
    double M = (((2.0 * PI) / P) * (DAYS - JD)) + M0;
    
    return M;
}
double Astrodynamics::Epoch_to_JD(int backwards_date, double UT) {
    int Y = int(backwards_date/10000);
    int M = int((backwards_date - (Y*10000))/100);
    int D = int((backwards_date) - (M*100) - (Y*10000));
    cout << "year: " << Y << "\tmonth: " << M << "\tday: " << D << endl;
    double JD = 367.0 * Y 
    - int(7.0 * (Y + int((M + 9.0)/12.0))/4.0) 
    - int(3.0 * (int((Y+(M-9)/7)/100.0) + 1.0)/4.0) 
    + int((275 * M) / 9) + D + 1721028.5 + UT/24.0;
    
    return JD;

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
    double dot(vector<double> a, vector<double> b);
    
    double line(double m, double b, double x);
    double parabola(double a, double b, double c, double x);
    double filter(double x, double x_mean, double x_min, double x_max);
    double amplifier(double y, double x_mean, double x_min, double x_max);  

    double bilinear(vector<double> x, vector<double> y, vector<double> w, double xa, double ya);
    double distance(double x1, double y1, double x2, double y2); 
    
    vector<double> ofs_algorithm(vector<double> &y);
    
    vector<int> category(vector<double> v, vector<double> val);
    vector<int> unique(vector<int> list);
    vector<_matrix> covariance(vector<double> x, vector<double> y);
    vector<_matrix> inverse_2x2(vector<_matrix> A);
    vector<double> spearman(vector<double> x, vector<double> y);
    
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
double Analysis::dot(vector<double> a, vector<double> b) {
    double a1b1, a2b2, a3b3, a2, b2, xA, yA, zA, xB, yB, zB, angle;

    xA = a.at(0);
    yA = a.at(1);
    zA = a.at(2);
    xB = b.at(0);
    yB = b.at(1);
    zB = b.at(2);

    a1b1 = xA * xB;
    a2b2 = yA * yB;
    a3b3 = zA * zB;
    a2 = sqrt((xA*xA)+(yA*yA)+(zA*zA));
    b2 = sqrt((xB*xB)+(yB*yB)+(zB*zB));
    angle = acos((a1b1 + a2b2 + a3b3) / (a2*b2)); 
    
    return angle;  
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
double Analysis::bilinear(vector<double> x, vector<double> y, vector<double> w, double xa, double ya) {

    double x1, x2, x3;
    double y1, y2, y3;
    double w1, w2, w3;
    double det, a, b, c, ans;
    
    x1 = x.at(0);
    x2 = x.at(1);
    x3 = x.at(2);
    
    y1 = y.at(0);
    y2 = y.at(1);
    y3 = y.at(2);
    
    w1 = w.at(0);
    w2 = w.at(1);
    w3 = w.at(2);
    
    det = x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3;
    a = ((y2-y3)*w1+(y3-y1)*w2+(y1-y2)*w3) / det;
    b = ((x3-x2)*w1+(x1-x3)*w2+(x2-x1)*w3) / det;
    c = ((x2*y3-x3*y2)*w1+(x3*y1-x1*y3)*w2+(x1*y2-w2*y1)*w3) / det;
    ans = a*xa + b*ya + c;

    return ans;
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
vector<double> Analysis::ofs_algorithm(vector<double> &y) {
    vector<double> z, za, zb, zpi, g, b, pk_freq, pk_index, area;
    vector<double> zs;
    double delta, cn, cn2, cs, val_h, val_w, z_smooth;
    double val_1, val_2, val_3, val_4, n_avg, b_val;
    int n, pk, s_pk, pk1, pk2;
    double pi = 3.1415926535897932;
    bool srchpeak = true;
 
    for (int i = 0; i < y.size(); i++) {
        z.push_back(y.at(i));
    }

    n = z.size();

    for (int i = 0; i < n; i++) {
        zb.push_back(double(i/(n-1.0)));
        zpi.push_back(pi*zb.back());
        za.push_back(zb.at(i)*(z.at(n-1)-z.at(0)));
        g.push_back(z.at(i) - z.at(0) - za.at(i));           
    }

    area.push_back(0.0);
    pk = 0;
    
    for (int k = 0; k < n-1; k++) {
        b_val = 0.0;

        for (int i = 1; i < n-1; i++) {
            b_val = b_val + g.at(i)*sin(double(k+1) * zpi.at(i));
        }
        b_val = b_val * 2.0/double(n);
        b.push_back(b_val);

        if (k <= 2) {
            pk_freq.push_back(b.at(k));
            pk_index.push_back(k);
            pk++;            
        }
        else {      
            if ((abs(b.at(k-2)) < abs(b.at(k-1))) && (abs(b.at(k-1)) > abs(b.at(k)))) {
                pk_freq.push_back(b.at(k-1));
                pk_index.push_back(k-1); 
                val_h = (abs(pk_freq.at(pk))+abs(pk_freq.at(pk-1)));
                val_w = (pk_index.at(pk) - pk_index.at(pk-1));                  
                area.push_back(abs(0.5 * val_h * val_w));               
                pk1 = area.size() - 1;
                pk2 = area.size() - 2;
                
                if ((pk >= 3) && (srchpeak == true)) {
                    delta = abs(area.at(pk1)-area.at(pk2))/area.at(pk1);
                    if (delta < 0.10) {
                        s_pk = pk;
                        srchpeak = false;
                    }
                }
                pk++;
            }
        }
    }
    
    val_1 = 0.0;
    val_2 = 0.0;

    for (int s = s_pk-2; s <= s_pk+1; s++) {
        val_1 = val_1 + (abs(pk_freq.at(s)) * (1.0/((double) pow(pk_index.at(s),3.0))));
        val_2 = val_2 + pow((1.0/((double) pow(pk_index.at(s), 3.0))), 2.0);
    }
    
    cs = val_1/val_2;
    cn = 0.0;
    n_avg = double(n)/10.0;
    
    for (int s = (s_pk+2); s <= (s_pk+int(n_avg)); s++) {
        cn = cn + pow(pk_freq.at(s),2.0);
    } 
    
    cn = sqrt(cn/n_avg);
    cn2 = cn*cn;
    
    for (int i = 0; i < n; i++) {

        val_3 = 0.0;
        val_4 = 0.0;

        for (int k = 0; k < n-1; k++) {
            val_3 = pow((cs/pow(double(k+1),3.0)),2.0);
            val_4 = val_4 + (val_3/(val_3+cn2) * b.at(k) * sin(double(k+1) * zpi.at(i)));
        }
        z_smooth = double(val_4 + z.at(0) + za.at(i));
        zs.push_back(z_smooth);
    }

    return zs;    
}
vector<int> Analysis::category(vector<double> v, vector<double> vals) {
    // determines if a value can be categorized as ultra-high, high, low, ultra-low, or average
    // uses 0.674 of the std deviation to partition the data into thirds
    // the average group has 50% of the data within its partition.
    // the ultra group is split for the top 10% and bottom 10%
    // ultra-low = -2, low = -1, average = 0, high = 1, ultra-high = 2

    double stddev = sqrt(stddev2(v));
    double mean = average(v);
    double interval_50pct = 0.674490 * stddev;
    double interval_20pct = 1.281552 * stddev;
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
 * Net class
 */
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
            cout << "Made a Neuron!" << endl;
        }

        // Force the bias node's output to 1.0 (it was the last neuron pushed in this layer):
        m_layers.back().back().setOutputVal(1.0);
    }
}
void showVectorVals(string label, vector<double> &v)
{
    cout << label << " ";
    for (unsigned i = 0; i < v.size(); ++i) {
        cout << v[i] << " ";
    }

    cout << endl;
}

/*
 * Genetic Algorithm class
 */
/*
 * Genetic Algorithm class
 */
class GeneticAlgorithm {
public:
    
    // genetic algorithm functions
    void generate_optimization(int numAsteroids, int numAntennas, int mutate, int popSize, double maxCycle); 
    double calculateFitness(int antennaID, double powerLevel);    
    double optimizationError(double targetVal, double currentVal);     
    vector<_pop> create_population(int number_of_asteroids, int number_of_antennas, int population_size);
    vector<_pop> compete_population(vector<_pop> population);
    vector<_pop> select_population(vector<_pop> population);
    vector<_pop> evolve_population(vector<_pop> population, int mutate);
    vector<double> create_values(double min, double max, int values);    
    bitset<64> create_64b_gene(double v, double min, double max);
    double create_64b_value(bitset<64> g, double min, double max);
    vector<_gn> mutate_genome(vector<_gn> old_genome, int mutate);
    _gn2 cross_genomes(_gn2 &old_genome_pair);
    void output_genome(vector<_gn> genome);
    
private:
    
};

void GeneticAlgorithm::generate_optimization(int numAsteroids, int numAntennas, int mutate, int popSize, double maxCycle) {
    _pop pop;
    vector<_gn> genome;
    vector<_pop> population, population_competed, population_ranked;

    int gen = 0;
    int j = 0;
    int m = 0;
        
    population.clear();
    population = create_population(numAsteroids, numAntennas, popSize); 
    double score = 0.0;

    while (j < maxCycle) {
        
        gen++;
        
        
        m = round((double) mutate 
                * (1.0 - cos(3.141592654 
                * ((double) gen
                / (double) popSize))) 
                + 1.0);
        population_competed.clear();
        population_ranked.clear();
        population_competed = compete_population(population);
        population_ranked = select_population(population_competed);
        population = evolve_population(population_ranked, m);
        j++;
        double newScore = population.at(0).fitness / numAntennas;
        if (score < newScore) {
            score = newScore;
            cout << "Generation #: " << gen << " Fitness: " << score << endl;
            j=0;
        }
    }

    cout << "Mutations: " << m << "\tGeneration #: " << gen << "\tBest Fitness: " << score << endl;  

    pop = population.at(0);
    genome = pop.genome;
    int gene_groups = genome.size() / 2;
    _gn gene_1, gene_2;

    for (int i = 0; i < gene_groups; i++) {
        int gene_num = i*2;
        gene_1 = genome.at(gene_num+0);
        gene_2 = genome.at(gene_num+1);

        cout << "Antenna #: " << i << endl;
        cout << "targetID (1): " << round(create_64b_value(gene_1.gene, gene_1.min_value, gene_1.max_value)) << endl
            << "powerLevel (25.0%): " << create_64b_value(gene_2.gene, gene_2.min_value, gene_2.max_value) << endl << endl;

    }
}
double GeneticAlgorithm::calculateFitness(int targetID, double powerLevel) {
    
    double fitness = 1000.0 - max(
            optimizationError(1.0, targetID),
            optimizationError(25.0, powerLevel)
            ); 
    
    //cout << "fitness: " << fitness << endl;

    return fitness;
}
double GeneticAlgorithm::optimizationError(double targetVal, double currentVal) {
    double targError;
    
    targError = abs((targetVal-currentVal)/targetVal) * 100.0;
    return targError;

}
vector<_pop> GeneticAlgorithm::create_population(int number_of_asteroids, int number_of_antennas, int population_size) {
// generate a population using with parameters defined to the objective function / competition
    
    // for each radio telescope antenna in the array:
    //     1) primary asteroid target #
    //     2) emitter power level from 0 - 100% (0% is receiving mode)
    //     3) delta time for the current command state

    _gn gene1;
    _gn gene2;

    vector<_gn> genome;
    _pop pop;
    vector<_pop> population;    
    
    vector<double> targetID;
    double targetID_min = 0.0;
    double targetID_max = (double) number_of_asteroids;
   
    vector<double> powerLevel;
    double powerLevel_min = 0.0;
    double powerLevel_max = 100.0;

    for (int i = 0; i < population_size; i++) {
        genome.clear();

        targetID = create_values(targetID_min, targetID_max, number_of_antennas);
        powerLevel = create_values(powerLevel_min, powerLevel_max, number_of_antennas);
        
        for (int j = 0; j < number_of_antennas; j++) {
            gene1.gene_id = 0;
            gene1.gene = create_64b_gene(targetID.at(j), targetID_min, targetID_max);
            gene1.min_value = targetID_min;     
            gene1.max_value = targetID_max;
            genome.push_back(gene1);
            
            gene2.gene_id = 1;
            gene2.gene = create_64b_gene(powerLevel.at(j), powerLevel_min, powerLevel_max);
            gene2.min_value = powerLevel_min;       
            gene2.max_value = powerLevel_max;
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
vector<_pop> GeneticAlgorithm::compete_population(vector<_pop> population) {
// determine the fitness of each member of the population

    _gn gn;
    _pop pop;
    vector<_pop> new_population;
    double targetID, powerLevel;
    int num_gene_groups = population.at(0).genome.size() / 2;

    for (int i = 0; i < population.size(); i++) {
        pop = population.at(i);
        pop.fitness = 0.0;
        for (int j = 0; j < num_gene_groups; j++) {
            int gene_num = j*2;
            targetID = round(create_64b_value(
                    pop.genome.at(gene_num+0).gene,
                    pop.genome.at(gene_num+0).min_value,
                    pop.genome.at(gene_num+0).max_value
                ));
            powerLevel = create_64b_value(
                    pop.genome.at(gene_num+1).gene,
                    pop.genome.at(gene_num+1).min_value,
                    pop.genome.at(gene_num+1).max_value
                );

            pop.fitness = pop.fitness + calculateFitness(targetID, powerLevel);
        }

        new_population.push_back(pop);
    }    
    
    return new_population;

}
vector<_pop> GeneticAlgorithm::select_population(vector<_pop> population) {
// select the primary candidates for evolution
    
    _pop pop;
    vector<object_tuple> test_tuple;    
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
    for(vector<object_tuple>::iterator iter = test_tuple.begin(); iter != test_tuple.end(); iter++){
        ranked_popID.push_back(get<0>(*iter)); 
        ranked_fitness.push_back(get<1>(*iter));
    }     
    
    for (int j = 0; j < ranked_popID.size(); j++) {
        new_population.push_back(population.at(ranked_popID.at(j)));
    }
    
    return new_population;
}
vector<_pop> GeneticAlgorithm::evolve_population(vector<_pop> population, int mutate) {
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
vector<double> GeneticAlgorithm::create_values(double min, double max, int values) {
    vector<double> random_values;
    
    default_random_engine rand_engine;
    uniform_real_distribution<double> distribution(min, max);    
 
    for (int i = 0; i < values; i++) {
        random_values.push_back(distribution(rand_engine));
    }
    
    return random_values;
}
bitset<64> GeneticAlgorithm::create_64b_gene(double v, double min, double max) {    
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
double GeneticAlgorithm::create_64b_value(bitset<64> g, double min, double max) {
    double value;
    double v, v_max;
    bitset<64> b_64_max;
    
    b_64_max.set();
    v_max = (double) b_64_max.to_ullong();
    v = (double) g.to_ullong();
    value = ((v/v_max) * (max-min)) + min;
    return value;
}
vector<_gn> GeneticAlgorithm::mutate_genome(vector<_gn> old_genome, int mutate) {
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
_gn2 GeneticAlgorithm::cross_genomes(_gn2 &old_genome_pair) {

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
void GeneticAlgorithm::output_genome(vector<_gn> genome) {
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

/*
 * AsteroidDetector class
 */
class AsteroidDetector {
public:
    int trainingData(int width, int height, vector<int> imageData_1, vector<string> header_1, vector<double> wcs_1, vector<int> imageData_2, vector<string> header_2, vector<double> wcs_2, vector<int> imageData_3, vector<string> header_3, vector<double> wcs_3, vector<int> imageData_4, vector<string> header_4, vector<double> wcs_4, vector<string> detections);
    int testingData(string imageID, int width, int height, vector<int> imageData_1, vector<string> header_1, vector<double> wcs_1, vector<int> imageData_2, vector<string> header_2, vector<double> wcs_2, vector<int> imageData_3, vector<string> header_3, vector<double> wcs_3, vector<int> imageData_4, vector<string> header_4, vector<double> wcs_4);
    vector<string> getAnswer();
    
    void initDetections(vector<string> &detections);
    void initImages(int width, int height, vector<int> &image, vector<string> &header, vector<double> &wcs, bool isTraining);
    void initSet(bool isTraining);
    
    void trainNN_Detect(void);
    void trainNN_NEO(void);
    void execNN_Detect(void);
    void execNN_NEO(void);
    void configNN_Detect(void);
    void configNN_NEO(void);
    void initNN_Detect_inputs(bool isTraining);
    void initNN_NEO_inputs(bool isTraining);
    void initNN_Detect_targets(void);
    void initNN_NEO_targets(void);
    void cycleNN_Detect(void);
    void cycleNN_NEO(void);
    void generateNN_Detect(void);
    void generateNN_NEO(void);
    void generateResults(void);
    
    // utility functions
    vector<double> convertRADEC2XY(_image I, double RA, double DEC);
    vector<double> convertXY2RADDEC(_image I, double X, double Y);
    
    // initialize vectors
    vector<_detection> tD; // collection of detections in the training data
    vector<_image> tI; // collection of 4 images in the training data
    vector<_imageSet> tS; // collection of vector<_image> for training
    
    vector<_detection> xD; // collection of detections in the test data
    vector<_image> xI; // collection of 4 images in the test data
    vector<_imageSet> xS; // collection of vector<_image> for testing
    
    vector<string> the_big_answer; // a vector<string> set of detections
    vector<double> inputVals, targetVals, resultVals;
    vector<unsigned> topology_NEO, topology_DETECT;
    
    // initialize classes
    Analysis A;
    Net Net_NEO;
    Net Net_DETECT;

    // initialize booleans
    bool startTraining = false;
    bool startTesting = false;
    bool continueTraining = true;

    // initialize integers

    // initialize double
    double training_start_time = 0.0;
    double testing_start_time = 0.0;
    double get_answer_start_time = 0.0;
    double training_elapsed_time = 60.0 * 1000000.0; // 60 seconds
    double testing_elapsed_time = 0.0;
    
private:
    
};
void AsteroidDetector::initDetections(vector<string> &detections) { // used for training only
    // convert vector<string> detections into vector<_detection>
    _detection d;    
    string line;
    
    xD.clear();
    // convert each string to a _detection
    for (int i = 0; i < detections.size(); i++) {
        line = detections.at(i);
        istringstream ss(line);
        if (!(ss >> d.uniqueID >> d.frameNumber >> d.RA >> d.DEC >> d.X >> d.Y >> d.magnitude >> d.NEO)) {
            break;
        }
        ss.str("");
        tD.push_back(d);

    }
}
void AsteroidDetector::initImages(int width, int height, vector<int> &image, vector<string> &header, vector<double> &wcs, bool isTraining) { // used for training only
    _image i;
    
    if (isTraining) {
        i.width = width;
        i.height = height;
        i.imageData = image;
        i.header = header;
        i.wcs = wcs;

        tI.push_back(i);
    } else {
        i.width = width;
        i.height = height;
        i.imageData = image;
        i.header = header;
        i.wcs = wcs;

        xI.push_back(i);
    }
}
void AsteroidDetector::initSet(bool isTraining) { // used for training only
    _imageSet s;
    
    if (isTraining) { 
        s.imageID = tS.size() + 1;
        s.images = tI;
        s.detections = tD;

        tS.push_back(s);
    } else {
        s.imageID = xS.size() + 1;
        s.images = xI;
        s.detections = xD;
        
        xS.push_back(s);
    }
}
void AsteroidDetector::trainNN_Detect(void){ 
    initNN_Detect_inputs(true);
    initNN_Detect_targets();
    configNN_Detect();
    cycleNN_Detect();
}
void AsteroidDetector::trainNN_NEO(void){ 
    initNN_NEO_inputs(true);
    initNN_NEO_targets();
    configNN_NEO();
    cycleNN_NEO();
}
void AsteroidDetector::execNN_Detect(void){ 
    initNN_Detect_targets();
    generateNN_Detect();
}
void AsteroidDetector::execNN_NEO(void){ 
    initNN_NEO_targets();
    generateNN_NEO();
}
void AsteroidDetector::configNN_Detect(void){ 
}
void AsteroidDetector::configNN_NEO(void){ 
}
void AsteroidDetector::initNN_Detect_inputs(bool isTraining){ 
}
void AsteroidDetector::initNN_NEO_inputs(bool isTraining){ 
}
void AsteroidDetector::initNN_Detect_targets(void){ 
}
void AsteroidDetector::initNN_NEO_targets(void){ 
}
void AsteroidDetector::cycleNN_Detect(void){ 
}
void AsteroidDetector::cycleNN_NEO(void){ 
}
void AsteroidDetector::generateNN_Detect(void){ 
}
void AsteroidDetector::generateNN_NEO(void){ 
}
void AsteroidDetector::generateResults(void){ 
}

/*
 * Utility functions for AsteroidDetector
 */
vector<double> AsteroidDetector::convertRADEC2XY(_image FITS, double RA, double DEC){
    vector<double> XY;
    double dR = RA - FITS.wcs.at(2);
    double dD = DEC - FITS.wcs.at(3);
    double dY = (dR * FITS.wcs.at(6) - dD * FITS.wcs.at(4))/(FITS.wcs.at(6)*FITS.wcs.at(5) - FITS.wcs.at(4)*FITS.wcs.at(7));
    double dX = (dR - dY * FITS.wcs.at(5)) / FITS.wcs.at(4);
    XY.push_back(dX + FITS.wcs.at(0));
    XY.push_back(dY + FITS.wcs.at(1));
    return XY;
}
vector<double> AsteroidDetector::convertXY2RADDEC(_image FITS, double X, double Y){
    vector<double> RD;
    double dX = X - FITS.wcs.at(0);
    double dY = Y - FITS.wcs.at(1);
    RD.push_back(dX * FITS.wcs.at(4) + dY * FITS.wcs.at(5) + FITS.wcs.at(2));
    RD.push_back(dX * FITS.wcs.at(6) + dY * FITS.wcs.at(7) + FITS.wcs.at(3));
    return RD;
}
/*
 * Firstly, your trainingData method will be called with FITS images, headers 
 * and known detected objects. Your method will be called 100 times, one 
 * time for every set in the training data. You can use this method to train 
 * your algorithm on the provided data if you want to. If your trainingData 
 * method returns the value 1, no more training data will be passed to your 
 * algorithm and the testing phase will begin.
 */
int AsteroidDetector::trainingData(int width, int height, 
        vector<int> imageData_1,vector<string> header_1, vector<double> wcs_1, 
        vector<int> imageData_2, vector<string> header_2, vector<double> wcs_2, 
        vector<int> imageData_3, vector<string> header_3, vector<double> wcs_3, 
        vector<int> imageData_4, vector<string> header_4, vector<double> wcs_4, 
        vector<string> detections) {

    /* 
     * train a NN to classify if a detected object is a NEO {inputs: velocity vector, magnitude.  output:  NEO (yes/no)}
     * train a NN to determine if an image has a detection in it.
     * run statistics on each detection to build a magnitude model (magnitude: continuous)
     * run statistics on each detection to build a detection model (detection: yes / no)
     * convert x,y to RA and DEC
     */
    
    // start the training clock
    if (startTraining == false) {
        training_start_time = clock(); 
        cout << "trainingData start time: " << (training_start_time / 1000000.0) << " sec" << endl;         
        startTraining = true;
    }   

    int total_bits_1 = A.max_int(imageData_1);
    
    // initialize the data vectors with image data and detections information
    initDetections(detections); // convert detections string to vector<_detections> -> tD
    initImages(width, height, imageData_1, header_1, wcs_1, true); // convert image vectors to vector<_image> -> tI
    initImages(width, height, imageData_2, header_2, wcs_2, true); // convert image vectors to vector<_image> -> tI
    initImages(width, height, imageData_3, header_3, wcs_3, true); // convert image vectors to vector<_image> -> tI
    initImages(width, height, imageData_4, header_4, wcs_4, true); // convert image vectors to vector<_image> -> tI
    initSet(true); // convert the set of images to vector<_imageSet> -> tS
    
    // convert the images to 32px x 32px tiles
    
    // perform training on the NN
    trainNN_Detect();
    trainNN_NEO();
    
    // check to see if the allocated training time has been completed
    if (clock() < (training_start_time + training_elapsed_time)) {
        continueTraining = true;
    } else {
        continueTraining = false;
        cout << "trainingData end time: " << (clock() / 1000000.0) << " sec" << endl;
        cout << "==========================================================================================================" << endl;                
    }
    
    // return a value of 1 if the training data is complete -> proceed to testing.
    if (continueTraining) { return 0; } else { return 1; }
}
/*
 * Secondly, your testingData method will be called 20 times with different 
 * image data than provided during training. The testingData method can return 
 * any value, it does not matter what you return.
 */
int AsteroidDetector::testingData(string imageID, int width, int height, 
        vector<int> imageData_1, vector<string> header_1, vector<double> wcs_1, 
        vector<int> imageData_2, vector<string> header_2, vector<double> wcs_2, 
        vector<int> imageData_3, vector<string> header_3, vector<double> wcs_3, 
        vector<int> imageData_4, vector<string> header_4, vector<double> wcs_4) {
    
    // start the testing clock
    if (startTesting == false) {
        testing_start_time = clock(); 
        cout << "testingData start time: " << (testing_start_time / 1000000.0) << " sec" << endl;         
        startTesting = true;
    }  
    
    // initialize the data vectors with image data -> xI
    initImages(width, height, imageData_1, header_1, wcs_1, false); // convert image vectors to vector<_image> -> xI
    initImages(width, height, imageData_2, header_2, wcs_2, false); // convert image vectors to vector<_image> -> xI
    initImages(width, height, imageData_3, header_3, wcs_3, false); // convert image vectors to vector<_image> -> xI
    initImages(width, height, imageData_4, header_4, wcs_4, false); // convert image vectors to vector<_image> -> xI
  
    // convert the images to 32px x 32px tiles
    
    
    // evaluate the images using the trained NN
    execNN_Detect();
    execNN_NEO();
    
    // transfer the data into a vector<_imageSet>
    initSet(false); // generate a set of images -> xS
    
    // generate a list of detections -> xD
    
    // stop the testing clock
    testing_elapsed_time = (clock() - testing_start_time) / 1000000.0;
    cout << "testingData end time: " << (clock() / 1000000.0) << " sec" << endl;
    cout << "testingData elapsed time: " << testing_elapsed_time << " sec" << endl;
    cout << "==========================================================================================================" << endl;                
    
    return _AsteroidDetector_testing_; 
}
/*
 * Finally, your getAnswer method will be called. This method should return 
 * a list of all the objects that your algorithm detected for each set provided 
 * to your algorithm through the testingData method. Your may not return 
 * more than 100000 detections. Each element in your return should contain 
 * the following information in space delimited format:

    1. imageID - the imageID associated with the image where the object was detected
    2. RA_1 - right ascension of object in decimal hours in frame 1.
    3. DEC_1 - declination in decimal degrees in frame 1.
    4. RA_2 - right ascension of object in decimal hours in frame 2.
    5. DEC_2 - declination in decimal degrees in frame 2.
    6. RA_3 - right ascension of object in decimal hours in frame 3.
    7. DEC_3 - declination in decimal degrees in frame 3.
    8. RA_4 - right ascension of object in decimal hours in frame 4.
    9. DEC_4 - declination in decimal degrees in frame 4.
    10. NEO - this value should be 1 if your algorithm think the object is 
              a Near Earth Object (NEO), 0 otherwise.
 */

vector<string> AsteroidDetector::getAnswer() {
    _detection d;
    _imageSet t;
    
    ostringstream targetLine;
    int detect = 0;
    
    // iterate thru the detections found in the testing to generate the answer
    for (int i = 0; i < xS.size(); i++) {
        t = xS.at(i);
        
        targetLine.str("");
        targetLine << t.imageID << " ";
        for (int j = 0; j < t.detections.size(); j++) {
            d.RA = t.detections.at(j).RA;
            d.DEC = t.detections.at(j).DEC;
            d.NEO = t.detections.at(j).NEO;
            
            // load the output into the_big_answer as a vector<string>
            targetLine << d.RA << " " << d.DEC << " "; 
            detect++;

        }
        targetLine << d.NEO;
        
        if (detect < 100000) { // return up to 100000 identified objects
            the_big_answer.push_back(targetLine.str());
        }
    }
    
    // output statistics on the results
    cout << "Found " << the_big_answer.size() << " Targets" << endl;
    cout << "==========================================================================================================" << endl;                
  
    return the_big_answer;
}
string clearSpaces(string str) {
    /*
     *     
     */
     
    std::string::iterator new_end = std::unique(str.begin(), str.end(), BothAreSpaces);
    str.erase(new_end, str.end());
    return str;
}

int main(int argc, char** argv) {
    
    string readLine;
    vector<string> dataNames, knownAsteroids;
    vector<_top> obj;
    _top o;
    
    Astrodynamics AD;
    
    string dataPath = "/Users/scott/Arctria/ARC-32/AsteroidDetector/data/supplemental/";
    string imagePath = "/Users/scott/Arctria/ARC-32/AsteroidDetector/data/traindata/";
    string trainingList = "/Users/scott/Arctria/ARC-32/AsteroidDetector/data/supplemental/traindata.txt";
    
    string Amors = "/users/scott/Arctria/ARC-32/AsteroidDetector/data/supplemental/Amors.txt";
    string Apollos = "/users/scott/Arctria/ARC-32/AsteroidDetector/data/supplemental/Apollos.txt";
    string Atens = "/users/scott/Arctria/ARC-32/AsteroidDetector/data/supplemental/Atens.txt";
    
    ifstream trainingFile (trainingList);
    ifstream amoFile (Amors);
    ifstream apolloFile (Apollos);
    ifstream atenFile (Atens);

    while (getline(trainingFile, readLine)) {
        dataNames.push_back(readLine);
    }
    
    cout << "==========================================================================================================" << endl;                
  
    cout << "Training Files: " << dataNames.size() << endl;
    /*for (int i = 0; i < dataNames.size(); i++) {
        cout << i << ": " << dataNames.at(i) << endl;
    }*/
    cout << "==========================================================================================================" << endl;                
    cout << "Asteroids: " << endl;
    
    int n = 0;
    string year, name, details;
    
    double q, Q;
    while (getline(amoFile, readLine)) {

    
/*
 * Designation (and name)     Prov. Des.     q       Q     EMoid     H     Epoch      M    Peri. Node   Incl.  e      a      Opps.     Ref.      Designation (and name)                  Discovery date, site and discoverer(s)
 */ 
        n++;
        if (n > 2) {
            istringstream iss(clearSpaces(readLine));
            iss >> o.YODA >> o.NAME >> o.q >> o.Q >> o.EMOID >> o.AVM >> o.EPOCH >> o.TMA >> o.APG >> o.RAAN >> o.INC >> o.ECC >> o.SMA >> details;
            obj.push_back(o);
            knownAsteroids.push_back(o.YODA + " " + o.NAME);
            /*cout << knownAsteroids.size() << ": " << knownAsteroids.back() 
                    << "\t\ta: " << o.SMA 
                    << "\te: " << o.ECC
                    << "\ti: " << o.INC
                    << "\tw: " << o.APG
                    << "\tAN: " << o.RAAN
                    << "\tMA: " << o.TMA
                    << endl;*/
        }
    }
    n = 0;
    while (getline(apolloFile, readLine)) {
        n++;
        if (n > 2) {
            istringstream iss(clearSpaces(readLine));
            iss >> o.YODA >> o.NAME >> o.q >> o.Q >> o.EMOID >> o.AVM >> o.EPOCH >> o.TMA >> o.APG >> o.RAAN >> o.INC >> o.ECC >> o.SMA >> details;
            obj.push_back(o);
            knownAsteroids.push_back(o.YODA + " " + o.NAME);
            /*cout << knownAsteroids.size() << ": " << knownAsteroids.back()
                    << "\t\ta: " << o.SMA 
                    << "\te: " << o.ECC
                    << "\ti: " << o.INC
                    << "\tw: " << o.APG
                    << "\tAN: " << o.RAAN
                    << "\tMA: " << o.TMA
                    << endl;*/
        }
    }
    
    n = 0;
    
    while (getline(atenFile, readLine)) {
        n++;
        if (n > 2) {
            istringstream iss(clearSpaces(readLine));
            iss >> o.YODA >> o.NAME >> o.q >> o.Q >> o.EMOID >> o.AVM >> o.EPOCH >> o.TMA >> o.APG >> o.RAAN >> o.INC >> o.ECC >> o.SMA >> details;
            obj.push_back(o);
            knownAsteroids.push_back(o.YODA + " " + o.NAME);
            /*cout << knownAsteroids.size() << ": " << knownAsteroids.back() 
                    << "\t\ta: " << o.SMA 
                    << "\te: " << o.ECC
                    << "\ti: " << o.INC
                    << "\tw: " << o.APG
                    << "\tAN: " << o.RAAN
                    << "\tMA: " << o.TMA
                    << endl;*/
        }
    }

    cout << "Known Asteroids: " << knownAsteroids.size() << endl;
    
    cout << "Julian Date: " << AD.Epoch_to_JD(20141225,2119.0);
    cout << "nu: " << AD.trueAnomaly(245.0, 0.95) << endl;
    
    return 0;
}
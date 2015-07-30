/* Analysis v10 [Competition Sensitive]
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
int main(int argc, char** argv) {
    Analysis A;
    vector<double> dataset;
    vector<double> vals;
    vector<int> cats;
    
    for (int i = 0; i < 100; i++) {
        dataset.push_back((double) i );
    }
    
    for (double val = 0; val < 100.0; val=val+0.1) {
        vals.push_back(val);
    }
    cats = A.category(dataset,vals);
    
    for (int i = 0; i < vals.size(); i++) {
        cout << "Category for " << vals.at(i) << " = " << cats.at(i) << endl;
    }

    return 0;
}
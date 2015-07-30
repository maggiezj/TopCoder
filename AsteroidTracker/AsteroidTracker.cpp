/* AsteroidTracker v10 [Competition Sensitive]
 * NASA Asteroid Grand Challenge | NASA | Planetary Resources | TopCoder
 * RigelFive - Hudson Ohio USA
 * 25 July 2014 - 9 August 2014
 */

#define _If_Nacho_Libre_Can_Wake_Up_At_Five_A_M_To_Make_The_Soup_ 0;
#define _It_Is_The_Besssssssst_ 1;

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

bool sort_object_score (const object_tuple &lhs, const object_tuple &rhs){
  return get<1>(lhs) > get<1>(rhs);  // sorts for largest numbers at beginning of list
}

struct _traj { // trajectory of an asteroid
    double t;
    double x;
    double y;
    double z;
};
struct _signal {
    int antennaID;
    double arrivalTime; // time of arrival for signal from asteroid
    double signalPower;
};
struct _i { // asteroid information
     int asteroidID;
     double scienceScoreMultiplier;
     double reflectivityMultiplier; 
     double imageInformation; 
     double trajectoryInformation; 
     double timeAppear;

     vector<_traj> trajectory;
     vector<_signal> signalStrength;    
};
struct _j { // antenna information
    // static parameters
    int antennaID;
    double Xpos;
    double Ypos;
    
    // dynamic parameters v. time
    double power; // dynamic power setting
    double Xpoint; // dynamic X pointing coordinate to target asteroid
    double Ypoint; // dynamic Y pointing coordinate
    double Zpoint; // dynamic Z pointing coordinate
    int targetAsteroidID; // dynamic condition as to which asteroid the antenna is pointing at
    double antennaMinDistance;
    double antennaMaxDistance;
    double nextCommandAvailabilityTime; 
    bool relocating;
};
struct _gain {
    double minDistanceGain;
    double maxDistanceGain;
};
struct _link {
    int antennaID;
    int asteroidID;
    double reflectivity;
    double distance;
    double travelTime;
    double figureOfMerit;
    bool visibleAtBeamStart;
    bool visibleAtBeamEnd;
};
struct _command {
    int commandTime;
    bool pointDelta;
    bool powerDelta;
    int antennaID;
    int targetAsteroidID;
    double powerLevel;
    string commandNote;
};

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
/*
 * AsteroidTracker Class
 */
class AsteroidTracker {
public:
    // primary functions and methods
    int initialize(vector<double> antennaPositions, double peakGain, vector<double> minDistanceGain, vector<double> maxDistanceGain);
    int asteroidAppearance(int asteroidIndex, double scienceScoreMultiplier, double reflectivityMultiplier, double initialImageInformation, double initialTrajectoryInformation, vector<double> trajectory);
    string nextCommand(double currentTime);
    
    // functions to calculate the fitness of the objective function
    double distanceBetweenAntennas(int antennaID_1, int antennaID_2); //verified
    _traj getAsteroidPosition(vector<_traj> asteroidTrajectory, double viewTime); //    
    double distanceToAsteroid(_traj asteroidPosition); //
    double signalTravelTime(_traj asteroidPosition);
    vector<int> getSubarray(int asteroidID); // verified       
    bool isAsteroidVisible(vector<_traj> asteroidTrajectory, double viewTime);
    bool isAntennaRelocating(int antennaID, double viewTime); // verified    
    bool checkAntennaBeamIntersection(int antennaID_T, vector<_traj> targetAsteroidTraj, double viewTime); // verified - but no events???    
    bool checkTargetProximity(int asteroidID_A, int asteroidID_B, double viewTime); // verified
    bool checkLinkProximity(vector<_link> rankedLinkList, double viewTime);
    vector<_link> generateLinks(double currentTime);
    vector<_link> rankLinks(vector<_link> firstLinkList);
    vector<int> generateTargetLinks(vector<_link> rankedList);
    double angleBetweenAsteroids(_traj asteroidPosition_a, _traj asteroidPosition_b);
    double calcSlewAngle(vector<_traj> asteroidTraj_a, vector<_traj> asteroidTraj_b, double viewTime); // verified    
    vector<double> generateTimesToRelocate(vector<int> targetIDList, double viewTime);
    void relocateAntennas(vector<int> targetID, double timeToSlew, double viewTime);
    double activateAntennas(double viewTime);
    void monitorAsteroids(double beamTime, double viewTime);
    void retireAntennas(void);

// static constants
    double SPEED_OF_LIGHT = 299792458.0; // speed of light in meters per second
    double SIMULATION_TIME = 604800.0; // 7 days
    double T_MIN = 0.1;
    double CRITICAL_TRACKING_ANGLE = 0.000290888208665721596153948461414768785573811981423620907407; // 1 arc minute
    double ANTENNA_SAFE_RADIUS = 11.0; // 11 meters
    double Q_LOST = 124648.8515328064383958974924385634806736622104388179847093028; // half-life: 24h
    double Q_TRAJECTORY = 6000.0;
    double Q_IMAGE = 1000000.0;
    double IMAGE_SCORE_MULTIPLIER = 30.0;
    double TRAJECTORY_SCORE_MULTIPLIER = 0.000099206349206349206349206349206349206349206349206349206349;
    double MAX_TRANSMITTING_POWER = 40000; // 40 kW
    double BACKGROUND_NOISE = 1.0e-30;
    double SHIELDING_MULTIPLIER = 1.0e-27;
    double RELOCATE_SPEED = 0.005235987755982988730771072305465838140328615665625176333333; // [rad/s] 10 min to relocate 180°
    double RELOCATE_POWER = 5000; // 5 kW
    
    double LATITUDE = 28.524812;
    double LONGITUDE = 0.0;
    double MEAN_RADIUS_EARTH = 6367.4447; // average radius of earth in km (Wolfram Alpha)
    
    double PEAK_GAIN = 0;
    double ANTENNA_MIN_DISTANCE = 0;
    double ANTENNA_MAX_DISTANCE = 0;

    // primary data vectors
    vector<_i> asteroids;
    vector<_j> antennas;
    vector<_gain> antennaGain; 
    vector<_command> commandList;
    
    // Simulator parameters
    int numberOfAntennas = 0;
    int numberOfAsteroids = 0; 
    double energySpent = 0.0;
    double viewTime = 0.0;
    double nextCommandTime = 0.0;
    
    bool commandDeckLoaded = false;
    bool antennasRetired = false;
    int commandNumber = 0;
    
private:
    
};
// Functions to determine the fitness of the objective function
double AsteroidTracker::distanceBetweenAntennas(int antennaID_1, int antennaID_2) {
    double dx, dy, dist;
    
    dx = antennas.at(antennaID_1).Xpos - antennas.at(antennaID_2).Xpos;
    dy = antennas.at(antennaID_1).Ypos - antennas.at(antennaID_2).Ypos;
    dist = sqrt((dx*dx)+(dy*dy));
    
    return dist;
}
_traj AsteroidTracker::getAsteroidPosition(vector<_traj> asteroidTrajectory, double viewTime) {
    _traj pos;

    double t1, x1, y1, z1, t2, x2, y2, z2;
    
    for (int i = 1; i < asteroidTrajectory.size(); i++) {
        t1 = asteroidTrajectory.at(i-1).t;
        t2 = asteroidTrajectory.at(i).t;
        
        if ((viewTime < t2) && (viewTime > t1)) {
            x1 = asteroidTrajectory.at(i-1).x;
            y1 = asteroidTrajectory.at(i-1).y;
            z1 = asteroidTrajectory.at(i-1).z;
            
            x2 = asteroidTrajectory.at(i).x;
            y2 = asteroidTrajectory.at(i).y;
            z2 = asteroidTrajectory.at(i).z;
            
            // the trajectory data points between the time point viewTime has been found
            pos.t = viewTime;
            pos.x = ((x2-x1)/(t2-t1))*viewTime+(x1-(((x2-x1)/(t2-t1))*t1));
            pos.y = ((y2-y1)/(t2-t1))*viewTime+(y1-(((y2-y1)/(t2-t1))*t1));
            pos.z = ((z2-z1)/(t2-t1))*viewTime+(z1-(((z2-z1)/(t2-t1))*t1));           
        } else if (viewTime == t1) {
            x1 = asteroidTrajectory.at(i-1).x;
            y1 = asteroidTrajectory.at(i-1).y;
            z1 = asteroidTrajectory.at(i-1).z;
            pos.x = x1;
            pos.y = y1;
            pos.z = z1;
            pos.t = viewTime;
        } else if (viewTime == t2) {
            x2 = asteroidTrajectory.at(i).x;
            y2 = asteroidTrajectory.at(i).y;
            z2 = asteroidTrajectory.at(i).z;
            pos.x = x2;
            pos.y = y2;
            pos.z = z2;
            pos.t = viewTime;
        }
    }
    
    return pos;
}
double AsteroidTracker::distanceToAsteroid(_traj asteroidPosition) { 
    _traj pos = asteroidPosition;
    return sqrt((pos.x*pos.x)+(pos.y*pos.y)+(pos.z*pos.z));
}
double AsteroidTracker::signalTravelTime(_traj asteroidPosition) {
    double travelTime, travelDist;
    
    travelDist = distanceToAsteroid(asteroidPosition);
    travelTime = (2.0 * travelDist)/SPEED_OF_LIGHT;
    
    return travelTime;
}
vector<int> AsteroidTracker::getSubarray(int asteroidID) {
    /*
     * For each asteroid, generate the list of antennas tracking it
     */ 
    
    _i asteroid;
    _j antenna;
    
    vector<int> tracked_by_antenna;

    // load the asteroid of interest to see if there are antennas pointed at it.
    asteroid = asteroids.at(asteroidID);

    for (int j = 0; j < antennas.size(); j++) {
        antenna = antennas.at(j);

        // if the asteroid is being tracked by an antenna, return the antennaID in the array.
        if (antenna.antennaID == asteroidID) {
            tracked_by_antenna.push_back(antenna.antennaID);
        }
    }

    return tracked_by_antenna;
}
bool AsteroidTracker::isAsteroidVisible(vector<_traj> asteroidTrajectory, double viewTime) {
    bool visible;
    _traj traj1, traj2;
    double z1, z2, t1, t2, t_min, t_max;
    
    t_min = asteroidTrajectory.at(0).t;
    t_max = asteroidTrajectory.back().t;
    
    visible = false;  // default to false
    
    for (int i = 1; i < asteroidTrajectory.size(); i++) {
        traj1 = asteroidTrajectory.at(i-1);
        traj2 = asteroidTrajectory.at(i);
        
        t1 = traj1.t;
        t2 = traj2.t;
        z1 = traj1.z;
        z2 = traj2.z;
        
        if ((viewTime <= t2) && (viewTime >= t1)) {
            if ((z1 > 0) && (z2 > 0)) {
                visible = true;
            } else if ((z1 < 0) && (z2 < 0)) {
                visible = false;
            } else { // where the data shows that the asteroid has crossed the horizon (assumption is over the past three minutes)
                if (abs(t1-t2) < 180.1) {
                    // if one of the z coordinates is negative, then check to see if there is a new day
                    // if the data is less than 180 seconds apart, then calculate the precise moment that the 
                    // asteroid appears on the horizon...

                    double timeOnHorizon =  -(z1 - (((z2-z1)/(t2-t1)) * t1))/((z2-z1)/(t2-t1));
                    if (z2 > 0) { // ascending
                        if (viewTime > timeOnHorizon) {
                            visible = true;
                        }
                    } else if (z2 < 0) { // descending
                        if (viewTime < timeOnHorizon) {
                            visible = true;
                        }
                    }
                } else {
                    visible = false;
                }
            }
        }
    }

    return visible;
}
bool AsteroidTracker::isAntennaRelocating(int antennaID, double viewTime) {
    bool relocating;
    
    _j antenna;
    
    antenna = antennas.at(antennaID);
    
    if (viewTime <= antenna.nextCommandAvailabilityTime)  { 
        // just read the current state if the command has been set for an amount of time.
        
        relocating = antenna.relocating;
    } else { // update the current state to non-relocating with availability for the next command to now.
        antenna.relocating = false;
        antenna.nextCommandAvailabilityTime = viewTime;
        antennas.at(antennaID) = antenna;
    }    
    
    return relocating;
}
bool AsteroidTracker::checkAntennaBeamIntersection(int antennaID_T, vector<_traj> targetAsteroidTraj, double viewTime) {
    /*
     * Check if the beam from one antenna intersects any other antenna
     */

    _j antenna_R; // receiving antenna
    _j antenna_T; // transmitting antenna
    _traj asteroidPos;
    bool beam_too_close = false;
    
    double beamAngle, beamDistance;
    double a1b1, a2b2, a3b3, a2, b2, dX_a, dY_a, dX_b, dY_b, dZ_b;
    
    double PI = 3.1415926535897932384626433832795028841971693993751058;
        
    antenna_T = antennas.at(antennaID_T);

    for (int antennaID_r = 0; antennaID_r < antennas.size(); antennaID_r++) {
        antenna_R = antennas.at(antennaID_r);
        if (antennaID_r == antennaID_T) {
            continue;
        } else {
            asteroidPos = getAsteroidPosition(targetAsteroidTraj, viewTime);

            dX_a = antenna_T.Xpos - antenna_R.Xpos;
            dY_a = antenna_T.Ypos - antenna_R.Ypos;
            dX_b = asteroidPos.x;
            dY_b = asteroidPos.y;
            dZ_b = asteroidPos.z;

            a1b1 = dX_a * dX_b;
            a2b2 = dY_a * dY_b;
            a3b3 = 0;
            a2 = sqrt((dX_a*dX_a)+(dY_a*dY_a));
            b2 = sqrt((dX_b*dX_b)+(dY_b*dY_a)+(dZ_b*dZ_b));
            beamAngle = acos((a1b1 + a2b2 + a3b3) / (a2*b2)); 

            beamDistance = a2*sin(beamAngle);
            if (beamDistance < ANTENNA_SAFE_RADIUS) {
                beam_too_close = true;
            }
        }
    }
    
    return beam_too_close;
}
bool AsteroidTracker::checkTargetProximity(int asteroidID_A, int asteroidID_B, double viewTime) {
    /*
     * Check the target proximity
     */
    _traj pos_a, pos_b;
    _i asteroid_a, asteroid_b;
    
    bool asteroids_too_close = false;
    
    asteroid_a = asteroids.at(asteroidID_A);
    asteroid_b = asteroids.at(asteroidID_B);

    if (isAsteroidVisible(asteroid_a.trajectory, viewTime)
        && isAsteroidVisible(asteroid_b.trajectory, viewTime)
            ) {
        pos_a = getAsteroidPosition(asteroid_a.trajectory, viewTime);
        pos_b = getAsteroidPosition(asteroid_b.trajectory, viewTime);
        double angle = angleBetweenAsteroids(pos_a, pos_b);
    
        if (angle < CRITICAL_TRACKING_ANGLE) {
            asteroids_too_close = true;
        }
    }

    return asteroids_too_close;
}
bool AsteroidTracker::checkLinkProximity(vector<_link> rankedLinkList, double viewTime) {    
    vector<int> asteroidList;
    _link link;
    bool valueInList = false;
    bool linksCompatible = true;

        
    cout << "Size of rankedLinkedList: " << rankedLinkList.size() << endl;
    if (rankedLinkList.size() > 0) {
        
    
        asteroidList.push_back(rankedLinkList.at(0).asteroidID);
    
        for (int i = 0; i < rankedLinkList.size(); i++) {
            link = rankedLinkList.at(i);
            valueInList = false;
            for (int j = 0; j < asteroidList.size(); j++) {
                if (asteroidList.at(j) == link.asteroidID) {
                valueInList = true;
                }
            }
            if (!valueInList) {
            asteroidList.push_back(link.asteroidID);
            }
        }

    
        if (asteroidList.size() > 1) {
            for (int i = 0; i < asteroidList.size(); i++) {
                for (int j = i+1; j < asteroidList.size(); j++) {
                    if (checkTargetProximity(asteroidList.at(i),asteroidList.at(j), viewTime)) {
                        linksCompatible = false;
                    }
                }
            }
        }
    }
    
    return linksCompatible;
}
vector<_link> AsteroidTracker::generateLinks(double currentTime) {

    _i asteroid;
    _link link;
    vector<_link> links;
    _traj pos;
    
    links.clear();
    
    // generate a list of all possible links
    for (int antennaID = 0; antennaID < numberOfAntennas; antennaID++) {
        for (int asteroidID = 0; asteroidID < numberOfAsteroids; asteroidID++) {
            asteroid = asteroids.at(asteroidID);
            link.antennaID = antennaID;
            link.asteroidID = asteroidID;
            link.reflectivity = asteroid.reflectivityMultiplier;
            pos = getAsteroidPosition(asteroid.trajectory, currentTime);
            link.travelTime = signalTravelTime(pos);
            link.distance = distanceToAsteroid(pos);
            link.visibleAtBeamStart = isAsteroidVisible(asteroid.trajectory, currentTime);
            link.visibleAtBeamEnd = isAsteroidVisible(asteroid.trajectory, currentTime+link.travelTime);
            
            if ((link.visibleAtBeamStart)  // asteroid is visible at the beam start
                    && (link.visibleAtBeamEnd) // asteroid is visible at the beam end
                    && (!antennas.at(antennaID).relocating) // antenna is not relocating
                    && (!checkAntennaBeamIntersection(antennaID, asteroid.trajectory, currentTime)) // antenna is not interfering now
                    && (!checkAntennaBeamIntersection(antennaID, asteroid.trajectory, currentTime+link.travelTime)) // antenna is not interfering later
                    ) //
                        {   link.figureOfMerit = 
                            (numberOfAntennas 
                            * MAX_TRANSMITTING_POWER
                            * link.reflectivity 
                            * link.travelTime)
                            / link.distance;
            } else {
                link.figureOfMerit = 0.0;
            }
            links.push_back(link);   
        }
    }
    return links;
}
vector<_link> AsteroidTracker::rankLinks(vector<_link> firstLinkList) {
    vector<_link> rankedLinkList;
    vector<int> rankedID;
    _link link;

    vector<object_tuple> test_tuple;    
    
    for (int i = 0; i < firstLinkList.size(); i++) {
            link = firstLinkList.at(i);
            test_tuple.push_back(make_tuple(i, link.figureOfMerit));    
    }
    
    // sort the firstLinkList to find the highest values for figureOfMerit
    sort(test_tuple.begin(),test_tuple.end(),sort_object_score);    

        // transfer the ranked population information
    for(vector<object_tuple>::iterator iter = test_tuple.begin(); iter != test_tuple.end(); iter++){
        rankedID.push_back(get<0>(*iter));
    }  
    
    for (int i = 0; i < rankedID.size(); i++) {
        link = firstLinkList.at(rankedID.at(i));
        if (link.figureOfMerit > 0.0) {
            rankedLinkList.push_back(link);
        }

    }
    
    for (int i = 0; i < rankedLinkList.size(); i++) {
        link = rankedLinkList.at(i);
        /*cout << "antennaID: " << link.antennaID 
                << "\tasteroidID: " << link.asteroidID 
                << "\tFOM: " << link.figureOfMerit << endl;*/

    }

    
    return rankedLinkList;   
}
vector<int> AsteroidTracker::generateTargetLinks(vector<_link> rankedList) {
    vector<int> targetID;
    
    _link link;

    
    // set the default condition to turn off the antenna if there is no viable target.
    for (int antennaID = 0; antennaID < numberOfAntennas; antennaID++) {
        targetID.push_back(-1);
    }

    // displace the off condition with target conditions to best FOM targets.
    for (int antennaID = 0; antennaID < numberOfAntennas; antennaID++) {
        double maxFOM = 0.0;
        double maxFOM_index = 0;
        
        for (int i = 0; i < rankedList.size(); i++) {
            link = rankedList.at(i);
            if (link.antennaID == antennaID) {
                if (link.figureOfMerit > maxFOM) {
                    maxFOM = link.figureOfMerit;
                    maxFOM_index = i;
                }
            }
        }
        if (maxFOM > 1e-99) {
            targetID.at(antennaID) = rankedList.at(maxFOM_index).asteroidID;
            /*cout << "antennaID: " << antennaID << "\ttargetAsteroidID: " 
            << targetID.at(antennaID) << "\tfigureOfMerit: " 
            << rankedList.at(maxFOM_index).figureOfMerit
            << endl;*/
        } else {
            targetID.at(antennaID) = -1;
            /*cout << "antennaID: " << antennaID << "\ttargetAsteroidID: " 
            << targetID.at(antennaID) << "\tfigureOfMerit: 0" 
            << endl;*/
        }
 
    }
    
    return targetID;

}
double AsteroidTracker::angleBetweenAsteroids(_traj asteroidPosition_a, _traj asteroidPosition_b) {
    double angle;
    _traj pos_a, pos_b;
    pos_a = asteroidPosition_a;
    pos_b = asteroidPosition_b;
    
    double a1b1, a2b2, a3b3, a2, b2;
    
    a1b1 = pos_a.x * pos_b.x;
    a2b2 = pos_a.y * pos_b.y;
    a3b3 = pos_a.z * pos_b.z;

    a2 = sqrt((pos_a.x*pos_a.x)+(pos_a.y*pos_a.y)+(pos_a.z*pos_a.z));
    b2 = sqrt((pos_b.x*pos_b.x)+(pos_b.y*pos_b.y)+(pos_b.z*pos_b.z));
    
    /*cout << "a1b1: " << a1b1 << " a2b2: " 
            << a2b2 << " a3b3: " << a3b3 << endl
            << "a2: " << a2 << " b2: " << b2 << endl
            << " xA: " << pos_a.x << " yA: " << pos_a.y << " zA: " << pos_a.z << endl
            << " xB: " << pos_b.x << " yA: " << pos_b.y << " zA: " << pos_b.z << endl
            << endl;*/
   
    angle = acos((a1b1 + a2b2 + a3b3) / (a2*b2));
    
    return angle;
}
double AsteroidTracker::calcSlewAngle(vector<_traj> asteroidTraj_a, vector<_traj> asteroidTraj_b, double viewTime) {
    
    double PI = 3.1415926535897932384626433832795028841971693993751058;
    double slewAngle, slewTime;
    
    double xA, yA, zA, xB, yB, zB;
    
    _traj position_A;
    _traj position_B;
    
    position_A = getAsteroidPosition(asteroidTraj_a, viewTime);
    position_B = getAsteroidPosition(asteroidTraj_b, viewTime);
    
    xA = position_A.x;
    yA = position_A.y;
    zA = position_A.z;
    xB = position_B.x;
    yB = position_B.y;
    zB = position_B.z;
    
    //cout << "xA: " << xA << "\tyA: " << yA << "\tzA: " << zA;
    //cout << "\txB: " << xB << "\tyB: " << yB << "\tzB: " << zB << endl;
    
    double a1b1, a2b2, a3b3, a2, b2;

    a1b1 = xA * xB;
    a2b2 = yA * yB;
    a3b3 = zA * zB;
    a2 = sqrt((xA*xA)+(yA*yA)+(zA*zA));
    b2 = sqrt((xB*xB)+(yB*yB)+(zB*zB));
    slewAngle = abs(acos((a1b1 + a2b2 + a3b3) / (a2*b2))); 
    
    // check to see if the slew angle is greater than 180 deg
    if (slewAngle > PI) {
        slewAngle = (2*PI) - slewAngle;
    }
    
    // Estimate the extra amount of rotation needed given the sky rotation
    slewTime = slewAngle / RELOCATE_SPEED;
    //cout << "Slew Time: " << slewTime << endl;
    
    position_A = getAsteroidPosition(asteroidTraj_a, viewTime);
    position_B = getAsteroidPosition(asteroidTraj_b, viewTime+slewTime);
    
    xA = position_A.x;
    yA = position_A.y;
    zA = position_A.z;
    xB = position_B.x;
    yB = position_B.y;
    zB = position_B.z;

    a1b1 = xA * xB;
    a2b2 = yA * yB;
    a3b3 = zA * zB;
    a2 = sqrt((xA*xA)+(yA*yA)+(zA*zA));
    b2 = sqrt((xB*xB)+(yB*yB)+(zB*zB));
    slewAngle = abs(acos((a1b1 + a2b2 + a3b3) / (a2*b2))); 
    
    // check to see if the slew angle is greater than 180 deg
    if (slewAngle > PI) {
        slewAngle = (2*PI) - slewAngle;
    }    
    
    return slewAngle;
}
vector<double> AsteroidTracker::generateTimesToRelocate(vector<int> targetIDList, double viewTime) {
    vector<double> timesToRelocate;
    
    _j antenna;
    _i asteroid_a, asteroid_b; 
    
    double slewTime;
    int targetID, newTargetID;
    
    double PI = 3.1415926535897932384626433832795028841971693993751058;

    for (int i = 0; i < numberOfAntennas; i++) {
        antenna = antennas.at(i);
        if (antenna.targetAsteroidID == -1) {
            slewTime = (PI/4) / RELOCATE_SPEED;
            timesToRelocate.push_back(slewTime);
        } else {
            targetID = antenna.targetAsteroidID;
            newTargetID = targetIDList.at(i);
            if (newTargetID != -1) {
                asteroid_a = asteroids.at(targetID);
                asteroid_b = asteroids.at(newTargetID);
            
                slewTime = calcSlewAngle(asteroid_a.trajectory, asteroid_b.trajectory, viewTime)
                    / RELOCATE_SPEED;
                timesToRelocate.push_back(slewTime);
            } else {
                slewTime = (PI/4) / RELOCATE_SPEED;
                timesToRelocate.push_back(slewTime);
            }
            
        }
    }
    
    return timesToRelocate;
}
void AsteroidTracker::relocateAntennas(vector<int> targetID, double timeToSlew, double viewTime) {
    
    _i asteroid;
    _j antenna;
    _traj pos;
    _command cmd;
    ostringstream commandLine;
    int ID;
    
    //cout << "targetID size: " << targetID.size() << "\ttimeToSlew: " 
    //        << timeToSlew << "\tviewTime: " << viewTime << endl;
    
    for (int antennaID = 0; antennaID < numberOfAntennas; antennaID++) {
        
        ID = targetID.at(antennaID);

        if (ID != -1) {
            asteroid = asteroids.at(ID);
            antenna.targetAsteroidID = targetID.at(antennaID);
            antenna.relocating = true;
            antenna.power = 0.0;
            pos = getAsteroidPosition(asteroid.trajectory, viewTime+timeToSlew);
            antenna.Xpoint = pos.x;
            antenna.Ypoint = pos.y;
            antenna.Zpoint = pos.z;
            antenna.nextCommandAvailabilityTime = viewTime + timeToSlew;
            antennas.at(antennaID) = antenna;   // added to v3

            cmd.antennaID = antennaID;
            cmd.commandTime = viewTime;
            cmd.pointDelta = true;
            cmd.powerDelta = false;
            cmd.powerLevel = 0.0 * MAX_TRANSMITTING_POWER;
            cmd.targetAsteroidID = targetID.at(antennaID);

            commandLine.str("");
            commandLine 
                    << cmd.commandTime << " " << cmd.antennaID 
                    << " R " << cmd.targetAsteroidID;

            cmd.commandNote = commandLine.str();
            commandList.push_back(cmd);
        } else {
            int prevID = antennas.at(antennaID).targetAsteroidID;
            if (prevID != -1) {
                antenna.targetAsteroidID = ID;
                antenna.relocating = true;
                antenna.power = 0.0;
                antenna.Xpoint = 0;
                antenna.Ypoint = 0;
                antenna.Zpoint = 1;
                antenna.nextCommandAvailabilityTime = viewTime + timeToSlew;
                antennas.at(antennaID) = antenna;   // added to v3

                cmd.antennaID = antennaID;
                cmd.commandTime = viewTime;
                cmd.pointDelta = true;
                cmd.powerDelta = false;
                cmd.powerLevel = 0.0 * MAX_TRANSMITTING_POWER;
                cmd.targetAsteroidID = ID;
                commandLine.str("");
                commandLine 
                        << cmd.commandTime << " " << cmd.antennaID 
                        << " R " << cmd.targetAsteroidID;

                cmd.commandNote = commandLine.str();
                commandList.push_back(cmd);
            } else {
                // skip null commands ... duplicates of relocating -1 to -1 
            }
        }

        // cout << "Command " << antennaID << ": " << cmd.commandNote << endl;
    }
    
    nextCommandTime = viewTime + timeToSlew;
    
}
double AsteroidTracker::activateAntennas(double viewTime) {
    _i asteroid;
    _j antenna;
    _traj pos;
    _command cmd;
    ostringstream commandLine;
    double beamTime;
    double powerLevel;
    double maxBeamTime = 0;
    
    for (int antennaID = 0; antennaID < numberOfAntennas; antennaID++) {
        
        if (antenna.targetAsteroidID == -1) {
            powerLevel = 0.0;
        } else {
            powerLevel = 1.0;
        }
        //antenna.targetAsteroidID = targetID.at(antennaID);
        antenna.relocating = false;
        antenna.power = powerLevel * MAX_TRANSMITTING_POWER;
        pos = getAsteroidPosition(asteroid.trajectory, viewTime);
        antenna.Xpoint = pos.x;
        antenna.Ypoint = pos.y;
        antenna.Zpoint = pos.z; 
        antennas.at(antennaID) = antenna;   // added to v3
        
        beamTime = signalTravelTime(pos);
        if (beamTime > maxBeamTime) {
            maxBeamTime = beamTime;
        }
        antenna.nextCommandAvailabilityTime = viewTime + beamTime;
        cmd.antennaID = antennaID;
        cmd.commandTime = viewTime;
        cmd.pointDelta = false;
        cmd.powerDelta = true;
        cmd.powerLevel = antenna.power;
        cmd.targetAsteroidID = antenna.targetAsteroidID;
        
        commandLine.str("");
        commandLine 
                << cmd.commandTime << " " << cmd.antennaID 
                << " P " << cmd.powerLevel;

        cmd.commandNote = commandLine.str();
        commandList.push_back(cmd);

        // cout << "Command " << antennaID << ": " << cmd.commandNote << endl;
    }
    
    nextCommandTime = viewTime + maxBeamTime;
    return maxBeamTime;
}
void AsteroidTracker::monitorAsteroids(double beamTime, double viewTime) {
    _i asteroid;
    _j antenna;
    _traj pos;
    _command cmd;
    ostringstream commandLine;
    double powerLevel;
    
    for (int antennaID = 0; antennaID < numberOfAntennas; antennaID++) {
        
        if (antenna.targetAsteroidID == -1) {
            powerLevel = 0.0;
        } else {
            powerLevel = 0.0;
        }
        //antenna.targetAsteroidID = targetID.at(antennaID);
        antenna.relocating = false;
        antenna.power = powerLevel * MAX_TRANSMITTING_POWER;
        pos = getAsteroidPosition(asteroid.trajectory, viewTime);
        antenna.Xpoint = pos.x;
        antenna.Ypoint = pos.y;
        antenna.Zpoint = pos.z;
        antenna.nextCommandAvailabilityTime = viewTime + beamTime;
        antennas.at(antennaID) = antenna;   // added to v3     
        
        cmd.antennaID = antennaID;
        cmd.commandTime = viewTime;
        cmd.pointDelta = false;
        cmd.powerDelta = true;
        cmd.powerLevel = antenna.power;
        cmd.targetAsteroidID = antenna.targetAsteroidID;
        
        commandLine.str("");
        commandLine 
                << cmd.commandTime << " " << cmd.antennaID 
                << " P " << cmd.powerLevel;

        cmd.commandNote = commandLine.str();
        commandList.push_back(cmd);

        // cout << "Command " << antennaID << ": " << cmd.commandNote << endl;
    }
    
    nextCommandTime = viewTime + beamTime;
}
void AsteroidTracker::retireAntennas(void) {

    _j antenna;
    _command cmd;
    
    ostringstream commandLine;
    
    double powerLevel;
    
    for (int antennaID = 0; antennaID < numberOfAntennas; antennaID++) {
        
        if (antenna.targetAsteroidID == -1) {
            powerLevel = 0.0;
        } else {
            powerLevel = 0.0;
        }
        antenna.targetAsteroidID = -1;
        antenna.relocating = false;
        antenna.power = powerLevel * MAX_TRANSMITTING_POWER;
        
        antenna.Xpoint = 0;
        antenna.Ypoint = 0;
        antenna.Zpoint = 1;
        antenna.nextCommandAvailabilityTime = SIMULATION_TIME + 1;
        antennas.at(antennaID) = antenna; // added to v3        
        
        cmd.antennaID = antennaID;
        cmd.commandTime = SIMULATION_TIME + 1;
        cmd.pointDelta = false;
        cmd.powerDelta = true;
        cmd.powerLevel = antenna.power;
        cmd.targetAsteroidID = antenna.targetAsteroidID;
        
        commandLine.str("");
        commandLine 
                << cmd.commandTime << " " << cmd.antennaID 
                << " P " << cmd.powerLevel;

        cmd.commandNote = commandLine.str();
        commandList.push_back(cmd);

        // cout << "Command " << antennaID << ": " << cmd.commandNote << endl;
    }
    
    nextCommandTime = SIMULATION_TIME;
}
// parent functions in AsteroidTracker
int AsteroidTracker::initialize(vector<double> antennaPositions, double peakGain, vector<double> minDistanceGain, vector<double> maxDistanceGain) {
    // initialize the antenna data vectors & parameters
    _j antenna;
    _gain gain;
    
    asteroids.clear();    
    antennas.clear();
    antennaGain.clear();
    commandList.clear();
    
    PEAK_GAIN = peakGain;
    
    commandNumber = 0;
    antennasRetired = false;
    commandDeckLoaded = false;
    
    // build the antenna data vectors
    for (int i = 0; i < antennaPositions.size(); i=i+2) {
        antenna.antennaID = i/2;
        antenna.Xpos = antennaPositions.at(i);
        antenna.Ypos = antennaPositions.at(i+1);
        antenna.targetAsteroidID = -1;
        antenna.Xpoint = 0.0;
        antenna.Ypoint = 0.0;
        antenna.Zpoint = 0.0;
        antenna.power = 0.0;
        antenna.nextCommandAvailabilityTime = 0.0;
        antenna.antennaMaxDistance = 0.0;
        antenna.antennaMinDistance = 0.0;
        antenna.relocating = false;
        antennas.push_back(antenna);
        //cout << "Antenna: " << antenna.antennaID << "\tXpos: " << antenna.Xpos
        //        << "\tYpos: " << antenna.Ypos << endl;
    }
    numberOfAntennas = antennas.size();
    
    // build the antennaDistanceGain min/max vectors
    for (int i = 0; i < minDistanceGain.size(); i++) {
        gain.minDistanceGain = minDistanceGain.at(i);
        gain.maxDistanceGain = maxDistanceGain.at(i);
        antennaGain.push_back(gain);
    }
    // cout << endl << "Total antenna Distance Gain Points: " << antennaGain.size() << endl << endl;
    
    // update the values for antennaMaxDistance & antennaMinDistance
    for (int i = 0; i < antennas.size(); i++) {
        double minDist = 1e20;
        double maxDist = 0.0;
        for (int j = 0; j < antennas.size(); j++) {
            if (i != j) {   
                double dist = distanceBetweenAntennas(i, j);
                if (dist > maxDist) {
                    maxDist = dist;
                }
                if (dist < minDist) {
                    minDist = dist;
                }
             }
        }
        antennas.at(i).antennaMaxDistance = maxDist;
        antennas.at(i).antennaMinDistance = minDist;
        // cout << "Antenna: " << i << "\tMax Dist: " << maxDist << "\tMin Dist: " << minDist << endl;
    }

    return _If_Nacho_Libre_Can_Wake_Up_At_Five_A_M_To_Make_The_Soup_;
}
int AsteroidTracker::asteroidAppearance(int asteroidIndex, double scienceScoreMultiplier, double reflectivityMultiplier, double initialImageInformation, double initialTrajectoryInformation, vector<double> trajectory) {
    
    _i asteroid;
    vector<_traj> traj;
    _traj pos;
    
    asteroid.asteroidID = asteroidIndex;
    asteroid.scienceScoreMultiplier = scienceScoreMultiplier;
    asteroid.reflectivityMultiplier = reflectivityMultiplier;
    asteroid.imageInformation = initialImageInformation;
    asteroid.trajectoryInformation = initialTrajectoryInformation;
    traj.clear();
    for (int i = 0; i < trajectory.size(); i=i+4) {
        pos.t = trajectory.at(i);
        pos.x = trajectory.at(i+1);
        pos.y = trajectory.at(i+2);
        pos.z = trajectory.at(i+3);
        traj.push_back(pos);
    }

    asteroid.trajectory = traj;
    asteroids.push_back(asteroid);
    numberOfAsteroids = asteroids.size();

    return _It_Is_The_Besssssssst_;
}
string AsteroidTracker::nextCommand(double currentTime) {    
    Analysis a;
    string the_big_answer;
    vector<_link> links, rankedLinks;
    vector<int> targetList;
    vector<double> timesToRelocate;
    int commandListSize;
    
    nextCommandTime = currentTime;
    while ((!commandDeckLoaded) && (nextCommandTime < SIMULATION_TIME)) {
        links.clear();
        rankedLinks.clear();
        targetList.clear();
        timesToRelocate.clear();
        
        links = generateLinks(nextCommandTime);
        rankedLinks = rankLinks(links);
        targetList = generateTargetLinks(rankedLinks);
        timesToRelocate = generateTimesToRelocate(targetList, nextCommandTime);
        double maxRelocateTime = a.max(timesToRelocate);
        relocateAntennas(targetList, maxRelocateTime, nextCommandTime);
        double beamTime = activateAntennas(nextCommandTime);
        monitorAsteroids(beamTime, nextCommandTime);
    }
    commandDeckLoaded = true;
    
    if (!antennasRetired) { // set a trap to prevent the deck from being overloaded
        retireAntennas();
        antennasRetired = true;
        commandListSize = commandList.size();
    }

    // fill the output returned with the command note.
    if (commandDeckLoaded) { // Added to v2
        if (commandNumber < commandListSize) {
            the_big_answer = commandList.at(commandNumber).commandNote;
        }
        commandNumber++;
    }

    // evaluate all links to determine the best candidate

    cout << "Number of Antennas: " << numberOfAntennas << endl;    
    cout << "Number of Asteroids: " << numberOfAsteroids << endl;
    cout << "Number of Commands: " << commandList.size() << endl;
    
    return the_big_answer;
}
int main(int argc, char** argv) {
    AsteroidTracker AT;
       
    string readLine;
    vector<string> testFiles;
    
    int numAntennas, numAsteroids, numSightings, numGainData;
    double x, y, z, t, val, val1, val2, val3, val4, peakGain;
    vector<double> antennaPosition;
    vector<double> minDistanceGain, maxDistanceGain;
    vector<double> scienceScoreMultiplier, reflectivityMultiplier, initialImageInformation, initialTrajectoryInformation;
    vector<double> asteroidTrajectory;
    
    // read data from test case file
    string pathToTestData = "/Users/scott/Arctria/ARC-32/AsteroidTracker/data/testcase";
    string testDataFileExtension = ".txt";
    int numTestCases = 1;//110;
    
    for (int testCase = 1; testCase <= numTestCases; testCase++) {
        ostringstream os;
        os << pathToTestData << testCase << testDataFileExtension;
        string testFile = os.str();
        testFiles.push_back(testFile);
    }

    for (int test = 0; test < testFiles.size(); test++) {
        antennaPosition.clear();
        minDistanceGain.clear();
        maxDistanceGain.clear();
        scienceScoreMultiplier.clear();
        reflectivityMultiplier.clear();
        initialImageInformation.clear();
        initialTrajectoryInformation.clear();
        asteroidTrajectory.clear();

        cout << "======================================================================================" << endl;         
        cout << "Opening test data file: " << testFiles.at(test) << endl;
        ifstream testDataFile (testFiles.at(test));

        // number of antennas
        getline(testDataFile, readLine);
        istringstream iss_1(readLine);
        iss_1 >> numAntennas;
        //cout << "Number of Antennas: " << numAntennas << endl;
            
        // antenna positions: [x, y] x number of antennas
        for (int i = 0; i < numAntennas; i++) {
            getline(testDataFile, readLine);
            istringstream iss_2(readLine);
            iss_2 >> x >> y;
            antennaPosition.push_back(x);
            antennaPosition.push_back(y);
        }

        // peak gain
        getline(testDataFile, readLine);
        istringstream iss_3(readLine);
        iss_3 >> peakGain;

        // length of vector for minDistanceGain & maxDistanceGain
        getline(testDataFile, readLine);
        istringstream iss_4(readLine);
        iss_4 >> numGainData;
       
        // minDistanceGain (x length)
        for (int i = 0; i < numGainData; i++) {
            getline(testDataFile, readLine, ' ');
            istringstream iss_5(readLine);
            iss_5 >> val;
            minDistanceGain.push_back(val);
        }

        // maxDistanceGain (x length)
        for (int i = 0; i < numGainData; i++) {
            getline(testDataFile, readLine, ' ');
            istringstream iss_6(readLine);
            iss_6 >> val;
            maxDistanceGain.push_back(val);
        }

        getline(testDataFile, readLine);    // need to add an extra readLine             
            
        // number of asteroids
        getline(testDataFile, readLine);
        istringstream iss_7(readLine);
        iss_7 >> numAsteroids;
        
        AT.initialize(antennaPosition, peakGain, minDistanceGain, maxDistanceGain);

        // scienceScoreMultiplier (x number of asteroids)
        // reflectivityMultiplier (x number of asteroids)
        // initialImageInformation (x number of asteroids)
        // initialTrajectoryInformation (x number of asteroids)
        // number of asteroid sightings   

        for (int i = 0; i < numAsteroids; i++) {

            getline(testDataFile, readLine);
            istringstream iss_8(readLine);
            iss_8 >> val1 >> val2 >> val3 >> val4 >> numSightings;
            scienceScoreMultiplier.push_back(val1);
            reflectivityMultiplier.push_back(val2);
            initialImageInformation.push_back(val3);
            initialTrajectoryInformation.push_back(val4);
            asteroidTrajectory.clear();

            // asteroid trajectory: [time, x, y, z] (x number of asteroid sightings )
            for (int j = 0; j < numSightings; j++) {
                getline(testDataFile, readLine);
                istringstream iss_9(readLine);
                iss_9 >> t >> x >> y >> z;
                asteroidTrajectory.push_back(t);
                asteroidTrajectory.push_back(x);
                asteroidTrajectory.push_back(y);
                asteroidTrajectory.push_back(z);
            }
            
            AT.asteroidAppearance(i, scienceScoreMultiplier.at(i), 
            reflectivityMultiplier.at(i), initialImageInformation.at(i), 
            initialTrajectoryInformation.at(i),asteroidTrajectory);            
        }        
   
        testDataFile.close();
        
        AT.nextCommand(0.0);      
    }

    return 0;
}
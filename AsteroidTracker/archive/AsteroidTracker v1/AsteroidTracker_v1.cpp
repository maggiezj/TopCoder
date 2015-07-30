/* AsteroidTracker [Competition Sensitive]
 * NASA Asteroid Grand Challenge | NASA | Planetary Resources | TopCoder
 * RigelFive - Hudson Ohio USA
 * 25 July 2014 - 8 August 2014
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
  return get<1>(lhs) > get<1>(rhs); 
}

struct _traj { // trajectory of an asteroid
    int asteroid_id;
    double t;
    double x;
    double y;
    double z;
    double distance;
    double travelTime;
    bool visible; 
    bool ascending;
};
struct _signal {
    int antenna_id;
    double arrivalTime; // time of arrival for signal from asteroid
    double signalPower;
};
struct _i { // asteroid information
     int asteroidIndex;
     double scienceScoreMultiplier;
     double reflectivityMultiplier; 
     double imageInformation; 
     double trajectoryInformation; 
     double trajInfoScore;
     double timeAppear;

     vector<_traj> trajectory;
     vector<_signal> signalStrength;    
};
struct _j { // antenna information
    // static parameters
    int antenna_id;
    double Xpos;
    double Ypos;
    
    // dynamic parameters v. time
    double t; // time
    int state; // dynamic states:  0 = antenna off, 1 = antenna transit, 2 = antenna emitting, 3 = antenna receiving 
    double power; // dynamic power setting
    double Xpoint; // dynamic X pointing coordinate
    double Ypoint; // dynamic Y pointing coordinate
    double Zpoint; // dynamic Z pointing coordinate
    int asteroid_id; // dynamic condition as to which asteroid the antenna is pointing at
    bool transmitting;
    bool relocating;
    double antennaMinDistance;
    double antennaMaxDistance;
    double nextCommandAvailabilityTime; 
};
struct _gain {
    int antenna_id;
    double minDistanceGain;
    double maxDistanceGain;
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

class AsteroidTracker {
public:
    // primary functions and methods
    int initialize(vector<double> antennaPositions, double peakGain, vector<double> minDistanceGain, vector<double> maxDistanceGain);
    int asteroidAppearance(int asteroidIndex, double scienceScoreMultiplier, double reflectivityMultiplier, double initialImageInformation, double initialTrajectoryInformation, vector<double> trajectory);
    string nextCommand(double currentTime);
    
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
  
    // objective functions adapted from Simulator provided in the AsteroidTracker problem statement
    vector<double> getAsteroidPositions(int asteroidIndex); // done - modified
    double getScore(void); // done
    void updateOutputSignalPower(void); // done
    void updateRelocatingAntennasDirections(double deltaTime); // done
    double updateScore(double deltaTime);
    void updateTrackingAntennasDirections(void); // done
    
    // methods to determine fitness (in some cases, adapted from Asteroid Tracker Simulator)
    double distanceBetweenAntennas(int antennaID_1, int antennaID_2); // verified
    double distanceToAsteroid(double x, double y, double z); // verified
    double signalTravelTime(double distToAsteroid); // verified
    double timeOfAppearance(vector<_traj> asteroidTrajectory); // verified
    vector<int> getSubarray(int asteroidID); // verified   
    _traj getAsteroidPosition(int asteroidID, double viewTime); //
    bool isAsteroidVisible(int asteroidID, double viewTime); // verified
    bool isAntennaRelocating(int antennaID, double viewTime); // verified
    bool isAntennaTransmitting(int antennaID, double viewTime); // verified
    bool isAscendingFromHorizon(int asteroidID, double viewTime); //
    bool checkTargetProximity(int asteroidID_A, int asteroidID_B, double viewTime); // verified
    bool checkAntennaBeamIntersection(int transmittingAntennaID, double viewTime); // verified - but no events???
    bool checkMaxPowerLimit(int antennaID, double viewTime); // verified
    double calcSlewAngle(int targetAsteroidID_A, int targetAsteroidID_B, double viewTime); // verified
    double calcPowerTransmitted(int asteroidID, double viewTime);
    double effectiveSignalPowerReturn(int asteroidIndex, double viewTime); // 
    double interpolateGain(double distanceBetweenAntennas, double angle); // 
    double inducedNoise(int receivingAntennaIndex, double viewTime); // 
    double calcTotalNoise(int asteroidID, double viewTime); // 
    double calcInformationRate(int asteroidID, double viewTime); // 
    double calcTrajectoryInformation(int asteroidID, double viewTime, double deltaTime); //
    
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
    vector<_gain> antennaDistanceGain;
    
    vector<double> antennaMinDistance, antennaMaxDistance;     
    
    // Simulator parameters
    int numberOfAntennas = 0;
    int numberOfAsteroids = 0; 
    
    double energySpent = 0.0;
    double simTime = 0.0;
    double score = 0.0;
    
    /* Note: 
     * outputPowerQueue -> ArrayDeque(Collection<E> c): 
     * constructs a deque containing the elements of the specified collection, 
     * in the order they are returned by the collection's iterator.  
     * The first element returned by the collection iterator becomes the first 
     * element or the front of the deque.
     */

private:
    
};

void AsteroidTracker::generate_optimization(int numAsteroids, int numAntennas, int mutate, int popSize, double maxCycle) {
    _pop pop;
    vector<_gn> genome;
    vector<_pop> population, population_competed, population_ranked;

    int gen = 0;
    int j = 0;
        
    population.clear();
    population = create_population(numAsteroids, numAntennas, popSize); 
    double score = 0.0;

    while (j < maxCycle) {
        
        gen++;
        population_competed.clear();
        population_ranked.clear();
        population_competed = compete_population(population);
        population_ranked = select_population(population_competed);
        population = evolve_population(population_ranked, mutate);
        j++;
        double newScore = population.at(0).fitness / numAntennas;
        if (score < newScore) {
            score = newScore;
            cout << "Generation #: " << gen << " Fitness: " << score << endl;
            j=0;
        }
    }

    cout << "Mutations: " << mutate << "\tGeneration #: " << gen << "\tBest Fitness: " << score << endl;  

    pop = population.at(0);
    genome = pop.genome;
    int gene_groups = genome.size() / 2;
    _gn gene_1, gene_2;

    for (int i = 0; i < gene_groups; i++) {
        int gene_num = i*2;
        gene_1 = genome.at(gene_num+0);
        gene_2 = genome.at(gene_num+1);

        cout << "Antenna #: " << i << endl;
        cout << "targetID: " << round(create_64b_value(gene_1.gene, gene_1.min_value, gene_1.max_value)) << endl
            << "powerLevel: " << round(create_64b_value(gene_2.gene, gene_2.min_value, gene_2.max_value)) << endl << endl;

    }
}
double AsteroidTracker::calculateFitness(int targetID, double powerLevel) {
    /*  
     * 1. antenna distances - distanceBetweenAntennas (antenna struct _j)
     * 2. distance to asteroid - distanceToAsteroid (asteroid struct _i)
     * 3. signal travel time (s) - signalTravelTime (asteroid struct _i)
     * 4. time of asteroid appearance (timeAppear):  - timeOfAppearance (asteroid struct _i)
     * 5. subarray (antennaIDs targeting asteroidID) - getSubarray
     * 6. VISIBILITY CONSTRAINT - isAsteroidVisible
     * 7. ANTENNA RELOCATING CONSTRAINT - isAntennaRelocating
     * 8. IS TRANSMITTING CONSTRAINT - isAntennaTransmitting
     * 9. TARGET PROXIMITY CONSTRAINT - checkTargetProximity
     * 10. NEAR-FIELD CONSTRAINT - checkAntennaBeamIntersection
     * 11. MAXIMUM TRANSMITTING POWER CONSTRAINT - checkMaxPowerLimit
     * 12. slewAngle & slewTime - calcSlewAngle
     * 13. transmittedPower - calcPowerTransmitted
     * 14. effectiveSignalPowerReturn - effectiveSignalPowerReturn
     * 15. interpolateGain  - interpolateGain
     * 16. inducedNoise - inducedNoise
     * 17. totalNoise - calcTotalNoise
     * 18. informationRate - calcInformationRate
     * 19. trajectoryInformation - calcTrajectoryInformation
     * 20. energySpent - [*] need to calculate
     * 21. imageKnowledgeScore - [*] need to calculate
     * 22. trajectoryKnowledgeScore - [*] need to calculate
     * 23. finalScore - [*] need to calculate
     */ 
    
    /*double fitness = 1000.0 - max(
            optimizationError(1.0, targetID),
            optimizationError(25.0, powerLevel)
            );  */
    _i asteroid;
    int asteroidID;
    double fitness = 0.0;

    if (targetID > asteroids.back().asteroidIndex) {
        asteroidID = asteroids.back().asteroidIndex;
    } else {
        asteroidID = (int) round(targetID);
    }

    asteroid = asteroids.at(asteroidID);
    double viewTime = simTime;
    double D = getAsteroidPosition(asteroidID, viewTime).distance;
    double truth = 0.0;
    if (isAscendingFromHorizon(asteroidID,viewTime)
            && isAsteroidVisible(asteroidID, viewTime)
            && (D > 1000.0)) {
        truth = 1.0;
        fitness = truth 
            //    * pow(antennas.size(),2.0)
                * powerLevel * MAX_TRANSMITTING_POWER * 180.0
                * ((asteroid.trajectoryInformation * TRAJECTORY_SCORE_MULTIPLIER
                + (asteroid.imageInformation * IMAGE_SCORE_MULTIPLIER))
                * asteroid.scienceScoreMultiplier 
                * PEAK_GAIN 
                * asteroid.reflectivityMultiplier)
                / (D*D);
    }
    // cout << "asteroidID: " << asteroidID << "\tMax: " << asteroids.size() << "\tDist: " << D << endl;

    //cout << "fitness: " << fitness << endl;

    return fitness;
}
double AsteroidTracker::optimizationError(double targetVal, double currentVal) {
    double targError;
    
    targError = abs((targetVal-currentVal)/targetVal) * 100.0;
    return targError;

}
vector<_pop> AsteroidTracker::create_population(int number_of_asteroids, int number_of_antennas, int population_size) {
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
vector<_pop> AsteroidTracker::compete_population(vector<_pop> population) {
// determine the fitness of each member of the population

    _gn gn;
    _pop pop;
    vector<_pop> new_population;
    double targetID, powerLevel, deltaTime;
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
            powerLevel = round(create_64b_value(
                    pop.genome.at(gene_num+1).gene,
                    pop.genome.at(gene_num+1).min_value,
                    pop.genome.at(gene_num+1).max_value
                ));

            pop.fitness = pop.fitness + calculateFitness(targetID, powerLevel);
        }

        new_population.push_back(pop);
    }    
    
    return new_population;

}
vector<_pop> AsteroidTracker::select_population(vector<_pop> population) {
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
    sort(test_tuple.begin(),test_tuple.end(),sort_object_score); 
    
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
vector<_pop> AsteroidTracker::evolve_population(vector<_pop> population, int mutate) {
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

vector<double> AsteroidTracker::create_values(double min, double max, int values) {
    vector<double> random_values;
    
    default_random_engine rand_engine;
    uniform_real_distribution<double> distribution(min, max);    
 
    for (int i = 0; i < values; i++) {
        random_values.push_back(distribution(rand_engine));
    }
    
    return random_values;
}
bitset<64> AsteroidTracker::create_64b_gene(double v, double min, double max) {    
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
double AsteroidTracker::create_64b_value(bitset<64> g, double min, double max) {
    double value;
    double v, v_max;
    bitset<64> b_64_max;
    
    b_64_max.set();
    v_max = (double) b_64_max.to_ullong();
    v = (double) g.to_ullong();
    value = ((v/v_max) * (max-min)) + min;
    return value;
}

vector<_gn> AsteroidTracker::mutate_genome(vector<_gn> old_genome, int mutate) {
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
_gn2 AsteroidTracker::cross_genomes(_gn2 &old_genome_pair) {

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
void AsteroidTracker::output_genome(vector<_gn> genome) {
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

vector<double> AsteroidTracker::getAsteroidPositions(int asteroidIndex) {
    /*
     * Get current asteroid position, given an asteroid index
     */    
    _i asteroid;
    _traj trajectory;
    vector<double> asteroid_position;
    
    asteroid = asteroids.at(asteroidIndex);    
    trajectory = asteroid.trajectory.back();

    double x = trajectory.x;
    double y = trajectory.y;
    double z = trajectory.z;
    double distance = trajectory.distance;
    
    asteroid_position.push_back(x);
    asteroid_position.push_back(y);
    asteroid_position.push_back(z);
    asteroid_position.push_back(distance);
  
    return asteroid_position;
}
double AsteroidTracker::getScore(void) {
    /*
     * Return the current total score
     */
    
    _i asteroid;
    double score = 0.0;
    double imageInfoScore, trajInfoScore;
    int id;

    for (int i = 0; i < numberOfAsteroids; i++) {
        asteroid = asteroids.at(i);
        
        id = asteroid.asteroidIndex;
        imageInfoScore = tanh(asteroid.imageInformation / Q_IMAGE);
        trajInfoScore = asteroid.trajInfoScore;
        
        score = score 
                + asteroid.scienceScoreMultiplier * ((IMAGE_SCORE_MULTIPLIER * imageInfoScore)
                + (TRAJECTORY_SCORE_MULTIPLIER * trajInfoScore));
    }
    
    score = score - (energySpent * 1.0e-9);
    return score;

}

void AsteroidTracker::updateOutputSignalPower(void) {
    /*
     * Update current transmitted signal output power
     */
    
    _i asteroid;
    _j antenna;
    _signal signal;
    vector<_signal> signals;
    
    bool visible;
    vector<int> trackingAntennas;
    double amplitudeSum, combinedOutputPower, R, travelTime;
    int antennaID;
    
    for (int asteroidID = 0; asteroidID < asteroids.size(); asteroidID++) {
        asteroid = asteroids.at(asteroidID);
        visible = asteroid.trajectory.at(0).visible;
        if (visible == false) {
            
        } else {
            trackingAntennas = getSubarray(asteroidID);
            if (trackingAntennas.size() > 0) {
                amplitudeSum = 0.0;
                for (int linkID = 0; linkID < trackingAntennas.size(); linkID++) {
                    antennaID = trackingAntennas.at(linkID);
                    amplitudeSum += sqrt(antennas.at(antennaID).power);
                }
                combinedOutputPower = amplitudeSum * amplitudeSum;
                R = asteroid.trajectory.back().distance; // last data point in trajectory vector defines current distance
                travelTime = 2.0 * R / SPEED_OF_LIGHT;
                signals = asteroid.signalStrength;
                for (int linkID = 0; linkID < trackingAntennas.size(); linkID++) {
                    signal.antenna_id = trackingAntennas.at(linkID);
                    signal.signalPower = combinedOutputPower;
                    signal.arrivalTime = simTime + travelTime;
                    signals.push_back(signal);
                    // need to add a tuple sort for power v. time for the signals on the asteroid
                    // make sure the vector is not being overloaded - do garbage collection
    
                }
                asteroid.signalStrength = signals;
                asteroids.at(asteroidID) = asteroid;
            }
        }
    }
}
void AsteroidTracker::updateRelocatingAntennasDirections(double deltaTime) {
    /*
     * Update relocating antennas directions given elapsed time 
     */
    
    _j antenna;
    
    double slewAngle;
    double ant_X, ant_Y, ant_Z, ast_X, ast_Y, ast_Z;
    double rot_X, rot_Y, rot_Z, mag_rot, uX, uY, uZ;
    double cosA, sinA, oneMinus_cosA;
    double m00, m01, m02, m10, m11, m12, m20, m21, m22;   
    double newX, newY, newZ, mag_newVector;
    
    vector<double> asteroid_position;
    
    if (deltaTime == 0.0) {
        return;
    }

    for (int antennaID = 0; antennaID < antennas.size(); antennaID++) {
        antenna = antennas.at(antennaID);
        
        if ((antenna.asteroid_id == -1) || (antenna.relocating == false)) {
            
        } else {
            // initialize constants for this problem
            slewAngle = deltaTime * RELOCATE_SPEED;
            cosA = cos(slewAngle);
            sinA = sin(slewAngle);
            oneMinus_cosA = 1.0 - cosA;
            
            // cross product of the antenna pointing vector x asteroid vector
            ant_X = antenna.Xpoint;
            ant_Y = antenna.Ypoint;
            ant_Z = antenna.Zpoint;
            
            asteroid_position = getAsteroidPositions(antenna.asteroid_id);
            ast_X = asteroid_position.at(0);
            ast_Y = asteroid_position.at(1);
            ast_Z = asteroid_position.at(2);
            
            // determine the axis of rotation for the antenna relocation
            rot_X = (ant_Y*ast_Z) - (ant_Z*ast_Y); // u2v3 - u3v2
            rot_Y = (ant_Z*ast_X) - (ant_X*ast_Z); // u3v1 - u1v3
            rot_Z = (ant_X*ast_Y) - (ant_Y*ast_X); // u1v2 - u2v1
            
            mag_rot = sqrt((rot_X*rot_X) + (rot_Y*rot_Y) + (rot_Z*rot_Z));
            
            // normalize the antenna rotation vector
            uX = rot_X / mag_rot;
            uY = rot_Y / mag_rot;
            uZ = rot_Z / mag_rot;

            // rotation matrix for rotation around rotationAxis with specified angle
            m00 = cosA + uX * uX * oneMinus_cosA;
            m01 = uX * uY * oneMinus_cosA - uZ * sinA;
            m02 = uX * uZ * oneMinus_cosA + uY * sinA;
            m10 = uX * uY * oneMinus_cosA + uZ * sinA;
            m11 = cosA + uY * uY * oneMinus_cosA;
            m12 = uY * uZ * oneMinus_cosA - uX * sinA;
            m20 = uX * uZ * oneMinus_cosA - uY * sinA;
            m21 = uY * uZ * oneMinus_cosA + uX * sinA;
            m22 = cosA + uZ * uZ * oneMinus_cosA;

            // calculate the new antenna pointing vector
            newX = m00 * ant_X + m01 * ant_Y + m02 * ant_Z;
            newY = m10 * ant_X + m11 * ant_Y + m12 * ant_Z;
            newZ = m20 * ant_X + m21 * ant_Y + m22 * ant_Z;  
            
            mag_newVector = sqrt((newX*newX) + (newY*newY) + (newZ*newZ));
     
            // calculate the normalized values of the new antenna pointing vector
            antenna.Xpoint = newX / mag_newVector;
            antenna.Ypoint = newY / mag_newVector;
            antenna.Zpoint = newZ / mag_newVector; 
     
            // ??? need to determine if the antenna has already arrived
            // at the desired angle and is no longer 'relocating'
            
            antennas.at(antennaID) = antenna; // update the relocating antenna pointing vector
        }
    }
}
double AsteroidTracker::updateScore(double deltaTime) {
    /*
     * Update image/trajectory information, energy spent and scores, given elapsed time
     */
    
    _i asteroid;
    _j antenna;
    
    double decayMultiplier, informationRate, eSPR, inducedNoiseSum, noise, gain;
    
    vector<int> antennaLinks;
    int N_receive;
    
    if (deltaTime == 0.0) {
        return score;
    } else {
        for (int antennaID = 0; antennaID < antennas.size(); antennaID++) {
            antenna = antennas.at(antennaID);
            if (antenna.asteroid_id == -1) {
                continue;
            }
            if (antenna.relocating == true) {
                energySpent += deltaTime * RELOCATE_POWER;
            } else {
                energySpent += deltaTime * antenna.power;
            }

        }
    }
    
    decayMultiplier = exp(-deltaTime / Q_LOST);
    score = 0.0;
    
    for (int asteroidID = 0; asteroidID < asteroids.size(); asteroidID++) {
        asteroid = asteroids.at(asteroidID);
        antennaLinks = getSubarray(asteroidID);
        N_receive = antennaLinks.size();
        informationRate = 0.0;
        if (N_receive > 0.0) {
            eSPR = effectiveSignalPowerReturn(asteroidID, N_receive);
            inducedNoiseSum = 0.0;
            for (int antennaID = 0; antennaID < antennas.size(); antennaID++) {
                antenna = antennas.at(antennaID);
                if (antenna.asteroid_id == asteroidID) {
                    inducedNoiseSum += inducedNoise(antennaID, 0); //     added a viewTime argument to inducedNoise...     
                }
            }           

            noise = (SHIELDING_MULTIPLIER * inducedNoiseSum) 
                    + (N_receive * BACKGROUND_NOISE);
            
            gain = eSPR / noise;
            if (gain > 1.0) {
                informationRate = 10.0 * log10(gain);
            }
        }
        
        asteroid.imageInformation += deltaTime * informationRate;

        double prevTrajInfo = asteroid.trajectoryInformation;
        asteroid.trajectoryInformation = 
                decayMultiplier * prevTrajInfo
                + (1.0 - decayMultiplier) * Q_LOST * informationRate;
        asteroid.trajInfoScore += 
                deltaTime * (tanh(prevTrajInfo / Q_TRAJECTORY)
                + tanh(asteroid.trajectoryInformation / Q_TRAJECTORY)) / 2.0;

        asteroids.at(asteroidID) = asteroid;
        score += asteroid.trajInfoScore;
    }
    return score;
}
void AsteroidTracker::updateTrackingAntennasDirections(void) {
    /*
     * Update tracking antennas directions 
     */
    
    _i asteroid;
    _j antenna;
    _traj traj;
    int asteroidIndex;    
    
    for (int i = 0; i < antennas.size(); i++) {
        antenna = antennas.at(i);
        asteroidIndex = antenna.asteroid_id;
        
        if ((asteroidIndex == -1) || (antenna.relocating)) {
            
        } else {
            asteroid = asteroids.at(asteroidIndex);
            traj = asteroid.trajectory.back();
            antenna.Xpoint = (traj.x/traj.distance);
            antenna.Ypoint = (traj.y/traj.distance);
            antenna.Zpoint = (traj.z/traj.distance);
            antennas.at(i) = antenna;
        }
    }
}


// AsteroidTracker methods to determine the objective function fitness
double AsteroidTracker::distanceBetweenAntennas(int antennaID_1, int antennaID_2) {
    double dx, dy, dist;
    
    dx = antennas.at(antennaID_1).Xpos - antennas.at(antennaID_2).Xpos;
    dy = antennas.at(antennaID_1).Ypos - antennas.at(antennaID_2).Ypos;
    dist = sqrt((dx*dx)+(dy*dy));
    
    return dist;
}
double AsteroidTracker::distanceToAsteroid(double x, double y, double z) {
    double distance;
    
    distance = sqrt((x*x)+(y*y)+(z*z));
    
    return distance;
}
double AsteroidTracker::signalTravelTime(double distToAsteroid) {
    double travelTime;
    
    travelTime = (2.0 * distToAsteroid) / SPEED_OF_LIGHT;
    return travelTime;
}
double AsteroidTracker::timeOfAppearance(vector<_traj> asteroidTrajectory) {
    double timeAppear;

    int appearIndex;
    _traj trajAsteroid;
    
    timeAppear = 1e20;
    
    for (int i = 0; i < asteroidTrajectory.size(); i++) {
        trajAsteroid = asteroidTrajectory.at(i);
        
        //cout << "t: " << trajAsteroid.t << "\tx: " << trajAsteroid.x << "\ty: " << trajAsteroid.y << "\tz: " << trajAsteroid.z << endl;

        if (trajAsteroid.t < timeAppear) {
            timeAppear = trajAsteroid.t;
            appearIndex = i;
        }
    }
    
    //cout << "Time of Appearance: " << timeAppear << "\tIndex: " << appearIndex  << endl;

    return timeAppear;
}
vector<int> AsteroidTracker::getSubarray(int asteroidID) {
    /*
     * For each asteroid, generate the list of antennas tracking it
     * 
     * Note:
     * Modified java test code to simplify the approach.  Returning a vector<int>
     * of antennas pointed at a specific asteroid.
     */ 
    
    _i asteroid;
    _j antenna;
    
    vector<int> tracked_by_antenna;

    // load the asteroid of interest to see if there are antennas pointed at it.
    asteroid = asteroids.at(asteroidID);

    for (int j = 0; j < antennas.size(); j++) {
        antenna = antennas.at(j);

        // if the asteroid is being tracked by an antenna, return the antennaID in the array.
        if (antenna.asteroid_id == asteroidID) {
            tracked_by_antenna.push_back(antenna.antenna_id);
        }
    }

    return tracked_by_antenna;
}
_traj AsteroidTracker::getAsteroidPosition(int asteroidID, double viewTime) {
    _traj pos;
    
    vector<_traj> asteroidTrajectory;
    double t1, x1, y1, z1, t2, x2, y2, z2;
    
    _i asteroid;
    
    asteroid = asteroids.at(asteroidID);
    asteroidTrajectory = asteroid.trajectory;
    
    for (int i = 1; i < asteroidTrajectory.size(); i++) {
        t1 = asteroidTrajectory.at(i-1).t;
        t2 = asteroidTrajectory.at(i).t;
        
        if ((t2 >= viewTime) && (t1 <= viewTime)) {
            x1 = asteroidTrajectory.at(i-1).x;
            y1 = asteroidTrajectory.at(i-1).y;
            z1 = asteroidTrajectory.at(i-1).z;
            
            x2 = asteroidTrajectory.at(i).x;
            y2 = asteroidTrajectory.at(i).y;
            z2 = asteroidTrajectory.at(i).z;
            // the trajectory data points between the time point viewTime has been found
            pos.asteroid_id = asteroidID;
            pos.t = viewTime;
            pos.x = ((x2-x1)/(t2-t1))*viewTime+(x1-(((x2-x1)/(t2-t1))*t1));
            pos.y = ((y2-y1)/(t2-t1))*viewTime+(y1-(((y2-y1)/(t2-t1))*t1));
            pos.z = ((z2-z1)/(t2-t1))*viewTime+(z1-(((z2-z1)/(t2-t1))*t1));
            pos.distance = distanceToAsteroid(pos.x, pos.y, pos.z);
            pos.visible = isAsteroidVisible(asteroidID, viewTime);
            pos.travelTime = signalTravelTime(pos.distance);            
        }
    }
    
    return pos;

}
bool AsteroidTracker::isAsteroidVisible(int asteroidID, double viewTime) {
    bool visible;
    _i asteroid;
    vector<_traj> asteroidTrajVector;
    _traj traj1, traj2;
    double z1, z2, t1, t2, t_min, t_max;
    
    asteroid = asteroids.at(asteroidID);
    asteroidTrajVector = asteroid.trajectory;
    
    t_min = asteroidTrajVector.at(0).t;
    t_max = asteroidTrajVector.back().t;
    
    visible = false;  // default to false
    
    for (int i = 1; i < asteroidTrajVector.size(); i++) {
        traj1 = asteroidTrajVector.at(i-1);
        traj2 = asteroidTrajVector.at(i);
        
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

bool AsteroidTracker::isAntennaTransmitting(int antennaID, double viewTime) {
    /*
     * Check to see if the antenna is transmitting (power > 0)
     */
    
    _j antenna;
    
    bool transmitting = false;    
    antenna = antennas.at(antennaID); 

    if (viewTime <= antenna.nextCommandAvailabilityTime)  { 
        // just read the current state if the command has been set for an amount of time.
        if (antenna.power > 0.0) {
            transmitting = true;
        }
    } else { // update the current state to non-relocating with availability for the next command to now.
        antenna.power = 0.0;
        antenna.nextCommandAvailabilityTime = viewTime;
        antennas.at(antennaID) = antenna;
    }    
   
    return transmitting;
}
bool AsteroidTracker::isAscendingFromHorizon(int asteroidID, double viewTime) {
    bool ascending = false;
    
    _i asteroid;
    _traj pos1, pos2;
    
    asteroid = asteroids.at(asteroidID);
    pos1 = getAsteroidPosition(asteroidID, viewTime+180);
    pos2 = getAsteroidPosition(asteroidID, viewTime);
    
    if (pos1.z > pos2.z) {
        ascending = true;
    }
    
    return ascending;
}
bool AsteroidTracker::checkTargetProximity(int asteroidID_A, int asteroidID_B, double viewTime) {
    /*
     * Check the target proximity
     */
    
    _i asteroid_a;
    _i asteroid_b;
    _traj traj_a;
    _traj traj_b;
    
    vector<int> subarrays_A;
    vector<int> subarrays_B;
    
    bool asteroids_too_close = false;
    
    double a1b1, a2b2, a3b3, a2, b2, angle;
    
    subarrays_A = getSubarray(asteroidID_A);
    subarrays_B = getSubarray(asteroidID_B);

    if ((subarrays_A.size() == 0) || (subarrays_B.size() == 0)) {
        
    } else {
        asteroid_a = asteroids.at(asteroidID_A);
        asteroid_b = asteroids.at(asteroidID_B);
        
        traj_a = getAsteroidPosition(asteroidID_A,viewTime);
        traj_b = getAsteroidPosition(asteroidID_B,viewTime);
        
        a1b1 = traj_a.x * traj_b.x;
        a2b2 = traj_a.y * traj_b.y;
        a3b3 = traj_a.z * traj_b.z;
        
        a2 = sqrt((traj_a.x*traj_a.x)+(traj_a.y*traj_a.y)+(traj_a.z*traj_a.z));
        b2 = sqrt((traj_b.x*traj_b.x)+(traj_b.y*traj_b.y)+(traj_b.z*traj_b.z));
        
        angle = acos((a1b1 + a2b2 + a3b3) / (a2*b2));
    }

    double minAngle = CRITICAL_TRACKING_ANGLE / min(subarrays_A.size(), subarrays_B.size());
    
    if (angle < minAngle) {
        asteroids_too_close = true;
    }

    return asteroids_too_close;
}
bool AsteroidTracker::checkAntennaBeamIntersection(int transmittingAntennaID, double viewTime) {
    /*
     * Check if the beam from one antenna intersects any other antenna
     */
    
    _i targetAsteroid_T;
    _j antenna_R; // receiving antenna
    _j antenna_T; // transmitting antenna
    _traj asteroidPos;
    bool beam_too_close = false;
    
    double beamAngle, beamDistance, kindaCloseAngle;
    double a1b1, a2b2, a3b3, a2, b2, angle, dX_a, dY_a, dX_b, dY_b, dZ_b;
    
    double PI = 3.1415926535897932384626433832795028841971693993751058;
        
    antenna_T = antennas.at(transmittingAntennaID);
    if (antenna_T.transmitting == false) {

    } else {

        for (int antennaID_r = 0; antennaID_r < antennas.size(); antennaID_r++) {
            antenna_R = antennas.at(antennaID_r);
            if (antennaID_r == transmittingAntennaID) {
                continue;
            } else {
                asteroidPos = getAsteroidPosition(antenna_T.asteroid_id, viewTime);
                   
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
                
                kindaCloseAngle = PI/4.0; // meh, let's say about 45 deg
                
                // if the beam is kinda close, check to see if the beam is too close
                if (abs(beamAngle) < kindaCloseAngle) {
                    beamDistance = a2*sin(beamAngle);
                    if (beamDistance < ANTENNA_SAFE_RADIUS) {
                        beam_too_close = true;
                    }
                }
            }
        }
    }
    
    return beam_too_close;
}
bool AsteroidTracker::checkMaxPowerLimit(int antennaID, double viewTime) {
    _j antenna;
    
    bool powerLimit = false;
    
    antenna = antennas.at(antennaID);
    if ((antenna.power > MAX_TRANSMITTING_POWER) || (antenna.power < 0.0)) {
        powerLimit = true;
    }
        
    return powerLimit;
}
double AsteroidTracker::calcSlewAngle(int targetAsteroidID_A, int targetAsteroidID_B, double viewTime) {
    
    double PI = 3.1415926535897932384626433832795028841971693993751058;
    double slewAngle, slewTime;
    
    double xA, yA, zA, xB, yB, zB;
    
    _i asteroid_A;
    _i asteroid_B;
    _traj position_A;
    _traj position_B;
    
    asteroid_A = asteroids.at(targetAsteroidID_A);
    asteroid_B = asteroids.at(targetAsteroidID_B);
    
    position_A = getAsteroidPosition(targetAsteroidID_A, viewTime);
    position_B = getAsteroidPosition(targetAsteroidID_B, viewTime);
    
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
    
    position_A = getAsteroidPosition(targetAsteroidID_A, viewTime);
    position_B = getAsteroidPosition(targetAsteroidID_B, viewTime+slewTime); 
    
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
double AsteroidTracker::calcPowerTransmitted(int asteroidID, double viewTime) {
    double powerTransmitted = 0;
    double powerSum = 0;
    vector<int> antennasTargetingAsteroid;
    
    _i asteroid;
    _signal signalForAntennaX; // individual signal at the asteroid for all antennas (unorganized)
    vector<_signal> signalStrength;
    vector<_signal> signalsFromAntennaX; // signals returned from the asteroid (organized))
    double t1, t2, pwr1, pwr2, pwr;    
    
    antennasTargetingAsteroid = getSubarray(asteroidID);
    asteroid = asteroids.at(asteroidID);
    signalStrength = asteroid.signalStrength;
    
    // build a vector from the unorganized deck of transmitted signal data sent to the asteroid at time {sendTime}
    if (antennasTargetingAsteroid.size() > 0) {
        for (int antID = 0; antID < antennasTargetingAsteroid.size(); antID++) {
            
            signalsFromAntennaX.clear();
            
            for (int j = 0; j < signalStrength.size(); j++) {
                signalForAntennaX = signalStrength.at(j);
                if (antennasTargetingAsteroid.at(antID) == signalForAntennaX.antenna_id) {
                    signalsFromAntennaX.push_back(signalForAntennaX);
                }
            }
            if (signalsFromAntennaX.size() > 1) {
                pwr = 0.0;
                for (int timeIndex = 1; timeIndex < signalsFromAntennaX.size(); timeIndex++) {
                    t1 = signalsFromAntennaX.at(timeIndex-1).arrivalTime;
                    t2 = signalsFromAntennaX.at(timeIndex).arrivalTime;
                    pwr1 = signalsFromAntennaX.at(timeIndex-1).signalPower;
                    pwr2 = signalsFromAntennaX.at(timeIndex).signalPower;
                    if ((viewTime >= t1) && (viewTime <= t2)) {
                        pwr = ((pwr2-pwr1)/(t2-t1))*viewTime+(pwr1-(((pwr2-pwr1)/(t2-t1))*t1));
                    } else {
                        pwr = 0.0;
                    }
                }
            }
            powerSum += sqrt(pwr);
        }
    }
    
    powerTransmitted = powerSum * powerSum;

    return powerTransmitted;
}
double AsteroidTracker::effectiveSignalPowerReturn(int asteroidIndex, double viewTime) {
    /*
     * Get the current signal return for an asteroid, given the number of antennas tracking it 
     */
    
    _i asteroid;
    
    asteroid = asteroids.at(asteroidIndex);
    
    int N_receive = getSubarray(asteroidIndex).size();

    double R = getAsteroidPosition(asteroidIndex, viewTime).distance;
    double R2 = R * R;
    double R4 = R2 * R2;

    double val_1 = tanh(asteroid.trajectoryInformation / Q_TRAJECTORY);
    double val_2 = T_MIN;
    
    double T_multiplier = max(val_1,val_2);
    
    asteroid = asteroids.at(asteroidIndex);
    
    // calc total quadratic power delivered to the asteroid
    double rSS = calcPowerTransmitted(asteroidIndex, viewTime);
    // rSS is the sum of transmitted power returned from the asteroid
    
    double ESPR = PEAK_GAIN
                  * asteroid.reflectivityMultiplier
                  * T_multiplier 
                  * rSS // receivedsignalStrength 
                  * (N_receive * N_receive) // needs to be calculated
                  / R4;
    
    return ESPR;
}
double AsteroidTracker::interpolateGain(double distanceBetweenAntennas, double angle) {
    /*
     * Interpolate the gain with bilinear interpolation
     */
    
    _gain gain_1;
    _gain gain_2;
    
    double PI = 3.1415926535897932384626433832795028841971693993751058;

    int len = antennaDistanceGain.size();
    double u = (len - 1.0) * angle / PI;
    int index = (int) floor(u);
    u = u - index;
    if (index == (len - 1)) {
        index = index - 1;
        u++;
    }
    
    gain_1 = antennaDistanceGain.at(index);
    gain_2 = antennaDistanceGain.at(index+1);
    double v = (distanceBetweenAntennas - ANTENNA_MIN_DISTANCE) / (ANTENNA_MAX_DISTANCE - ANTENNA_MIN_DISTANCE);
    
    double p0 = ((1.0 - v) * gain_1.minDistanceGain) + (v * gain_1.maxDistanceGain);
    double p1 = ((1.0 - v) * gain_2.minDistanceGain) + (v * gain_2.maxDistanceGain);
    
    return ((1.0 - u) * p0) + (u * p1);      
}
double AsteroidTracker::inducedNoise(int receivingAntennaIndex, double viewTime) {
    /*
     * The induced noise for a given antenna index
     */
    
    _i asteroid;
    _j antenna_t;
    _traj asteroidPos;
    
    double x_r, y_r, x_t, y_t, dX, dY;
    double D, scalarProduct, angle;
    
    double inducedNoiseSum = 0.0;
    x_r = antennas.at(receivingAntennaIndex).Xpos;
    y_r = antennas.at(receivingAntennaIndex).Ypos;
    for (int antennaID_r = 0; antennaID_r < antennas.size(); antennaID_r++) { // receiving antennas
        for (int antennaID_t = antennaID_r; antennaID_t < antennas.size(); antennaID_t++) { // transmitting antennas
            antenna_t = antennas.at(antennaID_t);
            asteroid = asteroids.at(antenna_t.asteroid_id);
            if (antenna_t.transmitting == false) {
                // do nothing if antenna_t is not transmitting
            } else {
                x_t = antenna_t.Xpos;
                y_t = antenna_t.Ypos;
                
                dX = x_t-x_r;
                dY = y_t-y_r;
                
                asteroidPos = getAsteroidPosition(asteroid.asteroidIndex, viewTime);
                
                D = sqrt((dX*dX)+(dY*dY));
                scalarProduct = (dX * asteroidPos.x) + (dY * asteroidPos.y);
                angle = acos(scalarProduct / D);
                
                inducedNoiseSum += ((antenna_t.power * interpolateGain(D, angle)) / (D*D));

            }
        }
    }
    
    return inducedNoiseSum;
}
double AsteroidTracker::calcTotalNoise(int asteroidID, double viewTime) {
    
    _j antenna;
    
    vector<int> antennasTargetingAsteroidID;
    
    antennasTargetingAsteroidID = getSubarray(asteroidID);
    int N_receive = antennasTargetingAsteroidID.size();
    double inducedNoiseSum = 0.0;
    
    for (int antennaID = 0; antennaID < antennas.size(); antennaID++) {
        antenna = antennas.at(antennaID);
        if (antenna.power == 0.0) {
            inducedNoiseSum += inducedNoise(antennaID, viewTime);
        }
    }

    double totalNoise = (SHIELDING_MULTIPLIER * inducedNoiseSum) 
        + (N_receive * BACKGROUND_NOISE);
    
    return totalNoise;
}
double AsteroidTracker::calcInformationRate(int asteroidID, double viewTime) {

    double totalGain = effectiveSignalPowerReturn(asteroidID, viewTime);
    double totalNoise = calcTotalNoise(asteroidID, viewTime);
    double val_1 = 0.0;
    double val_2 = 10 * log10(totalGain / totalNoise);
    
    double informationRate = max(val_1,val_2);
    return informationRate;
}
double AsteroidTracker::calcTrajectoryInformation(int asteroidID, double viewTime, double deltaTime) {
    _i asteroid;
    
    asteroid = asteroids.at(asteroidID);
    
    double prevTrajInfo = asteroid.trajectoryInformation;   
    double decayMultiplier = exp(-deltaTime / Q_LOST);
    double informationRate = calcInformationRate(asteroidID, viewTime);
    double trajInfo;
    
    trajInfo = 
        decayMultiplier * prevTrajInfo
        + (1.0 - decayMultiplier) * Q_LOST * informationRate;
    
    asteroid.trajectoryInformation = trajInfo;
    asteroids.at(asteroidID) = asteroid;
    
    return trajInfo;
}

// parent functions in AsteroidTracker
int AsteroidTracker::initialize(vector<double> antennaPositions, double peakGain, vector<double> minDistanceGain, vector<double> maxDistanceGain) {
    
    _j antenna;
    _gain gain;
    double dist;    
    
    numberOfAntennas = 0;

    asteroids.clear();    
    antennas.clear();
    antennaDistanceGain.clear();
    
    PEAK_GAIN = peakGain;

    for (int i = 0; i < antennaPositions.size(); i=i+2) {

        antenna.antenna_id = i / 2;
        antenna.Xpos = antennaPositions.at(i);
        antenna.Ypos = antennaPositions.at(i+1);
        antenna.t = 0.0;
        antenna.Xpoint = 0.0;
        antenna.Ypoint = 0.0;
        antenna.Zpoint = 1.0;
        antenna.power = 0.0;
        antenna.state = 0;
        antenna.asteroid_id = -1;
        antenna.transmitting = false;
        antenna.relocating = false;
        antenna.nextCommandAvailabilityTime = 0.0;
        
        antennas.push_back(antenna);     
    }
    
    for (int i = 0; i < minDistanceGain.size(); i++) {
        gain.antenna_id = i;
        gain.minDistanceGain = minDistanceGain.at(i);
        gain.maxDistanceGain = maxDistanceGain.at(i);
        antennaDistanceGain.push_back(gain);
    }
    
    for (int i = 0; i < antennas.size(); i++) {
        double minDist = 1e20;
        double maxDist = 0.0;
        for (int j = 0; j < antennas.size(); j++) {
            if (i != j) {   
                dist = distanceBetweenAntennas(i, j);
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
        //cout << "Antenna: " << i << "\tMax Dist: " << maxDist << "\tMin Dist: " << minDist << endl;
    }
 
    numberOfAntennas = antennas.size();
    
    // calculate antenna min & max distance    
    
    return _If_Nacho_Libre_Can_Wake_Up_At_Five_A_M_To_Make_The_Soup_;
}
int AsteroidTracker::asteroidAppearance(int asteroidIndex, double scienceScoreMultiplier, double reflectivityMultiplier, double initialImageInformation, double initialTrajectoryInformation, vector<double> trajectory) {
    _i asteroid;
    _j antenna;
    _traj traj;
    _signal signal;
    vector<_signal> signalStrength;
    vector<_traj> T;
   
    asteroid.asteroidIndex = asteroidIndex;
    asteroid.scienceScoreMultiplier = scienceScoreMultiplier;
    asteroid.reflectivityMultiplier = reflectivityMultiplier;
    asteroid.imageInformation = initialImageInformation;
    asteroid.trajectoryInformation = initialTrajectoryInformation;
    asteroid.trajInfoScore = 0.0;

    for (int i = 0; i < trajectory.size(); i=i+4) {
        traj.t = trajectory.at(i);
        traj.x = trajectory.at(i+1);
        traj.y = trajectory.at(i+2);
        traj.z = trajectory.at(i+3);
        traj.distance = distanceToAsteroid(traj.x, traj.y, traj.z);
        traj.travelTime = signalTravelTime(traj.distance);
        traj.ascending = false;

        
        if (traj.z > 0.0) {
            traj.visible = true;
        } else {
            traj.visible = false;
        }
        
        T.push_back(traj);
    }
    
    asteroid.timeAppear = timeOfAppearance(T);

    // cout << "Distance to Asteroid (last value):  " << T.back().distance << "\tSignal Travel Time (last value): " << T.back().travelTime << endl;
    // cout << "Asteroid ID: " << asteroid.asteroidIndex << "\tTime of Appearance: " << asteroid.timeAppear << endl;
    
    asteroid.trajectory = T;
    
    for (int antennaID = 0; antennaID < numberOfAntennas; antennaID++) {
        signal.antenna_id = antennaID;
        signal.arrivalTime = 0.0;
        signal.signalPower = 0.0;
        signalStrength.push_back(signal);
    }
    
    asteroid.signalStrength = signalStrength;
    
    _traj pos1, pos2;
    for (int i = 1; i < T.size(); i++) {
        pos1 = T.at(i-1);
        pos2 = T.at(i);
        
        if (pos2.z > pos1.z) {
            asteroid.trajectory.at(i).ascending = true;
        } else {
            asteroid.trajectory.at(i).ascending = false;
        }
    }

    asteroids.push_back(asteroid);
    numberOfAsteroids = asteroids.size(); 
  
    vector<double> targetID = create_values(0, (double) numberOfAsteroids, numberOfAntennas);
    for (int antennaID = 0; antennaID < antennas.size(); antennaID++) {
        antenna = antennas.at(antennaID);
        antenna.asteroid_id = (int) targetID.at(antennaID);
        antennas.at(antennaID) = antenna;
    }

    return _It_Is_The_Besssssssst_;
}
string AsteroidTracker::nextCommand(double currentTime) {    
    string the_big_answer;
    ostringstream powerCommand, relocateCommand;
    
    simTime = currentTime;
    
    generate_optimization(numberOfAsteroids, numberOfAntennas, 10, 100, 1000);
    /*
    double non_vis = 0;
    double vis = 0;
    
    Test to see what the % visibility of the asteroids are throughout the 
    // simulation time ... 20-42% (winter in the northern hemisphere?)
    for (double viewTime = 0.0; viewTime < SIMULATION_TIME; viewTime=viewTime+180.0) {
        for (int asteroidID = 0; asteroidID < asteroids.size(); asteroidID++) {
            if (isAsteroidVisible(asteroidID, viewTime)) {
                vis++;
                //cout << "Asteroid #: " << asteroidIndex << "\t is above the horizon at time = " << viewTime << endl;
            } else {
                non_vis++;
                //cout << "Asteroid #: " << asteroidID << "\t is NOT above the horizon at time = " << viewTime << endl;
            }
        }
    }

    // double pctVisible = 100.0 * (vis / (vis + non_vis));
     */

    // test isAntennaRelocating
    //cout << "Is antenna #3 relocating? " << isAntennaRelocating(2, 0.0) << endl;
   
    /*
    // test checkTargetProximity.  Found that asteroids are too close 86-98% of the time?!
    int tooClose = 0;   
    int notSoClose = 0;
    double viewTime = 352000;
    
    for (int i = 0; i < numberOfAsteroids; i++) {
        for (int j = i+1; j < numberOfAsteroids; j++) {
            if (checkTargetProximity(i,j,viewTime)) {
                //cout << "Asteroids #" << i << " and #" << j << " are too close at time: " << viewTime << endl;
                tooClose++;
            }
            else {
                notSoClose++;
            }
        }
    }

    double pctTooClose = 100 * ((double) tooClose / ((double) tooClose + (double) notSoClose)); */
    

    
    //double viewTime = 52000;
 
    
    /*
    // Test the function for isAntennaTransmitting
    cout << "Antennas Not Transmitting: ";
    for (int antennaID = 0; antennaID < numberOfAntennas; antennaID++) {
            bool isTransmitting = isAntennaTransmitting(antennaID, viewTime);

            if (!isTransmitting) {
                cout << " " << antennaID;
            }
    }*/
    
    /*
    int beamHit = 0;
    int no_beamHit = 0;
    int excessPower = 0;
    int no_excessPower = 0;
    
    for (int viewTime = 0; viewTime < SIMULATION_TIME; viewTime=viewTime+180) {
        antennas.at(0).power = 1.0 * MAX_TRANSMITTING_POWER;
        antennas.at(0).nextCommandAvailabilityTime = viewTime + 180; 
        antennas.at(0).asteroid_id = 0;
        
        for (int antennaID_t = 0; antennaID_t < numberOfAntennas; antennaID_t++) {
            bool beamInterference = checkAntennaBeamIntersection(antennaID_t, viewTime);
            bool powerLimitWarning = checkMaxPowerLimit(antennaID_t, viewTime);
            
            if (beamInterference) {
                beamHit++;
            } else {
                no_beamHit++;
            }
            
            if (powerLimitWarning) {
                excessPower++;
            } else {
                no_excessPower++;
            }
        }
    }
    
    double pctBeamHit = 100 * ((double) beamHit / ((double) beamHit + (double) no_beamHit));
    double pctExcessPower = 100 * ((double) excessPower / ((double) excessPower + (double) no_excessPower));
    
    double viewTime = 322000;
    double minTime = 1e20;
    double maxTime = 0;
    for (int asteroidID_1 = 0; asteroidID_1 < numberOfAsteroids; asteroidID_1++) {
        for (int asteroidID_2 = asteroidID_1+1; asteroidID_2 < numberOfAsteroids; asteroidID_2++) {
            if ((isAsteroidVisible(asteroidID_1,viewTime)) && (isAsteroidVisible(asteroidID_2,viewTime))) {
                double slewAngle = calcSlewAngle(asteroidID_1, asteroidID_2, viewTime);
                //cout << "SlewAngle: " << slewAngle << endl;
                double slewTime = slewAngle / RELOCATE_SPEED;
                if (slewTime > maxTime) {
                    maxTime = slewTime;
                }
                if (slewTime < minTime) {
                    minTime = slewTime;
                }   
            }
        }
    }*/
    
    cout << endl;

    cout << "Number of Antennas: " << numberOfAntennas << endl
            << "Number of Asteroids: " << numberOfAsteroids << endl
            << "peak Gain: " << PEAK_GAIN << endl;
            // << "Beam Hits: " << beamHit << "\t\t% Beam Interference: " << pctBeamHit << endl
            // << "Power Events: " << excessPower << "\t\t% Excess Power Events: " << pctExcessPower << endl
            // << "Min Slew Time: " << minTime << "\tMax Slew Time: " << maxTime << endl;
            // << "Too Close Events: " << tooClose << "\t% Too Close: " << pctTooClose << endl;
            // << "% Visibility: " << pctVisible << "%" << endl;
    
    // Determine next command
    bool powerCMD = false;
    
    // Generate parameters for next command
    int antennaIndex = 0;
    double powerLevel = 1.00 * MAX_TRANSMITTING_POWER;
    int newAsteroidIndex = -1;
    double deltaTime = 60.0;
    double commandTime = currentTime + deltaTime;
    
    if (powerCMD == true) {
        // "P" Change Power Command
        powerCommand << "<" << commandTime << "> <" << antennaIndex << "> P <" << powerLevel << ">";
        the_big_answer = powerCommand.str();
        cout << the_big_answer << endl;
    } else if (powerCMD == false) {
        // "R" Relocate Antenna Pointing Position        
        relocateCommand << "<" << commandTime << "> <" << antennaIndex << "> R <" << newAsteroidIndex << ">";
        the_big_answer = relocateCommand.str();
        cout << the_big_answer << endl;
    }
    
    // Transmit next command
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
    string pathToTestData = "/Users/scott/Arctria/ARC-36/data/testcase";
    string testDataFileExtension = ".txt";
    int numTestCases = 1;//10;
    
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
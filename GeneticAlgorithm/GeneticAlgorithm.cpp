/* GeneticAlgorithm [Competition Sensitive]
 * RigelFive - Hudson Ohio USA
 * 25 July 2014 - 10 August 2014
 */

#define _Genetic_Algorithm_ 0;

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

int main(int argc, char** argv) {
    GeneticAlgorithm GA;
    
    int numberOfAsteroids = 50;
    int numberOfAntennas = 50;
    int mutate = 10; // number of bits to randomly change in the genome
    int populationSize = 100; // number of genomes in the population
    int maxCycles = 1000; // number of cycles allowed with no improvement to population occurs
    
    GA.generate_optimization(numberOfAsteroids, numberOfAntennas, mutate, populationSize, maxCycles);

    return _Genetic_Algorithm_;
}
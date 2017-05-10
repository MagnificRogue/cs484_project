#include "GeneticAlgorithm.h"
#include "Constants.h"
#include <stdlib.h>     
#include <random>
#include <cmath>
#include "Simulate.cpp"
#include<iostream>
using namespace std;

/*
 * Function: initializePopulation
 * Param: Candidate* population
 * Purpose: Allocate memory for the population and initialize 
 *          members of the population with random variables
 */
void initializePopulation(Candidate *&population, bool varyNumOfLinks, bool varyWeightOfLinks, bool varyLengthOfLinks)
{
  std::default_random_engine generator;
  std::uniform_real_distribution<double> torqueDistribution(TORQUE_LOWER_BOUND,TORQUE_UPPER_BOUND);
  std::uniform_real_distribution<double> varyLinksDistribution(NUM_OF_LINKS_LOWER_BOUND,NUM_OF_LINKS_UPPER_BOUND +1 );
  std::uniform_real_distribution<double> varyWeightDistribution(WEIGHT_OF_LINK_LOWER_BOUND,WEIGHT_OF_LINK_UPPER_BOUND);
  std::uniform_real_distribution<double> varyLengthDistribution(LENGTH_OF_LINK_LOWER_BOUND,LENGTH_OF_LINK_UPPER_BOUND);

  population = (Candidate* ) malloc(POPULATION_SIZE * sizeof(*population));
   
 #pragma omp parallel for 
 for(int i=0; i < POPULATION_SIZE; i++){
    
    int numOfLinks = varyNumOfLinks ? int(varyLinksDistribution(generator)) : DEFAULT_NUM_OF_LINKS;

    population[i].numberOfLinks = numOfLinks;
    for(int j=0; j < numOfLinks; j++) {
      population[i].torque[j] = torqueDistribution(generator);
      population[i].weight[j] = varyWeightOfLinks? varyWeightDistribution(generator) : DEFAULT_WEIGHT_OF_LINK;
      population[i].length[j] = varyLengthOfLinks? varyLengthDistribution(generator) : DEFAULT_LENGTH_OF_LINK;
    }
  }
}

/*
 * Function: evaluatePopulation
 * Param: Candidate* population
 * Purpose: The purpose of this function is to assign a fitness to
 *          each member of the population dependant on their variables
 *
 */
void evaluatePopulation(Candidate *&population)
{
  #pragma omp parallel for 
  for(int i=0; i < POPULATION_SIZE; i++){
    population[i].fitness =  simulate(population[i]); 
  }
}

/*
 * Function: selectFitterIndividuals 
 * Param: Candidate* population
 * Purpose: The purpose of this function is to change population
 *          Such that the members of the population are more fit.
 *
 *          The algorithm that we have implemented is one of proportonal
 *          selection. That is, there is a random chance of selecting
 *          every member of the previous population and the probability
 *          of selecting that member is proportional to how fit they are.
 */
void selectFitterIndividuals(Candidate *&population)
{
  Candidate* newPopulation = (Candidate* ) malloc(POPULATION_SIZE * sizeof(*newPopulation));
  double totalFitness = 0;

  #pragma omp parallel for reduction(+:totalFitness)
  for(int i=0; i < POPULATION_SIZE; i++) {
    totalFitness += population[i].fitness; 
  }

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  
  #pragma omp parallel for
  for(int i=0; i < POPULATION_SIZE; i++) {
    double newPopulationTarget = distribution(generator) * totalFitness;
    for(int j=0; j < POPULATION_SIZE; j++) {
      newPopulationTarget -= population[j].fitness;

      if(newPopulationTarget <=0) {
        newPopulation[i] = population[j];
        break;
      }
    }

    if(newPopulation > 0) {
      newPopulation[i] = population[POPULATION_SIZE -1]; 
    }
  }

  free(population);
  population = newPopulation;
}

/*
 * Function: matePopulation
 * Param: Candidate* population
 * Purpose: The purpose of this function is to mate members of the population
 *          dependent on the Crossover rate.
 */
void matePopulation(Candidate *&population)
{

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  #pragma omp parallel for
  for(int i=0; i < POPULATION_SIZE; i += 2){
    if(distribution(generator) <= CROSSOVER_RATE) {
      Candidate parentA = population[i];
      Candidate parentB = population[i+1];

      int totalLength = parentA.numberOfLinks + parentB.numberOfLinks;
      int minLength = parentA.numberOfLinks < parentB.numberOfLinks ? parentA.numberOfLinks : parentB.numberOfLinks;
      int maxLength = parentA.numberOfLinks > parentB.numberOfLinks ? parentA.numberOfLinks : parentB.numberOfLinks;

      std::uniform_real_distribution<double> newNumOfLinks(minLength, maxLength + 0.1);

      Candidate childA;
      Candidate childB;

      childA.numberOfLinks = int(newNumOfLinks(generator));
      childB.numberOfLinks = int(newNumOfLinks(generator));

      double * weights = new double[totalLength];
      copy(parentA.weight, parentA.weight + parentA.numberOfLinks, weights);
      copy(parentB.weight, parentB.weight + parentB.numberOfLinks, weights + parentA.numberOfLinks);

      double * lengths = new double[totalLength];
      copy(parentA.length, parentA.length+ parentA.numberOfLinks, lengths);
      copy(parentB.length, parentB.length + parentB.numberOfLinks, lengths+ parentA.numberOfLinks);

      double * torques = new double[totalLength];
      copy(parentA.torque, parentA.torque + parentA.numberOfLinks, torques );
      copy(parentB.torque , parentB.torque + parentB.numberOfLinks, torques + parentA.numberOfLinks);

      for(int j=0; j < childA.numberOfLinks; j++){
        int torqueSeed = rand() % totalLength;
        int weightSeed = rand() % totalLength;
        int lengthSeed = rand() % totalLength;
        
        childA.torque[j] = torques[torqueSeed];
        childA.weight[j] = weights[weightSeed];
        childA.length[j] = lengths[lengthSeed];
      }

      for(int j=0; j < childB.numberOfLinks; j++){
        int torqueSeed = rand() % totalLength;
        int weightSeed = rand() % totalLength;
        int lengthSeed = rand() % totalLength;
        
        childB.torque[j] = torques[torqueSeed];
        childB.weight[j] = weights[weightSeed];
        childB.length[j] = lengths[lengthSeed];
      }
      
      population[i] = childA;
      population[i+1] = childB;
    }
  }
}

/*
 * Function: mutateIndividuals
 * Param: Candidate* population
 * Purpose: The purpose of this function is to apply a bitwise operation
 *          to mutate bits of each member of the population, just in case
 *          our fitness selection and mating accidentally removes
 *          information
 */
void mutateIndividuals(Candidate *&population)
{
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  std::uniform_real_distribution<double> torqueGenerator(TORQUE_LOWER_BOUND,TORQUE_UPPER_BOUND);

  // for each member in the population
  #pragma omp parallel for
  for(int i=0; i < POPULATION_SIZE; i++) {
    // for each random variable in the member
    for(int j=0; j < population[i].numberOfLinks-1; j++){
      if(distribution(generator) <= MUTATION_RATE) {
        population[i].torque[j] = (population[i].torque[j] + torqueGenerator(generator))/2;
      }
    }
  }
}

#include "GeneticAlgorithm.h"
#include <stdlib.h>     
#include <random>
#include <cmath>
#include "Simulate.cpp"
#include<iostream>
using namespace std;

//const double MUTATION_RATE = .1;
//const double CROSSOVER_RATE = .6;
//const int POPULATION_SIZE = 2000;
//const int TORQUE_BOUND = 100;

/*
 * Function: initializePopulation
 * Param: Candidate *&population
 * Purpose: Allocate memory for the population and initialize 
 *          members of the population with random variables
 */
void initializePopulation(Candidate *&population)
{
  population = (Candidate* ) malloc(POPULATION_SIZE * sizeof(*population));
  
  #pragma omp parallel for 
  for(int i=0; i < POPULATION_SIZE; i++){
    population[i].fitness = 0;
    for(int j=0; j< 4; j++) {
      population[i].torque[j] =  double(rand() % TORQUE_BOUND); 
    }
  }
}

/*
 * Function: evaluatePopulation
 * Param: Candidate *&population
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
 * Param: Candidate *&population
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
  

  #pragma omp parallel for private(generator)
 for(int i=0; i < POPULATION_SIZE; i++) {

    double newPopulationTarget = distribution(generator) * totalFitness;
    cout << " Population Target: " << newPopulationTarget << endl;
    for(int j=0; j < POPULATION_SIZE; j++) {
      newPopulationTarget -= population[j].fitness;

      if(newPopulationTarget <=0) {
        cout << " Found it at " << i << endl;
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
 * Param: Candidate *&population
 * Purpose: The purpose of this function is to mate members of the population
 *          dependent on the Crossover rate.
 */
void matePopulation(Candidate *&population)
{

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  #pragma omp parallel for private(generator) 
  for(int i=0; i < POPULATION_SIZE; i += 2){
    if(distribution(generator) <= CROSSOVER_RATE) {

      // for each random variable being crossed over
      for(int j=0; j < 4; j++) {
        // a bitwise crossover point, * 8 because sizeof returns how many bytes the datatype takes up
        int crossOverPoint = rand() % (sizeof(double) * 8);
        
        long lowMask = (long) pow(2,crossOverPoint) - 1;
        long highMask = lowMask ^ (-1); 

        double tmp = population[i].torque[j];
        
        population[i].torque[j] = double((long(population[i].torque[j]) & highMask) | (long(population[i+1].torque[j]) & lowMask));
        population[i+1].torque[j] = double((long(tmp) & highMask) | (long(population[i+1].torque[j]) & lowMask));
      }
    }
  }
}

/*
 * Function: mutateIndividuals
 * Param: Candidate *&population
 * Purpose: The purpose of this function is to apply a bitwise operation
 *          to mutate bits of each member of the population, just in case
 *          our fitness selection and mating accidentally removes
 *          information
 */
void mutateIndividuals(Candidate *&population)
{
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  // for each member in the population
  #pragma omp parallel for private(generator)
  for(int i=0; i < POPULATION_SIZE; i++) {
    // for each random variable in the member
    for(int j=0; j < 4; j++){
      //for each bit in that variable 
      for(int k=0; j < sizeof(double) * 8; k++) {
        if(distribution(generator) <= MUTATION_RATE){
          long mask = (long) pow(2, k);
          population[i].torque[j] = (double)(long(population[i].torque[j]) ^ mask);
        } 
      }
    }
  }
}

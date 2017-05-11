#include "Candidate.h"
#include "GeneticAlgorithm.h"
#include <iostream>
#include <random>
using namespace std;
int main(int argc, char** argv) {

  srand(123456789);
  Candidate *population;
  
  // Step 1, initialize the population
  initializePopulation(population, true, true, true);

  // Step 2, evaluate the population
  evaluatePopulation(population);

  double maxFitness = 0;
  Candidate optimalSnake;
  unsigned int iteration = 0;

  // Step 3 While Termination condition not met do
  while( iteration < 10000) { 
    iteration++;
    // Step 4 Select fitter individuals for reproduction
    selectFitterIndividuals(population);

    // Step 5 Recombine Individuals
    matePopulation(population);

    // Step 6  Mutate individuals
    mutateIndividuals(population);

    // Step 7 Evaluate the fitness of the modified individuals
    evaluatePopulation(population);

    double sum = 0;

    for(int i=0; i < POPULATION_SIZE; i++) {
      sum += population[i].fitness;
      if (population[i].fitness > maxFitness){
        maxFitness = population[i].fitness;
        optimalSnake = population[i];
      }
    }

    
  free(population);

  return 0;
}

#include "Candidate.h"
#include "GeneticAlgorithm.h"
#include <iostream>
#include <random>

int main(int argc, char** argv) {

  srand(123456789);
  Candidate *population;
  
  // Step 1, initialize the population
  initializePopulation(population);

  // Step 2, evaluate the population
  evaluatePopulation(population);

  unsigned int iteration = 0;
  // Step 3 While Termination condition not met do
  while( iteration < 10000) { 
    std::cout << "Iteration " << iteration << std::endl;
    iteration++;
    // Step 4 Select fitter individuals for reproduction
    selectFitterIndividuals(population);

    // Step 5 Recombine Individuals
    matePopulation(population);

    // Step 6  Mutate individuals
    mutateIndividuals(population);

    // Step 7 Evaluate the fitness of the modified individuals
    evaluatePopulation(population);
  }
  
  free(population);

  return 0;
}

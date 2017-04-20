#include "Candidate.h"
#include "Simulate.cpp"
#include "GeneticAlgorithm.h"
#include <iostream>

void evaluatePopulation(Candidate* population, unsigned long long &currentFitness);
void selectFitterIndividuals(Candidate* population);
void matePopulation(Candidate* population);
void mutateIndividuals(Candidate* population);

int main(int argc, char** argv) {

  srand(time(0));
  Candidate *population;
  
  // Step 1, initialize the population
  initializePopulation(population);

  // Step 2, evaluate the population
  evaluatePopulation(population);

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
  }

  return 0;
}

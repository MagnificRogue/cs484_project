#ifndef GENETICALGORITHM_H
#define GENETICALGORITHM_H

#include "Candidate.h"

const double MUTATION_RATE = .1;
const double CROSSOVER_RATE = .6;
const int POPULATION_SIZE = 2000;
const int TORQUE_BOUND = 100;

void initializePopulation(Candidate *&population);
void evaluatePopulation(Candidate *&population);
void selectFitterIndividuals(Candidate *&population);
void matePopulation(Candidate *&population);
void mutateIndividuals(Candidate *&population);

#endif

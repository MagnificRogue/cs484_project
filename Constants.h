#ifndef CONSTANTS_H
#define CONSTANTS_H

// Genetic Algorithm Constants
const double MUTATION_RATE = .1;
const double CROSSOVER_RATE = .6;
const int POPULATION_SIZE = 2000;


// Snake constants
const double TORQUE_UPPER_BOUND = 100;
const double TORQUE_LOWER_BOUND = .5;

const double pi = 3.141592653;

const int  DEFAULT_NUM_OF_LINKS = 5;
const int NUM_OF_LINKS_LOWER_BOUND = 2;
const int NUM_OF_LINKS_UPPER_BOUND = 10;

const double DEFAULT_LENGTH_OF_LINK = 1.0;
const double LENGTH_OF_LINK_LOWER_BOUND = .5; 
const double LENGTH_OF_LINK_UPPER_BOUND = 50;

const double DEFAULT_WEIGHT_OF_LINK = 1.0;
const double WEIGHT_OF_LINK_LOWER_BOUND = .5;
const double WEIGHT_OF_LINK_UPPER_BOUND = 50;

#endif

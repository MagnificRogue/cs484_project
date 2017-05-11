#ifndef CONSTANTS_H
#define CONSTANTS_H

// Genetic Algorithm Constants
// Snake constants

const double MUTATION_RATE = .2; // Rate to mutate the torque
const double CROSSOVER_RATE = .6; // Rate to mate members of the population
const double EXCHANGE_RATE = .2; // Rate at which each island sends its best snakes
const int POPULATION_SIZE = 8; // number of candidates in the population
const int TERMINATION = 2; // Number of generations to simulate
const int NUM_OF_TRAVELING_SNAKES = 2;

const double BATTLE_ROYAL_RATE = .2;

// Snake constants
const double TORQUE_UPPER_BOUND = 20;
const double TORQUE_LOWER_BOUND = 1;

const double pi = 3.141592653;

const int DEFAULT_NUM_OF_LINKS = 5;
const int NUM_OF_LINKS_LOWER_BOUND = 2;
const int NUM_OF_LINKS_UPPER_BOUND = 10;

const double DEFAULT_LENGTH_OF_LINK = 1.0;
const double LENGTH_OF_LINK_LOWER_BOUND = .5; 
const double LENGTH_OF_LINK_UPPER_BOUND = 3;

const double DEFAULT_WEIGHT_OF_LINK = 1.0;
const double WEIGHT_OF_LINK_LOWER_BOUND = .5;
const double WEIGHT_OF_LINK_UPPER_BOUND = 3;

#endif

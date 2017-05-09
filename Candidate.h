#ifndef CANDIDATE_H
#define CANDIDATE_H

#include "Constants.h"

struct Candidate 
{
  int numberOfLinks;
  double torque[NUM_OF_LINKS_UPPER_BOUND-1];
  double weight[NUM_OF_LINKS_UPPER_BOUND];
  double length[NUM_OF_LINKS_UPPER_BOUND];
  double fitness;
};




#endif

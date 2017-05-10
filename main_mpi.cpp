#include "Candidate.h"
#include "GeneticAlgorithm.h"
#include <iostream>
#include <random>
#include <mpi.h>

using namespace std;

int main(int argc, char** argv) {
  int size, rank;

  // Initialization for MPI 
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Using least significant bits of rank
  // to determine varation parameters (island model)
  unsigned bitwiseRank = rank;
  bool varyNumOfLinks = (1 & bitwiseRank);
  bool varyWeightOfLinks = (1 & bitwiseRank >> 1);
  bool varyLengthOfLinks = (1 & bitwiseRank >> 2);



  srand(time(NULL) + rank);
  Candidate *population;
  
  // Step 1, initialize the population
  initializePopulation(population, varyNumOfLinks, varyWeightOfLinks, varyLengthOfLinks);


  /* build datatype describing structure */ 

   
  MPI_Datatype MPI_CANDIDATE_TYPE; 
  MPI_Datatype type[5] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; 
  int structLength[5] = {1, NUM_OF_LINKS_UPPER_BOUND-1, NUM_OF_LINKS_UPPER_BOUND, NUM_OF_LINKS_UPPER_BOUND, 1}; 
  MPI_Aint disp[5]; 
  int base; 

   
  /* compute displacements of structure components */ 
  // Refer to here as an example: http://mpi-forum.org/docs/mpi-1.1/mpi-11-html/node61.html
  MPI_Address(population, disp); 
  MPI_Address(population[0].torque,  disp + 1); 
  MPI_Address(population[0].weight , disp + 2); 
  MPI_Address(population[0].length,  disp + 3); 
  MPI_Address(population[0].fitness, disp + 4); 

  base = disp[0]; 
  for (i=0; i < 5; i++) 
    disp[i] -= base; 

  MPI_Type_struct(5, structLength, disp, type, &MPI_CANDIDATE_TYPE); 

    

  /* Now, lets create a series of windows that can hold 32 snakes from each other process */

  MPI_Win windows[size];
  Candidate *travelingSnakes;
  MPI_Alloc_mem(32*sizeof(Candidate)*size, MPI_INFO_NULL, &travelingSnakes); 
  
  for(int i=0; i < size; i++) {
    MPI_Win_create(travelingSnakes + i, 32*sizeof(Candidate), sizeof(Candidate), MPI_INFO_NULL, MPI_COMM_WORLD, windows+i);
  }

  // Step 2, evaluate the population
  evaluatePopulation(population);

  double maxFitness = 0;
  Candidate optimalSnake;
  unsigned int iteration = 0;
  // Step 3 While Termination condition not met do
  while( iteration < 1000) { 
    iteration++;
    // Step 4 Select fitter individuals for reproduction
    selectFitterIndividuals(population);

    // Step 5 Recombine Individuals
    matePopulation(population);

    // Step 6  Mutate individuals
    mutateIndividuals(population);

    // Step 7 Evaluate the fitness of the modified individuals
    evaluatePopulation(population);

    #pragma omp parallel for 
    for(int i=0; i < POPULATION_SIZE; i++) {
      if (population[i].fitness > maxFitness){
        #pragma omp critical
        {
          maxFitness = population[i].fitness;
          optimalSnake = population[i];
        }
      }
    }

  free(population); // free our population
  free(travelingSnakes); // free our visiting friends
  
  for(int i=0; i < size; i++)
    MPI_Win_free(windows + i); // free the window

  Candidate *travelingSnakes;

	MPI_Finalize();
  return 0;
}

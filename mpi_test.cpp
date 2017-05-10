#include "Candidate.h"
#include "GeneticAlgorithm.h"
#include <iostream>
#include <random>
#include "mpi.h"
#include <assert.h>

using namespace std;

int main(int argc, char** argv) {
  int size, rank;

  // Initialization for MPI 
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  Candidate* population;

  population = new Candidate;

  Candidate ct;
  population[0] =  ct;


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
  MPI_Address(&population[0].fitness, disp + 4); 

  base = disp[0]; 
  for (int i=0; i < 5; i++) 
    disp[i] -= base; 

  MPI_Type_struct(5, structLength, disp, type, &MPI_CANDIDATE_TYPE); 
  MPI_Type_commit(&MPI_CANDIDATE_TYPE);


  /* Now, lets create a series of windows that can hold NUM_OF_TRAVELING_SNAKES snakes from each other process */

  MPI_Win win;
  Candidate *travelingSnakes;
  MPI_Alloc_mem(NUM_OF_TRAVELING_SNAKES*sizeof(Candidate)*size, MPI_INFO_NULL, &travelingSnakes); 

  MPI_Win_create(travelingSnakes, NUM_OF_TRAVELING_SNAKES*sizeof(Candidate)*size, sizeof(Candidate), MPI_INFO_NULL, MPI_COMM_WORLD, &win);


  Candidate c;

  if(rank == 0 ) {
    c.numberOfLinks = 42;
    c.fitness = 1337;
  }
 
  MPI_Bcast(&c, 1, MPI_CANDIDATE_TYPE, 0, MPI_COMM_WORLD);
  
  assert(c.numberOfLinks == 42);
  assert(c.fitness == 1337);
  
  
 //Put example, Move data from origin to target
  
  // I put a candidate in the first position of my window
  if(rank == 1) {
    c.numberOfLinks = -1;
  }
    MPI_Win_fence(0, win);

    if(rank == 1){
      c.numberOfLinks = -1;
      
//   What, how many, What type, to who, displacement, count, type, which window
     MPI_Put(&c, 1, MPI_CANDIDATE_TYPE, 0,0, 1, MPI_CANDIDATE_TYPE, win); 

  }

     MPI_Win_fence(0, win);

  if(rank==0){
       assert(travelingSnakes[0].numberOfLinks == -1);    
  } 
  
  
  free(population); // free our population
  free(travelingSnakes); // free our visiting friends
  
  MPI_Win_free(&win);// free the window

	MPI_Finalize();
  return 0;
}

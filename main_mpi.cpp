#include "Candidate.h"
#include "GeneticAlgorithm.h"
#include <iostream>
#include <random>
#include "mpi.h"
#include <algorithm>
#include <fstream>
#include <string> 

using namespace std;

int main(int argc, char** argv) {
  int size, rank;

  // Initialization for MPI 
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  // open a file with this process' rank to pipe output
  ofstream outFile;
  outFile.open(std::to_string(rank) + "-out.txt");


  // Using least significant bits of rank
  // to determine varation parameters (island model)
  unsigned bitwiseRank = rank;
  bool varyNumOfLinks = (1 & bitwiseRank);
  bool varyWeightOfLinks = (1 & bitwiseRank >> 1);
  bool varyLengthOfLinks = (1 & bitwiseRank >> 2);

  srand(time(NULL) + rank);
  Candidate *population;

  // a uniform distribution generator
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);

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
  MPI_Address(&population[0].fitness, disp + 4); 

  base = disp[0]; 
  for (int i=0; i < 5; i++) 
    disp[i] -= base; 

  MPI_Type_struct(5, structLength, disp, type, &MPI_CANDIDATE_TYPE); 
  MPI_Type_commit(&MPI_CANDIDATE_TYPE);

  /* Now, lets create a window that can hold NUM OF TRAVELING SNAKES snakes from each other process */

  MPI_Win win;
  Candidate *travelingSnakes;
  MPI_Alloc_mem(NUM_OF_TRAVELING_SNAKES*sizeof(Candidate)*size, MPI_INFO_NULL, &travelingSnakes); 
  
  //for our "flag" conditions, set our traveling snakes to have a fitness of -1
  for(int i=0; i < NUM_OF_TRAVELING_SNAKES; i++)
    travelingSnakes[i].fitness = -1;

  MPI_Win_create(travelingSnakes, NUM_OF_TRAVELING_SNAKES*sizeof(Candidate)*size, sizeof(Candidate), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

  // Step 2, evaluate the population
  evaluatePopulation(population);

  double maxFitness = 0;
  Candidate optimalSnake;
  unsigned int iteration = 0;
  // Step 3 While Termination condition not met do
  while( iteration < TERMINATION) { 
    iteration++;
    // Step 4 Select fitter individuals for reproduction
    selectFitterIndividuals(population);

    // Step 5 Recombine Individuals
    matePopulation(population);

    // Step 6  Mutate individuals
    mutateIndividuals(population);

    // Step 7 Evaluate the fitness of the modified individuals
    evaluatePopulation(population);


    // Step 8: BATTLE ROYALE!
    if(distribution(generator) <= BATTLE_ROYAL_RATE) {
      //First, allocate space to hold NUM_OF_TRAVELING_SNAKES 
      Candidate toSend[NUM_OF_TRAVELING_SNAKES]; 
      //Second, pick a random number of processes that this
      // process is going to send to
      int numberToSend = (rand() % 3) + 1;
      
      // now pick out the ranks
      int ranks[numberToSend];
      
      for(int i=0; i < numberToSend; i++){
        int sendId = rand() % size;

        while(sendId == rank)
          sendId = rand() % size;

        ranks[i] = sendId;
      }
      
      //Now, fill the snakes we're sending with our best and brightest
      // first, sort our array
      sort(population, population + POPULATION_SIZE, [](const Candidate &a, const Candidate &b) {
        return a.fitness < b.fitness;
      });
      
      //second, fill our our toSend array
      for(int i=0; i < NUM_OF_TRAVELING_SNAKES; i++)
        toSend[i] = population[i];

      // third, fill the windows at my rank * NUM_OF_TRAVELING_SNAKES offset
      MPI_Win_fence(0, win);
      
      for(int i=0; i < numberToSend; i++) {
        MPI_Put(toSend, NUM_OF_TRAVELING_SNAKES, MPI_CANDIDATE_TYPE, ranks[i], rank*NUM_OF_TRAVELING_SNAKES, NUM_OF_TRAVELING_SNAKES, MPI_CANDIDATE_TYPE, win); 
      }
     MPI_Win_fence(0, win);


     // and, lastly, lets randomize our population as to not affect future generations
     std::random_shuffle(population, population + POPULATION_SIZE);
    }

    // now, conversely, lets  check to see if we have any snakes for ourself and battle royal them!
    MPI_Win_fence(0, win);
    for(int i=0; i < size * NUM_OF_TRAVELING_SNAKES; i++){

      // ughoh! we have a visitor
      if(travelingSnakes[i].fitness != -1){
        Candidate visitor = travelingSnakes[i]; // grab that snake
        travelingSnakes[i].fitness = -1; // reset this position
        
        int targetIndex = rand() % POPULATION_SIZE;
        
        // if the visitor is fitter, kill that unlucky snake
        if(visitor.fitness > population[targetIndex].fitness)
          population[targetIndex] = visitor;
      }
    } 


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

    outFile << "Iteration: " << iteration << "\n";
    outFile << "Fitness: " << maxFitness;
    outFile << "\t Number of links: " << optimalSnake.numberOfLinks << "\n";
    for(int i=0; i < optimalSnake.numberOfLinks; i++){
      outFile << "\t\tTorque[" << i << "]: " << optimalSnake.torque[i] << "\n"; 
      outFile << "\t\tWeight[" << i << "]: " << optimalSnake.weight[i] << "\n";
      outFile << "\t\tLength[" << i << "]: " << optimalSnake.length[i] << "\n";
    }
    outFile << "\n";
  }

  free(population); // free our population
  MPI_Win_free(&win); // free the window
  free(travelingSnakes); // free our visiting friends
  
  outFile.close();

	MPI_Finalize();
  return 0;
}

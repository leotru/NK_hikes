//****************************************************************************
//
// nk_walk.cpp version 3.2  2021-02-04:
//
// Main routine to simulate adative walks on the Stuart Kauffman's
// NK-fitness landscape model, with point mutations (as usual) and 
// inversions (as new).
//
// Help: run the command line ./nk_walk -help  to print out the required 
// input parameters. 
//
// Why reinvent the wheel? 
// This code uses the one developed by Wim Hordijk (in its version of 
// August 23, 2010): "which uses an idea and some code from Terry Jones 
// to calculate the fitness values. Instead of storing a table with 
// N*(2^K) fitness contributions, a unique (long) integer is deterministically 
// calculated for each gene and every possible neighborhood configuration. 
// This unique integer is then used as the seed for drawing a random number
// which is the fitness contribution of that gene given tha neighborhood 
// configuration."
// 
// Our nk_walk.cpp history:
// v.0.0 2010-08-23: original code by Wim Hordijk.
// v.1.0 2019-09-12: first versions with the inversions mutations operation.
//                   Done by Paul Banse (paul.banse@inria.fr) and
//                   Leo Trujillo (leonardo.trujillo@inria.fr)
// v.1.1 2020-02-28: final version to perfom exhaustive exploratory expirical 
//                   analysis. Assembled to the "AWNK numerical project", with
//                   the scripts to run the simulations and data analysis.
//                   Done by Leo.
// v.2.0 2020-08-06: calibrations of the numerical experiments and benchmarks. 
//                   Exhaustive empirical analysis (first report 2020-10-15).
// v.3.0 2020-11-24: some modifications and the correction of the bias towards 
//                   inversions of size two. Done by
//                   Guillaume Beslon (guillaume.beslon@inria.fr)
// v.3.1 2021-01-15: benchmarks and test of new parameters. Reorganization
//                   of the code and standarization with the other components
//                   of  AWNK. Algorithms + Data Structures = Programs. So this
//                   version was incorporated in the program project 
//                   simNK v.0.1
// v.3.2 2021-02-04: implemnts the option to read an initial genome from a 
//                   file containing a genome string of size N (for a binary 
//                   alphabet A = {0,1}).
//****************************************************************************
// PROGRAM PARAMETERS (Global variables):
//
// Command-line option variables:
// N:           Length of the genomes.
// K:           Number of epistatic interactions.
// A:           Alphabet size.
//              Remark: the inversion operation implemented in the present 
//                      version of the code are consistent for a binary
//                      alphabet A = {0,1}. 
// epi:         Type of epistatic interactions (adjacent or random).
// nk_seed:     Seed value for the random number generator of the landscape.
// wlk_seed:    Seed value for the random number generator of the random walk.
// M:           Fraction of inversions:
//                                    M = 0: point mutations, 
//                                    M = 1: inversions, 
//                                    0.0 < M < 1.0: both types of mutations.
// T:           Number of steps to perform in the adaptive walk.
// wikde_type:  Initial genome to be read from a given file containing a genome
//              string of size N (for a binary alphabet A = {0,1}).
//
// Pointers:
// nk:       A pointer to a NK-landscape.
// rnd:      A pointer to a random number generator.
//****************************************************************************
#include "NK.h"
#include "Random.h"
#include <string.h>
#include <stdio.h>
#include <iostream>
using namespace std;

// ---
// Global variables and pointers.
// ---
int N, K, A, epi, nk_seed, wlk_seed, T;
double M;

static char*  wild_type;

Random  *rnd;
NK_Landscape *nk;

// ---
// Function prototypes.
// ---
int GetArguments (int argc, char **argv);
double AdaptiveWalk(char* wild_type);

//****************************************************************************
// main: Main routine of the program.
// 
// Parameters:
//   - argc: The number of arguments to the program.
//   - argv: A pointer to the list of arguments.
// 
// Returns:
//   If everything went fine: 0,
//   otherwise:               1.
//****************************************************************************
int main (int argc, char **argv)
{
  int status;
  double fitness;

  FILE * summary;

  status = 0;

  // ---
  // Get and check the arguments.
  // ---
  if (GetArguments (argc, argv) == -1)
  {
    status = 1;
    goto End_of_Routine;
  }

  // ---
  // Create an NK-landscape.
  // ---
  nk = new NK_Landscape (N, K, epi, A, nk_seed);
  //nk -> Test();
  if (!nk -> init_OK)
  {
    status = 1;
    cerr << "Could not create NK-landscape." << endl;
    goto End_of_Routine;
  }

  // ---
  // Create a random number generator for the walk.
  // ---
  rnd = new Random();
  rnd -> SetSeed(wlk_seed);

  // ---
  // Perform a random walk.
  // ---
  //if (AdaptiveWalk(wild_type) == -1)
  //{
  //  status = 1;
  //  goto End_of_Routine;
  //}
  //else{
  fitness = AdaptiveWalk(wild_type);
  //}
  
  // ---  
  // Clean up after use...
  // ---
  delete nk;
  delete rnd;

  // ---
  // Save the final fitness to comma-separated values file.
  // Format: N K M fitness
  // ---
  //summary = fopen("../../experiment_results.csv","a");
  summary = fopen("final_fitness.csv","a");
  fprintf(summary,"%d;%d;%lf;%.12lf \n",N,K,M,fitness);
  fclose(summary);
  
  End_of_Routine:

  // ---
  // Return the status.
  // ---
  return (status);
}

//****************************************************************************
// GetArguments: Get and check the command line arguments.
//
// Parameters:
//   - argc: The number of arguments to the program.
//   - argv: A pointer to the list of arguments.
// 
// Returns:
//   If everything went fine:  0.
//   Otherwise:               -1.
//****************************************************************************
int GetArguments (int argc, char **argv)
{
  int status, i;

  status = 0;

  // ---
  // Set defaults.
  // ---
  N = -1;
  K = -1;
  A = 2;
  M = 0;
  epi = NK_Landscape::ADJ;
  nk_seed = -1;
  wlk_seed = -1;
  T = 10000;
  wild_type = nullptr;

  // ---
  // Get and check all arguments.
  // ---
  i = 1;
  while (i < argc)
  {
    if (strcmp (argv[i], "-n") == 0)
    {
      if ((sscanf (argv[++i], "%d", &N) != 1) || (N < 2))
      {
	      status = -1;
	      cerr << "Invalid value for N: " <<  argv[i]  <<  endl;
	      goto End_of_Routine;
      }
      i++;
    }
    else if (strcmp (argv[i], "-k") == 0)
    {
      if ((sscanf (argv[++i], "%d", &K) != 1) || (K < 0) || (K >= N))
      {
	      status = -1;
	      cerr << "Invalid value for K: " << argv[i] << endl;
	      goto End_of_Routine;
      }
      i++;
    }
    else if (strcmp (argv[i], "-a") == 0)
    {
      if ((sscanf (argv[++i], "%d", &A) != 1) || (A < 2))
      {
	      status = -1;
	      cerr << "Invalid value for A: " << argv[i] << endl;
	      goto End_of_Routine;
      }
      i++;
    }
    else if (strcmp (argv[i], "-m") == 0)
    {
      if ((sscanf (argv[++i], "%lf", &M) < 0) || (M > 1))
      {
	      status = -1;
	      cerr << "Invalid value for M: " << argv[i] << endl;
	      goto End_of_Routine;
      }
      i++;
    }
    else if (strcmp (argv[i], "-e") == 0)
    {
      if (strcmp (argv[++i], "ADJ") == 0)
      {
	      epi = NK_Landscape::ADJ;
	      //printf("adj (%d)\n",epi);
      }
      else if (strcmp (argv[i], "RND") == 0)
      {
	      epi = NK_Landscape::RND;
	      //printf("rnd (%d)\n",epi);
      }
      else
      {
	      status = -1;
	      cerr << "Unknown type of epistatic interactions: " << argv[i] << endl;
	      goto End_of_Routine;
      }
      i++;
    }
    else if (strcmp (argv[i], "-snk") == 0)
    {
      if (sscanf (argv[++i], "%d", &nk_seed) != 1)
      {
	      status = -1;
	      cerr << "Invalid value for nk_seed: " << argv[i] << endl;
	      goto End_of_Routine;
      }
      i++;
    }
    else if (strcmp (argv[i], "-swlk") == 0)
    {
      if (sscanf (argv[++i], "%d", &wlk_seed) != 1)
      {
	      status = -1;
	      cerr << "Invalid value for wlk_seed: " << argv[i] << endl;
	      goto End_of_Routine;
      }
      i++;
    }
    //-------------------
    // NEW 
    // Revisar ya que cuando los nombres no coinciden da: Segmentation fault: 11
    else if (strcmp (argv[i], "-f") == 0)
    {
      wild_type = new char[strlen(argv[++i]) + 1];
      if ((sscanf (argv[i], "%s", wild_type) != 1))
      {
	      status = -1;
	      cerr << "Invalid file name: " << argv[i] << endl;
	      goto End_of_Routine;
      }
      i++;
    }
    //-------------------
    else if (strcmp (argv[i], "-t") == 0)
    {
      if ((sscanf (argv[++i], "%d", &T) != 1) || (T < 1))
      {
	      status = -1;
	      cerr << "Invalid value for T: " << argv[i] << endl;
	      goto End_of_Routine;
      }
      i++;
    }
    else if (strcmp (argv[i], "-help") == 0)
    {
      cout //<< endl
      << "-----------------------------------------------------------------------" << endl
      << argv[0] 
      << " -n <N> -k <K> [-a <A>] [-e <epi>] [-snk <nk_seed>]" << endl
      << " [-swlk <wlk_seed>] [-t <T>] [-m <M>] [-f <wild_type>] [-help]" << endl
	    << endl
	    << " N:           Length of the genomes." << endl
	    << " K:           Number of epistatic interactions (0 <= K < N)." << endl
	    //<< endl
      << " A:           Alphabet size (default = 2)." << endl
	    << " epi:         Type of epistatic interactions: " << endl
      << "                 ADJ: adjacent interactions (default),"<< endl
      << "                 RND: random interactions." << endl
	    << " nk_seed:     PRNG's seed for the landscape (default = -1)." << endl
	    << " wlk_seed:    PRNG's seed for the random walk (default = -1)." << endl
	    << " T:           Number of steps in the random walk (default = 10000)." << endl
	    << " M:           Fraction of inversions: "<< endl 
      << "                 M = 0: point mutations, "<< endl 
      << "                 M = 1: inversions, "<< endl 
      << "                 0.0 < M < 1.0: both types of mutations."<< endl 
      //<< "              Inversions of length 1 are allowed (point mutations)."<< endl
	    //<< endl
      << " wild_type:  (a file name) If the initial genome is to be read from " << endl
      << "              a file, otherwise a random genome is generated." << endl
      << endl
	    << " help:        Print out this message and exit." << endl
      << "-----------------------------------------------------------------------" << endl;
      
      status = -1;
      goto End_of_Routine;
    }
    else
    {
      status = -1;
      cerr << "Unknow option " << argv[i] << endl;
      goto End_of_Routine;
    }
  }

  // ---
  // Make sure the user has set at least the N and K values.
  // ---
  if ((N < 0) || (K < 0))
  {
    status = -1;
    cerr << "Expecting at least the -n and -k options..." << endl;
    goto End_of_Routine;
  }

  End_of_Routine:

  // ---
  // Return the status.
  // ---
  return (status);
}

//****************************************************************************
// AdaptiveWalk: Perform an adaptive walk on the landscape and print out the
//               generated time series of fitness values.
//****************************************************************************
//
//--------------------------------------
// The adaptive walk algoritm:
//--------------------------------------
// Input: genome x and fx <-- Fitness(x)
//    y <-- Mutate(x)
//    fy <-- Fitness(y)
//    if fx < fy
//        x <-- y
//        fx <-- fy
//    end if
// Output: genome x and fitness fx
//-------------------------------------- 
//
//****************************************************************************
double AdaptiveWalk(char* wild_type)
{
  int status, i, j, l, loc, ini_loc, end_loc, point_loc, val, mut_type;
  int *genome, *x, *y, *comp;
  double fx, fy, mut;

  status = 0;
    
  genome = new int[N];
  x = new int[N];
  y = new int[N];
  comp = new int[N];
    
  // ---
  // Start by creating a random genome.
  // ---
  //for (i = 0; i < N; i++)  genome[i] = rnd -> Unif(A);
  //for (i = 0; i < N; i++)  genome[i] = 0;
  //-------------------
  // Start by reading a given genome (for example: genome.csv)
  if (wild_type != nullptr)
  {

  FILE* genome_file = fopen(wild_type, "r");

      int i =0;
      do 
      {
        char c = fgetc (genome_file);
        genome[i] = (int) (c - '0');
        i++;
      } while (i != N);
  }
  else{
    for (i = 0; i < N; i++)  genome[i] = rnd -> Unif(A);
  }
  //-------------------   
  
  for (loc = 0; loc < N; loc++) x[loc] = genome[loc];
  fx = nk -> Fitness(x);

  // ---
  // Print in terminal.
  // ---
  //sprintf("genome;generation;fitness;ini_loc;end_loc;type\n");
  // To save the actual genotype uncomment the next line.
  for (loc = 0; loc < N; loc++) cout << x[loc]; cout << ";";
  printf("%d;%.12lf;%d;%d;%d\n",0,fx,-1,-1,0);
    
  // ---
  // Perform T-1 mutation/selection steps.
  // ---
  for (i = 1; i < T; i++) 
  {
    
    // ---
    // Mutation.
    // ---
    mut = rnd -> Unif(100000000) / 100000000.0;
      
    if (mut > M)
	    mut_type = 0;
    else
	    mut_type = 1;
    
    if (mut_type == 1)
    {
      // ---
      // Mutate the genome with inversions.
      // ---
      ini_loc = rnd -> Unif(N);
      end_loc = rnd -> Unif(N);

      if (ini_loc < end_loc)
      {
        j = end_loc;
        for(loc = ini_loc; loc <= end_loc; loc ++)
        {
          comp[loc] = (genome[j] + 1) % 2;
          j--;
        }
        for(loc = ini_loc; loc <= end_loc; loc ++) genome[loc] = comp[loc];
      }
      else{
        j = end_loc + N;
        for(loc = ini_loc; loc <= end_loc + N; loc ++)
        {
          l=(loc % N);
          comp[l] = (genome[j % N] + 1) % 2;
          j--;
        }
        for(loc = ini_loc; loc <= end_loc + N; loc ++) genome[loc % N] = comp[loc % N];
      }
      
      // ---
      // Calculate the fitness of the inversion mutation.
      // ---
      for (loc = 0; loc < N; loc++) y[loc] = genome[loc];
      fy = nk -> Fitness(y);
      }
    else
    {
      // ---
      // Mutate the genome with point mutations.
      // ---
      loc = rnd -> Unif(N);
      point_loc = loc;
      val = genome[loc];
      while (val == genome[loc])
      {
        val = rnd -> Unif(A); 
      }
      genome[loc] = val;
      
      // ---
      // Calculate the fitness of the point mutation.
      // ---
      for (loc = 0; loc < N; loc++) y[loc] = genome[loc];
      fy = nk -> Fitness(y);

    }

    // ---
    // Selection.
    // ---
    if(fy > fx)
    {
      for (loc = 0; loc < N; loc++) x[loc] = y[loc];
      for (loc = 0; loc < N; loc++) genome[loc] = y[loc];
      fx = fy;

		  // ---
      // Print in terminal when the genotype is selected.
      // ---
      // To save the actual genotype uncomment the next line.
      for (loc = 0; loc < N; loc++) cout << x[loc];cout << ";";
		  printf("%d;%.12lf", i, fx);
		  
      if(mut_type == 1)
      {
		    printf(";%d;%d;%d\n", ini_loc, end_loc,mut_type);
		  }
		  else
      {
		    printf(";%d;%d;%d\n", point_loc, point_loc,mut_type);
		  }
    }
    else
    {
      for (loc = 0; loc < N; loc++) genome[loc] = x[loc];
	    // To save all the data per time-step, uncomment the next line.
      //for (loc = 0; loc < N; loc++) cout << x[loc];cout << "  ";
	    //printf("%d;%.12lf",i,fx);
		  //printf(";-1;-1;0\n");
    }
  }

  End_of_Routine:

  // ---
  // Return the final fitness value.
  // ---
  return (fx);
}

//****************************************************************************
// Bye!
//****************************************************************************

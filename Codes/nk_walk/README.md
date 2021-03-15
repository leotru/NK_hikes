 nk_walk.cpp version 3.2  2021-02-04:

 Main routine to simulate adative walks on the Stuart Kauffman's
 NK-fitness landscape model, with point mutations (as usual) and 
 inversions (as new).

 Help: run the command line ./nk_walk -help  to print out the required 
 input parameters. 

 Why reinvent the wheel? 
 This code uses the one developed by Wim Hordijk (in its version of 
 August 23, 2010): "which uses an idea and some code from Terry Jones 
 to calculate the fitness values. Instead of storing a table with 
 N*(2^K) fitness contributions, a unique (long) integer is deterministically 
 calculated for each gene and every possible neighborhood configuration. 
 This unique integer is then used as the seed for drawing a random number
 which is the fitness contribution of that gene given tha neighborhood 
 configuration."
 
http://www.cs.unibo.it/~fioretti/CODE/NK/ 
https://github.com/terrycojones/nk-landscapes)


 Our nk_walk.cpp history:
 v.0.0 2010-08-23: original code by Wim Hordijk.
 v.1.0 2019-09-12: first versions with the inversions mutations operation.
                   Done by Paul Banse (paul.banse@inria.fr) and
                   Leo Trujillo (leonardo.trujillo@inria.fr)
 v.1.1 2020-02-28: final version to perfom exhaustive exploratory expirical 
                   analysis. Assembled to the "AWNK numerical project", with
                   the scripts to run the simulations and data analysis.
                   Done by Leo.
 v.2.0 2020-08-06: calibrations of the numerical experiments and benchmarks. 
                   Exhaustive empirical analysis (first report 2020-10-15).
 v.3.0 2020-11-24: some modifications and the correction of the bias towards 
                   inversions of size two. Done by
                   Guillaume Beslon (guillaume.beslon@inria.fr)
 v.3.1 2021-01-15: benchmarks and test of new parameters. Reorganization
                   of the code and standarization with the other components
                   of  AWNK. Algorithms + Data Structures = Programs. So this
                   version was incorporated in the program project 
                   simNK v.0.1
 v.3.2 2021-02-04: implemnts the option to read an initial genome from a 
                   file containing a genome string of size N (for a binary 
                   alphabet A = {0,1}).

 PROGRAM PARAMETERS (Global variables):

 Command-line option variables:
 N:           Length of the genomes.
 K:           Number of epistatic interactions.
 A:           Alphabet size.
              Remark: the inversion operation implemented in the present 
                      version of the code are consistent for a binary
                      alphabet A = {0,1}. 
 epi:         Type of epistatic interactions (adjacent or random).
 nk_seed:     Seed value for the random number generator of the landscape.
 wlk_seed:    Seed value for the random number generator of the random walk.
 M:           Fraction of inversions:
                                    M = 0: point mutations, 
                                    M = 1: inversions, 
                                    0.0 < M < 1.0: both types of mutations.
 T:           Number of steps to perform in the adaptive walk.
 wikde_type:  Initial genome to be read from a given file containing a genome
              string of size N (for a binary alphabet A = {0,1}).

 Pointers:
 nk:       A pointer to a NK-landscape.
 rnd:      A pointer to a random number generator.


:.
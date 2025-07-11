/*******************************************************
 * This file is intended to test the implementation
 * of the Sorting by Reversals method.
 * 
 * It creates random permutations of a given size,
 * specified by the user in the input, and it checks
 * if the solutions are correct. 
 * 
 * Specifically, it checks if a sequence of reversals
 * provided by the method transforms the input permutation
 * into the identity permutation and if the number of 
 * reversals corresponds to the minimum number expected.
 * 
 * In addition to check the correctness of the method, it also 
 * saves many stats about the permutations (number of breakpoints,
 * number of cycles, etc.) and the time needed to run the method.
 * 
 * Compile:
 * g++ tests_RandomPermutations.cpp sortByReversals.cpp findComponents_Bader2001.cpp sortOrientedByReversals_Tannier2007.cpp solveUnoriented_HannenhalliPevzner1999.cpp genome.cpp -o gensort_random
 * 
 * Run:
 * ./gensort_random 42 10 3 0 1
 * 
 *******************************************************/

#include <iostream>
#include <vector>
#include <cstdlib>   // exit
#include <cmath>	 // abs
#include <chrono>

#include "sortByReversals.hpp"

/*******************************************************
 * Some basic tests.
*******************************************************/

SortByReversals testCase_generalSort(GenomeMultichrom<int>& genome_A, GenomeMultichrom<int>& genome_B, std::mt19937& rng, const bool debug, const bool printSol){
	
	SortByReversals sortGenome(genome_A,genome_B,debug);
	
	// Sort genome and measure time.
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	sortGenome.sort(rng);
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	// [Check 1] Check if the sorting scenario ends with the identity permutation (i.e. no breakpoints).
	bool correctSolution = (printSol ? sortGenome.printSolution() : sortGenome.checkSolution());
	if(!correctSolution){
		std::cout << "ERROR! Problem during the sorting. Program is aborting." << std::endl;
		exit(1);
	}
	// [Check 2] Check if the number of reversals is minimal.
	correctSolution = sortGenome.printStats();
	if(!correctSolution){
		std::cout << "ERROR! Problem during the sorting. Program is aborting." << std::endl;
		exit(1);
	}
	std::cout << "- Running time: " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin)).count() << " [µs]." << std::endl;	
	return sortGenome;
}

int main(int argc, char* argv[]) {
	
	if (argc < 2) {
		std::cout << "ERROR! No extra command line argument passed other than program name. Program is aborting." << std::endl;
		exit(1);
	} else if (argc != 6) {
		std::cout << "ERROR! Wrong number of arguments. Expected: ./program [seed] [#genes] [#tests] [debug] [print solution]. Program is aborting." << std::endl;
		exit(1);
	}

	// Parameters.
	int const seed      = std::stoi(argv[1]);        // input (command line argument): seed random number generator.
	int const nbgenes   = std::stoi(argv[2]);        // input (command line argument): number of genes.
	int const nbtests   = std::stoi(argv[3]);        // input (command line argument): how many random permutations will be generated.
	bool const debug    = (std::stoi(argv[4]) == 1); // input (command line argument): debug flag (1: debug is ON; 0 (or any other integer for now): debug is OFF).
	bool const printSol = (std::stoi(argv[5]) == 1); // input (command line argument): flag for printing solution (1: solution is going to be printed; 0 (or any other integer for now): solution is NOT printed).
	int const nbchrom   = 1; // For now only unichromosomal genomes.
	

	std::cout << "**********************************************" << std::endl;
	std::cout << "* Sort by reversals with random permutations *" << std::endl;
	std::cout << "**********************************************" << std::endl << std::endl;
	std::cout << "- Seed to reproduce tests: "    << seed    << std::endl;
	std::cout << "- Genome size (nb. of genes): " << nbgenes << std::endl;
	std::cout << "- Nb. random permutations: "    << nbtests << std::endl;

	// Create a random number generator
	std::mt19937 rng(seed); // Seed the generator
	// Distribution for the probability of negative signs in genes.
	std::uniform_real_distribution<> distr(0.0, 1.0);

	// Generate random tests.
	for(int i=0; i<nbtests; ++i){

		// Probability for a gene to have negative sign.
		double const probRev = distr(rng);
		std::cout << std::endl << std::endl << "Test " << (i+1) << ": Random genome (n=" << nbgenes << ";p_rev=" << probRev << ")" << std::endl;

		// Create random genomes using the random constructor.
		GenomeMultichrom<int> genome_A(rng, nbgenes, nbchrom, probRev);
		GenomeMultichrom<int> genome_B(rng, nbgenes, nbchrom, probRev, genome_A.gene_labels_map);

		// Sort genomes using the minimum number of reversals.
		SortByReversals solution = testCase_generalSort(genome_A, genome_B, rng, debug, printSol);

		// Save stats.
		// Save genomes.
	} 
	
	std::cout << std::endl << "Bye bye" << std::endl;
	return 0;
}

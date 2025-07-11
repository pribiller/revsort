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

SortByReversals testCase_generalSort(GenomeMultichrom<int>& genome_A, GenomeMultichrom<int>& genome_B, std::mt19937& rng, int verbose){
	
	const bool debug = (verbose > 2);
	SortByReversals sortGenome(genome_A,genome_B,debug);

	// Sort genome and measure time.
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	sortGenome.sort(rng);
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	// [Check 1] Check if the sorting scenario ends with the identity permutation (i.e. no breakpoints).
	bool correctSolution = ((verbose > 1) ? sortGenome.printSolution() : sortGenome.checkSolution());
	if(!correctSolution){
		std::cout << "ERROR! Problem during the sorting. Program is aborting." << std::endl;
		exit(1);
	}
	// [Check 2] Check if the number of reversals is minimal.
	correctSolution = ((verbose > 0) ? sortGenome.printStats() : sortGenome.checkDistance());
	if(verbose > 0){std::cout << "- Running time: " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin)).count() << " [µs]." << std::endl;}
	if(!correctSolution){
		std::cout << "ERROR! Problem during the sorting. Program is aborting." << std::endl;
		exit(1);
	}
	return sortGenome;
}

void printHelper() {
	std::cout << "Use: ./gensort_random [seed] [genes] [tests] [verbose]" << std::endl;
	std::cout << "Example: ./gensort_random 42 10 3 2" << std::endl;
	std::cout << "Parameters:" << std::endl;
	std::cout << "- [seed]: any positive integer number greater than 0;" << std::endl;
	std::cout << "- [genes]: the number of genes. In the example above, all random genomes will have 1 chromosome and 10 genes;" << std::endl;
	std::cout << "- [tests]: the number of tests. Each test will create a pair of random genomes and sort them. In the example above, 3 scenarios will be randomly generated and solved;" << std::endl;
	std::cout << "- [verbose]: a number between 0 and 4 that controls the verbosity level for the output messages." << std::endl;
	std::cout << "  - [verbose=0][for cluster]: only basic stats are printed;" << std::endl;
	std::cout << "  - [verbose=1][for user]: basic messages;" << std::endl;
	std::cout << "  - [verbose=2][for user]: basic messages + replay of the sorting scenario (not recommended if the genome is too big);" << std::endl;
	std::cout << "  - [verbose=3][for developer]: very detailed messages, useful for debugging;" << std::endl;
}

int main(int argc, char* argv[]) {
	
	if (argc < 2) {
		std::cout << "ERROR! No extra command line argument passed other than program name. Program is aborting." << std::endl;
		printHelper();
		exit(1);
	} else if (argc != 5) {
		std::cout << "ERROR! Wrong number of arguments. Program is aborting." << std::endl;
		printHelper();
		exit(1);
	}

	// Parameters.
	int const seed    = std::stoi(argv[1]); // input (command line argument): seed random number generator.
	int const nbgenes = std::stoi(argv[2]); // input (command line argument): number of genes.
	int const nbtests = std::stoi(argv[3]); // input (command line argument): how many random permutations will be generated.
	int const verbose = std::stoi(argv[4]); // input (command line argument): verbose mode (0: [cluster] only basic stats are printed; 1: [user] basic messages; 2: [user] basic messages+replay of the sorting scenario (not recommended if the genome is too big); 3: [developer] very detailed debugging messages.
	int const nbchrom   = 1; // For now only unichromosomal genomes.
	
	if(verbose > 0) {
		std::cout << "**********************************************" << std::endl;
		std::cout << "* Sort by reversals with random permutations *" << std::endl;
		std::cout << "**********************************************" << std::endl << std::endl;
		std::cout << "- Seed to reproduce tests: "    << seed    << std::endl;
		std::cout << "- Genome size (nb. of genes): " << nbgenes << std::endl;
		std::cout << "- Nb. random permutations: "    << nbtests << std::endl;
	}

	// Create a random number generator
	std::mt19937 rng(seed); // Seed the generator
	// Distribution for the probability of negative signs in genes.
	std::uniform_real_distribution<> distr(0.0, 1.0);

	// Generate random tests.
	for(int i=0; i<nbtests; ++i){

		// Probability for a gene to have negative sign.
		double const probRev = distr(rng);
		if(verbose > 0) {
			std::cout << std::endl << std::endl << "Test " << (i+1) << ": Random genome (n=" << nbgenes << ";p_rev=" << probRev << ")" << std::endl;
		}

		// Create random genomes using the random constructor.
		GenomeMultichrom<int> genome_A(rng, nbgenes, nbchrom, probRev);
		GenomeMultichrom<int> genome_B(rng, nbgenes, nbchrom, probRev, genome_A.gene_labels_map);

		// Sort genomes using the minimum number of reversals.
		SortByReversals solution = testCase_generalSort(genome_A, genome_B, rng, verbose);

		// Save stats.
		// Save genomes.
	} 
	
	std::cout << std::endl << "Bye bye" << std::endl;
	return 0;
}

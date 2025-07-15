/*******************************************************
 * Original author: Priscila Biller
 * Created: July/2025
 * License: GPL v3
 * 
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
 * g++ tests_RandomPermutations.cpp sortByReversals.cpp findComponents_Bader2001.cpp sortOrientedByReversals_Tannier2007.cpp solveUnoriented_HannenhalliPevzner1999.cpp genome.cpp -o revsort_random
 * 
 * Run:
 * ./revsort_random 42 10 4 3 2
 * 
 *******************************************************/

#include <iostream>
#include <vector>
#include <cstdlib>   // exit
#include <cmath>	 // abs

#include "sortByReversals.hpp"

/*******************************************************
 * Some basic tests.
*******************************************************/

SortByReversals testCase_generalSort(GenomeMultichrom<int>& genome_A, GenomeMultichrom<int>& genome_B, std::mt19937& rng, int verbose){
	const bool debug = (verbose > 2);
	SortByReversals sortGenome(genome_A,genome_B,debug);
	// Sort genome and measure time.
	sortGenome.sort(rng);
	// [Check 1] Check if the sorting scenario ends with the identity permutation (i.e. no breakpoints).
	bool correctSolution = ((verbose > 1) ? sortGenome.printSolution() : sortGenome.checkSolution());
	if(!correctSolution){
		std::cout << "ERROR! Problem during the sorting. Program is aborting." << std::endl;
		exit(1);
	}
	// [Check 2] Check if the number of reversals is minimal.
	correctSolution = ((verbose > 0) ? sortGenome.printStats() : sortGenome.checkDistance());
	if(!correctSolution){
		std::cout << "ERROR! Problem during the sorting. Program is aborting." << std::endl;
		exit(1);
	}
	return sortGenome;
}

GenomeMultichrom<int> createRandomGenome(GenomeMultichrom<int>& genome_A, const int nb_reversals, const double probRev, std::mt19937& rng){
	if(nb_reversals > 0){
		GenomePermutation<BlockSimple> genome_B(genome_A.getExtendedGenome());
		applyRandomReversals(genome_B, nb_reversals, rng);
		// Convert labels in the random genome.
		std::vector<int> genome_B_perm = genome_B.getUnextendedPerm();
		std::vector<int> genome_B_labels;
		std::vector<bool> genome_B_signs;
		genome_B_labels.reserve(genome_B_perm.size());
		genome_B_signs.reserve(genome_B_perm.size());
		// std::cout << " - Genome B: ";
		for(int const &gene_id : genome_B_perm){
			genome_B_signs.emplace_back((gene_id < 0));
			genome_B_labels.emplace_back(genome_A.gene_labels_map.getLabel(gene_id));
			// std::cout << gene_id << " [" << genome_A.gene_labels_map.getLabel(gene_id) << "] ";
		}
		// std::cout << std::endl;
		return GenomeMultichrom<int>(genome_B_labels, genome_B_signs, genome_A.gene_labels_map);
	} else {
		return GenomeMultichrom<int>(rng, genome_A.n, genome_A.m, probRev, genome_A.gene_labels_map);
	}
}

void printHelper() {
	std::cout << "Use: ./gensort_random [seed] [genes] [reversals] [tests] [verbose]" << std::endl;
	std::cout << "Example: ./gensort_random 42 10 4 3 2" << std::endl;
	std::cout << "Parameters:" << std::endl;
	std::cout << "- [seed]: any positive integer number greater than 0;" << std::endl;
	std::cout << "- [genes]: the number of genes. In the example above, all random genomes will have 1 chromosome and 10 genes;" << std::endl;
	std::cout << "- [reversals]: the number of reversals separating two genomes. If this value is equal to 0, two random independent permutations will be created. If this value is greater than 0, the second genome is obtained from the first genome by applying the specified number of random reversals. In the example above, all pairs of random genomes will be separated by 4 reversals. Notice that, if the number of reversals is too high, the reversal distance (minimum number of reversals that separates two genomes) will be lower than the actual number of reversals;" << std::endl;
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
	} else if (argc != 6) {
		std::cout << "ERROR! Wrong number of arguments. Program is aborting." << std::endl;
		printHelper();
		exit(1);
	}

	// Parameters.
	int const seed    = std::stoi(argv[1]); // input (command line argument): seed random number generator.
	int const nbgenes = std::stoi(argv[2]); // input (command line argument): number of genes.
	int const nbreversals = std::stoi(argv[3]); // input (command line argument): number of reversals.
	int const nbtests = std::stoi(argv[4]); // input (command line argument): how many random permutations will be generated.
	int const verbose = std::stoi(argv[5]); // input (command line argument): verbose mode (0: [cluster] only basic stats are printed; 1: [user] basic messages; 2: [user] basic messages+replay of the sorting scenario (not recommended if the genome is too big); 3: [developer] very detailed debugging messages.
	int const nbchrom   = 1; // For now only unichromosomal genomes.
	
	if(verbose > 0) {
		std::cout << "**********************************************" << std::endl;
		std::cout << "* Sort by reversals with random permutations *" << std::endl;
		std::cout << "**********************************************" << std::endl << std::endl;
		std::cout << "- Seed to reproduce tests: "    << seed    << std::endl;
		std::cout << "- Genome size (nb. of genes): " << nbgenes << std::endl;
		std::cout << "- Nb. reversals: " << nbreversals << ((nbreversals < 1) ? " (two random independent permutations will be created for each test)" : "") << std::endl;
		std::cout << "- Nb. random permutations: "    << nbtests << std::endl;
	}

	// Create a random number generator
	std::mt19937 rng(seed); // Seed the generator
	// Distribution for the probability of negative signs in genes.
	std::uniform_real_distribution<> distr(0.0, 1.0);

	if(verbose < 1) {std::cout << "test;t_µs;d_obs;d_exp;nb_breakpoints;nb_hurdles;cyc_nontriv;cyc_tot;comp_tot;comp_ori;comp_unori;" << std::endl;}

	// Generate random tests.
	for(int i=0; i<nbtests; ++i){

		// Probability for a gene to have negative sign.
		double const probRev = distr(rng);
		if(verbose > 0) {std::cout << std::endl << std::endl << "Test " << (i+1) << ": Random genome (n=" << nbgenes << ";p_rev=" << probRev << ")" << std::endl;} 

		// Create random genomes using the random constructor.
		GenomeMultichrom<int> genome_A(rng, nbgenes, nbchrom, probRev);
		GenomeMultichrom<int> genome_B = createRandomGenome(genome_A, nbreversals, probRev, rng);

		// Sort genomes using the minimum number of reversals.
		SortByReversals solution = testCase_generalSort(genome_A, genome_B, rng, verbose);
		if(verbose < 1) {std::cout << (i+1) << ";" << solution.t << ";" << solution.obs_distance << ";" << solution.exp_distance << ";" << solution.nb_breakpoints << ";" << solution.nb_hurdles << ";" << solution.nb_cycles_nontrivial << ";" << solution.nb_cycles << ";" << solution.nb_components << ";" << solution.nb_components_oriented << ";" << solution.nb_components_unoriented << ";" << std::endl;}
	} 
	if(verbose > 0) {std::cout << std::endl << "Bye bye" << std::endl;}
	return 0;
}

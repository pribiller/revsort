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
 * reversals corresponds exactly to (method "opt"), 
 * or at least (method "mcmc"), the minimum number expected.
 * 
 * In addition to check the correctness of the method, it also 
 * saves many stats about the permutations (number of breakpoints,
 * number of cycles, etc.) and the time needed to run the method.
 * 
 * Compile (the option -lboost_serialization *must* be in the end):
 * g++ -O3 -fopenmp run_RandomPermutations.cpp sortByReversals.cpp findComponents_Bader2001.cpp reversalMCMC_York2002.cpp sampleReversal_York2002.cpp sortOrientedByReversals_Tannier2007.cpp solveUnoriented_HannenhalliPevzner1999.cpp genome.cpp utils.cpp -o revsort_random -lboost_serialization
 * 
 * Run:
 * ./revsort_random opt 42 10 4 3 2
 * ./revsort_random mcmc.in 42 10 4 1 2
 * 
 *******************************************************/

#include <iostream>
#include <vector>
#include <cstdlib>   // exit
#include <cmath>	 // abs

#include "sortByReversals.hpp"
#include "reversalMCMC_York2002.hpp"

#include "utils.hpp"
#include "inputParameters.hpp"

/*******************************************************
 * Some basic tests.
*******************************************************/

SortByReversals testCase_generalSort(GenomeMultichrom<int>& genome_A, GenomeMultichrom<int>& genome_B, std::mt19937& rng, int debug){
	SortByReversals sortGenome(genome_A,genome_B,debug);
	// Sort genome and measure time.
	sortGenome.sort(rng);
	// [Check 1] Check if the sorting scenario ends with the identity permutation (i.e. no breakpoints).
	bool correctSolution = ((debug >= DEBUG_MEDIUM) ? sortGenome.printSolution() : sortGenome.checkSolution());
	if(!correctSolution){
		std::cout << "ERROR! Problem during the sorting. Program is aborting." << std::endl;
		exit(1);
	}
	// [Check 2] Check if the number of reversals is minimal.
	correctSolution = ((debug >= DEBUG_LOW) ? sortGenome.printStats() : sortGenome.checkDistance());
	if(!correctSolution){
		std::cout << "ERROR! Problem during the sorting. Program is aborting." << std::endl;
		exit(1);
	}
	return sortGenome;
}

ReversalMCMC testCase_reversalMCMC(GenomeMultichrom<int>& genome_A, GenomeMultichrom<int>& genome_B, std::mt19937& rng, int debug, McmcOptions* parameters){
	if(debug >= DEBUG_LOW){std::cout << "Starting " << parameters->print() << "..." << std::endl;}
	ReversalMCMC mcmc(genome_A,genome_B,rng,(*parameters),false);

	// Load last saved state if the option "resume_run" is 'true' and the checkpoint file exists.
	const std::string bkp_filename = mcmc.getBackupFilename();
	const std::string out_filename = mcmc.getSampleFilename();
	if((parameters->resume_run) && isValidPath(bkp_filename)){
		mcmc.loadState(bkp_filename);
		mcmc.saveSamples(out_filename); // Save sampled results so far.
	}
	
	// Update parameters if needed (not all parameters can be modified).
	mcmc.updateParameters((*parameters));
	
	mcmc.run();
	return mcmc;

}

GenomeMultichrom<int> createRandomGenome(GenomeMultichrom<int>& genome_A, const int nb_reversals, const double probRev, std::mt19937& rng, const int debug){
	if(nb_reversals > 0){
		GenomePermutation<BlockSimple> genome_B(genome_A.getExtendedGenome());
		std::vector<Reversal> reversals = applyRandomReversals(genome_B, nb_reversals, rng, debug);
		// Convert labels in the random genome.
		std::vector<int> genome_B_perm = genome_B.getUnextendedPerm();
		std::vector<int> genome_B_labels;
		std::vector<bool> genome_B_signs;
		genome_B_labels.reserve(genome_B_perm.size());
		genome_B_signs.reserve(genome_B_perm.size());
		// std::cout << " - Genome B: ";
		for(int const &gene_id : genome_B_perm){
			genome_B_signs.emplace_back((genome_A.gene_labels_map.getAdjustedId(gene_id) < 0));
			genome_B_labels.emplace_back(genome_A.gene_labels_map.getLabel(gene_id));
			// std::cout << gene_id << " [" << genome_A.gene_labels_map.getLabel(gene_id) << "] ";
		}
		// std::cout << std::endl;

		// Convert labels of reversal history.
		GenomeMultichrom<int> genome_B_rdm(genome_B_labels, genome_B_signs, genome_A.gene_labels_map);

		if(debug >= DEBUG_LOW){
			std::cout << "Random genomes" << std::endl;

			std::cout << "Genome A: ";
			for(int i : genome_A.getExtendedGenome()){std::cout << i << " ";}
			std::cout << std::endl;

			std::cout << "Genome B: ";
			for(int i : genome_B_rdm.getExtendedGenome()){std::cout << i << " ";}
			std::cout << std::endl;
		}

		// std::cout << "Reversal history" << std::endl;
		// RandomReversalScenario rdmReversalScenario;
		// for (Reversal& rev : reversals){
		// 	std::cout << " Reversal (" << rev.g_beg << ", " << rev.g_end << "] becomes ";
		// 	std::cout << " (" << rev.g_beg << ", " << rev.g_beg_next << "]" << std::endl;
		// }

		return genome_B_rdm;
	} else {
		return GenomeMultichrom<int>(rng, genome_A.n, genome_A.m, probRev, genome_A.gene_labels_map);
	}
}

void printHelper() {
	std::cout << "Use: ./revsort_random [method] [seed] [genes] [reversals] [tests] [verbose]" << std::endl;
	std::cout << "Example: ./revsort_random opt 42 10 4 3 2" << std::endl;
	std::cout << "Parameters:" << std::endl;
	std::cout << "- [method]: which method will be used to find an evolutionary path with reversals connecting the two genomes. Three options are available:" << std::endl;
	std::cout << "  - [method=opt]:  Implementation based on Tannier et al. (2007). It computes a scenario with the minimum number of inversions." << std::endl;
	std::cout << "  - [method=mcmc]: Implementation based on York, Durrett, and Nielsen. (2002). It samples reversal scenarios with possibly different number of inversions. The average size of the sampled paths should correspond to the expected number of inversions." << std::endl;
	std::cout << "  - [method=file.in]: The path to a file containing details of the selected method, which allows for more flexibility in setting parameter values specific to each method (see examples in the repository)." << std::endl;
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

RevMethodOptions* checkParameters(std::string& method){
	InputPars inputPars = InputPars();
	if (!inputPars.isValidMethod(method)){
		std::cout << "ERROR! Invalid method (Found: '" << method << "'). Valid options: 'opt'; 'mcmc'. Program is aborting." << std::endl;
		printHelper();
		exit(1);
	}
	return inputPars.getMethodPars(method);
}

int main(int argc, char* argv[]) {
	
	if (argc < 2) {
		std::cout << "ERROR! No extra command line argument passed other than program name. Program is aborting." << std::endl;
		printHelper();
		exit(1);
	} else if (argc != 7) {
		std::cout << "ERROR! Wrong number of arguments. Program is aborting." << std::endl;
		printHelper();
		exit(1);
	}

	// Parameters.
	std::string method  = argv[1];            // input (command line argument): method to find evolutionary path (optimal or MCMC).
	int const seed      = std::stoi(argv[2]); // input (command line argument): seed random number generator.
	int const nbgenes   = std::stoi(argv[3]); // input (command line argument): number of genes.
	int const nbreversals = std::stoi(argv[4]); // input (command line argument): number of reversals.
	int const nbtests   = std::stoi(argv[5]); // input (command line argument): how many random permutations will be generated.
	int const debug     = std::stoi(argv[6]); // input (command line argument): verbose mode (0: [cluster] only basic stats are printed; 1: [user] basic messages; 2: [user] basic messages+replay of the sorting scenario (not recommended if the genome is too big); 3: [developer] very detailed debugging messages.
	int const nbchrom   = 1; // For now only unichromosomal genomes.

	RevMethodOptions* methodPars = checkParameters(method);
	
	if(debug >= DEBUG_LOW) {
		std::cout << "**********************************************" << std::endl;
		std::cout << "* Sort by reversals with random permutations *" << std::endl;
		std::cout << "**********************************************" << std::endl << std::endl;
		std::cout << "- Method: "                     << methodPars->print() << std::endl;
		std::cout << "- Seed to reproduce tests: "    << seed    << std::endl;
		std::cout << "- Genome size (nb. of genes): " << nbgenes << std::endl;
		std::cout << "- Nb. reversals: " << nbreversals << ((nbreversals < 1) ? " (two random independent permutations will be created for each test)" : "") << std::endl;
		std::cout << "- Nb. random permutations: "    << nbtests << std::endl;
	}

	// Create a random number generator
	std::mt19937 rng(seed); // Seed the generator
	// Distribution for the probability of negative signs in genes.
	std::uniform_real_distribution<> distr(0.0, 1.0);

	if(debug == DEBUG_OFF) {std::cout << "test;t_µs;d_obs;d_exp;nb_breakpoints;nb_hurdles;cyc_nontriv;cyc_tot;comp_tot;comp_ori;comp_unori;" << std::endl;}

	// Generate random tests.
	for(int i=0; i<nbtests; ++i){

		// Probability for a gene to have negative sign.
		double const probRev = distr(rng);
		if(debug >= DEBUG_LOW) {std::cout << std::endl << std::endl << "Test " << (i+1) << ": Random genome (n=" << nbgenes << ";p_rev=" << probRev << ")" << std::endl;} 

		// Create random genomes using the random constructor.
		GenomeMultichrom<int> genome_A(rng, nbgenes, nbchrom, probRev);
		GenomeMultichrom<int> genome_B = createRandomGenome(genome_A, nbreversals, probRev, rng, debug);

		// Sort genomes using the expected number of reversals.
		if(methodPars->method == RevMethodType::MCMC){
			ReversalMCMC solution = testCase_reversalMCMC(genome_A, genome_B, rng, debug, dynamic_cast<McmcOptions*>(methodPars));
		// Sort genomes using the minimum number of reversals.
		} else {
			SortByReversals solution = testCase_generalSort(genome_A, genome_B, rng, debug);
			if(debug == DEBUG_OFF) {std::cout << (i+1) << ";" << solution.t << ";" << solution.obs_distance << ";" << solution.exp_distance << ";" << solution.nb_breakpoints << ";" << solution.nb_hurdles << ";" << solution.nb_cycles_nontrivial << ";" << solution.nb_cycles << ";" << solution.nb_components << ";" << solution.nb_components_oriented << ";" << solution.nb_components_unoriented << ";" << std::endl;}
		}

	} 
	if(debug >= DEBUG_LOW) {std::cout << std::endl << "Bye bye" << std::endl;}
	return 0;
}

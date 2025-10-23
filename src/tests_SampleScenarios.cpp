/*******************************************************
 * Original author: Priscila Biller
 * Created: July/2025
 * License: GPL v3
 * 
 * This file contains many small examples of permutations
 * coming from different papers (the paper's reference is 
 * indicated in the name of each function).
 * 
 * These small examples are used for testing the 
 * implementation of the sorting by reversals method.
 * 
 * Compile:
 * g++ -fopenmp tests_SampleScenarios.cpp sortByReversals.cpp reversalMCMC_York2002.cpp sampleReversal_York2002.cpp findComponents_Bader2001.cpp sortOrientedByReversals_Tannier2007.cpp solveUnoriented_HannenhalliPevzner1999.cpp genome.cpp -o revsampler
 * 
 * Run:
 * ./revsampler 42 1
 * 
 *******************************************************/

#include <iostream>
#include <unordered_map>
#include <vector>
#include <stack>
#include <list>
#include <cstdlib>   // exit
#include <cmath>	 // abs
#include <memory>	 // shared_ptr

#include "sampleReversal_York2002.hpp"
#include "reversalMCMC_York2002.hpp"

#include "sortByReversals.hpp"

/*******************************************************
 * Some basic tests.
*******************************************************/

void testCase_reversalMCMC(GenomeMultichrom<int>& genome_A, GenomeMultichrom<int>& genome_B, std::mt19937& rng, const bool debug){

	std::cout << "Starting MCMC..." << std::endl;

	const int nb_chains=10; //10
	const double rev_mean_range=10;
	const int max_steps=100000;
	const int pre_burnin_steps=10;

	// If the probability of bad reversals is too low, then the proposal ratio of better paths decreases:
	// - Proposal ratio increases as P(old|new) / P(new|old).
	// If the old path is composed of one or more bad reversals, and the probability of bad reversals is very low,
	// then the probability of coming back from a better new path to the old path reduces.
	
	// My values.
	const double p_good=1.0; 
	const double p_neutralgood=0.0276; 
	const double p_neutral=0.0276;
	const double p_bad=0.001;
	
	// Larget (2004)'s values.
	// They converge too fast to the minimum number.
	// const double p_good=1.0; 
	// const double p_neutralgood=0.0276; 
	// const double p_neutral=0.0276;
	// const double p_bad=0.0024;

	// Miklos (2003) values (Table 1).
	// const double p_good=1.0; 
	// const double p_neutralgood=0.006; 
	// const double p_neutral=0.003;
	// const double p_bad=0.001;

	// York, Durrett and Nielsen (2002)'s values.
	// They yeld larger paths, most of the time rejected.
	// const double p_good=1.0; 
	// const double p_neutralgood=0.035; 
	// const double p_neutral=0.030;
	// const double p_bad=0.015;

	std::vector<double> probs = {p_good, p_neutralgood, p_neutral, p_bad};

	const double p_stop=0.999; //0.99;

	ReversalMCMC mcmc(genome_A,genome_B,rng,nb_chains,rev_mean_range,max_steps,pre_burnin_steps,probs,p_stop,false);
	mcmc.run();
}

void testCase_sampleModifiedScenario(GenomeMultichrom<int>& genome_A, GenomeMultichrom<int>& genome_B, std::mt19937& rng, const bool debug){
	
	SortByReversals sortGenome(genome_A,genome_B,false);
	//sortGenome.printInputGenomes();

	sortGenome.sort(rng);

	const double p_good=1.0; 
	const double p_neutralgood=0.025; 
	const double p_neutral=0.020;
	const double p_bad=0.015; 
	std::vector<double> probs = {p_good, p_neutralgood, p_neutral, p_bad};

	const double p_stop=0.99; //0.99;
	RandomReversalScenario sampler(probs,p_stop,false);

	std::cout << "Sampling scenario..." << std::endl;
	std::vector<ReversalRandom> reversals = sampler.sampleScenario(genome_B,rng);
	std::cout << " - Sampled scenario = " << reversals.size() << " reversals." << std::endl;

	bool correctSolution = sortGenome.printStats();
		
	std::cout << "\nSampling modified scenario..." << std::endl;
	ProposalReversalScenario proposalHist = sampler.sampleModifiedScenario(genome_B, reversals, rng);
	
	// Compute proposal ratio.
	std::cout << "\nComputing proposal ratio..." << std::endl;
	double propRatio = proposalHist.getProposalRatio(sampler.rev_weights, sampler.p_stop);

	// Compute posterior ratio.
	std::cout << "\nComputing posterior ratio..." << std::endl;

	double rev_mean  = sortGenome.obs_distance;
	double postRatio = proposalHist.getPosteriorRatio(rev_mean);
	std::cout << "\nPosterior ratio = " << postRatio << std::endl;

	// Compute acceptance probability.
	std::cout << "\nComputing acceptance probability..." << std::endl;
	double acceptanceProb = proposalHist.getAcceptanceProb(rev_mean, sampler.rev_weights, sampler.p_stop);
	std::cout << "\nAcceptance prob. = " << acceptanceProb << std::endl;

	std::cout << "\nSampling new mean..." << std::endl;
	// rev_mean_range how far from the current value the new sampled value can be.
	const double rev_mean_range = 10.0;
	ProposalReversalMean proposalMean(rng, rev_mean, rev_mean_range);
	double rev_mean_prop = proposalMean.sampleReversalMean(rev_mean);
	double acceptanceProb_mean = proposalMean.getAcceptanceProb(proposalHist.L_new, rev_mean, rev_mean_prop);
	std::cout << "\nCur. mean= " << rev_mean << "; New mean=" << rev_mean_prop << "; Acceptance prob.=" << acceptanceProb_mean << std::endl;
}

void testCase_sampleReversalScenario(GenomeMultichrom<int>& genome_A, GenomeMultichrom<int>& genome_B, std::mt19937& rng, const bool debug){

	SortByReversals sortGenome(genome_A,genome_B,false);
	sortGenome.sort(rng);


	const double p_good=1.0; 
	const double p_neutralgood=0.025; 
	const double p_neutral=0.020;
	const double p_bad=0.015; 
	std::vector<double> probs = {p_good, p_neutralgood, p_neutral, p_bad};

	const double p_stop=0.99; //0.99;
	RandomReversalScenario sampler(probs,p_stop,true);

	std::cout << "Sampling scenario..." << std::endl;
	std::vector<ReversalRandom> reversals = sampler.sampleScenario(genome_B,rng);
	std::cout << " - Sampled scenario = " << reversals.size() << " reversals." << std::endl;

	bool correctSolution = sortGenome.printStats();
}

void testCase_sampleOneReversal(GenomeMultichrom<int>& genome_A, GenomeMultichrom<int>& genome_B, std::mt19937& rng, const bool debug){

	GenomePermutation<BlockSimple> genperm(genome_B.getExtendedGenome());

	ReversalSampler sampler(genperm,debug);
	sampler.countReversals();

	std::cout << "Sampling reversal..." << std::endl;
	const int nb_reversals = 100;
	for (int i = 0; i < nb_reversals; ++i) {
		ReversalRandom reversal = sampler.sampleReversal(rng,false);
		std::cout << " Type: " << reversal.type << "; Cycle: " << reversal.cycle << " Reversal (" << reversal.g_beg << ", " << reversal.g_end << "] " << std::endl;
	}
}

/*******************************************************
* Tests for the overall algorithm (find components, 
* clear hurdles, sort connected components, etc.).
*******************************************************/

// Example used in the paper from Garg et al. (2019).
// {-2,5,4,-1,3,6,9,-7,-8}
// Reversal distance = 5 reversals. 
// Details for Reversal distance (d) computation: 10 breakpoints (b); 5 cycles(c); 0 hurdles(h): d = b-c+h (+1 if fortress).
void testCase_Garg2019(std::mt19937& rng, bool debug){
	std::vector<int>  genome_multichrom_A  = {1, 2, 3, 4, 5, 6, 7, 8, 9};
	std::vector<bool> genome_orientation_A = {false, false, false, false, false, false, false, false, false};

	std::vector<int>  genome_multichrom_B  = {2, 5, 4, 1, 3, 6, 9, 7, 8};
	std::vector<bool> genome_orientation_B = {true, false, false, true, false, false, false, true, true};

	GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
	GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);

	std::cout << "\n\nTest: Example from Garg et al.(2019)\n";
	testCase_reversalMCMC(genome_A, genome_B, rng, debug);
}

// Example used in the paper from Bader et al. (2001).
// (+3, +9, −7, +5, −10, +8, +4, −6, +11, +2, +1)
// This is the only case where the lower bound for the 
// reversal distance does not match the value found:
// - Reversal distance (d) = expected: 7 reversals (or 8 if there is a fortress); found: 8 reversals.
// The graph has 1 unoriented component composed of **two** intertwining cycles.
// As the definition of hurdle seems to be a cycle in an unoriented component with 
// all its elements appearing in a consecutive order, apparently 
// in this case there are no hurdles (h=0). 
// However, the program finds h=1. Notice that, if h=2, the expected and observed values would match.
void testCase_Bader2001(std::mt19937& rng, bool debug){
	std::vector<int>  genome_multichrom_A  = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
	std::vector<bool> genome_orientation_A = {false, false, false, false, false, false, false, false, false, false, false};

	std::vector<int>  genome_multichrom_B  = {3, 9, 7, 5, 10, 8, 4, 6, 11, 2, 1};
	std::vector<bool> genome_orientation_B = {false, false, true, false, true, false, false, true, false, false, false};

	GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
	GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);

	std::cout << "\n\nTest: Example from Bader et al.(2001)\n";
	testCase_reversalMCMC(genome_A, genome_B, rng, debug);
}

// Example used in the paper from Tannier et al. (2007) (Figure 4).
// {0, -1, 3, 2, 4}
void testCase_Tannier2007_Figure4(std::mt19937& rng, bool debug){
	std::vector<int>  genome_multichrom_A  = {1, 2, 3, 4, 5};
	std::vector<bool> genome_orientation_A = {false, false, false, false, false};

	std::vector<int>  genome_multichrom_B  = {1, 2, 4, 3, 5};
	std::vector<bool> genome_orientation_B = {false, true, false, false, false};

	GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
	GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);

	std::cout << "\n\nTest: Example from Tannier et al. (2007) (Figure 4)\n";
	testCase_reversalMCMC(genome_A, genome_B, rng, debug);
}

// Example used in the paper from Hannehalli and Pevzner (1999) (Figure 4(a)).
// {+5, +7, +6, +8, +1, +3, +2, +4}
void testCase_Hannehalli1999_Fig4a(std::mt19937& rng, bool debug){
	std::vector<int>  genome_multichrom_A  = {1, 2, 3, 4, 5, 6, 7, 8};
	std::vector<bool> genome_orientation_A = {false, false, false, false, false, false, false, false};

	std::vector<int>  genome_multichrom_B  = {5, 7, 6, 8, 1, 3, 2, 4};
	std::vector<bool> genome_orientation_B = {false, false, false, false, false, false, false, false};

	GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
	GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);

	std::cout << "\n\nTest: Example from Hannehalli and Pevzner (1999) - Figure 4(a)\n";
	testCase_reversalMCMC(genome_A, genome_B, rng, debug);
}

// Example used in the paper from Hannehalli and Pevzner (1999) (Figure 4(b)).
// {+2, +4, +3, +5, +7, +6, +8, +1}
void testCase_Hannehalli1999_Fig4b(std::mt19937& rng, bool debug){
	std::vector<int>  genome_multichrom_A  = {1, 2, 3, 4, 5, 6, 7, 8};
	std::vector<bool> genome_orientation_A = {false, false, false, false, false, false, false, false};

	std::vector<int>  genome_multichrom_B  = {2, 4, 3, 5, 7, 6, 8, 1};
	std::vector<bool> genome_orientation_B = {false, false, false, false, false, false, false, false};

	GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
	GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);

	std::cout << "\n\nTest: Example from Hannehalli and Pevzner (1999) - Figure 4(b)\n";
	testCase_reversalMCMC(genome_A, genome_B, rng, debug);
}

// Example used in the book ``Mathematics of Evolution and Phylogeny`` (2005) (Section 10.4.2).
// {0, 2, 1, 3, 5, 7, 6, 8, 9, 4, 10}
void testCase_Bergeron2005_Sec10_4_2(std::mt19937& rng, bool debug){
	std::vector<int>  genome_multichrom_A  = {1, 2, 3, 4, 5, 6, 7, 8, 9};
	std::vector<bool> genome_orientation_A = {false, false, false, false, false, false, false, false, false};

	std::vector<int>  genome_multichrom_B  = {2, 1, 3, 5, 7, 6, 8, 9, 4};
	std::vector<bool> genome_orientation_B = {false, false, false, false, false, false, false, false, false};

	GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
	GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);

	std::cout << "\n\nTest: Example from the book 'Mathematics of Evolution and Phylogeny' (2005) (Section 10.4.2)\n";
	testCase_reversalMCMC(genome_A, genome_B, rng, debug);
}

// Example used in the book ``Mathematics of Evolution and Phylogeny`` (2005) (Figure 10.6).
// P_2 = {0, -3, 1, 2, 4, 6, 5, 7, -15, -13, -14, -12, -10, -11, -9, 8, 16}
// d(P_2) = 13.
void testCase_Bergeron2005_Fig10_6(std::mt19937& rng, bool debug){
	std::vector<int>  genome_multichrom_A  = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
	std::vector<bool> genome_orientation_A = {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false};

	std::vector<int>  genome_multichrom_B  = {3, 1, 2, 4, 6, 5, 7, 15, 13, 14, 12, 10, 11, 9, 8};
	std::vector<bool> genome_orientation_B = {true, false, false, false, false, false, false, true, true, true, true, true, true, true, false};

	GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
	GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);

	std::cout << "\n\nTest: Example from the book 'Mathematics of Evolution and Phylogeny' (2005) (permutation P_2, Figure 10.6)\n";
	testCase_reversalMCMC(genome_A, genome_B, rng, debug);
}

int main(int argc, char* argv[]) {
	
	if (argc < 2) {
		std::cout << "ERROR! No extra command line argument passed other than program name. Program is aborting." << std::endl;
		exit(1);
	} else if (argc != 3) {
		std::cout << "ERROR! Wrong number of arguments. Expected: ./program [seed] [debug]. Program is aborting." << std::endl;
		exit(1);
	}

	// Parameters specified by the user.
	int const seed   = std::stoi(argv[1]); // input (command line argument): seed random number generator.
	bool const debug = (std::stoi(argv[2]) > 0);

	// Create a random number generator
	std::cout << "- Seed to reproduce tests: " << seed << std::endl;
	std::mt19937 rng(seed); // Seed the generator

	// Examples from various papers.
	// testCase_Garg2019(rng,debug);               // --> work: OK
	// testCase_Bader2001(rng,debug);              // --> work: OK* (check the header of the function)
	// testCase_Hannehalli1999_Fig4a(rng,debug);   // --> work: OK
	// testCase_Hannehalli1999_Fig4b(rng,debug);   // --> work: OK
	// testCase_Bergeron2005_Sec10_4_2(rng,debug); // --> work: OK
	testCase_Bergeron2005_Fig10_6(rng,debug);   // --> work: OK
	// testCase_Tannier2007_Figure4(rng,debug);    // --> work: OK
	
	std::cout << std::endl << "Bye bye" << std::endl;
	return 0;
}

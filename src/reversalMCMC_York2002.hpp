/*******************************************************
 * Original author: Priscila Biller
 * Created: October/2025
 * License: GPL v3
 * 
 * This class implements a method to sample a random
 * sorting scenario by reversals given two genomes,
 * i.e., it finds a sequence of reversals to transform 
 * one of the input genomes into the other genome.
 * Note that the sequence of reversals might not be optimal,
 * specially if the two given genomes are far apart.
 * 
 * The given scenario is just one of the many possible
 * scenarios that could transform one genome
 * into another.
 * 
 * To find a random scenario, a method similar to York's
 * work is implemented. [see York et al. (2002); the 
 * implemented method is also similar to 
 * Larget et al. (2004), Larget et al. (2005)] 
 * 
 * The method implements a MCMCMC (also know as 
 * MCMC with Parallel Tempering), in which each
 * several chains are ran in parallel and swaped
 * every once in a while. The chains have different
 * 'temperatures', a parameter that controls how
 * flat or rugged is the solution space. Flatter 
 * landscapes allow for bigger jumps in the solution space.
 * 
 * To find a random reversal history, at each step
 * a random reversal is sampled. Reversals are categorized
 * in 4 types: "good reversals", "neutral good reversals",
 * "neutral reversals", and "bad reversals". These types are 
 * defined based on how close the genome gets to the target 
 * genome after they are applied.
 * 
 * References
 * ----------
 * 
 * For additional information on similar methods that inspired 
 * this implementation, please check:
 * 
 * - Thomas York, Richard Durrett, and Rasmus Nielsen. "Bayesian estimation 
 * of the number of inversions in the history of two chromosomes". Journal of 
 * Computational Biology (2002), 9(6), 805-818.
 * 
 * - Bret Larget et al. "A Bayesian analysis of metazoan mitochondrial 
 * genome arrangements." Molecular Biology and Evolution 22.3 (2005): 486-495.
 * 
 * - Bret Larget et al. "Bayesian phylogenetic inference from animal mitochondrial 
 * genome arrangements." Journal of the Royal Statistical Society Series B: 
 * Statistical Methodology 64.4 (2002): 681-693.
 * 
 * - Miklós, István, and Aaron E. Darling. "Efficient sampling of parsimonious 
 * inversion histories with application to genome rearrangement in Yersinia." 
 * Genome biology and evolution 1 (2009): 153-164.
 * 
 * Notice that the method implemented here is similar to, but not exactly 
 * the same as the ones mentioned above. It should be seen more like a "mix"
 * of the mentioned methods: 
 * 
 * 1) The generation of a reversal history is similar to Miklos et al. (2009);
 * 
 * 2) The proposal probability, acceptance probability, etc., are computed in a 
 * similar way to the method of York et al. (2002).
 * 
 * 3) The parallel tempering is similar to Larget et al. (2005).
 * 
 *******************************************************/

#pragma once // It avoids class redefinition.

#include <iostream>
#include <unordered_map>
#include <algorithm> // set_difference, sort
#include <vector>
#include <stack>
#include <list>
#include <cstdlib>   // exit
#include <cmath>	 // abs
#include <numeric>   // iota
#include <memory>	 // shared_ptr
#include <chrono>
#include <utility>   // swap 

#include "genome.hpp"
#include "reversal.hpp"
#include "findComponents_Bader2001.hpp"
#include "solveUnoriented_HannenhalliPevzner1999.hpp"

#include "sortByReversals.hpp"
#include "sampleReversal_York2002.hpp"
#include "inputParameters.hpp"

class ProposalReversalMean {
public:
	std::mt19937& rng;

	double rev_mean_min;
	double rev_mean_range; // Proposed values are uniformly sampled from the range [cur-range, cur+range].

	double acceptanceProb{0.0};
	double rev_mean_new{0.0};

	// The occurrence of inversions is a Poisson process with unknown mean lambda.
	ProposalReversalMean(std::mt19937& rng_, double rev_mean_min_, double rev_mean_range_):rng(rng_),rev_mean_min(rev_mean_min_),rev_mean_range(rev_mean_range_){}
	ProposalReversalMean(std::mt19937& rng_):rng(rng_){}

	double sampleReversalMean(const double rev_mean_cur);
	double getAcceptanceProb(const int rev_path_len, const double rev_mean_cur, const double rev_mean_prop);
	inline double proposeNewValue(const int rev_path_len, const double rev_mean_cur) {return getAcceptanceProb(rev_path_len, rev_mean_cur, sampleReversalMean(rev_mean_cur));}
};

class ProposalReversalScenario {
public:
	std::vector<ReversalRandom> currentReversalScenario;
	std::vector<ReversalRandom> proposedReversalScenario;

	int path_beg; // Starting position of the replaced path (same position both in the current and proposed scenarios).

	int L_cur; // Length of the current reversal scenario (total number of reversals).
	int l_cur; // Length of the replaced sub-path in the current reversal scenario (number of reversals replaced).

	int L_new; // Length of the proposed reversal scenario (total number of reversals).
	int l_new; // Length of the new sub-path in the proposed reversal scenario (number of new reversals).

	int N; // Number of markers in the unextended uncapped permutation.

	// Sampling subpath parameters.
	// q(l) ~ 1 - tanh(epsilon*(l/(alpha*N)-1))
	// - alpha: lengths small than N*alpha are roughly equally likely represented.
	//          For example, alpha=0.65, means that there is more or less the same 
	//          chance to sample a size between 0 and 65% of the total size of the 
	//          path (N). From this point on, the probability drops to almost 0.
	double alpha{0.65}; 
	double epsilon{8};

	double proposalRatio{0.0};  // Updated in getProposalRatio.
	double posteriorRatio{0.0}; // Updated in getPosteriorRatio.

	bool debug{false};

	ProposalReversalScenario(std::vector<ReversalRandom>& currentReversalScenario_, std::vector<ReversalRandom>& proposedReversalScenario_, 
		const int path_beg_, const int L_cur_, const int l_cur_, const int L_new_, const int l_new_, const int N_):currentReversalScenario(currentReversalScenario_),proposedReversalScenario(proposedReversalScenario_),path_beg(path_beg_),L_cur(L_cur_),l_cur(l_cur_),L_new(L_new_),l_new(l_new_),N(N_){
	}

	// Probability of sampling a subpath of size l from a path of total size L (size = number of reversals?).
	inline double q_L(const int L, const int l) const {return 1.0-std::tanh(epsilon*(l/(alpha*L)-1));}

	// Probability of sampling a subpath of size l from a path of total size L such that the path starts at position j.
	inline double q_lj(const int L, const int l, const int j) const {return q_L(L,l)/(L+1-l);}

	// Probability of proposing a particular path with l inversions.
	std::vector<double> q_path_factors(const std::vector<ReversalRandom>& path, const std::vector<double>& rev_weights, const double p_stop);

	// Proposal ratio (Hastings ratio).
	double getProposalRatio(const std::vector<double>& rev_weights, const double p_stop);

	// Posterior ratio.
	double getPosteriorRatio(const double rev_mean);

	// Acceptance probability.
	double getAcceptanceProb(const double rev_mean, const std::vector<double>& rev_weights, const double p_stop);

};

class RandomReversalScenario {
// private:

public:

	std::vector<double> rev_weights;  // Weight of each type of reversal.
	double p_stop{0.99};

	bool debug{false};
		
	RandomReversalScenario(const std::vector<double>& rev_weights_, const double p_stop, const bool debug):rev_weights(rev_weights_),p_stop(p_stop),debug(debug){
	}

	std::vector<ReversalRandom> getSubpath(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, const int pos_beg, const int pos_end);

	std::vector<ReversalRandom> sampleScenario(GenomeMultichrom<int>& genome_B, std::mt19937& rng);
	ProposalReversalScenario sampleModifiedScenario(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, std::mt19937& rng);

	std::vector<ReversalRandom> updateReversalScenario(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom> reversals, std::vector<ReversalRandom> reversals_new, const int pos_beg, const int pos_end);

	//int samplePathLength(std::mt19937& rng, const int N, const double alpha=0.65, const double epsilon=8);
	int samplePathLength(std::mt19937& rng, const int N, const double alpha=0.85, const double epsilon=8);
	int samplePathStart(std::mt19937& rng, const int N, const int l);

	std::pair<std::vector<int>,std::vector<int>> getPathEnds(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, const int pos_beg, const int pos_end);
	GenomeMultichrom<int> getGenomes(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, const int pos_beg, const int pos_end);
	int getGeneLabelDefault(GenomeMultichrom<int>& genome, int gene_lbl_ext);
	int getGeneLabelSubpath(GenomePermutation<BlockSimple>& genperm, GenomeMultichrom<int>& genome, int gene_lbl_ext);
	int getGeneExtremity(const int gene_id, const GenomePermutation<BlockSimple>& genperm) const;
	
	void printGenome(std::vector<int> perm);
	void printSortingScenario(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, int pos_beg, int pos_end);

};

class ReversalMCMC {
public:

	int nb_chains;

	GenomeMultichrom<int>& genome_B;
	std::mt19937& rng;
	bool debug;

	int rev_dist{-1}; // It makes sure that estimated nb. of reversals is equal or above the reversal distance.
	double rev_mean_range; // Parameter used to propose new values for the mean based on the current value.
	const std::vector<double>& rev_weights;  // Weight of each type of reversal.
	double p_stop{0.99};

	// It stores the current state of all chains.
	std::vector<std::vector<ReversalRandom>> currentState_revHists;
	std::vector<double> currentState_revMeans;

	// It stores the history of a chain, i.e., all the states visited 
	// so far. More specifically, given a chain i, the length of all reversal paths 
	// visited by chain i are the keys of the map stored at position i,
	// whereas the frequency of the path lengths are the values of map i.
	// To save space, the order in which the states are visited is not saved.
	std::vector<std::unordered_map<int,int>> hist_chains;

	// Convergence measures.
	double max_steps;
	double pre_burnin_steps;
	double cur_step{0}; // It stores the current step.
	std::vector<double> rev_path_avgsize; // Position i stores the average path size in chain i so far.

	// It stores sampled values.
	std::vector<std::vector<ReversalRandom>> sampled_revHists;
	std::vector<double> sampled_revMeans;

	// Distribution used for computing acceptance probabilities.
	std::uniform_real_distribution<double> distr_accept;

	ProposalReversalMean proposalMean;

	// The occurrence of inversions is a Poisson process with unknown mean lambda.
	ReversalMCMC(GenomeMultichrom<int>& genome_A_, GenomeMultichrom<int>& genome_B_, std::mt19937& rng_, const int nb_chains_, const double rev_mean_range_, const int max_steps_, const int pre_burnin_steps_, std::vector<double>& rev_weights_, double p_stop_, const bool debug_):nb_chains(nb_chains_),genome_B(genome_B_),rng(rng_),max_steps(max_steps_),pre_burnin_steps(pre_burnin_steps_),debug(debug_),rev_mean_range(rev_mean_range_),currentState_revHists(nb_chains_),currentState_revMeans(nb_chains_),distr_accept(0.0, 1.0),proposalMean(rng_),rev_path_avgsize(nb_chains_),rev_weights(rev_weights_),p_stop(p_stop_),hist_chains(nb_chains_){

		// Compute reversal distance.
		SortByReversals sortGenome(genome_A_,genome_B_,false);
		sortGenome.sort(rng);
		rev_dist = sortGenome.obs_distance;

		// Initialize proposal mean.
		proposalMean.rev_mean_min   = rev_dist;
		proposalMean.rev_mean_range = rev_mean_range;

		// Initialize reversal histories/means.
		initializeChains();
	}

	// The occurrence of inversions is a Poisson process with unknown mean lambda.
	ReversalMCMC(GenomeMultichrom<int>& genome_A_, GenomeMultichrom<int>& genome_B_, std::mt19937& rng_, const McmcOptions& parameters, const bool debug_):nb_chains(parameters.nb_chains),genome_B(genome_B_),rng(rng_),max_steps(parameters.max_steps),pre_burnin_steps(parameters.pre_burnin_steps),debug(debug_),rev_mean_range(parameters.rev_mean_range),currentState_revHists(parameters.nb_chains),currentState_revMeans(parameters.nb_chains),distr_accept(0.0, 1.0),proposalMean(rng_),rev_path_avgsize(parameters.nb_chains),rev_weights(parameters.probs),p_stop(parameters.p_stop),hist_chains(parameters.nb_chains){

		// Compute reversal distance.
		SortByReversals sortGenome(genome_A_,genome_B_,false);
		sortGenome.sort(rng);
		rev_dist = sortGenome.obs_distance;

		// Initialize proposal mean.
		proposalMean.rev_mean_min   = rev_dist;
		proposalMean.rev_mean_range = rev_mean_range;
		
		std::cout << "- Reversal distance = " << rev_dist << std::endl;
		
		// Initialize reversal histories/means.
		initializeChains();
	}

	void initializeChains();
	std::string runSingleChain(const int chainIdx);
	void run();

	double computeWithinChainVariance();
	double computeBetweenChainVariance();
};

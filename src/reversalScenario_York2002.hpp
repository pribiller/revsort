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

#include "sampleReversal_York2002.hpp"

class ProposalReversalMean {
public:

	double rev_mean_min;
	double rev_mean_range; // Proposed values are uniformly sampled from the range [cur-range, cur+range].

	// The occurrence of inversions is a Poisson process with unknown mean lambda.
	ProposalReversalMean(double rev_mean_min_, double rev_mean_range_):rev_mean_min(rev_mean_min_),rev_mean_range(rev_mean_range_){}

	double sampleReversalMean(std::mt19937& rng, const double rev_mean_cur) const;

	double getAcceptanceProb(const int rev_path_len, const double rev_mean_cur, const double rev_mean_prop) const;
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

	bool debug{false};

	ProposalReversalScenario(std::vector<ReversalRandom>& currentReversalScenario_, std::vector<ReversalRandom>& proposedReversalScenario_, 
		const int path_beg_, const int L_cur_, const int l_cur_, const int L_new_, const int l_new_, const int N_):currentReversalScenario(currentReversalScenario_),proposedReversalScenario(proposedReversalScenario_),path_beg(path_beg_),L_cur(L_cur_),l_cur(l_cur_),L_new(L_new_),l_new(l_new_),N(N_){
	}

	// Probability of sampling a subpath of size l from a path of total size L (size = number of reversals?).
	inline double q_L(const int L, const int l) const {return 1.0-std::tanh(epsilon*(l/(alpha*L)-1));}

	// Probability of sampling a subpath of size l from a path of total size L such that the path starts at position j.
	inline double q_lj(const int L, const int l, const int j) const {return q_L(L,l)/(L+1-l);}

	// Probability of proposing a particular path with l inversions.
	std::vector<double> q_path_factors(const std::vector<ReversalRandom>& path, const std::vector<float>& rev_weights, const double p_stop);

	// Proposal ratio (Hastings ratio).
	double getProposalRatio(const std::vector<float>& rev_weights, const double p_stop);

	// Posterior ratio.
	double getPosteriorRatio(const double rev_mean);

	// Acceptance probability.
	double getAcceptanceProb(const double rev_mean, const std::vector<float>& rev_weights, const double p_stop);

};

class RandomReversalScenario {
// private:

public:

	std::vector<float> rev_weights;  // Weight of each type of reversal.
	float p_stop{0.99};

	bool debug{false};

	RandomReversalScenario(const bool debug=false, const float p_good=1.0, const float p_neutralgood=0.025, const float p_neutral=0.020, const float p_bad=0.015, const float p_stop=0.99):rev_weights(ReversalType_COUNT, 0),p_stop(p_stop),debug(debug){
		rev_weights = {p_good, p_neutralgood, p_neutral, p_bad};
	}
	
	std::vector<ReversalRandom> getSubpath(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, const int pos_beg, const int pos_end);

	std::vector<ReversalRandom> sampleScenario(GenomeMultichrom<int>& genome_B, std::mt19937& rng);
	ProposalReversalScenario sampleModifiedScenario(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, std::mt19937& rng);

	std::vector<ReversalRandom> updateReversalScenario(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom> reversals, std::vector<ReversalRandom> reversals_new, const int pos_beg, const int pos_end);

	int samplePathLength(std::mt19937& rng, const int N, const float alpha=0.65, const float epsilon=8);
	int samplePathStart(std::mt19937& rng, const int N, const int l);

	std::pair<std::vector<int>,std::vector<int>> getPathEnds(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, const int pos_beg, const int pos_end);
	GenomeMultichrom<int> getGenomes(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, const int pos_beg, const int pos_end);
	int getGeneLabelDefault(GenomeMultichrom<int>& genome, int gene_lbl_ext);
	int getGeneLabelSubpath(GenomePermutation<BlockSimple>& genperm, GenomeMultichrom<int>& genome, int gene_lbl_ext);
	int getGeneExtremity(const int gene_id, const GenomePermutation<BlockSimple>& genperm) const;
	
	void printGenome(std::vector<int> perm);
	void printSortingScenario(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, int pos_beg, int pos_end);

};

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

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/optional.hpp>

#include "utils.hpp"
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

	double acceptanceProb{0.0};
	double rev_mean_new{0.0};

	// The occurrence of inversions is a Poisson process with unknown mean lambda.
	ProposalReversalMean(std::mt19937& rng_, double rev_mean_min_):rng(rng_),rev_mean_min(rev_mean_min_){}
	ProposalReversalMean(std::mt19937& rng_):rng(rng_){}

	// Serialization with Boost.
	template<class Archive>
	void serialize(Archive& ar, const unsigned int version) {
		ar & rev_mean_min;
		ar & acceptanceProb;
		ar & rev_mean_new;
		ar & rng;
	}

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

	std::string id_run{"revMCMC"};
	int nb_chains;

	GenomeMultichrom<int>& genome_B;
	std::mt19937& rng;
	bool debug;

	int rev_dist{-1}; // It makes sure that estimated nb. of reversals is equal or above the reversal distance.
	std::vector<double>& rev_weights; // Weight of each type of reversal.
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

	// Convergence evaluation: Same method as in York, Durrett, and Nielsen (2002).
	// We use the method of Gelman and Rubin (1992) to decide when the Markov chain has converged.
	// If this flag is true,  sampling starts after convergence criteria is satisfied.
	// If this flag is false, sampling starts after the end of the pre-burn-in phase (mode used for tests).
	bool check_convergence{true}; // Criteria: sqrt(R) <= 1.1

	// Method of Gelman and Rubin (1992).
	double W{1000}; // within-chain variance.
	double B{1000}; // between chain variance.
	double R{1000}; // convergence measure. Initialize w/a value higher than 1.1.

	// Test parameters.
	// In case the flag 'check_convergence' is false,
	// then these parameters are used instead. 
	// Sampling starts directly after pre_burning_steps.
	int max_steps{100000};
	int pre_burnin_steps{1000};

	// Sampling parameters.
	int sample_interval{1000}; // Based on Larget et al. (2004)
	int sample_amount{10000};  // 10 million steps. Based on Miklos and Darling (2009).

	// Steps in which the whole state of the MCMC is saved.
	// This allows the program to recover in case of crash, time out, 
	// or any other problem.
	int backup_interval{500};
	int print_interval{50};

	// It stores sampled values.
	std::vector<std::vector<ReversalRandom>> sampled_revHists;
	std::vector<double> sampled_revMeans;

	double cur_step{0}; // It stores the current step.
	std::vector<double> rev_path_avgsize; // Position i stores the average path size in chain i so far.

	// Distribution used for computing acceptance probabilities.
	std::uniform_real_distribution<double> distr_accept;

	ProposalReversalMean proposalMean;

	ReversalMCMC(GenomeMultichrom<int>& genome_B_, std::mt19937& rng_, std::vector<double>& rev_weights_):genome_B(genome_B_),rng(rng_),rev_weights(rev_weights_),proposalMean(rng_){}

	// The occurrence of inversions is a Poisson process with unknown mean lambda.
	ReversalMCMC(GenomeMultichrom<int>& genome_A_, GenomeMultichrom<int>& genome_B_, std::mt19937& rng_, 
		const std::string id_run_, const int nb_chains_, const bool check_convergence_, const int max_steps_, const int pre_burnin_steps_, 
		const int sample_interval_, const int sample_amount_, int backup_interval_, int print_interval_,
		std::vector<double>& rev_weights_, double p_stop_, const bool debug_):genome_B(genome_B_),rng(rng_),
		id_run(id_run_),nb_chains(nb_chains_),check_convergence(check_convergence_),max_steps(max_steps_),pre_burnin_steps(pre_burnin_steps_),
		sample_interval(sample_interval_),sample_amount(sample_amount_),backup_interval(backup_interval_),print_interval(print_interval_),
		rev_weights(rev_weights_),p_stop(p_stop_),debug(debug_),
		currentState_revHists(nb_chains_),currentState_revMeans(nb_chains_),rev_path_avgsize(nb_chains_),hist_chains(nb_chains_),
		distr_accept(0.0, 1.0),proposalMean(rng_){

		// Compute reversal distance.
		SortByReversals sortGenome(genome_A_,genome_B_,false);
		sortGenome.sort(rng);
		rev_dist = sortGenome.obs_distance;

		// Initialize proposal mean.
		proposalMean.rev_mean_min   = rev_dist;

		// Initialize reversal histories/means.
		initializeChains();
	}

	// The occurrence of inversions is a Poisson process with unknown mean lambda.
	ReversalMCMC(GenomeMultichrom<int>& genome_A_, GenomeMultichrom<int>& genome_B_, std::mt19937& rng_, 
		McmcOptions& parameters, const bool debug_):genome_B(genome_B_),rng(rng_),
		id_run(parameters.id_run),nb_chains(parameters.nb_chains),check_convergence(parameters.check_convergence),
		max_steps(parameters.max_steps),pre_burnin_steps(parameters.pre_burnin_steps),
		sample_interval(parameters.sample_interval),sample_amount(parameters.sample_amount),
		backup_interval(parameters.backup_interval),print_interval(parameters.print_interval),
		rev_weights(parameters.probs),p_stop(parameters.p_stop),debug(debug_),
		currentState_revHists(parameters.nb_chains),currentState_revMeans(parameters.nb_chains),
		rev_path_avgsize(parameters.nb_chains),hist_chains(parameters.nb_chains),
		distr_accept(0.0, 1.0),proposalMean(rng_){

		// Compute reversal distance.
		SortByReversals sortGenome(genome_A_,genome_B_,false);
		sortGenome.sort(rng);
		rev_dist = sortGenome.obs_distance;

		// Initialize proposal mean.
		proposalMean.rev_mean_min   = rev_dist;
		
		std::cout << "- Reversal distance = " << rev_dist << std::endl;

		// Initialize reversal histories/means.
		initializeChains();
	}

	// Serialization with Boost.
	template<class Archive>
	void serialize(Archive& ar, const unsigned int version) {
		ar & id_run;
		ar & nb_chains;
		ar & genome_B;
		ar & rng;
		ar & debug;
		ar & rev_dist;
		ar & rev_weights;
		ar & p_stop;
		ar & currentState_revHists;
		ar & currentState_revMeans;
		ar & hist_chains;
		ar & check_convergence;
		ar & W;
		ar & B;
		ar & R;
		ar & max_steps;
		ar & pre_burnin_steps;
		ar & sample_interval;
		ar & sample_amount;
		ar & backup_interval;
		ar & print_interval;
		ar & sampled_revHists;
		ar & sampled_revMeans;
		ar & cur_step;
		ar & rev_path_avgsize;
		ar & proposalMean;
		// Special to serialize std::uniform_real_distribution<double>.
		if (typename Archive::is_saving()){
			double lower_bound = distr_accept.a(); // Get lower bound
			double upper_bound = distr_accept.b(); // Get upper bound
			ar & lower_bound;
			ar & upper_bound;
		} else {
			double lower_bound;
			double upper_bound;
			ar & lower_bound;
			ar & upper_bound;
			// Reconstruct the distribution
			distr_accept = std::uniform_real_distribution<double>(lower_bound, upper_bound);
		}
	}

	void initializeIdRun(); // TODO: Currently not used.
	void initializeChains();
	std::string runSingleChain(const int chainIdx);
	void run();

	// Save the state of the object to a file
	void saveState(const std::string& filename) const {
		std::ofstream ofs(filename);
		boost::archive::binary_oarchive oa(ofs);
		oa << *this; // Serialize the object
	}

	// Load the state of the object from a file
	void loadState(const std::string& filename) {
		std::ifstream ifs(filename);
		boost::archive::binary_iarchive ia(ifs);
		ia >> *this; // Deserialize the object
	}

	double computeWithinChainVariance();
	double computeBetweenChainVariance();
};

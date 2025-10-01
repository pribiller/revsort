/*******************************************************
 * Original author: Priscila Biller
 * Created: October/2025
 * License: GPL v3
 * 
 * This class implements a method to sample one random reversal. 
 * Reversals are categorized in 4 types: "good reversals", 
 * "neutral good reversals", "neutral reversals", and "bad reversals". 
 * These types are defined based on how close the genome gets to the target 
 * genome after they are applied:
 * 
 * 1) Good reversals
 * -----------------
 * 
 * Good reversals (also called "proper reversals" in Larget et al.) 
 * are any reversals that could belong to a parsimonious scenario.
 * Basically, there are two types, depending on whether the cycle 
 * belongs to an oriented or unoriented component in the Breakpoint Graph:
 * 
 * 1.1) Oriented components: any reversal that affects two "oriented" 
 *      black edges belonging to the same cycle in the breakpoint graph. 
 *      It results in one cycle being split into two, which increases 
 *      the number of cycles by 1.
 * 
 * 1.2) Unoriented components: any reversal that merges two non-adjacent 
 *      "hurdles". Although the number of components decreases by 1, they 
 *      are also considered "good", as they optimally transform unoriented 
 *      components into oriented components.
 *      Here are few exceptional cases that should be considered:
 *      - If the number of hurdles is equal to 1, then any reversal inside
 *        the component is considered "good";
 *      - If the number of hurdles is equal to 2 or 3, then any reversal 
 *        involving adjacent hurdles is considered "good".
 * 
 * 2) Neutral good reversals
 * -------------------------
 * 
 * Neutral good reversals is a new category in this implementation.
 * They occur only in unoriented components, and can be one of those:
 * 
 * 1.1) Hurdles (if the number of hurdles is >3): any reversal that 
 *      affects two "oriented" black edges belonging to the same hurdle 
 *      or adjacent hurdles. Although they do not occur in parsimonious 
 *      scenarios, these reversals transform unoriented components into 
 *      oriented ones, helping to approach the target in future steps.
 * 
 * 1.2) Unoriented components (non-trivial, non-hurdle): any reversal 
 *      that affects two "oriented" black edges belonging to the same 
 *      cycle. Although they might not occur in parsimonious scenarios, 
 *      these reversals transform unoriented components into oriented ones, 
 *      helping to approach the target in future steps.
 *      
 * Initially, Larget et al. considered these reversals as "neutral"
 * and in a later work from the same group these reversals were 
 * categorized as "good".
 * 
 * 2) Neutral reversals
 * --------------------
 * 
 * Similar to neutral good reversals, but they occur in oriented components.
 * These reversals affect two black edges belonging to the same cycle (in an 
 * oriented component), which are not "good reversals".
 * A neutral reversal does not change the number of cycles, and does not help
 * to approach the target in future steps.
 * 
 * 3) Bad reversals
 * ----------------
 * 
 * Bad reversals are the reversals that affect two black edges in different 
 * cycles and that are not encompassed in the previous categories 
 * (e.g. non-adjacent hurdles). These reversals move the current genome 
 * away from the target genome, and decrease the number of cycles by 1.
 * 
 * References
 * ----------
 * 
 * For additional information on similar methods that inspired 
 * this implementation, please check:
 * 
 * - Bret Larget et al. "A Bayesian analysis of metazoan mitochondrial 
 * genome arrangements." Molecular Biology and Evolution 22.3 (2005): 486-495.
 * 
 * - Bret Larget et al. "Bayesian phylogenetic inference from animal mitochondrial 
 * genome arrangements." Journal of the Royal Statistical Society Series B: 
 * Statistical Methodology 64.4 (2002): 681-693.
 * 
 * Less related, but also nice:
 * 
 * - Miklós, István, and Aaron E. Darling. "Efficient sampling of parsimonious 
 * inversion histories with application to genome rearrangement in Yersinia." 
 * Genome biology and evolution 1 (2009): 153-164.
 * 
 * Notice that the method implemented here is similar to, but not exactly 
 * the same as the ones mentioned above. For example, reversals that solve 
 * hurdles in an optimal way are considered "good" here, but they would be 
 * considered "bad" in the previous methods.
 * 
 *******************************************************/

#pragma once // It avoids class redefinition.

#include <iostream>
#include <unordered_map>
#include <algorithm> // set_difference
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


// Types of reversals
enum ReversalType : int {
	GOOD = 0,         // GOOD is assigned the value 0
	              // Reversals that increase the number of cycles by 1.
	NEUTRAL_GOOD, // NEUTRAL_GOOD is assigned the value 1
				  // Reversals that do not change the number of cycles, but help to approach the target in future steps.
	NEUTRAL,      // NEUTRAL is assigned the value 2
				  // Reversals that do not change the number of cycles, and do not help to approach the target in future steps.
	BAD,          // BAD is assigned the value 3
				  // Reversals that decrease the number of cycles by 1.
	// Count
};
const int ReversalType_COUNT = 4;

// Types of cycles
enum CycleType : int {
	ORIENTED,     // ORIENTED is assigned the value 0
	UNORIENTED,   // UNORIENTED is assigned the value 1
	TRIVIAL,      // TRIVIAL is assigned the value 2
	HURDLE,       // HURDLE is assigned the value 3
	// Count
};
const int CycleType_COUNT = 4;

/* This class store information about the reversals in a cycle.*/
class CycleCounters {
public:
	std::vector<int> counts; // How many of each type of reversal the cycle has.
	std::vector<std::vector<std::pair<int,int>>> reversals; // List reversals of a certain type in the cycle (if the number is not too big).
	CycleCounters():counts(ReversalType_COUNT, 0),reversals(ReversalType_COUNT){}
	CycleCounters(std::vector<int>& counts_):counts(ReversalType_COUNT, 0),reversals(ReversalType_COUNT){
		for (int i = 0; i < ReversalType_COUNT; ++i) {counts[i] = counts_[i];}
	}
	void clearCounters() {
		for (int i = 0; i < ReversalType_COUNT; ++i) {
			counts[i] = 0;
			reversals[i].clear();
		}
	}
};

/* This class store information about general properties of a cycle.*/
class CycleInfo {
public:
	CycleType type;     // Type of cycle (oriented, hurdle, trivial, etc.)
	int size;           // Nb. of black edges in each cycle.
	int hurdle_idx{-1}; // Index of the hurdle that cycle belongs to.
	std::pair<int,int> hurdles_adj{-1,-1}; // Indices of cycles that are adjacent hurdles (in case the cycle belongs to a hurdle).

	CycleInfo(){}
	CycleInfo(CycleType type_, int size_):type(type_),size(size_){}
	CycleInfo(CycleType type_, int size_, int hurdle_idx_, std::pair<int,int> hurdles_adj_):type(type_),size(size_),hurdle_idx(hurdle_idx_),hurdles_adj(hurdles_adj_){}
};

class SampleRandomReversal {
// private:

public:

	GenomePermutation<BlockSimple> genperm; // Unsigned extended permutation (it assumes that one of the permutations is the identity).
	ConnectedComponents comps;              // It stores the current state of the components.
	std::vector<int> hurdles;               // indices of all roots (cycles) of components that are hurdles.
	std::vector<CycleInfo> cycles_info;     // General properties of each cycle (nb. of black edges, if it is oriented or not, etc.).

	std::vector<CycleCounters> revcounters; // Counts the number of reversals in each category (good, neutral, bad...) for *each cycle*.
	
	std::vector<float> revprobs;  // Probability of each type of reversal.
	std::vector<float> revtotals; // Total count of each type of reversal.

	//std::vector<Reversal> reversals; // It saves all sampled reversals in the solution.

	float p_stop{0.99};

	int nb_cycles{0};
	int nb_gene_ext{0};
	int nb_genes{0};
	int nb_hurdles{0};

	bool debug{false};

	SampleRandomReversal(GenomePermutation<BlockSimple>& genperm, bool debug=false, float p_good=0.97, float p_neutralgood=0.92, float p_neutral=0.92, float p_stop=0.99):genperm(genperm),revprobs(ReversalType_COUNT, 0),revtotals(ReversalType_COUNT, 0.0),p_stop(p_stop),debug(debug){
		initializeComponents();
		initializeCounts();
		std::vector<float> probs = {p_good, p_neutralgood, p_neutral, static_cast<float>(1.0-p_neutral)};
		// For testing:
		// std::vector<float> probs = {0.01, 0.01, 0.9, 0.3};
		initializeRevProbs(probs);
	}

	void initializeComponents();
	void updateComponents();
	void getCyclesInfo();

	void initializeCounts();
	void initializeRevProbs(std::vector<float>& probs);

	// Computes a list of good/proper reversals given one oriented component.
	std::vector<std::pair<int,int>> getGoodReversalsOriented(Cycle& oriented_component);

	float countReversals();
	float countReversalsOrientedComponents();
	float countReversalsUnorientedComponents();

	float countReversalsCycle_trivial(    CycleCounters& counters);
	float countReversalsCycle_unoriented( CycleCounters& counters, const int cycle_size);
	float countReversalsCycle_manyHurdles(CycleCounters& counters, const int cycle_size, const int hurdles_total_size, const int hurdle_cur_size, const int hurdles_adj_size);
	float countReversalsCycle_fewHurdles( CycleCounters& counters, const int cycle_size, const int hurdles_total_size, const int hurdle_cur_size);
	float countReversalsCycle_oneHurdle(  CycleCounters& counters, const int cycle_size, const int hurdles_total_size, const int hurdle_cur_size);
	
	int getSizesAdjacentHurdles(std::vector<int>& sizes_hurdles, const int hurdle_idx);

	std::vector<int> getGeneExtremitiesComponent(const int comp_idx);
	std::vector<int> getAdjacentHurdleExtremities(const int cycle_idx);
	std::vector<int> getAllHurdleExtremities();

	ReversalType sampleReversalType(std::mt19937& rng);
	int  sampleCycle(ReversalType revtype, std::mt19937& rng);
	std::pair<int,int> sampleReversal(std::mt19937& rng, bool updateComps=false);
	std::pair<int,int> sampleReversalFromCycle(const int cycle_idx, const ReversalType revtype, std::mt19937& rng);

	std::pair<int,int> sampleOriented(const int cycle_idx, const ReversalType revtype, std::mt19937& rng);
	std::pair<int,int> sampleUnoriented(const int cycle_idx, const ReversalType revtype, std::mt19937& rng);
	std::pair<int,int> sampleTrivial(const int cycle_idx, const ReversalType revtype, std::mt19937& rng);
	std::pair<int,int> sampleHurdle(const int cycle_idx, const ReversalType revtype, std::mt19937& rng);

	std::vector<std::pair<int,int>> getReversalsFromCycle(const int cycle_idx);
	std::vector<int> getGeneExtNotInCycle(const int cycle_idx);

	std::vector<std::pair<int,int>> makeAllCombinations(std::vector<int> elems);
	std::vector<std::pair<int,int>> makeAllCombinations(std::vector<int> elems1, std::vector<int> elems2);
	
};

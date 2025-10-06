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
 * Similar to the approach implemented here, Miklos and Darling (2009) 
 * consider "hurdle-cutting" and "hurdle-merging" operations as "good 
 * operations", since these operations also compose sorting scenarios.
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
 * the same as the ones mentioned above. For example, reversals that solve 
 * hurdles in an optimal way are considered "good" here, but they would be 
 * considered "bad" in the previous methods.
 * 
 *******************************************************/

#include <iostream>
#include <unordered_map>
#include <vector>
#include <stack>
#include <list>
#include <cstdlib> // exit
#include <cmath>	 // abs
#include <memory>	 // shared_ptr

#include "sampleReversal_York2002.hpp"

void ReversalSampler::initializeComponents(){

	comps = ConnectedComponents(genperm.getUnsignedExtendedPerm(),debug);

	// Update stats.
	nb_cycles    = comps.forest.size();
	nb_gene_ext  = comps.perm.size();
	nb_genes     = nb_gene_ext/2;

	// Find hurdles.
	UnorientedComponents comps_unoriented = UnorientedComponents(genperm, comps);
	hurdles = comps_unoriented.findHurdles(); // indices of all roots (cycles) of components that are hurdles.
	nb_hurdles = hurdles.size();

	// Find information about each cycle (nb. of black edges, type of component that the cycle belongs to).
	getCyclesInfo();
}

void ReversalSampler::getCyclesInfo(){

	// Clean data structure.
	cycles_info.clear();
	cycles_info.resize(nb_cycles);

	// Number of black edges in each cycle (cycle size).
	std::vector<int> cycles_sizes = comps.getCycleSizes(true);
	
	int cur_hurdle_idx = (nb_hurdles > 0) ? 0 : -1;
	int cur_hurdle     = (nb_hurdles > 0) ? hurdles[0] : -1;

	// Find type of each cycle (oriented, unoriented, trivial, or hurdle).
	for (int const& root_idx : comps.rootList) {

		std::vector<int> cycles = comps.getCycles(comps.forest[root_idx]);
		std::pair<int, int> hurdles_adj{-1, -1};
		int hurdle_idx = -1;

		// Get type of component (all cycles in the component have the same type).
		CycleType type = CycleType::ORIENTED;
		if(!comps.forest[root_idx].oriented){
			if(comps.forest[root_idx].genes.size() == 2) {	
				type = CycleType::TRIVIAL;
			} else if (root_idx == cur_hurdle) {
				type = CycleType::HURDLE;
				hurdle_idx = root_idx;
				// Get adjacent hurdles.
				if(hurdles.size() > 1){
					int prev_hurdle = (cur_hurdle_idx == 0) ? hurdles[hurdles.size()-1] : hurdles[cur_hurdle_idx-1];
					int next_hurdle = (cur_hurdle_idx == (hurdles.size()-1)) ? hurdles[0] : hurdles[cur_hurdle_idx+1];
					hurdles_adj = std::make_pair(prev_hurdle,next_hurdle);
				}
				// Update hurdle.
				cur_hurdle_idx += 1;
				cur_hurdle = (nb_hurdles > cur_hurdle_idx) ? hurdles[cur_hurdle_idx] : -1;
			} else {
				type = CycleType::UNORIENTED;
			}
		}

		// Update info of each cycle (cycle type and cycle size).
		for (int const& cycle_idx : cycles) {
			cycles_info[cycle_idx].type = type;
			cycles_info[cycle_idx].size = cycles_sizes[cycle_idx];
			cycles_info[cycle_idx].hurdles_adj = hurdles_adj;
			cycles_info[cycle_idx].hurdle_idx  = hurdle_idx;
		}
	}
}

void ReversalSampler::updateComponents(){
	// Update connected components.
	initializeComponents();
	// Clear reversal counters.
	for (int cycle_idx = 0; cycle_idx < nb_cycles; ++cycle_idx) {rev_counters[cycle_idx].clearCounters();}
	for (int revtype_idx = 0; revtype_idx < ReversalType_COUNT; ++revtype_idx) {rev_totals[revtype_idx] = 0;}
}

void ReversalSampler::initializeCounts(){
	// Initialize structures that store the counts for good, neutral good, neutral, and bad reversals.
	// The counts store information for *each cycle* (not per component).
	rev_counters.resize(nb_cycles);
}

void ReversalSampler::initializeRevProbs(std::vector<float>& probs){
	rev_weights_total = 0;
	for (int i=0; i<ReversalType_COUNT; ++i) {
		rev_weights[i]     = probs[i];
		rev_weights_total += probs[i];
	}
	// Cumulative weights.
	rev_weights_cum[0] = rev_weights[0];
	for (int i=1; i<ReversalType_COUNT; ++i) {
		rev_weights_cum[i] = rev_weights_cum[i-1] + rev_weights[i];
	}
}

void ReversalSampler::updateRevProbs() {
	// Make sure that reversal types without any 
	// reversals associated have 0 chance to be sampled.
	float rev_weight_cur = (rev_totals[0] == 0) ? 0.0 : rev_weights[0];
	rev_weights_cum[0]   = rev_weight_cur;
	rev_weights_total    = rev_weight_cur;
	for (int i=1; i<ReversalType_COUNT; ++i) {
		rev_weight_cur     = (rev_totals[i] == 0) ? 0.0 : rev_weights[i];
		rev_weights_cum[i] = rev_weights_cum[i-1] + rev_weight_cur;
	}
}

// It returns a list of pairs of gene extremities. 
// Each pair corresponds to one good reversal.
// Each pair contains two gene extremities from the 
// unsigned extended permutation representation.
std::vector<std::pair<int,int>> ReversalSampler::getGoodReversalsOriented(Cycle& oriented_component){
	
	std::vector<std::pair<int,int>> goodReversals;

	// Permutation **must** start at 1: [1 2 .. gene]
	std::unordered_map<int,std::pair<int,int>> newlabels_map;
	std::vector<int> perm = genperm.getExtendedPerm(oriented_component.genes,newlabels_map);
	const int n = perm.size();

	// Compute position of each gene in this sub-permutation.
	std::vector<int> genepos(n);
	for (int pos=0; pos<n; pos++) {
		const int g_id  = perm[pos];
		const int g_idx = std::abs(g_id)-1;
		const int sign  = (g_id > 0 ? 1 : -1);
		genepos[g_idx] = pos*sign;
	}

	// Find proper/good reversals.
	for (int g_idx=0; g_idx<(n-1); g_idx++) {
		// Good reversal: genes have different signs.
		const bool isGoodReversal = ((genepos[g_idx] < 0) != (genepos[g_idx+1] < 0));
		if(isGoodReversal){
			const int gene_id   = g_idx+1;
			int uns_ext_1 = newlabels_map[gene_id].second;
			int uns_ext_2 = newlabels_map[gene_id+1].first;

			// Get the start extremity of the black edges involved in the reversal.
			if(comps.idxs[uns_ext_1] % 2 != 0) {uns_ext_1 = comps.perm[comps.idxs[uns_ext_1]-1];}
			if(comps.idxs[uns_ext_2] % 2 != 0) {uns_ext_2 = comps.perm[comps.idxs[uns_ext_2]-1];}

			// Make sure that lowest position always comes first.
			if(comps.idxs[uns_ext_1] > comps.idxs[uns_ext_2]){std::swap(uns_ext_1,uns_ext_2);}

			// Extremities should be in order. Check if they already exist in the list.
			std::pair<int,int> rev(uns_ext_1,uns_ext_2);
			if (std::find(goodReversals.begin(), goodReversals.end(), rev) == goodReversals.end()){
				goodReversals.emplace_back(rev);
				std::cout << "Good reversal (gene " << (g_idx+1) << "): [" << uns_ext_1 << ", " << uns_ext_2 << "] - ADDED" << std::endl;
			} else {
				std::cout << "Good reversal (gene " << (g_idx+1) << "): [" << uns_ext_1 << ", " << uns_ext_2 << "] - IGNORED" << std::endl;
			}
		}
	}
	return goodReversals;
}

void ReversalSampler::debugPrintRevTotals(){
	std::cout << " good = " << rev_totals[ReversalType::GOOD] << "; neutral good = " << rev_totals[ReversalType::NEUTRAL_GOOD] << "; neutral = " << rev_totals[ReversalType::NEUTRAL] << "; bad = " << rev_totals[ReversalType::BAD] << std::endl;
}

float ReversalSampler::countReversalsOrientedComponents(){
	// Find good, neutral, and bad reversals.
	float obs_nb_reversals = 0;
	for (int const& root_idx : comps.rootList) {
		// Component is oriented.
		if(comps.forest[root_idx].oriented){

			// Compute good reversals.
			std::vector<std::pair<int,int>> goodReversals = getGoodReversalsOriented(comps.forest[root_idx]);

			// Update counts of good reversals per cycle.
			for (const std::pair<int, int>& goodReversal : goodReversals) {

				// Get the index of the cycle that contains the gene extremity.
				int const cycle_idx = comps.getCycleIdx(goodReversal.first);
				rev_counters[cycle_idx].counts[ReversalType::GOOD] += 1;
				rev_counters[cycle_idx].reversals[ReversalType::GOOD].emplace_back(goodReversal);
				
				// Check if values make sense (a good reversal should always happen inside a single cycle).
				if(comps.getCycleIdx(goodReversal.first) != comps.getCycleIdx(goodReversal.second)) {
					std::cout << "ERROR! A reversal identified as 'good' affects two different cycles (cycle indices " << comps.getCycleIdx(goodReversal.first) << " and "<< comps.getCycleIdx(goodReversal.second) << ").\nProgram is aborting." << std::endl;
					exit(1);
				}
			}
			// Get all cycles belonging to the component.
			std::vector<int> cycles = comps.getCycles(comps.forest[root_idx]);

			// Compute neutral good, neutral and bad reversals.
			for (int const& cycle_idx : cycles) {
				const int cycle_size      = cycles_info[cycle_idx].size;
				
				rev_counters[cycle_idx].counts[ReversalType::NEUTRAL_GOOD] = 0;
				rev_counters[cycle_idx].counts[ReversalType::NEUTRAL]      = ((cycle_size*(cycle_size-1))/2) - rev_counters[cycle_idx].counts[ReversalType::GOOD];
				rev_counters[cycle_idx].counts[ReversalType::BAD]          = cycle_size*(nb_genes-cycle_size);

				rev_totals[ReversalType::GOOD]         += rev_counters[cycle_idx].counts[ReversalType::GOOD];
				rev_totals[ReversalType::NEUTRAL_GOOD] += rev_counters[cycle_idx].counts[ReversalType::NEUTRAL_GOOD];
				rev_totals[ReversalType::NEUTRAL]      += rev_counters[cycle_idx].counts[ReversalType::NEUTRAL];
				rev_totals[ReversalType::BAD]          += rev_counters[cycle_idx].counts[ReversalType::BAD]/2.0;

				// Update count of observed reversals.
				// As bad reversals are counted twice in two different cycles, their value is divided by 2.
				obs_nb_reversals += (  rev_counters[cycle_idx].counts[ReversalType::GOOD]
											+ rev_counters[cycle_idx].counts[ReversalType::NEUTRAL_GOOD] 
											+ rev_counters[cycle_idx].counts[ReversalType::NEUTRAL] 
											+ rev_counters[cycle_idx].counts[ReversalType::BAD]/2.0);
			}
		}
	}
	return obs_nb_reversals;
}

// Trivial component (i.e. 2 gene extremities that are already sorted).
// There are no neutral or good reversals in this case (the component is already sorted).
float ReversalSampler::countReversalsCycle_trivial(CycleCounters& counters) {
	counters.counts[ReversalType::GOOD]         = 0;
	counters.counts[ReversalType::NEUTRAL_GOOD] = 0;
	counters.counts[ReversalType::NEUTRAL]      = 0;
	counters.counts[ReversalType::BAD]          = (nb_genes-1);

	rev_totals[ReversalType::GOOD]         += 0;
	rev_totals[ReversalType::NEUTRAL_GOOD] += 0;
	rev_totals[ReversalType::NEUTRAL]      += 0;
	rev_totals[ReversalType::BAD]          += counters.counts[ReversalType::BAD]/2.0;

	// Return count of observed reversals.
	// As bad reversals are counted twice in two different cycles, their value is divided by 2.
	return (counters.counts[ReversalType::GOOD] 
			+ counters.counts[ReversalType::NEUTRAL_GOOD] 
			+ counters.counts[ReversalType::NEUTRAL] 
			+ counters.counts[ReversalType::BAD]/2.0);
}

// "Normal" unoriented component (not a hurdle, not trivial).
float ReversalSampler::countReversalsCycle_unoriented(CycleCounters& counters, const int cycle_size) {
	counters.counts[ReversalType::GOOD]         = 0;
	counters.counts[ReversalType::NEUTRAL_GOOD] = ((cycle_size*(cycle_size-1))/2);
	counters.counts[ReversalType::NEUTRAL]      = 0;
	counters.counts[ReversalType::BAD]          = cycle_size*(nb_genes-cycle_size);

	rev_totals[ReversalType::GOOD]         += 0;
	rev_totals[ReversalType::NEUTRAL_GOOD] += counters.counts[ReversalType::NEUTRAL_GOOD];
	rev_totals[ReversalType::NEUTRAL]      += 0;
	rev_totals[ReversalType::BAD]          += counters.counts[ReversalType::BAD]/2.0;

	return (counters.counts[ReversalType::GOOD] 
			+ counters.counts[ReversalType::NEUTRAL_GOOD] 
			+ counters.counts[ReversalType::NEUTRAL] 
			+ counters.counts[ReversalType::BAD]/2.0);
}

// Hurdles - General case (i.e. #hurdles > 3).
// Good reversals: a black edge from the cycle and a black edge from 
//                 any other cycle belonging to a non-adjacent hurdle.
// Good neutral reversals: 
//    (1) Two black edges from the cycle; 
//    (2) A black edge from the cycle and another black edge belonging to its own hurdle or an adjacent hurdle.
float ReversalSampler::countReversalsCycle_manyHurdles(CycleCounters& counters, const int cycle_size, 
							const int hurdles_total_size, const int hurdle_cur_size, const int hurdles_adj_size) {

	const int black_edges_same_cycle = (cycle_size*(cycle_size-1))/2;
	const int black_edges_own_hurdle = cycle_size*(hurdle_cur_size-cycle_size);
	const int black_edges_adj_hurdle = cycle_size*hurdles_adj_size;

	counters.counts[ReversalType::GOOD]         = cycle_size*(hurdles_total_size-(hurdles_adj_size+hurdle_cur_size));
	counters.counts[ReversalType::NEUTRAL_GOOD] = black_edges_same_cycle + black_edges_own_hurdle + black_edges_adj_hurdle;
	counters.counts[ReversalType::NEUTRAL]      = 0;
	counters.counts[ReversalType::BAD]          = cycle_size*(nb_genes-hurdles_total_size);

	rev_totals[ReversalType::GOOD]         += counters.counts[ReversalType::GOOD]/2.0;
	rev_totals[ReversalType::NEUTRAL_GOOD] += (black_edges_same_cycle + (black_edges_own_hurdle + black_edges_adj_hurdle)/2.0);
	rev_totals[ReversalType::NEUTRAL]      += 0;
	rev_totals[ReversalType::BAD]          += counters.counts[ReversalType::BAD]/2.0;

	return (counters.counts[ReversalType::GOOD]/2.0 
			+ black_edges_same_cycle + (black_edges_own_hurdle + black_edges_adj_hurdle)/2.0 
			+ counters.counts[ReversalType::NEUTRAL] 
			+ counters.counts[ReversalType::BAD]/2.0);
}

// Hurdle - Base case.
// 2 or 3 hurdles. All hurdles are adjacent to each other.
// Good reversals: a black edge from the cycle and a black edge from the other hurdle.
// Good neutral reversals: 
//    (1) Two black edges from the cycle; 
//    (2) A black edge from the cycle and another black edge belonging to its own hurdle.
float ReversalSampler::countReversalsCycle_fewHurdles(CycleCounters& counters, const int cycle_size, 
							const int hurdles_total_size, const int hurdle_cur_size) {

	const int black_edges_same_cycle = (cycle_size*(cycle_size-1))/2;
	const int black_edges_own_hurdle = cycle_size*(hurdle_cur_size-cycle_size);

	counters.counts[ReversalType::GOOD]         = cycle_size*(hurdles_total_size-hurdle_cur_size);
	counters.counts[ReversalType::NEUTRAL_GOOD] = black_edges_same_cycle + black_edges_own_hurdle;
	counters.counts[ReversalType::NEUTRAL]      = 0;
	counters.counts[ReversalType::BAD]          = cycle_size*(nb_genes-hurdles_total_size);

	rev_totals[ReversalType::GOOD]         += counters.counts[ReversalType::GOOD]/2.0;
	rev_totals[ReversalType::NEUTRAL_GOOD] += (black_edges_same_cycle + black_edges_own_hurdle/2.0);
	rev_totals[ReversalType::NEUTRAL]      += 0;
	rev_totals[ReversalType::BAD]          += counters.counts[ReversalType::BAD]/2.0;

	return (counters.counts[ReversalType::GOOD]/2.0 
			+ black_edges_same_cycle + black_edges_own_hurdle/2.0 
			+ counters.counts[ReversalType::NEUTRAL] 
			+ counters.counts[ReversalType::BAD]/2.0);
}

// Hurdle - Base case.
// A single hurdle.
// Good reversals:
//    (1) Two black edges from the cycle; 
//    (2) A black edge from the cycle and another black edge belonging to its own hurdle.
float ReversalSampler::countReversalsCycle_oneHurdle(CycleCounters& counters, const int cycle_size, 
							const int hurdles_total_size, const int hurdle_cur_size) {

	const int black_edges_same_cycle = (cycle_size*(cycle_size-1))/2;
	const int black_edges_own_hurdle = cycle_size*(hurdle_cur_size-cycle_size);

	counters.counts[ReversalType::GOOD]         = black_edges_same_cycle + black_edges_own_hurdle;
	counters.counts[ReversalType::NEUTRAL_GOOD] = 0;
	counters.counts[ReversalType::NEUTRAL]      = 0;
	counters.counts[ReversalType::BAD]          = cycle_size*(nb_genes-hurdles_total_size);

	rev_totals[ReversalType::GOOD]         += black_edges_same_cycle + black_edges_own_hurdle/2.0;
	rev_totals[ReversalType::NEUTRAL_GOOD] += 0;
	rev_totals[ReversalType::NEUTRAL]      += 0;
	rev_totals[ReversalType::BAD]          += counters.counts[ReversalType::BAD]/2.0;

	return (black_edges_same_cycle + black_edges_own_hurdle/2.0 
			+ counters.counts[ReversalType::NEUTRAL_GOOD]
			+ counters.counts[ReversalType::NEUTRAL] 
			+ counters.counts[ReversalType::BAD]/2.0);
}

int ReversalSampler::getSizesAdjacentHurdles(std::vector<int>& sizes_hurdles, const int hurdle_idx) {
	// Next hurdle.
	int hurdles_adj_size = (sizes_hurdles.size() > (hurdle_idx+1)) ? sizes_hurdles[hurdle_idx+1] : sizes_hurdles[0];
	// Previous hurdle.
	hurdles_adj_size = (-1 < (hurdle_idx-1)) ? sizes_hurdles[hurdle_idx-1] : sizes_hurdles[sizes_hurdles.size()-1];
	return hurdles_adj_size;
}

float ReversalSampler::countReversalsUnorientedComponents(){

	// Find total edges in hurdles.
	std::vector<int> sizes_hurdles;
	int hurdles_total_size = 0;
	if(nb_hurdles > 0){
		for (int root_idx : hurdles) {
			int hurdle_size = 0;
			std::vector<int> cycles = comps.getCycles(comps.forest[root_idx]);
			for (int const& cycle_idx : cycles) {hurdle_size += cycles_info[cycle_idx].size;}
			sizes_hurdles.emplace_back(hurdle_size);
			hurdles_total_size += hurdle_size;
		}
	}

	float obs_nb_reversals = 0;

	int next_hurdle_idx = (nb_hurdles > 0) ? 0 : -1;
	int next_hurdle = (nb_hurdles > 0) ? hurdles[0] : -1;

	// Find good, neutral, and bad reversals.
	for (int const& root_idx : comps.rootList) {

		// Skip if component is not unoriented.
		if(comps.forest[root_idx].oriented){continue;}

		// Check if component is trivial.
		if(comps.forest[root_idx].genes.size() == 2) {
			obs_nb_reversals += countReversalsCycle_trivial(rev_counters[root_idx]);

		// Check if component is a hurdle.
		} else if (root_idx == next_hurdle) {

			std::vector<int> cycles   = comps.getCycles(comps.forest[root_idx]);
			const int hurdle_cur_size = sizes_hurdles[next_hurdle_idx];

			// Hurdle - General case (#hurdles > 3).
			if(nb_hurdles > 3){ 
				// Get sizes of adjacent hurdles.
				int hurdles_adj_size = getSizesAdjacentHurdles(sizes_hurdles, next_hurdle_idx);
				for (int const& cycle_idx : cycles) {
					obs_nb_reversals += countReversalsCycle_manyHurdles(rev_counters[cycle_idx], cycles_info[cycle_idx].size, hurdles_total_size, hurdle_cur_size, hurdles_adj_size);
				}

			// Hurdle - Base case (#hurdles = 2 or 3).
			} else if(nb_hurdles > 1){
				for (int const& cycle_idx : cycles) {
					obs_nb_reversals += countReversalsCycle_fewHurdles(rev_counters[cycle_idx], cycles_info[cycle_idx].size, hurdles_total_size, hurdle_cur_size);
				}

			// Hurdle - Base case (#hurdles = 1).
			} else {
				for (int const& cycle_idx : cycles) {
					obs_nb_reversals += countReversalsCycle_oneHurdle(rev_counters[cycle_idx], cycles_info[cycle_idx].size, hurdles_total_size, hurdle_cur_size);
				}
			}

			// Update hurdle.
			next_hurdle_idx += 1;
			next_hurdle = (nb_hurdles > next_hurdle_idx) ? hurdles[next_hurdle_idx] : -1;

		// "Normal" unoriented component.
		} else {
			std::vector<int> cycles = comps.getCycles(comps.forest[root_idx]);
			for (int const& cycle_idx : cycles) {
				obs_nb_reversals += countReversalsCycle_unoriented(rev_counters[cycle_idx], cycles_info[cycle_idx].size);
			}
		}
	}
	return obs_nb_reversals;
}

float ReversalSampler::countReversals(){

	const int exp_nb_reversals = (nb_genes*(nb_genes-1))/2;

	float obs_nb_reversals = countReversalsOrientedComponents();
	obs_nb_reversals    += countReversalsUnorientedComponents();
	
	std::cout << " - Total nb. of reversals : Obs=" << obs_nb_reversals << "; Exp=" << exp_nb_reversals << std::endl;
	if(static_cast<int>(obs_nb_reversals) != exp_nb_reversals){
		std::cout << "ERROR! The number of observed reversals differs from the expected number of reversals (Obs.=" << obs_nb_reversals << "; Exp.="<< exp_nb_reversals << ").\nProgram is aborting." << std::endl;
		exit(1);
	}
	return obs_nb_reversals;
}

ReversalType ReversalSampler::sampleReversalType(std::mt19937& rng){
	ReversalType revtype = ReversalType::GOOD;
	// Make sure that reversal types without any 
	// reversals associated have 0 chance to be sampled.
	updateRevProbs();
	// Sample reversal type.
	std::uniform_real_distribution distr(0.0, static_cast<double>(rev_weights_total)); // [0, rev_weights_total)
	float rdmval = distr(rng);
	// Find the index corresponding to the random value
	for (int revtype_idx = 0; revtype_idx < ReversalType_COUNT; ++revtype_idx) {
		if (rdmval < rev_weights_cum[revtype_idx]) {
			revtype = static_cast<ReversalType>(revtype_idx);
			break;
		}
	}
	return revtype;
}

int ReversalSampler::sampleCycle(ReversalType revtype, std::mt19937& rng){
	const int weight_tot = rev_totals[revtype];
	// Sample weight.
	std::uniform_int_distribution distr(0, weight_tot - 1);
	const int weight_rdm = distr(rng);
	// Sample cycle based on the random weight.
	int cycle_idx_sampled = 0;
	int weight_cum = 0;
	for (int cycle_idx = 0; cycle_idx < nb_cycles; ++cycle_idx) {
		weight_cum += rev_counters[cycle_idx].counts[revtype];
		if(weight_rdm < weight_cum){
			cycle_idx_sampled = cycle_idx;
			break;
		}
	}
	return cycle_idx_sampled;
}


std::vector<std::pair<int,int>> ReversalSampler::makeAllCombinations(std::vector<int> elems){
	std::vector<std::pair<int,int>> allpairs;
	for (int elem_i = 0; elem_i < elems.size(); ++elem_i) {
		for (int elem_j = elem_i+1; elem_j < elems.size(); ++elem_j) {
			// Extremities are sorted by the position they appear in the unsigned permutation.
			if(comps.idxs[elems[elem_i]]>comps.idxs[elems[elem_j]]) {
				allpairs.emplace_back(std::make_pair(elems[elem_j],elems[elem_i]));	
			} else {
				allpairs.emplace_back(std::make_pair(elems[elem_i],elems[elem_j]));
			}
		}
	}
	return allpairs;
}

std::vector<std::pair<int,int>> ReversalSampler::makeAllCombinations(std::vector<int> elems1, std::vector<int> elems2){
	std::vector<std::pair<int,int>> allpairs;
	for (int elem_i = 0; elem_i < elems1.size(); ++elem_i) {
		for (int elem_j = 0; elem_j < elems2.size(); ++elem_j) {
			// Extremities are sorted by the position they appear in the unsigned permutation.
			if(comps.idxs[elems1[elem_i]]>comps.idxs[elems2[elem_j]]) {
				allpairs.emplace_back(std::make_pair(elems2[elem_j],elems1[elem_i]));	
			} else {
				allpairs.emplace_back(std::make_pair(elems1[elem_i],elems2[elem_j]));
			}
		}
	}
	return allpairs;
}

std::vector<std::pair<int,int>> ReversalSampler::getReversalsFromCycle(const int cycle_idx){
	std::vector<int> genes = comps.getGeneExtremities(comps.forest[cycle_idx], true);
	return makeAllCombinations(genes);
}

std::vector<int> ReversalSampler::getGeneExtNotInCycle(const int cycle_idx){

	// All gene extremities.
	std::vector<int> gene_exts_all = comps.getGeneExtremities(true);

	// Gene extremities in the cycle.
	std::vector<int> gene_exts_cyc = comps.getGeneExtremities(comps.forest[cycle_idx], true);

	std::vector<int> gene_exts_diff;
	std::set_difference(gene_exts_all.begin(), gene_exts_all.end(), gene_exts_cyc.begin(), gene_exts_cyc.end(),
						std::inserter(gene_exts_diff, gene_exts_diff.begin()));
	return gene_exts_diff;
}

std::pair<int,int> ReversalSampler::sampleOriented(const int cycle_idx, const ReversalType revtype, std::mt19937& rng){

	// rev_counters[cycle_idx].counts[ReversalType::NEUTRAL_GOOD] = 0;
	// rev_counters[cycle_idx].counts[ReversalType::NEUTRAL]      = ((cycle_size*(cycle_size-1))/2) - rev_counters[cycle_idx].counts[ReversalType::GOOD];
	// rev_counters[cycle_idx].counts[ReversalType::BAD]          = cycle_size*(nb_genes-cycle_size);
	//std::vector<std::pair<int,int>>& goodReversals = rev_counters[cycle_idx].reversals[ReversalType::GOOD];

	std::pair<int, int> reversal{-1, -1};
	switch(revtype) {
		case ReversalType::GOOD: {
			std::cout << " --- Good" << std::endl;
			std::vector<std::pair<int,int>>& goodReversals = rev_counters[cycle_idx].reversals[ReversalType::GOOD];
			std::uniform_int_distribution distr(0, static_cast<int>(goodReversals.size()) - 1);
			const int rdmidx = distr(rng);
			reversal = goodReversals[rdmidx];
			break;
		}
		case ReversalType::NEUTRAL_GOOD: {
			std::cout << "ERROR! Oriented cycles do not have 'neutral good' reversals." << std::endl;
			break;
		}
		case ReversalType::NEUTRAL: {
			std::cout << " --- Neutral" << std::endl;

			std::vector<std::pair<int,int>>  allReversals  = getReversalsFromCycle(cycle_idx);
			std::vector<std::pair<int,int>>& goodReversals = rev_counters[cycle_idx].reversals[ReversalType::GOOD];

			std::cout << "\t\t Good reversals (" << goodReversals.size() << ") = " ;
			for (std::pair<int, int>& r : goodReversals) {
				std::cout << " (" << r.first << ", " << r.second << ")";
			}
			std::cout << std::endl;

			std::vector<std::pair<int,int>> neutralReversals;
			std::set_difference(allReversals.begin(), allReversals.end(), goodReversals.begin(), goodReversals.end(),
								std::inserter(neutralReversals, neutralReversals.begin()));

			std::uniform_int_distribution distr(0, static_cast<int>(neutralReversals.size()) - 1);
			const int rdmidx = distr(rng);
			reversal = neutralReversals[rdmidx];
			break;
		}
		case ReversalType::BAD: {
			std::cout << " --- Bad" << std::endl;
			
			std::vector<std::pair<int,int>> badReversals = makeAllCombinations(comps.getGeneExtremities(comps.forest[cycle_idx], true), getGeneExtNotInCycle(cycle_idx));
			std::uniform_int_distribution distr(0, static_cast<int>(badReversals.size()) - 1);
			const int rdmidx = distr(rng);
			reversal = badReversals[rdmidx];
			break;
		}
	}

	return reversal;
}

std::pair<int,int> ReversalSampler::sampleUnoriented(const int cycle_idx, const ReversalType revtype, std::mt19937& rng){

	// UNORIENTED
	// counters.counts[ReversalType::GOOD]         = 0;
	// counters.counts[ReversalType::NEUTRAL_GOOD] = ((cycle_size*(cycle_size-1))/2);
	// counters.counts[ReversalType::NEUTRAL]      = 0;
	// counters.counts[ReversalType::BAD]          = cycle_size*(nb_genes-cycle_size);

	std::pair<int, int> reversal{-1, -1};
	switch(revtype) {
		case ReversalType::GOOD: {
			std::cout << "ERROR! Unoriented cycles do not have 'good' reversals." << std::endl;
			break;
		}
		case ReversalType::NEUTRAL_GOOD: {
			std::cout << "\t It is neutral good" << std::endl;
			std::vector<std::pair<int,int>> allReversals = getReversalsFromCycle(cycle_idx);
			std::uniform_int_distribution distr(0, static_cast<int>(allReversals.size()) - 1);
			const int rdmidx = distr(rng);
			reversal = allReversals[rdmidx];
			break;
		}
		case ReversalType::NEUTRAL: {
			std::cout << "ERROR! Unoriented cycles do not have 'neutral' reversals." << std::endl;
			break;
		}
		case ReversalType::BAD: {
			std::cout << " --- Bad " << std::endl;
			std::vector<std::pair<int,int>> badReversals = makeAllCombinations(comps.getGeneExtremities(comps.forest[cycle_idx], true), getGeneExtNotInCycle(cycle_idx));
			std::uniform_int_distribution distr(0, static_cast<int>(badReversals.size()) - 1);
			const int rdmidx = distr(rng);
			reversal = badReversals[rdmidx];
			break;
		}
	}
	return reversal;
}

std::pair<int,int> ReversalSampler::sampleTrivial(const int cycle_idx, const ReversalType revtype, std::mt19937& rng){

	// TRIVIAL
	// counters.counts[ReversalType::GOOD]         = 0;
	// counters.counts[ReversalType::NEUTRAL_GOOD] = 0;
	// counters.counts[ReversalType::NEUTRAL]      = 0;
	// counters.counts[ReversalType::BAD]          = (nb_genes-1);

	std::pair<int, int> reversal{-1, -1};
	switch(revtype) {
		case ReversalType::GOOD: {
			std::cout << "ERROR! Trivial cycles do not have 'good' reversals." << std::endl;
			break;
		}
		case ReversalType::NEUTRAL_GOOD: {
			std::cout << "ERROR! Trivial cycles do not have 'neutral good' reversals." << std::endl;
			break;
		}
		case ReversalType::NEUTRAL: {
			std::cout << "ERROR! Trivial cycles do not have 'neutral' reversals." << std::endl;
			break;
		}
		case ReversalType::BAD: {
			std::cout << " --- Bad" << std::endl;
			std::vector<std::pair<int,int>> badReversals = makeAllCombinations(comps.getGeneExtremities(comps.forest[cycle_idx], true), getGeneExtNotInCycle(cycle_idx));
			std::uniform_int_distribution distr(0, static_cast<int>(badReversals.size()) - 1);
			const int rdmidx = distr(rng);
			reversal = badReversals[rdmidx];
			break;
		}
	}
	return reversal;
}

std::vector<int> ReversalSampler::getGeneExtremitiesComponent(const int comp_idx) {
	std::vector<int> all_comp_exts = comps.forest[comp_idx].genes;
	std::vector<int> sel_comp_exts;
	for(int g : all_comp_exts){if(comps.idxs[g] % 2 == 0){sel_comp_exts.emplace_back(g);}}
	return sel_comp_exts;
}

std::vector<int> ReversalSampler::getAdjacentHurdleExtremities(const int cycle_idx){
	std::vector<int> adjhurdle_exts;
	if(hurdles.size() == 1) {return adjhurdle_exts;} 

	// At least two hurdles.
	int adjhurdle_idx = cycles_info[cycle_idx].hurdles_adj.first;
	adjhurdle_exts    = getGeneExtremitiesComponent(adjhurdle_idx);

	// At least three hurdles.
	if(hurdles.size() > 2) {
		adjhurdle_idx = cycles_info[cycle_idx].hurdles_adj.second;
		std::vector<int> adjhurdle_exts_extra = getGeneExtremitiesComponent(adjhurdle_idx);
		adjhurdle_exts.insert(adjhurdle_exts.end(), adjhurdle_exts_extra.begin(), adjhurdle_exts_extra.end());
	}
	return adjhurdle_exts;
}

std::vector<int> ReversalSampler::getAllHurdleExtremities(){
	std::vector<int> hurdle_exts_all;
	for (int hurdle_idx : hurdles) {
		std::vector<int> hurdle_exts = getGeneExtremitiesComponent(hurdle_idx);
		hurdle_exts_all.insert(hurdle_exts_all.end(), hurdle_exts.begin(), hurdle_exts.end());
	}
	return hurdle_exts_all;
}

std::pair<int,int> ReversalSampler::sampleHurdle(const int cycle_idx, const ReversalType revtype, std::mt19937& rng){

	// MANY HURDLES
	// counters.counts[ReversalType::GOOD]         = cycle_size*(hurdles_total_size-(hurdles_adj_size+hurdle_cur_size));
	// counters.counts[ReversalType::NEUTRAL_GOOD] = black_edges_same_cycle + black_edges_own_hurdle + black_edges_adj_hurdle;
	// counters.counts[ReversalType::NEUTRAL]      = 0;
	// counters.counts[ReversalType::BAD]          = cycle_size*(nb_genes-hurdles_total_size);

	// FEW HURDLES
	// counters.counts[ReversalType::GOOD]         = cycle_size*(hurdles_total_size-hurdle_cur_size);
	// counters.counts[ReversalType::NEUTRAL_GOOD] = black_edges_same_cycle + black_edges_own_hurdle;
	// counters.counts[ReversalType::NEUTRAL]      = 0;
	// counters.counts[ReversalType::BAD]          = cycle_size*(nb_genes-hurdles_total_size);

	// ONE HURDLE
	// counters.counts[ReversalType::GOOD]         = black_edges_same_cycle + black_edges_own_hurdle;
	// counters.counts[ReversalType::NEUTRAL_GOOD] = 0;
	// counters.counts[ReversalType::NEUTRAL]      = 0;
	// counters.counts[ReversalType::BAD]          = cycle_size*(nb_genes-hurdles_total_size);

	std::vector<int> cycle_exts      = comps.getGeneExtremities(comps.forest[cycle_idx], true);
	std::vector<int> hurdle_exts_own = getGeneExtremitiesComponent(cycles_info[cycle_idx].hurdle_idx);
	std::vector<int> hurdle_exts_adj = getAdjacentHurdleExtremities(cycle_idx);
	std::vector<int> hurdle_exts_all = getAllHurdleExtremities();

	std::pair<int, int> reversal{-1, -1};
	switch(revtype) {
		case ReversalType::GOOD: {
			std::cout << " --- Good ";
			std::vector<std::pair<int,int>> goodReversals;

			// Case: Many hurdles (h > 3).
			if(hurdles.size() > 3){
				std::cout << " --- Many hurdles ";

				std::vector<int> hurdle_exts_notadj;
				std::set_difference(hurdle_exts_all.begin(), hurdle_exts_all.end(), hurdle_exts_adj.begin(), hurdle_exts_adj.end(),
									std::inserter(hurdle_exts_notadj, hurdle_exts_notadj.begin()));
				goodReversals = makeAllCombinations(cycle_exts, hurdle_exts_notadj);

			// Case: Few hurdles (h=2 or 3).
			} else if(hurdles.size() > 1){ 
				std::cout << " --- Few hurdles ";

				goodReversals = makeAllCombinations(cycle_exts, hurdle_exts_adj);

			// Case: Single hurdle (h=1).
			} else {
				std::cout << " --- Single hurdle ";

				std::vector<int> hurdle_exts_own_diff;
				std::set_difference(hurdle_exts_own.begin(), hurdle_exts_own.end(), cycle_exts.begin(), cycle_exts.end(),
									std::inserter(hurdle_exts_own_diff, hurdle_exts_own_diff.begin()));

				goodReversals = makeAllCombinations(cycle_exts, hurdle_exts_own_diff);
				std::vector<std::pair<int,int>> cycReversals = getReversalsFromCycle(cycle_idx);
				goodReversals.insert(goodReversals.end(), cycReversals.begin(), cycReversals.end());
			}

			std::cout << std::endl;

			std::uniform_int_distribution distr(0, static_cast<int>(goodReversals.size()) - 1);
			const int rdmidx = distr(rng);
			reversal = goodReversals[rdmidx];
			break;
		}
		case ReversalType::NEUTRAL_GOOD: {
			std::cout << " --- Neutral good " << std::endl;
			
			// Case: many hurdles or few hurdles.
			if(hurdles.size() > 1){
				std::vector<int> hurdle_exts_own_diff;
				std::set_difference(hurdle_exts_own.begin(), hurdle_exts_own.end(), cycle_exts.begin(), cycle_exts.end(),
									std::inserter(hurdle_exts_own_diff, hurdle_exts_own_diff.begin()));
				std::vector<std::pair<int,int>> neutralGoodReversals = makeAllCombinations(cycle_exts, hurdle_exts_own_diff);
				
				std::vector<std::pair<int,int>> cycReversals = getReversalsFromCycle(cycle_idx);
				neutralGoodReversals.insert(neutralGoodReversals.end(), cycReversals.begin(), cycReversals.end());

				std::uniform_int_distribution distr(0, static_cast<int>(neutralGoodReversals.size()) - 1);
				const int rdmidx = distr(rng);
				reversal = neutralGoodReversals[rdmidx];

			// Case: Single hurdle (h=1).
			} else {
				std::cout << "ERROR! Cycles in a single hurdle do not have 'neutral good' reversals." << std::endl;
			}
			break;
		}
		case ReversalType::NEUTRAL: {
			std::cout << "ERROR! Cycles in hurdles do not have 'neutral' reversals." << std::endl;
			break;
		}
		case ReversalType::BAD: {
			std::cout << " --- Bad " << std::endl;

			// All gene extremities.
			std::vector<int> gene_exts_all = comps.getGeneExtremities(true);

			// Gene extremities not in hurdles.
			std::vector<int> exts_not_in_hurdle;
			std::set_difference(gene_exts_all.begin(), gene_exts_all.end(), hurdle_exts_all.begin(), hurdle_exts_all.end(),
								std::inserter(exts_not_in_hurdle, exts_not_in_hurdle.begin()));

			std::vector<std::pair<int,int>> badReversals = makeAllCombinations(cycle_exts, exts_not_in_hurdle);
			std::uniform_int_distribution distr(0, static_cast<int>(badReversals.size()) - 1);
			const int rdmidx = distr(rng);
			reversal = badReversals[rdmidx];

			break;
		}
	}
	return reversal;
}

std::pair<int,int> ReversalSampler::sampleReversalFromCycle(const int cycle_idx, const ReversalType revtype, std::mt19937& rng){
	std::pair<int, int> reversal{-1, -1};

	// Depending on the type of cycle and the type of reversal, some categories were pre-computed
	// while counting the number of reversals in a category. As it is prohibitive to do this for
	// all categories, not all of them have a list of reversals pre-computed.
	switch(cycles_info[cycle_idx].type) {

		case CycleType::ORIENTED: {
			std::cout << "It is oriented ";
			reversal = sampleOriented(cycle_idx, revtype, rng);
			break;
		}
		case CycleType::UNORIENTED: {
			std::cout << "It is UNoriented ";
			reversal = sampleUnoriented(cycle_idx, revtype, rng);
			break;
		}
		case CycleType::TRIVIAL: {
			std::cout << "It is trivial ";
			reversal = sampleTrivial(cycle_idx, revtype, rng);
			break;
		}
		case CycleType::HURDLE: {
			std::cout << "It is hurdle ";
			reversal = sampleHurdle(cycle_idx, revtype, rng);
			break;
		}
	}
	return reversal;
}

ReversalRandom ReversalSampler::sampleReversal(std::mt19937& rng, bool updateComps){
	// Update components if needed.
	if(updateComps){
		updateComponents();
		countReversals();
	}
	// Sample reversal type.
	ReversalType revtype = sampleReversalType(rng);
	std::cout << "\n > Sampled type " << revtype << " : nb. rev=" << rev_totals[revtype] << "; prob.=" << rev_weights[revtype] << std::endl;
	
	// Sample a cycle based on how many reversals of the chosen reversal type each cycle has.
	int cycle_idx = sampleCycle(revtype, rng);

	// Sample a reversal given the chosen cycle and the chosen type.
	std::pair<int,int> blackEdges = sampleReversalFromCycle(cycle_idx, revtype, rng);
	ReversalRandom reversal(comps.getGene(blackEdges.first), comps.getGene(blackEdges.second), // start of black edges
							comps.getGene(comps.perm[comps.idxs[blackEdges.first] +1]), // end point of black edge
							comps.getGene(comps.perm[comps.idxs[blackEdges.second]+1]), // end point of black edge
							revtype, cycles_info[cycle_idx].type, rev_totals);

	std::cout << "- Sampled reversal: exts=(" << blackEdges.first << ", " << blackEdges.second << "); genes: (" << reversal.g_beg << ", " << reversal.g_end << "]" << std::endl;
	return reversal;
}

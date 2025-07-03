/*******************************************************
 * This class implements the algorithm described in 
 * Hannenhalli and Pevzner (1999) to transform unoriented 
 * components into oriented components with a minimum
 * number of reversals.
 * 
 * The idea consists of first sorting the unoriented components 
 * by their rightmost element, and then merge non-consecutive 
 * hurdles H_1 and H_2 by applying a reversal that affects a
 * gene in H_1 and a gene in H_2.
 * 
 * Once all components are oriented, they can be sorted
 * using, for example, the algorithm from Tannier et al. (2007). 
 * 
 * For more details on the algorithm to transform unoriented
 * connected components of the overlap graph into oriented ones, 
 * please check:
 * 
 * Hannenhalli and Pevzner, 1999: "Transforming cabbage into 
 * turnip: polynomial algorithm for sorting signed 
 * permutations by reversals" (Lemma 18 and Lemma 19).
 * 
 * Haim Kaplan et al., 1997: "Faster and simpler algorithm for 
 * sorting signed permutations by reversals" (Lemma 5.1);
 * 
 *******************************************************/

#include <iostream>
#include <unordered_map>
#include <vector>
#include <stack>
#include <list>
#include <cstdlib>   // exit
#include <cmath>	 // abs, trunc, ceil
#include <memory>	 // shared_ptr
#include <utility>   // move, pair
#include <algorithm> // sort, swap

#include "findComponents_Bader2001.hpp"
#include "reversal.hpp"

Reversal getReversal(int gen_ext_beg, int gen_ext_end, const ConnectedComponents& comps) {

	std::cout << "1-> gen_ext_beg=" << gen_ext_beg << " [" << comps.idxs[gen_ext_beg] << "] " << " gen_ext_end=" << gen_ext_end << " [" << comps.idxs[gen_ext_end] << "] " << std::endl;

	// Check if start appears before the end.
	if (comps.idxs[gen_ext_beg] > comps.idxs[gen_ext_end]){
		std::swap(gen_ext_beg, gen_ext_end);
		gen_ext_beg = comps.getBlackEdge(gen_ext_beg, true);
		gen_ext_end = comps.getBlackEdge(gen_ext_end, false);
	}

	std::cout << "2-> gen_ext_beg=" << gen_ext_beg << " [" << comps.idxs[gen_ext_beg] << "] " << " gen_ext_end=" << gen_ext_end << " [" << comps.idxs[gen_ext_end] << "] " << std::endl;

	// Get genes from gene adjacencies affected by reversal.
	const int gen_beg = comps.getGene(gen_ext_beg);
	const int gen_end = comps.getGene(gen_ext_end);
	const int gen_beg_next = comps.getGene(comps.perm[comps.idxs[gen_ext_beg]+1]);
	const int gen_end_next = comps.getGene((((comps.idxs[gen_ext_end]+1) < comps.perm.size()) ? comps.perm[comps.idxs[gen_ext_end]+1] : comps.perm.size()));
	return Reversal(gen_beg, gen_end, gen_beg_next, gen_end_next);
}

// Get a **random** black edge, so different parsimonious scenarios can be sampled.
std::vector<Reversal> clearUnorientedComponents(std::vector<int> extendedPerm, ConnectedComponents& comps, std::mt19937& rng) {

	// List of reversals to transform all unoriented components into oriented ones.
	std::vector<Reversal> reversals;

	std::vector<int> sortedUnorientedComps;
	sortedUnorientedComps.reserve(comps.rootList.size());

	// Make a list containing only the unoriented components.
	for (int const& root_idx : comps.rootList) {
		if(!comps.forest[root_idx].oriented){sortedUnorientedComps.emplace_back(root_idx);}
	}

	// Sort unoriented components by the index of their rightmost element.
	auto compareRightmostIdx = [&comps](int a, int b) {return comps.forest[a].max < comps.forest[b].max;};
	std::sort(sortedUnorientedComps.begin(), sortedUnorientedComps.end(), compareRightmostIdx);

	// Print unoriented components sorted by rightmost element.
	std::cout << "\nUnoriented components = " << sortedUnorientedComps.size() << std::endl;
	for (int const& root_idx : sortedUnorientedComps) {
		std::cout << ">Component max=" << comps.forest[root_idx].max << std::endl;
		comps.printComponent(comps.forest[root_idx], "", comps.perm.size(), comps.forest);
	}

	// Apply minimum number of reversals to transform unoriented 
	// components into oriented components.
	int comp_2 = (int)std::trunc(sortedUnorientedComps.size()/2);
	for (int comp_1=0; comp_1 < comp_2; ++comp_1) {

		// Get one black edge in each component.
		int gen_ext_beg = comps.getRandomBlackEdge(comps.forest[sortedUnorientedComps[comp_1]], rng, true);
		int gen_ext_end = comps.getRandomBlackEdge(comps.forest[sortedUnorientedComps[comp_1+comp_2]], rng, false);

		Reversal rev = getReversal(gen_ext_beg, gen_ext_end, comps);
		reversals.emplace_back(rev);

		// Merge two non-consecutive unoriented components with one reversal.
		std::cout << " [Merge components] Reversal: (" << rev.g_beg << " " << rev.g_end << "]" << std::endl;
	}

	// Check if there is an odd number of components.
	// If so, "cut" the last component, by applying one 
	// reversal that affects two black edges of the last component,
	// to make it oriented (see Proposition 10.18 in the book
	// "Mathematics of Evolution and Phylogeny" from Olivier Gascuel).
	if (sortedUnorientedComps.size() % 2 == 1) {

		const int idx_last = sortedUnorientedComps.size()-1;

		// Get two black edges from the same component.
		// This is done in the case where the last component in the list  
		// remains untouched by the previous merging operations. 
		// This happens when the list of unoriented components has an odd number of elements.
		int gen_ext_beg = comps.getRandomBlackEdge(comps.forest[sortedUnorientedComps[idx_last]], rng,  true);
		int gen_ext_end = gen_ext_beg;
		while(gen_ext_end == gen_ext_beg){
			gen_ext_end = comps.getRandomBlackEdge(comps.forest[sortedUnorientedComps[idx_last]], rng,  true);
		}
		gen_ext_end  = comps.getBlackEdge(gen_ext_end, false);

		Reversal rev = getReversal(gen_ext_beg, gen_ext_end, comps);
		reversals.emplace_back(rev);

		// Merge two non-consecutive unoriented components with one reversal.
		std::cout << " [Split last component] Reversal: (" << rev.g_beg << " " << rev.g_end << "]" << std::endl;
	}
	
	return reversals;
}


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
#include <vector>
#include <list>
#include <cstdlib>   // exit
#include <cmath>	 // abs, trunc, ceil
#include <memory>	 // shared_ptr
#include <utility>   // move, pair
#include <algorithm> // sort, swap

#include "solveUnoriented_HannenhalliPevzner1999.hpp"

/* It returns a sequence of elements that appear in unoriented 
components, sorted by their indices in the permutation. For more details,
see the definition of ``CR`` in Kaplan, Shamir, and Tarjan (1997). */
std::vector<std::pair<int,int>> UnorientedComponents::getSortedUnorientedElements() {
	// Make a list containing only elements from unoriented components.
	std::vector<std::pair<int,int>> sortedUnorientedElems;
	for (int const& root_idx : comps.rootList) {
		if((!comps.forest[root_idx].oriented) && (comps.forest[root_idx].genes.size() > 2)){ // 2 genes = trivial component (already sorted).
			for (int const& gene : comps.forest[root_idx].genes) {
				sortedUnorientedElems.emplace_back(std::make_pair(gene,root_idx));
			}
		}
	}
	// Sort elements from unoriented components by their index.
	auto compareIdx = [&](std::pair<int,int> a, std::pair<int,int> b) {return comps.idxs[a.first] < comps.idxs[b.first];};
	std::sort(sortedUnorientedElems.begin(), sortedUnorientedElems.end(), compareIdx);
	return sortedUnorientedElems;
}

Reversal UnorientedComponents::getReversal(int gen_ext_beg, int gen_ext_end) {
	// Check if start appears before the end.
	if (comps.idxs[gen_ext_beg] > comps.idxs[gen_ext_end]){
		std::swap(gen_ext_beg, gen_ext_end);
	}
	// Get genes from gene adjacencies affected by reversal.
	const int gen_beg = comps.getGene(gen_ext_beg);
	const int gen_end = comps.getGene(gen_ext_end);
	const int gen_beg_next = comps.getGene(comps.perm[comps.idxs[gen_ext_beg]+1]);
	const int gen_end_next = comps.getGene((((comps.idxs[gen_ext_end]+1) < comps.perm.size()) ? comps.perm[comps.idxs[gen_ext_end]+1] : comps.perm.size()));
	return Reversal(gen_beg, gen_end, gen_beg_next, gen_end_next);
}

/* A hurdle is an unoriented component whose elements 
appear contiguous in a sequence considering only elements
from unoriented components. 

They are important because they cannot become oriented by 
a reversal that affects breakpoints from other 
unoriented components. 

Notice that unoriented components that are not a hurdle 
could become oriented by reversals applied to other unoriented 
components.
*/
std::vector<int> UnorientedComponents::findHurdles() {
	// Sort elements from unoriented components by their index.
	// The sorted list has as elements a pair of ints, where the first 
	// int refers to the gene value, whereas the second int refers to 
	// the component where the gene is (more precisely, it is the index of the tree's root).
	std::vector<std::pair<int,int>> sortedUnorientedElems = getSortedUnorientedElements();
	// Print sorted unoriented elements .
	if(debug){
		std::cout << std::endl << "Unoriented elements = " << sortedUnorientedElems.size() << std::endl;
		for (std::pair<int, int> const& elem : sortedUnorientedElems) {
			std::cout << elem.first << "[" << elem.second << "] ";
		}
		std::cout << std::endl;
	}
	// Find hurdles.
	std::vector<int> hurdles;
	if (sortedUnorientedElems.size() == 0) {
		if(debug){std::cout << std::endl << "Hurdles = 0 (no unoriented components)" << std::endl;}
		return hurdles;
	} 
	int current_idx =  0;
	int count_last  = -1;
	const int root_idx_last = sortedUnorientedElems[sortedUnorientedElems.size()-1].second;
	while(current_idx < sortedUnorientedElems.size()){
		int root_idx = sortedUnorientedElems[current_idx].second;

		// Check if the first and last elements are in the same component.
		// This is done to determine if the unoriented component with the rightmost element
		// is a hurdle, as this component may wrap around to the beginning of the list.
		if(count_last < 0){
			++count_last; 
			while(root_idx == root_idx_last){
				++count_last; ++current_idx;
				root_idx = ((current_idx < sortedUnorientedElems.size()) ? sortedUnorientedElems[current_idx].second : -1);
			}
			if (root_idx < 0){break;}
		}

		// Except by the ``last component wrapping around at the beginning`` case,
		// all other components are checked directly: if the unoriented component 
		// is a hurdle, then the rightmost element will be located at:
		// [leftmost element's index] + [number of elements in the component] - 1.
		const int last_elem_expected = comps.perm[comps.forest[root_idx].max];
		const int last_elem_idx      = current_idx+comps.forest[root_idx].genes.size()-1;
		const int last_elem_found    = ((last_elem_idx < sortedUnorientedElems.size()) ? sortedUnorientedElems[last_elem_idx].first : -1);

		// Current unoriented component is a hurdle.
		if(last_elem_expected == last_elem_found){
			hurdles.emplace_back(root_idx);
			current_idx = last_elem_idx+1; // Advance to the next component.

		// Current unoriented component is NOT a hurdle.
		} else {
			// Advance to the next component.
			const int root_idx_current = root_idx;
			while(root_idx == root_idx_current){
				++current_idx;
				root_idx = ((current_idx < sortedUnorientedElems.size()) ? sortedUnorientedElems[current_idx].second : -1);
			}
		}
	}
	// Check if the last unoriented component is a hurdle.
	if(count_last > 0){
		if(count_last < comps.forest[root_idx_last].genes.size()){
			current_idx  = sortedUnorientedElems.size()-1;
			int root_idx = sortedUnorientedElems[current_idx].second;
			while(root_idx == root_idx_last) {
				--current_idx; ++count_last;
				root_idx = (current_idx < 0) ? -1 : sortedUnorientedElems[current_idx].second;
			}			
		}
		if(count_last==comps.forest[root_idx_last].genes.size()){
			hurdles.emplace_back(root_idx_last);
		} 
	}
	// Print hurdles sorted by rightmost element.
	if(debug){
		std::cout << std::endl << "Hurdles = " << hurdles.size() << std::endl;
		for (int const& root_idx : hurdles) {
			std::cout << ">Component max=" << comps.forest[root_idx].max << std::endl;
			comps.printComponent(comps.forest[root_idx], "", comps.perm.size(), comps.forest);
		}
	}
	return hurdles;
}

// Get a **random** black edge, so different parsimonious scenarios can be sampled.
std::vector<Reversal> UnorientedComponents::clearUnorientedComponents(std::mt19937& rng) {

	if(debug) genperm.printBlocks();

	// List of reversals to transform *all* unoriented components into oriented ones.
	std::vector<Reversal> reversals;
	std::vector<int> hurdles = findHurdles();
	nb_hurdles = hurdles.size();
	if(nb_hurdles == 0){return reversals;}
	
	// Apply minimum number of reversals to transform unoriented 
	// components into oriented components.
	int comp_2 = (int)std::trunc(hurdles.size()/2.0);
	for (int comp_1=0; comp_1 < comp_2; ++comp_1) {

		// Get one black edge in each component.
		int gen_ext_beg = comps.getRandomBlackEdge(comps.forest[hurdles[comp_1]], rng, true);
		int gen_ext_end = comps.getRandomBlackEdge(comps.forest[hurdles[comp_1+comp_2]], rng, true);

		// Merge two non-consecutive unoriented components with one reversal.
		Reversal rev = getReversal(gen_ext_beg, gen_ext_end);
		if(debug) std::cout << "[Merge components] Reversal: (" << rev.g_beg << " " << rev.g_end << "]" << std::endl;
		applyReversal(genperm, rev.g_beg, rev.g_end);
		if(debug) genperm.printBlocks();
		// Labels correspond to extended signed permutation.
		reversals.emplace_back(rev); 
	}

	// Check if there is an odd number of components.
	// If so, "cut" the last component, by applying one 
	// reversal that affects two black edges of the last component,
	// to make it oriented (see Proposition 10.18 in the book
	// "Mathematics of Evolution and Phylogeny" from Olivier Gascuel).
	if (hurdles.size() % 2 == 1) {

		const int idx_last = hurdles.size()-1;

		// Get two black edges from the same component.
		// This is done in the case where the last component in the list  
		// remains untouched by the previous merging operations. 
		// This happens when the list of unoriented components has an odd number of elements.
		int gen_ext_beg = comps.getRandomBlackEdge(comps.forest[hurdles[idx_last]], rng,  true);
		int gen_ext_end = gen_ext_beg;
		while(gen_ext_end == gen_ext_beg){
			gen_ext_end = comps.getRandomBlackEdge(comps.forest[hurdles[idx_last]], rng,  true);
		}
		// Merge two non-consecutive unoriented components with one reversal.
		Reversal rev = getReversal(gen_ext_beg, gen_ext_end);
		if(debug) std::cout << "[Split last component] Reversal: (" << rev.g_beg << " " << rev.g_end << "]" << std::endl;
		applyReversal(genperm, rev.g_beg, rev.g_end);
		genperm.printBlocks();
		// Labels correspond to extended signed permutation.
		reversals.emplace_back(rev);
	}
	return reversals;
}

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
#include <cmath>	 // abs
#include <memory>	 // shared_ptr
#include <utility>   // move, pair
#include <algorithm> // sort

#include "findComponents_Bader2001.hpp"

void clearUnorientedComponents(ConnectedComponents& comps) {
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
	std::cout << "\n Unoriented components = " << sortedUnorientedComps.size() << std::endl;
	for (int const& root_idx : sortedUnorientedComps) {
		std::cout << "> Component max=" << comps.forest[root_idx].max << std::endl;
		comps.printComponent(comps.forest[root_idx], "", comps.perm.size(), comps.forest);
	}

	
}

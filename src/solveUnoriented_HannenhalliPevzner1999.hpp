/*******************************************************
 * Original author: Priscila Biller
 * Created: July/2025
 * License: GPL v3
 * 
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

#pragma once // It avoids class redefinition.

#include <iostream>
#include <vector>
#include <list>
#include <cstdlib>   // exit
#include <cmath>	 // abs, trunc, ceil
#include <memory>	 // shared_ptr
#include <utility>   // move, pair
#include <algorithm> // sort, swap

#include "findComponents_Bader2001.hpp"
#include "reversal.hpp"
#include "utils.hpp"

class BlockSimple : public BlockBase<BlockSimple> {
public:
    // Inherit constructors from BlockBase
    using BlockBase<BlockSimple>::BlockBase;
};

class UnorientedComponents {
protected:
	
	/* It returns a sequence of elements that appear in unoriented 
	components, sorted by their indices in the permutation. For more details,
	see the definition of ``CR`` in Kaplan, Shamir, and Tarjan (1997). */
	std::vector<std::pair<int,int>> getSortedUnorientedElements();

	Reversal getReversal(int gen_ext_beg, int gen_ext_end);

	int getGeneExtIndex(const int gene_ext) const;
	int getRandomBlackEdge(const Cycle& comp, std::mt19937& rng, const bool beg_ext) const;

public:
	ConnectedComponents& comps;
	GenomePermutation<BlockSimple>& genperm;
	int debug{DEBUG_OFF};

	int nb_hurdles{0};

	UnorientedComponents(GenomePermutation<BlockSimple>& genperm, ConnectedComponents& comps, int debug=DEBUG_OFF):genperm(genperm),comps(comps),debug(debug){

	}

	/* Get a **random** black edge, so different parsimonious scenarios can be sampled.*/
	std::vector<Reversal> clearUnorientedComponents(std::mt19937& rng);

	/* A hurdle is an unoriented component whose elements 
	appear contiguous in a sequence considering only elements
	from unoriented components. 

	They are important because they cannot become oriented by 
	a reversal that affects breakpoints from other 
	unoriented components. 

	Notice that unoriented components that are not a hurdle 
	could become oriented by reversals applied to other unoriented 
	components. */
	std::vector<int> findHurdles();

};

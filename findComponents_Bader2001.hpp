/*******************************************************
 * This class implements the data structure described in 
 * Bader, Moret, and Yan (2001), called ``overlap forest``.
 * Given two genomes, the overlap forest is a way to 
 * find the connected components of the overlap graph 
 * (i.e. cycles in the cycle graph that partially overlap) 
 * in linear time.
 * 
 * Finding the connected components of the overlap graph
 * in an efficient way is an important step in solving
 * several genome rearrangements problems including,
 * for example, the sorting by reversals problem.
 * 
 * Each of these connected components (if they are oriented)
 * can be handled independently, facilitating the 
 * computation of a solution.
 * 
 * For more details on the data structures and the algorithm 
 * to find connected components of the overlap graph, 
 * please check:
 * 
 * David Bader et al., 2001: "A Linear-Time Algorithm for 
 * Computing Inversion Distance between Signed Permutations 
 * with an Experimental Study";
 * 
 * For extra references with slightly modified approaches:
 * 
 * Garg et al., 2019: "Sorting by Reversals: A faster Approach 
 * for Building Overlap Forest";
 * 
 * Olivier Gascuel, 2005: "Mathematics of evolution and phylogeny" 
 * (book's chapter "The inversion distance problem", section 10.5.
 * 
 *******************************************************/

#pragma once // It avoids class redefinition.

#include <iostream>
#include <unordered_map>
#include <vector>
#include <stack>
#include <list>
#include <cstdlib>   // exit
#include <cmath>	 // abs
#include <memory>	 // shared_ptr
#include <utility>   // move, pair
#include <random>

/*******************************************************
 *  Auxiliary structures
*******************************************************/

// Color of edges in the Breakpoint Graph (BG).
enum class BgColor {
	GRAY,   // GRAY is assigned the value 0
	BLACK   // BLACK is assigned the value 1
};

/* This class store information about a cycle in a breakpoint graph.*/
class Cycle {
public:
	int parent{-1};
	int id;
	std::list<int>::iterator root; // It keeps a valid iterator if this cycle is a root in the overlap forest.
	std::vector<int> children;     // It can store genes (indices between 0 and 2n+1) or cycles (index of the cycle + [2n+1]).
	std::vector<int> genes;        // It store all genes (indices between 0 and 2n+1).
	int min; // Smallest index of an element belonging to the forest rooted by this cycle.
	int max; // Biggest index of an element belonging to the forest rooted by this cycle.
	bool oriented{false}; // true if the cycle contains two genes with flipped signs; false otherwise.
	Cycle(const int id, const int min_idx, const int max_idx):id(id),min(min_idx),max(max_idx){
	}
	std::string printCycle() const {
		//std::vector<int> children
		std::string children_str = "";
		for(const int& c: children) {children_str += (std::to_string(c) + " ");}
		return "C_" + std::to_string(id) + ": parent=" + std::to_string(parent) + ";min=" + std::to_string(min) + ";max=" + std::to_string(max) + ";children=" + children_str;
	}
};

/*******************************************************
 *  Connected component data structure.
*******************************************************/

/* Implements algorithm from Bader et al. (2001):
"A Linear-Time Algorithm for Computing Inversion Distance between Signed Permutations with an Experimental Study".
*/
class ConnectedComponents {
protected:
	void findConnectedComponents(); // TODO: Implement a version for multichromosomal genomes.
	
public:
	
	std::vector<int> perm;   // pos i stores a gene extremity.
	std::vector<int> idxs;   // pos i stores the current position of gene extremity i in the permutation.
	std::vector<int> cycles; // pos i stores the cycle id that gene extremity i belongs.

	std::vector<Cycle> forest; // Store all cycles in the breakpoint graph, one element per cycle.
	std::list<int> rootList;   // Store all roots in the overlap (forest) graph, one element per root.

	ConnectedComponents(const std::vector<int> unsignedExtPerm):perm(unsignedExtPerm),idxs(perm.size()),cycles(perm.size(),-1){
		printUnsignedExtPerm();
		findConnectedComponents();
		printComponents();
	}

	void printComponent(const Cycle& comp, std::string indent, const int& n, const std::vector<Cycle>& forest) const;

	void printUnsignedExtPerm() const;
	void printComponents() const;

	int getBlackEdge(const int gene_ext, const bool beg_ext) const;
	int getRandomBlackEdge(const Cycle& comp, std::mt19937& rng, const bool beg_ext) const;

	inline int getGene(const int gene_ext) const {return (int)(std::ceil(gene_ext/2.0)+1);}
};

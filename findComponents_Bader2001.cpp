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
#include <algorithm> // shuffle

#include "findComponents_Bader2001.hpp"

/*******************************************************
 *  Connected component data structure.
*******************************************************/

void ConnectedComponents::printComponent(const Cycle& comp, std::string indent, const int& n, const std::vector<Cycle>& forest) const {
	const int comp_idx = comp.id-n;
	std::cout << indent << comp_idx << " = ";
	for(const int& c: comp.children) {if(c < n){std::cout << c << " ";}}
	std::cout << " [min=" << comp.min << ";max=" << comp.max << ";oriented=" << (comp.oriented ? "T" : "F") << "]" << std::endl;
	indent += "   ";
	for(const int& c: comp.children) {if(c >= n){printComponent(forest[c-n], indent, n, forest);}}
}
	
void ConnectedComponents::findConnectedComponents() {

	std::stack<int> stack; // Store all cycles currently "active", one element per cycle.

	for (int i=0; i<perm.size(); ++i){idxs[perm[i]] = i;}
	// Find components.
	int cycle_id = perm.size();
	for (int i=0; i<perm.size(); i+=2){ // Loops every 2 elements.
		
		// std::cout << "- Element " << i << ": " << perm[i] << std::endl;

		// Case 1: element i has no cycle. 
		if(cycles[i] < 0) {
			// The id of a cycle corresponds to the index of the cycle in the forest, plus perm.size().
			const int cycle_idx = cycle_id-perm.size();

			// Each node in the overlap forest is a cycle.
			// The list of roots saves the indices of all 
			// cycles (in the vector ``forest``) that 
			// are roots in the forest (i.e. they have no parent).
			rootList.emplace_back(cycle_idx);

			// Create a new cycle containing i.
			forest.emplace_back(cycle_id,i,i);
			Cycle& new_cycle = forest[cycle_idx];
			new_cycle.root   = std::prev(rootList.end());

			BgColor edge_color = BgColor::GRAY;
			int current = perm[i];
			// Finds the cycle by traversing the alternating gray and black edges from the cycle.
			// std::cout << "\tNew cycle : ";
			do {
				// std::cout << " " << current;
				// Map element to the new cycle.
				cycles[idxs[current]] = cycle_idx;
				new_cycle.children.emplace_back(current);

				// Update maximum index if needed.
				if(new_cycle.max < idxs[current]) {new_cycle.max = idxs[current];}

				// Case 1.1 : Traversing a gray edge.
				if(edge_color == BgColor::GRAY){
					// A gray edge connecting genes i and i+1 always starts at 
					// an even value (2*i) and points to an odd value (2*(i+1)-1).
					// We use the parity to infer at which side of the edge we are now,
					// and where we should go (if we are at the head, we go to the tail, and vice-versa).
					const bool idx_cur_parity  = (idxs[current] % 2);
					current += (((current % 2)==0) ? 1 : -1);
					const bool idx_next_parity = (idxs[current] % 2);

					edge_color = BgColor::BLACK;
					// Check if gray edge is oriented (i.e. if connects two genes with flipped signs).
					// An edge is oriented if the index of both of its extremities have the same parity.
					if(idx_cur_parity == idx_next_parity){new_cycle.oriented = true;}
				// Case 1.2 : Traversing a black edge.
				} else {
					// A black edge points to the consecutive element in the current permutation.
					// Depending on the orientation of the gene (which can be inferred by the parity
					// of the index), we move left or right in the current permutation.
					current = (((idxs[current] % 2)==0) ? perm[idxs[current]+1] : perm[idxs[current]-1]);
					edge_color = BgColor::GRAY;
				}
				// std::cout << "(next:" << current << ")";
			} while ((idxs[current]-i) != 0);
			// std::cout << std::endl;
			// Add cycle to the stack.
			// Do not add only if the max index comes just next the current index.
			if(new_cycle.min+1 < new_cycle.max){
				stack.push(cycle_idx);	
			}
			++cycle_id;

		// Case 2: Cycle of element i was already created.
		} else {
			// std::cout << "\tExistent cycle with element "  << perm[i] << ": ";

			// Find the cycle that is the root of the 
			// connected component containing element i.
			int root_idx = cycles[i]; // Index of cycle i in the list ``forest``.
			while(forest[root_idx].parent > 0){root_idx = forest[root_idx].parent;}
			// std::cout << " " << forest[root_idx].printCycle() << std::endl;

			// Merge cycles that overlap.
			while(forest[stack.top()].min > forest[root_idx].min){
				// Update max if needed.
				if(forest[root_idx].max < forest[stack.top()].max){
					forest[root_idx].max = forest[stack.top()].max;
				}
				// Update parity if needed.
				if(!forest[root_idx].oriented){forest[root_idx].oriented = forest[stack.top()].oriented;}
				// Update roots of the forest.
				forest[stack.top()].parent = root_idx;
				forest[root_idx].children.emplace_back(forest[stack.top()].id);
				// Top of the stack is not a root anymore. Delete root.
				rootList.erase(forest[stack.top()].root);
				// Remove top of the stack.
				stack.pop();
			}
			if(i == (forest[root_idx].max-1)){stack.pop();}
		}
	}
}

void ConnectedComponents::printUnsignedExtPerm() const {
	std::cout << "Unsigned extended permutation:\n";
	for (int i=0; i<perm.size(); ++i){std::cout << perm[i] << " ";}
	std::cout << std::endl;
}

void ConnectedComponents::printComponents() const {
	// Print components.
	std::cout << "Components:" << std::endl;
	for (int const& root_idx : rootList) {
		printComponent(forest[root_idx], "", perm.size(), forest);
	}
	std::cout << std::endl;
}

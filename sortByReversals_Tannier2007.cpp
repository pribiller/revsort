/*******************************************************
 * One genome is represented as one permutation. 
 * The permutation is given as a vector of ints,
 * where each int corresponds to a gene, and its sign
 * indicates the direction its transcription occurs.
 * 
 * This class stores this permutation as blocks
 * of size Θ(√n×log(n)), with the additional property
 * that each block has a balanced binary tree associated
 * with, where safe reversals can be easily retrieved.
 * 
 * For more details on the data structures and 
 * the algorithm to sort by reversals, please check:
 * 
 * Eric Tannier et al., 2007: "Advances on sorting by reversals";
 * 
 *******************************************************/

#include <iostream>
#include <list>
#include <vector>
#include <numeric> // iota
#include <cmath> 
#include <string> // convert input args (string to int)
#include <utility> // move
#include <algorithm> // find_if, reverse
#include <unordered_map>
#include <unordered_set>
#include "genomePermutation.hpp"
//#include <memory> // shared_ptr
//#include <chrono>

#include <boost/intrusive/rbtree.hpp>

// It makes code a bit clear: 
// instead of std::cout you can type directly cout.
//using namespace std;

/*******************************************************
 * Data structures to implement balanced 
 * binary trees for blocks.
 *******************************************************/

/* Each node corresponds to an arc v_i. */
                    // Base hook with default options (i.e., optimized for speed).
template <typename BlockT>
class Node : public boost::intrusive::set_base_hook<> {
public:
	const Gene<BlockT>& gene;
	const Gene<BlockT>& gene_next;
	int orientedTotal{0}; // the number of oriented arcs in the subtree rooted at node i.
	bool reversed{false}; // nodes should be ordered backwards with respect to the original order, and all nodes change orientation.

	Node(const Gene<BlockT>& g, const Gene<BlockT>& g_next):gene(g),gene_next(g_next){ 
	}

	// Returns if the arc i is oriented or not (oriented: sign of i and (i+1) are flipped).
	inline bool getOrientation() const {
		return (gene.reversed != gene_next.reversed);
	}

	// Returns the position of node i's successor in the current permutation.
	inline int getPosNext() const {
		return gene_next.block->genePosAbs(gene_next);
	}

	// Overload the '<' operator.
	friend bool operator<(const Node &a, const Node &b){
		return a.getPosNext() < b.getPosNext();
	}
	// Overload the '>' operator.
	friend bool operator>(const Node &a, const Node &b){
		return a.getPosNext() > b.getPosNext();
	}
	// Overload the '==' operator.
	friend bool operator==(const Node &a, const Node &b){
		return a.getPosNext() == b.getPosNext();
	}
};

class BlockTree : public BlockBase<BlockTree> {
public:
	boost::intrusive::rbtree<Node<BlockTree>> tree; // Balanced binary tree with segment elements as nodes.	

    // Inherit constructors from BlockBase
    using BlockBase<BlockTree>::BlockBase;
    
	void makeTree(std::vector<Node<BlockTree>>& nodes){
		tree.clear();
		for (std::list<Gene<BlockTree>>::iterator gene = permutationSegment.begin(); gene != permutationSegment.end(); ++gene){
			if(gene->id < nodes.size()){ // The last gene does not have an arc.
				tree.insert_equal(nodes[gene->id-1]);				
			}
		}
	}
	void updateTree(){}
};

/*******************************************************
 * Class to find a sorting scenario between two genomes.
 *******************************************************/

class GenomeSort {
private:

	// Permutation stored as blocks, and each block 
	// has an additional structure to quickly access
	// safe reversals (balanced binary tree with arcs as nodes).
	GenomePermutation<BlockTree> genperm;

	// References of nodes in balanced binary trees for each block (size=n-1 arcs).
	std::vector<Node<BlockTree>> nodes;

	/* Method used only during construction of the object, 
	to initialize tree of arcs for each block. */
	void initializeTrees(){
		for(BlockTree &b : genperm.blockList) {b.makeTree(nodes);}
	}
	
	/* Method used only during construction of the object, 
	to initialize references to nodes. */
	void initializeNodes(){
		nodes.reserve(genperm.n-1);
		for (int g=0; g<(genperm.n-1); g++) {nodes.emplace_back((*genperm.genes[g]), (*genperm.genes[g+1]));}
	}
	
public:
	// Parameterized constructor.
	GenomeSort(const std::vector<int>& perm): genperm(perm){
		// Initialize map of arcs.
		initializeNodes();
		// Initialize trees.
		initializeTrees();
		// Print blocks.
		//std::cout << genperm.printBlocks() << std::endl;
	}

	void applyReversal(int g_beg, int g_end){

		g_beg = std::abs(g_beg);
		g_end = std::abs(g_end);

		std::string applyReversal_str = "<applyReversal (" + std::to_string(g_beg) + ", " + std::to_string(g_end) + "] >";
		std::cout << applyReversal_str << " Before splitting blocks: " << genperm.printBlocks("\n\t") << std::endl;

		// (1) Split at most two blocks so that the endpoints of the reversal correspond to endpoints of blocks;
		genperm.splitBlock(g_beg);
		genperm.splitBlock(g_end);
		// (1.1) After splitting blocks, we must reconstruct from scratch the associated trees.

		std::cout << applyReversal_str << " After splitting blocks: " << genperm.printBlocks("\n\t") << std::endl;

		// (2) Flip ``reversed`` flag of each block between the endpoints of the reversal;
		std::list<BlockTree>::iterator reversal_beg  = std::next(genperm.getBlock(g_beg)); // this block is the first reversed block.
		std::list<BlockTree>::iterator reversal_last = genperm.getBlock(g_end);  // this block is the last reversed block.
		std::list<BlockTree>::iterator reversal_end  = std::next(reversal_last); // this block will not be reversed.
		for (std::list<BlockTree>::iterator b = reversal_beg; b != reversal_end; ++b) { b->reversed = !(b->reversed);}

		// (3) Reverse the order of the blocks between the endpoints of the reversal;
		const int g_after_beg = reversal_beg->permutationSegment.front().id;
		const int g_after_end = reversal_end->permutationSegment.front().id;
		std::reverse(reversal_beg, reversal_end); // Reverses the order of the elements in the range [first, last).
		// (3.1) Update positions of reversed blocks.
		int blockPos = genperm.getBlock(g_beg)->pos + genperm.getBlock(g_beg)->permutationSegment.size();
		for (std::list<BlockTree>::iterator b = genperm.getBlock(g_end); b != std::next(genperm.getBlock(g_after_beg)); ++b) { 
			b->pos = blockPos;
			blockPos += b->permutationSegment.size();
		}

		std::cout << applyReversal_str << " After reversing blocks: " << genperm.printBlocks("\n\t") << std::endl;

		// (4) Concatenate and split blocks in such a way that the size of each block 
		// lies within the interval [½×√(n×log(n)), 2×√(n×log(n))].
		// It must check all 4 blocks that were split.
		genperm.balanceBlock(g_beg);
		genperm.balanceBlock(g_end);
		genperm.balanceBlock(g_after_beg);
		genperm.balanceBlock(g_after_end);
	}

	/* Temporary function to handle memory of boost trees. */
	void clearTrees(){
		for(auto &b : genperm.blockList) {b.tree.clear();}
	}
};

int main(int argc, char* argv[]) {

	int const n = std::stoi(argv[1]); // input (command line argument): number of genes

	// Start vector with Identity permutation (for testing).
	std::vector<int> perm(n);
	std::iota(perm.begin(), perm.end(), 1);

	std::cout << "Hello World! :-D\nVector with " << n << " elements!" << std::endl;
	
	for(int i=0; i<n; i++){
		std::cout << "\t pos. " << i << " = " << perm[i] << std::endl;
	}
	
	GenomeSort genomeSort = GenomeSort(perm);
	// genomeSort.applyReversal(2, 9);   // after rev: 1 2 -9 -8 -7 -6   -5 -4 -3 10 11 .. 20
	// genomeSort.applyReversal(-6, -5); // after rev: 1 2 -9 -8 -7 -6    5 -4 -3 10 11 .. 20
	// genomeSort.applyReversal(-6, 10); // after rev: 1 2 -9 -8 -7 -6  -10  3  4 -5 11 .. 20

	genomeSort.clearTrees();
	std::cout << "Bye bye\n";
	return 0;
}

/*******************************************************
 * One genome is represented as one permutation. 
 * The permutation is given as a vector of integers (int),
 * where each int corresponds to a gene, and its sign
 * indicates the direction its transcription occurs.
 * 
 * This class stores this permutation as blocks
 * of size Θ(√(n×log(n))), with the additional property
 * that each block has a balanced binary tree associated
 * with, where safe reversals can be easily retrieved.
 * 
 * For more details on the data structures and 
 * the algorithm to sort by reversals, please check:
 * 
 * Eric Tannier et al., 2007: "Advances on sorting by reversals";
 * 
 *******************************************************/

#pragma once // It avoids class redefinition.

#include <iostream>
#include <list>
#include <vector>
#include <numeric>   // iota
#include <cmath> 
#include <string>    // convert input args (string to int)
#include <utility>   // move, pair
#include <algorithm> // find_if, reverse
#include <cstdlib>   // exit
#include <deque>
#include <iterator>  // reverse_iterator

#include "genomePermutation_KaplanVerbin2003.hpp"
#include "reversalBST_KaplanVerbin2003.hpp"
#include "reversal.hpp"

/*******************************************************
 * Data structures to implement balanced 
 * binary trees for blocks.
 *******************************************************/

class BlockTree : public BlockBase<BlockTree> {
private:
	// Auxiliary structures used during tree update phase during a reversal.

	// ``t2`` contains all arcs (i, i+1) whose ``remote`` element (i+1) is part of the reversal.
	RedBlackTree<BlockTree> t2;
	// ``t3`` contains all arcs (i, i+1) whose ``remote`` element (i+1) appears after reversal.
	RedBlackTree<BlockTree> t3;

	void printAuxTrees(const int g_beg, const int g_end) const {
		std::cout << "T1 (x <= " << g_beg << ")" << std::endl;
		tree.printTree();
		std::cout << "T2 (" << g_beg << " < x <=" << g_end << ")" << std::endl;
		t2.printTree();
		std::cout << "T3 (x >" << g_end << ")" << std::endl;
		t3.printTree();
	}

public:
	RedBlackTree<BlockTree> tree; // Balanced binary tree with segment elements as nodes.	

    // Inherit constructors from BlockBase
    using BlockBase<BlockTree>::BlockBase;

	void makeTree(std::vector<Node<BlockTree>>& nodes){
		// std::cout << "Make tree from scratch: " << printBlock() << "\n";
		tree.cleanTree();
		for (std::list<Gene<BlockTree>>::iterator gene = permutationSegment.begin(); gene != permutationSegment.end(); ++gene){
			// The last gene does not have an arc.
			// Node list: from 0 until n-2; Gene ids: from 1 until n (last gene (id=n) is ignored).
			if(gene->id <= nodes.size()){tree.insert(&nodes[gene->id-1]);}
		}
		// tree.printTree();
		status = OK;
	}

	// Save state of the tree just before the reversal (g_beg, g_end]
	// After moving blocks around, the tree is going to be in an
	// inconsistent state, i.e., the order of the nodes might not be 
	// correct anymore, so this step needs to be done before the reversal.
	// Split tree into three trees:
	// The middle tree (t2) will contain all pairs whose ``remote`` element 
	// is part of the reversal, and the other two trees (tree and t3) will 
	// contain the pairs whose ``remote`` element is out of the reversal.
	void splitTree(const int g_beg, const int g_end) {
		t2.cleanTree();
		t3.cleanTree();

		// std::cout << "[Before reversal] Split tree before applying reversal (" << g_beg << "," << g_end << "]:" << printBlock() << std::endl;
		// tree.printTree();
		tree.split(g_beg, t2); // at this point: ``tree`` keeps elements before reversal (x_next <= g_beg); ``t2`` has elements after the start of the reversal (x_next > g_beg).

		// std::cout << "[Before reversal]  (" << g_beg << "," << g_end << "]: T1, (T2 & T3):" << printBlock() << std::endl;
		// tree.printTree();
		// t2.printTree();

		t2.split(g_end, t3);   // at this point: ``t2`` has elements inside the reversal (g_beg < x_next <= g_end); ``t3`` has elements after the end of the reversal (x_next > g_end).

		// std::cout << "[Before reversal] Split tree - Resulting trees:" << std::endl;
		// printAuxTrees(g_beg, g_end);
	}

	// Update tree after reversal (g_beg, g_end].
	void mergeTrees(const int g_beg, const int g_end){

		// If T is a tree corresponding to a block outside the reversal, 
		// flip the ﬂag at the root of T2 , and concatenate T1 , T2 and T3.
		// This flip will affect the direction *and* orientation of T2.

		// If T is a tree corresponding to a block inside the reversal,
		// this flip will affect the direction of T2, but the orientation 
		// of T1 and T3 (the ``reversed`` flag of the block is taken into
		// account to compute the orientation).
		if(t2.root != nullptr) {t2.root->reversed = !t2.root->reversed;}

		// std::cout << "[Merge trees] Before join t1 + t2 | Block: " << printBlock() << std::endl;
		// printAuxTrees(g_beg,g_end);

		// Concatenate t1, t2, t3.
		// Concatenate - Step 1: t1 = t1 + t2
		if (tree.root == nullptr) {
			tree.root = t2.root;
		} else {
			// TODO: Depending on tree sizes, get max of one or the min from the other.
			Node<BlockTree>* max_t1 = tree.getGlobalMax();
			// std::cout << "[Merge trees] Max from T1: (" << max_t1->printNode() << ")" << std::endl;
			tree.remove(max_t1);
			// std::cout << "[Merge trees] T1 (x <= " << g_beg << ") after removing (" << max_t1->printNode() << ")" << std::endl;
			// tree.printTree();
			tree.join(max_t1,t2);
		}

		// std::cout << "[Merge trees] After join t1 + t2 | Block: " << printBlock() << std::endl;
		// tree.printTree();

		// Concatenate - Step 2: t1 = t1 + t2 + t3
		if (tree.root == nullptr) {
			tree.root = t3.root;
		} else {
			// TODO: Depending on tree sizes, get max of one or the min from the other.
			// std::cout << "[Merge trees] join t1 + t2 + t3 | Block: " << printBlock() << std::endl;
			// std::cout << "T1 (x <= " << g_beg << ")" << std::endl;
			// tree.printTree();
			Node<BlockTree>* max_t1 = tree.getGlobalMax();
			// std::cout << "[Merge trees] Max from T1: (" << max_t1->printNode() << ")" << std::endl;
			tree.remove(max_t1);
			// std::cout << "[Merge trees] T1 (x <= " << g_beg << ") after removing (" << max_t1->printNodeDetailed() << ")" << std::endl;
			// tree.printTree();
			tree.join(max_t1,t3);
		}
		status = OK;

		// std::cout << "[Merge trees] Resulting tree | Block: " << printBlock() << std::endl;
		// tree.printTree();
		// exit(0);
	}
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
	void initializeTrees();
	
	/* Method used only during construction of the object, 
	to initialize references to nodes. */
	void initializeNodes();
	
	void printNodesBlockDetails();

	// Split the tree from each block into 3 trees t1, t2, t3:
	// - ``t1`` contains all nodes whose arcs (i, i+1) have the ``remote`` 
	// element (i+1) before the reversal (g_beg, g_end];
	// - ``t2`` contains all nodes whose arcs (i, i+1) have the ``remote`` 
	// element (i+1) inside the reversal (g_beg, g_end];
	// - ``t3`` contains all nodes whose arcs (i, i+1) have the ``remote`` 
	// element (i+1) after the reversal (g_beg, g_end];
	// Blocks that were SPLIT at the start of the reversal will be ignored,
	// since their trees are made from scratch after the reversal.
	void updateTreesBeforeReversal(const int g_beg, const int g_end);

	// Merge t1, t2, t3 into a single tree.
	// The ``reversed`` flag of each tree will depend whether
	// the respective block is inside the reversal or not.
	void updateTreesAfterReversal(const int g_beg, const int g_end);

	// Update the trees that were affected by balancing blocks operation.
	// These trees will be made from scratch.
	void updateTreesEndReversal();
	
	void printTrees();
	
	bool hasUnused();

	Node<BlockTree>* getUnusedOriented();

	/* Given an arc, there are four cases to consider:
		Reversal (i, i+1] (normal case). 2 cases:
			1 2 4 -3 5 -> 1 2 3 -4 5
			1 -3 4 2 5 -> 1 -3 -2 -4 5
		Reversal [i, i+1) (positions need to be adjusted). 2 cases:
			1 -2 4 3 5 -> 1 -4 2 3 5
			1 3 4 -2 5 -> 1 -4 -3 -2 5
	*/
	Reversal getReversal(Node<BlockTree>* node) const;

	std::vector<Node<BlockTree>*> getNewSolvedAdjacencies(Reversal rev);

public:
	// Parameterized constructor.
	GenomeSort(const std::vector<int>& perm): genperm(perm){
		// Initialize map of arcs.
		initializeNodes();
		// Initialize trees.
		initializeTrees();
		// Print blocks.
		// std::cout << genperm.printBlocks() << std::endl;
		printTrees();
	}

	// Apply reversal (g_beg, g_end].
	void applyReversal(int g_beg, int g_end);

	/* Sort by reversals.
	Transform permutation A into B by reversal operations.
	- Constraint: all components in the overlap graph of the 
				  permutations A and B must be **oriented**.
	- Constraint: this method treats **a single component**, so 
				  each component in the overlap graph should 
				  be treated separately.
	*/
	std::deque<Reversal> sortByReversals();
};

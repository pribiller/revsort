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

// It makes code a bit clear: 
// instead of std::cout you can type directly cout.
//using namespace std;

/*******************************************************
 * Data structures to keep information on 
 * reversal operation.
 *******************************************************/
class Reversal {
public:
	int g_beg;  // Gene i.
	int g_end;  // Gene i+1.
	int g_beg_next; // Current gene after gene i.
	int g_end_next; // Current gene after gene i+1.
	Reversal(const int start, const int next_ideal, const int next_cur, const int next_end):g_beg(start),g_end(next_ideal),g_beg_next(next_cur),g_end_next(next_end){

	}
};

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
	void initializeTrees(){
		for(BlockTree &b : genperm.blockList) {b.makeTree(nodes);}
	}
	
	/* Method used only during construction of the object, 
	to initialize references to nodes. */
	void initializeNodes(){
		nodes.reserve(genperm.n-1);
		for (int g=0; g<(genperm.n-1); g++) {
			nodes.emplace_back((*genperm.genes[g]), (*genperm.genes[g+1]));
		}
		std::cout << "Nodes:" << std::endl;
		for(Node<BlockTree> &n : nodes) {std::cout << " \t" << n.printNode() << std::endl;}
	}
	
	void printNodesBlockDetails(){
		std::cout << "Nodes:" << std::endl;
		for(Node<BlockTree> &n : nodes) {
			std::cout << " \t" << n.printNodeDetailed() << std::endl;
			std::cout << " \t Block g (" << n.gene.id << ")   = " << n.gene.block->printBlock() << std::endl;
			std::cout << " \t Block g+1 (" << n.gene_next.id << ") = " << n.gene_next.block->printBlock() << std::endl;
		}
	}

	// Split the tree from each block into 3 trees t1, t2, t3:
	// - ``t1`` contains all nodes whose arcs (i, i+1) have the ``remote`` 
	// element (i+1) before the reversal (g_beg, g_end];
	// - ``t2`` contains all nodes whose arcs (i, i+1) have the ``remote`` 
	// element (i+1) inside the reversal (g_beg, g_end];
	// - ``t3`` contains all nodes whose arcs (i, i+1) have the ``remote`` 
	// element (i+1) after the reversal (g_beg, g_end];
	// Blocks that were SPLIT at the start of the reversal will be ignored,
	// since their trees are made from scratch after the reversal.
	void updateTreesBeforeReversal(const int g_beg, const int g_end){
		for(BlockTree &b : genperm.blockList) {
			if(b.status == OK) {b.splitTree(g_beg, g_end);}
		}
	}

	// Merge t1, t2, t3 into a single tree.
	// The ``reversed`` flag of each tree will depend whether
	// the respective block is inside the reversal or not.
	void updateTreesAfterReversal(const int g_beg, const int g_end){

		// std::cout << "After Reversal (" << g_beg << " " << g_end << "]" << std::endl;
		// printNodesBlockDetails();

		for(BlockTree &b : genperm.blockList) {
			if ((b.status.find(SPLIT) != std::string::npos) || (b.status.find(CONC) != std::string::npos)) {
				b.makeTree(nodes); // Make tree from scratch.
			} else {
				b.mergeTrees(g_beg, g_end);	
			}
		}
	}

	// Update the trees that were affected by balancing blocks operation.
	// These trees will be made from scratch.
	void updateTreesEndReversal(){
		for(BlockTree &b : genperm.blockList) {
			if ((b.status.find(SPLIT) != std::string::npos) || (b.status.find(CONC) != std::string::npos)) {
				b.makeTree(nodes); // Make tree from scratch.
			}
		}
		// std::cout << "End Reversal " << std::endl;
		// printNodesBlockDetails();
	}
	
	void printTrees() {
		for(BlockTree &b : genperm.blockList) {std::cout << "\nBlock: " << b.printBlock() << std::endl; b.tree.printTree();}
	}
	
	bool hasUnused() {
		for(BlockTree &b : genperm.blockList) {
			if((b.tree.root != nullptr) && (b.tree.root->unused_tot > 0)){return true;}
		}
		return false;
	}

	Node<BlockTree>* getUnusedOriented() {
		Node<BlockTree>* node = nullptr;
		for(BlockTree &b : genperm.blockList) {
			node = b.tree.getUnusedOriented();
			if(node != nullptr){break;}
		}
		return node;
	}

	/* Given an arc, there are four cases to consider:
		Reversal (i, i+1] (normal case). 2 cases:
			1 2 4 -3 5 -> 1 2 3 -4 5
			1 -3 4 2 5 -> 1 -3 -2 -4 5
		Reversal [i, i+1) (positions need to be adjusted). 2 cases:
			1 -2 4 3 5 -> 1 -4 2 3 5
			1 3 4 -2 5 -> 1 -4 -3 -2 5
	*/
	Reversal getReversal(Node<BlockTree>* node) const {
		const int pos_beg = node->gene.block->genePosAbs(node->gene);
		const int pos_end = node->gene_next.block->genePosAbs(node->gene_next);
		std::list<Gene<BlockTree>>::iterator g_beg = (pos_beg < pos_end) ? genperm.genes[node->gene.id-1] : genperm.genes[node->gene_next.id-1];
		std::list<Gene<BlockTree>>::iterator g_end = (pos_beg < pos_end) ? genperm.genes[node->gene_next.id-1] : genperm.genes[node->gene.id-1];
		// Check if we have the case that needs adjustment: Reversal [i, i+1)
		if (((g_beg->id < g_end->id) && (g_beg->block->reversed != g_beg->reversed)) || ((g_beg->id > g_end->id) && (g_beg->block->reversed == g_beg->reversed))) {
			if (genperm.isFirstGene((*g_beg))){
				std::cout << "ERROR! Reversal " << node->printNode() << " involves the element in the first position of the current permutation. This element is not expected to be involved in any reversal, as it should be in the correct place from the start (the given permutation should be ``expanded`` with the elements 1 and n+2, both never touched by reversals)." << std::endl;
				exit(0);				
			}
			g_beg = g_beg->block->getPrevGene(g_beg);
			g_end = g_end->block->getPrevGene(g_end);
		}
		// Information needed to apply a reversal and also to undo it later if needed.
		const int g_beg_id  = g_beg->id;
		const int g_end_id  = g_end->id;
		std::list<Gene<BlockTree>>::iterator gene_next = g_beg->block->getNextGene(g_beg);
		const int g_beg_next_id = gene_next->id;
		const bool g_end_isLast = genperm.isLastGene((*g_end));
		if(!g_end_isLast){gene_next = g_beg->block->getNextGene(g_end);}
		const int g_end_next_id = (!g_end_isLast) ? gene_next->id : -1;
		// std::cout << "\n\n[getReversal] " << g_beg_id << " " << g_end_id << " " << g_beg_next_id << " " << g_end_next_id << " " << std::endl;
		return Reversal(g_beg_id, g_end_id, g_beg_next_id, g_end_next_id); // std::pair{g_beg, g_end};
	}

	std::vector<Node<BlockTree>*> getNewSolvedAdjacencies(Reversal rev){
		std::vector<int> candidates = {rev.g_beg, rev.g_end, rev.g_beg_next, rev.g_end_next};
		std::vector<Node<BlockTree>*> solvedAdj;
		//solvedAdj.resize(3); // Maximum 2 adjacencies can be solved.
		for(int const &g_id : candidates) {
			if ((g_id > 0) && (g_id < genperm.n)){
				Node<BlockTree>& node = nodes[g_id-1];
				if((node.unused) && node.isSorted()){solvedAdj.emplace_back(&nodes[g_id-1]);} //&nodes[g_id-1]
			}
		}
		return solvedAdj;
	}

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
	void applyReversal(int g_beg, int g_end) {

		g_beg = std::abs(g_beg);
		g_end = std::abs(g_end);

		std::string applyReversal_str = "\n<applyReversal (" + std::to_string(g_beg) + ", " + std::to_string(g_end) + "] >";
		std::cout << applyReversal_str << " Before splitting blocks: " << genperm.printBlocks("\n\t") << std::endl;
		
		// (1) Split at most two blocks so that the endpoints of the reversal correspond to endpoints of blocks;
		genperm.splitBlock(g_beg);
		genperm.splitBlock(g_end);
		// (1.1) After splitting blocks, we must reconstruct from scratch the associated trees.
		std::cout << applyReversal_str << " After splitting blocks: " << genperm.printBlocks("\n\t") << std::endl;

		// (2) Flip ``reversed`` flag of each block between the endpoints of the reversal;
		updateTreesBeforeReversal(g_beg, g_end);
		std::list<BlockTree>::iterator reversal_beg  = std::next(genperm.getBlock(g_beg)); // this block is the first reversed block.
		std::list<BlockTree>::iterator reversal_last = genperm.getBlock(g_end);  // this block is the last reversed block.
		std::list<BlockTree>::iterator reversal_end  = std::next(reversal_last); // this block will not be reversed.
		for (std::list<BlockTree>::iterator b = reversal_beg; b != reversal_end; ++b) { b->reversed = !(b->reversed); b->status += MUT;}
		
		// (3) Reverse the order of the blocks between the endpoints of the reversal;
		const int g_after_beg = reversal_beg->permutationSegment.front().id;
		const int g_after_end = reversal_end->permutationSegment.front().id;
		
		std::list<BlockTree> temp;
		temp.splice(temp.begin(), genperm.blockList, reversal_beg, reversal_end);
		temp.reverse();
		genperm.blockList.splice(reversal_end, temp);

		// (3.1) Update positions of reversed blocks.
		int blockPos = genperm.getBlock(g_beg)->pos + genperm.getBlock(g_beg)->permutationSegment.size();
		for (std::list<BlockTree>::iterator b = genperm.getBlock(g_end); b != std::next(genperm.getBlock(g_after_beg)); ++b) { 
			b->pos = blockPos;
			blockPos += b->permutationSegment.size();
		}

		// (3.2) After reversing the order and ﬂags of the blocks inside the reversal, 
		// every tree, even those outside the reversal, are processed separately.
		updateTreesAfterReversal(g_beg, g_end);

		std::cout << applyReversal_str << " After reversing blocks: " << genperm.printBlocks("\n\t") << std::endl;

		// (4) Concatenate and split blocks in such a way that the size of each block 
		// lies within the interval [½×√(n×log(n)), 2×√(n×log(n))].
		// It must check all 4 blocks that were split.
		genperm.balanceBlock(g_beg);
		genperm.balanceBlock(g_end);
		genperm.balanceBlock(g_after_beg);
		genperm.balanceBlock(g_after_end);
		// (4.1) After concatenating and splitting blocks, reconstruct from scratch the associated trees.
		updateTreesEndReversal();

		std::cout << applyReversal_str << " End of reversal: " << genperm.printBlocks("\n\t") << std::endl;
		printTrees();
	}

	/* Sort by reversals.
	Transform permutation A into B by reversal operations.
	- Constraint: all components in the overlap graph of the 
				  permutations A and B must be **oriented**.
	- Constraint: this method treats **a single component**, so 
				  each component in the overlap graph should 
				  be treated separately.
	*/
	std::deque<Reversal> sortByReversals(){

		std::cout << "\n\nSort by reversals: " << genperm.printBlocks("\n\t") << std::endl;
		printTrees();

		std::deque<Reversal> s1{};
		std::deque<Reversal> s2{};

		Node<BlockTree>* safeReversal = getUnusedOriented();
		while(safeReversal != nullptr){

			Reversal rev = getReversal(safeReversal);
			applyReversal(rev.g_beg, rev.g_end);

			// Update status of arc from ``used`` to ``unused``.
			// Check if the other adjacency affected by the reversal was also ``solved``.
			//safeReversal->gene.block->tree.setUsed(safeReversal);
			std::cout << "\n\nMark used genes" << std::endl;
			std::vector<Node<BlockTree>*> solved = getNewSolvedAdjacencies(rev);
			for (Node<BlockTree>* node : solved) {
				node->gene.block->tree.setUsed(node);
				std::cout << "\t> Mark gene " << node->gene.id << "=USED" << std::endl;
			}
			std::cout << genperm.printBlocks("\n\t") << std::endl;
			printTrees();

			// s1 <- s1 + [new arc].
			s1.push_back(rev);

			// Check if the first element s2 is oriented.
			// If not oriented, go one step back.
			// Apparently this case only happens when the algorithm is about to end.

			// Get another unused oriented arc.
			safeReversal = getUnusedOriented();

			// If there are still ``unused`` arcs.
			if(hasUnused()){
				// while current permutation does not have unused oriented arc.
				while(safeReversal == nullptr){
					// Undo reversal.
					applyReversal(rev.g_beg, rev.g_beg_next);
					rev = s1.back();
					s2.push_front(rev); // s2 <- [new arc] + s2
					s1.pop_back();      // s1 <- s1 - [new arc]
					safeReversal = getUnusedOriented();
				}
			}
		}
		// Merge s2 into s1.
		s1.insert(s1.end(), s2.begin(), s2.end());
		return s1;
	}
};


void testCase1(int const n) {
	// Start vector with Identity permutation (for testing).
	std::vector<int> perm(n);
	std::iota(perm.begin(), perm.end(), 1);	
	GenomeSort genomeSort = GenomeSort(perm);
	genomeSort.applyReversal(2, 9);   // after rev: 1 2 -9 -8 -7 -6   -5 -4 -3 10 11 .. 20
	genomeSort.applyReversal(-6, -5); // after rev: 1 2 -9 -8 -7 -6    5 -4 -3 10 11 .. 20
	genomeSort.applyReversal(-6, 10); // after rev: 1 2 -9 -8 -7 -6  -10  3  4 -5 11 .. 20
}

void testCase2() {
	// Permutation **must** start at 1: [1 2 .. gene]
	std::vector<int> perm{1, -2, 4, 3, 5}; // {0, -1, 3, 2, 4}
	GenomeSort genomeSort = GenomeSort(perm);
	std::deque<Reversal> allrev = genomeSort.sortByReversals();
	std::cout << "Sorting by reversals---Solution" << std::endl;
	for(Reversal const &rev : allrev) {
		std::cout << "(gene " << rev.g_beg << ", gene " << rev.g_end << "]" << std::endl;
	}
}

int main(int argc, char* argv[]) {
	int const n = std::stoi(argv[1]); // input (command line argument): number of genes
	
	//testCase1(n);
	testCase2();

	std::cout << "Bye bye\n";
	return 0;
}

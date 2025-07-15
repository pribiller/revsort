/*******************************************************
 * Original author: Priscila Biller
 * Created: July/2025
 * License: GPL v3
 * 
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
#include <unordered_set>
#include <cmath> 
#include <string>    // convert input args (string to int)
#include <utility>   // move, pair
#include <algorithm> // find_if, reverse
#include <cstdlib>   // exit
#include <deque>
#include <iterator>  // reverse_iterator

#include "sortOrientedByReversals_Tannier2007.hpp"

/*******************************************************
 * Class to find a sorting scenario between two genomes.
 *******************************************************/

/* Method used only during construction of the object, 
to initialize tree of arcs for each block. */
void GenomeSort::initializeTrees(){
	for(BlockTree &b : genperm.blockList) {b.makeTree(nodes);}
}
	
/* Method used only during construction of the object, 
to initialize references to nodes. */
void GenomeSort::initializeNodes(){
	nodes.reserve(genperm.n-1);
	for (int g=0; g<(genperm.n-1); g++) {
		nodes.emplace_back((*genperm.genes[g]), (*genperm.genes[g+1]));
	}
	if (debug) {std::cout << "Nodes:" << std::endl;}
	if (debug) {for(Node<BlockTree> &n : nodes) {std::cout << " \t" << n.printNode() << std::endl;}}
}
	
void GenomeSort::printNodesBlockDetails(){
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
void GenomeSort::updateTreesBeforeReversal(const int g_beg, const int g_end){
	const int g_beg_pos = genperm.getBlock(g_beg)->genePosAbs((*genperm.genes[g_beg-1]));
	const int g_beg_end = genperm.getBlock(g_end)->genePosAbs((*genperm.genes[g_end-1]));
	for(BlockTree &b : genperm.blockList) {
		if(b.status == OK) {b.splitTree(g_beg, g_end, g_beg_pos, g_beg_end);}
	}
}

// Merge t1, t2, t3 into a single tree.
// The ``reversed`` flag of each tree will depend whether
// the respective block is inside the reversal or not.
void GenomeSort::updateTreesAfterReversal(const int g_beg, const int g_end){
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
void GenomeSort::updateTreesEndReversal(){
	for(BlockTree &b : genperm.blockList) {
		if ((b.status.find(SPLIT) != std::string::npos) || (b.status.find(CONC) != std::string::npos)) {
			b.makeTree(nodes); // Make tree from scratch.
		}
	}
	// std::cout << "End Reversal " << std::endl;
	// printNodesBlockDetails();
}
	
void GenomeSort::printTrees() {
	for(BlockTree &b : genperm.blockList) {std::cout << "\nBlock: " << b.printBlock() << std::endl; b.tree.printTree();}
}

bool GenomeSort::hasUnused() {
	for(BlockTree &b : genperm.blockList) {
		if((b.tree.root != nullptr) && (b.tree.root->unused_tot > 0)){return true;}
	}
	return false;
}

Node<BlockTree>* GenomeSort::getUnusedOriented() {
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
Reversal GenomeSort::getReversal(Node<BlockTree>* node) const {
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
	const int g_beg_id  = (g_beg->block->reversed == g_beg->reversed) ? g_beg->id : -g_beg->id;
	const int g_end_id  = (g_end->block->reversed == g_end->reversed) ? g_end->id : -g_end->id;
	std::list<Gene<BlockTree>>::iterator gene_next = g_beg->block->getNextGene(g_beg);
	const int g_beg_next_id = (gene_next->block->reversed == gene_next->reversed) ? gene_next->id : -gene_next->id;
	const bool g_end_isLast = genperm.isLastGene((*g_end));
	if(!g_end_isLast){gene_next = g_end->block->getNextGene(g_end);}
	const int g_end_next_id = (!g_end_isLast) ? ((gene_next->block->reversed == gene_next->reversed) ? gene_next->id : -gene_next->id) : -1;
	// std::cout << "\n\n[getReversal] " << g_beg_id << " " << g_end_id << " " << g_beg_next_id << " " << g_end_next_id << " " << std::endl;
	return Reversal(node->gene.id,g_beg_id, g_end_id, g_beg_next_id, g_end_next_id); // std::pair{g_beg, g_end};
}

std::vector<Node<BlockTree>*> GenomeSort::getNewSolvedAdjacencies(Reversal rev){
	std::unordered_set<int> candidates = {std::abs(rev.g_beg), std::abs(rev.g_end), std::abs(rev.g_beg_next), std::abs(rev.g_end_next)};
	std::vector<Node<BlockTree>*> solvedAdj;
	for(int const &g_id : candidates) {
		// std::cout << "\t[getNewSolvedAdjacencies] Candidate : " << g_id << std::endl;
		if ((g_id > 0) && (g_id < genperm.n)){
			Node<BlockTree>& node = nodes[g_id-1];
			// std::cout << "\t\t[getNewSolvedAdjacencies] unused : " << node.unused << " / isSorted : "  << node.isSorted() << std::endl;
			if((node.unused) && node.isSorted()){solvedAdj.emplace_back(&nodes[g_id-1]);} //&nodes[g_id-1]
		}
	}
	return solvedAdj;
}

// Apply reversal (g_beg, g_end].
void GenomeSort::applyReversal(int g_beg, int g_end) {

	g_beg = std::abs(g_beg);
	g_end = std::abs(g_end);

	std::string applyReversal_str = "\n<applyReversal (" + std::to_string(g_beg) + ", " + std::to_string(g_end) + "] >";
	if (debug) {std::cout << applyReversal_str << " Before splitting blocks: " << genperm.printBlocks("\n\t") << std::endl;}
	
	// (1) Split at most two blocks so that the endpoints of the reversal correspond to endpoints of blocks;
	genperm.splitBlock(g_beg);
	genperm.splitBlock(g_end);
	// (1.1) After splitting blocks, we must reconstruct from scratch the associated trees.
	if (debug) {std::cout << applyReversal_str << " After splitting blocks: " << genperm.printBlocks("\n\t") << std::endl;}

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

	if (debug) {std::cout << applyReversal_str << " After reversing blocks: " << genperm.printBlocks("\n\t") << std::endl;}

	// (4) Concatenate and split blocks in such a way that the size of each block 
	// lies within the interval [½×√(n×log(n)), 2×√(n×log(n))].
	// It must check all 4 blocks that were split.
	genperm.balanceBlock(g_beg);
	genperm.balanceBlock(g_end);
	genperm.balanceBlock(g_after_beg);
	genperm.balanceBlock(g_after_end);
	// (4.1) After concatenating and splitting blocks, reconstruct from scratch the associated trees.
	updateTreesEndReversal();

	if (debug) {
		std::cout << applyReversal_str << " End of reversal: " << genperm.printBlocks("\n\t") << std::endl;
		printTrees();
	}
}

/* Sort by reversals.
Transform permutation A into B by reversal operations.
- Constraint: all components in the overlap graph of the 
			  permutations A and B must be **oriented**.
- Constraint: this method treats **a single component**, so 
			  each component in the overlap graph should 
			  be treated separately.
*/
std::deque<Reversal> GenomeSort::sortByReversals(){

	if (debug) {
		std::cout << "\n\nSort by reversals: " << genperm.printBlocks("\n\t") << std::endl;
		printTrees();
	}

	std::deque<Reversal> s1{};
	std::deque<Reversal> s2{};

	Node<BlockTree>* safeReversal = getUnusedOriented();
	while(safeReversal != nullptr){

		Reversal rev = getReversal(safeReversal);
		applyReversal(rev.g_beg, rev.g_end);

		// Update status of arc from ``used`` to ``unused``.
		// Check if the other adjacency affected by the reversal was also ``solved``.
		//safeReversal->gene.block->tree.setUsed(safeReversal);
		if (debug) {std::cout << "\n\nMark used genes" << std::endl;}
		std::vector<Node<BlockTree>*> solved = getNewSolvedAdjacencies(rev);
		for (Node<BlockTree>* node : solved) {
			node->gene.block->tree.setUsed(node);
			if (debug) {std::cout << "\t> Mark gene " << node->gene.id << "=USED" << std::endl;}
		}
		if (debug) {std::cout << genperm.printBlocks("\n\t") << std::endl;}
		if (debug) {printTrees();}

		// s1 <- s1 + [new arc].
		s1.push_back(rev);
		
		// Check if the first element of s2 is sorted (unoriented **AND** in the correct place).
		// If sorted, either the last element of s1 or the first element of s2 should be popped. 
		// Supposedly, it makes no difference.
		while((!s2.empty()) && (nodes[s2.front().g_arc-1].isSorted())){
			if (debug) {std::cout << "\nRemoving reversal (" << s2.front().g_arc << ") from s2 (reversal is equivalente to last reversal applied)..." << std::endl;}
			// s2 <- -[new arc] +s2 
			s2.pop_front();
		}

		// Get another unused oriented arc.
		safeReversal = getUnusedOriented();

		// If there are still ``unused`` arcs.
		if(hasUnused()){
			// while current permutation does not have unused oriented arc.
			while(safeReversal == nullptr){
				if (debug) {std::cout << std::endl << "Undoing reversal (no safe reversals)... " << rev.printReversal() << std::endl;}
				// Undo reversal.
				rev = s1.back();
				applyReversal(rev.g_beg, rev.g_beg_next);
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

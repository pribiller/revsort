/*******************************************************
 * This class provides some basic utilities for the 
 * sorting by reversals problem.
 * 
 * In particular, it implements a basic function to apply
 * a reversal in a genome represented as a list of blocks 
 * (Kaplan and Verbin's representation) in case no other 
 * additional structures need to be maintained.
 * 
 * For more details on how to represent a genome as 
 * a list of blocks, please check:
 * 
 * - Kaplan and Verbin, 2003: "Efficient Data Structures and a 
 * New Randomized Approach for Sorting Signed Permutations by Reversals"
 * 
 *******************************************************/

#pragma once // It avoids class redefinition.

#include <iostream>
#include <vector>
#include <cstdlib>   // exit
#include <cmath>	 // abs
#include <memory>	 // shared_ptr
#include <utility>   // move, pair

#include "genomePermutation_KaplanVerbin2003.hpp"

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

// Apply reversal (g_beg, g_end].
template <typename BlockT>
void applyReversal(GenomePermutation<BlockTree>& genperm, int g_beg, int g_end) {

	g_beg = std::abs(g_beg);
	g_end = std::abs(g_end);

	std::string applyReversal_str = "\n<applyReversal (" + std::to_string(g_beg) + ", " + std::to_string(g_end) + "] >";
	std::cout << applyReversal_str << " Before splitting blocks: " << genperm.printBlocks("\n\t") << std::endl;
	
	// (1) Split at most two blocks so that the endpoints of the reversal correspond to endpoints of blocks;
	genperm.splitBlock(g_beg);
	genperm.splitBlock(g_end);

	std::cout << applyReversal_str << " After splitting blocks: " << genperm.printBlocks("\n\t") << std::endl;

	// (2) Flip ``reversed`` flag of each block between the endpoints of the reversal;
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

	std::cout << applyReversal_str << " After reversing blocks: " << genperm.printBlocks("\n\t") << std::endl;

	// (4) Concatenate and split blocks in such a way that the size of each block 
	// lies within the interval [½×√(n×log(n)), 2×√(n×log(n))].
	// It must check all 4 blocks that were split.
	genperm.balanceBlock(g_beg);
	genperm.balanceBlock(g_end);
	genperm.balanceBlock(g_after_beg);
	genperm.balanceBlock(g_after_end);
	
	std::cout << applyReversal_str << " End of reversal: " << genperm.printBlocks("\n\t") << std::endl;
}

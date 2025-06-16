/*******************************************************
 * One genome is represented as one permutation. 
 * The permutation is given as a vector of ints,
 * where each int corresponds to a gene, and its sign
 * indicates the direction its transcription occurs.
 * 
 * This class stores this permutation as blocks
 * of size Θ(√n×log(n)), and it makes sure that
 * all blocks will always have their size between 
 * [½×√(n×log(n)), 2×√(n×log(n))], regardless if
 * they are split, concatenated, etc.
 * 
 * For more details:
 * - Eric Tannier et al., 2007: "Advances on sorting by reversals";
 * - Kaplan and Verbin, 2003: "Efficient Data Structures and a 
 * New Randomized Approach for Sorting Signed Permutations by Reversals"
 *******************************************************/

#include <iostream>
#include <list>
#include <vector>
#include <numeric>   // iota
#include <cmath> 
#include <string>    // convert input args (string to int)
#include <utility>   // move
#include <algorithm> // find_if, reverse
#include "genomePermutation.hpp" 

/*******************************************************
 * Implementation of Block class.
 *******************************************************/

std::string Block::printBlock() const {
	std::string block_str = "[ ";
	for(Gene const &g : permutationSegment) {block_str +=  (std::to_string(g.id) + "(" + std::to_string(genePosAbs(g)) + ") ");}
	block_str += ("].rev=" + std::to_string(reversed));
	return block_str;
}

/*******************************************************
 * Implementation of GenomePermutation class.
 *******************************************************/

std::string GenomePermutation::printBlocks(const std::string& sep) const {
	std::string blockList_str;
	for(auto &b : blockList) {blockList_str += (sep + b.printBlock());}
	return blockList_str;
}

/* Adds a new block after the specified position. 
WARNING! ``permSegment`` will become unusable at the end of this method.
*/
std::list<Block>::iterator GenomePermutation::createNewBlock(const int pos, std::list<Gene>& permSegment, std::list<Block>::iterator position){
	// Create a new block containing all genes from ``permSegment``.
	// This new block will be inserted after the block given in the input.
	std::list<Block>::iterator new_block;
	if (position != blockList.end()) {
		blockList.emplace(std::next(position), ++maxBlockId, pos, permSegment);
		// Get an iterator to the new block added to the list.
		new_block = std::next(position);
	} else {
		blockList.emplace_back(++maxBlockId, pos, permSegment); // It return a reference to the block (Block& b)
		// Get an iterator to the new block added to the list.
		new_block = std::prev(blockList.end());
	}
	// Updates the reference block for affected genes.
	for (Gene& g : new_block->permutationSegment) {g.block = new_block;}
	return new_block;
}


void GenomePermutation::splitBlock(const int gene){
	int g = std::abs(gene);

	// Retrieve block where gene is currently located.
	std::list<Gene>::iterator  g_it = genes[g-1];
	std::list<Block>::iterator b_it = g_it->block;

	// Creates an iterator for a potential new block.
	std::list<Block>::iterator new_block; // new_block = b_it;

	// Split depends on the block's flag ``reversed``.

	// Case 1: block->reversed = true. before: A = B + C, rev(A)=true / after: C + B, rev(B) = rev(C) = true.
	if (b_it->reversed){

		// Check if ``gene`` is not the first element in the list (nothing to split).
		if (g_it != b_it->permutationSegment.begin()) {

			// Move genes that appear before `gene` to a new list.
			std::list<Gene> permSegment;
			permSegment.splice(permSegment.begin(), b_it->permutationSegment, b_it->permutationSegment.begin(), g_it);
			// Update relative positions of genes.
			int relPos = 0;
			for (Gene& g : b_it->permutationSegment) {g.pos = relPos++;}

			// Create a new block containing all genes that appear before the gene specified in the input.
			// ``permSegment`` becomes unusable after this point.
			new_block = createNewBlock(b_it->pos + b_it->permutationSegment.size(), permSegment, b_it);
			new_block->reversed = b_it->reversed;
		}

	// Case 2: block->reversed = false. before: A = B + C, rev(A)=false / after: B + C, rev(B) = rev(C) = false.
	} else {
		// Check if ``gene`` is not the last element in the list.
		if ((g_it != b_it->permutationSegment.end()) && (std::next(g_it) != b_it->permutationSegment.end())) {

			// Move genes that appear after `gene` to a new list.
			std::list<Gene> permSegment;
			permSegment.splice(permSegment.begin(), b_it->permutationSegment, std::next(g_it), b_it->permutationSegment.end());

			// Create a new block containing all genes that appear after the gene specified in the input.
			// ``permSegment`` becomes unusable after this point.
			new_block = createNewBlock(b_it->pos+b_it->permutationSegment.size(), permSegment, b_it);
		}
	}
	// Return new block.
	// return new_block;
}

void GenomePermutation::reverseBlock(std::list<Block>::iterator& b){
	// (1) Flip ``reversed`` flag of the block.
	b->reversed = !(b->reversed);
	// (2) Change the order of the elements.
	b->permutationSegment.reverse();
	// (3) Change the sign of the elements and their relative position.
	int relPos = 0;
	for (Gene& g : b->permutationSegment) {g.reversed = !g.reversed; g.pos = relPos++;}
}

/* Block b1 = b1 + b2. Block b2: deleted. */
void GenomePermutation::concatenateBlocks(std::list<Block>::iterator& b1, std::list<Block>::iterator& b2){
	// Update the reference block for affected genes.
	for (Gene& g : b2->permutationSegment) {genes[g.id-1]->block = b1;}
	// Flag ``reversed``: 4 cases to consider.
	// Case 1 and Case 2: Blocks have ``reversed`` flag with the same value.
	if(b1->reversed == b2->reversed){
		// Case 1: b1.reversed = true; b2.reversed = true --> new_b1 = b2 + b1; new_b1.reversed = true;
		if (b1->reversed) {
			// Update relative positions of genes in b1.
			for (Gene& g : b1->permutationSegment) {g.pos += b2->permutationSegment.size();}
			// Move items of ``b2`` to the start of ``b1`` (emptying ``b2`` at the same time).
			(b1->permutationSegment).splice((b1->permutationSegment).begin(), b2->permutationSegment);
		// Case 2: b1.reversed = false; b2.reversed = false --> new_b1 = b1 + b2; new_b1.reversed = false;
		} else {
			// Update relative positions of genes in b2.
			for (Gene& g : b2->permutationSegment) {g.pos += b1->permutationSegment.size();}
			// Move items of ``b2`` to the end of ``b1`` (emptying ``b2`` at the same time).
			(b1->permutationSegment).splice((b1->permutationSegment).end(), b2->permutationSegment);
		}
	// Case 3 and Case 4: Blocks have ``reversed`` flag with a different value.
	} else {
		// Case 3: b1.reversed = true; b2.reversed = false.
		if (b1->reversed) {
			// Case 3.1: b1.reversed = true; b2.reversed = false; |b1| < |b2| --> new_b1 = flip(b1) + b2; new_b1.reversed = false;
			if(b1->permutationSegment.size() < b2->permutationSegment.size()){
				reverseBlock(b1);
				// Update relative positions of genes in b2.
				for (Gene& g : b2->permutationSegment) {g.pos += b1->permutationSegment.size();}
				// Move items of ``b2`` to the end of ``b1`` (emptying ``b2`` at the same time).
				(b1->permutationSegment).splice((b1->permutationSegment).end(), b2->permutationSegment);
			// Case 3.2: b1.reversed = true; b2.reversed = false; |b1| >= |b2| --> new_b1 = flip(b2) + b1; new_b1.reversed = true;
			} else {
				reverseBlock(b2);
				// Update relative positions of genes in b1.
				for (Gene& g : b1->permutationSegment) {g.pos += b2->permutationSegment.size();}
				// Move items of ``b2`` to the start of ``b1`` (emptying ``b2`` at the same time).
				(b1->permutationSegment).splice((b1->permutationSegment).begin(), b2->permutationSegment);
			}
		// Case 4: b1.reversed = false; b2.reversed = true.
		} else {
			// Case 4.1: b1.reversed = false; b2.reversed = true; |b1| < |b2| --> new_b1 = b2 + flip(b1); new_b1.reversed = true;
			if(b1->permutationSegment.size() < b2->permutationSegment.size()){
				reverseBlock(b1);
				// Update relative positions of genes in b1.
				for (Gene& g : b1->permutationSegment) {g.pos += b2->permutationSegment.size();}
				// Move items of ``b2`` to the start of ``b1`` (emptying ``b2`` at the same time).
				(b1->permutationSegment).splice((b1->permutationSegment).begin(), b2->permutationSegment);
			// Case 4.2: b1.reversed = false; b2.reversed = true; |b1| >= |b2| --> new_b1 = b1 + flip(b2); new_b1.reversed = false;
			} else {
				reverseBlock(b2);
				// Update relative positions of genes in b2.
				for (Gene& g : b2->permutationSegment) {g.pos += b1->permutationSegment.size();}
				// Move items of ``b2`` to the end of ``b1`` (emptying ``b2`` at the same time).
				(b1->permutationSegment).splice((b1->permutationSegment).end(), b2->permutationSegment);
			}
		}
	}
	// Delete block ``b2``.
	blockList.erase(b2);
}

/* Makes sure that the size of the block where gene g belongs
lies within the interval [½×√(n×log(n)), 2×√(n×log(n))]. */
void GenomePermutation::balanceBlock(const int gene){

	std::list<Block>::iterator b = getBlock(std::abs(gene));

	// If the block is too small, then concatenate the block with a neighboring block.
	while(b->permutationSegment.size() < minBlockSize){

		std::cout << "\t [balanceBlock] Block needs to increase in size: cur.size=" <<  b->permutationSegment.size() << " / minBlockSize=" << minBlockSize << " / " << b->printBlock() << std::endl;

		std::list<Block>::iterator b_prev = std::prev(b);
		std::list<Block>::iterator b_next = std::next(b);
		const int b_prev_size = (b != blockList.begin())    ? b_prev->permutationSegment.size() : (n+1);
		const int b_next_size = (b_next != blockList.end()) ? b_next->permutationSegment.size() : (n+1);
		// Concatenate blocks.
		if(b_prev_size < b_next_size){
			std::cout << "\t [balanceBlock] Concatenate blocks: " << b_prev->printBlock() << " and " << b->printBlock() << std::endl;
			concatenateBlocks(b_prev, b);
			b = b_prev; // b was erased; update reference.
		} else {
			std::cout << "\t [balanceBlock] Concatenate blocks: " << b->printBlock() << " and " << b_next->printBlock() << std::endl;
			concatenateBlocks(b, b_next);
		}
		std::cout << "\t [balanceBlock] After concatenation: " << printBlocks() << std::endl;
	}
	
	// If the block is too big, then split the block into two blocks.
	if(b->permutationSegment.size() > maxBlockSize) {
		std::cout << "\t [balanceBlock] Block needs to reduce in size: cur.size=" <<  b->permutationSegment.size() << " / maxBlockSize=" << maxBlockSize << " / " << b->printBlock() << std::endl;

		std::list<Gene>::iterator g_mid_it = b->permutationSegment.begin();
		std::advance(g_mid_it, std::floor(b->permutationSegment.size()/2));
		const int g_mid = g_mid_it->id;
		splitBlock(g_mid);
	}
}

/* Method used only during construction of the object, 
to initialize references to genes. */
void GenomePermutation::initializeGenes(){
	genes.resize(n);
	for(auto &b : blockList) {
		for (std::list<Gene>::iterator gene_it = b.permutationSegment.begin(); gene_it != b.permutationSegment.end(); ++gene_it) {
			int g = gene_it->id; // Access the first element of the iterator.
			genes[g-1] = gene_it;
		}
	}
}

/* Method used only during construction of the object, 
to initialize list of blocks. */
void GenomePermutation::initializeBlocks(const std::vector<int>& perm){
	// Split permutation into blocks.
	int blockSize = std::round(std::sqrt(n*std::log(n)));		
	int const numChunks		= n / blockSize; // number of blocks.
	int const baseChunkSize = n / numChunks;
	int const remainder		= n % numChunks;

	// Create a list of blocks with size Θ(√n×log(n)).
	int pos = 0;
	for (int i=0; i<numChunks; i++) {
		// Take a segment of the permutation.
		int const currentChunkSize = baseChunkSize + (i < remainder ? 1 : 0);

		std::list<Gene> permSegment;
		for (int size = 0; size < currentChunkSize; ++size, ++pos) {
			permSegment.emplace_back(std::abs(perm[pos]), size, (perm[pos] < 0));
		}

		// Make a block with the permutation segment.
		// ``permSegment`` becomes unusable after this point.
		createNewBlock(pos-currentChunkSize+1,permSegment, blockList.end());
	}
}

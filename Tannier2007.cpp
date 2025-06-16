#include <iostream>
#include <list>
#include <vector>
#include <numeric> // iota
#include <cmath> 
#include <string> // convert input args (string to int)
#include <utility> // move
#include <algorithm> // find_if, reverse
#include "genome.hpp" 
//#include <memory> // shared_ptr
//#include <chrono>

#include <boost/intrusive/rbtree.hpp>

// It makes code a bit clear: 
// instead of std::cout you can type directly cout.
//using namespace std;

/*******************************************************
 * Auxiliary data structures.
 *******************************************************/
inline int Node::getPosNext() const {
	return gene_next.block->genePosAbs(gene_next);
}

/*******************************************************
 * Class to find a sorting scenario between two genomes.
 *******************************************************/

class GenomeSort {
private:
	std::list<Block> blockList;
	std::vector<std::list<Gene>::iterator> genes; // references of genes (size=n genes).
	std::vector<Node> nodes; // references of nodes in balanced binary trees for each block (size=n-1 arcs).
	int n; // number of genes

	// Invariant: at *all* points during the sorting algorithm, 
	// the block sizes are between [½×√(n×log(n)), 2×√(n×log(n))].
	int minBlockSize;
	int maxBlockSize;

	/* Adds a new block after the specified position. 
	WARNING! ``permSegment`` will become unusable at the end of this method.
	*/
	std::list<Block>::iterator createNewBlock(const int pos, std::list<Gene>& permSegment, std::list<Block>::iterator position){
		// Create a new block containing all genes from ``permSegment``.
		// This new block will be inserted after the block given in the input.
		std::list<Block>::iterator new_block;
		if (position != blockList.end()) {
			blockList.emplace(std::next(position), pos, permSegment);
			// Get an iterator to the new block added to the list.
			new_block = std::next(position);
		} else {
			blockList.emplace_back(pos, permSegment); // It return a reference to the block (Block& b)
			// Get an iterator to the new block added to the list.
			new_block = std::prev(blockList.end());
		}
		// Updates the reference block for affected genes.
		for (Gene& g : new_block->permutationSegment) {g.block = new_block;}
		return new_block;
	}

	void splitBlock(const int gene){
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

	void reverseBlock(std::list<Block>::iterator& b){
		// (1) Flip ``reversed`` flag of the block.
		b->reversed = !(b->reversed);
		// (2) Change the order of the elements.
		b->permutationSegment.reverse();
		// (3) Change the sign of the elements and their relative position.
		int relPos = 0;
		for (Gene& g : b->permutationSegment) {g.reversed = !g.reversed; g.pos = relPos++;}
	}

	/* Block b1 = b1 + b2. Block b2: deleted. */
	void concatenateBlocks(std::list<Block>::iterator& b1, std::list<Block>::iterator& b2){
		// Update the reference block for affected genes.
		for (Gene& g : b2->permutationSegment) {g.block = b1;}
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
	void balanceBlock(const int gene){

		std::list<Block>::iterator b = genes[std::abs(gene)-1]->block;

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

	void createTree(Block& block){
		block.tree.clear();
		for (std::list<Gene>::iterator gene = block.permutationSegment.begin(); gene != block.permutationSegment.end(); ++gene){
			if(gene->id < nodes.size()){ // The last gene does not have an arc.
				std::cout << "Insert node = " << gene->id << "\n";
				block.tree.insert_equal(nodes[gene->id-1]);				
			}
		}
	}

	/* Method used only during construction of the object, 
	to initialize tree of arcs for each block. */
	void initializeTrees(){
		for(auto &b : blockList) {createTree(b);}
	}

	/* Method used only during construction of the object, 
	to initialize references to nodes. */
	void initializeNodes(){
		nodes.reserve(n-1);
		for (int g=0; g<(n-1); g++) {nodes.emplace_back((*genes[g]), (*genes[g+1]));}
	}
	
	/* Method used only during construction of the object, 
	to initialize references to genes. */
	void initializeGenes(){
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
	void initializeBlocks(const std::vector<int>& perm){
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

public:
	// Parameterized constructor.
	GenomeSort(const std::vector<int>& perm):n(perm.size()),minBlockSize(std::floor(0.5*(std::sqrt(n*std::log(n))))),maxBlockSize(std::ceil(2*(std::sqrt(n*std::log(n))))){
		std::cout << "min block size=" << minBlockSize << " max block size=" << maxBlockSize << "\n";
		// Split permutation into blocks.
		// Create a list of blocks with size Θ(√n×log(n)). 
		initializeBlocks(perm);
		// Initialize map of genes.
		initializeGenes();
		// Initialize map of arcs.
		initializeNodes();
		// Initialize trees.
		initializeTrees();
		// Print blocks.
		std::cout << printBlocks() << std::endl;
	}

	std::string printBlocks(const std::string& sep=" ") const {
		std::string blockList_str;
		for(auto &b : blockList) {blockList_str += (sep + b.printBlock());}
		return blockList_str;
	}
		
	void applyReversal(int g_beg, int g_end){

		g_beg = std::abs(g_beg);
		g_end = std::abs(g_end);

		std::string applyReversal_str = "<applyReversal (" + std::to_string(g_beg) + ", " + std::to_string(g_end) + "] >";
		
		std::cout << applyReversal_str << " Before splitting blocks: " << printBlocks("\n\t") << std::endl;

		// (1) Split at most two blocks so that the endpoints of the reversal correspond to endpoints of blocks;
		splitBlock(g_beg);
		splitBlock(g_end);

		std::cout << applyReversal_str << " After splitting blocks: " << printBlocks("\n\t") << std::endl;

		// (2) Flip ``reversed`` flag of each block between the endpoints of the reversal;
		std::list<Block>::iterator reversal_beg  = std::next(genes[g_beg-1]->block); // this block is the first reversed block.
		std::list<Block>::iterator reversal_last = genes[g_end-1]->block;    // this block is the last reversed block.
		std::list<Block>::iterator reversal_end  = std::next(reversal_last); // this block will not be reversed.
		for (std::list<Block>::iterator b = reversal_beg; b != reversal_end; ++b) { b->reversed = !(b->reversed);}

		// (3) Reverse the order of the blocks between the endpoints of the reversal;
		const int g_after_beg = reversal_beg->permutationSegment.front().id;
		const int g_after_end = reversal_end->permutationSegment.front().id;
		std::reverse(reversal_beg, reversal_end); // Reverses the order of the elements in the range [first, last).
		// (3.1) Update positions of reversed blocks.
		int blockPos = (genes[g_beg-1]->block)->pos + (genes[g_beg-1]->block)->permutationSegment.size();
		for (std::list<Block>::iterator b = genes[g_end-1]->block; b != std::next(genes[g_after_beg-1]->block); ++b) { 
			b->pos = blockPos;
			blockPos += b->permutationSegment.size();
		}

		std::cout << applyReversal_str << " After reversing blocks: " << printBlocks("\n\t") << std::endl;

		// (4) Concatenate and split blocks in such a way that the size of each block 
		// lies within the interval [½×√(n×log(n)), 2×√(n×log(n))].
		// It must check all 4 blocks that were split.
		balanceBlock(g_beg);
		balanceBlock(g_end);
		balanceBlock(g_after_beg);
		balanceBlock(g_after_end);
	}

	/* Temporary function to handle memory of boost trees. */
	void clearTrees(){
		for(auto &b : blockList) {b.tree.clear();}
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
	genomeSort.applyReversal(2, 9);   // after rev: 1 2 -9 -8 -7 -6   -5 -4 -3 10 11 .. 20
	genomeSort.applyReversal(-6, -5); // after rev: 1 2 -9 -8 -7 -6    5 -4 -3 10 11 .. 20
	genomeSort.applyReversal(-6, 10); // after rev: 1 2 -9 -8 -7 -6  -10  3  4 -5 11 .. 20

	genomeSort.clearTrees();
	std::cout << "Bye bye\n";
	return 0;
}

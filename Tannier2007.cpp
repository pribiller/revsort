#include <iostream>
#include <list>
#include <vector>
#include <numeric> // iota
#include <cmath> 
#include <string> // convert input args (string to int)
#include <utility> // move
#include <algorithm> // find_if, reverse
//#include <memory> // shared_ptr
//#include <chrono>


// It makes code a bit clear: 
// instead of std::cout you can type directly cout.
//using namespace std;

/*******************************************************
 * Auxiliary data structures.
 *******************************************************/

class Block {
public:

	int blockId;
	std::list<int> permutationSegment;
	bool reversed{false}; // raised if the block should be read in the reverse order, changing the sign of the elements.

	/* p is a **rvalue reference** (&&) to a ``permutation segment`` created externally 
	(i.e., a list of ints with type std::list<int>). The parameter ``p`` itself is a **lvalue**.
	If the attribute ``permutationSegment`` is constructed using ``permutationSegment(p)``, 
	then the list in ``permutationSegment`` is **copy constructed** from p (not very efficient, references to genes need to be updated). 
	In order to **move constructed** instead (no extra memory allocated, references are kept intact), 
	we use instead ``permutationSegment(std::move(p))``.
	After this operation, the parameter ``p`` becomes in an unspecified state and should no longer be used.
	*/
	Block(int id, std::list<int>& p_tmp):blockId(id),permutationSegment(std::move(p_tmp)){ 
	}

	void printBlock(){
		for(auto const &g : permutationSegment) {std::cout << g << " ";}
		std::cout.put('\n');
	}
};

/*******************************************************
 * Class to find a sorting scenario between two genomes.
 *******************************************************/

class GenomeSort {
private:
	std::list<Block> blockList;
	std::vector<std::list<Block>::iterator> genes_to_blocks; // map between genes and blocks.
	std::vector<std::list<int>::iterator> genes; // references of genes.
	int blockId_max{0};
	int n; // number of genes

	// Invariant: at *all* points during the sorting algorithm, 
	// the block sizes are between [½×√(n×log(n)), 2×√(n×log(n))].
	int minBlockSize;
	int maxBlockSize;

	/* Adds a new block after the specified position. 
	WARNING! ``permSegment`` will become unusable at the end of this method.
	*/
	std::list<Block>::iterator createNewBlock(std::list<int>& permSegment, std::list<Block>::iterator position){
		// Create a new block containing all genes from ``permSegment``.
		// This new block will be inserted after the block given in the input.
		std::list<Block>::iterator new_block;
		if (position != blockList.end()) {
			blockList.emplace(std::next(position), ++blockId_max, permSegment);
			// Get an iterator to the new block added to the list.
			new_block = std::next(position);
		} else {
			blockList.emplace_back(++blockId_max, permSegment); // It return a reference to the block (Block& b)
			// Get an iterator to the new block added to the list.
			new_block = std::prev(blockList.end());
		}
		// Updates the reference block for affected genes.
		for (const int& g : new_block->permutationSegment) {genes_to_blocks[std::abs(g)] = new_block;}
		return new_block;
	}

	std::list<Block>::iterator splitBlock(const int gene){
		int g = std::abs(gene);

		// Retrieve block where gene is currently located.
		std::list<Block>::iterator b_it = genes_to_blocks[g];
		std::list<int>::iterator g_it   = genes[g];

		// Creates an iterator for a potential new block.
		std::list<Block>::iterator new_block;

		// Check if ``gene`` is not the last element in the list.
		if ((g_it != b_it->permutationSegment.end()) && (std::next(g_it) != b_it->permutationSegment.end())) {

			// Move genes that appear after `gene` to a new list.
			std::list<int> permSegment;
			permSegment.splice(permSegment.begin(), b_it->permutationSegment, std::next(g_it), b_it->permutationSegment.end());

			// Create a new block containing all genes that appear after the gene specified in the input.
			// ``permSegment`` becomes unusable after this point.
			// Position of the new block will depend on the flag ``reversed`` of current block.
			// reversed=true -> The new block will appear *before* the current block in the block list.
			if(b_it->reversed){
				new_block = createNewBlock(permSegment, ((b_it != blockList.begin()) ? std::prev(b_it) : blockList.begin()));
			// reversed=false -> The new block will appear *next* the current block in the block list.
			} else {
				new_block = createNewBlock(permSegment, b_it);
			}
			new_block->reversed = b_it->reversed;
		// If ``gene`` is the last element in the list, there is nothing to be split.
		} else {
			new_block = b_it;
		}
		// Return new block.
		return new_block;
	}

	void reverseBlock(std::list<Block>::iterator& b){
		// (1) Flip ``reversed`` flag of the block.
		b->reversed = !(b->reversed);
		// (2) Change the order of the elements.
		b->permutationSegment.reverse();
		// (3) Change the sign of the elements.
		for (int& g : b->permutationSegment) {g = -g;}
	}

	/* Block b1 = b1 + b2. Block b2: deleted. */
	void concatenateBlocks(std::list<Block>::iterator& b1, std::list<Block>::iterator& b2){
		// Update the reference block for affected genes.
		for (const int& g : b2->permutationSegment) {genes_to_blocks[std::abs(g)] = b1;}
		// Flag ``reversed``: 4 cases to consider.
		// Case 1 and Case 2: Blocks have ``reversed`` flag with the same value.
		if(b1->reversed == b2->reversed){
			// Case 1: b1.reversed = true; b2.reversed = true --> new_b1 = b2 + b1; new_b1.reversed = true;
			if (b1->reversed) {
				// Move items of ``b2`` to the start of ``b1`` (emptying ``b2`` at the same time).
				(b1->permutationSegment).splice((b1->permutationSegment).begin(), b2->permutationSegment);
			// Case 2: b1.reversed = false; b2.reversed = false --> new_b1 = b1 + b2; new_b1.reversed = false;
			} else {
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
					// Move items of ``b2`` to the end of ``b1`` (emptying ``b2`` at the same time).
					(b1->permutationSegment).splice((b1->permutationSegment).end(), b2->permutationSegment);
				// Case 3.2: b1.reversed = true; b2.reversed = false; |b1| >= |b2| --> new_b1 = flip(b2) + b1; new_b1.reversed = true;
				} else {
					reverseBlock(b2);
					// Move items of ``b2`` to the start of ``b1`` (emptying ``b2`` at the same time).
					(b1->permutationSegment).splice((b1->permutationSegment).begin(), b2->permutationSegment);
				}
			// Case 4: b1.reversed = false; b2.reversed = true.
			} else {
				// Case 4.1: b1.reversed = false; b2.reversed = true; |b1| < |b2| --> new_b1 = b2 + flip(b1); new_b1.reversed = true;
				if(b1->permutationSegment.size() < b2->permutationSegment.size()){
					reverseBlock(b1);
					// Move items of ``b2`` to the start of ``b1`` (emptying ``b2`` at the same time).
					(b1->permutationSegment).splice((b1->permutationSegment).begin(), b2->permutationSegment);
				// Case 4.2: b1.reversed = false; b2.reversed = true; |b1| >= |b2| --> new_b1 = b1 + flip(b2); new_b1.reversed = false;
				} else {
					reverseBlock(b2);
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

		std::list<Block>::iterator b = genes_to_blocks[std::abs(gene)];

		// If the block is too small, then concatenate the block with a neighboring block.
		while(b->permutationSegment.size() < minBlockSize){

			std::cout << "Block needs to increase in size: cur.size=" <<  b->permutationSegment.size() << " / minBlockSize=" << minBlockSize << " / ";
			b->printBlock();

			std::list<Block>::iterator b_prev = std::prev(b);
			std::list<Block>::iterator b_next = std::next(b);
			const int b_prev_size = (b != blockList.begin())    ? b_prev->permutationSegment.size() : (n+1);
			const int b_next_size = (b_next != blockList.end()) ? b_next->permutationSegment.size() : (n+1);
			// Concatenate blocks.
			if(b_prev_size < b_next_size){
				std::cout << "Concatenate blocks: " << b_prev_size << ", " <<  b_next_size <<"\n";
				b_prev->printBlock();
				std::cout << "and \n";
				b->printBlock();
				concatenateBlocks(b_prev, b);
				b = b_prev; // b was erased; update reference.
			} else {
				std::cout << "Concatenate blocks: " << b_prev_size << ", " <<  b_next_size <<"\n";
				b->printBlock();
				std::cout << "and \n";
				b_next->printBlock();
				concatenateBlocks(b, b_next);
			}
			std::cout << "After concatenation\n";
			printBlocks();
		}
		
		// If the block is too big, then split the block into two blocks.
		if(b->permutationSegment.size() > maxBlockSize){

		}
		//
	}

	/* Method used only during construction of the object, 
	to initialize references to genes. */
	void initializeGenes(){
		genes.resize(n);
		for(auto &b : blockList) {
			for (auto gene_it = b.permutationSegment.begin(); gene_it != b.permutationSegment.end(); ++gene_it) {
				int g = std::abs(*gene_it); // Access the first element of the iterator.
				genes[g] = gene_it;
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
		genes_to_blocks.resize(n);
		int begIdx = 0;
		for (int i=0; i<numChunks; i++) {
			// Take a segment of the permutation.
			int const currentChunkSize = baseChunkSize + (i < remainder ? 1 : 0);
			std::list<int> permSegment( perm.begin() + begIdx, 
									perm.begin() + begIdx + currentChunkSize);
			// Make a block with the permutation segment.
			// ``permSegment`` becomes unusable after this point.
			createNewBlock(permSegment, blockList.end());

			// Update position in the permutation.
			begIdx += currentChunkSize;
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
		// Print blocks.
		printBlocks();
	}

	void printBlocks(){
		for(auto &b : blockList) {b.printBlock();}
		std::cout.put('\n');
	}
		
	void applyReversal(const int g_beg, const int g_end){
		
		std::cout << "Before splitting blocks (" << g_beg << ", " << g_end << "): \n";
		printBlocks();


		// (1) Split at most two blocks so that the endpoints of the reversal correspond to endpoints of blocks;
		splitBlock(g_beg);
		splitBlock(g_end);

		std::cout << "After splitting blocks: \n";
		printBlocks();

		// (2) Flip ``reversed`` flag of each block between the endpoints of the reversal;
		std::list<Block>::iterator reversal_beg = std::next(genes_to_blocks[std::abs(g_beg)]);
		std::list<Block>::iterator reversal_end = std::next(genes_to_blocks[std::abs(g_end)]);
		for (auto b = reversal_beg; b != reversal_end; ++b) { b->reversed = !(b->reversed);}

		// (3) Reverse the order of the blocks between the endpoints of the reversal;
		const int g_after_beg = reversal_beg->permutationSegment.front();
		const int g_after_end = reversal_end->permutationSegment.front();
		std::reverse(reversal_beg, reversal_end); // Reverses the order of the elements in the range [first, last).

		std::cout << "After reversing blocks: \n";
		printBlocks();

		// (4) Concatenate and split blocks in such a way that the size of each block 
		// lies within the interval [½×√(n×log(n)), 2×√(n×log(n))].
		// It must check all 4 blocks that were split.
		balanceBlock(g_beg);
		balanceBlock(g_end);
		balanceBlock(g_after_beg);
		balanceBlock(g_after_end);
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
	genomeSort.applyReversal(2, 9);
	// genomeSort.applyReversal(4, 5);
	// genomeSort.applyReversal(4, 10);

	return 0;
}

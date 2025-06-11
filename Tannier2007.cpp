#include <iostream>
#include <list>
#include <vector>
#include <numeric> // iota
#include <cmath> 
#include <string> // convert input args (string to int)
#include <algorithm> // find_if
//#include <memory> // shared_ptr
#include <chrono>
#include <utility> // move

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
	Block(int id, std::list<int>&& p_tmp):blockId(id),permutationSegment(std::move(p_tmp)){ 
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
	std::list<Block> blocks;
	std::vector<std::list<Block>::iterator> genes_to_blocks; // map between genes and blocks.
	std::vector<std::list<int>::iterator> genes; // references of genes.
	int blockId_max{0};
	int n; // number of genes

	/* Adds a new block after the specified position. 
	WARNING! ``permSegment`` will become unusable at the end of this method.
	*/
	void createNewBlock(std::list<int>&& permSegment, std::list<Block>::iterator position){
		// Create a new block containing all genes from ``permSegment``.
		// This new block will be inserted after the block given in the input.
		std::list<Block>::iterator new_block;
		if (position != blocks.end()) {
			blocks.emplace(std::next(position), ++blockId_max, std::move(permSegment));
			// Get an iterator to the new block added to the list.
			new_block = std::next(position);
		} else {
			blocks.emplace_back(++blockId_max, std::move(permSegment)); // It return a reference to the block (Block& b)
			// Get an iterator to the new block added to the list.
			new_block = std::prev(blocks.end());
		}
		// Updates the reference block for affected genes.
	    for (const int& g : new_block->permutationSegment) {genes_to_blocks[g] = new_block;}

	    std::cout << "New block=";
		new_block->printBlock();
	}

	std::list<Block>::iterator splitBlock(int gene){
		int g = std::abs(gene);

		// Retrieve block where gene is currently located.
		std::list<Block>::iterator b_it = genes_to_blocks[g];
		std::list<int>::iterator g_it   = genes[g];

		// Print debug info on the retrieved block.
		std::cout << "Gene 1: " << g << "; Block 1: ";
		b_it->printBlock();
		std::cout << "Gene 1 (from ref): " << (*g_it) << ";\n";
		std::cout << "Gene 1 (next): " << (*std::next(g_it)) << ";\n";

		// Creates an iterator for a potential new block.
		std::list<Block>::iterator b_new_it = b_it;

		// Check if gene is not the last element in the list.
		if (g_it != b_it->permutationSegment.end()) {

			std::list<int>::iterator g_next = std::next(g_it);
			
			std::cout << "Gene next (before splice): " << *(std::next(g_next)) << ";\n";

			// Move genes that appear after `gene` to a new list.
			std::list<int> permSegment;
			permSegment.splice(permSegment.begin(), b_it->permutationSegment, std::next(g_it), b_it->permutationSegment.end());

			std::cout << "Gene next (after splice): " << *(std::next(g_next)) << ";\n";

			// Create a new block containing all genes that appear after the gene specified in the input.
			// This new block will appear *next* the current block in the block list.
			// ``permSegment`` becomes unusable after this point.
			createNewBlock(std::move(permSegment), b_it);

			std::cout << "Gene next (after new block): " << *(std::next(g_next)) << ";\n";
		}
		// Return new block.
		return b_new_it;
	}

	/* Method used only during construction of the object, 
	to initialize references to genes. */
	void initializeGenes(){
		genes.resize(n);
		for(auto &b : blocks) {
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
		int const numChunks	    = n / blockSize; // number of blocks.
		int const baseChunkSize = n / numChunks;
		int const remainder	    = n % numChunks;

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
			createNewBlock(std::move(permSegment), blocks.end());

			// Update position in the permutation.
			begIdx += currentChunkSize;
		}
	}

public:
	// Parameterized constructor.
	GenomeSort(const std::vector<int>& perm):n(perm.size()){
		// Split permutation into blocks.
		// Create a list of blocks with size Θ(√n×log(n)). 
		initializeBlocks(perm);
		// Initialize map of genes.
		initializeGenes();
		// Print blocks.
		printBlocks();
	}

	void printBlocks(){
		for(auto &b : blocks) {b.printBlock();}
		std::cout.put('\n');
	}
		
	void applyReversal(int g1, int g2){
		splitBlock(g1);
		printBlocks();
		splitBlock(g2);
		printBlocks();
		// Reverse segment.
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

	return 0;
}

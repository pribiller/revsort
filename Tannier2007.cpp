#include <iostream>
#include <list>
#include <vector>
#include <numeric> // iota
#include <cmath> 
#include <string> // convert input args (string to int)
#include <algorithm> // find_if
//#include <memory> // shared_ptr
#include <chrono>

// It makes code a bit clear: 
// instead of std::cout you can type directly cout.
//using namespace std;

/***********************
 * Utils
 ***********************/

// Splits a vector into chunks of size approximately [chunkSize].
std::vector<std::vector<int>> splitVector(const std::vector<int>& original, int idealChunkSize) {
	std::vector<std::vector<int>> result;

	int const numChunks	 = original.size() / idealChunkSize; // number of chunks.
	int const baseChunkSize = original.size() / numChunks;
	int const remainder	 = original.size() % numChunks;

	// At each step creates a new sub-list.
	int begIdx = 0;
	for (int i=0; i<numChunks; i++) {
		int const currentChunkSize = baseChunkSize + (i < remainder ? 1 : 0);
		// Create a sub-vector from the current position to the next N elements.
		std::vector<int> subVector( original.begin() + begIdx, 
									original.begin() + begIdx + currentChunkSize);
		// Add the new list at the end.
		result.push_back(subVector);
		begIdx += currentChunkSize;
	}
	return result;
}

// Splits a vector into chunks of size approximately [chunkSize].
std::vector<std::list<int>> splitVectorIntoLists(const std::vector<int>& original, int idealChunkSize) {
	std::vector<std::list<int>> result;

	int const numChunks	 = original.size() / idealChunkSize; // number of chunks.
	int const baseChunkSize = original.size() / numChunks;
	int const remainder	 = original.size() % numChunks;

	// At each step creates a new sub-list.
	int begIdx = 0;
	for (int i=0; i<numChunks; i++) {
		int const currentChunkSize = baseChunkSize + (i < remainder ? 1 : 0);
		// Create a sub-vector from the current position to the next N elements.
		std::list<int> subList( original.begin() + begIdx, 
								original.begin() + begIdx + currentChunkSize);
		// Add the new list at the end.
		result.push_back(subList);
		begIdx += currentChunkSize;
	}
	return result;
}

/*******************************************************
 * Auxiliary data structures.
 *******************************************************/

class Block {
public:

	int blockId;
	std::list<int> block;
	bool reversed{false}; // raised if the block should be read in the reverse order, changing the sign of the elements.

	Block(int id, const std::list<int>& b):blockId(id),block(b){
	}
	
	void printBlock(){
		for(auto const &g : block) {std::cout << g << " ";}
		std::cout.put('\n');
	}
};

/*******************************************************
 * Class to find a sorting scenario between two genomes.
 *******************************************************/

class GenomeSort {
private:
	std::list<Block> blocks;
	int blockId_max;
	std::vector<std::list<Block>::iterator> genes_to_blocks; // map between genes and blocks.
	std::vector<std::list<int>::iterator> genes; // references of genes.

	std::list<Block>::iterator getBlock(int gene){
		return genes_to_blocks[std::abs(gene)];
	}

	std::list<Block>::iterator splitBlock(int gene){

		int g = std::abs(gene);

		// Retrieve block where gene is currently located.
		std::list<Block>::iterator b_it = genes_to_blocks[g];
		std::list<int>::iterator g_it   = genes[g];

		// Print debug info on the retrieved block.
		std::cout << "Gene 1: " << g << "; Block 1: ";
		(*b_it).printBlock();
		std::cout << "Gene 1: " << (*g_it) << "\n";
		std::cout << "Gene 2: " << (*std::next(g_it)) << "\n";

		// Creates an iterator for a potential new block.
		std::list<Block>::iterator b_new_it = b_it;

		// Check if gene is not the last element in the list.
		if (g_it != (*b_it).block.end()) {

			std::cout << "Not the last element\n";
			
			// Move genes that appear after `gene` to a new list.
			std::list<int> new_b;
			new_b.splice(new_b.begin(), (*b_it).block, std::next(g_it), (*b_it).block.end());

			std::cout << "Splice done\n";
			(*b_it).printBlock();
			// Create a new block containing all genes that appear after the gene specified in the input.
			// This new block will appear next the current block in the block list.

			// Updates the reference block for affected genes.
		}
		// Return new block.
		return b_new_it;
	}

public:
	// Parameterized constructor.
	GenomeSort(const std::vector<int>& perm){

		// Split permutation into blocks of size Θ(√n×log(n)).
		int n = perm.size();
		int blockSize = std::round(std::sqrt(n*std::log(n)));
		std::vector<std::list<int>> permList = splitVectorIntoLists(perm, blockSize);

		std::cout << "n=" << perm.size() << " -> Block size " << blockSize << "." << std::endl;

		// Create a map between genes and blocks.
		// WARNING! It assumes that (absolute) elements in the permutation [perm] are between 0 and n-1.
		genes.resize(n);
		genes_to_blocks.resize(n);

		// Create a list of blocks.
		for (int i = 0; i < permList.size(); ++i) {
			Block& b = blocks.emplace_back(i, permList[i]); // It return a reference to the block (Block&)
			// Get an iterator to the last block added to the list.
			std::list<Block>::iterator new_block = std::prev(blocks.end());
			// Update map between genes and blocks.
			for (auto gene_it = b.block.begin(); gene_it != b.block.end(); ++gene_it) {
				int g = std::abs(*gene_it); // Access the first element of the iterator.
				genes_to_blocks[g] = new_block;
				genes[g] = gene_it;
			}
		}
		blockId_max = permList.size();
		
		// Print blocks.
		printBlocks();
	}

	void printBlocks(){
		for(auto &b : blocks) {b.printBlock();}
		std::cout.put('\n');
	}
		
	void applyReversal(int g1, int g2){
		splitBlock(g1);
		//splitBlock(g2);
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

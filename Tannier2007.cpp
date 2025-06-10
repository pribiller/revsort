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

/*******************************************************
 * Auxiliary data structures.
 *******************************************************/

class Block {
public:

	int blockId;
	std::vector<int> block;
	bool reversed; // raised if the block should be read in the reverse order, changing the sign of the elements.

	Block(int id, std::vector<int> b){
		blockId  = id;
		block	 = b;
		reversed = false;
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

	void splitBlock(std::list<Block>::iterator b, int g){
		
		b->block;
		
		return ;
	}

public:
	// Parameterized constructor.
	GenomeSort(const std::vector<int>& perm){

		// Split permutation into blocks of size Θ(√n×log(n)).
		int n = perm.size();
		int blockSize = std::round(std::sqrt(n*std::log(n)));
		std::vector<std::vector<int>> permList = splitVector(perm, blockSize);

		std::cout << "n=" << perm.size() << " -> Block size " << blockSize << "." << std::endl;

		// Create a map between genes and blocks.
		// WARNING! It assumes that (absolute) elements in the permutation [perm] are between 0 and n-1.
		genes_to_blocks.resize(n);

		// Create a list of blocks.
		for (int i = 0; i < permList.size(); ++i) {
			blocks.emplace_back(i, permList[i]); // It return a reference to the block (Block&)
			// Get an iterator to the last block added to the list.
			std::list<Block>::iterator new_block = std::prev(blocks.end());
			// Update map between genes and blocks.
			for(auto const &g : permList[i]) {genes_to_blocks[std::abs(g)] = new_block;}
		}
		blockId_max = permList.size();

		// Print blocks.
		printBlocks();
	}

	void printBlocks(){
		for(auto &b : blocks) {b.printBlock();}
		std::cout.put('\n');
	}

	std::list<Block>::iterator getBlock(int gene){
		return genes_to_blocks[std::abs(gene)];
	}
	
	void applyReversal(int g1, int g2){

		// Retrieve blocks where genes are.
		std::list<Block>::iterator b1 = getBlock(g1);
		std::list<Block>::iterator b2 = getBlock(g2);

		std::cout << "Gene 1: " << g1 << "; Block 1: ";
		(*b1).printBlock();
		std::cout << "Gene 2: " << g2 << "; Block 2: ";
		(*b2).printBlock();

		// Split blocks.

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

	return 0;
}

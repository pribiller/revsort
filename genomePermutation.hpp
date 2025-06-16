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

/*******************************************************
 * Auxiliary data structures.
 *******************************************************/

template <typename BlockT>
class BlockBase; // forward declaration

template <typename BlockT>
class Gene {
public:
	int id;  // Absolute value of the gene.
	int pos; // Relative position inside the block.
	bool reversed{false};              // Orientation of the gene (+ or -).
	typename std::list<BlockT>::iterator block; // Block where gene currently is.
	Gene(const int id_, const int pos_, const bool sign_):id(id_),pos(pos_),reversed(sign_){
	}
};

template <typename BlockT>
class BlockBase {
public:
	int id;
	std::list<Gene<BlockT>> permutationSegment;
	bool reversed{false}; // Raised if the block should be read in the reverse order, changing the sign of the elements.
	int pos;              // Starting index of the permutation segment in the current permutation.

	/* p is a **rvalue reference** (&&) to a ``permutation segment`` created externally 
	(i.e., a list of ints with type std::list<int>). The parameter ``p`` itself is a **lvalue**.
	If the attribute ``permutationSegment`` is constructed using ``permutationSegment(p)``, 
	then the list in ``permutationSegment`` is **copy constructed** from p (not very efficient, references to genes need to be updated). 
	In order to **move constructed** instead (no extra memory allocated, references are kept intact), 
	we use instead ``permutationSegment(std::move(p))``.
	After this operation, the parameter ``p`` becomes in an unspecified state and should no longer be used.
	*/
	BlockBase(const int id, const int posBeg, std::list<Gene<BlockT>>& p_tmp):id(id),pos(posBeg),permutationSegment(std::move(p_tmp)){ 
		// Update relative position of genes in the new block.
		int relPos = 0;
		for(Gene<BlockT> &g : permutationSegment) {g.pos=relPos++;}
	}
	
	/* It computes the absolute position of the gene based on its 
	relative position and the absolute starting position of the block.
	*/
	inline int genePosAbs(Gene<BlockT> const &g) const {
		return (reversed ? (pos + permutationSegment.size() - 1 -g.pos) : (pos + g.pos));
	}

	std::string printBlock() const;

	// TODO: updateBlock method that is called during splitTree?

};

/***********************************************************
 * Genome represented as a permutation divided into blocks.
 ***********************************************************/
template <typename BlockT>
class GenomePermutation {
protected:
	std::list<BlockT> blockList;
	std::vector<typename std::list<Gene<BlockT>>::iterator> genes;  // references of genes (size=n genes).
	int n; // number of genes
		
	int maxBlockId{0};

	// Invariant: at *all* points during the sorting algorithm, 
	// the block sizes are between [½×√(n×log(n)), 2×√(n×log(n))].
	int minBlockSize;
	int maxBlockSize;

	/* Adds a new block after the specified position. 
	WARNING! ``permSegment`` will become unusable at the end of this method.
	*/
	typename std::list<BlockT>::iterator createNewBlock(const int pos, std::list<Gene<BlockT>>& permSegment, typename std::list<BlockT>::iterator position);

	void splitBlock(const int gene);

	void reverseBlock(typename std::list<BlockT>::iterator& b);

	/* Block b1 = b1 + b2. Block b2: deleted. */
	void concatenateBlocks(typename std::list<BlockT>::iterator& b1, typename std::list<BlockT>::iterator& b2);

	/* Makes sure that the size of the block where gene g belongs
	lies within the interval [½×√(n×log(n)), 2×√(n×log(n))]. */
	void balanceBlock(const int gene);

	/* Method used only during construction of the object, 
	to initialize references to genes. */
	void initializeGenes();

	/* Method used only during construction of the object, 
	to initialize list of blocks. */
	void initializeBlocks(const std::vector<int>& perm);

	inline typename std::list<BlockT>::iterator& getBlock(const int gene){return genes[gene-1]->block;} 

public:
	// Parameterized constructor.
	GenomePermutation(const std::vector<int>& perm):n(perm.size()),minBlockSize(std::floor(0.5*(std::sqrt(n*std::log(n))))),maxBlockSize(std::ceil(2*(std::sqrt(n*std::log(n))))){
		// Split permutation into blocks.
		// Create a list of blocks with size Θ(√n×log(n)). 
		initializeBlocks(perm);
		// Initialize map of genes.
		initializeGenes();
	}

	std::string printBlocks(const std::string& sep=" ") const;
};

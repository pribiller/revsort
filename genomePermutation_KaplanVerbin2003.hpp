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

#pragma once // It avoids class redefinition.

#include <iostream>
#include <list>
#include <vector>
#include <unordered_map>
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

// Status of blocks.
const std::string OK    =  ""; // No events happened.
const std::string SPLIT = "S"; // Block split.
const std::string CONC  = "C"; // Block concatenation.
const std::string MUT   = "M"; // Block inside mutation.

template <typename BlockT>
class BlockBase {
public:
	int id;
	std::list<Gene<BlockT>> permutationSegment;
	bool reversed{false};   // Raised if the block should be read in the reverse order, changing the sign of the elements.
	int pos;                // Starting index of the permutation segment in the current permutation.
	std::string status{OK}; // Status of a block. It keeps track if the block was modified, or affected by an external event.

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

	// Get next gene in the current permutation.
	// Behavior is undefined if given gene is last gene.
	typename std::list<Gene<BlockT>>::iterator getNextGene(typename std::list<Gene<BlockT>>::iterator const &g_it) {
		// Check if ``gene`` is the last element of the block.
		const bool isLastElement = ((reversed  && (g_it == g_it->block->permutationSegment.begin()))
								|| (!reversed && ((g_it == g_it->block->permutationSegment.end()) || (std::next(g_it) == g_it->block->permutationSegment.end()))));

		// Get the block of the next gene depending if the 
		// given gene is the last element of the block or not.
		typename std::list<BlockT>::iterator b_it = isLastElement ? std::next(g_it->block) : g_it->block;
		if (b_it->reversed) {
			return (isLastElement ? (--b_it->permutationSegment.end()) : std::prev(g_it));
		} else {
			return (isLastElement ? b_it->permutationSegment.begin() : std::next(g_it));
		}
	}

	// Get previous gene in the current permutation.
	// Behavior is undefined if given gene is first gene.
	typename std::list<Gene<BlockT>>::iterator getPrevGene(typename std::list<Gene<BlockT>>::iterator const &g_it) {
		// Check if ``gene`` is the first element of the block.
		const bool isFirstElement = ((reversed  && ((g_it == g_it->block->permutationSegment.end()) || (std::next(g_it) == g_it->block->permutationSegment.end()))))
								 || (!reversed  && (g_it == g_it->block->permutationSegment.begin()));

		// Get the block of the next gene depending if the 
		// given gene is the last element of the block or not.
		typename std::list<BlockT>::iterator b_it = isFirstElement ? std::prev(g_it->block) : g_it->block;
		if (b_it->reversed) {
			return (isFirstElement ? b_it->permutationSegment.begin() : std::next(g_it));
		} else {
			return (isFirstElement ? (--b_it->permutationSegment.end()) : std::prev(g_it));
		}
	}

	std::string printBlock() const;
};

/***********************************************************
 * Genome represented as a permutation divided into blocks.
 ***********************************************************/
template <typename BlockT>
class GenomePermutation {
protected:	
	int maxBlockId{0};

	// Invariant: at *all* points during the sorting algorithm, 
	// the block sizes are between [½×√(n×log(n)), 2×√(n×log(n))].
	int minBlockSize{0};
	int maxBlockSize{0};

	/* Adds a new block after the specified position. 
	WARNING! ``permSegment`` will become unusable at the end of this method.
	*/
	typename std::list<BlockT>::iterator createNewBlock(const int pos, std::list<Gene<BlockT>>& permSegment, typename std::list<BlockT>::iterator position);

	void reverseBlock(typename std::list<BlockT>::iterator& b);

	/* Block b1 = b1 + b2. Block b2: deleted. */
	void concatenateBlocks(typename std::list<BlockT>::iterator& b1, typename std::list<BlockT>::iterator& b2);

	/* Method used only during construction of the object, 
	to initialize references to genes. */
	void initializeGenes();

	/* Method used only during construction of the object, 
	to initialize list of blocks. */
	void initializeBlocks(const std::vector<int>& perm);

	/* Auxiliary method for ``getExtendedPerm``. 
	It sorts gene extremities by their position in the permutation.
	This function affects the vector passed as input.
	It returns a mapping between the gene extremity and its position in the sorted vector. */
	std::unordered_map<int, int> sortGeneExtremities(std::vector<int>& gene_extremities);

public:
	int n; // number of genes
	
	std::list<BlockT> blockList;
	std::vector<typename std::list<Gene<BlockT>>::iterator> genes;  // references of genes (size=n genes).

	// Parameterized constructor.
	GenomePermutation(const std::vector<int>& perm):n(perm.size()){
		
		minBlockSize = std::floor(0.5*(std::sqrt(n*std::log(n))));
		maxBlockSize = std::ceil(2*(std::sqrt(n*std::log(n))));

		// Split permutation into blocks.
		// Create a list of blocks with size Θ(√n×log(n)). 
		initializeBlocks(perm);
		// Initialize map of genes.
		initializeGenes();
	}

	void splitBlock(const int gene);

	std::string printBlocks(const std::string& sep=" ") const;

	std::vector<int> getExtendedPerm() const;
	std::vector<int> getExtendedPerm(std::vector<int>& gene_extremities);
	std::vector<int> getUnextendedPerm() const;
	std::vector<int> getUnsignedExtendedPerm() const;

	inline typename std::list<BlockT>::iterator& getBlock(const int gene){return genes[gene-1]->block;} 

	/* Makes sure that the size of the block where gene g belongs
	lies within the interval [½×√(n×log(n)), 2×√(n×log(n))]. */
	void balanceBlock(const int gene);

	inline bool isFirstGene(Gene<BlockT> const &gene) const {
		return (gene.block->genePosAbs(gene) == 1);
	}

	inline bool isLastGene(Gene<BlockT> const &gene) const {
		return (gene.block->genePosAbs(gene) == n);
	}
};

/*******************************************************/
/*******************************************************/
/*******************************************************/

/*******************************************************
 * Implementation of Block class.
 *******************************************************/

template <typename BlockT>
std::string BlockBase<BlockT>::printBlock() const {
	std::string block_str  = "[ ";
	const std::string status_str = ((status == "") ? "" : (" status: " + status));
	for(Gene<BlockT> const &g : permutationSegment) {block_str +=  ((g.reversed ? "-" : "+") + std::to_string(g.id) + "(" + std::to_string(genePosAbs(g)) + ") ");}
	block_str += ("].rev=" + std::to_string(reversed)) + status_str;
	return block_str;
}

/*******************************************************
 * Implementation of GenomePermutation class.
 *******************************************************/
template <typename BlockT>
std::string GenomePermutation<BlockT>::printBlocks(const std::string& sep) const {
	std::string blockList_str;
	for(auto &b : blockList) {blockList_str += (sep + b.printBlock());}
	return blockList_str;
}

template <typename BlockT>
std::vector<int> GenomePermutation<BlockT>::getExtendedPerm() const {
	std::vector<int> perm(n);
	for(auto &gene : genes) {
		const int pos = gene->block->genePosAbs((*gene))-1;
		perm[pos] = gene->id;
		if (gene->reversed != gene->block->reversed){perm[pos] = -perm[pos];}
	}
	return perm;
}

template <typename BlockT>
std::unordered_map<int, int> GenomePermutation<BlockT>::sortGeneExtremities(std::vector<int>& gene_extremities) {
	// Find positions of gene extremities in the current permutation.
	// - Labels in the current permutation: 1, 2, .., n (1 and n are extended);
	// - Labels in the ``gene extremities`` representation: 0, 1, 2, 3, 4, .., 2*n, 2*n+1 (0 and 2n+1 are extended).
	std::unordered_map<int, int> gene_to_pos_map;
	for(const int &g_ext : gene_extremities) {
		const int gene_id  = (int)std::trunc((g_ext+1)/2.0)+1;
		Gene<BlockT>& gene = (*genes[gene_id-1]);
		const int gene_pos = gene.block->genePosAbs(gene);
		int g_ext_pos;
		// Gene sign is ``+``, order of gene extremities: 2i-1, 2i
		if (gene.reversed == gene.block->reversed) {
			g_ext_pos = (((g_ext % 2) == 0) ? (2*gene_pos) : 2*gene_pos-1);
		// Gene sign is ``-``, order of gene extremities: 2i, 2i-1
		} else {
			g_ext_pos = (((g_ext % 2) == 0) ? (2*gene_pos-1) : 2*gene_pos);
		}
		gene_to_pos_map[g_ext] = g_ext_pos;
	}
	// Sort gene extremities by their position in the current permutation.
	auto comparePos = [&gene_to_pos_map](int a, int b) {return gene_to_pos_map[a] < gene_to_pos_map[b];};
	std::sort(gene_extremities.begin(), gene_extremities.end(), comparePos);
	// Update the map between gene extremities and their positions.
	int cur_idx = 0;
	for(const int g_ext: gene_extremities) {gene_to_pos_map[g_ext] = cur_idx++;}
	return gene_to_pos_map;
}

// It creates an extended signed ``sub-permutation``, composed 
// only of gene extremities specified in the input.
template <typename BlockT>
std::vector<int> GenomePermutation<BlockT>::getExtendedPerm(std::vector<int>& gene_extremities) {
	std::vector<int> perm(gene_extremities.size()/2+1);
	// Sort gene extremities by their position in the current permutation.
	std::unordered_map<int,int> gene_to_pos_map = sortGeneExtremities(gene_extremities);

	std::cout << "Sorted gene extremities: " << std::endl;
	for(const int& g: gene_extremities){std::cout << g << " ";}
	std::cout << std::endl;
	
	// Re-label gene extremities.
	int g_ext = gene_extremities[0];
	int g_ext_pos = 0;
	int g_label = 1; // New label for gene.
	int g_pos   = 0;
	for (int i=0; i<perm.size()-1; i++) {
		perm[g_pos] = g_label;
		// std::cout << "i=" << i << ": g_pos=" << g_pos << " / g_label=" << g_label << std::endl;
		// Move through a gray edge.
		// A gray edge has its start in an even value and its end in an odd value.	
		// std::cout << "i=" << i << ": Moves through gray edge: " << g_ext << "[" << g_ext_pos << "]";
		g_ext     = (((g_ext % 2) == 0) ? (g_ext+1) : (g_ext-1));
		g_ext_pos = gene_to_pos_map[g_ext];
		// std::cout << " -> " << g_ext << "[" << g_ext_pos << "]" << std::endl;
		// Update new gene label.
		g_label = std::abs(g_label)+1;
		if((g_ext_pos % 2) == 0){g_label=-g_label;}
		// Update gene extremity.
		const bool lastGeneExt = ((g_ext_pos == gene_extremities.size()-1) && ((g_ext_pos % 2) == 1));
		if (!lastGeneExt) {
			g_ext = ((g_ext_pos % 2) == 0) ? gene_extremities[g_ext_pos-1] : gene_extremities[g_ext_pos+1];
			g_ext_pos = gene_to_pos_map[g_ext];
		}
		// Update new gene position.
		g_pos = (int)std::trunc((g_ext_pos+1)/2.0);
	}
	// Re-label last gene.
	perm[g_pos] = g_label;
	return perm;
}

template <typename BlockT>
std::vector<int> GenomePermutation<BlockT>::getUnextendedPerm() const {
	// Removes the extended elements ``1`` and ``n``.
	std::vector<int> perm(n-2); 
	// Check if extended elements (``1`` and ``n``) are in their positions.
	typename std::list<Gene<BlockT>>::iterator first = genes[0];
	typename std::list<Gene<BlockT>>::iterator last  = genes[genes.size()-1];
	if ((first->block->genePosAbs((*first)) != 1) || (last->block->genePosAbs((*last)) != n)){
		std::cout << "ERROR! Extended elements were moved:\n- Position of 1 (should be 1): " << first->block->genePosAbs((*first)) << "\n- Position of n (should be " << n << "): " << last->block->genePosAbs((*last)) << std::endl;
		return perm;
	}
	for(auto &gene : genes) {
		const int pos = gene->block->genePosAbs((*gene))-2;
		if((pos >= 0) && (pos < perm.size())){
			perm[pos] = gene->id-1;
			if (gene->reversed != gene->block->reversed){perm[pos] = -perm[pos];}
		}
	}
	return perm;
}

// Creates an extended unsigned permutation, where a gene i (1 <= i <= n)
// is represented by its gene extremities i and i+1.
// Two additional gene extremities are added: 0 and 2n+1.
// Thus, the total set of extremities are {0, 1, 2, ..., 2n, 2n+1}.
template <typename BlockT>
std::vector<int> GenomePermutation<BlockT>::getUnsignedExtendedPerm() const {
	std::vector<int> unsignedExtPerm(2*n-2);
	// Extended elements.
	unsignedExtPerm[0]     = 0;
	unsignedExtPerm[2*n-3] = 2*n-3;
	// Main elements.
	int idx_perm = 1;
	for(auto &gene : genes) {
		// Skip ``extended`` elements (1 and n).
		if((gene->id > 1) && (gene->id < n)){
			const int pos_2 = 2*(gene->block->genePosAbs((*gene))-1);
			const int pos_1 = pos_2-1;
			const int ext_2 = (gene->id-1)*2;
			const int ext_1 = ext_2-1;
			// Gene sign is "+"
			if (gene->reversed == gene->block->reversed){
				unsignedExtPerm[pos_1] = ext_1;
				unsignedExtPerm[pos_2] = ext_2;
			// Gene sign is "-"
			} else {
				unsignedExtPerm[pos_2] = ext_1;
				unsignedExtPerm[pos_1] = ext_2;
			}
		}
	}
	return unsignedExtPerm;
}

/* Adds a new block after the specified position. 
WARNING! ``permSegment`` will become unusable at the end of this method.
*/
template <typename BlockT>
typename std::list<BlockT>::iterator GenomePermutation<BlockT>::createNewBlock(const int pos, std::list<Gene<BlockT>>& permSegment, typename std::list<BlockT>::iterator position){
	// Create a new block containing all genes from ``permSegment``.
	// This new block will be inserted after the block given in the input.
	typename std::list<BlockT>::iterator new_block;
	if (position != blockList.end()) {
		new_block = blockList.emplace(std::next(position), ++maxBlockId, pos, permSegment);
		// Get an iterator to the new block added to the list.
		//new_block = std::next(position);
	} else {
		blockList.emplace_back(++maxBlockId, pos, permSegment); // It return a reference to the block (Block& b)
		// Get an iterator to the new block added to the list.
		new_block = std::prev(blockList.end());
	}
	// Updates the reference block for affected genes.
	for (Gene<BlockT>& g : new_block->permutationSegment) {g.block = new_block;}
	return new_block;
}

/* Block is split into [1 2 .. gene] [gene+1 .. n]*/
// replace by pair (or tuple)
template <typename BlockT>
void GenomePermutation<BlockT>::splitBlock(const int gene){
	int g = std::abs(gene);

	// Retrieve block where gene is currently located.
	typename std::list<Gene<BlockT>>::iterator g_it = genes[g-1];
	typename std::list<BlockT>::iterator b_it = g_it->block;

	// Creates an iterator for a potential new block.
	typename std::list<BlockT>::iterator new_block; // new_block = b_it;

	// Split depends on the block's flag ``reversed``.

	// Case 1: block->reversed = true. before: A = B + C, rev(A)=true / after: C + B, rev(B) = rev(C) = true.
	if (b_it->reversed){

		// Check if ``gene`` is not the first element in the list (nothing to split).
		if (g_it != b_it->permutationSegment.begin()) {

			// Move genes that appear before `gene` to a new list.
			std::list<Gene<BlockT>> permSegment;
			permSegment.splice(permSegment.begin(), b_it->permutationSegment, b_it->permutationSegment.begin(), g_it);
			// Update relative positions of genes.
			int relPos = 0;
			for (Gene<BlockT>& g : b_it->permutationSegment) {g.pos = relPos++;}

			// Create a new block containing all genes that appear before the gene specified in the input.
			// ``permSegment`` becomes unusable after this point.
			new_block = createNewBlock(b_it->pos + b_it->permutationSegment.size(), permSegment, b_it);
			new_block->reversed = b_it->reversed;

			// Update status of blocks.
			b_it->status      += SPLIT;
			new_block->status += SPLIT;
		}

	// Case 2: block->reversed = false. before: A = B + C, rev(A)=false / after: B + C, rev(B) = rev(C) = false.
	} else {
		// Check if ``gene`` is not the last element in the list.
		if ((g_it != b_it->permutationSegment.end()) && (std::next(g_it) != b_it->permutationSegment.end())) {

			// Move genes that appear after `gene` to a new list.
			std::list<Gene<BlockT>> permSegment;
			permSegment.splice(permSegment.begin(), b_it->permutationSegment, std::next(g_it), b_it->permutationSegment.end());

			// Create a new block containing all genes that appear after the gene specified in the input.
			// ``permSegment`` becomes unusable after this point.
			new_block = createNewBlock(b_it->pos+b_it->permutationSegment.size(), permSegment, b_it);

			// Update status of blocks.
			b_it->status      += SPLIT;
			new_block->status += SPLIT;
		}
	}
	// Return affected blocks.
	// return {(*b_it), (*new_block)};
}

template <typename BlockT>
void GenomePermutation<BlockT>::reverseBlock(typename std::list<BlockT>::iterator& b){
	// (1) Flip ``reversed`` flag of the block.
	b->reversed = !(b->reversed);
	// (2) Change the order of the elements.
	b->permutationSegment.reverse();
	// (3) Change the sign of the elements and their relative position.
	int relPos = 0;
	for (Gene<BlockT>& g : b->permutationSegment) {g.reversed = !g.reversed; g.pos = relPos++;}
}

/* Block b1 = b1 + b2. Block b2: deleted. */
template <typename BlockT>
void GenomePermutation<BlockT>::concatenateBlocks(typename std::list<BlockT>::iterator& b1, typename std::list<BlockT>::iterator& b2){
	// Update the reference block for affected genes.
	for (Gene<BlockT>& g : b2->permutationSegment) {getBlock(g.id) = b1;}
	// Flag ``reversed``: 4 cases to consider.
	// Case 1 and Case 2: Blocks have ``reversed`` flag with the same value.
	if(b1->reversed == b2->reversed){
		// Case 1: b1.reversed = true; b2.reversed = true --> new_b1 = b2 + b1; new_b1.reversed = true;
		if (b1->reversed) {
			// Update relative positions of genes in b1.
			for (Gene<BlockT>& g : b1->permutationSegment) {g.pos += b2->permutationSegment.size();}
			// Move items of ``b2`` to the start of ``b1`` (emptying ``b2`` at the same time).
			(b1->permutationSegment).splice((b1->permutationSegment).begin(), b2->permutationSegment);
		// Case 2: b1.reversed = false; b2.reversed = false --> new_b1 = b1 + b2; new_b1.reversed = false;
		} else {
			// Update relative positions of genes in b2.
			for (Gene<BlockT>& g : b2->permutationSegment) {g.pos += b1->permutationSegment.size();}
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
				for (Gene<BlockT>& g : b2->permutationSegment) {g.pos += b1->permutationSegment.size();}
				// Move items of ``b2`` to the end of ``b1`` (emptying ``b2`` at the same time).
				(b1->permutationSegment).splice((b1->permutationSegment).end(), b2->permutationSegment);
			// Case 3.2: b1.reversed = true; b2.reversed = false; |b1| >= |b2| --> new_b1 = flip(b2) + b1; new_b1.reversed = true;
			} else {
				reverseBlock(b2);
				// Update relative positions of genes in b1.
				for (Gene<BlockT>& g : b1->permutationSegment) {g.pos += b2->permutationSegment.size();}
				// Move items of ``b2`` to the start of ``b1`` (emptying ``b2`` at the same time).
				(b1->permutationSegment).splice((b1->permutationSegment).begin(), b2->permutationSegment);
			}
		// Case 4: b1.reversed = false; b2.reversed = true.
		} else {
			// Case 4.1: b1.reversed = false; b2.reversed = true; |b1| < |b2| --> new_b1 = b2 + flip(b1); new_b1.reversed = true;
			if(b1->permutationSegment.size() < b2->permutationSegment.size()){
				reverseBlock(b1);
				// Update relative positions of genes in b1.
				for (Gene<BlockT>& g : b1->permutationSegment) {g.pos += b2->permutationSegment.size();}
				// Move items of ``b2`` to the start of ``b1`` (emptying ``b2`` at the same time).
				(b1->permutationSegment).splice((b1->permutationSegment).begin(), b2->permutationSegment);
			// Case 4.2: b1.reversed = false; b2.reversed = true; |b1| >= |b2| --> new_b1 = b1 + flip(b2); new_b1.reversed = false;
			} else {
				reverseBlock(b2);
				// Update relative positions of genes in b2.
				for (Gene<BlockT>& g : b2->permutationSegment) {g.pos += b1->permutationSegment.size();}
				// Move items of ``b2`` to the end of ``b1`` (emptying ``b2`` at the same time).
				(b1->permutationSegment).splice((b1->permutationSegment).end(), b2->permutationSegment);
			}
		}
	}
	// Update status of the blocks.
	b1->status += CONC;
	b2->status += CONC;
	// Delete block ``b2``.
	blockList.erase(b2);
}

/* Makes sure that the size of the block where gene g belongs
lies within the interval [½×√(n×log(n)), 2×√(n×log(n))]. */
template <typename BlockT>
void GenomePermutation<BlockT>::balanceBlock(const int gene){

	typename std::list<BlockT>::iterator b = getBlock(std::abs(gene));

	// If the block is too small, then concatenate the block with a neighboring block.
	while(b->permutationSegment.size() < minBlockSize){

		std::cout << "\t [balanceBlock] Block needs to increase in size: cur.size=" <<  b->permutationSegment.size() << " / minBlockSize=" << minBlockSize << " / " << b->printBlock() << std::endl;

		typename std::list<BlockT>::iterator b_prev = std::prev(b);
		typename std::list<BlockT>::iterator b_next = std::next(b);
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

		typename std::list<Gene<BlockT>>::iterator g_mid_it = b->permutationSegment.begin();
		std::advance(g_mid_it, std::floor(b->permutationSegment.size()/2));
		const int g_mid = g_mid_it->id;
		splitBlock(g_mid);
	}
}

/* Method used only during construction of the object, 
to initialize references to genes. */
template <typename BlockT>
void GenomePermutation<BlockT>::initializeGenes(){
	genes.resize(n);
	for(auto &b : blockList) {
		for (typename std::list<Gene<BlockT>>::iterator gene_it = b.permutationSegment.begin(); gene_it != b.permutationSegment.end(); ++gene_it) {
			int g = gene_it->id; // Access the first element of the iterator.
			genes[g-1] = gene_it;
		}
	}
}

/* Method used only during construction of the object, 
to initialize list of blocks. */
template <typename BlockT>
void GenomePermutation<BlockT>::initializeBlocks(const std::vector<int>& perm){
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

		std::list<Gene<BlockT>> permSegment;
		for (int size = 0; size < currentChunkSize; ++size, ++pos) {
			permSegment.emplace_back(std::abs(perm[pos]), size, (perm[pos] < 0));
		}

		// Make a block with the permutation segment.
		// ``permSegment`` becomes unusable after this point.
		createNewBlock(pos-currentChunkSize+1,permSegment, blockList.end());
	}
}

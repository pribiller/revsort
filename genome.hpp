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

#include <boost/intrusive/rbtree.hpp>

/*******************************************************
 * Auxiliary data structures.
 *******************************************************/

class Block; // forward declaration

class Gene {
public:
	int id;  // Absolute value of the gene.
	int pos; // Relative position inside the block.
	bool reversed{false};             // Orientation of the gene (+ or -).
	std::list<Block>::iterator block; // Block where gene currently is.
	Gene(const int id_, const int pos_, const bool sign_):id(id_),pos(pos_),reversed(sign_){
	}
};

/* Each node corresponds to an arc v_i. */
                    // Base hook with default options (i.e., optimized for speed).
class Node : public boost::intrusive::set_base_hook<> {
public:
	const Gene& gene;
	const Gene& gene_next;
	int orientedTotal{0}; // the number of oriented arcs in the subtree rooted at node i.
	bool reversed{false}; // nodes should be ordered backwards with respect to the original order, and all nodes change orientation.

	Node(const Gene& g, const Gene& g_next):gene(g),gene_next(g_next){ 
	}

	// Returns the position of node i's successor in the current permutation.
	inline int getPosNext() const;

	// Returns if the arc i is oriented or not (oriented: sign of i and (i+1) are flipped).
	inline bool getOrientation() const {
		return (gene.reversed != gene_next.reversed);
	}

	// Overload the '<' operator.
	friend bool operator<(const Node &a, const Node &b){
		std::cout << "Comp <\n";
		return a.getPosNext() < b.getPosNext();
	}
	// Overload the '>' operator.
	friend bool operator>(const Node &a, const Node &b){
		std::cout << "Comp >\n";
		return a.getPosNext() > b.getPosNext();
	}
	// Overload the '==' operator.
	friend bool operator==(const Node &a, const Node &b){
		std::cout << "Comp ==\n";
		return a.getPosNext() == b.getPosNext();
	}
};

class Block {
public:
	
	std::list<Gene> permutationSegment;
	bool reversed{false};                // Raised if the block should be read in the reverse order, changing the sign of the elements.
	int pos;                             // Starting index of the permutation segment in the current permutation.
	boost::intrusive::rbtree<Node> tree; // Balanced binary tree with segment elements as nodes.

	/* p is a **rvalue reference** (&&) to a ``permutation segment`` created externally 
	(i.e., a list of ints with type std::list<int>). The parameter ``p`` itself is a **lvalue**.
	If the attribute ``permutationSegment`` is constructed using ``permutationSegment(p)``, 
	then the list in ``permutationSegment`` is **copy constructed** from p (not very efficient, references to genes need to be updated). 
	In order to **move constructed** instead (no extra memory allocated, references are kept intact), 
	we use instead ``permutationSegment(std::move(p))``.
	After this operation, the parameter ``p`` becomes in an unspecified state and should no longer be used.
	*/
	Block(const int posBeg, std::list<Gene>& p_tmp):pos(posBeg),permutationSegment(std::move(p_tmp)){ 
		// Update relative position of genes in the new block.
		int relPos = 0;
		for(Gene &g : permutationSegment) {g.pos=relPos++;}
	}
	
	/* It computes the absolute position of the gene based on its 
	relative position and the absolute starting position of the block.
	*/
	inline int genePosAbs(Gene const &g) const {
		return (reversed ? (pos + permutationSegment.size() - 1 -g.pos) : (pos + g.pos));
	}

	std::string printBlock() const {
		std::string block_str = "[ ";
		for(Gene const &g : permutationSegment) {block_str +=  (std::to_string(g.id) + "(" + std::to_string(genePosAbs(g)) + ") ");}
		block_str += ("].rev=" + std::to_string(reversed));
		return block_str;
	}

};

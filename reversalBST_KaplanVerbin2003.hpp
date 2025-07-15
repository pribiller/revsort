/*******************************************************
 * Original author: Priscila Biller
 * Created: July/2025
 * License: GPL v3
 *  
 * This class implements a Balanced Search Tree (BST)
 * whose nodes are 'arcs', a concept related to how
 * reversals are represented (see Kaplan and Verbin (2003) 
 * for more details).
 * 
 * Each node i stores: 
 * (1) a flag to indicate if the arc (i, i+1) is ``oriented``;
 * (2) a flag to indicate if the arc (i, i+1) is ``unused`` (∈ V, see Tannier et al. 2007);
 * (3) the total number of nodes in the subtree that are ``unused`` (∈ V);
 * (4) the total number of nodes in the subtree that are ``unused`` and ``oriented``;
 * 
 * The BST tree is implemented as a Red-black tree. Thus, in addition to 
 * the properties listed above, the nodes also store a ``color`` attribute 
 * (red or black), and its ``black height``, which is the maximum number 
 * of black nodes on a path from the root node to a leaf node. 
 * Each node also keeps pointers to the nodes holding the minimum and maximum 
 * values found in the subtree rooted by the node.
 * 
 * The following operations are implemented:
 * (1) Split of a tree given a value k; -> O(log n)
 * (2) Join two trees T1 and T2, such that all values of T1 are smaller than T2; -> O(log n)
 * (3) Insert a node; -> O(log n)
 * (4) Remove a node; -> O(log n)
 * 
 * All operations are implemented in such way to keep all flags 
 * and counters updated at all points.
 * 
 * Red-black tree implementation is based on:
 * - Cormen's book: "Introduction to Algorithms", chapter 13 ("Red-black trees").
 * - An implementation of Cormen's pseudocode in C++ is provided in: 
 * https://www.geeksforgeeks.org/cpp/red-black-tree-in-cpp/
 * 
 * For more details:
 * - Eric Tannier et al., 2007: "Advances on sorting by reversals";
 * - Kaplan and Verbin, 2003: "Efficient Data Structures and a 
 * New Randomized Approach for Sorting Signed Permutations by Reversals"
 * 
 *******************************************************/

#pragma once // It avoids class redefinition.

#include <iostream>
#include <list>
#include <vector>
#include <numeric>   // iota
#include <cmath> 
#include <string>	// convert input args (string to int)
#include <utility>   // move
#include <algorithm> // find_if, reverse

#include "genomePermutation_KaplanVerbin2003.hpp"

/*******************************************************
 * Data structures to implement balanced 
 * binary trees for blocks.
 *******************************************************/

// Colors of the nodes in the Red Black Tree (RB).
enum class RbColor {
	RED,   // RED is assigned the value 0
	BLACK  // BLACK is assigned the value 1
};

/* Each node corresponds to an arc v_i. */
template <typename BlockT>
class Node {
public:

	/* General tree properties. */
	Node* parent{nullptr};
	Node* left{nullptr};
	Node* right{nullptr};

	/* Red-black tree properties. */

	// Color: black; red;
	RbColor color{RbColor::RED};  // Initially RED so number of black nodes in the tree is unchanged.
	// maximum number of black nodes on a path from the root node to a leaf node.
	// Here ``black height`` takes into account the own node's color, i.e.,
	// if the node itself is black, then its blackHeight is >= 1.
	int blackHeight{0}; 

	// Keep minimum and maximum values found in the subtree.
	Node* min{nullptr};
	Node* max{nullptr};

	/* Genome properties. */

	const Gene<BlockT>& gene;
	const Gene<BlockT>& gene_next;

	// If reversed=true, then nodes should be ordered backwards 
	// with respect to the original order, and all nodes change orientation.
	bool reversed{false};

	// If unused=true, then the arc represented by the node has not been
	// included yet in the sorting scenario.
	bool unused{true};

	// Total number of nodes in the subtree that have flag unused=true;
	int unused_tot{0};

	// Total number of nodes in the subtree that have flag 
	// unused=true and oriented=true;
	// If the ``cummulated`` reversed flag from the root until this node 
	// is equal to ``true``, then this amount should be read as:
	// unused=true and oriented=false;
	int unused_oriented_tot{0};

	Node(const Gene<BlockT>& g, const Gene<BlockT>& g_next):gene(g),gene_next(g_next),min(this),max(this){ 
		const bool sign_g	   = (g.reversed != g.block->reversed);
		const bool sign_g_next = (gene_next.reversed != gene_next.block->reversed);
	}

	bool isSorted() const {
		const int pos_g = gene.block->genePosAbs(gene);
		const int pos_g_next   = gene_next.block->genePosAbs(gene_next);
		const bool sign_g	   = (gene.reversed != gene.block->reversed);
		const bool sign_g_next = (gene_next.reversed != gene_next.block->reversed);
		return ((sign_g == sign_g_next) && (((!sign_g) && ((pos_g_next-pos_g) == 1)) || ((sign_g) && ((pos_g-pos_g_next) == 1))));
	}

	// ``true`` if signs of genes (i) and (i+1) are flipped in the 
	// current permutation; false otherwise.
	inline bool getOrientation() const {
		const bool sign_g	   = (gene.reversed != gene.block->reversed);
		const bool sign_g_next = (gene_next.reversed != gene_next.block->reversed);
		return (sign_g != sign_g_next);
	}

	// Returns the position of node i's successor in the current permutation.
	inline int getPosNext() const {
		return gene_next.block->genePosAbs(gene_next);
	}

	inline std::string printNode() const {
		return std::to_string(gene.id) + ","  + std::to_string(gene_next.id) + "[" + std::to_string(getPosNext()) +"]";
	}

	inline std::string printNodeDetailed() const {
		return std::to_string(gene.id) + ","  + std::to_string(gene_next.id) + "[" + std::to_string(getPosNext()) +"]" + " (" + ((color == RbColor::RED) ? "RED" : "BLACK") + ";BH=" + std::to_string(blackHeight) + ";UN=" + std::to_string(unused_tot) + ";OR=" + std::to_string(unused_oriented_tot) + ";MIN=" + min->printNode() + ";MAX=" + max->printNode() + ";REV=" + (reversed ? "T" : "F") + ";PARENT=" + ((parent == nullptr)? "-" : std::to_string(parent->gene.id))+ ")";
	}

	// Overload the '<' operator.
	friend bool operator<(const Node &a, const Node &b){
		return a.getPosNext() < b.getPosNext();
	}
	// Overload the '>' operator.
	friend bool operator>(const Node &a, const Node &b){
		return a.getPosNext() > b.getPosNext();
	}
	// Overload the '==' operator.
	friend bool operator==(const Node &a, const Node &b){
		return a.getPosNext() == b.getPosNext();
	}

	// Some properties related to the arc are **not** cleared: 
	// pointers to the genes (i and i+1) and flag ``unused``.
	void cleanNode() {
		parent = nullptr;
		left   = nullptr;
		right  = nullptr;

		min    = this;
		max    = this;

		color  = RbColor::RED;
		blackHeight = 0;

		unused_tot  = 0;
		unused_oriented_tot = 0;

		reversed = false;
	}

	void swapChildren() {
		Node* temp = left;
		left  = right;
		right = temp;
	}

	// Utility function specific to Genomic context:
	// Pushing down the reversed ﬂags before applying Rotation.
	// A rotation on nodes with reversed flags turned off 
	// correctly maintains the structure.
	void clearReversedFlag(){
		if(reversed) {
			// Exchange children.
			// TODO: Can I use std::swap?
			swapChildren();
			// Flip the reversed flag in each of them.
			// As the ``reversed`` flag is also flipped for the parent node, 
			// the current reversed state for the children is not modified.
			if (right != nullptr) {right->reversed = !right->reversed;}
			if (left  != nullptr) {left->reversed  = !left->reversed;}
			// Flipping the sign of the element at the node and update its counts.
			reversed = !reversed;
			unused_oriented_tot = unused_tot - unused_oriented_tot;
			// Flip min and max nodes.
			Node* tmp = min;
			min = max;
			max = tmp;
		}
	}
};

/********************************************************
 * Implementation of Red-black trees, with additional 
 * adaptations for the sorting by reversal problem.
 * 
 * Adaptations are detailed in Kaplan and Verbin (2003).
 * 
 ********************************************************/
template <typename BlockT>
class RedBlackTree {
private:

	// Utility function: Left Rotation
	/* After rotateLeft :
	- Right child of node changes; 
		- In more details:
			- Right child becomes the node's parent;
			- New right child is the left grandchild of the right child.
	- Left child of node remains the same;
	- In a Red-Black tree, every path from a given node to any of its leaf 
		nodes goes through the same number of black nodes.
	- Thus, the rotation does not change the ``black height`` of a node,
		as its left child does not change.
	*/
	void rotateLeft(Node<BlockT>*& node) {
		
		// WARNING - Be careful with edge cases:
		// node has reversed=True, it has a right child but no left child.
		// after clear reversed flag, now node has a left child but NO RIGHT CHILD to do the rotation.

		// Clear reversed flags (turn them off),
		// to make sure that rotation is valid.
		node->clearReversedFlag();
		node->right->clearReversedFlag();
		if (node->right->left != nullptr) {node->right->left->clearReversedFlag();}

		Node<BlockT>* child = node->right;
		node->right = child->left;
		if (node->right != nullptr)
			node->right->parent = node;
		child->parent = node->parent;
		if (node->parent == nullptr)
			root = child;
		else if (node == node->parent->left)
			node->parent->left = child;
		else
			node->parent->right = child;
		child->left = node;
		node->parent = child;

		// Update min and max values.
		node->max  = (node->right != nullptr) ? node->right->max : node;
		child->min = node->min;

		// Update counts of rotated nodes.
		// node, node->right (child), and node->right->left (node's new right) have the same ``reversed`` flag (=false).
		const int unused_node_nochild = node->unused_tot - child->unused_tot;
		const int unused_oriented_node_nochild = node->unused_oriented_tot - child->unused_oriented_tot;
		const int unused_new_right = ((node->right == nullptr) ? 0 : node->right->unused_tot);
		const int unused_oriented_new_right = ((node->right == nullptr) ? 0 : node->right->unused_oriented_tot);
		// The order the update is done is important (bottom-up direction).
		// The new count for node adds the counts of the new right child and remove the counts of the old right child.
		node->unused_tot = unused_node_nochild + unused_new_right;
		node->unused_oriented_tot = unused_oriented_node_nochild + unused_oriented_new_right;
		// The new count for child adds the counts of the new left child and remove the counts of the old left child.
		child->unused_tot = child->unused_tot + node->unused_tot -unused_new_right;
		child->unused_oriented_tot = child->unused_oriented_tot + node->unused_oriented_tot-unused_oriented_new_right;
	}

	// Utility function: Right Rotation
	/* After rotateRight :
	- Left child of node changes; 
		- In more details:
			- Left child becomes the node's parent;
			- New left child is the right grandchild of the left child.
	- Right child of node remains the same;
	- In a Red-Black tree, every path from a given node to any of its leaf 
		nodes goes through the same number of black nodes.
	- Thus, the rotation does not change the ``black height`` of a node,
		as its right child does not change.
	*/
	void rotateRight(Node<BlockT>*& node)
	{
		// Clear reversed flags (turn them off),
		// to make sure that rotation is valid.
		node->clearReversedFlag();
		node->left->clearReversedFlag();
		if (node->left->right != nullptr) {node->left->right->clearReversedFlag();}

		Node<BlockT>* child = node->left;
		node->left = child->right;
		if (node->left != nullptr)
			node->left->parent = node;
		child->parent = node->parent;
		if (node->parent == nullptr)
			root = child;
		else if (node == node->parent->left)
			node->parent->left = child;
		else
			node->parent->right = child;
		child->right = node;
		node->parent = child;

		// Update min and max values.
		node->min  = (node->left != nullptr) ? node->left->min : node;
		child->max = node->max;

		// Update counts of rotated nodes.
		// node, node->left (child), and node->left->right (node's new left) have the same ``reversed`` flag (=false).
		const int unused_node_nochild = node->unused_tot - child->unused_tot;
		const int unused_oriented_node_nochild = node->unused_oriented_tot - child->unused_oriented_tot;
		const int unused_new_left = ((node->left == nullptr) ? 0 : node->left->unused_tot);
		const int unused_oriented_new_left = ((node->left == nullptr) ? 0 : node->left->unused_oriented_tot);
		// The order the update is done is important (bottom-up direction).
		// The new count for node adds the counts of the new left child and remove the counts of the old left child.
		node->unused_tot = unused_node_nochild + unused_new_left;
		node->unused_oriented_tot = unused_oriented_node_nochild + unused_oriented_new_left;
		// The new count for child adds the counts of the new right child and remove the counts of the old right child.
		child->unused_tot = child->unused_tot + node->unused_tot -unused_new_left;
		child->unused_oriented_tot = child->unused_oriented_tot + node->unused_oriented_tot -unused_oriented_new_left;
	}

	// Utility function: Fixing insertion violation.
	void fixInsert(Node<BlockT>* node) {
		Node<BlockT>* parent	  = nullptr;
		Node<BlockT>* grandparent = nullptr;
		while ((node != root) && (node->color == RbColor::RED) && (node->parent->color == RbColor::RED)) {
			parent	  = node->parent;
			grandparent = parent->parent;

			// Clear reversed flags of grandparent, parent, and node.
			// It must be in this order since reversed flags are pushed down in the tree.
			grandparent->clearReversedFlag();
			parent->clearReversedFlag();
			node->clearReversedFlag();

			if (parent == grandparent->left) {
				Node<BlockT>* uncle = grandparent->right;
				if ((uncle != nullptr) && (uncle->color == RbColor::RED)) {
					grandparent->color = RbColor::RED;
					parent->color = RbColor::BLACK;
					uncle->color  = RbColor::BLACK;
					// Update the ``blackHeight`` of parent and uncle nodes.
					parent->blackHeight += 1;
					uncle->blackHeight  += 1;
					// Repeat the ``fixing insertion violation`` process for the grandparent node.
					node = grandparent;
				} else {
					if (node == parent->right) {
						rotateLeft(parent);
						node   = parent;
						parent = node->parent;
					}
					rotateRight(grandparent);
					// At this point, grandparent = BLACK; parent = RED.
					std::swap(parent->color, grandparent->color);
					// Update the black heights of the swaped nodes:
					parent->blackHeight += 1;
					grandparent->blackHeight -= 1;
					node = parent;
				}
			// Symmetrical case (parent == grandparent->right).
			} else {
				Node<BlockT>* uncle = grandparent->left;
				if ((uncle != nullptr) && (uncle->color == RbColor::RED)) {
					grandparent->color = RbColor::RED;
					parent->color = RbColor::BLACK;
					uncle->color  = RbColor::BLACK;
					// Update the ``blackHeight`` of parent and uncle nodes.
					parent->blackHeight += 1;
					uncle->blackHeight  += 1;
					// Repeat the ``fixing insertion violation`` process for the grandparent node.
					node = grandparent;
				} else {
					if (node == parent->left) {
						rotateRight(parent);
						node = parent;
						parent = node->parent;
					}
					rotateLeft(grandparent);
					// At this point, grandparent = BLACK; parent = RED.
					std::swap(parent->color, grandparent->color);
					// Update the black heights of the swaped nodes:
					parent->blackHeight += 1;
					grandparent->blackHeight -= 1;
					node = parent;
				}
			}
		}
		// Check if the root is BLACK. Root must always be black (following Cormen's defition).
		if(root->color == RbColor::RED){
			root->color = RbColor::BLACK;
			root->blackHeight += 1;
		}
	}

	// Utility function: Fixing Deletion Violation
	void fixDelete(Node<BlockT>* node, Node<BlockT>* parent) {
		// Edge case: node is the new root.
		if((node != nullptr) && (node == root)){
			if (node->color == RbColor::RED){
				node->blackHeight = 0; // It will be incremented later in the end of the function.
			} else {
				node->blackHeight = 1; // It will not be incremented later in the end of the function.
			}
		}
		// Normal case: node is not a root.
		while ((node != root) && ((node == nullptr) || (node->color == RbColor::BLACK))) {
			// Clear reversed flags of parent, node, and its sibling.
			if(parent != nullptr){
				parent->clearReversedFlag();
				if(parent->left  != nullptr){parent->left->clearReversedFlag();}
				if(parent->right != nullptr){parent->right->clearReversedFlag();}
			} else if (node != nullptr) {
				node->clearReversedFlag();
			}
			if (node == parent->left) {
				Node<BlockT>* sibling = parent->right;
				if (sibling->color == RbColor::RED) {
					sibling->color = RbColor::BLACK;
					parent->color  = RbColor::RED;
					// Update the ``blackHeight`` of parent and sibling nodes.
					sibling->blackHeight += 1;
					parent->blackHeight -= 1;
					rotateLeft(parent);
					sibling = parent->right;
				}
				// Sibling can be either BLACK or RED (both its children are BLACK).
				if (((sibling->left  == nullptr) || (sibling->left->color  == RbColor::BLACK))
				&&  ((sibling->right == nullptr) || (sibling->right->color == RbColor::BLACK))) {
					// Update ``blackHeight`` if needed.
					if (sibling->color == RbColor::BLACK) {
						sibling->blackHeight -= 1;
						parent->blackHeight  -= 1;
					}
					sibling->color = RbColor::RED;
					node   = parent;
					parent = node->parent;
				// Sibling is BLACK (at least one of its children is RED).
				} else {
					// Right is BLACK. Left is RED.
					if ((sibling->right == nullptr) || (sibling->right->color == RbColor::BLACK)) {
						if (sibling->left != nullptr){
							sibling->left->blackHeight += 1;
							sibling->left->color = RbColor::BLACK;
						}
						sibling->blackHeight -= 1;
						sibling->color = RbColor::RED;
						rotateRight(sibling);
						sibling = parent->right;
					}
					// Update ``blackHeight`` of parent and sibling.
					// In this case, parent is going to become the left child of ``sibling`` (via rotateLeft).
					// Sibling color at the end is going to be parent's color.
					// Parent color at the end is going to be BLACK.
					if(parent->color == RbColor::BLACK){ // Parent is BLACK.
						// Sibling is going to become parent of its black parent node. 
						// Thus, black height of sibling increases by 1.
						sibling->blackHeight += 1;
						if(sibling->color != parent->color){ // Sibling: RED -> BLACK. Parent does not lose any black child.
							sibling->blackHeight += 1;
						} else { // Sibling is BLACK. Parent loses one black height.
							parent->blackHeight -= 1;
						}
					}
					sibling->color = parent->color;
					parent->color  = RbColor::BLACK;
					if (sibling->right != nullptr){
						if(sibling->right->color == RbColor::RED){sibling->right->blackHeight += 1;}
						sibling->right->color = RbColor::BLACK;
					}
					rotateLeft(parent);
					node = root;
				}
			} else {
				Node<BlockT>* sibling = parent->left;
				if (sibling->color == RbColor::RED) {
					sibling->color = RbColor::BLACK;
					parent->color  = RbColor::RED;
					// Update the ``blackHeight`` of parent and sibling nodes.
					sibling->blackHeight += 1;
					parent->blackHeight  -= 1;
					rotateRight(parent);
					sibling = parent->left;
				}
				// Sibling can be either BLACK or RED (both its children are BLACK).
				if ( ((sibling->left  == nullptr) || (sibling->left->color  == RbColor::BLACK))
				  && ((sibling->right == nullptr) || (sibling->right->color == RbColor::BLACK))) {
					// Update ``blackHeight`` if needed.
					if (sibling->color == RbColor::BLACK) {
						sibling->blackHeight -= 1;
						parent->blackHeight  -= 1;
					}
					sibling->color = RbColor::RED;
					node   = parent;
					parent = node->parent;
				// Sibling is BLACK (at least one of its children is RED).
				} else {
					// Left is BLACK. Right is RED.
					if ((sibling->left == nullptr) || (sibling->left->color == RbColor::BLACK)) {
						if (sibling->right != nullptr){
							sibling->right->blackHeight += 1;
							sibling->right->color = RbColor::BLACK;
						}
						sibling->blackHeight -= 1;
						sibling->color = RbColor::RED;
						rotateLeft(sibling);
						sibling = parent->left;
					}
					// Update ``blackHeight`` of parent and sibling.
					// In this case, parent is going to become the right child of ``sibling`` (via rotateRight).
					// Sibling color at the end is going to be parent's color.
					// Parent color at the end is going to be BLACK.
					if(parent->color == RbColor::BLACK){ // Parent is BLACK.
						// Sibling is going to become parent of its black parent node. 
						// Thus, black height of sibling increases by 1.
						sibling->blackHeight += 1;
						if(sibling->color != parent->color){ // Sibling: RED -> BLACK. Parent does not lose any black child.
							sibling->blackHeight += 1;
						} else { // Sibling is BLACK. Parent loses one black height.
							parent->blackHeight -= 1;
						}
					}
					sibling->color = parent->color;
					parent->color  = RbColor::BLACK;
					if (sibling->left != nullptr){
						if(sibling->left->color == RbColor::RED){sibling->left->blackHeight += 1;}
						sibling->left->color = RbColor::BLACK;
					}
					rotateRight(parent);
					node = root;
				}
			}
		}
		if(node != nullptr){
			if (node->color == RbColor::RED) {node->blackHeight += 1;}
			node->color = RbColor::BLACK;
		}
	}

	// Utility function: Transplant nodes in Red-Black Tree
	void transplant(Node<BlockT>*& u, Node<BlockT>*& v){

		// Clear ``reversed`` flags before transplanting.
		if (u->parent != nullptr){
			u->parent->clearReversedFlag();
		}
		u->clearReversedFlag();

		// TODO: Update any counts? (blackHeight, total unused, etc.?)

		// Transplant: ``v`` takes the place of ``u``.
		if (u->parent == nullptr){
			root = v;
		} else if (u == u->parent->left) {
			u->parent->left = v;
		} else {
			u->parent->right = v;
		}
		if (v != nullptr)
			v->parent = u->parent;
	}

	// Utility function: Helper to print Red-Black Tree
	void printHelper(Node<BlockT>* root, std::string indent, bool last) const {
		if (root != nullptr) {
			std::cout << indent;
			if (last) {
				std::cout << "R----";
				indent += "   ";
			}
			else {
				std::cout << "L----";
				indent += "|  ";
			}
			std::cout << root->printNodeDetailed() << std::endl;
			printHelper(root->left, indent, false);
			printHelper(root->right, indent, true);
		}
	}

	// Case: t1 < k < t2; t1.blackHeight >= t2.blackHeight.
	// The tree with smaller black height will always end up 
	// as a subtree of the tree with bigger black height.
	void joinRight(RedBlackTree& t1, Node<BlockT>*& node, RedBlackTree& t2){
		// Clear reversed flags (turn them off),
		// in order to simplify update of counts.
		t1.root->clearReversedFlag();
		node->clearReversedFlag();
		t2.root->clearReversedFlag();

		bool const reversed_block = node->gene.block->reversed;
		node->unused_tot = (node->unused) ? 1 : 0;
		node->unused_oriented_tot = ((node->unused) && (node->getOrientation() != reversed_block)) ? 1 : 0;

		// Update the counts of all nodes from the root until insertion position.
		const int unused_tot_upd = node->unused_tot + t2.root->unused_tot;
		const int unused_oriented_tot_upd = node->unused_oriented_tot + t2.root->unused_oriented_tot;
		Node<BlockT>* current = t1.root;
		while((current->blackHeight != t2.root->blackHeight) || (current->color != RbColor::BLACK)){
			// std::cout << "[Join trees] current="<< current->printNodeDetailed() << std::endl;
			// Update counts.
			current->unused_tot += unused_tot_upd;
			current->unused_oriented_tot += unused_oriented_tot_upd;
			current->max = t2.root->max;
			// Move to the next node in the path.
			current = current->right;
			current->clearReversedFlag();
		}
		// At this point: t1 and t2 have the same height; t1 and t2 roots are black.
		// std::cout << "[Join trees] Same height="<< current->printNodeDetailed() << std::endl;

		// Insert node.
		Node<BlockT>* parent = current->parent;
		if (parent == nullptr) {
			t1.root = node;
		} else {
			parent->right = node;
		}
		node->parent = parent;  // red (potential conflict with parent; fixed later in fixInsert)
		node->left   = current; // black
		node->right  = t2.root; // black
		node->left->parent  = node;
		node->right->parent = node;
		node->blackHeight = current->blackHeight; // New node is always (initially) a RED node.
		node->unused_tot += (node->left->unused_tot + node->right->unused_tot);
		// node, t1 and t2 have the same ``reversed`` flag (= false).
		node->unused_oriented_tot += (node->left->unused_oriented_tot + node->right->unused_oriented_tot);
		node->min    = node->left->min;
		node->max    = node->right->max;
		// Fix tree balance if needed.
		t1.fixInsert(node);
	}

	// Case: t1 < k < t2; t1.blackHeight < t2.blackHeight.
	// The tree with smaller black height will always end up 
	// as a subtree of the tree with bigger black height.
	void joinLeft(RedBlackTree& t1, Node<BlockT>*& node, RedBlackTree& t2){
		// Clear reversed flags (turn them off),
		// in order to simplify update of counts.
		t2.root->clearReversedFlag();
		node->clearReversedFlag();
		t1.root->clearReversedFlag();

		bool const reversed_block = node->gene.block->reversed;
		node->unused_tot = (node->unused) ? 1 : 0;
		node->unused_oriented_tot = ((node->unused) && (node->getOrientation() != reversed_block)) ? 1 : 0;

		// Update in the counts of all nodes in the way from the root until insertion position.
		const int unused_tot_upd = node->unused_tot + t1.root->unused_tot;
		const int unused_oriented_tot_upd = node->unused_oriented_tot + t1.root->unused_oriented_tot;
		Node<BlockT>* current = t2.root;
		while((current->blackHeight != t1.root->blackHeight) || (current->color != RbColor::BLACK)){
			// Update counts.
			current->unused_tot += unused_tot_upd;
			current->unused_oriented_tot += unused_oriented_tot_upd;
			current->min = t1.root->min;
			current = current->left;
			current->clearReversedFlag();
		}
		// At this point: t1 and t2 have the same height; t1 and t2 roots are black.

		// Insert node.
		Node<BlockT>* parent = current->parent;
		if (parent == nullptr) {
			t2.root = node;
		} else {
			parent->left = node;
		}
		node->parent = parent;   // red (potential conflict with parent; fixed later in fixInsert)
		node->left   = t1.root;  // black
		node->right  = current;  // black
		node->left->parent  = node;
		node->right->parent = node;
		node->blackHeight = current->blackHeight; // New node is always (initially) a RED node.
		node->unused_tot += (node->left->unused_tot + node->right->unused_tot);
		// node, t1 and t2 have the same ``reversed`` flag (= false).
		node->unused_oriented_tot += (node->left->unused_oriented_tot + node->right->unused_oriented_tot);
		node->min    = node->left->min;
		node->max    = node->right->max;
		// Fix tree balance if needed. It also updates the counts of nodes.
		t2.fixInsert(node);
	}

public:
	Node<BlockT>* root{nullptr}; // Root of the Red-Black Tree

	RedBlackTree() {

	}

	void cleanTree(){
		root = nullptr;
	}

	inline Node<BlockT>* getGlobalMin() const {
		return (root->reversed) ? root->max : root->min;
	}

	inline Node<BlockT>* getGlobalMax() const {
		return (root->reversed) ? root->min : root->max;
	}

	void reroot(Node<BlockT>* new_root){
		root = new_root;
		if(root != nullptr){
			if(root->parent != nullptr){
				// WARNING! If new_root was part of some previous tree,
				// this other tree will become invalid (unbalanced, outdated, etc.).
				if (root->parent->right == root){
					root->parent->right = nullptr;
				} else if (root->parent->left == root){
					root->parent->left = nullptr;
				}
				root->parent = nullptr;
			}
			if (root->color == RbColor::RED) {
				root->color = RbColor::BLACK;
				root->blackHeight += 1;
			}
			root->clearReversedFlag();
		}
	}

	// Join acts on two red-black trees t1 and t2 and a key k, 
	// where t1 < k < t2, i.e. all keys in t1 are less than k (node), and 
	// all keys in t2 are greater than k. 
	// WARNING: The data structure of both trees will be modified.
	//			Specifically, node pointers will be modified by this function.
	void join(Node<BlockT>* node, RedBlackTree& another_tree){
		if (root == nullptr) {
			// std::cout << "[Join trees] Current tree is empty - middle node: " << node->printNodeDetailed() << " " << std::endl;
			another_tree.insert(node);
			root = another_tree.root;
			// std::cout << "[Join trees] Current tree is empty - result" << std::endl;
			// printTree();
		} else if (another_tree.root == nullptr) {
			insert(node);
		} else {
			node->cleanNode();
			// std::cout << "[Join trees]  Middle (k) =" << node->gene.id << " "<< std::endl;
			// this tree < k < another tree.
			if((*root) < (*node)){
				// std::cout << "[Join trees] this tree < k < another tree"<< std::endl;
				// std::cout << "[Join trees] this tree"<< std::endl;
				// printTree();
				// std::cout << "[Join trees] another tree"<< std::endl;
				// another_tree.printTree();

				// Add nodes from t2 to t1.
				if (root->blackHeight > another_tree.root->blackHeight) {
					joinRight(*this, node, another_tree);
				// Add nodes from t1 to t2.
				} else {
					joinLeft(*this, node, another_tree);
					root = another_tree.root;
				}
			// another tree < k < this tree
			} else {
				// std::cout << "[Join trees] another tree < k < this tree"<< std::endl;
				// std::cout << "[Join trees] this tree"<< std::endl;
				// printTree();
				// std::cout << "[Join trees] another tree"<< std::endl;
				// another_tree.printTree();

				// Add nodes from t1 to t2.
				if (root->blackHeight > another_tree.root->blackHeight) {
					joinLeft(another_tree, node, *this);
				// Add nodes from t2 to t1.
				} else {
					joinRight(another_tree, node, *this);
					root = another_tree.root;
				}
			}
		}
	}

	// Split the current tree into two subtrees: 
	// current tree (with all elements <= key); 
	// another_tree (with all elements > key).
	void split(int key, RedBlackTree& another_tree) {
		if(root == nullptr){return;} // Nothing to split.
		// Check extreme cases.
		int const min_key = (root->reversed) ? root->max->getPosNext() : root->min->getPosNext();
		int const max_key = (root->reversed) ? root->min->getPosNext() : root->max->getPosNext();
		if(min_key > key){
			// std::cout << "[split] Min is bigger than key : min_key:" << min_key << " < key:" << key << std::endl;
			another_tree.reroot(root);
			cleanTree();
			return;
		} else if(max_key <= key){
			// std::cout << "[split] Max is smaller than key : min_key" << max_key << " <= key:" << key << std::endl;
			return;
		}
		// std::cout << "[split] No extreme cases" << std::endl;
		Node<BlockT>* current = root;
		Node<BlockT>* parent  = nullptr;
		RedBlackTree tmp_t{};
		// Clean both trees.
		cleanTree();
		another_tree.cleanTree();
		// Find position to insert node.
		while (current != nullptr){
			// Make sure that the ``reversed`` flag is set to false 
			// in all the way between the root and the insertion point.
			// std::cout << "[split] Current node: " << current->printNodeDetailed() << std::endl;
			// std::cout << "[split] key=" << key << " current=" << current->getPosNext() << std::endl;
			current->clearReversedFlag();
			if(key < current->getPosNext()) {
				parent  = current;
				current = current->left;
				// All nodes in the right subtree (+ the current node) 
				// are bigger than key, but smaller than all values 
				// inserted in another_tree so far.
				tmp_t.reroot(parent->right);
				if(tmp_t.root == nullptr) {
					another_tree.insert(parent);
				} else {
					// Parent is the lowest value in ``tmp``.
					tmp_t.insert(parent);
					// Extract the maximum node from the right tree.
					// Both trees have ``reversed`` flag false.
					Node<BlockT>* k = tmp_t.root->max;
					tmp_t.remove(k);
					another_tree.join(k, tmp_t);
				}
				// std::cout << "[split] [key<current] After joining" << std::endl;
				// another_tree.printTree();

			} else if(key > current->getPosNext()) {
				parent  = current;
				current = current->right;

				// All nodes in the left subtree (+ the current node) 
				// are smaller than key, but bigger than all values 
				// inserted in another_tree so far.
				tmp_t.reroot(parent->left);
				if(tmp_t.root == nullptr) {
					insert(parent);
				} else {
					// Parent is the biggest value in ``tmp``.
					tmp_t.insert(parent);
					// Extract the minimum node from the left tree.
					// Both trees have ``reversed`` flag false.
					Node<BlockT>* k = tmp_t.root->min;
					tmp_t.remove(k);
					join(k, tmp_t);
				}
				// std::cout << "[split] [key>current] After joining" << std::endl;
				// printTree();

			// key == current.
			// ``Shortcut case``: loop does not need to go until a leaf.
			} else {
				parent  = current;
				current = current->right;

				// All nodes in the left subtree are smaller than key,
				// but bigger than all values inserted in ``tree`` so far.
				tmp_t.reroot(parent->left);
				if(tmp_t.root == nullptr) {
					insert(parent);
				} else {
					// Parent is the biggest value in the ``tmp``.
					tmp_t.insert(parent);
					// Extract the minimum node from the left tree.
					// Both trees have ``reversed`` flag false.
					Node<BlockT>* k = tmp_t.root->min;
					tmp_t.remove(k);
					join(k, tmp_t);
				}				
				// std::cout << "[split] [key=current] After joining" << std::endl;
				// printTree();
				// All nodes in the right subtree are bigger than key,
				// but smaller than all values inserted in another_tree so far.
				if(current != nullptr){
					tmp_t.reroot(current);
					// Extract the minimum or maximum node from the right tree.
					// Both trees have ``reversed`` flag false.
					Node<BlockT>* k = ((another_tree.root == nullptr) || (*(another_tree.root->min)) > (*(tmp_t.root->max))) ? tmp_t.root->max : tmp_t.root->min;
					tmp_t.remove(k);
					another_tree.join(k, tmp_t);
				}
				break;

				/*************************************/
				// Alternative implementation of this case.
				// Here, ``key`` value is adjusted to force
				// the loop go until one of the leaves.
				// WARNING: key needs to be declared as a double.

				// parent  = current;
				// current = current->right;
				
				// reroot(parent->left);
				// insert(parent);

				// // ``key`` value is adjusted to force
				// // the loop go until one of the leaves.
				// key += 0.1;
			}
		}
	}

	// Public function: Remove a value from the Red-Black Tree.
	// ``node`` is the node to be removed.
	void remove(Node<BlockT>* node) {
		Node<BlockT>* z = node;
		Node<BlockT>* x = nullptr;
		Node<BlockT>* y = nullptr;
		if (z == nullptr) {return;}
		// Make sure the ``reversed`` flags are off so comparison is valid.
		Node<BlockT>* current = root;
		while(current != node){
			current->clearReversedFlag();
			// Traverse tree in the correct order (reversed flag=false).
			if ((*node) < (*current)){
				current = current->left;
			} else {
				current = current->right;
			}
		}
		node->clearReversedFlag();
		// Check if the deleted node is the min or max of its parent.
		// In this scenario, the node has at most one child.
		Node<BlockT>* new_min = nullptr;
		Node<BlockT>* new_max = nullptr;
		if((z->parent != nullptr) && ((z == z->parent->min) || (z == z->parent->max))){
			Node<BlockT>* child = (z->left != nullptr) ? z->left : z->right;
			// Make sure the ``reversed`` flags are off so comparison is valid.
			if (child != nullptr) {child->clearReversedFlag();}
			// Finds new min/max.
			if(z == z->parent->min){
				if ((child != nullptr) && ((*(child->min)) < (*(z->parent)))){
					new_min = child->min;
				} else {
					new_min = z->parent;
				}
			}
			if(z == z->parent->max){
				if ((child != nullptr) && ((*(child->max)) > (*(z->parent)))){
					new_max = child->max;
				} else {
					new_max = z->parent;
				}
			}
		}
		// Update min/max and other counts in the path between the root and the deleted node.
		// Make sure that all ``reversed`` flags in the path 
		// between the root and the deleted node are ``off``.
		current = root;
		bool const reversed_block = node->gene.block->reversed;
		const int unused_upd  = ((node->unused) ? 1 : 0);
		const int unused_oriented_upd = ((node->unused) && (node->getOrientation() != reversed_block)) ? 1 : 0;
		while(current != node){
			if(current->max == node) {current->max = new_max;}
			if(current->min == node) {current->min = new_min;}
			current->unused_tot -= unused_upd;
			current->unused_oriented_tot -= unused_oriented_upd;
			// Traverse tree in the correct order (reversed flag=false).
			if ((*node) < (*current)){
				current = current->left;
			} else {
				current = current->right;
			}
		}
		y = z;
		RbColor yOriginalColor = y->color;
		if (z->left == nullptr) {
			x = z->right;
			transplant(z, z->right);
		} else if (z->right == nullptr) {
			x = z->left;
			transplant(z, z->left);
		} else {
			y = z->right->min;
			yOriginalColor = y->color;
			x = y->right;
			if (y->parent == z) {
				if (x != nullptr) {x->parent = y;}
			} else {
				transplant(y, y->right);
				y->right = z->right;
				y->right->parent = y;
			}
			transplant(z, y);
			y->left = z->left;
			y->left->parent = y;
			y->color = z->color;
		}
		if (yOriginalColor == RbColor::BLACK) {
			// Check if tree became empty.
			if ((x == root) && (x == nullptr)){
				cleanTree();
			} else {
				Node<BlockT>* parent = (x == nullptr) ? z->parent : x->parent;
				fixDelete(x, parent);
			}
		}
	}

	// Public function: Insert a node into the Red-Black Tree.
	// Important properties that need to be updated: 
	// pointers to the genes (i and i+1); flag ``unused``.
	void insert(Node<BlockT>* node) {

		// Re-set properties related to the tree.
		node->cleanNode();
		Node<BlockT>* current = root;
		Node<BlockT>* parent  = nullptr;

		// ``reversed`` flag of the block should be taken 
		// into account in the orientation of a node 
		// (but **not** on how the tree is traversed).
		bool const reversed_block = node->gene.block->reversed;
		// Find position to include node in the tree.
		// Position is always as a leaf (parent->left/right=NULL).
		// As it goes through the tree, update the counts of nodes.
		while (current != nullptr) {

			// Make sure that all ``reversed`` flags from root 
			// until the inserted position are turned off.
			current->clearReversedFlag();

			parent = current;
			// Update counts of the current node.
			if(node->unused){
				current->unused_tot += 1;
				// Update number of unused oriented arcs depending on 
				// the current state of flag ``reversed``.
				// If reversed = false, ``unused_oriented_tot`` keeps all 
				// arcs that are unused and oriented.
				// If reversed = true, ``unused_oriented_tot`` keeps all 
				// arcs that are unused and **not** oriented.
				if(node->getOrientation() != reversed_block) {current->unused_oriented_tot += 1;} 
			}
			// Traverse the tree in the correct order (Left;Root;Right).
			// Cumulated ``reversed`` flag from root until current node is always false.
			if ((*node) < (*current)){
				if((*node) < (*(current->min))){current->min = node;}
				current = current->left;
			} else {
				if((*node) > (*(current->max))){current->max = node;}
				current = current->right;
			}
		}
		// Insert node.
		node->parent	  = parent;
		node->blackHeight = 0; // New node is always (initially) a RED node.
		node->unused_tot  = (node->unused) ? 1 : 0;
		// It takes into account the ``reversed`` flag of the block.
		node->unused_oriented_tot = ((node->unused) && (node->getOrientation() != reversed_block)) ? 1 : 0;

		const bool sign_g	  = (node->gene.reversed != node->gene.block->reversed);
		const bool sign_g_next = (node->gene_next.reversed != node->gene_next.block->reversed);
		
		// std::cout << "[Insert] Block g   = " << node->gene.block->printBlock() << std::endl;
		// std::cout << "[Insert] Block g+1 = " << node->gene_next.block->printBlock() << std::endl;
		// std::cout << "[Insert] unused=" << node->unused_tot << " oriented=" << node->unused_oriented_tot << " " << node->getOrientation() << " / " << reversed_block << " / " << sign_g << " / " << sign_g_next << std::endl;

		if (parent == nullptr) {
			root = node;
		} else if ((*node) > (*parent)) {
			parent->right = node;
		} else {
			parent->left  = node;
		}
		// Fix tree balance if needed. It also updates the counts of nodes.
		fixInsert(node);
	}

	// Returns an unused oriented arc from the tree (if it exists).
	Node<BlockT>* getUnusedOriented(){
		if(root == nullptr){return nullptr;}
		Node<BlockT>* current = root;
		current->clearReversedFlag();
		// ``reversed`` flag of the block should be taken 
		// into account in the orientation of a node 
		// (but **not** on how the tree is traversed).
		bool const reversed_block = root->gene.block->reversed;
		// Find an oriented unused arc in the tree.
		int unused_oriented_tot  = reversed_block ? (current->unused_tot-current->unused_oriented_tot) : current->unused_oriented_tot;
		bool curIsUnusedOriented = ((current->unused) && (current->getOrientation()));
		while ((current != nullptr) && (!curIsUnusedOriented) && (unused_oriented_tot > 0)) {
			// Compute amount of unused arcs in the right and left child.
			int unused_oriented_left = 0;			
			if(current->left != nullptr){
				current->left->clearReversedFlag();
				unused_oriented_left = reversed_block ? (current->left->unused_tot - current->left->unused_oriented_tot) : current->left->unused_oriented_tot;
			}
			int unused_oriented_right = 0;			
			if(current->right != nullptr){
				current->right->clearReversedFlag();
				unused_oriented_right = reversed_block ? (current->right->unused_tot - current->right->unused_oriented_tot) : current->right->unused_oriented_tot;
			}
			// Move to the child with the highest amount of unused oriented arcs.
			if(unused_oriented_left > unused_oriented_right){
				current = current->left;
			} else {
				current = current->right;
			}
			// Update flag and count.
			if(current != nullptr){
				unused_oriented_tot = reversed_block ? (current->unused_tot-current->unused_oriented_tot) : current->unused_oriented_tot;
			} else {
				unused_oriented_tot = 0;
			}
			curIsUnusedOriented = ((current != nullptr) && (current->unused) && (current->getOrientation()));
		}
		// Check if current node represents an arc that is unused and oriented.
		return curIsUnusedOriented ? current : nullptr;
	}

	// Keep node in the tree, but update its status to ``used``.
	void setUsed(Node<BlockT>* node){

		if ((node == nullptr) || (!node->unused)){return;}

		// ``reversed`` flag of the block should be taken 
		// into account in the orientation of a node 
		// (but **not** on how the tree is traversed).
		bool const reversed_block = root->gene.block->reversed;
		const int unused_upd  = 1;
		const int unused_oriented_upd = reversed_block ? 1 : 0; // The current status of the node is ``unoriented``. Thus, if the block is reversed, the node is being counted as ``oriented``.

		// Update node status.
		node->unused = false;
		node->unused_tot -= unused_upd;
		node->unused_oriented_tot -= ((node->unused_tot == 0) ? node->unused_oriented_tot : unused_oriented_upd);

		// Update counts in the path between the root and the node (node's ancestors).
		// Make sure that all ``reversed`` flags in the path 
		// between the root and the node are ``off``.
		Node<BlockT>* current = root;
		while(current != node){
			current->clearReversedFlag();

			current->unused_tot -= unused_upd;
			current->unused_oriented_tot -= unused_oriented_upd;

			// Traverse tree in the correct order (reversed flag=false).
			if ((*node) < (*current)){
				current = current->left;
			} else {
				current = current->right;
			}
		}
	}

	void printTree() const {
		if (root == nullptr)
			std::cout << "Tree is empty." << std::endl;
		else {
			std::cout << "Red-Black Tree:" << std::endl;
			printHelper(root, "", true);
		}
	}
};

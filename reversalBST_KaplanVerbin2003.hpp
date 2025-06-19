/*******************************************************
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
 * 
 * The following operations are implemented:
 * (1) Split of a tree given a value k; -> O(log n)
 * (2) Join two trees T1 and T2, such that all values of T1 are smaller than T2; -> O(log n)
 * (3) Insert a node; -> O(log n)
 * 
 * All operations are implemented in such way to keep all flags 
 * and counters updated at all points.
 * 
 * Red-black tree implementation is based on:
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

#include "genomePermutation.hpp"

/*******************************************************
 * Data structures to implement balanced 
 * binary trees for blocks.
 *******************************************************/

enum Color {
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
	Color color{RED};  // Initially RED so number of black nodes in the tree is unchanged.
	// maximum number of black nodes on a path from the root node to a leaf node.
	// Here ``black height`` takes into account the own node's color, i.e.,
	// if the node itself is black, then its blackHeight is >= 1.
	int blackHeight{0}; 

	/* Genome properties. */

	const Gene<BlockT>& gene;
	const Gene<BlockT>& gene_next;

	// ``true`` if signs of genes (i) and (i+1) are flipped in the 
	// current permutation; false otherwise.	
	bool oriented{false}; 

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

	Node(const Gene<BlockT>& g, const Gene<BlockT>& g_next):gene(g),gene_next(g_next){ 
		const bool sign_g	  = (g.reversed != g.block->reversed);
		const bool sign_g_next = (gene_next.reversed != gene_next.block->reversed);
		oriented = getOrientation();
	}

	bool getOrientation(){
		const bool sign_g	  = (gene.reversed != gene.block->reversed);
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

	// Properties related to the arc itself are **not** touched by this function: 
	// pointers to the genes (i and i+1) and flag ``unused``.
	void cleanNode(Node<BlockT>*& node) {
		node->parent = nullptr;
		node->left   = nullptr;
		node->right  = nullptr;

		node->color		 = RED;
		node->blackHeight = 0;

		node->unused_tot = 0;
		node->unused_oriented_tot = 0;

		node->reversed = false;
		node->oriented = node->getOrientation();
	}

	void swapChildren(Node<BlockT>*& node) {
		if (node != nullptr) {
			Node<BlockT>* temp = node->left;
			node->left  = node->right;
			node->right = temp;
		}
	}

	// Utility function specific to Genomic context:
	// Pushing down the reversed ﬂags before applying Rotation.
	// A rotation on nodes with reversed flags turned off 
	// correctly maintains the structure.
	void clearReversedFlag(Node<BlockT>*& node){
		if((node != nullptr) && (node->reversed)){
			// Exchange children.
			// TODO: Can I use std::swap?
			swapChildren(node);
			// Flip the reversed flag in each of them.
			// As the ``reversed`` flag is also flipped for the parent node, 
			// the current reversed state for the children is not modified.
			if (node->right != nullptr) {node->right->reversed = !node->right->reversed;}
			if (node->left  != nullptr) {node->left->reversed  = !node->left->reversed;}
			// Flipping the sign of the element at the node and update its counts.
			node->reversed = !node->reversed;
			node->unused_oriented_tot = node->unused_tot - node->unused_oriented_tot;
		}
	}

	// Return the ``cummulated`` reversed flag from the root until this node.
	bool getReversedFlag(Node<BlockT>*& node){
		bool reversed = (node == root) ? root->reversed : false;
		while (node != root) {
			reversed = (reversed == node->reversed);
			node = node->parent;
		}
		return reversed;
	}

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
	void rotateLeft(Node<BlockT>*& node)
	{
		// Clear reversed flags (turn them off),
		// to make sure that rotation is valid.
		clearReversedFlag(node);
		clearReversedFlag(node->right);
		if (node->right->left != nullptr) {clearReversedFlag(node->right->left);}

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
		clearReversedFlag(node);
		clearReversedFlag(node->left);
		if (node->left->right != nullptr) {clearReversedFlag(node->left->right);}

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
		while ((node != root) && (node->color == RED) && (node->parent->color == RED)) {
			parent	  = node->parent;
			grandparent = parent->parent;
			if (parent == grandparent->left) {
				Node<BlockT>* uncle = grandparent->right;
				if ((uncle != nullptr) && (uncle->color == RED)) {
					grandparent->color = RED;
					parent->color = BLACK;
					uncle->color  = BLACK;
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
				if ((uncle != nullptr) && (uncle->color == RED)) {
					grandparent->color = RED;
					parent->color = BLACK;
					uncle->color  = BLACK;
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
		if(root->color == RED){
			root->color = BLACK;
			root->blackHeight += 1;
		}
	}

	// Utility function: Helper to print Red-Black Tree
	void printHelper(Node<BlockT>* root, std::string indent, bool last) {
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
			std::string sColor = (root->color == RED) ? "RED" : "BLACK";
			std::cout << root->printNode() << "(" << sColor << ";BH=" << root->blackHeight << ";UN=" << root->unused_tot << ";OR=" << root->unused_oriented_tot << ")" << std::endl;
			printHelper(root->left, indent, false);
			printHelper(root->right, indent, true);
		}
	}

public:
	Node<BlockT>* root{nullptr}; // Root of the Red-Black Tree

	RedBlackTree() {

	}

	// Public function: Insert a node into the Red-Black Tree.
	// Important properties that need to be updated: 
	// pointers to the genes (i and i+1); flag ``unused``.
	void insert(Node<BlockT>* node) {

		// Re-set properties related to the tree.
		cleanNode(node);

		Node<BlockT>* current = root;
		Node<BlockT>* parent  = nullptr;
		bool reversed = false;

		// Find position to include node in the tree.
		// Position is always as a leaf (parent->left/right=NULL).
		// As it goes through the tree, update the counts of nodes.
		while (current != nullptr) {
			parent = current;
			// Update global/accumulated ``reversed`` flag. 
			reversed = (reversed != current->reversed);
			// Update counts of the current node.
			if(node->unused){
				current->unused_tot += 1;
				// Update number of unused oriented arcs depending on 
				// the current state of flag ``reversed``.
				// If reversed = false, ``unused_oriented_tot`` keeps all 
				// arcs that are unused and oriented.
				// If reversed = true, ``unused_oriented_tot`` keeps all 
				// arcs that are unused and **not** oriented.
				if(node->oriented != reversed) {current->unused_oriented_tot += 1;} 
			}
			// Traverse the tree in an inverted order (Right;Root;Left).
			if(reversed){
				if ((*node) < (*current)){
					current = current->right;
				} else {
					current = current->left;
				}				
			// Traverse the tree in the correct order (Left;Root;Right).
			} else{
				if ((*node) < (*current)){
					current = current->left;
				} else {
					current = current->right;
				}
			}
		}

		// Insert node.
		node->parent	  = parent;
		node->blackHeight = 0; // New node is always (initially) a RED node.
		node->unused_tot  = (node->unused) ? 1 : 0;
		// It takes into account current state of ``reversed`` flag.
		node->unused_oriented_tot = ((node->unused) && (node->oriented != reversed)) ? 1 : 0;

		if (parent == nullptr) {
			root = node;
		} else if (((*node) < (*parent)) == reversed) {
			parent->right = node;
		} else {
			parent->left  = node;
		}

		// Fix tree balance if needed. It also updates the counts of nodes.
		fixInsert(node);
	}

	void printTree() {
		if (root == nullptr)
			std::cout << "Tree is empty." << std::endl;
		else {
			std::cout << "Red-Black Tree:" << std::endl;
			printHelper(root, "", true);
		}
	}
};

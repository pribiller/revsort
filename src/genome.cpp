/*******************************************************
 * Original author: Priscila Biller
 * Created: July/2025
 * License: GPL v3
 * 
 * This class provides some basic data structures to convert 
 * a multichromosomal (or unichromosomal) genome with generic 
 * labels to the internal representation used by the methods  
 * to sort genomes.
 * 
 * In the internal representation, a multichromosomal genome 
 * is represented as a list of permutations, with one 
 * permutation per chromosome.
 * 
 * The permutation of a chromosome is represented by a 
 * vector of integers (int), where each value corresponds 
 * to a gene, and its sign indicates the direction its 
 * transcription occurs.
 * 
 * Among other features, this class makes it easier to 
 * convert one of the genomes to an identity permutation, 
 * while making the labels of the other genomes consistent 
 * with these new labels.
 * 
 * This class can also create random genomes with arbitrary
 * amount of genes and chromosomes, which can be useful 
 * for testing.
 * 
 *******************************************************/

#include <iostream>
#include <unordered_map>
#include <vector>
#include <cstdlib>   // exit
#include <cmath>	 // abs
#include <memory>	 // shared_ptr
#include <utility>   // move, pair
#include <random>
#include <algorithm> // shuffle

#include "genome.hpp"

/*******************************************************
 *  Auxiliary functions to create random genomes.
*******************************************************/

// Creates a random permutation of [begElem]..n+[begElem].
std::vector<int> createRandomPermutation(std::mt19937& rng, const int n, int begElem) {
	// Create a vector with integers from 1 to n
	std::vector<int> permutation(n,begElem);
	for (int i = 0; i < n; ++i) {permutation[i] += i;}
	// Shuffle the vector to create a random permutation.
	std::shuffle(permutation.begin(), permutation.end(), rng);
	return permutation;
}

// Creates a vector of boolean with chance [probRev] of being true and 1-[probRev] of being false.
std::vector<bool> createRandomSigns(std::mt19937& rng, const int n, const double probRev) {
	if ((probRev < 0) || (probRev > 1)){
		std::cout << "ERROR! Probability of a reversed gene must be a value in the range [0, 1]. Value found=" << probRev << "\nProgram is aborting." << std::endl;
		exit(1);
	}
	std::uniform_real_distribution<> distr(0.0, 1.0); // Distribution for probability
	// Create a vector with integers from 1 to n
	std::vector<bool> signs(n);
	for (int i = 0; i < n; ++i) {signs[i] = (distr(rng) < probRev);}
	return signs;
}

// Split n into m bins.
std::vector<int> createRandomPartition(std::mt19937& rng, int n, const int m) {
	if(n < m) {
		std::cout << "ERROR! There are more partitions than elements. Values found: elements(n)=" << n << "; partitions(m)=" << m << "\nProgram is aborting." << std::endl;
		exit(1);
	}
	// Every partition has at least one element.
	std::vector<int> partitions(m, 1);
	n = n - m;
	// Distribute remaining elements in a random way.
	if(n > 0){
		// Uniform distribution on the closed interval [1,m].
		std::uniform_int_distribution distr(0, m-1); 
		for (int i = 0; i < n; ++i) {partitions[distr(rng)] += 1;}
	}
	return partitions;
}

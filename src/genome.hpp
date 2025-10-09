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

#pragma once // It avoids class redefinition.

#include <iostream>
#include <unordered_map>
#include <vector>
#include <cstdlib>   // exit
#include <cmath>	 // abs
#include <memory>	 // shared_ptr
#include <utility>   // move, pair
#include <random>
#include <algorithm> // shuffle

/*******************************************************
 *  Auxiliary functions to create random genomes.
*******************************************************/

// Creates a random permutation of [begElem]..n+[begElem].
std::vector<int> createRandomPermutation(std::mt19937& rng, const int n, int begElem=1);

// Creates a vector of boolean with chance [probRev] of being true and 1-[probRev] of being false.
std::vector<bool> createRandomSigns(std::mt19937& rng, const int n, const double probRev);

// Split n into m bins.
std::vector<int> createRandomPartition(std::mt19937& rng, int n, const int m);

/*******************************************************
 *  Auxiliary functions to create genomes from permutations.
*******************************************************/

// Get the gene order and the gene signs from a given permutation.
std::pair<std::vector<int>,std::vector<bool>> permToGenome(std::vector<int> permutation);

/*******************************************************
 *  Genome data structure.
*******************************************************/

// Map between the labels used in the input and the labels
// used internally by the class. The class relabel the genes
// such that one of the genomes has an ordered sequence 2, 3, .., n+1.
// Labels 1 and n+2 are reserved to extended caps used later.
template <typename GeneLabelT>
class GeneLabelMap {
public:
	// Given a custom label to a gene (string, int, etc.), it 
	// returns a signed integer. The integer value denotes the 
	// gene id used internally. The sign of the integer, if negative (positive),
	// means that originally the gene had negative (positive) sign.
	std::shared_ptr<std::unordered_map<GeneLabelT,int>> labels;
	std::shared_ptr<std::vector<GeneLabelT>> ids_to_labels;

	GeneLabelMap(std::unordered_map<GeneLabelT, int>* gene_labels, std::vector<GeneLabelT>* gene_ids_to_labels):labels(gene_labels),ids_to_labels(gene_ids_to_labels){

	}
	int add(const GeneLabelT& gene_label, const bool reversed){
		const int gene_id = (reversed ? -(ids_to_labels->size()+1) : (ids_to_labels->size()+1));
		(*labels)[gene_label] = gene_id;
		(*ids_to_labels).emplace_back(gene_label);
		return gene_id;
	}
	inline bool labelExists(const GeneLabelT& gene_label) const {
		return (labels->find(gene_label) != labels->end());
	}
	inline GeneLabelT getLabel(const int gene_id) const {
		return (*ids_to_labels)[std::abs(gene_id)-1];
	}
	inline int getId(const GeneLabelT& gene_label) const {
		return (*labels)[gene_label];
	}
	// Depending on the current orientation of the gene,
	// check if it corresponds to the ideal orientation.
	// Even if the gene is reversed, if the original orientation
	// is also reversed, it returns a positive integer.
	int getAdjustedId(const GeneLabelT& gene_label, bool const reversed) const {
		if (labels->find(gene_label) != labels->end()) {
			int gene_id = getId(gene_label);
			// Sign: [true]  = the reference gene and the query gene have different orientations. 
			//       [false] = the reference gene and the query gene have the same orientation. 
			const bool sign = ((gene_id < 0) != reversed);
			gene_id = std::abs(gene_id);
			return (sign ? -gene_id : gene_id);
		} else {
			return 0;
		}
	}
	int getAdjustedId(const int gene_id) const {
		GeneLabelT gene_label = getLabel(gene_id);
		return getAdjustedId(gene_label, (gene_id < 0));
	}
	std::string getLabelStr(const int gene_id) const {
		GeneLabelT gene_label = getLabel(gene_id);
		return (((getAdjustedId(gene_label, (gene_id < 0)) < 0) ? "-" : "+") + std::to_string(gene_label));
	}
};

template <typename GeneLabelT>
class GenomeMultichrom {
public:

	// Map between the labels used in the input and the labels
	// used internally by the class. The class relabel the genes
	// such that one of the genomes has an ordered sequence 2, 3, .., n+1.
	// Labels 1 and n+2 are reserved to extended caps used later.
	GeneLabelMap<GeneLabelT> gene_labels_map;

	std::vector<std::vector<int>> genome;

	int n{0}; // Number of genes.
	int m{0}; // Number of chromosomes.

	// Create an internal representation of a multichromosomal genome.
	// A genome with m chromosomes is represented by a vector containing
	// m vectors. Each chromosome is a vector.
	// The orientation of the genes (true if reversed, false otherwise)
	// is indicated in the vector ``genome_orientation``. 
	// The orientation of the gene located in genome_multichrom[i][j] 
	// can be found in genome_orientation[i][j].
	// In the internal representation, genes are represented as integers, 
	// appearing in increasing order: 1, ..., n.
	GenomeMultichrom(const std::vector<std::vector<GeneLabelT>>& genome_multichrom, const std::vector<std::vector<bool>>& genome_orientation):genome(genome_multichrom.size()),gene_labels_map(new typename std::unordered_map<GeneLabelT, int>,new typename std::vector<GeneLabelT>){
		// Compute the number of genes, and check if inputs make sense.
		initializeGenome(genome_multichrom, genome_orientation);
		// Create chromosomes and save mapping between gene labels and internal gene ids.
		createPermutation(genome_multichrom, genome_orientation, false);
	}

	// Create genome **without** update mapping between gene labels and internal ids.
	GenomeMultichrom(const std::vector<std::vector<GeneLabelT>>& genome_multichrom, const std::vector<std::vector<bool>>& genome_orientation, GeneLabelMap<GeneLabelT>& gene_labels_map_):genome(genome_multichrom.size()),gene_labels_map(gene_labels_map_){
		// Compute the number of genes, and check if inputs make sense.
		initializeGenome(genome_multichrom, genome_orientation);
		// Create chromosomes **without** creating a new mapping between gene labels and internal gene ids.
		createPermutation(genome_multichrom, genome_orientation, true);
	}

	// A simpler alternative to constructor if genome is unichromosomal.
	GenomeMultichrom(const std::vector<GeneLabelT>& genome_unichrom, const std::vector<bool>& genome_orientation_unichrom):genome(1),gene_labels_map(new typename std::unordered_map<GeneLabelT, int>,new typename std::vector<GeneLabelT>){
		// Make genome ``multichromosomal``.
		const std::vector<std::vector<GeneLabelT>> genome_multichrom{genome_unichrom};
		const std::vector<std::vector<bool>> genome_orientation{genome_orientation_unichrom};
		// Compute the number of genes, and check if inputs make sense.
		initializeGenome(genome_multichrom, genome_orientation);
		// Create chromosomes and save mapping between gene labels and internal gene ids.
		createPermutation(genome_multichrom, genome_orientation, false);
	}

	// A simpler alternative to constructor if genome is unichromosomal.
	GenomeMultichrom(const std::vector<GeneLabelT>& genome_unichrom, const std::vector<bool>& genome_orientation_unichrom, GeneLabelMap<GeneLabelT>& gene_labels_map_):genome(1),gene_labels_map(gene_labels_map_){
		// Make genome ``multichromosomal``.
		const std::vector<std::vector<GeneLabelT>> genome_multichrom{genome_unichrom};
		const std::vector<std::vector<bool>> genome_orientation{genome_orientation_unichrom};
		// Compute the number of genes, and check if inputs make sense.
		initializeGenome(genome_multichrom, genome_orientation);
		// Create chromosomes **without** creating a new mapping between gene labels and internal gene ids.
		createPermutation(genome_multichrom, genome_orientation, true);
	}

	// Create a random genome with n genes and m chromosomes.
	GenomeMultichrom(std::mt19937& rng, int n, int m, double probRev):n(n),m(m),genome(m),gene_labels_map(new typename std::unordered_map<GeneLabelT, int>,new typename std::vector<GeneLabelT>){
		// Initialize a random genome.
		std::pair<std::vector<std::vector<int>>,std::vector<std::vector<bool>>> rdmgen = initializeRandomGenome(rng, n, m, probRev);
		// Create chromosomes and save mapping between gene labels and internal gene ids.
		createPermutation(rdmgen.first, rdmgen.second, false);
	}

	// Create a random genome with n genes and m chromosomes.
	GenomeMultichrom(std::mt19937& rng, int n, int m, double probRev, GeneLabelMap<GeneLabelT>& gene_labels_map_):n(n),m(m),genome(m),gene_labels_map(gene_labels_map_){
		// Initialize a random genome.
		std::pair<std::vector<std::vector<int>>,std::vector<std::vector<bool>>> rdmgen = initializeRandomGenome(rng, n, m, probRev);
		// Create chromosomes **without** creating a new mapping between gene labels and internal gene ids.
		createPermutation(rdmgen.first, rdmgen.second, true);
	}

	// Creates an extended signed permutation, where two 
	// additional gene extremities are added: 1 and n+2.
	// The label of gene i becomes i+1.
	std::vector<int> getExtendedGenome(){
		std::vector<int> extendedPerm(n+2);
		// Extended elements.
		extendedPerm[0]   = 1;
		extendedPerm[n+1] = n+2;
		// Main elements.
		int idx_perm = 1;
		for(std::vector<int> const &chrom : genome) {
			for(int const &g_id : chrom) {
				extendedPerm[idx_perm++] = (g_id < 0) ? (g_id-1) : (g_id+1);
			}
		}
		return extendedPerm;
	}

	// Creates an extended unsigned permutation, where a gene i (1 <= i <= n)
	// is represented by its gene extremities i and i+1.
	// Two additional gene extremities are added: 0 and 2n+1.
	// Thus, the total set of extremities are {0, 1, 2, ..., 2n, 2n+1}.
	std::vector<int> getUnsignedExtendedPerm() {
		std::vector<int> unsignedExtPerm(2*(n+1));
		// Extended elements.
		unsignedExtPerm[0]     = 0;
		unsignedExtPerm[2*n+1] = 2*n+1;
		// Main elements.
		int idx_perm = 1;
		for(std::vector<int> const &chrom : genome) {
			for(int const &g_id : chrom) {
				if(g_id < 0) {
					unsignedExtPerm[idx_perm++] = 2*std::abs(g_id);
					unsignedExtPerm[idx_perm++] = 2*std::abs(g_id)-1;
				} else {
					unsignedExtPerm[idx_perm++] = 2*std::abs(g_id)-1;
					unsignedExtPerm[idx_perm++] = 2*std::abs(g_id);
				}

			}
		}
		return unsignedExtPerm;
	}

	std::pair<std::vector<std::vector<int>>, std::vector<std::vector<bool>>> initializeRandomGenome(std::mt19937& rng, int n, int m, const double probRev) {
		// Initialize a random genome.
		std::vector<std::vector<int>>  genome_multichrom;
		std::vector<std::vector<bool>> genome_orientation;
		std::vector<int> chrom_sizes = createRandomPartition(rng, n, m);
		int begElem=1; m=0;
		for(int const &chrom_size : chrom_sizes) {
			(genome[m++]).resize(chrom_size);
			genome_multichrom.emplace_back(createRandomPermutation(rng, chrom_size, begElem));
			genome_orientation.emplace_back(createRandomSigns(rng, chrom_size, probRev));
			begElem += chrom_size;
		}
		return std::make_pair(genome_multichrom,genome_orientation);
	}

	// Compute the number of genes, and check if inputs make sense.
	void initializeGenome(const std::vector<std::vector<GeneLabelT>>& genome_multichrom, const std::vector<std::vector<bool>>& genome_orientation){
		m = 0; n = 0;
		for(std::vector<GeneLabelT> const &chrom : genome_multichrom) {
			n += chrom.size();
			(genome[m]).resize(chrom.size());
			if(chrom.size() != genome_orientation[m].size()){
				std::cout << "ERROR! The vector of gene labels and the vector of gene orientations given as input have a different number of elements in the chromosome occupying the " << (m+1) << " position:\n\tNumber of gene labels at chromosome " << (m+1) << ": " << chrom.size() << ";\n\tNumber of gene orientations at chromosome " << (m+1) << ": " << genome_orientation[m].size() << ".\nProgram is aborting." << std::endl;
				exit(1);
			}
			++m;
		}
	}
	
	void createPermutation(const std::vector<std::vector<GeneLabelT>>& genome_multichrom, const std::vector<std::vector<bool>>& genome_orientation, const bool labelShouldExist){
		// Create chromosomes and save mapping between gene labels and internal gene ids.
		if(!labelShouldExist){gene_labels_map.ids_to_labels->reserve(n);}
		int chrom_idx = 0;
		for(std::vector<int> const &chrom : genome_multichrom) {createChrom(chrom_idx,chrom,genome_orientation[chrom_idx],labelShouldExist);++chrom_idx;}
	}

	int getGeneId(GeneLabelT const &geneLabel, bool const &geneSign){
		if(!gene_labels_map.labelExists(geneLabel)) {
			std::cout << "ERROR! The gene with label '" << geneLabel << "' was not found in another genome. For now, all genomes must have the same gene content, and missing genes are not supported. Program is aborting." << std::endl;
			exit(1);			
		} 
		return gene_labels_map.getAdjustedId(geneLabel, geneSign); // Retrieve ID corresponding to the label.
	}
	
	int createGeneId(GeneLabelT const &geneLabel, bool const &geneSign) {
		if(gene_labels_map.labelExists(geneLabel)) {
			std::cout << "ERROR! The gene with label '" << geneLabel << "' appears more than once in the genome. For now, genomes with duplicated genes are not supported. Program is aborting." << std::endl;
			exit(1);
		}
		return std::abs(gene_labels_map.add(geneLabel, geneSign)); // Create a new map entry.
	}

	// Create chromosome and update mapping between gene labels and internal ids.
	void createChrom(int const chrom_idx, std::vector<int> const &chrom, std::vector<bool> const &chrom_orientation, const bool labelShouldExist) {
		int gene_idx = 0;
		for(GeneLabelT const &gene : chrom) {
			const int gene_id = labelShouldExist ? getGeneId(gene, chrom_orientation[gene_idx]) : createGeneId(gene, chrom_orientation[gene_idx]);
			genome[chrom_idx][gene_idx] = gene_id;
			++gene_idx;
		}
	}

	void printGenome(){
		int chr_idx = 0;
		for(std::vector<int> const &chrom : genome) {
			std::cout << "> Chrom. " << ++chr_idx << std::endl;
			for(int const &g_id : chrom) {
				std::cout << ((g_id < 0) ? "-" : "+") << std::abs(g_id) << " ";
			}
			std::cout << std::endl;
		}
	}

	void printOriginalGenome(){
		int chr_idx = 0;
		for(std::vector<int> const &chrom : genome) {
			std::cout << "> Chrom. " << ++chr_idx << std::endl;
			for(int const &g_id : chrom) {
				const GeneLabelT gene_label = gene_labels_map.getLabel(std::abs(g_id));
				const int g_id_ideal = gene_labels_map.getId(gene_label);
				const bool reversed  = ((g_id_ideal < 0) != (g_id < 0));
				std::cout << (reversed ? "-" : "+") << gene_label << " ";
			}
			std::cout << std::endl;
		}
	}

};

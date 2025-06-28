/*******************************************************
 * One genome is represented as a list of permutations,
 * one per chromosome. 
 * 
 * Each permutation is represented as a vector of integers (int),
 * where each int corresponds to a gene, and its sign
 * indicates the direction its transcription occurs.
 * 
 * This class implements the data structure described in 
 * Bader, Moret, and Yan (2001), called ``overlap forest``.
 * Given two genomes represented as described earlier, 
 * the overlap forest is a way to find the connected 
 * components of the overlap graph (i.e. cycles in the
 * cycle graph that partially overlap) in linear time.
 * 
 * Finding the connected components of the overlap graph
 * in an efficient way is an important step in solving
 * several genome rearrangements problems including,
 * for example, the sorting by reversals problem.
 * 
 * Each of these connected components (if they are oriented)
 * can be handled independently, facilitating the 
 * computation of a solution.
 * 
 * For more details on the data structures and the algorithm 
 * to find connected components of the overlap graph, 
 * please check:
 * 
 * David Bader et al., 2001: "A Linear-Time Algorithm for 
 * Computing Inversion Distance between Signed Permutations 
 * with an Experimental Study";
 * 
 *******************************************************/

#include <iostream>
#include <unordered_map>
#include <vector>
#include <cstdlib>   // exit
#include <cmath>     // abs
#include <memory>    // shared_ptr
#include <utility>   // move, pair

// Map between the labels used in the input and the labels
// used internally by the class. The class relabel the genes
// such that one of the genomes has an ordered sequence 2, 3, .., n+1.
// Labels 1 and n+2 are reserved to extended caps used later.

template <typename GeneLabelT>
class GeneLabelMapKey {
public:
	GeneLabelT label;
	bool reversed{false};
	GeneLabelMapKey(const GeneLabelT& label, const bool reversed):label(label),reversed(reversed){
	}
	// Overload the '<' operator.
	friend bool operator<(const GeneLabelMapKey &a, const GeneLabelMapKey &b){
		return a.label < b.label;
	}
	// Overload the '>' operator.
	friend bool operator>(const GeneLabelMapKey &a, const GeneLabelMapKey &b){
		return a.label > b.label;
	}
	// Overload the '==' operator.
	friend bool operator==(const GeneLabelMapKey &a, const GeneLabelMapKey &b){
		return a.label == b.label;
	}
};

// Specialization of std::hash for GeneLabelMapKey.
template <typename GeneLabelT>
struct std::hash<GeneLabelMapKey<GeneLabelT>> {
	std::size_t operator () (const GeneLabelMapKey<GeneLabelT>& key) const {
		return std::hash<GeneLabelT>{}(key.label); 
    }
};

template <typename GeneLabelT>
class GeneLabelMap {
public:
	std::shared_ptr<std::unordered_map<GeneLabelMapKey<GeneLabelT>,int>> labels;
	std::shared_ptr<std::vector<GeneLabelMapKey<GeneLabelT>>> ids_to_labels;
	GeneLabelMap(std::unordered_map<GeneLabelMapKey<GeneLabelT>, int>* gene_labels, std::vector<GeneLabelMapKey<GeneLabelT>>* gene_ids_to_labels):labels(gene_labels),ids_to_labels(gene_ids_to_labels){

	}
	int add(const GeneLabelT& gene_label, const bool sign){
		GeneLabelMapKey<GeneLabelT> key = GeneLabelMapKey<GeneLabelT>(gene_label,sign);
		const int gene_id = ids_to_labels->size()+1;
		(*labels)[key] = gene_id;
		(*ids_to_labels).emplace_back(key);
		return gene_id;
	}
	inline bool labelExists(const GeneLabelT& gene_label) const {
		GeneLabelMapKey<GeneLabelT> key{gene_label, true};
		return (labels->find(key) != labels->end());
	}
	int getId(const GeneLabelT& gene_label, bool const reversed) const {
		GeneLabelMapKey<GeneLabelT> key{gene_label, reversed};
		if (labels->find(key) != labels->end()) {
			int gene_id = (*labels)[key];
			const bool reversed = (((gene_id < 0) && reversed) || ((gene_id > 0) && !reversed));
			return (reversed ? -gene_id : gene_id);
		} else {
			return 0;
		}
	}
	GeneLabelMapKey<GeneLabelT> getGene(const int gene_id) const {
		GeneLabelMapKey<GeneLabelT> key = (*ids_to_labels)[std::abs(gene_id)-1];
		const bool reversed = (((gene_id < 0) && key.reversed) || ((gene_id > 0) && !key.reversed));
		return GeneLabelMapKey<GeneLabelT>(key.label,reversed);
	}
	std::string getLabel(const int gene_id) const {
		GeneLabelMapKey<GeneLabelT> key = getGene(gene_id);
		return (key.reversed ? "-" : "+") + std::to_string(key.label);
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
	// A genome with m chromosomes is represented by a vector with
	// m vectors. Each chromosome is a vector.
	// The orientation of the genes (true if reversed, false otherwise)
	// is indicated in the vector ``genome_orientation``. 
	// The orientation of the gene located in genome_multichrom[i][j] 
	// can be found in genome_orientation[i][j].
	// In the internal representation, genes are represented as integers, 
	// appearing in increasing order: 1, ..., n.
	GenomeMultichrom(const std::vector<std::vector<GeneLabelT>>& genome_multichrom, const std::vector<std::vector<bool>>& genome_orientation):genome(genome_multichrom.size()),gene_labels_map(new typename std::unordered_map<GeneLabelMapKey<GeneLabelT>, int>,new typename std::vector<GeneLabelMapKey<GeneLabelT>>){

		// Compute the number of genes, and check if inputs make sense.
		initializeGenome(genome_multichrom, genome_orientation);
		// Create chromosomes and save mapping between gene labels and internal gene ids.
		gene_labels_map.ids_to_labels->reserve(n); 
		int chrom_idx = 0;
		for(std::vector<int> const &chrom : genome_multichrom) {createChrom(chrom_idx,chrom,genome_orientation[chrom_idx],false);++chrom_idx;}
	}

	// Create genome **without** update mapping between gene labels and internal ids.
	GenomeMultichrom(const std::vector<std::vector<int>>& genome_multichrom, const std::vector<std::vector<bool>>& genome_orientation, GeneLabelMap<GeneLabelT>& gene_labels_map_):genome(genome_multichrom.size()),gene_labels_map(gene_labels_map_){
		// Compute the number of genes, and check if inputs make sense.
		initializeGenome(genome_multichrom, genome_orientation);
		// Create chromosomes **without** creating a new mapping between gene labels and internal gene ids.
		int chrom_idx = 0;
		for(std::vector<int> const &chrom : genome_multichrom) {createChrom(chrom_idx,chrom,genome_orientation[chrom_idx],true);++chrom_idx;}
	}

	// GenomeMultichrom(const std::vector<std::vector<int>>& genome_unichrom){
	// }

	// GenomeMultichrom(const std::vector<std::vector<int>>& genome_unichrom, map){
	// }

	// // Create a random genome with n genes and m chromosomes.
	// GenomeMultichrom(int n, int m){

	// }


	// void getExtendedGenome(){

	// }

	// void getUnsignedExtendedPerm(){

	// }
	
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
	
	int getGeneId(GeneLabelT const &geneLabel, bool const &geneSign){
		if(!gene_labels_map.labelExists(geneLabel)) {
			std::cout << "ERROR! The gene with label '" << geneLabel << "' was not found in another genome. For now, all genomes must have the same gene content, and missing genes are not supported. Program is aborting." << std::endl;
			exit(1);			
		} 
		return gene_labels_map.getId(geneLabel, geneSign); // Retrieve ID correspoding to the label.
	}
	
	int createGeneId(GeneLabelT const &geneLabel, bool const &geneSign) {
		if(gene_labels_map.labelExists(geneLabel)) {
			std::cout << "ERROR! The gene with label '" << geneLabel << "' appears more than once in the genome. For now, genomes with duplicated genes are not supported. Program is aborting." << std::endl;
			exit(1);
		}
		return gene_labels_map.add(geneLabel, geneSign); // Create a new map entry.
	}

	// Create chromosome and update mapping between gene labels and internal ids.
	void createChrom(int const chrom_idx, std::vector<int> const &chrom, std::vector<bool> const &chrom_orientation, const bool labelShouldExist) {
		int gene_idx = 0;
		for(GeneLabelT const &gene : chrom) {
			const int gene_id = labelShouldExist ? getGeneId(gene, chrom_orientation[gene_idx]) : createGeneId(gene, chrom_orientation[gene_idx]);
			std::cout << gene << "=" << gene_id << std::endl;
			genome[chrom_idx][gene_idx] = gene_id;
			++gene_idx;
		}
	}

	void printGenome(){
		int chr_idx = 0;
		for(std::vector<int> const &chrom : genome) {
			std::cout << "> Chrom. " << ++chr_idx << std::endl;
			for(int const &g_id : chrom) {
				std::cout << g_id << " ";
			}
			std::cout << std::endl;
		}
	}

	void printOriginalGenome(){
		int chr_idx = 0;
		for(std::vector<int> const &chrom : genome) {
			std::cout << "> Chrom. " << ++chr_idx << std::endl;
			for(int const &g_id : chrom) {
				std::cout << gene_labels_map.getLabel(std::abs(g_id)) << " ";
			}
			std::cout << std::endl;
		}
	}

};


int main(int argc, char* argv[]) {
	//int const n = std::stoi(argv[1]); // input (command line argument): number of genes
	
    std::vector<std::vector<int>>  genome_multichrom_A  = {{5, 4, 3, 2, 1}};
    std::vector<std::vector<bool>> genome_orientation_A = {{false, true, false, true, false}};

    std::vector<std::vector<int>>  genome_multichrom_B  = {{1, 3, 5, 2, 4}};
    std::vector<std::vector<bool>> genome_orientation_B = {{true, true, true, true, true}};

	GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
	GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);

    genome_A.printGenome();
    genome_A.printOriginalGenome();

    genome_B.printGenome();
    genome_B.printOriginalGenome();

	std::cout << "Bye bye\n";
	return 0;
}

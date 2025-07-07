/*******************************************************
 * This class implements a method to sort two genomes
 * by reversals, i.e., it finds a minimum number of 
 * reversals to transform one of the input genomes 
 * into the other genome.
 * 
 * The given scenario is just one of the many possible
 * minimum scenarios that could transform one genome
 * into another.
 * 
 * To sort two genomes, this method builds up on the work
 * of many people:
 * 
 * 1) To find the connected components of the overlap 
 * graph in an efficient way, the algorithm from 
 * Bader et al. (2001) was implemented. Each component
 * can be either ``oriented`` or ``unoriented``.
 * An unoriented component means that all genes have 
 * the same sign.
 * 
 * 2) To transform unoriented components into oriented
 * components using the minimum number of reversals
 * needed, the method from Hannehalli and Pevzner 
 * (1999) was implemented. In the literature, this method 
 * is also called ``clear hurdles``.
 * 
 * 3) To sort an oriented component using the minimum
 * number of reversals needed, the method from Tannier et al.
 * (2007) was implemented. This method runs in Θ(√(n×log(n))),
 * and it was the most efficient way to sort two genomes 
 * by reversals until very recently, when two independent 
 * groups proposed two new methods, both inspired by 
 * Tannier's method, that runs in Θ(n×log(n)).
 * 
 * For more details on each step, please check the 
 * header of the files: 
 * 
 * 1) This step is implemented in ``findComponents_Bader2001.hpp``;
 * 2) This step is implemented in ``clearUnoriented_HannenhalliPevzner1999.hpp``;
 * 3) This step is implemented in ``sortOrientedByReversals_Tannier2007.hpp``.
 * 
 *******************************************************/

#include <iostream>
#include <unordered_map>
#include <vector>
#include <stack>
#include <list>
#include <cstdlib>   // exit
#include <cmath>	 // abs
#include <memory>	 // shared_ptr

#include "genome.hpp"
#include "reversal.hpp"
#include "findComponents_Bader2001.hpp"
#include "sortOrientedByReversals_Tannier2007.hpp"
#include "solveUnoriented_HannenhalliPevzner1999.hpp"


/*******************************************************
 * Auxiliary functions.
*******************************************************/

inline int convertLabel(const int g_label, std::unordered_map<int,std::pair<int,int>>& labels_map, const bool rightmost){
	return (rightmost == (g_label > 0)) ? (int)(std::trunc((labels_map[std::abs(g_label)].second+1)/2)+1) : (int)(std::trunc((labels_map[std::abs(g_label)].first+1)/2)+1);
}

inline Reversal transformLabels(Reversal rev, std::unordered_map<int,std::pair<int,int>>& labels_map){
	return Reversal(convertLabel(rev.g_beg,labels_map,true), convertLabel(rev.g_end,labels_map,true), convertLabel(rev.g_beg_next,labels_map,false), convertLabel(rev.g_end_next,labels_map,false));
}

void printGenome(std::vector<int> perm){
	for (int const& gene : perm) {
		std::cout << gene << " ";
	}
	std::cout << std::endl;
}

/*******************************************************
 * Some basic tests.
*******************************************************/

/*******************************************************
 * Tests for genome.hpp
*******************************************************/

void testCase_MakeRandomPerm(std::mt19937& rng, int n, int m, double probRev=0.3) {
	GenomeMultichrom<int> genome_A(rng, n, m, probRev);
	GenomeMultichrom<int> genome_B(rng, n, m, probRev, genome_A.gene_labels_map);

	std::cout << "\n\nTest: Random genome constructor (n=" << n << ";m=" << m  << ";p_rev=" << probRev << ")\n";
	std::cout << "Genome A -- Original:\n";
	genome_A.printOriginalGenome();
	std::cout << "Genome A -- Internal representation:\n";
	genome_A.printGenome();
	
	std::cout << "Genome B -- Original:\n";
	genome_B.printOriginalGenome();
	std::cout << "Genome B -- Internal representation:\n";
	genome_B.printGenome();
}

void testCase_MakeMultichromGenome(){
	std::vector<std::vector<int>>  genome_multichrom_A  = {{5, 4, 3, 2, 1}};
	std::vector<std::vector<bool>> genome_orientation_A = {{false, true, false, true, false}};

	std::vector<std::vector<int>>  genome_multichrom_B  = {{1, 3, 5, 2, 4}};
	std::vector<std::vector<bool>> genome_orientation_B = {{true, true, true, true, true}};

	GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
	GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);

	std::cout << "\n\nTest: Multichromosomal genome constructor\n";
	std::cout << "Genome A -- Original:\n";
	genome_A.printOriginalGenome();
	std::cout << "Genome A -- Internal representation:\n";
	genome_A.printGenome();
	
	std::cout << "Genome B -- Original:\n";
	genome_B.printOriginalGenome();
	std::cout << "Genome B -- Internal representation:\n";
	genome_B.printGenome();
}

void testCase_MakeUnichromGenome(){
	std::vector<int>  genome_multichrom_A  = {5, 4, 3, 2, 1};
	std::vector<bool> genome_orientation_A = {false, true, false, true, false};

	std::vector<int>  genome_multichrom_B  = {1, 3, 5, 2, 4};
	std::vector<bool> genome_orientation_B = {true, true, true, true, true};

	GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
	GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);

	std::cout << "\n\nTest: Unichromosomal genome constructor\n";
	std::cout << "Genome A -- Original:\n";
	genome_A.printOriginalGenome();
	std::cout << "Genome A -- Internal representation:\n";
	genome_A.printGenome();
	
	std::cout << "Genome B -- Original:\n";
	genome_B.printOriginalGenome();
	std::cout << "Genome B -- Internal representation:\n";
	genome_B.printGenome();
	
}

/*******************************************************
 * Tests for findComponents_Bader2001.hpp
*******************************************************/

// Example used in the paper from Garg et al. (2019).
// {-2,5,4,-1,3,6,9,-7,-8}
void testCase_Garg2019(){
	std::vector<int>  genome_multichrom_A  = {1, 2, 3, 4, 5, 6, 7, 8, 9};
	std::vector<bool> genome_orientation_A = {false, false, false, false, false, false, false, false, false};

	std::vector<int>  genome_multichrom_B  = {2, 5, 4, 1, 3, 6, 9, 7, 8};
	std::vector<bool> genome_orientation_B = {true, false, false, true, false, false, false, true, true};

	GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
	GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);

	std::cout << "\n\nTest: Example from Garg et al.(2019)\n";
	std::cout << "Genome A -- Original:\n";
	genome_A.printOriginalGenome();
	std::cout << "Genome A -- Internal representation:\n";
	genome_A.printGenome();
	
	std::cout << "Genome B -- Original:\n";
	genome_B.printOriginalGenome();
	std::cout << "Genome B -- Internal representation:\n";
	genome_B.printGenome();

	ConnectedComponents comps = ConnectedComponents(genome_B.getUnsignedExtendedPerm());
}

// Example used in the paper from Bader et al. (2001).
// (+3, +9, −7, +5, −10, +8, +4, −6, +11, +2, +1)
void testCase_Bader2001(){
	std::vector<int>  genome_multichrom_A  = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
	std::vector<bool> genome_orientation_A = {false, false, false, false, false, false, false, false, false, false, false};

	std::vector<int>  genome_multichrom_B  = {3, 9, 7, 5, 10, 8, 4, 6, 11, 2, 1};
	std::vector<bool> genome_orientation_B = {false, false, true, false, true, false, false, true, false, false, false};

	GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
	GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);

	std::cout << "\n\nTest: Example from Garg et al.(2019)\n";
	std::cout << "Genome A -- Original:\n";
	genome_A.printOriginalGenome();
	std::cout << "Genome A -- Internal representation:\n";
	genome_A.printGenome();
	
	std::cout << "Genome B -- Original:\n";
	genome_B.printOriginalGenome();
	std::cout << "Genome B -- Internal representation:\n";
	genome_B.printGenome();

	ConnectedComponents comps = ConnectedComponents(genome_B.getUnsignedExtendedPerm());
}

/*******************************************************
* Tests for sortOrientedByReversals_Tannier2007.hpp
*******************************************************/

void testCase_SortOrientedComponent(int const n) {
	// Start vector with Identity permutation (for testing).
	std::vector<int> perm(n);
	std::iota(perm.begin(), perm.end(), 1);	
	GenomeSort genomeSort = GenomeSort(perm);
	genomeSort.applyReversal(2, 9);   // after rev: 1 2 -9 -8 -7 -6   -5 -4 -3 10 11 .. 20
	genomeSort.applyReversal(-6, -5); // after rev: 1 2 -9 -8 -7 -6    5 -4 -3 10 11 .. 20
	genomeSort.applyReversal(-6, 10); // after rev: 1 2 -9 -8 -7 -6  -10  3  4 -5 11 .. 20
}

// Example used in the paper from Tannier et al. (2007) (Figure 4).
// {0, -1, 3, 2, 4}
void testCase_Tannier2007() {
	// Permutation **must** start at 1: [1 2 .. gene]
	std::vector<int> perm{1, -2, 4, 3, 5};
	GenomeSort genomeSort = GenomeSort(perm);
	std::deque<Reversal> allrev = genomeSort.sortByReversals();
	std::cout << "Sorting by reversals---Solution" << std::endl;
	for(Reversal const &rev : allrev) {
		std::cout << "(gene " << rev.g_beg << ", gene " << rev.g_end << "]" << std::endl;
	}
}

/*******************************************************
* Tests for solveUnoriented_HannenhalliPevzner1999.hpp
*******************************************************/

// Example used in the paper from Hannehalli and Pevzner (1999) (Figure 4(a)).
// {+5, +7, +6, +8, +1, +3, +2, +4}
void testCase_Hannehalli1999_Fig4a(std::mt19937& rng){
	std::vector<int>  genome_multichrom_A  = {1, 2, 3, 4, 5, 6, 7, 8};
	std::vector<bool> genome_orientation_A = {false, false, false, false, false, false, false, false};

	std::vector<int>  genome_multichrom_B  = {5, 7, 6, 8, 1, 3, 2, 4};
	std::vector<bool> genome_orientation_B = {false, false, false, false, false, false, false, false};

	GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
	GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);

	std::cout << "\n\nTest: Example from Hannehalli and Pevzner (1999) - Figure 4(a)\n";
	std::cout << "Genome A -- Original:\n";
	genome_A.printOriginalGenome();
	std::cout << "Genome A -- Internal representation:\n";
	genome_A.printGenome();
	
	std::cout << "Genome B -- Original:\n";
	genome_B.printOriginalGenome();
	std::cout << "Genome B -- Internal representation:\n";
	genome_B.printGenome();

	//////////////////////////////////////////
	// Part I: Unoriented components -> oriented components.
	GenomePermutation<BlockSimple> genperm(genome_B.getExtendedGenome());
	// Find connected components.
	ConnectedComponents comps = ConnectedComponents(genperm.getUnsignedExtendedPerm());
	// Transform unoriented components into oriented components using the minimum number of reversals.
	UnorientedComponents comps_unoriented = UnorientedComponents(genperm, comps);
	std::vector<Reversal> reversals = comps_unoriented.clearUnorientedComponents(rng);

	//////////////////////////////////////////

	std::cout << "Genome B -- Oriented unextended:\n";
	printGenome(genperm.getUnextendedPerm());

	std::cout << "Genome B -- Oriented extended:\n";
	printGenome(genperm.getExtendedPerm());

	std::cout << "Genome B -- Unsigned extended:\n";
	printGenome(genperm.getUnsignedExtendedPerm());

	/////////////////////////////////////////////
	// Part II: Sort oriented components.
	// Find connected components again (they should be all oriented now).
	comps = ConnectedComponents(genperm.getUnsignedExtendedPerm());
	// Sort each connected component separately.
	std::cout << "\nSort connected components by reversals" << std::endl;
	for(const int& root_idx: comps.rootList){
		comps.printComponent(comps.forest[root_idx], "", comps.perm.size(), comps.forest);
		std::unordered_map<int,std::pair<int,int>> newlabels_map;
		std::vector<int> perm = genperm.getExtendedPerm(comps.forest[root_idx].genes,newlabels_map);
		std::cout << "Extended permutation: " << std::endl;
		for(const int& g: perm){std::cout << g << "[" << newlabels_map[std::abs(g)].first << "," << newlabels_map[std::abs(g)].second << "] ";}
		std::cout << std::endl;
		// Sort component.
		std::cout << "Sort permutation..." << std::endl;
		GenomeSort genomeSort = GenomeSort(perm);
		std::deque<Reversal> reversalsPerComp = genomeSort.sortByReversals();
		// Apply reversals to the permutation.
		std::cout << "Save reversals..." << std::endl;	
		// printGenome(genperm.getExtendedPerm());
		for(Reversal const &rev : reversalsPerComp) {
			Reversal rev_ = transformLabels(rev, newlabels_map);
			std::cout << "(" << rev.g_beg << "{" << rev_.g_beg << "}, " << rev.g_end << "{" << rev_.g_end << "}]" << std::endl;
			applyReversal(genperm, rev_.g_beg, rev_.g_end);
			reversals.emplace_back(rev_);
			// printGenome(genperm.getExtendedPerm());
		}
		genperm.clearBlockStatus();
	}

	std::cout << "\nSorting by reversals - Solution" << std::endl;
	GenomePermutation<BlockSimple> genperm_final(genome_B.getExtendedGenome());
	printGenome(genperm_final.getExtendedPerm());
	for(Reversal const &rev : reversals) {
		std::cout << "(" << rev.g_beg << "," << rev.g_end << "]" << std::endl;
		applyReversal(genperm_final, rev.g_beg, rev.g_end);
		printGenome(genperm_final.getExtendedPerm());
	}
}

// Example used in the paper from Hannehalli and Pevzner (1999) (Figure 4(b)).
// {+2, +4, +3, +5, +7, +6, +8, +1}
void testCase_Hannehalli1999_Fig4b(std::mt19937& rng){
	std::vector<int>  genome_multichrom_A  = {1, 2, 3, 4, 5, 6, 7, 8};
	std::vector<bool> genome_orientation_A = {false, false, false, false, false, false, false, false};

	std::vector<int>  genome_multichrom_B  = {2, 4, 3, 5, 7, 6, 8, 1};
	std::vector<bool> genome_orientation_B = {false, false, false, false, false, false, false, false};

	GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
	GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);

	std::cout << "\n\nTest: Example from Hannehalli and Pevzner (1999) - Figure 4(b)\n";
	std::cout << "Genome A -- Original:\n";
	genome_A.printOriginalGenome();
	std::cout << "Genome A -- Internal representation:\n";
	genome_A.printGenome();
	
	std::cout << "Genome B -- Original:\n";
	genome_B.printOriginalGenome();
	std::cout << "Genome B -- Internal representation:\n";
	genome_B.printGenome();

	// Find connected components.
	// ConnectedComponents comps = ConnectedComponents(genome_B.getUnsignedExtendedPerm());

	// Transform unoriented components into oriented components using the minimum number of reversals.
	// UnorientedComponents comps_unoriented = UnorientedComponents(genome_B.getExtendedGenome(), comps);
	// std::vector<Reversal> reversals = comps_unoriented.clearUnorientedComponents(rng);
}

// Example used in the book ``Mathematics of Evolution and Phylogeny`` (2005) (Section 10.4.2).
// {0, 2, 1, 3, 5, 7, 6, 8, 9, 4, 10}
void testCase_Bergeron2005(std::mt19937& rng){
	std::vector<int>  genome_multichrom_A  = {1, 2, 3, 4, 5, 6, 7, 8, 9};
	std::vector<bool> genome_orientation_A = {false, false, false, false, false, false, false, false, false};

	std::vector<int>  genome_multichrom_B  = {2, 1, 3, 5, 7, 6, 8, 9, 4};
	std::vector<bool> genome_orientation_B = {false, false, false, false, false, false, false, false, false};

	GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
	GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);

	std::cout << "\n\nTest: Example from the book 'Mathematics of Evolution and Phylogeny' (2005)\n";
	std::cout << "Genome A -- Original:\n";
	genome_A.printOriginalGenome();
	std::cout << "Genome A -- Internal representation:\n";
	genome_A.printGenome();
	
	std::cout << "Genome B -- Original:\n";
	genome_B.printOriginalGenome();
	std::cout << "Genome B -- Internal representation:\n";
	genome_B.printGenome();

	// Find connected components.
	// ConnectedComponents comps = ConnectedComponents(genome_B.getUnsignedExtendedPerm());

	// Transform unoriented components into oriented components using the minimum number of reversals.
	// UnorientedComponents comps_unoriented = UnorientedComponents(genome_B.getExtendedGenome(), comps);
	// std::vector<Reversal> reversals = comps_unoriented.clearUnorientedComponents(rng);
}

int main(int argc, char* argv[]) {
	
	if (argc < 2) {
		std::cout << "ERROR! No extra command line argument passed other than program name. Program is aborting." << std::endl;
		exit(1);
	} else if (argc != 5) {
		std::cout << "ERROR! Wrong number of arguments. Expected: ./program [seed] [#genes] [#chrom] [probRev]. Program is aborting." << std::endl;
		exit(1);
	}

	int const seed       = std::stoi(argv[1]); // input (command line argument): seed random number generator.
	int const nbgenes    = std::stoi(argv[2]); // input (command line argument): number of genes.
	int const nbchrom    = std::stoi(argv[3]); // input (command line argument): number of chromosomes.
	double const probRev = std::stod(argv[4]); // input (command line argument): probability of gene w/negative sign.

	// Create a random number generator
	std::cout << "- Seed to reproduce tests: " << seed << std::endl;
	std::mt19937 rng(seed); // Seed the generator

	// testCase_MakeRandomPerm(rng,nbgenes,nbchrom,probRev);
	// testCase_MakeMultichromGenome();
	// testCase_MakeUnichromGenome();
	// testCase_Garg2019();
	// testCase_Bader2001();
	testCase_Hannehalli1999_Fig4a(rng);
	// testCase_Hannehalli1999_Fig4b(rng);
	// testCase_Bergeron2005(rng);

	// testCase_SortOrientedComponent(20);
	// testCase_Tannier2007();

	std::cout << "Bye bye\n";
	return 0;
}

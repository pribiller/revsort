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

#include "sortByReversals.hpp"

void SortByReversals::printGenome(std::vector<int> perm){
	for (int const& gene : perm) {
		std::cout << gene << " ";
	}
	std::cout << std::endl;
}

void SortByReversals::printInputGenomes(){
	std::cout << "Genome A -- Original:" << std::endl;
	genome_A.printOriginalGenome();
	std::cout << "Genome A -- Internal representation:" << std::endl;
	genome_A.printGenome();
	
	std::cout << "Genome B -- Original:" << std::endl;
	genome_B.printOriginalGenome();
	std::cout << "Genome B -- Internal representation:" << std::endl;
	genome_B.printGenome();
}

void SortByReversals::printSolution(){
	std::cout << "\nSorting by reversals - Solution" << std::endl;
	GenomePermutation<BlockSimple> genperm_final(genome_B.getExtendedGenome());
	printGenome(genperm_final.getExtendedPerm());
	for(Reversal const &rev : reversals) {
		std::cout << "(" << rev.g_beg << "," << rev.g_end << "]" << std::endl;
		applyReversal(genperm_final, rev.g_beg, rev.g_end);
		printGenome(genperm_final.getExtendedPerm());
	}
}

std::vector<Reversal> SortByReversals::sort(std::mt19937& rng){
	// Part I: Unoriented components -> oriented components.
	GenomePermutation<BlockSimple> genperm(genome_B.getExtendedGenome());
	// Find connected components.
	if(debug){std::cout << "Find connected components..." << std::endl << std::endl;}
	ConnectedComponents comps = ConnectedComponents(genperm.getUnsignedExtendedPerm(),debug);
	// Transform unoriented components into oriented components using the minimum number of reversals.
	if(debug){std::cout << "Transform unoriented components into oriented components..." << std::endl;}
	UnorientedComponents comps_unoriented = UnorientedComponents(genperm, comps);
	comps_unoriented.debug = debug;
	reversals = comps_unoriented.clearUnorientedComponents(rng);

	// For debugging.
	if(debug){
		std::cout << "Genome B -- Oriented unextended:\n";
		printGenome(genperm.getUnextendedPerm());
		std::cout << "Genome B -- Oriented extended:\n";
		printGenome(genperm.getExtendedPerm());
		std::cout << "Genome B -- Unsigned extended:\n";
		printGenome(genperm.getUnsignedExtendedPerm());
	}

	// Part II: Sort oriented components.
	// Find connected components again (they should be all oriented now).
	comps = ConnectedComponents(genperm.getUnsignedExtendedPerm(),debug);
	// Sort each connected component separately.
	if(debug){std::cout << "Sort connected components by reversals" << std::endl;}
	for(const int& root_idx: comps.rootList){
		comps.printComponent(comps.forest[root_idx], "", comps.perm.size(), comps.forest);
		std::unordered_map<int,std::pair<int,int>> newlabels_map;
		// Permutation **must** start at 1: [1 2 .. gene]
		std::vector<int> perm = genperm.getExtendedPerm(comps.forest[root_idx].genes,newlabels_map);
		if(debug){
			std::cout << "\n\nExtended permutation: " << std::endl;
			for(const int& g: perm){std::cout << g << "[" << newlabels_map[std::abs(g)].first << "," << newlabels_map[std::abs(g)].second << "] ";}
			std::cout << std::endl;
		}
		// Sort component.
		if(debug){std::cout << "Sort permutation..." << std::endl;}
		GenomeSort genomeSort = GenomeSort(perm);
		genomeSort.debug = debug;
		std::deque<Reversal> reversalsPerComp = genomeSort.sortByReversals();
		// Apply reversals to the permutation.
		if(debug){std::cout << "Save reversals..." << std::endl;}
		// printGenome(genperm.getExtendedPerm());
		for(Reversal const &rev : reversalsPerComp) {
			Reversal rev_ = convertLabels(rev, newlabels_map);
			// std::cout << "(" << rev.g_beg << "{" << rev_.g_beg << "}, " << rev.g_end << "{" << rev_.g_end << "}]" << std::endl;
			applyReversal(genperm, rev_.g_beg, rev_.g_end);
			reversals.emplace_back(rev_);
			// printGenome(genperm.getExtendedPerm());
		}
		genperm.clearBlockStatus();
	}

	if(debug){printSolution();}
	return reversals;
}

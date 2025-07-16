/*******************************************************
 * Original author: Priscila Biller
 * Created: July/2025
 * License: GPL v3
 * 
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

bool SortByReversals::printSolution(){
	bool debug_prev = debug;
	debug = true;
	bool is_correct = checkSolution();
	debug = debug_prev;
	return is_correct;
}

bool SortByReversals::printStats(){
	std::cout << std::endl << "Some stats: " << std::endl;
	std::cout << "- Running time: " << t << " [µs]." << std::endl;
	std::cout << "- Breakpoints (b) = " << nb_breakpoints << std::endl;
	std::cout << "- Non-trivial cycles (c) = " << nb_cycles_nontrivial << "; total (including trivial cycles): " << nb_cycles << std::endl;
	std::cout << "- Hurdles (h) = " << nb_hurdles << std::endl;
	std::cout << "- Components = total: " << nb_components << "; oriented: " << nb_components_oriented << "; unoriented: " << nb_components_unoriented << "." << std::endl;
	std::cout << "- Reversal distance (d) = expected: " << exp_distance << " reversals (or " << (exp_distance+1) << " if there is a fortress); found: " << obs_distance << " reversals." << std::endl;
	const bool is_correct = ((obs_distance-exp_distance) < 2);
	return is_correct;
}

// Check if final permutation is equal to the identity permutation.
bool SortByReversals::checkSolution(){
	if(debug){std::cout << "\nSorting by reversals - Solution (" << reversals.size() << " reversals)" << std::endl;}
	GenomePermutation<BlockSimple> genperm_final(genome_B.getExtendedGenome());
	if(debug){printGenome(genperm_final.getExtendedPerm());}
	for(Reversal const &rev : reversals) {
		if(debug){std::cout << "(" << rev.g_beg << "," << rev.g_end << "]" << std::endl;}
		applyReversal(genperm_final, rev.g_beg, rev.g_end);
		if(debug){printGenome(genperm_final.getExtendedPerm());}
	}
	return (genperm_final.getBreakpoints() == 0);
}

// Check if number of reversals matches the minimum number expected.
bool SortByReversals::checkDistance(){
	int const exp_distance = nb_breakpoints - nb_cycles_nontrivial + nb_hurdles;
	int const obs_distance = reversals.size();
	const bool is_correct = ((obs_distance-exp_distance) < 2);
	return is_correct;
}

// TODO: This function is very similar to another function: GenomeSort::getReversal.
//       Maybe refactoring both in the future?
Reversal SortByReversals::convertLabels(const Reversal& rev, std::unordered_map<int,std::pair<int,int>>& labels_map, const GenomePermutation<BlockSimple>& genperm){
	// std::cout << "[convertLabels] Arc: " << rev.g_arc << std::endl;
	// Check orientation of genes 1 and 2.
	bool g1_rev   = genperm.isReversed(labels_map[std::abs(rev.g_arc)].first, labels_map[std::abs(rev.g_arc)].second);
	bool g2_rev   = genperm.isReversed(labels_map[std::abs(rev.g_arc+1)].first, labels_map[std::abs(rev.g_arc+1)].second);
	// std::cout << "[convertLabels] Is arc reversed? " << (g1_rev ? "T" : "F") << " - Extremities: "<< labels_map[std::abs(rev.g_arc)].first << " and " << labels_map[std::abs(rev.g_arc)].second << std::endl;
	// std::cout << "[convertLabels] Is arc->next reversed? " << (g2_rev ? "T" : "F") << " - Extremities: "<< labels_map[std::abs(rev.g_arc+1)].first << " and " << labels_map[std::abs(rev.g_arc+1)].second << std::endl;	
	if(g1_rev==g2_rev){
		std::cout << "ERROR! Reversal in genes with same orientation: " << rev.printReversal() << ". Program is aborting."<< std::endl;
		exit(1);
	}
	// Check if gene g2 appears before gene g1 in the permutation.
	bool g1g2_rev = genperm.isReversed(labels_map[std::abs(rev.g_arc)].first, labels_map[std::abs(rev.g_arc+1)].first);
	int g_beg;
	int g_end;
	// Reversal (i, i+1] (normal case).
	if(!g1_rev){
		// std::cout << "[convertLabels] Reversal (i, i+1] (normal case): ";
		// Reversal (g2, g1].
		if (g1g2_rev) {
			g_beg = convertLabel(rev.g_arc+1, labels_map, !g2_rev); // Current rightmost extremity depends on the gene orientation.
			g_end = convertLabel(rev.g_arc,   labels_map, !g1_rev); // Current rightmost extremity depends on the gene orientation.
		// Reversal (g1, g2].
		} else {
			g_beg = convertLabel(rev.g_arc,  labels_map, !g1_rev); // Current rightmost extremity depends on the gene orientation.
			g_end = convertLabel(rev.g_arc+1,labels_map, !g2_rev); // Current rightmost extremity depends on the gene orientation.
		}
	// Reversal [i, i+1) (case needs adjustment).
	} else {
		// std::cout << "[convertLabels] Reversal [i, i+1) (case needs adjustment): ";
		// Convert labels.
		g_beg = convertLabel(rev.g_arc,  labels_map, g1_rev); // Current leftmost extremity depends on the gene orientation.
		g_end = convertLabel(rev.g_arc+1,labels_map, g2_rev); // Current leftmost extremity depends on the gene orientation.
		// Reversal [g1, g2).
		g_beg = genperm.genes[g_beg-1]->block->getPrevGene(genperm.genes[g_beg-1])->id;
		g_end = genperm.genes[g_end-1]->block->getPrevGene(genperm.genes[g_end-1])->id;
		// Reversal [g2, g1).
		if (g1g2_rev) {std::swap(g_beg,g_end);}
	}
	// std::cout << "(" << g_beg << "," << g_end << "]" << std::endl;
	return Reversal(rev.g_arc, g_beg, g_end, -1, -1);
}

/* Sort a unichromosomal genome by reversals.
   It returns ``true`` if the input permutation is 
   properly sorted in the end, and ``false`` if 
   something went wrong.
*/
void SortByReversals::sort(std::mt19937& rng){
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	// Part I: Unoriented components -> oriented components.
	GenomePermutation<BlockSimple> genperm(genome_B.getExtendedGenome());
	nb_breakpoints = genperm.getBreakpoints();
	// Find connected components.
	if(debug){std::cout << "Find connected components..." << std::endl << std::endl;}
	ConnectedComponents comps = ConnectedComponents(genperm.getUnsignedExtendedPerm(),debug);
	nb_components=comps.rootList.size();
	nb_cycles=comps.forest.size();
	// Stats on components.
	nb_cycles_nontrivial     = nb_cycles;
	nb_components_oriented   = 0;
	nb_components_unoriented = 0;
	for (int const& root_idx : comps.rootList) {
		if(comps.forest[root_idx].oriented){
			++nb_components_oriented;
		} else {
			++nb_components_unoriented;
		}
		if(comps.forest[root_idx].genes.size() == 2){--nb_cycles_nontrivial;}
	}
	// Transform unoriented components into oriented components using the minimum number of reversals.
	if(debug){std::cout << "Transform unoriented components into oriented components..." << std::endl;}
	UnorientedComponents comps_unoriented = UnorientedComponents(genperm, comps);
	comps_unoriented.debug = debug;
	reversals  = comps_unoriented.clearUnorientedComponents(rng);
	nb_hurdles = comps_unoriented.nb_hurdles;
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
		if(debug){comps.printComponent(comps.forest[root_idx], "", comps.perm.size(), comps.forest);}
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
			Reversal rev_ = convertLabels(rev, newlabels_map, genperm);
			// std::cout << "(" << rev.g_beg << "{" << rev_.g_beg << "}, " << rev.g_end << "{" << rev_.g_end << "}]" << std::endl;
			applyReversal(genperm, rev_.g_beg, rev_.g_end);
			reversals.emplace_back(rev_);
			// printGenome(genperm.getExtendedPerm());
		}
		genperm.clearBlockStatus();
	}
	exp_distance = nb_breakpoints - nb_cycles_nontrivial + nb_hurdles;
	obs_distance = reversals.size();
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	t = (std::chrono::duration_cast<std::chrono::microseconds>(end - begin)).count();
}

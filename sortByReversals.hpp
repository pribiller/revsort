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

#pragma once // It avoids class redefinition.

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

class SortByReversals {
private:

	// Input : gene A with (unsigned) extremities x_1, y_2
	// Output: gene X or gene Y, depending on which extremity is selected.
	inline int convertLabel(const int g_label, std::unordered_map<int,std::pair<int,int>>& labels_map, const bool rightmost){
		return (rightmost == (g_label > 0)) ? (int)(std::trunc((labels_map[std::abs(g_label)].second+1)/2)+1) : (int)(std::trunc((labels_map[std::abs(g_label)].first+1)/2)+1);
	}

	// inline Reversal convertLabels(Reversal rev, std::unordered_map<int,std::pair<int,int>>& labels_map){
	// 	return Reversal(convertLabel(rev.g_beg,labels_map,true), convertLabel(rev.g_end,labels_map,true), convertLabel(rev.g_beg_next,labels_map,false), convertLabel(rev.g_end_next,labels_map,false));
	// }

	// TODO: This function is very similar to another function: GenomeSort::getReversal.
	//       Maybe refactoring both in the future?
	Reversal convertLabels(const Reversal& rev, std::unordered_map<int,std::pair<int,int>>& labels_map, const GenomePermutation<BlockSimple>& genperm);

public:

	GenomeMultichrom<int>& genome_A;
	GenomeMultichrom<int>& genome_B;

	std::vector<Reversal> reversals; // It saves all reversals in the solution.

	bool debug{false};

	// Some stats.
	int nb_breakpoints{0};
	int nb_cycles{0};
	int nb_cycles_nontrivial{0};
	int nb_components{0};
	int nb_components_oriented{0};
	int nb_components_unoriented{0};
	int nb_hurdles{0};

	SortByReversals(GenomeMultichrom<int>& genome_A, GenomeMultichrom<int>& genome_B, bool debug=false):genome_A(genome_A),genome_B(genome_B),debug(debug){
		if(debug){printInputGenomes();}
	}

	void printGenome(std::vector<int> perm);
	void printInputGenomes();
	bool printSolution();
	bool printStats();
	
	// It ``replays`` the sorting scenario, starting with the input 
	// permutation and applying each reversal from the list ``reversals``.
	// If the final permutation is equal to the identity permutation,
	// it returns true; otherwise false.
	bool checkSolution();
	// Check if number of reversals matches the minimum number expected.
	bool checkDistance();

	void sort(std::mt19937& rng);

};

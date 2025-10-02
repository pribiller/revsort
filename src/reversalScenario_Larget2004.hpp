/*******************************************************
 * Original author: Priscila Biller
 * Created: October/2025
 * License: GPL v3
 * 
 * This class implements a method to sample a random
 * sorting scenario by reversals given two genomes,
 * i.e., it finds a sequence of reversals to transform 
 * one of the input genomes into the other genome.
 * Note that the sequence of reversals might not be optimal,
 * specially if the two given genomes are far apart.
 * 
 * The given scenario is just one of the many possible
 * scenarios that could transform one genome
 * into another.
 * 
 * To find a random scenario, a method similar to Larget's
 * work [Larget et al. (2004), Larget et al. (2005)] 
 * is implemented.
 * 
 * The method implements a MCMCMC (also know as 
 * MCMC with Parallel Tempering), in which each
 * several chains are ran in parallel and swaped
 * every once in a while. The chains have different
 * 'temperatures', a parameter that controls how
 * flat or rugged is the solution space. Flatter 
 * landscapes allow for bigger jumps in the solution space.
 * 
 * To find a random reversal history, at each step
 * a random reversal is sampled. Reversals are categorized
 * in 4 types: "good reversals", "neutral good reversals",
 * "neutral reversals", and "bad reversals". These types are 
 * defined based on how close the genome gets to the target 
 * genome after they are applied.
 * 
 * References
 * ----------
 * 
 * For additional information on similar methods that inspired 
 * this implementation, please check:
 * 
 * - Bret Larget et al. "A Bayesian analysis of metazoan mitochondrial 
 * genome arrangements." Molecular Biology and Evolution 22.3 (2005): 486-495.
 * 
 * - Bret Larget et al. "Bayesian phylogenetic inference from animal mitochondrial 
 * genome arrangements." Journal of the Royal Statistical Society Series B: 
 * Statistical Methodology 64.4 (2002): 681-693.
 * 
 * Less related, but also nice:
 * 
 * - Miklós, István, and Aaron E. Darling. "Efficient sampling of parsimonious 
 * inversion histories with application to genome rearrangement in Yersinia." 
 * Genome biology and evolution 1 (2009): 153-164.
 * 
 * Notice that the method implemented here is similar to, but not exactly 
 * the same as the ones mentioned above. For example, reversals that solve 
 * hurdles in an optimal way are considered "good" here, but they would be 
 * considered "bad" in the previous methods.
 * 
 *******************************************************/

#pragma once // It avoids class redefinition.

#include <iostream>
#include <unordered_map>
#include <algorithm> // set_difference
#include <vector>
#include <stack>
#include <list>
#include <cstdlib>   // exit
#include <cmath>	 // abs
#include <numeric>   // iota
#include <memory>	 // shared_ptr
#include <chrono>
#include <utility>   // swap 

#include "genome.hpp"
#include "reversal.hpp"
#include "findComponents_Bader2001.hpp"
#include "solveUnoriented_HannenhalliPevzner1999.hpp"

#include "sampleReversal_Larget2004.hpp"

class RandomReversalScenario {
// private:

public:

	GenomePermutation<BlockSimple>& genperm; // Unsigned extended permutation (it assumes that one of the permutations is the identity).
	std::vector<ReversalRandom> reversals;  // It saves all sampled reversals in the solution.

	std::vector<float> revprobs;  // Probability of each type of reversal.
	float p_stop{0.99};

	bool debug{false};

	RandomReversalScenario(GenomePermutation<BlockSimple>& genperm_, const bool debug=false, const float p_good=0.97, const float p_neutralgood=0.92, const float p_neutral=0.92, const float p_stop=0.99):genperm(genperm_),revprobs(ReversalType_COUNT, 0),p_stop(p_stop),debug(debug){
		revprobs = {p_good, p_neutralgood, p_neutral, static_cast<float>(1.0-p_neutral)};
	}
	
	std::vector<ReversalRandom> sampleScenario(std::mt19937& rng);
	std::vector<ReversalRandom> sampleModifiedScenario(std::mt19937& rng, int updpos_beg=-1, int updpos_end=-1);

	void printGenome(std::vector<int> perm);
};

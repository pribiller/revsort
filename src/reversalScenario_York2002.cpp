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
 * To find a random scenario, a method similar to York's
 * work is implemented. [see York et al. (2002); the 
 * implemented method is also similar to 
 * Larget et al. (2004), Larget et al. (2005)] 
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
 * - Thomas York, Richard Durrett, and Rasmus Nielsen. "Bayesian estimation 
 * of the number of inversions in the history of two chromosomes". Journal of 
 * Computational Biology (2002), 9(6), 805-818.
 * 
 * - Bret Larget et al. "A Bayesian analysis of metazoan mitochondrial 
 * genome arrangements." Molecular Biology and Evolution 22.3 (2005): 486-495.
 * 
 * - Bret Larget et al. "Bayesian phylogenetic inference from animal mitochondrial 
 * genome arrangements." Journal of the Royal Statistical Society Series B: 
 * Statistical Methodology 64.4 (2002): 681-693.
 * 
 * - Miklós, István, and Aaron E. Darling. "Efficient sampling of parsimonious 
 * inversion histories with application to genome rearrangement in Yersinia." 
 * Genome biology and evolution 1 (2009): 153-164.
 * 
 * Notice that the method implemented here is similar to, but not exactly 
 * the same as the ones mentioned above. It should be seen more like a "mix"
 * of the mentioned methods: 
 * 
 * 1) The generation of a reversal history is similar to Miklos et al. (2009);
 * 
 * 2) The proposal probability, acceptance probability, etc., are computed in a 
 * similar way to the method of York et al. (2002).
 * 
 * 3) The parallel tempering is similar to Larget et al. (2005).
 * 
 *******************************************************/

#include <iostream>
#include <unordered_map>
#include <vector>
#include <stack>
#include <list>
#include <cstdlib> // exit
#include <cmath>	 // abs
#include <memory>	 // shared_ptr

#include "reversalScenario_York2002.hpp"

void RandomReversalScenario::printGenome(std::vector<int> perm){
	for (int const& gene : perm) {
		std::cout << gene << " ";
	}
	std::cout << std::endl;
}

std::vector<ReversalRandom> RandomReversalScenario::sampleScenario(std::mt19937& rng){

	std::uniform_real_distribution distr(0.0, 1.0);
	int prob_rdm = distr(rng);

	std::vector<ReversalRandom> reversals;
	while((genperm.getBreakpoints() > 0) || (p_stop < prob_rdm)) { // not identity and p_stop.
		std::cout << "\n\n\n------------------------------------------------\n\n - Current nb. of breakpoints = " << genperm.getBreakpoints() << std::endl;

		// Sample a random reversal.
		ReversalSampler sampler(genperm,false);
		sampler.debug      = true;
		ReversalRandom rev = sampler.sampleReversal(rng,true);
		//if(debug){printGenome(genperm.getExtendedPerm());}

		// Apply reversal.
		//genperm.debug = false;
		applyReversal(genperm, rev.g_beg, rev.g_end);
		genperm.clearBlockStatus();
		
		// Add reversal to list of reversals.
		reversals.emplace_back(rev);
		// Sample probability of stopping if genome is sorted.
		if(genperm.getBreakpoints() == 0){prob_rdm = distr(rng);}
	}
	if(debug){printGenome(genperm.getExtendedPerm());}
	return reversals;
}

std::vector<ReversalRandom> RandomReversalScenario::sampleModifiedScenario(std::mt19937& rng, int updpos_beg, int updpos_end){
	std::vector<ReversalRandom> reversals;
	return reversals;
}


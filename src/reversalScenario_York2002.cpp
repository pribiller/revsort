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
#include <cmath>	 // abs, tanh
#include <memory>	 // shared_ptr

#include "reversalScenario_York2002.hpp"

void RandomReversalScenario::printGenome(std::vector<int> perm){
	for (int const& gene : perm) {
		std::cout << gene << " ";
	}
	std::cout << std::endl;
}

std::vector<ReversalRandom> RandomReversalScenario::sampleScenario(GenomeMultichrom<int>& genome_B, std::mt19937& rng){

	// Unsigned extended permutation (it assumes that one of the permutations is the identity).
	GenomePermutation<BlockSimple> genperm(genome_B.getExtendedGenome());
	std::uniform_real_distribution distr(0.0, 1.0);
	int prob_rdm = distr(rng);

	std::vector<ReversalRandom> reversals;
	while((genperm.getBreakpoints() > 0) || (p_stop < prob_rdm)) { // not identity and p_stop.
		if(debug){std::cout << "\n\n\n------------------------------------------------\n\n - Current nb. of breakpoints = " << genperm.getBreakpoints() << std::endl;}

		// Sample a random reversal.
		ReversalSampler sampler(genperm,false);
		//sampler.debug      = true;
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

int RandomReversalScenario::samplePathLength(std::mt19937& rng, const int N, const float alpha, const float epsilon){
	// q(l) ~ 1 - tanh(epsilon*(l/(alpha*N)-1))
	// - alpha: lengths small than N*alpha are roughly equally likely represented.
	//          For example, alpha=0.65, means that there is more or less the same 
	//          chance to sample a size between 0 and 65% of the total size of the 
	//          path (N). From this point on, the probability drops to almost 0.

	// Compute CDF of q(l).
	std::vector<double> lprobs_cum(N,0.0);
	double l_total = 0;
	int l = 1;
	double prob_l = 1-std::tanh(epsilon*(l/(alpha*N)-1));
	lprobs_cum[0] = prob_l;	
	for(l=2; l<=N; ++l){ // l: 1..N; l_idx: 0..N-1
		prob_l = 1-std::tanh(epsilon*(l/(alpha*N)-1));
		lprobs_cum[l-1] = lprobs_cum[l-2] + prob_l;
		l_total += prob_l;
	}
	// Sample l.
	std::uniform_real_distribution distr(0.0, l_total); // [0, rev_weights_total)
	float rdmval = distr(rng);
	int l_chosen = 1;
	// Find the index corresponding to the random value
	for (int l_idx = 0; l_idx < N; ++l_idx) {
		if (rdmval < lprobs_cum[l_idx]) {
			l_chosen = (l_idx+1);
			break;
		}
	}
	return l_chosen;
}

int RandomReversalScenario::samplePathStart(std::mt19937& rng, const int N, const int l) {
	// Sample j (the start of a path of size l, given that total path size is N).
	std::uniform_int_distribution distr(0, N-l);
	return distr(rng);
}

std::pair<std::vector<int>,std::vector<int>> RandomReversalScenario::getPathEnds(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, const int pos_beg, const int pos_end) {
	// Unsigned extended permutation (it assumes that one of the permutations is the identity).
	GenomePermutation<BlockSimple> genperm(genome_B.getExtendedGenome());
	std::vector<int> perm_A;
	std::vector<int> perm_B;
	for(int rev_idx=0; rev_idx<(pos_end+1); ++rev_idx){
		// Before applying reversal, check if the current genome is one of the path ends.
		if(rev_idx == pos_beg){
			perm_A = genperm.getUnextendedPerm();
		}
		if(rev_idx == pos_end){
			perm_B = genperm.getUnextendedPerm();
			break;
		}
		// Apply reversal.
		ReversalRandom rev = reversals[rev_idx];
		applyReversal(genperm, rev.g_beg, rev.g_end);
	}
	return std::make_pair(perm_A, perm_B);
}

GenomeMultichrom<int> RandomReversalScenario::getGenomes(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, const int pos_beg, const int pos_end) {
	
	std::pair<std::vector<int>,std::vector<int>> pathEnds_perm = getPathEnds(genome_B, reversals, pos_beg, pos_end);

	std::pair<std::vector<int>,std::vector<bool>> genomeInfo_beg = permToGenome(pathEnds_perm.first);
	std::pair<std::vector<int>,std::vector<bool>> genomeInfo_end = permToGenome(pathEnds_perm.second); // new 'identity' genome

	GenomeMultichrom<int> genome_id(genomeInfo_end.first, genomeInfo_end.second);
	GenomeMultichrom<int> genome_other(genomeInfo_beg.first, genomeInfo_beg.second, genome_id.gene_labels_map);

	return genome_other;
}

int RandomReversalScenario::getGeneLabelDefault(GenomeMultichrom<int>& genome, int gene_lbl_ext) {

	// Check for caps or invalid (negative) genes.
	if((gene_lbl_ext < 2) || (gene_lbl_ext == genome.n+2)) {return gene_lbl_ext;}
	
	// gene_labels_map stores the mapping between the labels of the *unextended version*
	// of the new genome and the *unextended version* of the old genome 
	// In the *unextended version*, genes go between 1 and n, with gene i in index i-1.

	// On the other hand, labels in the list of reversals are represented in 
	// the *extended version*, with genes between 1 and n+2, where gene 1 and 
	// gene n+2 are "caps". For example, gene 2 in the extended version corresponds 
	// to gene 1 in the unextended version and so on.

	// Thus, to get the corresponding label of gene i in the extended version:
	// 1) Gene i in the new extended version corresponds to gene i-1 in the unextended version;
	// 2) Gene i-1 occupies the position i-2 in gene_labels_map (done in the function getLabel).
	// 3) Gene linked to mapped position i-2 is the default unextended gene label. It needs to be converted to extended version.
	// int gene_label_unext_default = std::stoi(genome.gene_labels_map.getLabelStr(gene_lbl_ext-1));
	// return (gene_label_unext_default < 0) ? (gene_label_unext_default-1) : (gene_label_unext_default+1);
	return genome.gene_labels_map.getLabel(gene_lbl_ext-1)+1;
}

std::vector<ReversalRandom> RandomReversalScenario::updateReversalScenario(GenomeMultichrom<int>& genome, std::vector<ReversalRandom> reversals, std::vector<ReversalRandom> reversals_new, const int pos_beg, const int pos_end) {
	// Make sure that new reversals and old reversals are based on the same gene labels.
	for (ReversalRandom& rev : reversals_new){
		//std::cout << " Reversal (" << rev.g_beg << ", " << rev.g_end << ") becomes ";
		rev.g_beg = getGeneLabelDefault(genome, rev.g_beg);
		rev.g_end = getGeneLabelDefault(genome, rev.g_end);
		rev.g_beg_next = getGeneLabelDefault(genome, rev.g_beg_next);
		rev.g_end_next = getGeneLabelDefault(genome, rev.g_end_next);
		//std::cout << " (" << rev.g_beg << ", " << rev.g_end << ")" << std::endl;
	}
	// Remove old path.
	reversals.erase(reversals.begin() + pos_beg, reversals.begin() + pos_end);
	// Add new path.
	reversals.insert(reversals.begin() + pos_beg, reversals_new.begin(), reversals_new.end());
	return reversals;
}

std::vector<ReversalRandom> RandomReversalScenario::sampleModifiedScenario(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, std::mt19937& rng){

	// Choosing a section of the path to replace.
	const int N = reversals.size();
	const int l = samplePathLength(rng, N);
	const int j = samplePathStart(rng, N, l);

	std::cout << " Total path size (nb. of reversals) = " << N << "; length modified path = " << l << "; Start modified path: " << j << std::endl;

	// Get genomes at these positions in the path.
	GenomeMultichrom<int> genome_B_new = getGenomes(genome_B, reversals, j, j+l);

	// Generate a new subpath.
	std::vector<ReversalRandom> reversals_new = sampleScenario(genome_B_new, rng);

	// // Integrate new reversals to the current path.
	std::vector<ReversalRandom> reversals_upd = updateReversalScenario(genome_B_new, reversals, reversals_new, j, j+l);

	std::cout << "\nOld scenario: " << std::endl;
	printSortingScenario(genome_B, reversals, j, j+l);

	std::cout << "\nNew scenario " << std::endl;
	printSortingScenario(genome_B, reversals_upd, j, j+l);

	// std::vector<ReversalRandom> reversals_upd;
	return reversals_upd;
}

void RandomReversalScenario::printSortingScenario(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, int pos_beg, int pos_end) {
	// Unsigned extended permutation (it assumes that one of the permutations is the identity).
	GenomePermutation<BlockSimple> genperm(genome_B.getExtendedGenome());

	std::cout << " Genome [start]" << std::endl;
	std::vector<int> genome =  genperm.getExtendedPerm();
	for (int gene : genome) {std::cout << gene << " ";}
	// Check if the current genome is one of the path ends.
	if(0 == pos_beg){std::cout << " [A] ";}
	if(0 == pos_end){std::cout << " [B] ";}
	std::cout << std::endl;

	for(int rev_idx=0; rev_idx<reversals.size(); ++rev_idx){		
		
		// Apply reversal.
		ReversalRandom rev = reversals[rev_idx];
		applyReversal(genperm, rev.g_beg, rev.g_end);

		// Print genome after reversal.
		std::cout << " Genome [" << (rev_idx+1) << "] after reversal (" << rev.g_beg << ", " << rev.g_end << "]: " << std::endl;
		std::vector<int> genome =  genperm.getExtendedPerm();
		for (int gene : genome) {std::cout << gene << " ";}
		// Check if the current genome is one of the path ends.
		if((rev_idx-1) == pos_beg){std::cout << " [A] ";}
		if((rev_idx-1) == pos_end){std::cout << " [B] ";}
		std::cout << std::endl;

	}
}

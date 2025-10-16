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


//////////////////////////////////////////////////////////////////////
// MCMC Metropolis-Hastings
//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
// Mean number of reversals update.
//////////////////////////////////////////////////////////////////////

double ProposalReversalMean::sampleReversalMean(std::mt19937& rng, const double rev_mean_cur) const{
	const double rev_mean_lb = std::max(rev_mean_min, rev_mean_cur-rev_mean_range/2);
	const double rev_mean_ub = rev_mean_lb + rev_mean_range;

	// Proposed values are uniformly sampled from the range [cur-range, cur+range].
	std::uniform_real_distribution<double> distrib(rev_mean_lb, rev_mean_ub);
	return distrib(rng);
}

double ProposalReversalMean::getAcceptanceProb(const int rev_path_len, const double rev_mean_cur, const double rev_mean_prop) const{
	// P(lambda | X, D) ~ P(X|lambda) P(lambda) ~ e^(-lambda) * lambda^L_x * P(lambda)
	// P(lambda) = 1/(lambda_max-lambda_min+1)
	return std::exp(-rev_mean_prop+rev_mean_cur)*(std::pow(rev_mean_prop/rev_mean_cur, rev_path_len));
}

//////////////////////////////////////////////////////////////////////
// Inversion history update.
//////////////////////////////////////////////////////////////////////

// Probability of proposing a particular path with l inversions.
// The probability of a path with l inversions will be the product 
// of l+1 factors, one for each inversion and a factor of 1-q_stop
// for actually stopping upon reaching the target genome 
// [from York, Durrett, and Nielsen (2002)].
std::vector<double> ProposalReversalScenario::q_path_factors(const std::vector<ReversalRandom>& path, const std::vector<float>& rev_weights, const double p_stop){
	std::vector<double> factors;
	factors.reserve(path.size());
	for (ReversalRandom rev : path){
		// Compute probability of choosing a particular inversion.
		// Reversal types without any reversals associated with have weight = 0.
		double w_total = 0.0;
		for (int revtype_idx = 0; revtype_idx < ReversalType_COUNT; ++revtype_idx) {
			if(rev.rev_totals[revtype_idx] > 0){w_total += rev_weights[revtype_idx];}
		}
		const double w_revtype = rev_weights[rev.type];

		// Number of inversions of a certain type.
		// For now all inversions have the same chance to be sampled. 
		// For example, there is no bias depending on the reversal location or its size.
		const int N_revtype = rev.rev_totals[rev.type];

		factors.emplace_back((w_total/w_revtype)*(N_revtype/1.0));

		std::cout << " - Factors : Type=" << rev.type << "; W=" << w_total << "; N_revtype=" << N_revtype  << "; prod=" << (w_total/w_revtype)*(N_revtype/1.0) << std::endl;

	}
	return factors;
}

// Each proposal ratio (also called a Hastings ratio) is the probability of
// proposing the original state from proposed state divided by
// the probability of proposing the proposed state from the
// original state.
double ProposalReversalScenario::getProposalRatio(const std::vector<float>& rev_weights, const double p_stop){

	// Proposal probability [from York, Durrett, and Nielsen (2002)].
	// q(Y|X) = q_L(l,j) * q_new
	const double q_lj_cur = q_lj(L_cur, l_cur, path_beg);

	// Probability of getting current state (X) from the new state (Y).
	// q(X|Y) = q_L'(l',j) * q_old
	const double q_lj_new =q_lj(L_new, l_new, path_beg);

	// Looks only at the subpath that was modified by the proposal. 
	// The factors from the common parts cancel out.
	std::vector<ReversalRandom> p_cur( currentReversalScenario.begin() + path_beg, currentReversalScenario.begin()  + (path_beg + l_cur));
	std::vector<ReversalRandom> p_new(proposedReversalScenario.begin() + path_beg, proposedReversalScenario.begin() + (path_beg + l_new));

	std::cout << "\n\n - Current Factors " << std::endl;
	std::vector<double> q_cur = q_path_factors(p_cur, rev_weights, p_stop);

	std::cout << "\n\n - Proposed Factors " << std::endl;
	std::vector<double> q_new = q_path_factors(p_new, rev_weights, p_stop);
	std::cout << "\n\n";
	
	std::sort(q_cur.begin(), q_cur.end(), std::greater<double>());
	std::sort(q_new.begin(), q_new.end(), std::greater<double>());

	// Compute the proposal ratio (Hastings ratio).
	// Proposal ratio: q(X|Y) / q(Y|X)
	// q(X|Y) = q_L'(l',j) * q_old
	// q(Y|X) = q_L(l,j)   * q_new
	double proposalRatio = q_lj_new/q_lj_cur; // q_L'(l',j) / q_L(l,j)
	std::cout << " - q_L'(l',j)=" << q_lj_new << "; q_L(l,j)=" << q_lj_cur << "; ratio=" << proposalRatio << std::endl;
	int factors_max = std::max(q_cur.size(), q_new.size());
	int factors_idx = 0;
	while(factors_idx < factors_max){
		// Instead of doing (1/q_old) / (1/q_new) = q_new / q_old
		const double num = (factors_idx < q_new.size()) ? q_new[factors_idx] : 1.0;
		const double den = (factors_idx < q_cur.size()) ? q_cur[factors_idx] : 1.0;
		proposalRatio *= (num/den);
		std::cout << " - Factors : cur=" << q_cur[factors_idx] << "; prop=" << q_new[factors_idx] << std::endl;
		++factors_idx;
	}
	std::cout << " - Proposal ratio : " << proposalRatio << std::endl;
	return proposalRatio;
}

// The posterior ratio is the ratio of the posterior probability
// of the proposed state over that of the current state. 
double ProposalReversalScenario::getPosteriorRatio(const double rev_mean){
	
	// Simple case: both paths have the same size.
	// For now there is no bias towards small reversals or specific regions for example.
	if (L_cur == L_new) {return 1.0;}
	
	// Number of possible reversals that can be done in a certain step (including all types: good, bad, neutral).
	const int nb_rev_step = N*(N+1)/2;
	
	// Compute the term [lambda^(-L_cur + L_new) * L_cur! / L_new!]
	// avoiding the calculation of factorials which can become large.
	const int L_ub = std::max(L_cur, L_new);
	const int L_lb = std::min(L_cur, L_new);
	double factorial = 1.0;
	for(int factor = L_lb; factor < L_ub; ++factor){
		if(L_cur > L_new) {
			factorial *= ((factor+1.0)/rev_mean);
		} else {
			factorial *= (rev_mean/(factor+1.0));
		}
	}

	// P(X_new|lambda) / P(X_cur|lambda)
	// P(X|lambda) is defined in York, Durrett, and Nielsen (2002)
	return std::pow(nb_rev_step, -L_new+L_cur) * factorial;
}

// The acceptance probability is the minimum of one and
// the product of the posterior ratio and the proposal ratio.
double ProposalReversalScenario::getAcceptanceProb(const double rev_mean, const std::vector<float>& rev_weights, const double p_stop){
	return std::min(1.0, getProposalRatio(rev_weights, p_stop)*getPosteriorRatio(rev_mean));
}

//////////////////////////////////////////////////////////////////////
// Proposal distribution
//////////////////////////////////////////////////////////////////////

void RandomReversalScenario::printGenome(std::vector<int> perm){
	for (int const& gene : perm) {
		std::cout << gene << " ";
	}
	std::cout << std::endl;
}


int RandomReversalScenario::getGeneExtremity(const int gene_id, const GenomePermutation<BlockSimple>& genperm) const {
	std::pair<int,int> gene_exts = genperm.geneToUnsExt(gene_id);
	return (genperm.isReversed(gene_id)) ? gene_exts.first : gene_exts.second;
}

std::vector<ReversalRandom> RandomReversalScenario::getSubpath(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, const int pos_beg, const int pos_end){

	// Unsigned extended permutation (it assumes that one of the permutations is the identity).
	GenomePermutation<BlockSimple> genperm(genome_B.getExtendedGenome());

	std::vector<ReversalRandom> rev_subpath(reversals.begin() + pos_beg, reversals.begin() + pos_end);
	for(ReversalRandom& rev : rev_subpath){

		ReversalRandom rev_ori(rev);

		// Convert gene labels to the gene labels used in the subpath.
		//ReversalRandom rev; 
		rev.g_beg = getGeneLabelSubpath(genperm, genome_B, rev.g_beg);;
		rev.g_end = getGeneLabelSubpath(genperm, genome_B, rev.g_end);
		rev.g_beg_next = getGeneLabelSubpath(genperm, genome_B, rev.g_beg_next);
		rev.g_end_next = getGeneLabelSubpath(genperm, genome_B, rev.g_end_next);

		// Get the type of reversal and probability in the subpath.
		std::pair<int,int> rev_extremities = std::make_pair(getGeneExtremity(rev.g_beg, genperm),getGeneExtremity(rev.g_end, genperm));
		ReversalSampler sampler(genperm,debug);
		sampler.updateComponents();
		sampler.countReversals();

		rev.rev_totals = sampler.rev_totals;
		rev.type  = sampler.getReversalType(rev_extremities); // ReversalType type;
		rev.cycle = sampler.getCycleType(rev_extremities); // CycleType cycle;

		// std::cout << "- ORI.rev.= (" << rev_ori.g_beg << ", " << rev_ori.g_end << "), type=" << rev_ori.type << ", cyctype=" << rev_ori.cycle << " // ";
		// std::cout << "MOD.rev.= (" << rev.g_beg << ", " << rev.g_end << "), type=" << rev.type << ", cyctype=" << rev.cycle << ", exts=" << rev_extremities.first << ", " << rev_extremities.second << ";" << std::endl;
		// for (int revtype_idx = 0; revtype_idx < ReversalType_COUNT; ++revtype_idx) {
		// 	std::cout << "\t type=" << revtype_idx <<  "; count=" << rev.rev_totals[revtype_idx] << std::endl;
		// }

		// Apply reversal.
		genperm.debug = false;
		applyReversal(genperm, rev.g_beg, rev.g_end);
		genperm.clearBlockStatus();

	}
	return rev_subpath;
}


std::vector<ReversalRandom> RandomReversalScenario::sampleScenario(GenomeMultichrom<int>& genome_B, std::mt19937& rng){

	// Unsigned extended permutation (it assumes that one of the permutations is the identity).
	GenomePermutation<BlockSimple> genperm(genome_B.getExtendedGenome());
	std::uniform_real_distribution distr(0.0, 1.0);
	double prob_rdm = distr(rng);

	std::vector<ReversalRandom> reversals;
	while((genperm.getBreakpoints() > 0) || (p_stop < prob_rdm)) { // not identity and p_stop.
		if(debug){std::cout << "\n\n\n------------------------------------------------\n\n - Current nb. of breakpoints = " << genperm.getBreakpoints() << std::endl;}

		// Sample a random reversal.
		ReversalSampler sampler(genperm,debug);
		//sampler.debug      = true;
		ReversalRandom rev = sampler.sampleReversal(rng,true);
		//if(debug){printGenome(genperm.getExtendedPerm());}

		// std::cout << "- NEW.rev.= (" << rev.g_beg << ", " << rev.g_end << "), type=" << rev.type << ", cyctype=" << rev.cycle << ";" << std::endl;
		// for (int revtype_idx = 0; revtype_idx < ReversalType_COUNT; ++revtype_idx) {
		// 	std::cout << "\t type=" << revtype_idx <<  "; count=" << rev.rev_totals[revtype_idx] << std::endl;
		// }

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
	double prob_l = 1.0-std::tanh(epsilon*(l/(alpha*N)-1));
	lprobs_cum[0] = prob_l;	
	for(l=2; l<=N; ++l){ // l: 1..N; l_idx: 0..N-1
		prob_l = 1.0-std::tanh(epsilon*(l/(alpha*N)-1));
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

int RandomReversalScenario::getGeneLabelSubpath(GenomePermutation<BlockSimple>& genperm, GenomeMultichrom<int>& genome, int gene_lbl_ext) {
	// Check for caps or invalid (negative) genes.
	if((gene_lbl_ext < 2) || (gene_lbl_ext == genome.n+2)) {
		return gene_lbl_ext;
	// Return extended label in the new representation (see comments in the 
	// function ``getGeneLabelDefault''; it explains how the mapping works).
	} else {
		return std::abs(genome.gene_labels_map.getId(gene_lbl_ext-1))+1;	
	}
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

ProposalReversalScenario RandomReversalScenario::sampleModifiedScenario(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, std::mt19937& rng){

	// Choosing a section of the path to replace.
	const int N = reversals.size();
	const int l = samplePathLength(rng, N);
	const int j = samplePathStart(rng, N, l);

	std::cout << " Total path size (nb. of reversals) = " << N << "; length modified path = " << l << "; Start modified path: " << j << std::endl;

	// Get genomes at these positions in the path.
	GenomeMultichrom<int> genome_B_new = getGenomes(genome_B, reversals, j, j+l);

	// Get current subpath.
	//debug = false;
	std::vector<ReversalRandom> reversals_sub_cur = getSubpath(genome_B_new, reversals, j, j+l);

	// Integrate probabilities of the subpath to the current path.
	std::vector<ReversalRandom> reversals_cur = updateReversalScenario(genome_B_new, reversals, reversals_sub_cur, j, j+l);

	// Generate a new subpath.
	std::vector<ReversalRandom> reversals_sub_new = sampleScenario(genome_B_new, rng);
	const int l_new = reversals_sub_new.size();

	// Integrate new reversals to the current path.
	std::vector<ReversalRandom> reversals_new = updateReversalScenario(genome_B_new, reversals, reversals_sub_new, j, j+l);
	const int N_new = reversals_new.size();

	std::cout << "\nOld scenario: " << std::endl;
	printSortingScenario(genome_B, reversals, j, j+l);

	std::cout << "\nNew scenario " << std::endl;
	printSortingScenario(genome_B, reversals_new, j, j+l_new);

	// std::vector<ReversalRandom> reversals_new;
	return ProposalReversalScenario(reversals_cur, reversals_new, j, N, l, N_new, l_new, genome_B.n);
}

void RandomReversalScenario::printSortingScenario(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, int pos_beg, int pos_end) {
	// Unsigned extended permutation (it assumes that one of the permutations is the identity).
	GenomePermutation<BlockSimple> genperm(genome_B.getExtendedGenome());

	std::cout << " Genome [start]" << std::endl;
	std::vector<int> genome =  genperm.getExtendedPerm();
	for (int gene : genome) {std::cout << gene << " ";}
	// Check if the current genome is one of the path ends.
	if(0 == pos_beg){std::cout << " [A; pos_beg=" << pos_beg << "] ";}
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
		if((rev_idx+1) == pos_beg){std::cout << " [A; pos_beg=" << pos_beg << "] ";}
		if((rev_idx+1) == pos_end){std::cout << " [B] ";}
		std::cout << std::endl;

	}
}

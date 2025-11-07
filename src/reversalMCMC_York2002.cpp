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
 * in 3 types: "good reversals", "neutral reversals", and 
 * "bad reversals". These types are defined based on how 
 * close the genome gets to the target genome after they are applied.
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
 * similar way to the method of York et al. (2002) [also similar to Miklos (2003)].
 * 
 * 3) The parallel tempering is similar to Larget et al. (2005), but a 
 * clearer explanation is provided in Altekar (2004) (Mr.Bayes paper).
 * 
 *******************************************************/

#include <iostream>
#include <unordered_map>
#include <vector>
#include <stack>
#include <list>
#include <cstdlib> // exit
#include <cmath>   // abs, tanh, log
#include <memory>  // shared_ptr
#include <iomanip> // put_time
#include <ctime>
#include <sstream>
#include <boost/lexical_cast.hpp> // lexical_cast<string>

#include "reversalMCMC_York2002.hpp"


//////////////////////////////////////////////////////////////////////
// MCMC Metropolis-Hastings
//////////////////////////////////////////////////////////////////////

void ReversalMCMC::initializeIdRun(){
	std::time_t t = std::time(nullptr); 
	std::tm tm = *std::localtime(&t);
	std::ostringstream oss;
	oss << std::put_time(&tm, "%Y%m%d%H%M%S");
	id_run = "revMCMC_" + oss.str();
}

void ReversalMCMC::initializeChains(){
	
	std::cout << "Initializing chains..." << std::endl;
	cur_step = 1;
	
	// Varying the probability of sampling neutral reversals allows 
	// the sampling of shorter and longer paths.
	// - p_neutral used in initial paths in York, Durrett, and Nielsen:
	// [0.025 (short initial inversions paths), 0.7 (long initial inversion paths)]
	RandomReversalScenario sampler_std(rev_weights,p_stop,false);
	const double p_neutral     = sampler_std.rev_weights[ReversalType::GOOD]/2.0;
	const double p_neutral_min = sampler_std.rev_weights[ReversalType::BAD];
	const double p_delta = (nb_chains > 1) ? (p_neutral-p_neutral_min)/(nb_chains-1) : 0.0;

	#pragma omp parallel for collapse(2)
	for (int chain_idx=0; chain_idx<nb_chains; ++chain_idx){
		for (int temp_idx=0; temp_idx<nb_temperatures; ++temp_idx){

			const int idx     = temp_idx*nb_chains + chain_idx;
			const double temp = 1.0/(1.0+temp_idx*delta_temp);

			RandomReversalScenario sampler(rev_weights,p_stop,false);
			sampler.rev_weights[ReversalType::NEUTRAL] = p_neutral-chain_idx*p_delta;
			currentState_revHists[idx] = sampler.sampleScenario(genome_B,rng);
			currentState_revMeans[idx] = proposalMean.sampleReversalMean(currentState_revHists[idx].size());
			rev_path_avgsize[idx]      = currentState_revHists[idx].size();
			std::cout << "- Chain " << idx << ": p_neutral=" << sampler.rev_weights[ReversalType::NEUTRAL] << "; size=" << currentState_revHists[idx].size() << "; size_min=" << rev_dist << std::endl;
		}
	}
}

std::string ReversalMCMC::runSingleChain(const int chainIdx, const double chainTemp){

	std::vector<ReversalRandom>& reversals = currentState_revHists[chainIdx];
	double rev_mean = currentState_revMeans[chainIdx];
	// Sample a modified reversal history.
	if(debug){std::cout << "\nSampling modified scenario (chainTemp=" << chainTemp << ")..." << std::endl;}
	RandomReversalScenario sampler(rev_weights, p_stop, debug);
	ProposalReversalScenario proposalHist = sampler.sampleModifiedScenario(genome_B, reversals, rng);

	std::string status =  "  L_cur=" + std::to_string(proposalHist.L_cur)
						+ "  l_cur=" + std::to_string(proposalHist.l_cur)
						+ "  L_new=" + std::to_string(proposalHist.L_new)
						+ "  l_new=" + std::to_string(proposalHist.l_new);

	// Compute the acceptance probability of the modified reversal history.
	if(debug){std::cout << "\nComputing acceptance probability..." << std::endl;}
	double acceptanceProb = proposalHist.getAcceptanceProb(rev_mean, sampler.rev_weights, sampler.p_stop, chainTemp);
	if(debug){std::cout << "Acceptance prob. = " << acceptanceProb << "; Posterior ratio = " << proposalHist.posteriorRatio << "; Proposal ratio = " << proposalHist.proposalRatio << std::endl;}

	status = status + "; Rev.Types: Cur= ";
	for (int revtype_idx = 0; revtype_idx < ReversalType_COUNT; ++revtype_idx) {
		status = status + std::to_string(revtype_idx) + ": " + std::to_string(proposalHist.P_cur_revtypes[revtype_idx]) + ", ";
	}
	status = status + "New= ";
	for (int revtype_idx = 0; revtype_idx < ReversalType_COUNT; ++revtype_idx) {
		status = status + std::to_string(revtype_idx) + ": " + std::to_string(proposalHist.P_new_revtypes[revtype_idx]) + ", ";
	}

	// Check if current reversal history should be updated or not based on the acceptance probability.
	const double val = distr_accept(rng);
	if(val <= acceptanceProb){
		if(debug){std::cout << "\nNew scenario: ACCEPTED" << std::endl;}
		currentState_revHists[chainIdx] = proposalHist.proposedReversalScenario;
		status = status + " [ACCEPTED] : " + std::to_string(currentState_revHists[chainIdx].size()) + "; rdm=" + std::to_string(val) + "; accepProb=" + boost::lexical_cast<std::string>(acceptanceProb) + "; postRatio=" + boost::lexical_cast<std::string>(proposalHist.posteriorRatio) + "; propRatio=" + boost::lexical_cast<std::string>(proposalHist.proposalRatio);
	} else {
		if(debug){std::cout << "\nNew scenario: REJECTED" << std::endl;}
		status = status + " [REJECTED] : " + std::to_string(currentState_revHists[chainIdx].size()) + "; rdm=" + std::to_string(val) + "; accepProb=" + boost::lexical_cast<std::string>(acceptanceProb) + "; postRatio=" + boost::lexical_cast<std::string>(proposalHist.posteriorRatio) + "; propRatio=" + boost::lexical_cast<std::string>(proposalHist.proposalRatio);
	}

	// Sample an updated reversal mean.
	proposalMean.proposeNewValue(currentState_revHists[chainIdx].size(), currentState_revMeans[chainIdx]);
	if(debug){std::cout << "\nCur. mean= " << currentState_revMeans[chainIdx] << "; New mean=" << proposalMean.rev_mean_new << "; Acceptance prob.=" << proposalMean.acceptanceProb << std::endl;}
	if(distr_accept(rng) < proposalMean.acceptanceProb){
		if(debug){std::cout << "\nNew mean (λ=" << proposalMean.rev_mean_new << "): ACCEPTED..." << std::endl;}
		currentState_revMeans[chainIdx] = proposalMean.rev_mean_new;
	} else {
		if(debug){std::cout << "\nNew mean (λ=" << proposalMean.rev_mean_new << "): REJECTED..." << std::endl;}
	}
	return status;
}

double ReversalMCMC::computeWithinChainVariance(){
	// Within-chain variance:
	// W = 1/m * [ sum_j (1/(n-1)) * (L_ij - L_j)^2 ]
	double W = 0.0; 
	for (int idx=0; idx<nb_chains; ++idx){
		const double L_j = rev_path_avgsize[idx]; // Avg. path size in chain j.
		for (auto L_hist : hist_chains[idx]){
			const int L_ij = L_hist.first;  // Size of the reversal path.
			const int freq = L_hist.second; // Frequency that a paths of a certain size were visited.
			W += (freq*std::pow((L_ij-L_j), 2.0));
		}
	}
	return W/(nb_chains*(cur_step-1.0));
}

double ReversalMCMC::computeBetweenChainVariance(){
	// Average path size of all chains at all steps.
	double L = 1.0 * std::accumulate(rev_path_avgsize.begin(), rev_path_avgsize.end(), 0LL) / nb_chains;
	// std::cout << "\t Global average so far: L=" << L << std::endl;

	// Between chain variance:
	// B = 1/(m-1) * sum_j [L_j - L]^2
	double B = 0.0;
	for (int idx=0; idx<nb_chains; ++idx){
		const double L_j = rev_path_avgsize[idx]; // Avg. path size in chain j.
		// std::cout << "\t Average chain " <<  idx << ": L_j=" << L_j << std::endl;
		B += std::pow((L_j-L), 2.0);
	}
	return B / (nb_chains-1.0);
}

void ReversalMCMC::run(){
	
	// TODO: Change directory for the run.
	const std::string bkp_filename = id_run + ".bkp";

	bool sample_state     = false;
	int last_sampled_step = 0;
	std::uniform_int_distribution distr_sample(0, nb_chains-1);
	const int N = genome_B.n;

	while((cur_step < max_steps) && (sampled_revMeans.size() < sample_amount)){
		// Update current step.
		++cur_step;
		std::vector<std::string> status_all(nb_chains*nb_temperatures);

		// std::cout << "Running step " << cur_step << std::endl;

		// Run current step for all chains in parallel.
		// std::cout << "Updating chains..." << std::endl;
		#pragma omp parallel for collapse(2)
		for (int chain_idx=0; chain_idx<nb_chains; ++chain_idx){
			for (int temp_idx=0; temp_idx<nb_temperatures; ++temp_idx){

				const int idx     = temp_idx*nb_chains + chain_idx;
				const double temp = 1.0/(1.0+temp_idx*delta_temp);
				status_all[idx]   = runSingleChain(idx, temp);
				
				// Update average path size for each chain (L_j).
				if(cur_step > pre_burnin_steps) {
					const int steps_after_pre_burnin = cur_step-pre_burnin_steps;
					// Initial step.
					if(steps_after_pre_burnin == 1){
						rev_path_avgsize[idx] = currentState_revHists[idx].size();
					// General step.
					} else {
						rev_path_avgsize[idx] = ((steps_after_pre_burnin-1.0)/steps_after_pre_burnin)*rev_path_avgsize[idx] + static_cast<double>(currentState_revHists[idx].size())/steps_after_pre_burnin;
					}
				}
			}
		}
		
		// Swap chains.
		if((int(cur_step) % print_interval) == 0){
			std::cout << "\nStep " << cur_step << " - Swapping chains" << std::endl;
		} 
		// std::cout << "Swapping chains..." << std::endl;
		ProposalReversalScenario propSwap;
		for (int temp_idx=0; temp_idx<(nb_temperatures-1); ++temp_idx){
			for (int chain_idx=0; chain_idx<nb_chains; ++chain_idx){

				const int idx_cur  = temp_idx*nb_chains     + chain_idx;
				const int idx_heat = (temp_idx+1)*nb_chains + chain_idx;
				
				const double temp_cur  = 1.0/(1.0+temp_idx*delta_temp);
				const double temp_heat = 1.0/(1.0+(temp_idx+1)*delta_temp);

				const int L_cur  = currentState_revHists[idx_cur].size();
				const int L_heat = currentState_revHists[idx_heat].size();

				double rev_mean_cur  = currentState_revMeans[idx_cur];
				double rev_mean_heat = currentState_revMeans[idx_heat];
				
				// Check if states of chain i at temperature j and j+1 should be swaped.
				// Swap of states is accepted with probability: min(1,[(post_i^temp_j)*(post_j^temp_i)]/[(post_j^temp_j)*(post_i^temp_i)])
				const double acceptanceProb_swap = std::min(1.0, std::pow(propSwap.getPosteriorRatio(rev_mean_cur, L_cur, L_heat, N),temp_cur)*std::pow(propSwap.getPosteriorRatio(rev_mean_heat, L_heat, L_cur, N),temp_heat));

				const double val = distr_accept(rng);
				if(val <= acceptanceProb_swap){
					if((int(cur_step) % print_interval) == 0){
						std::cout << " [ACCEPTED] : Chains w/temps " << temp_idx << " (L=" << L_cur<< ") and " << (temp_idx+1) << " (L=" << L_heat << ") swaped; rdm=" << val << "; accepProb=" << boost::lexical_cast<std::string>(acceptanceProb_swap) << std::endl;
					}
					std::swap(currentState_revHists[idx_cur], currentState_revHists[idx_heat]);
					std::swap(currentState_revMeans[idx_cur], currentState_revMeans[idx_heat]);					
				} else {
					if((int(cur_step) % print_interval) == 0){
						std::cout << " [REJECTED] : Chains w/temps " << temp_idx << " (L=" << L_cur<< ") and " << (temp_idx+1) << " (L=" << L_heat << ") not swaped; rdm=" << val << "; accepProb=" << boost::lexical_cast<std::string>(acceptanceProb_swap) << std::endl;
					}
				}
			}
		}

		// Compute convergence measures if pre burn-in phase is over.
		if((cur_step > pre_burnin_steps) && (nb_chains > 1)){
			const int steps_after_pre_burnin = cur_step-pre_burnin_steps;
			// It stores the history of a chain (visited states so far).
			for (int idx=0; idx<nb_chains; ++idx){
				const int curSize = currentState_revHists[idx].size();
				hist_chains[idx][curSize] += 1;
			}
			// Measures used to assess convergence [York, Durrett, Nielsen (2002)].
			// Compute between chain variance and within-chain variance.
			W = computeWithinChainVariance();
			B = computeBetweenChainVariance();
			// Method of Gelman and Rubin (1992).
			R = std::sqrt((steps_after_pre_burnin-1.0)/steps_after_pre_burnin + B/W);
		}

		// Sample one of the current states.
		sample_state = (sample_state || ((check_convergence && (R <= 1.1)) || (!check_convergence && (cur_step > pre_burnin_steps))));
		if (sample_state){
			if((last_sampled_step % sample_interval) == 0){
				const int rdmidx = distr_sample(rng);
				sampled_revHists.emplace_back(currentState_revHists[rdmidx]);
				sampled_revMeans.emplace_back(currentState_revMeans[rdmidx]);
				// for(ReversalRandom rev : currentState_revHists[rdmidx]){
				// 	std::cout << "Reversal (" << rev.g_beg << ", " << rev.g_end << ")" << std::endl;
				// }
			}
			++last_sampled_step;
		}

		// Save a checkpoint, in case the run crashes later.
		if((int(cur_step) % backup_interval) == 0){saveState(bkp_filename);}

		// Print stats about the chains and their convergence.
		if((int(cur_step) % print_interval) == 0){
			std::cout << "\nStep " << cur_step << std::endl;
			// Print sorted average path sizes (one value per chain).
			std::vector<std::pair<int,double>> path_averages(nb_chains);
			for (int idx=0; idx<nb_chains; ++idx){
				path_averages[idx] = std::make_pair(idx, rev_path_avgsize[idx]);
			}
			// Sort the vector by descending order.
			std::sort(path_averages.begin(), path_averages.end(), [](const auto& a, const auto& b) {return a.second > b.second;});
			// Print values.
			double avg_revMean_cur = 0.0;
			double avg_revPath_cur = 0.0;
			for (const auto& p : path_averages) {
				const int idx = p.first;
				avg_revMean_cur += currentState_revMeans[idx];
				avg_revPath_cur += currentState_revHists[idx].size();
				std::cout << " [Chain " << idx << "] Avg. path size: " << rev_path_avgsize[idx] << "; cur path = " << currentState_revHists[idx].size() << " ; cur lambda = " << currentState_revMeans[idx] << " / " << status_all[idx] << std::endl;
			}
			std::cout << " Avg. lambda = " << (avg_revMean_cur/nb_chains) << "; Avg. path size=" << (avg_revPath_cur/nb_chains) << std::endl;
			if ((cur_step > pre_burnin_steps) && (nb_chains > 1)) {
				std::cout << "  R=" << R << "; W=" << W << "; B=" << B << std::endl;
			} else {
				std::cout << "  R=N/A; W=N/A; B=N/A " << ((nb_chains > 1) ? "[Pre burn in phase]" : "[Available only if nb. chains > 1]") << std::endl;
			}
		}
	}
	// Save the last state.
	saveState(bkp_filename);
}

//////////////////////////////////////////////////////////////////////
// Mean number of reversals update.
//////////////////////////////////////////////////////////////////////

double ProposalReversalMean::sampleReversalMean(const double rev_path_len_cur) {
	// Proposed values are sampled from the range [rev_mean_min, rev_mean_max].
	std::vector<double> rev_probs_cdf(rev_mean_max-rev_mean_min);
	double rev_probs_sum=0.0;
	for(int rev_samp=rev_mean_min; rev_samp<rev_mean_max; ++rev_samp){
		const double prob_rev_samp = -rev_samp + rev_path_len_cur*std::log(rev_samp);
		rev_probs_sum += prob_rev_samp;
		rev_probs_cdf[rev_samp-rev_mean_min] = rev_probs_sum;
	}	
	for(int rev_samp=rev_mean_min; rev_samp<rev_mean_max; ++rev_samp){
		rev_probs_cdf[rev_samp-rev_mean_min] = rev_probs_cdf[rev_samp-rev_mean_min]/rev_probs_sum;
	}
	std::uniform_real_distribution<double> distr(0.0, 1.0);
	double rdmval = distr(rng);
	for(int rev_samp=rev_mean_min; rev_samp<rev_mean_max; ++rev_samp){
		if (rev_probs_cdf[rev_samp-rev_mean_min] > rdmval){
			rev_mean_new = rev_samp;
			break;
		}

	}
	// Proposed values are uniformly sampled from the range [cur-range, cur+range].
	// std::uniform_real_distribution<double> distrib(rev_mean_lb, rev_mean_ub);
	//const double rev_mean_lb = rev_mean_min;
	//const double rev_mean_ub = rev_mean_min*2.0;
	// rev_mean_new = distrib(rng);
	return rev_mean_new;
}

double ProposalReversalMean::getAcceptanceProb(const int rev_path_len, const double rev_mean_cur, const double rev_mean_prop) {
	// P(lambda | X, D) ~ P(X|lambda) P(lambda) ~ e^(-lambda) * lambda^L_x * P(lambda)
	// P(lambda) = 1/(lambda_max-lambda_min+1)
	acceptanceProb = std::exp(-rev_mean_prop+rev_mean_cur)*(std::pow(rev_mean_prop/rev_mean_cur, rev_path_len));
	return acceptanceProb;
}

//////////////////////////////////////////////////////////////////////
// Inversion history update.
//////////////////////////////////////////////////////////////////////

// Probability of proposing a particular path with l inversions.
// The probability of a path with l inversions will be the product 
// of l+1 factors, one for each inversion and a factor of 1-q_stop
// for actually stopping upon reaching the target genome 
// [from York, Durrett, and Nielsen (2002)].
std::pair<std::vector<double>,std::vector<int>> ProposalReversalScenario::q_path_factors(const std::vector<ReversalRandom>& path, const std::vector<double>& rev_weights, const double p_stop){
	// std::cout << "Factors: " << std::endl;

	std::vector<double> factors;
	std::vector<int> path_rev_types(ReversalType_COUNT);
	factors.reserve(path.size());
	for (ReversalRandom rev : path){
		// Compute probability of choosing a particular inversion.
		// Reversal types without any reversals associated with have weight = 0.
		double w_total = 0.0;
		for (int revtype_idx = 0; revtype_idx < ReversalType_COUNT; ++revtype_idx) {
			if(rev.rev_totals[revtype_idx] > 0){w_total += rev_weights[revtype_idx];}
			//w_total += (rev.rev_totals[revtype_idx]*rev_weights[revtype_idx]);
		}
		const double w_revtype = rev_weights[rev.type];
		//const double w_revtype = rev.rev_totals[rev.type]*rev_weights[rev.type];
		path_rev_types[rev.type] += 1;

		// Number of inversions of a certain type.
		// For now all inversions have the same chance to be sampled. 
		// For example, there is no bias depending on the reversal location or its size.
		const int N_revtype = rev.rev_totals[rev.type];

		// std::cout << "\t Type: " << rev.type << "; Nb: " << N_revtype << "; w= " << w_revtype << " /" << w_total  << "(" << w_total/w_revtype << ")" << std::endl;
		if(w_revtype > 0.0){
			factors.emplace_back((w_total/w_revtype)*(N_revtype/1.0));
			// factors.emplace_back(w_total/w_revtype);

		// It is a reversal that is not allowed in the walk (e.g. a bad reversal), but it occurred in the path.
		} else {
			for (int revtype_idx = 0; revtype_idx < ReversalType_COUNT; ++revtype_idx) {
				std::cout << " type = " << revtype_idx << "; w=" << rev_weights[revtype_idx] << std::endl;
			}
			std::cout << "ERROR! There is a state in the chain with a reversal with weight = 0. Due to this, the state becomes unreachable. The program is being aborted.";
			exit(1);
		}

		if(debug){std::cout << " - Factors : Type=" << rev.type << "; W=" << w_total << "; W_type=" << w_revtype << "; N_revtype=" << N_revtype  << "; prod=" << (w_total/w_revtype) << std::endl;}
	}
	return std::make_pair(factors,path_rev_types);
}

// Each proposal ratio (also called a Hastings ratio) is the probability of
// proposing the original state from proposed state divided by
// the probability of proposing the proposed state from the
// original state.
double ProposalReversalScenario::getProposalRatio(const std::vector<double>& rev_weights, const double p_stop){
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

	if(debug){std::cout << "\n\n - Current Factors " << std::endl;}
	std::pair<std::vector<double>,std::vector<int>> factors_cur = q_path_factors(p_cur, rev_weights, p_stop);
	std::vector<double> q_cur = factors_cur.first;
	P_cur_revtypes = factors_cur.second;

	if(debug){std::cout << "\n\n - Proposed Factors " << std::endl;}
	std::pair<std::vector<double>,std::vector<int>> factors_new = q_path_factors(p_new, rev_weights, p_stop);
	std::vector<double> q_new = factors_new.first;
	P_new_revtypes = factors_new.second;

	if(debug){std::cout << "\n\n";}
	
	std::sort(q_cur.begin(), q_cur.end(), std::greater<double>());
	std::sort(q_new.begin(), q_new.end(), std::greater<double>());

	// debug=true;
	// Compute the proposal ratio (Hastings ratio).
	// Proposal ratio: q(X|Y) / q(Y|X)
	// q(X|Y) = q_L'(l',j) * q_old
	// q(Y|X) = q_L(l,j)   * q_new
	proposalRatio = q_lj_new/q_lj_cur; // q_L'(l',j) / q_L(l,j)
	if(debug){std::cout << " - q_L'(l',j)=" << q_lj_new << "; q_L(l,j)=" << q_lj_cur << "; ratio=" << proposalRatio  << " " << q_L(L_cur, l_cur)  << " " << q_L(L_new, l_new) << " " << std::endl;}
	int factors_max = std::max(q_cur.size(), q_new.size());
	int factors_idx = 0;
	while(factors_idx < factors_max){
		// Instead of doing (1/q_old) / (1/q_new) = q_new / q_old
		const double num = (factors_idx < q_new.size()) ? q_new[factors_idx] : 1.0;
		const double den = (factors_idx < q_cur.size()) ? q_cur[factors_idx] : 1.0;
		proposalRatio *= (num/den);
		if(debug){std::cout << " - Factors : cur=" << ((factors_idx < q_cur.size()) ? q_cur[factors_idx] : 0.0) << "; prop=" << ((factors_idx < q_new.size()) ? q_new[factors_idx] : 0.0) << std::endl;}
		++factors_idx;
	}
	if(debug){std::cout << " - Proposal ratio : " << proposalRatio << std::endl;}
	// debug=false;
	return proposalRatio;
}

// The posterior ratio is the ratio of the posterior probability
// of the proposed state over that of the current state. 
double ProposalReversalScenario::getPosteriorRatio(const double rev_mean){
	return getPosteriorRatio(rev_mean, L_cur, L_new, N);
}

// The posterior ratio is the ratio of the posterior probability
// of the proposed state over that of the current state. 
double ProposalReversalScenario::getPosteriorRatio(const double rev_mean, const int L_cur_, const int L_new_, const int N_){
	
	// Simple case: both paths have the same size.
	// For now there is no bias towards small reversals or specific regions for example.
	posteriorRatio = 1.0;
	if (L_cur_ == L_new_) {return 1.0;}
	
	// Number of possible reversals that can be done in a certain step (including all types: good, bad, neutral).
	const int nb_rev_step = N_*(N_+1)/2;
	
	// Compute the term [lambda^(-L_cur + L_new) * L_cur! / L_new!]
	// avoiding the calculation of factorials which can become large.
	const int L_ub = std::max(L_cur_, L_new_);
	const int L_lb = std::min(L_cur_, L_new_);
	double factorial = 1.0;
	for(int factor = L_lb; factor < L_ub; ++factor){
		if(L_cur_ > L_new_) {
			factorial *= ((factor+1.0)/rev_mean);
		} else {
			factorial *= (rev_mean/(factor+1.0));
		}
	}

	// P(X_new|lambda) / P(X_cur|lambda)
	// P(X|lambda) is defined in York, Durrett, and Nielsen (2002)
	posteriorRatio = std::pow(nb_rev_step, -L_new_+L_cur_) * factorial;
	if(debug){std::cout << " - Posterior ratio : " << posteriorRatio << "; Comb=" << std::pow(nb_rev_step, -L_new_+L_cur_) << "; Fact=" << factorial << std::endl;}
	return posteriorRatio;
}

// The acceptance probability is the minimum of one and
// the product of the posterior ratio and the proposal ratio.
double ProposalReversalScenario::getAcceptanceProb(const double rev_mean, const std::vector<double>& rev_weights, const double p_stop, const double chain_temp){
	double acceptanceProb = std::min(1.0, getProposalRatio(rev_weights, p_stop)*std::pow(getPosteriorRatio(rev_mean),chain_temp));
	// std::cout << " > Acceptance : " << acceptanceProb << "; Proposal ratio : " << proposalRatio << "; Posterior ratio : " << posteriorRatio << std::endl;
	return acceptanceProb;
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
		// debug=true;
		ReversalSampler sampler(genperm,debug);
		// debug=false;
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
	while((genperm.getBreakpoints() > 0) || ((p_stop < prob_rdm) && (rev_weights[ReversalType::BAD] > 0.0))) { // not identity and p_stop.

		if(debug){std::cout << "\n\n\n------------------------------------------------\n\n - Current nb. of breakpoints = " << genperm.getBreakpoints() << std::endl;}

		// Sample a random reversal.
		ReversalSampler sampler(genperm,rev_weights,debug);
		// sampler.debug      = true;
		ReversalRandom rev = sampler.sampleReversal(rng,true);
		//if(debug){printGenome(genperm.getExtendedPerm());}

		// std::cout << "- NEW.rev.= (" << rev.g_beg << ", " << rev.g_end << "), type=" << rev.type << ", cyctype=" << rev.cycle << ";" << std::endl;
		// for (int revtype_idx = 0; revtype_idx < ReversalType_COUNT; ++revtype_idx) {
		// 	std::cout << "\t type=" << revtype_idx <<  "; count=" << rev.rev_totals[revtype_idx] << std::endl;
		// }

		// Apply reversal.
		//genperm.debug = false;
		// const int nb_breakpoints_before = genperm.getBreakpoints();
		applyReversal(genperm, rev.g_beg, rev.g_end);
		genperm.clearBlockStatus();
		// const int nb_breakpoints_after  = genperm.getBreakpoints();
		
		// Add reversal to list of reversals.
		reversals.emplace_back(rev);
		// Sample probability of stopping if genome is sorted.
		if(genperm.getBreakpoints() == 0){prob_rdm = distr(rng);}
		
	}
	if(debug){printGenome(genperm.getExtendedPerm());}
	return reversals;
}

int RandomReversalScenario::samplePathLength(std::mt19937& rng, const int N, const double alpha, const double epsilon){
	// q(l) ~ 1 - tanh(epsilon*(l/(alpha*N)-1))
	// - alpha: lengths small than N*alpha are roughly equally likely represented.
	//          For example, alpha=0.65, means that there is more or less the same 
	//          chance to sample a size between 0 and 65% of the total size of the 
	//          path (N). From this point on, the probability drops to almost 0.
	// Compute CDF of q(l).
	std::vector<double> lprobs_cum(N,0.0);
	double l_total = 0;
	int l = 1;
	double prob_l = 1.0-std::tanh(epsilon*(static_cast<double>(l)/(alpha*static_cast<double>(N))-1.0));
	lprobs_cum[0] = prob_l;	
	for(l=2; l<=N; ++l){ // l: 1..N; l_idx: 0..N-1
		// By commenting this line, all path lenghts have the same probability of being sampled.
		prob_l = 1.0-std::tanh(epsilon*(static_cast<double>(l)/(alpha*static_cast<double>(N))-1.0));
		lprobs_cum[l-1] = lprobs_cum[l-2] + prob_l;
		l_total += prob_l;
	}
	// Sample l.
	std::uniform_real_distribution distr(0.0, l_total); // [0, rev_weights_total)
	double rdmval = distr(rng);
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
	
	if(debug) {std::cout << " Total path size (nb. of reversals) = " << N << "; length modified path = " << l << "; Start modified path: " << j << std::endl;}

	// Get genomes at these positions in the path.
	GenomeMultichrom<int> genome_B_new = getGenomes(genome_B, reversals, j, j+l);

	// Get current subpath.
	std::vector<ReversalRandom> reversals_sub_cur = getSubpath(genome_B_new, reversals, j, j+l);

	// Integrate probabilities of the subpath to the current path.
	std::vector<ReversalRandom> reversals_cur = updateReversalScenario(genome_B_new, reversals, reversals_sub_cur, j, j+l);

	// Generate a new subpath.
	std::vector<ReversalRandom> reversals_sub_new = sampleScenario(genome_B_new, rng);
	const int l_new = reversals_sub_new.size();

	// Integrate new reversals to the current path.
	std::vector<ReversalRandom> reversals_new = updateReversalScenario(genome_B_new, reversals, reversals_sub_new, j, j+l);
	const int N_new = reversals_new.size();
	if(debug) {
		std::cout << "\nOld scenario: " << std::endl;
		printSortingScenario(genome_B, reversals, j, j+l);

		std::cout << "\nNew scenario " << std::endl;
		printSortingScenario(genome_B, reversals_new, j, j+l_new);
	}
	// std::vector<ReversalRandom> reversals_new;
	return ProposalReversalScenario(reversals_cur, reversals_new, j, N, l, N_new, l_new, genome_B.n);
}

void RandomReversalScenario::printSortingScenario(GenomeMultichrom<int>& genome_B, std::vector<ReversalRandom>& reversals, int pos_beg, int pos_end) {
	// Unsigned extended permutation (it assumes that one of the permutations is the identity).
	GenomePermutation<BlockSimple> genperm(genome_B.getExtendedGenome());

	if(debug){
		std::vector<int> genome =  genperm.getExtendedPerm();
		std::cout << " Genome [start]" << std::endl;
		for (int gene : genome) {std::cout << gene << " ";}
		// Check if the current genome is one of the path ends.
		if(0 == pos_beg){std::cout << " [A; pos_beg=" << pos_beg << "] ";}
		if(0 == pos_end){std::cout << " [B] ";}
		std::cout << std::endl;
	}

	for(int rev_idx=0; rev_idx<reversals.size(); ++rev_idx){		
		
		// Apply reversal.
		ReversalRandom rev = reversals[rev_idx];
		applyReversal(genperm, rev.g_beg, rev.g_end);

		// Print genome after reversal.
		if(debug){
			std::cout << " Genome [" << (rev_idx+1) << "] after reversal (" << rev.g_beg << ", " << rev.g_end << "]: " << std::endl;
			std::vector<int> genome =  genperm.getExtendedPerm();
			for (int gene : genome) {std::cout << gene << " ";}
			// Check if the current genome is one of the path ends.
			if((rev_idx+1) == pos_beg){std::cout << " [A; pos_beg=" << pos_beg << "] ";}
			if((rev_idx+1) == pos_end){std::cout << " [B] ";}
			std::cout << std::endl;
		}
	}
}

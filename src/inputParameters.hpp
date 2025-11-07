/*******************************************************
 * Original author: Priscila Biller
 * Created: October/2025
 * License: GPL v3
 * 
 * This class provides some basic data structures to set
 * up the parameters for different methods to solve the
 * the problem of finding a reversal scenario given
 * two genomes.
 * 
 *******************************************************/

#pragma once // It avoids class redefinition.

#include <iostream>
#include <fstream>
#include <cstdlib>   // exit

#include "utils.hpp"

//////////////////////////////
// Methods to sort genomes
//////////////////////////////

enum RevMethodType : int {
	OPT = 0, // OPT is assigned the value 0
	         // Scenario with minimum number of inversions.
	MCMC     // MCMC is assigned the value 1
			 // Sampled scenario with minimum number of inversions or more.
	// Count
};
const int ReversalMethods_COUNT = 2;

class RevMethodOptions {
public:
	RevMethodType method;
	
	RevMethodOptions(){}

	virtual std::string print() {return "Undefined";}
};

class OptimalOptions : public RevMethodOptions {
public:
	OptimalOptions(){method = RevMethodType::OPT;}
	std::string print() override {return "Optimal";}
};

class McmcOptions : public RevMethodOptions {
public:

	std::string id_run{"revMCMC"};
	int nb_chains{1};

	// The MCMCMC method, also known as Parallel Tempering or MC^3, allows a 
	// faster convergence of a chain by "flattening" the landscape.
	// The "temperature" t of a chain is defined as a value between 0 and 1 
	// such that (PostRate(X))^t is the posterior rate of chain heated at t.
	// The cold chain is when t=1; this corresponds to the "normal" case. 
	// The extreme heated chain is when t=0, where the whole landscape becomes flat.
	int nb_temperatures{1}; 
	double delta_temp{0.1}; // Parameter used to compute the difference in temperature between each chain.

	// Convergence evaluation: Same method as in York, Durrett, and Nielsen (2002).
	// We use the method of Gelman and Rubin (1992) to decide when the Markov chain has converged.
	// If this flag is true,  sampling starts after convergence criteria is satisfied.
	// If this flag is false, sampling starts after the end of the pre-burn-in phase (mode used for tests).
	bool check_convergence{true}; // Criteria: sqrt(R) <= 1.1

	// Test parameters.
	// In case the flag 'check_convergence' is false,
	// then these parameters are used instead. 
	// Sampling starts directly after pre_burning_steps.
	int max_steps{100000};
	int pre_burnin_steps{1000};

	// Sampling parameters.
	int sample_interval{1000}; // Based on Larget et al. (2004)
	int sample_amount{10000};  // 10 million steps. Based on Miklos and Darling (2009).

	// Steps in which the whole state of the MCMC is saved.
	// This allows the program to recover in case of crash, time out, 
	// or any other problem.
	int backup_interval{500};

	// Steps in which stats about the chains are printed.
	int print_interval{50};

	// If the probability of bad reversals is too low, then the proposal ratio of better paths decreases:
	// - Proposal ratio increases as P(old|new) / P(new|old).
	// If the old path is composed of one or more bad reversals, and the probability of bad reversals is very low,
	// then the probability of coming back from a better new path to the old path reduces.

	// Probability values.
	double p_good{1.0};            // Based on Larget et al. (2004) value.
	double p_neutral{0.0276};      // Based on Larget et al. (2004) value.
	double p_bad{0.001};           // Based on Miklos (2003) value.
	double p_stop{0.999};

	std::vector<double> probs;
	
	McmcOptions(){method = RevMethodType::MCMC;}

	McmcOptions(const int nb_chains_, const int nb_temperatures_, const int delta_temp_, 
		const int max_steps_, const int pre_burnin_steps_, 
		const double p_good_, const double p_neutral_, const double p_bad_, 
		const double p_stop_):nb_chains(nb_chains_),nb_temperatures(nb_temperatures_),delta_temp(delta_temp_),
		max_steps(max_steps_),pre_burnin_steps(pre_burnin_steps_),
		p_good(p_good_),p_neutral(p_neutral_),p_bad(p_bad_),p_stop(p_stop_){
		method = RevMethodType::MCMC;
		probs  = {p_good, p_neutral, p_bad};
	}

	McmcOptions(std::unordered_map<std::string,std::string>& parvalues_map){
		method = RevMethodType::MCMC;

		if (parvalues_map.find("id_run") != parvalues_map.end()) {
			id_run = parvalues_map["id_run"];
		}
		if (parvalues_map.find("nb_chains") != parvalues_map.end()) {
			nb_chains = std::stoi(parvalues_map["nb_chains"]);
		}
		if (parvalues_map.find("nb_temperatures") != parvalues_map.end()) {
			nb_temperatures = std::stoi(parvalues_map["nb_temperatures"]);
		}
		if (parvalues_map.find("delta_temp") != parvalues_map.end()) {
			delta_temp = std::stod(parvalues_map["delta_temp"]);
		}
		if (parvalues_map.find("check_convergence") != parvalues_map.end()) {
			check_convergence = (parvalues_map["check_convergence"] == "true");
		}
		if (parvalues_map.find("max_steps") != parvalues_map.end()) {
			max_steps = std::stoi(parvalues_map["max_steps"]);
		}
		if (parvalues_map.find("pre_burnin_steps") != parvalues_map.end()) {
			pre_burnin_steps = std::stoi(parvalues_map["pre_burnin_steps"]);
		}
		if (parvalues_map.find("sample_interval") != parvalues_map.end()) {
			sample_interval = std::stoi(parvalues_map["sample_interval"]);
		}
		if (parvalues_map.find("sample_amount") != parvalues_map.end()) {
			sample_amount = std::stoi(parvalues_map["sample_amount"]);
		}
		if (parvalues_map.find("backup_interval") != parvalues_map.end()) {
			backup_interval = std::stoi(parvalues_map["backup_interval"]);
		}
		if (parvalues_map.find("print_interval") != parvalues_map.end()) {
			print_interval = std::stoi(parvalues_map["print_interval"]);
		}
		if (parvalues_map.find("p_good") != parvalues_map.end()) {
			p_good = std::stod(parvalues_map["p_good"]);
		}
		if (parvalues_map.find("p_neutral") != parvalues_map.end()) {
			p_neutral = std::stod(parvalues_map["p_neutral"]);
		}
		if (parvalues_map.find("p_bad") != parvalues_map.end()) {
			p_bad = std::stod(parvalues_map["p_bad"]);
		}
		if (parvalues_map.find("p_stop") != parvalues_map.end()) {
			p_stop = std::stod(parvalues_map["p_stop"]);
		}
		probs = {p_good, p_neutral, p_bad};
	}

	std::string print() override {
		return "MCMC; id_run=" + id_run
			+ ", nb_chains=" + std::to_string(nb_chains) 
			+ ", nb_temperatures=" + std::to_string(nb_temperatures) 
			+ ", delta_temp=" + std::to_string(delta_temp) 
			+ ", check_convergence=" + std::to_string(check_convergence)
			+ ", max_steps=" + std::to_string(max_steps) 
			+ ", pre_burnin_steps=" + std::to_string(pre_burnin_steps)
			+ ", sample_interval="  + std::to_string(sample_interval)
			+ ", sample_amount="    + std::to_string(sample_amount)
			+ ", backup_interval="  + std::to_string(backup_interval)
			+ ", print_interval="  + std::to_string(print_interval)
			+ ", p_good=" + std::to_string(p_good)
			+ ", p_neutral=" + std::to_string(p_neutral)
			+ ", p_bad="  + std::to_string(p_bad)
			+ ", p_stop=" + std::to_string(p_stop);
	}
};

//////////////////////////////
// Parameter parser.
//////////////////////////////

class InputPars {
public:

	InputPars(){}

	bool isValidMethod(std::string& method){
		bool validMethod = isValidPath(method);
		if (!validMethod) {
			// Convert method's name to lowercase.
			method      = lowercase(method);
			validMethod = ((method == "opt") || (method == "mcmc"));
		}
		return validMethod;
	}

	RevMethodOptions* initializeOptions(std::string& method){
		if(method == "opt"){
			return new OptimalOptions();
		} else if(method == "mcmc") {
			return new McmcOptions();
		} else {
			return new RevMethodOptions();
		}
	}

	RevMethodOptions* initializeOptions(std::string& method, std::unordered_map<std::string,std::string> parvalues_map){
		if(method == "opt"){
			return new OptimalOptions();
		} else if(method == "mcmc") {
			return new McmcOptions(parvalues_map);
		} else {
			return new RevMethodOptions();
		}
	}

	RevMethodOptions* getMethodPars(std::string& method){
		RevMethodOptions* parValues;
		if (isValidPath(method)){
			parValues = parseParameterFile(method);
		} else {
			method    = lowercase(method);
			parValues = initializeOptions(method);
		}
		return parValues;
	}

	RevMethodOptions* parseParameterFile(std::string& filename){
		std::unordered_map<std::string,std::string> parvalues_map = parseFilenameDict(filename, true);
		if (parvalues_map.find("method") == parvalues_map.end()) {
			std::cout << "ERROR! File with parameter values ('" << filename << "') is lacking the specification of a mandatory parameter: 'method=[opt|mcmc]'. Program is aborting." << std::endl;
			exit(1);
		}
		return initializeOptions(parvalues_map["method"], parvalues_map);
	}
};

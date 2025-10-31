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
	MCMC     // NEUTRAL_GOOD is assigned the value 1
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

	int nb_chains{1};
	int max_steps{100000};
	int pre_burnin_steps{1000};

	double rev_mean_range{10}; // TODO: currently not used.

	// If the probability of bad reversals is too low, then the proposal ratio of better paths decreases:
	// - Proposal ratio increases as P(old|new) / P(new|old).
	// If the old path is composed of one or more bad reversals, and the probability of bad reversals is very low,
	// then the probability of coming back from a better new path to the old path reduces.

	// Probability values.
	double p_good{1.0};            // Based on Larget et al. (2004) value.
	double p_neutralgood{0.0276};  // Based on Larget et al. (2004) value. [new category]
	double p_neutral{0.0276};      // Based on Larget et al. (2004) value.
	double p_bad{0.001};           // Based on Miklos (2003) value.
	double p_stop{0.999};

	std::vector<double> probs;
	
	McmcOptions(){method = RevMethodType::MCMC;}

	McmcOptions(const int nb_chains_, const double rev_mean_range_, const int max_steps_, const int pre_burnin_steps_, const double p_good_, const double p_neutralgood_, const double p_neutral_, const double p_bad_, const double p_stop_):nb_chains(nb_chains_),max_steps(max_steps_),pre_burnin_steps(pre_burnin_steps_),rev_mean_range(rev_mean_range_),p_good(p_good_),p_neutralgood(p_neutralgood_),p_neutral(p_neutral_),p_bad(p_bad_),p_stop(p_stop_){
		method = RevMethodType::MCMC;
		probs  = {p_good, p_neutralgood, p_neutral, p_bad};
	}

	McmcOptions(std::unordered_map<std::string,std::string>& parvalues_map){
		method = RevMethodType::MCMC;

		if (parvalues_map.find("nb_chains") != parvalues_map.end()) {
			nb_chains = std::stoi(parvalues_map["nb_chains"]);
		}
		if (parvalues_map.find("max_steps") != parvalues_map.end()) {
			max_steps = std::stoi(parvalues_map["max_steps"]);
		}
		if (parvalues_map.find("pre_burnin_steps") != parvalues_map.end()) {
			pre_burnin_steps = std::stoi(parvalues_map["pre_burnin_steps"]);
		}
		if (parvalues_map.find("rev_mean_range") != parvalues_map.end()) {
			rev_mean_range = std::stod(parvalues_map["rev_mean_range"]);
		}
		if (parvalues_map.find("p_good") != parvalues_map.end()) {
			p_good = std::stod(parvalues_map["p_good"]);
		}
		if (parvalues_map.find("p_neutralgood") != parvalues_map.end()) {
			p_neutralgood = std::stod(parvalues_map["p_neutralgood"]);
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
		probs = {p_good, p_neutralgood, p_neutral, p_bad};
	}

	std::string print() override {
		return "MCMC; nb_chains=" + std::to_string(nb_chains) 
			+ ", max_steps=" + std::to_string(max_steps) 
			+ ", pre_burnin_steps=" + std::to_string(pre_burnin_steps)
			+ ", rev_mean_range=" + std::to_string(rev_mean_range)
			+ ", p_good=" + std::to_string(p_good)
			+ ", p_neutralgood=" + std::to_string(p_neutralgood)
			+ ", p_neutral=" + std::to_string(p_neutral)
			+ ", p_bad=" + std::to_string(p_bad)
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

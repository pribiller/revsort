/*******************************************************
 * Original author: Priscila Biller
 * Created: September/2025
 * License: GPL v3
 * 
 * This file is intended to test the implementation
 * of the Sorting by Reversals method using one or more
 * genome pairs defined in a file.
 * 
 * Specifically, it checks if a sequence of reversals
 * provided by the method transforms the input permutation
 * into the identity permutation and if the number of 
 * reversals corresponds to the minimum number expected.
 * 
 * In addition to check the correctness of the method, it also 
 * saves many stats about the permutations (number of breakpoints,
 * number of cycles, etc.) and the time needed to run the method.
 * 
 * Input File - Format:
 * 
 * The input file defines the gene order of one or more pair of genomes.
 * 
 * Each pair of genomes requires 3 lines:
 *
 * - First line: identifier of the test. It can be, for example, a combination 
 * of the names of the genomes, or any other useful information that 
 * identifies the test. This information is used later in the output file.
 * 
 * - Second line: gene order of the first genome. 
 * 
 * Each gene is associated to an integer value, where the sign (positive or negative) 
 * designates the orientation of the gene in the DNA strand.
 * Although the first genome is considered as the reference genome,
 * its gene order does not need to be an identity permutation. 
 * The conversion to an identity permutation is handled internally.
 * 
 * - Third line: gene order of the second genome. 
 * 
 * - Example:
 * 
 * >TestcaseNamePair1
 * 1 -2 3 -4 5 -6
 * 1 -2 -3 4 5 6
 * >TestcaseNamePair2
 * 1 2 3 4 5 6
 * -2 -1 3 -5 4 6
 * 
 * Compile:
 * g++ run_InputFile.cpp sortByReversals.cpp findComponents_Bader2001.cpp sortOrientedByReversals_Tannier2007.cpp solveUnoriented_HannenhalliPevzner1999.cpp genome.cpp -o revsort_input
 * 
 * Run:
 * ./revsort_input 42 1 inputFile.txt
 * 
 *******************************************************/

#include <iostream>
#include <fstream>
#include <regex>
#include <unordered_map>
#include <vector>
#include <stack>
#include <list>
#include <cstdlib>   // exit
#include <cmath>	 // abs
#include <memory>	 // shared_ptr

#include "sortByReversals.hpp"

/*******************************************************
 * Some basic tests.
*******************************************************/

void testCase_generalSort(GenomeMultichrom<int>& genome_A, GenomeMultichrom<int>& genome_B, std::mt19937& rng, const bool debug){
	SortByReversals sortGenome(genome_A,genome_B,debug);
	sortGenome.sort(rng);
	// [Check 1] Check if the sorting scenario ends with the identity permutation (i.e. no breakpoints).
	bool correctSolution = sortGenome.printSolution();
	if(!correctSolution){
		std::cout << "ERROR! Problem during the sorting. Program is aborting." << std::endl;
		exit(1);
	}
	// [Check 2] Check if the number of reversals is minimal.
	correctSolution = sortGenome.printStats();
	if(!correctSolution){
		std::cout << "ERROR! Problem during the sorting. Program is aborting." << std::endl;
		exit(1);
	}
}

/*******************************************************
* Tests for the overall algorithm (find components, 
* clear hurdles, sort connected components, etc.).
*******************************************************/

// Parses pairs of gene orders in a file.
// It returns a vector where each item has the identifier of the pair, followed by the two genome objects.
std::vector<std::pair<std::string, std::pair<GenomeMultichrom<int>, GenomeMultichrom<int>>>> parseInputFile(std::string filename) {

	std::regex header_regex(">(\\S+)");
	std::regex gene_regex("([+-]?\\d+)\\S*");
	std::smatch matches;

	std::string line;
	std::ifstream inputFile(filename);

	std::vector<std::pair<std::string, std::pair<GenomeMultichrom<int>, GenomeMultichrom<int>>>> genomePairs;

	if (inputFile.is_open()) {

		std::string testId = "";
		bool readGenomeA   = false;

		std::vector<int>  genome_multichrom_A  = {};
		std::vector<bool> genome_orientation_A = {};
		std::vector<int>  genome_multichrom_B  = {};
		std::vector<bool> genome_orientation_B = {};

		while(getline(inputFile,line)) {

			// Header.
			if (std::regex_search(line, matches, header_regex)) {
				// Save the previous info.
				if(not testId.empty()){
					// Check if input genomes are properly defined.
					if (genome_multichrom_A.size() != genome_orientation_A.size()) {
						std::cout << "[" << testId << "] ERROR! Reference genome: The number of genes and signs must be the same (Found: genes= " << genome_multichrom_A.size() << "; signs= " << genome_orientation_A.size() << ")." << std::endl;
						exit(1);
					}
					if (genome_multichrom_B.size() != genome_orientation_B.size()) {
						std::cout << "[" << testId << "] ERROR! Query genome: The number of genes and signs must be the same (Found: genes= " << genome_multichrom_B.size() << "; signs= " << genome_orientation_B.size() << ")." << std::endl;
						exit(1);
					}
					if (genome_multichrom_A.size() != genome_multichrom_B.size()) {
						std::cout << "[" << testId << "] ERROR! Reference and query genomes must have the same number of genes (Found: genes ref.= " << genome_multichrom_A.size() << "; genes qry.= " << genome_multichrom_B.size() << ")." << std::endl;
						exit(1);
					}
					// Create genomes.
					GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
					GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);
					genomePairs.emplace_back(testId, std::make_pair(genome_A, genome_B));
				}

				// Initialize a new genome.		
				testId = matches[1];
				readGenomeA = true;
				genome_multichrom_A  = {};
				genome_orientation_A = {};
				genome_multichrom_B  = {};
				genome_orientation_B = {};
			// Content.
			} else {
				auto genes_begin = std::sregex_iterator(line.begin(), line.end(), gene_regex);
				auto genes_end   = std::sregex_iterator();
				for (std::sregex_iterator match = genes_begin; match != genes_end; ++match) {
					int gene      = std::stoi((*match).str());
					bool reversed = (gene < 0);
					gene = std::abs(gene);
					if(readGenomeA){
						genome_multichrom_A.push_back(gene);
						genome_orientation_A.push_back(reversed);
					} else {
						genome_multichrom_B.push_back(gene);
						genome_orientation_B.push_back(reversed);
					}
				}
				readGenomeA = false;
			}
		}

		// Save the last info.
		if(not testId.empty()){
			// Check if input genomes are properly defined.
			if (genome_multichrom_A.size() != genome_orientation_A.size()) {
				std::cout << "[" << testId << "] ERROR! Reference genome: The number of genes and signs must be the same (Found: genes= " << genome_multichrom_A.size() << "; signs= " << genome_orientation_A.size() << ")." << std::endl;
				exit(1);						
			}
			if (genome_multichrom_B.size() != genome_orientation_B.size()) {
				std::cout << "[" << testId << "] ERROR! Query genome: The number of genes and signs must be the same (Found: genes= " << genome_multichrom_B.size() << "; signs= " << genome_orientation_B.size() << ")." << std::endl;
				exit(1);						
			}
			if (genome_multichrom_A.size() != genome_multichrom_B.size()) {
				std::cout << "[" << testId << "] ERROR! Reference and query genomes must have the same number of genes (Found: genes ref.= " << genome_multichrom_A.size() << "; genes qry.= " << genome_multichrom_B.size() << ")." << std::endl;
				exit(1);						
			}
			GenomeMultichrom<int> genome_A(genome_multichrom_A, genome_orientation_A);
			GenomeMultichrom<int> genome_B(genome_multichrom_B, genome_orientation_B, genome_A.gene_labels_map);
			genomePairs.emplace_back(testId, std::make_pair(genome_A, genome_B));
		}
		inputFile.close();
	}
	return genomePairs;
}

int main(int argc, char* argv[]) {
	
	if (argc < 2) {
		std::cout << "ERROR! No extra command line argument passed other than program name. Program is aborting." << std::endl;
		exit(1);
	} else if (argc != 4) {
		std::cout << "ERROR! Wrong number of arguments. Expected: ./program [seed] [debug] [input file]. Program is aborting." << std::endl;
		exit(1);
	}

	// Parameters specified by the user.
	int const seed   = std::stoi(argv[1]); // input (command line argument): seed random number generator.
	bool const debug = (std::stoi(argv[2]) > 0);
	std::string const filename = argv[3];

	// Create a random number generator
	std::cout << "- Seed to reproduce tests: " << seed << std::endl;
	std::mt19937 rng(seed); // Seed the generator

	// Gene orders from input file.
	std::vector<std::pair<std::string, std::pair<GenomeMultichrom<int>, GenomeMultichrom<int>>>> genomePairs = parseInputFile(filename);
	for (const auto& info : genomePairs) {

		std::string testCase = info.first;
		GenomeMultichrom<int> genome_A = info.second.first;
		GenomeMultichrom<int> genome_B = info.second.second;

		std::cout << "\n\nTest : " << testCase << "\n";
		testCase_generalSort(genome_A, genome_B, rng, debug);
	}

	std::cout << std::endl << "Bye bye" << std::endl;
	return 0;
}

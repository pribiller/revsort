/*******************************************************
 * Original author: Priscila Biller
 * Created: October/2025
 * License: GPL v3
 * 
 * This class provides some functions to handle strings, 
 * files, etc.
 * 
 *******************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>   // exit
#include <vector>
#include <unordered_map>
#include <algorithm> // transform, find_if
#include <filesystem>

#include "utils.hpp"

//////////////////////////////
// String functions.
//////////////////////////////

std::vector<std::string> split(const std::string &str, char delim) {
	std::vector<std::string> result;
	std::stringstream ss(str);
	std::string item;
	while (getline(ss, item, delim)) {
	    result.push_back(item);
	}
	return result;
}

std::string lowercase(std::string str){
	std::transform(str.begin(), str.end(), str.begin(),[](unsigned char c){return std::tolower(c);});
	return str;
}

// Trim from the start (in place)
inline void ltrim(std::string &str) {
	str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](unsigned char ch) {
		return !std::isspace(ch);
	}));
}

// Trim from the end (in place)
inline void rtrim(std::string &str) {
	str.erase(std::find_if(str.rbegin(), str.rend(), [](unsigned char ch) {
		return !std::isspace(ch);
	}).base(), str.end());
}

// Trim from both ends (in place)
inline void trim(std::string &str) {
	rtrim(str);
	ltrim(str);
}

//////////////////////////////
// File functions.
//////////////////////////////

bool isValidPath(const std::string& path) {
	// Create a path object
	std::filesystem::path p(path);
	// Check if the path exists and is a regular file or directory
	return std::filesystem::exists(p);
}

std::unordered_map<std::string,std::string> parseFilenameDict(std::string& filename, bool makeLowercase){
	
	// Parse files whose lines are like "parname=parvalue".

	std::unordered_map<std::string,std::string> parvalues_map;
	std::string line;
	std::ifstream inputFile(filename);
	if (inputFile.is_open()) {
		while(getline(inputFile,line)) {
			// Skip blank lines.
			trim(line);
			if(line == ""){continue;}
			// Parse info from lines.
			std::vector<std::string> line_vals = split(line, '=');
			if(line_vals.size() != 2){
				std::cout << "ERROR! Problem to parse line '" << line << "' in the file '" << filename << "'.\nMake sure that lines follow the pattern: 'parname=parvalue'. Program is aborting." << std::endl;
				exit(1);
			}
			// Format and save info.
			std::string parname    = makeLowercase ? lowercase(line_vals[0]) : line_vals[0];
			std::string parvalue   = makeLowercase ? lowercase(line_vals[1]) : line_vals[1];
			trim(parname); trim(parvalue);
			parvalues_map[parname] = parvalue;
		}
	}
	return parvalues_map;
}


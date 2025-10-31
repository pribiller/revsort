/*******************************************************
 * Original author: Priscila Biller
 * Created: October/2025
 * License: GPL v3
 * 
 * This class provides some functions to handle strings, 
 * files, etc.
 * 
 *******************************************************/

#pragma once // It avoids class redefinition.

#include <unordered_map>
#include <vector>

//////////////////////////////
// String functions.
//////////////////////////////

std::vector<std::string> split(const std::string &str, char delim);

std::string lowercase(std::string str);

// Trim from the start (in place)
void ltrim(std::string &str);

// Trim from the end (in place)
void rtrim(std::string &str);

// Trim from both ends (in place)
void trim(std::string &str);

//////////////////////////////
// File functions.
//////////////////////////////

bool isValidPath(const std::string& path);

std::unordered_map<std::string,std::string> parseFilenameDict(std::string& filename, bool makeLowercase=false);

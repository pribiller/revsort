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

#include <boost/serialization/access.hpp>
#include <random>
#include <stdexcept>
#include <sstream>

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

//////////////////////////////////////////////////////
// Serialization with Boost of Rdm Number Generators.
//////////////////////////////////////////////////////

namespace boost { namespace serialization {
#define RNG_TPARAMS typename UIntType, size_t w, size_t n, size_t m, size_t r, UIntType a, size_t u, UIntType d, size_t s, UIntType b, size_t t, UIntType c, size_t l, UIntType f
#define RNG_TARGLIST UIntType, w, n, m, r, a, u, d, s, b, t, c, l, f

	template<typename Archive, RNG_TPARAMS>
	void load(Archive& ar, std::mersenne_twister_engine<RNG_TARGLIST>& rng, unsigned) {
		std::string text;
		ar & text;
		std::istringstream iss(text);
		if (!(iss >> rng)){
			throw std::invalid_argument("ERROR! Problem to load the random number generator.");
		}
	}

	template<typename Archive, RNG_TPARAMS>
	void save(Archive& ar, std::mersenne_twister_engine<RNG_TARGLIST> const& rng, unsigned) {
		std::ostringstream oss;
		if (!(oss << rng)){
			throw std::invalid_argument("ERROR! Problem to save the random number generator.");
		}
		std::string text = oss.str();
		ar & text;
	}

	template<typename Archive, RNG_TPARAMS>
	void serialize(Archive& ar, std::mersenne_twister_engine<RNG_TARGLIST>& rng, unsigned version) {
		if (typename Archive::is_saving()){
			save(ar, rng, version);
		} else {
			load(ar, rng, version);
		}
	}

#undef RNG_TPARAMS
#undef RNG_TARGLIST
} }

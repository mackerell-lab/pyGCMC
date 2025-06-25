// src/io/forcefield/PrmParserStructures.hpp

#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <tuple>

namespace pygcmc {
namespace io {

// Forward declarations - will be resolved by including ModelModule.hpp where needed

struct PrmParserStructures {
    
    // Helper utility functions for string conversion
    static double safe_stod(const std::string& str, const std::string& context);
    static int safe_stoi(const std::string& str, const std::string& context);
    
    // Line processing utilities
    static std::vector<std::string> tokenize(const std::string& line);
    static bool isCommentLine(const std::string& line);
    static std::string removeComments(const std::string& line);
    static std::string trim(const std::string& str);
    static std::string readContinuationLine(std::istream& input, std::string firstLine);
    
    // Section detection utilities
    static bool isAtomsSection(const std::string& line);
    static bool isBondsSection(const std::string& line);
    static bool isAnglesSection(const std::string& line);
    static bool isDihedralsSection(const std::string& line);
    static bool isImproperSection(const std::string& line);
    static bool isNonbondedSection(const std::string& line);
    static bool isNBFixSection(const std::string& line);
    
    // Parameter processing helpers
    static std::pair<std::string, std::string> make_type_pair(
        const std::string& type1, const std::string& type2);
    static std::tuple<std::string, std::string, std::string> make_type_triple(
        const std::string& type1, const std::string& type2, const std::string& type3);
    static std::tuple<std::string, std::string, std::string, std::string> make_type_quad(
        const std::string& type1, const std::string& type2, 
        const std::string& type3, const std::string& type4);
};

} // namespace io
} // namespace pygcmc
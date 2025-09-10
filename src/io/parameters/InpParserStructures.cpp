// src/io/parameters/InpParserStructures.cpp

#include "InpParserStructures.hpp"
#include <sstream>
#include <stdexcept>

namespace pygcmc {
namespace io {
namespace parameters {

std::vector<std::string> InpParserStructures::split(const std::string& str, char delim) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delim)) {
        if (!token.empty()) {
            tokens.push_back(token);
        }
    }
    return tokens;
}

std::string InpParserStructures::trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}

std::array<float, 3> InpParserStructures::parse_float_array(const std::string& str) {
    std::array<float, 3> result = {0.0f, 0.0f, 0.0f};
    std::istringstream iss(str);
    for (int i = 0; i < 3; ++i) {
        if (!(iss >> result[i])) {
            throw std::runtime_error("Failed to parse float array: " + str);
        }
    }
    return result;
}

std::vector<std::string> InpParserStructures::parse_string_vector(const std::string& str) {
    std::vector<std::string> result;
    std::istringstream iss(str);
    std::string item;
    while (iss >> item) {
        result.push_back(item);
    }
    return result;
}

std::vector<float> InpParserStructures::parse_float_vector(const std::string& str) {
    std::vector<float> result;
    std::istringstream iss(str);
    float value;
    while (iss >> value) {
        result.push_back(value);
    }
    return result;
}

std::vector<int> InpParserStructures::parse_int_vector(const std::string& str) {
    std::vector<int> result;
    std::istringstream iss(str);
    int value;
    while (iss >> value) {
        result.push_back(value);
    }
    return result;
}

} // namespace parameters
} // namespace io
} // namespace pygcmc
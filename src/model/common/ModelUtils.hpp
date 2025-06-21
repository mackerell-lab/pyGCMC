#pragma once

#ifndef PYGCMC_MODEL_UTILS_HPP
#define PYGCMC_MODEL_UTILS_HPP

#include <string>
#include <functional>
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <sstream>
#include <iomanip>
#include "ModelConstants.hpp"

namespace pygcmc {
namespace model {
namespace utils {

/**
 * @brief Hash utilities for model objects
 */
namespace hash {
    
    /**
     * @brief Hash function for string with case sensitivity option
     */
    inline std::size_t string_hash(const std::string& str, bool case_sensitive = true) {
        if (case_sensitive) {
            return std::hash<std::string>{}(str);
        } else {
            std::string lower_str = str;
            std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
            return std::hash<std::string>{}(lower_str);
        }
    }
    
    /**
     * @brief Hash function for coordinate array
     */
    inline std::size_t coordinate_hash(const std::array<double, 3>& coords, double precision = constants::EPSILON) {
        std::size_t seed = 0;
        for (double coord : coords) {
            // Round to specified precision for consistent hashing
            double rounded = std::round(coord / precision) * precision;
            seed ^= std::hash<double>{}(rounded) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
    
    /**
     * @brief Combine multiple hash values
     */
    template<typename... Args>
    inline std::size_t combine_hash(Args... args) {
        std::size_t seed = 0;
        ((seed ^= std::hash<Args>{}(args) + 0x9e3779b9 + (seed << 6) + (seed >> 2)), ...);
        return seed;
    }
    
} // namespace hash

/**
 * @brief Comparison utilities for model objects
 */
namespace compare {
    
    /**
     * @brief Compare floating point numbers with tolerance
     */
    inline bool almost_equal(double a, double b, double tolerance = constants::EPSILON) {
        return std::abs(a - b) < tolerance;
    }
    
    /**
     * @brief Compare coordinate arrays with tolerance
     */
    inline bool coordinates_equal(const std::array<double, 3>& a, 
                                 const std::array<double, 3>& b, 
                                 double tolerance = constants::EPSILON) {
        for (size_t i = 0; i < 3; ++i) {
            if (!almost_equal(a[i], b[i], tolerance)) {
                return false;
            }
        }
        return true;
    }
    
    /**
     * @brief String comparison with case sensitivity option
     */
    inline bool string_equal(const std::string& a, const std::string& b, bool case_sensitive = true) {
        if (case_sensitive) {
            return a == b;
        } else {
            return std::equal(a.begin(), a.end(), b.begin(), b.end(),
                            [](char c1, char c2) { return std::tolower(c1) == std::tolower(c2); });
        }
    }
    
} // namespace compare

/**
 * @brief ID generation utilities
 */
namespace id {
    
    /**
     * @brief Generate unique atom ID
     */
    inline std::string generate_atom_id(const std::string& type, int residue_num, 
                                       const std::string& segment_id) {
        std::ostringstream oss;
        oss << segment_id << ":" << residue_num << ":" << type;
        return oss.str();
    }
    
    /**
     * @brief Generate unique residue ID
     */
    inline std::string generate_residue_id(const std::string& resname, int residue_num,
                                          const std::string& segment_id, char chain = ' ') {
        std::ostringstream oss;
        oss << segment_id;
        if (chain != ' ') oss << ":" << chain;
        oss << ":" << residue_num << ":" << resname;
        return oss.str();
    }
    
    /**
     * @brief Generate random ID with prefix
     */
    inline std::string generate_random_id(const std::string& prefix = "obj") {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_int_distribution<> dis(10000, 99999);
        
        std::ostringstream oss;
        oss << prefix << "_" << dis(gen);
        return oss.str();
    }
    
} // namespace id

/**
 * @brief String formatting utilities
 */
namespace format {
    
    /**
     * @brief Format PDB atom name (4 characters, element right-aligned)
     */
    inline std::string format_pdb_atom_name(const std::string& name) {
        if (name.length() >= 4) return name.substr(0, 4);
        
        std::string result(4, ' ');
        
        // Check if first character is a digit (indicating branch)
        if (!name.empty() && std::isdigit(name[0])) {
            // Left justify if starts with digit
            for (size_t i = 0; i < name.length() && i < 4; ++i) {
                result[i] = name[i];
            }
        } else {
            // Right justify element symbol
            if (name.length() == 1) {
                result[1] = name[0];
            } else if (name.length() > 1) {
                result[0] = name[0];
                result[1] = name[1];
                // Add remaining characters
                for (size_t i = 2; i < name.length() && i < 4; ++i) {
                    result[i] = name[i];
                }
            }
        }
        
        return result;
    }
    
    /**
     * @brief Format coordinate with specified precision
     */
    inline std::string format_coordinate(double coord, int precision = 3) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(precision) << coord;
        return oss.str();
    }
    
    /**
     * @brief Trim whitespace from string
     */
    inline std::string trim(const std::string& str) {
        size_t start = str.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) return "";
        size_t end = str.find_last_not_of(" \t\r\n");
        return str.substr(start, end - start + 1);
    }
    
} // namespace format

/**
 * @brief Validation utilities
 */
namespace validate {
    
    /**
     * @brief Check if coordinate is valid
     */
    inline bool is_valid_coordinate(double coord) {
        return std::isfinite(coord) && 
               coord >= constants::MIN_VALID_COORDINATE && 
               coord <= constants::MAX_VALID_COORDINATE;
    }
    
    /**
     * @brief Check if mass is valid
     */
    inline bool is_valid_mass(double mass) {
        return std::isfinite(mass) && 
               mass >= constants::MIN_VALID_MASS && 
               mass <= constants::MAX_VALID_MASS;
    }
    
    /**
     * @brief Check if charge is valid
     */
    inline bool is_valid_charge(double charge) {
        return std::isfinite(charge) && 
               charge >= constants::MIN_VALID_CHARGE && 
               charge <= constants::MAX_VALID_CHARGE;
    }
    
    /**
     * @brief Check if name is valid (not empty, reasonable length)
     */
    inline bool is_valid_name(const std::string& name, size_t max_length = 8) {
        return !name.empty() && name.length() <= max_length && 
               name.find_first_not_of(" \t\r\n") != std::string::npos;
    }
    
} // namespace validate

/**
 * @brief Mathematical utilities
 */
namespace math {
    
    /**
     * @brief Calculate distance between two points
     */
    inline double distance(const std::array<double, 3>& a, const std::array<double, 3>& b) {
        double dx = a[0] - b[0];
        double dy = a[1] - b[1]; 
        double dz = a[2] - b[2];
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
    
    /**
     * @brief Calculate squared distance (avoids sqrt for performance)
     */
    inline double distance_squared(const std::array<double, 3>& a, const std::array<double, 3>& b) {
        double dx = a[0] - b[0];
        double dy = a[1] - b[1];
        double dz = a[2] - b[2];
        return dx*dx + dy*dy + dz*dz;
    }
    
    /**
     * @brief Normalize vector
     */
    inline std::array<double, 3> normalize(const std::array<double, 3>& vec) {
        double norm = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
        if (norm < constants::EPSILON) {
            return {0.0, 0.0, 0.0};
        }
        return {vec[0]/norm, vec[1]/norm, vec[2]/norm};
    }
    
} // namespace math

} // namespace utils
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_UTILS_HPP 
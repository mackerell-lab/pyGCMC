#pragma once

#include "MCMain.hpp"
#include <stdexcept>
#include <sstream>

namespace pygcmc {
namespace model {
namespace montecarlo {

/**
 * @brief Validates MCState consistency and fixes common issues
 * 
 * This is critical for preventing array bounds errors after Monte Carlo moves
 */
class MCStateValidator {
public:
    /**
     * @brief Validate and potentially fix MCState consistency
     * @param state The state to validate
     * @param autoFix If true, attempt to fix issues (e.g., activeCount > size)
     * @return true if state is valid (after any fixes), false otherwise
     */
    static bool validateAndFix(MCState& state, bool autoFix = true) {
        bool isValid = true;
        std::stringstream issues;
        
        // Check 1: activeAtomCount must not exceed atoms.size()
        if (state.activeAtomCount > static_cast<int>(state.atoms.size())) {
            issues << "activeAtomCount (" << state.activeAtomCount 
                   << ") > atoms.size() (" << state.atoms.size() << "); ";
            if (autoFix) {
                state.activeAtomCount = static_cast<int>(state.atoms.size());
                issues << "Fixed. ";
            } else {
                isValid = false;
            }
        }
        
        // Check 2: activeResidueCount must not exceed residues.size()
        if (state.activeResidueCount > static_cast<int>(state.residues.size())) {
            issues << "activeResidueCount (" << state.activeResidueCount 
                   << ") > residues.size() (" << state.residues.size() << "); ";
            if (autoFix) {
                state.activeResidueCount = static_cast<int>(state.residues.size());
                issues << "Fixed. ";
            } else {
                isValid = false;
            }
        }
        
        // Check 3: Each residue's atom range must be valid
        for (int r = 0; r < state.activeResidueCount; ++r) {
            const auto& res = state.residues[r];
            
            // Check atomStart is non-negative
            if (res.atomStart < 0) {
                issues << "Residue " << r << " has negative atomStart (" 
                       << res.atomStart << "); ";
                isValid = false;
            }
            
            // Check atom range doesn't exceed activeAtomCount
            if (res.atomStart + res.atomCount > state.activeAtomCount) {
                issues << "Residue " << r << " atom range [" << res.atomStart 
                       << ", " << (res.atomStart + res.atomCount) 
                       << ") exceeds activeAtomCount (" << state.activeAtomCount << "); ";
                if (autoFix && res.atomStart < state.activeAtomCount) {
                    // Truncate atomCount to fit within bounds
                    state.residues[r].atomCount = state.activeAtomCount - res.atomStart;
                    issues << "Fixed atomCount to " << state.residues[r].atomCount << ". ";
                } else {
                    isValid = false;
                }
            }
            
            // Check atom range doesn't exceed atoms.size()
            if (res.atomStart + res.atomCount > static_cast<int>(state.atoms.size())) {
                issues << "Residue " << r << " atom range exceeds atoms.size(); ";
                isValid = false;
            }
        }
        
        // Check 4: Ensure arrays are large enough for active counts
        if (state.atoms.size() < static_cast<size_t>(state.activeAtomCount)) {
            issues << "atoms.size() < activeAtomCount; ";
            if (autoFix) {
                state.atoms.resize(state.activeAtomCount);
                issues << "Resized atoms array. ";
            } else {
                isValid = false;
            }
        }
        
        if (state.residues.size() < static_cast<size_t>(state.activeResidueCount)) {
            issues << "residues.size() < activeResidueCount; ";
            if (autoFix) {
                state.residues.resize(state.activeResidueCount);
                issues << "Resized residues array. ";
            } else {
                isValid = false;
            }
        }
        
        // Log any issues found
        if (!issues.str().empty()) {
            if (isValid) {
                platform::log(platform::LogLevel::WARNING, 
                    "MCState validation found and fixed issues: ", issues.str());
            } else {
                platform::log(platform::LogLevel::ERROR, 
                    "MCState validation failed: ", issues.str());
            }
        }
        
        return isValid;
    }
    
    /**
     * @brief Throw exception if state is invalid
     */
    static void enforceValid(const MCState& state) {
        // Create a copy to test without modifying
        MCState testState = state;
        if (!validateAndFix(testState, false)) {
            throw std::runtime_error("MCState validation failed - inconsistent atom/residue indices");
        }
    }
};

} // namespace montecarlo
} // namespace model
} // namespace pygcmc
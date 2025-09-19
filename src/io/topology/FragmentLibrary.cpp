#include "FragmentLibrary.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

namespace pygcmc {
namespace io {
namespace topology {

void FragmentLibrary::addTemplate(const TemplateData& t) {
    byName_[t.name] = t;
    if (t.typeId >= 0) byType_[t.typeId] = t.name;
}

bool FragmentLibrary::loadFromITP(const std::string& path, const std::string& name, int typeId) {
    std::ifstream file(path);
    if (!file.good()) return false;
    
    TemplateData tmpl;
    tmpl.name = name;
    tmpl.typeId = typeId;
    
    std::string line, section;
    double totalMass = 0.0;
    double totalCharge = 0.0;
    std::vector<double> positions_x, positions_y, positions_z;
    
    while (std::getline(file, line)) {
        // Remove comments
        size_t comment = line.find(';');
        if (comment != std::string::npos) {
            line = line.substr(0, comment);
        }
        
        // Trim whitespace
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](char c) { return !std::isspace(c); }));
        line.erase(std::find_if(line.rbegin(), line.rend(), [](char c) { return !std::isspace(c); }).base(), line.end());
        
        // Check for section headers
        if (line.find("[ atoms ]") != std::string::npos) {
            section = "atoms";
            continue;
        } else if (line.find("[ bonds ]") != std::string::npos) {
            section = "bonds";
            continue;
        } else if (line.find("[") != std::string::npos && line.find("]") != std::string::npos) {
            // Other section, stop processing current section
            section = "";
            continue;
        }
        
        // Skip empty lines
        if (line.empty()) continue;
        
        // Parse atoms section
        if (section == "atoms") {
            std::istringstream iss(line);
            int idx;
            std::string atomType, resname, atomname;
            int resnr, cgnr;
            double charge, mass;
            
            // Format: nr type resnr residue atom cgnr charge mass
            if (iss >> idx >> atomType >> resnr >> resname >> atomname >> cgnr >> charge >> mass) {
                model::montecarlo::MCAtom atom;
                atom.charge = static_cast<float>(charge);
                atom.mass = static_cast<float>(mass);
                atom.type = 0;  // TODO: Map atomType string to type index
                atom.name = atomname;
                
                // Initialize position to zero (will be set from PDB if available)
                atom.x = 0.0f;
                atom.y = 0.0f;
                atom.z = 0.0f;
                atom.updatePosition();
                
                tmpl.atoms.push_back(atom);
                totalMass += mass;
                totalCharge += charge;
                
                // Store for radius calculation (if we had positions)
                positions_x.push_back(atom.position.x);
                positions_y.push_back(atom.position.y);
                positions_z.push_back(atom.position.z);
            }
        }
        // Note: bonds section parsing could be added here if needed for connectivity
        else if (section == "bonds") {
            // Could parse bond information if needed
            // Format: ai aj funct [parameters]
            // For now, we skip bond parsing as MCAtom doesn't store bond info
        }
    }
    
    // Calculate molecular weight and radius
    tmpl.molecularWeight = totalMass;
    
    // Calculate radius of gyration if we have atoms
    if (!tmpl.atoms.empty()) {
        // Since positions are not in ITP, use a default radius estimate
        // based on number of atoms (rough approximation)
        tmpl.radius = 0.15 * std::sqrt(static_cast<double>(tmpl.atoms.size()));
    }
    
    file.close();
    
    // Only add template if we successfully parsed atoms
    if (!tmpl.atoms.empty()) {
        addTemplate(tmpl);
        return true;
    }
    
    return false;
}

bool FragmentLibrary::loadFromDirectory(const std::string& /*dir*/) {
    // Optional: not implemented in stub
    return false;
}

const FragmentLibrary::TemplateData* FragmentLibrary::get(const std::string& name) const {
    auto it = byName_.find(name);
    return it == byName_.end() ? nullptr : &it->second;
}

const FragmentLibrary::TemplateData* FragmentLibrary::getByType(int typeId) const {
    auto it = byType_.find(typeId);
    if (it == byType_.end()) return nullptr;
    return get(it->second);
}

} // namespace topology
} // namespace io
} // namespace pygcmc
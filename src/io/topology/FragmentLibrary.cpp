#include "FragmentLibrary.hpp"
#include <array>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cmath>

namespace pygcmc {
namespace io {
namespace topology {

namespace {

constexpr double kAngstromToNm = 0.1;

struct PdbAtomCoord {
    std::string name;
    double xA{0.0};
    double yA{0.0};
    double zA{0.0};
};

std::string trimCopy(const std::string& s) {
    size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

std::vector<PdbAtomCoord> readPdbAtomCoords(const std::string& path) {
    std::ifstream in(path);
    if (!in.good()) return {};

    std::vector<PdbAtomCoord> atoms;
    std::string line;
    while (std::getline(in, line)) {
        if (!(line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0)) continue;
        if (line.size() < 54) continue;

        PdbAtomCoord a;
        a.name = trimCopy(line.substr(12, 4));
        try {
            a.xA = std::stod(line.substr(30, 8));
            a.yA = std::stod(line.substr(38, 8));
            a.zA = std::stod(line.substr(46, 8));
        } catch (const std::exception&) {
            continue;
        }
        atoms.push_back(std::move(a));
    }
    return atoms;
}

} // namespace

void FragmentLibrary::addTemplate(const TemplateData& t) {
    byName_[t.name] = t;
    if (t.typeId >= 0) byType_[t.typeId] = t.name;
}

bool FragmentLibrary::loadFromITP(const std::string& path, const std::string& name, int typeId,
                                  const std::string& coordinatePdbFile) {
    std::ifstream file(path);
    if (!file.good()) return false;
    
    TemplateData tmpl;
    tmpl.name = name;
    tmpl.typeId = typeId;
    
    std::string line, section;
    double totalMass = 0.0;
    double totalCharge = 0.0;
    
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
            double charge;
            double mass = 0.0;
            
            // Format: nr type resnr residue atom cgnr charge mass
            // Some GROMACS ITPs omit the mass column (mass can be derived from atomtypes);
            // accept both 7- and 8-column forms.
            if (iss >> idx >> atomType >> resnr >> resname >> atomname >> cgnr >> charge) {
                if (!(iss >> mass)) {
                    mass = 0.0;
                }
                model::montecarlo::MCAtom atom;
                atom.charge = static_cast<float>(charge);
                atom.mass = static_cast<float>(mass);
                // Placeholder until remapped to MCState atom type indices.
                atom.type = 0;
                atom.name = atomname;
                
                // Initialize position to zero (will be set from PDB if available)
                atom.x = 0.0f;
                atom.y = 0.0f;
                atom.z = 0.0f;
                atom.updatePosition();
                
                tmpl.atoms.push_back(atom);
                tmpl.atomTypeNames.push_back(atomType);
                totalMass += mass;
                totalCharge += charge;
            }
        }
        else if (section == "bonds") {
            std::istringstream iss(line);
            int ai = 0;
            int aj = 0;
            int funct = 0;
            if (!(iss >> ai >> aj)) {
                continue;
            }
            if (!(iss >> funct)) {
                funct = 0;
            }
            // ITP uses 1-based atom indices.
            ai -= 1;
            aj -= 1;
            if (ai < 0 || aj < 0) {
                continue;
            }
            tmpl.bonds.push_back(FragmentLibrary::TemplateData::Bond{ai, aj});
        }
    }

    // Validate any parsed bonds against the atom count (ignore out-of-range entries).
    if (!tmpl.bonds.empty() && !tmpl.atoms.empty()) {
        const int nAtoms = static_cast<int>(tmpl.atoms.size());
        std::vector<FragmentLibrary::TemplateData::Bond> filtered;
        filtered.reserve(tmpl.bonds.size());
        for (const auto& b : tmpl.bonds) {
            if (b.atom1 < 0 || b.atom2 < 0) {
                continue;
            }
            if (b.atom1 >= nAtoms || b.atom2 >= nAtoms) {
                continue;
            }
            filtered.push_back(b);
        }
        tmpl.bonds.swap(filtered);
    }
    
    // Calculate molecular weight and radius
    tmpl.molecularWeight = totalMass;
    
    bool coordsLoaded = false;
    if (!coordinatePdbFile.empty() && !tmpl.atoms.empty()) {
        const auto pdbAtoms = readPdbAtomCoords(coordinatePdbFile);
        if (!pdbAtoms.empty()) {
            // Map coordinates to ITP atoms. Prefer 1:1 ordering when sizes match.
            std::vector<std::array<double, 3>> coordsNm(tmpl.atoms.size(), {0.0, 0.0, 0.0});
            bool mapped = false;

            if (pdbAtoms.size() == tmpl.atoms.size()) {
                for (size_t i = 0; i < tmpl.atoms.size(); ++i) {
                    coordsNm[i] = {pdbAtoms[i].xA * kAngstromToNm,
                                   pdbAtoms[i].yA * kAngstromToNm,
                                   pdbAtoms[i].zA * kAngstromToNm};
                }
                mapped = true;
            } else {
                // Name-based sequential matching fallback.
                size_t p = 0;
                mapped = true;
                for (size_t i = 0; i < tmpl.atoms.size(); ++i) {
                    const std::string& want = tmpl.atoms[i].name;
                    while (p < pdbAtoms.size() && pdbAtoms[p].name != want) ++p;
                    if (p >= pdbAtoms.size()) {
                        mapped = false;
                        break;
                    }
                    coordsNm[i] = {pdbAtoms[p].xA * kAngstromToNm,
                                   pdbAtoms[p].yA * kAngstromToNm,
                                   pdbAtoms[p].zA * kAngstromToNm};
                    ++p;
                }
            }

            if (mapped) {
                // Center to COM in nm (use masses from ITP; fallback to geometric center)
                double mSum = 0.0;
                double cx = 0.0, cy = 0.0, cz = 0.0;
                for (size_t i = 0; i < tmpl.atoms.size(); ++i) {
                    const double m = static_cast<double>(tmpl.atoms[i].mass);
                    if (m > 0.0) {
                        mSum += m;
                        cx += m * coordsNm[i][0];
                        cy += m * coordsNm[i][1];
                        cz += m * coordsNm[i][2];
                    }
                }
                if (mSum > 0.0) {
                    cx /= mSum;
                    cy /= mSum;
                    cz /= mSum;
                } else {
                    for (const auto& c : coordsNm) {
                        cx += c[0];
                        cy += c[1];
                        cz += c[2];
                    }
                    cx /= static_cast<double>(coordsNm.size());
                    cy /= static_cast<double>(coordsNm.size());
                    cz /= static_cast<double>(coordsNm.size());
                }

                double maxDist = 0.0;
                for (size_t i = 0; i < tmpl.atoms.size(); ++i) {
                    const double x = coordsNm[i][0] - cx;
                    const double y = coordsNm[i][1] - cy;
                    const double z = coordsNm[i][2] - cz;
                    tmpl.atoms[i].x = static_cast<float>(x);
                    tmpl.atoms[i].y = static_cast<float>(y);
                    tmpl.atoms[i].z = static_cast<float>(z);
                    tmpl.atoms[i].updatePosition();
                    maxDist = std::max(maxDist, std::sqrt(x * x + y * y + z * z));
                }
                tmpl.radius = maxDist;
                coordsLoaded = true;
            }
        }
    }

    // Fallback radius estimate if coordinates were not loaded.
    if (!coordsLoaded && !tmpl.atoms.empty()) {
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

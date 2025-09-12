// src/io/TrajectoryWriter.cpp

#include "TrajectoryWriter.hpp"
#include "model/montecarlo/MCStructures.hpp"
#include <iomanip>
#include <ctime>
#include <stdexcept>
#include <map>
#include <cctype>
#include <set>
#include <algorithm>

namespace pygcmc {
namespace io {

namespace {
// Known two-letter element symbols (common subset sufficient for writer)
const std::set<std::string> kTwoLetterElements = {
    "He","Li","Be","Ne","Na","Mg","Al","Si","Cl","Ar","Ca","Sc","Ti","Cr","Mn","Fe",
    "Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Zr","Nb","Mo","Tc",
    "Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","Xe","Cs","Ba","La","Ce","Pr","Nd",
    "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","Re","Os","Ir",
    "Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np",
    "Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Y","V"
};

inline std::string normalizeElement(const std::string& raw) {
    if (raw.empty()) return "";
    // Extract leading letters only
    std::string letters;
    for (char c : raw) {
        if (std::isalpha(static_cast<unsigned char>(c))) letters.push_back(c);
        else break;
    }
    if (letters.empty()) return "";

    // Candidate of length >= 2
    if (letters.size() >= 2) {
        std::string cand;
        cand.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(letters[0]))));
        cand.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(letters[1]))));
        if (kTwoLetterElements.count(cand)) return cand;
    }
    // Fallback to first letter uppercased
    std::string one(1, static_cast<char>(std::toupper(static_cast<unsigned char>(letters[0]))));
    return one;
}

inline std::string inferElementSymbol(const std::string& typeName, const std::string& atomName) {
    // Prefer typeName when available, then atomName, else default to C
    std::string e = normalizeElement(typeName);
    if (e.empty()) e = normalizeElement(atomName);
    if (e.empty()) e = "C";
    return e;
}

inline std::string formatAtomNameForPDB(const std::string& element, const std::string& raw) {
    // Build a 4-char atom name following PDB alignment conventions.
    // If element is one-letter, name is right-justified in 4 columns; if two-letter, left-justified.
    std::string name = raw;
    if (name.size() > 4) name = name.substr(0, 4);
    // Strip spaces
    name.erase(std::remove_if(name.begin(), name.end(), ::isspace), name.end());
    std::ostringstream oss;
    if (element.size() == 1) {
        oss << std::right << std::setw(4) << name;
    } else {
        oss << std::left << std::setw(4) << name;
    }
    return oss.str();
}
} // anonymous namespace

TrajectoryWriter::TrajectoryWriter(const std::string& filename, Format format)
    : format_(format), filename_(filename), frameCount_(0), isOpen_(false) {
    
    file_.open(filename);
    if (!file_.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }
    isOpen_ = true;
    
    // Write header for PDB format
    if (format_ == Format::PDB) {
        file_ << "REMARK   PyGCMC Trajectory File\n";
        file_ << "REMARK   Created: " << std::time(nullptr) << "\n";
    }
}

TrajectoryWriter::~TrajectoryWriter() {
    close();
}

void TrajectoryWriter::close() {
    if (isOpen_ && file_.is_open()) {
        if (format_ == Format::PDB) {
            file_ << "END\n";
        }
        file_.close();
        isOpen_ = false;
    }
}

void TrajectoryWriter::writeFrame(const model::MCState& state, int frameNumber) {
    if (!isOpen_) return;
    
    int frame = (frameNumber >= 0) ? frameNumber : frameCount_;
    
    switch (format_) {
        case Format::PDB:
            writePDBFrame(state, frame);
            break;
        case Format::XYZ:
            writeXYZFrame(state, frame);
            break;
        case Format::DAT:
            writeDATFrame(state, frame);
            break;
        case Format::TOP:
            writeTOPData(state);
            break;
    }
    
    frameCount_++;
}

void TrajectoryWriter::writePDBFrame(const model::MCState& state, int frameNumber) {
    file_ << "MODEL " << std::setw(8) << frameNumber + 1 << "\n";
    file_ << "CRYST1" 
          << std::fixed << std::setprecision(3)
          << std::setw(9) << state.info.box[0] * 10.0  // nm to Angstrom
          << std::setw(9) << state.info.box[1] * 10.0
          << std::setw(9) << state.info.box[2] * 10.0
          << std::setw(7) << "90.00"
          << std::setw(7) << "90.00"
          << std::setw(7) << "90.00"
          << " P 1           1\n";
    
    int atomSerial = 1;
    for (size_t i = 0; i < state.residues.size(); ++i) {
        const auto& residue = state.residues[i];
        if (!residue.active) continue;
        
        std::string resName = !residue.resname.empty() ? residue.resname : std::string("MOL");
        int resSeq = (residue.resid > 0) ? residue.resid : static_cast<int>(i + 1);
        
        for (int j = 0; j < residue.atomCount; ++j) {
            int atomIdx = residue.atomStart + j;
            if (atomIdx >= static_cast<int>(state.atoms.size())) continue;
            
            const auto& atom = state.atoms[atomIdx];
            
            // Determine element and atom name
            std::string typeName = state.atomTypes.getTypeName(atom.type);
            std::string element = inferElementSymbol(typeName, atom.name);
            std::string baseName = !atom.name.empty() ? atom.name : (element + std::to_string(j + 1));
            std::string atomNamePDB = formatAtomNameForPDB(element, baseName);

            file_ << formatPDBAtom(atomSerial++, 
                                   atomNamePDB,
                                   resName,
                                   resSeq,
                                   atom.x * 10.0,  // nm to Angstrom
                                   atom.y * 10.0,
                                   atom.z * 10.0,
                                   1.0,  // occupancy
                                   0.0,  // temperature factor
                                   element);
        }
    }
    
    file_ << "ENDMDL\n";
}

void TrajectoryWriter::writeXYZFrame(const model::MCState& state, int frameNumber) {
    // Count actual active atoms that can be written
    int activeAtoms = 0;
    for (const auto& residue : state.residues) {
        if (residue.active) {
            // Only count atoms that are within valid range
            for (int j = 0; j < residue.atomCount; ++j) {
                int atomIdx = residue.atomStart + j;
                if (atomIdx < static_cast<int>(state.atoms.size())) {
                    activeAtoms++;
                }
            }
        }
    }
    
    file_ << activeAtoms << "\n";
    file_ << "Frame " << frameNumber << " Energy: ";
    
    // Calculate total energy
    double totalEnergy = 0.0;
    for (const auto& residue : state.residues) {
        if (residue.active) {
            totalEnergy += residue.energy_vdw + residue.energy_elec;
        }
    }
    file_ << std::fixed << std::setprecision(6) << totalEnergy << "\n";
    
    // Write atoms
    for (const auto& residue : state.residues) {
        if (!residue.active) continue;
        
        for (int j = 0; j < residue.atomCount; ++j) {
            int atomIdx = residue.atomStart + j;
            if (atomIdx >= static_cast<int>(state.atoms.size())) continue;
            
            const auto& atom = state.atoms[atomIdx];
            
            // Determine element using robust inference
            std::string typeName = state.atomTypes.getTypeName(atom.type);
            std::string element = inferElementSymbol(typeName, atom.name);
            
            file_ << std::setw(2) << element
                  << std::fixed << std::setprecision(8)
                  << std::setw(16) << atom.x * 10.0  // nm to Angstrom
                  << std::setw(16) << atom.y * 10.0
                  << std::setw(16) << atom.z * 10.0
                  << "\n";
        }
    }
}

void TrajectoryWriter::writeDATFrame(const model::MCState& state, int frameNumber) {
    // Write simple data format: frame, n_molecules, total_energy, box_volume
    int activeMolecules = 0;
    double totalEnergy = 0.0;
    
    for (const auto& residue : state.residues) {
        if (residue.active) {
            activeMolecules++;
            totalEnergy += residue.energy_vdw + residue.energy_elec;
        }
    }
    
    double volume = state.info.box[0] * state.info.box[1] * state.info.box[2];
    
    file_ << std::setw(10) << frameNumber
          << std::setw(10) << activeMolecules
          << std::fixed << std::setprecision(6)
          << std::setw(15) << totalEnergy
          << std::setw(15) << volume
          << "\n";
}

void TrajectoryWriter::writeTOPData(const model::MCState& state) {
    // Write topology information
    file_ << "# PyGCMC Topology File\n";
    file_ << "# Box: " << state.info.box[0] << " " << state.info.box[1] 
          << " " << state.info.box[2] << " nm\n";
    file_ << "# Temperature: " << 1.0 / (state.info.beta * model::MCInfo::BOLTZMANN) << " K\n";
    file_ << "\n";
    
    // Write atom types
    file_ << "[ atomtypes ]\n";
    file_ << "; type  mass  charge  sigma  epsilon\n";
    
    for (size_t i = 0; i < static_cast<size_t>(state.forcefield.numTotalTypes); ++i) {
        std::string typeName = state.atomTypes.getTypeName(static_cast<int>(i));
        if (typeName.empty()) {
            typeName = "T" + std::to_string(i);
        }
        
        // Get LJ parameters: prefer per-type if available, otherwise from diagonal of mixing matrix
        double sigma, eps;
        if (!state.forcefield.ljSigmaType.empty() && i < state.forcefield.ljSigmaType.size()) {
            sigma = state.forcefield.ljSigmaType[i];
            eps = state.forcefield.ljEpsType[i];
        } else if (!state.forcefield.ljSigma.empty()) {
            // Get from diagonal of mixing matrix (self-interaction)
            int idx = i * state.forcefield.numTotalTypes + i;
            sigma = state.forcefield.ljSigma[idx];
            eps = state.forcefield.ljEps[idx];
        } else {
            // Fallback defaults
            sigma = 0.3;
            eps = 0.5;
        }
        
        file_ << std::setw(6) << typeName
              << std::setw(8) << "12.01"  // Default mass
              << std::setw(8) << "0.0"    // Default charge
              << std::fixed << std::setprecision(4)
              << std::setw(8) << sigma
              << std::setw(8) << eps
              << "\n";
    }
    
    file_ << "\n[ molecules ]\n";
    file_ << "; molname  count\n";
    
    // Count molecule types
    std::map<std::string, int> molCounts;
    for (size_t i = 0; i < state.residues.size(); ++i) {
        if (state.residues[i].active) {
            std::string molName;
            if (!state.residues[i].resname.empty()) {
                molName = state.residues[i].resname;
            } else if (state.residues[i].type >= 0) {
                molName = "MOL" + std::to_string(state.residues[i].type);
            } else {
                molName = "MOL" + std::to_string(i);
            }
            molCounts[molName]++;
        }
    }
    
    for (const auto& [name, count] : molCounts) {
        file_ << std::setw(10) << name << std::setw(8) << count << "\n";
    }
}

std::string TrajectoryWriter::formatPDBAtom(int serial, const std::string& name,
                                           const std::string& resName, int resSeq,
                                           double x, double y, double z,
                                           double occupancy, double tempFactor,
                                           const std::string& element) {
    std::ostringstream oss;
    // PDB format: strict column alignment per PDB specification
    // Columns 1-6: "ATOM  " or "HETATM"
    // Columns 7-11: Atom serial number (right-aligned)
    // Column 12: blank
    // Columns 13-16: Atom name (left-aligned if 4 chars, else centered)
    // Column 17: Alternate location indicator  
    // Columns 18-20: Residue name (right-aligned)
    // Column 21: blank
    // Column 22: Chain identifier
    // Columns 23-26: Residue sequence number (right-aligned)
    // Column 27: Code for insertion of residues
    // Columns 28-30: blank
    // Columns 31-38: X coordinate (8.3 format)
    // Columns 39-46: Y coordinate (8.3 format)
    // Columns 47-54: Z coordinate (8.3 format)
    // Columns 55-60: Occupancy (6.2 format)
    // Columns 61-66: Temperature factor (6.2 format)
    // Columns 67-76: blank
    // Columns 77-78: Element symbol (right-aligned)
    
    oss << "ATOM  ";
    oss << std::right << std::setw(5) << serial;
    oss << " ";
    
    // Atom name: centered if less than 4 chars
    if (name.length() < 4) {
        oss << " " << std::left << std::setw(3) << name;
    } else {
        oss << std::left << std::setw(4) << name.substr(0, 4);
    }
    
    oss << " ";
    oss << std::right << std::setw(3) << resName.substr(0, 3);
    oss << " ";
    oss << " ";  // Chain ID (single space)
    oss << std::right << std::setw(4) << resSeq;
    oss << "    ";  // Insertion code + blanks
    
    oss << std::fixed << std::setprecision(3);
    oss << std::right << std::setw(8) << x;
    oss << std::right << std::setw(8) << y;
    oss << std::right << std::setw(8) << z;
    
    oss << std::fixed << std::setprecision(2);
    oss << std::right << std::setw(6) << occupancy;
    oss << std::right << std::setw(6) << tempFactor;
    
    oss << "          ";  // 10 blanks
    oss << std::right << std::setw(2) << element;
    oss << "\n";
    
    return oss.str();
}

void TrajectoryWriter::writeStatistics(const model::MCInfo::Statistics& stats) {
    if (!isOpen_ || format_ != Format::DAT) return;
    
    file_ << "# Statistics:\n";
    file_ << "# Total moves: " << stats.totalMoves << "\n";
    file_ << "# Accepted moves: " << stats.acceptedMoves << "\n";
    file_ << "# Acceptance rate: " 
          << (stats.totalMoves > 0 ? static_cast<double>(stats.acceptedMoves) / stats.totalMoves : 0.0)
          << "\n";
    file_ << "# Insertion attempts: " << stats.insertionAttempts << "\n";
    file_ << "# Accepted insertions: " << stats.acceptedInsertions << "\n";
    file_ << "# Deletion attempts: " << stats.deletionAttempts << "\n";
    file_ << "# Accepted deletions: " << stats.acceptedDeletions << "\n";
}

void TrajectoryWriter::writeTopology(const model::MCState& state) {
    if (format_ != Format::TOP) {
        format_ = Format::TOP;
    }
    writeTOPData(state);
}

// DataWriter implementation
DataWriter::DataWriter(const std::string& filename) : headerWritten_(false) {
    file_.open(filename);
    if (!file_.is_open()) {
        throw std::runtime_error("Failed to open data file: " + filename);
    }
}

DataWriter::~DataWriter() {
    close();
}

void DataWriter::close() {
    if (file_.is_open()) {
        file_.close();
    }
}

void DataWriter::writeHeader(const std::vector<std::string>& columns) {
    if (headerWritten_) return;
    
    file_ << "# ";
    for (size_t i = 0; i < columns.size(); ++i) {
        file_ << std::setw(14) << columns[i];
        if (i < columns.size() - 1) file_ << " ";
    }
    file_ << "\n";
    headerWritten_ = true;
}

void DataWriter::writeRow(const std::vector<double>& values) {
    for (size_t i = 0; i < values.size(); ++i) {
        file_ << std::fixed << std::setprecision(6)
              << std::setw(15) << values[i];
        if (i < values.size() - 1) file_ << " ";
    }
    file_ << "\n";
}

void DataWriter::writeComment(const std::string& comment) {
    file_ << "# " << comment << "\n";
}

} // namespace io
} // namespace pygcmc

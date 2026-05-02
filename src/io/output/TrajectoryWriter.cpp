#include "TrajectoryWriter.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <stdexcept>
#include <algorithm>
#include <set>
#include <cctype>
#include <map>

namespace pygcmc {
namespace io {
namespace output {

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
    if (e.empty()) e = "C";  // Default to carbon
    return e;
}

// Overload for single parameter (backward compatibility)
inline std::string inferElementSymbol(const std::string& atomName) {
    return inferElementSymbol("", atomName);
}

inline std::string formatAtomNameForPDB(const std::string& element, const std::string& raw) {
    // Build a 4-char atom name following PDB alignment conventions.
    // If element is one-letter, name is right-justified in 4 columns; if two-letter, left-justified.
    std::string name = raw;
    if (name.size() > 4) name = name.substr(0, 4);
    // Strip spaces
    name.erase(std::remove_if(name.begin(), name.end(), ::isspace), name.end());
    if (name.empty()) name = element;  // Use element as fallback

    std::ostringstream oss;
    if (element.size() == 1 && name.size() < 4) {
        // One-letter element: right-justify in 4 columns with leading space
        oss << " " << std::left << std::setw(3) << name;
    } else {
        // Two-letter element or 4-char name: left-justify
        oss << std::left << std::setw(4) << name;
    }
    return oss.str();
}
} // anonymous namespace

TrajectoryWriter::TrajectoryWriter(const Config& config)
    : config_(config), frameCount_(0), isOpen_(false) {
}

TrajectoryWriter::~TrajectoryWriter() {
    close();
}

bool TrajectoryWriter::open(const std::string& filename) {
    if (isOpen_) {
        close();
    }

    std::string fname = filename.empty() ? generateFilename(0) : filename;
    file_.open(fname);

    if (!file_.is_open()) {
        std::cerr << "Failed to open trajectory file: " << fname << std::endl;
        return false;
    }

    isOpen_ = true;
    currentFilename_ = fname;

    // Write header for PDB format
    if (config_.format == "pdb") {
        file_ << "REMARK   PyGCMC Trajectory File\n";
        std::time_t now = std::time(nullptr);
        file_ << "REMARK   Created: " << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S") << "\n";
    }

    return true;
}

void TrajectoryWriter::close() {
    if (isOpen_ && file_.is_open()) {
        if (config_.format == "pdb") {
            file_ << "END\n";
        }
        file_.close();
        isOpen_ = false;
    }
}

bool TrajectoryWriter::writeTrajectory(const model::montecarlo::MCState& state,
                                      int step,
                                      const std::string& filename) {
    // Handle file opening
    if (!isOpen_) {
        if (!open(filename.empty() ? generateFilename(step) : filename)) {
            return false;
        }
    } else if (!filename.empty() && filename != currentFilename_) {
        // Different filename requested, reopen
        close();
        if (!open(filename)) {
            return false;
        }
    }

    // Write based on format
    bool success = false;
    if (config_.format == "pdb") {
        success = writePDB(state, step, file_);
    } else if (config_.format == "xyz") {
        success = writeXYZ(state, step, file_);
    } else if (config_.format == "dat") {
        success = writeDAT(state, step, file_);
    } else if (config_.format == "top") {
        success = writeTOP(state, step, file_);
    }

    if (success) {
        frameCount_++;
        file_.flush();  // Ensure data is written
    }

    // If not writing continuous trajectory, close after each write
    if (!config_.continuousFile) {
        close();
    }

    return success;
}

bool TrajectoryWriter::writePDB(const model::montecarlo::MCState& state,
                               int /* step */,
                               std::ofstream& out) {
    // Write MODEL for multi-frame PDB
    if (config_.multiFrame) {
        out << "MODEL " << std::setw(8) << frameCount_ + 1 << "\n";
    }

    // Write crystal information
    std::vector<double> box;
    if (!state.periodicBox.empty() && state.periodicBox.size() >= 3) {
        box = {state.periodicBox[0], state.periodicBox[1], state.periodicBox[2]};
    } else {
        // info.box is a fixed array [3]
        box = {state.info.box[0], state.info.box[1], state.info.box[2]};
    }
    writePDBBox(out, box);

    // Write atoms
    int atomSerial = 1;

    // Check if we have residues to work with
    if (!state.residues.empty()) {
        // Write by residues
        for (size_t i = 0; i < state.residues.size(); ++i) {
            const auto& residue = state.residues[i];
            if (!residue.active) continue;

            std::string resName = !residue.resname.empty() ? residue.resname : "UNK";
            int resSeq = (residue.resid > 0) ? residue.resid : static_cast<int>(i + 1);

            for (int j = 0; j < residue.atomCount; ++j) {
                int atomIdx = residue.atomStart + j;
                if (atomIdx >= static_cast<int>(state.atoms.size())) continue;

                const auto& atom = state.atoms[atomIdx];
                // Try to get type name for better element inference
                std::string typeName = "";
                if (atom.type >= 0) {
                    typeName = state.atomTypes.getTypeName(atom.type);
                }
                std::string element = inferElementSymbol(typeName, atom.name);
                std::string atomNamePDB = formatAtomNameForPDB(element, atom.name);

                out << formatPDBAtom(atomSerial++, atomNamePDB, resName, resSeq,
                                   atom.x * 10.0,  // nm to Angstrom
                                   atom.y * 10.0,
                                   atom.z * 10.0,
                                   1.0,  // occupancy
                                   0.0,  // temperature factor
                                   element);
            }
        }
    } else {
        // Fallback: write atoms directly
        for (int i = 0; i < state.activeAtomCount; ++i) {
            const auto& atom = state.atoms[i];
            // Try to get type name for better element inference
            std::string typeName = "";
            if (atom.type >= 0) {
                typeName = state.atomTypes.getTypeName(atom.type);
            }
            std::string element = inferElementSymbol(typeName, atom.name);
            std::string atomNamePDB = formatAtomNameForPDB(element, atom.name);

            out << formatPDBAtom(atomSerial++, atomNamePDB, "UNK", i + 1,
                               atom.x * 10.0,  // nm to Angstrom
                               atom.y * 10.0,
                               atom.z * 10.0,
                               1.0,  // occupancy
                               0.0,  // temperature factor
                               element);
        }
    }

    if (config_.multiFrame) {
        out << "ENDMDL\n";
    }

    return true;
}

bool TrajectoryWriter::writeXYZ(const model::montecarlo::MCState& state,
                               int step,
                               std::ofstream& out) {
    // Count actual atoms
    int atomCount = 0;
    double totalEnergy = 0.0;

    if (!state.residues.empty()) {
        for (const auto& residue : state.residues) {
            if (residue.active) {
                for (int j = 0; j < residue.atomCount; ++j) {
                    int atomIdx = residue.atomStart + j;
                    if (atomIdx < static_cast<int>(state.atoms.size())) {
                        atomCount++;
                    }
                }
                totalEnergy += residue.energy_vdw + residue.energy_elec;
            }
        }
    } else {
        atomCount = state.activeAtomCount;
    }

    // Write XYZ header
    out << atomCount << "\n";
    out << "Frame " << step << " Energy: " << std::fixed << std::setprecision(6)
        << totalEnergy << "\n";

    // Write atoms
    if (!state.residues.empty()) {
        for (const auto& residue : state.residues) {
            if (!residue.active) continue;

            for (int j = 0; j < residue.atomCount; ++j) {
                int atomIdx = residue.atomStart + j;
                if (atomIdx >= static_cast<int>(state.atoms.size())) continue;

                const auto& atom = state.atoms[atomIdx];
                // Try to get type name for better element inference
                std::string typeName = "";
                if (atom.type >= 0) {
                    typeName = state.atomTypes.getTypeName(atom.type);
                }
                std::string element = inferElementSymbol(typeName, atom.name);

                out << std::setw(2) << element
                    << std::fixed << std::setprecision(8)
                    << std::setw(16) << atom.x * 10.0  // nm to Angstrom
                    << std::setw(16) << atom.y * 10.0
                    << std::setw(16) << atom.z * 10.0
                    << "\n";
            }
        }
    } else {
        for (int i = 0; i < state.activeAtomCount; ++i) {
            const auto& atom = state.atoms[i];
            // Try to get type name for better element inference
            std::string typeName = "";
            if (atom.type >= 0) {
                typeName = state.atomTypes.getTypeName(atom.type);
            }
            std::string element = inferElementSymbol(typeName, atom.name);

            out << std::setw(2) << element
                << std::fixed << std::setprecision(8)
                << std::setw(16) << atom.x * 10.0
                << std::setw(16) << atom.y * 10.0
                << std::setw(16) << atom.z * 10.0
                << "\n";
        }
    }

    return true;
}

bool TrajectoryWriter::writeDAT(const model::montecarlo::MCState& state,
                               int step,
                               std::ofstream& out) {
    // Write simple data format: step, n_molecules, total_energy, box_volume
    int activeMolecules = 0;
    double totalEnergy = 0.0;

    // Only access residues if they exist
    if (!state.residues.empty()) {
        for (const auto& residue : state.residues) {
            if (residue.active) {
                activeMolecules++;
                totalEnergy += residue.energy_vdw + residue.energy_elec;
            }
        }
    } else {
        // Fallback: count active atoms as molecules
        activeMolecules = state.activeAtomCount;
    }

    // Calculate volume from periodic box or info.box
    double volume = 0.0;
    if (!state.periodicBox.empty() && state.periodicBox.size() >= 3) {
        volume = state.periodicBox[0] * state.periodicBox[1] * state.periodicBox[2];
    } else {
        volume = state.info.box[0] * state.info.box[1] * state.info.box[2];
    }

    out << std::setw(10) << step
        << std::setw(10) << activeMolecules
        << std::fixed << std::setprecision(6)
        << std::setw(15) << totalEnergy
        << std::setw(15) << volume
        << "\n";

    return true;
}

bool TrajectoryWriter::writeTOP(const model::montecarlo::MCState& state,
                               int /* step */,
                               std::ofstream& out) {
    // Write topology information
    out << "# PyGCMC Topology File\n";

    // Write box information if available
    if (!state.periodicBox.empty() && state.periodicBox.size() >= 3) {
        out << "# Box: " << state.periodicBox[0] << " " << state.periodicBox[1]
            << " " << state.periodicBox[2] << " nm\n";
    } else {
        // info.box is a fixed array [3]
        out << "# Box: " << state.info.box[0] << " " << state.info.box[1]
            << " " << state.info.box[2] << " nm\n";
    }

    // Write temperature if beta is available
    if (state.info.beta > 0) {
        // Use the BOLTZMANN constant from MCInfo
        double temperature = 1.0 / (state.info.beta * model::montecarlo::MCInfo::BOLTZMANN);
        out << "# Temperature: " << std::fixed << std::setprecision(2)
            << temperature << " K\n";
    }
    out << "\n";

    // Write atom types section
    out << "[ atomtypes ]\n";
    out << "; type  mass  charge  sigma  epsilon\n";

    // Write each atom type
    for (int i = 0; i < state.forcefield.numTotalTypes; ++i) {
        std::string typeName = "T" + std::to_string(i);

        // Get LJ parameters
        double sigma = 0.3;  // Default
        double eps = 0.5;    // Default

        // Try to get actual parameters if available
        if (!state.forcefield.ljSigmaType.empty() && i < static_cast<int>(state.forcefield.ljSigmaType.size())) {
            sigma = state.forcefield.ljSigmaType[i];
            eps = state.forcefield.ljEpsType[i];
        } else if (!state.forcefield.ljSigma.empty()) {
            // Get from mixing matrix diagonal
            int idx = i * state.forcefield.numTotalTypes + i;
            if (idx < static_cast<int>(state.forcefield.ljSigma.size())) {
                sigma = state.forcefield.ljSigma[idx];
                eps = state.forcefield.ljEps[idx];
            }
        }

        out << std::setw(6) << typeName
            << std::setw(8) << "12.01"  // Default mass
            << std::setw(8) << "0.0"    // Default charge
            << std::fixed << std::setprecision(4)
            << std::setw(8) << sigma
            << std::setw(8) << eps
            << "\n";
    }

    // Write molecules section
    out << "\n[ molecules ]\n";
    out << "; molname  count\n";

    // Count molecule types
    std::map<std::string, int> molCounts;
    if (!state.residues.empty()) {
        for (size_t i = 0; i < state.residues.size(); ++i) {
            if (state.residues[i].active) {
                std::string molName = !state.residues[i].resname.empty() ?
                    state.residues[i].resname : ("MOL" + std::to_string(state.residues[i].type));
                molCounts[molName]++;
            }
        }
    } else {
        // Fallback if no residues
        molCounts["MOL"] = state.activeAtomCount;
    }

    for (const auto& [name, count] : molCounts) {
        out << std::setw(10) << name << std::setw(8) << count << "\n";
    }

    return true;
}

std::string TrajectoryWriter::generateFilename(int step) const {
    std::stringstream ss;
    ss << config_.prefix << "_";
    ss << std::setfill('0') << std::setw(8) << step;

    if (config_.format == "pdb") {
        ss << ".pdb";
    } else if (config_.format == "xyz") {
        ss << ".xyz";
    } else if (config_.format == "dat") {
        ss << ".dat";
    } else {
        ss << ".traj";
    }

    return ss.str();
}

void TrajectoryWriter::writePDBHeader(std::ofstream& out, int step) const {
    out << "REMARK GCMC Trajectory\n";
    out << "REMARK Step: " << step << "\n";

    std::time_t now = std::time(nullptr);
    out << "REMARK Generated: "
        << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S")
        << "\n";
}

void TrajectoryWriter::writePDBBox(std::ofstream& out,
                                  const std::vector<double>& box) const {
    if (box.size() >= 3) {
        out << "CRYST1"
            << std::fixed << std::setprecision(3)
            << std::setw(9) << box[0] * 10.0  // nm to Angstrom
            << std::setw(9) << box[1] * 10.0
            << std::setw(9) << box[2] * 10.0
            << std::setw(7) << "90.00"
            << std::setw(7) << "90.00"
            << std::setw(7) << "90.00"
            << " P 1           1\n";
    }
}

std::string TrajectoryWriter::formatPDBAtom(int serial,
                                           const std::string& atomName,
                                           const std::string& resName,
                                           int resSeq,
                                           double x, double y, double z,
                                           double occupancy,
                                           double tempFactor,
                                           const std::string& element) const {
    std::ostringstream oss;

    // PDB format: strict column alignment
    oss << "ATOM  ";
    oss << std::right << std::setw(5) << serial;
    oss << " ";

    // Atom name: centered if less than 4 chars
    if (atomName.length() < 4) {
        oss << " " << std::left << std::setw(3) << atomName;
    } else {
        oss << std::left << std::setw(4) << atomName.substr(0, 4);
    }

    oss << " ";
    oss << std::right << std::setw(3) << resName.substr(0, 3);
    oss << " ";
    oss << " ";  // Chain ID
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

} // namespace output
} // namespace io
} // namespace pygcmc

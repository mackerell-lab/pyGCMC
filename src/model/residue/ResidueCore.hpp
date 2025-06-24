#pragma once

#include <string>
#include <array>

namespace pygcmc {
namespace model {
namespace residue {

/**
 * @brief Secondary structure types following CHARMM conventions
 */
enum class SecondaryStructure {
    NONE,
    HELIX_ALPHA_RIGHT,    // Right-handed alpha (default)
    HELIX_OMEGA_RIGHT,    // Right-handed omega
    HELIX_PI_RIGHT,       // Right-handed pi
    HELIX_GAMMA_RIGHT,    // Right-handed gamma
    HELIX_310_RIGHT,      // Right-handed 3/10
    HELIX_ALPHA_LEFT,     // Left-handed alpha
    HELIX_OMEGA_LEFT,     // Left-handed omega
    HELIX_GAMMA_LEFT,     // Left-handed gamma
    HELIX_27_RIBBON,      // 2/7 ribbon/helix
    HELIX_POLYPROLINE,    // Polyproline
    SHEET_PARALLEL,       // Parallel beta sheet
    SHEET_ANTIPARALLEL    // Anti-parallel beta sheet
};

/**
 * @brief Sheet strand information
 */
struct SheetStrand {
    int strandNumber;     // Strand number in current sheet
    std::string sheetId;  // Sheet identifier
    int numStrands;       // Number of strands in current sheet
    int sense;            // 0: first, 1: parallel, -1: antiparallel
    // Hydrogen bond information
    std::string atom1;    // First atom name
    std::string atom2;    // Second atom name
    int resnum1;          // First residue number
    int resnum2;          // Second residue number
};

/**
 * @brief Disulfide bond information
 */
struct SSBond {
    int serialNumber;     // Bond serial number
    std::string partner;  // Partner residue identifier
    char partnerChain;    // Partner chain identifier
    int partnerResnum;    // Partner residue number
    char partnerInscode;  // Partner insertion code
    double bondLength;    // Length of the disulfide bond
};

/**
 * @brief Core residue data following CHARMM naming conventions
 */
struct ResidueData {
    // CHARMM standard fields
    std::string resname;     ///< Residue name (RESNAME)
    int ires;                ///< Residue number (IRES)
    std::string segid;       ///< Segment ID (SEGID)
    int iseg;                ///< Segment number (ISEG)
    int igro;                ///< Group number (IGRO)
    char chain;              ///< Chain identifier
    char inscode;            ///< Insertion code
    int move;                ///< Movement flag (MOVE)
    int ignore;              ///< Ignore flag (IGNORE)
    int constrain;           ///< Constraint flag (CONSTRAIN)
    bool hetatm;             ///< HETATM flag

    // Scalar properties (SCA1-9)
    std::array<double, 9> scalar;  ///< User-defined scalar properties

    // Secondary structure information
    SecondaryStructure secStruct{SecondaryStructure::NONE};
    SheetStrand sheetInfo;
    SSBond ssbond;

    // Center of mass calculation and storage
    std::array<double, 3> com{0.0, 0.0, 0.0};  // Center of mass

    // Default constructor
    ResidueData() : 
        resname(""),      // Residue name (e.g., ALA)
        ires(0),         // Residue number (IRES)
        segid(""),       // Segment ID (SEGID)
        iseg(0),         // Segment number (ISEG)
        igro(0),         // Group number (IGRO)
        chain(' '),      // Chain ID
        inscode(' '),    // Insertion code
        move(1),         // Movement flag (MOVE)
        ignore(0),       // Ignore flag (IGNORE)
        constrain(0),    // Constraint flag (CONSTRAIN)
        hetatm(false)    // HETATM flag
    {}
};

} // namespace residue
} // namespace model
} // namespace pygcmc
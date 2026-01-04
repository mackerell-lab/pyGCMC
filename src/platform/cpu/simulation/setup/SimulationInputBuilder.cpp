#include "SimulationInputBuilder.hpp"
#include "../../../../io/parameters/InpParserGCMC.hpp"
#include "../../../../io/structure/PdbParserMain.hpp"
#include "../../../../io/topology/topParserMain.hpp"
#include "../../../../io/topology/psfParserMain.hpp"
#include "../../../../io/forcefield/ItpNonbondedParser.hpp"
#include "../../../../io/forcefield/PrmParserMain.hpp"
#include "../../../../io/topology/FragmentLibrary.hpp"
#include "../../../../system/molecular/MolecularCombiner.hpp"
#include "../../../../system/montecarlo/MCInitializer.hpp"
#include "../../../../system/log/LogMain.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <filesystem>
#include <unordered_set>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {
namespace setup {

using namespace pygcmc::model;
using namespace pygcmc::io;

SimulationInputBuilder::SimulationInputBuilder(const Config& config)
    : config_(config) {
}

SimulationInputBuilder::Result SimulationInputBuilder::build() {
    Result result;

    auto pushUnique = [](std::vector<std::string>& v, const std::string& s) {
        if (std::find(v.begin(), v.end(), s) == v.end()) {
            v.push_back(s);
        }
    };

    // If the initial PDB contains additional fragment residues not present in the TOP,
    // we split the structure into:
    // - framework residues (must match TOP exactly) for MCInitializer
    // - extra residues (must match loaded fragment templates) which are added to MCState later
    // This supports gcmc_gpu/opencl style inputs where the TOP describes the protein/framework,
    // while fragment ITPs define the GCMC species (which may already appear in the initial PDB).
    std::vector<std::shared_ptr<Residue>> extraStructureResidues;
    
    // Step 1: Parse INP file
    log("Parsing INP file: " + config_.inpFile);
    result.parameters = parseINP(config_.inpFile);
    if (!result.parameters) {
        throw std::runtime_error("Failed to parse INP file");
    }
    
    // Determine base directory for resolving relative paths
    std::filesystem::path baseDir = std::filesystem::path(config_.inpFile).parent_path();
    if (baseDir.empty()) {
        baseDir = std::filesystem::current_path();
    }
    
    const auto& fileInfo = result.parameters->get_file_info();
    io::ItpNonbondedParser::Result itpNonbonded;
    bool itpNonbondedLoaded = false;
    
    // Step 2: Load PDB structure if available
    std::shared_ptr<Structure> structure;
    if (config_.loadStructure && !fileInfo.input_pdb_file.empty()) {
        std::string pdbPath = resolveFilePath(fileInfo.input_pdb_file, baseDir);
        log("Loading PDB structure: " + pdbPath);
        try {
            structure = loadPDB(pdbPath);
            if (structure) {
                result.structureLoaded = true;
                log("Loaded " + std::to_string(structure->get_atoms().size()) + " atoms from PDB");
            }

            // Heuristic warning: detect likely 10x unit mismatch between PDB CRYST1 (Å)
            // and INP box_size (interpreted by inp_units conversion).
            if (result.parameters) {
                auto& basic = result.parameters->get_basic_info();
                auto& warnings = basic.inp_warnings;
                auto pushWarning = [&](const std::string& code, const std::string& message) {
                    const auto it = std::find_if(warnings.begin(), warnings.end(), [&](const auto& w) { return w.code == code; });
                    if (it == warnings.end()) {
                        warnings.push_back({code, message});
                    }
                };

                std::array<float, 3> crystA = {0.0f, 0.0f, 0.0f};
                bool hasCryst1 = false;
                {
                    std::ifstream ifs(pdbPath);
                    std::string line;
                    // Scan a small prefix only; CRYST1 appears at the top in well-formed PDBs.
                    for (int i = 0; i < 200 && std::getline(ifs, line); ++i) {
                        if (line.rfind("CRYST1", 0) != 0) continue;
                        std::istringstream iss(line);
                        std::string tag;
                        double a = 0.0, b = 0.0, c = 0.0;
                        if (iss >> tag >> a >> b >> c) {
                            crystA = {static_cast<float>(a), static_cast<float>(b), static_cast<float>(c)};
                            hasCryst1 = true;
                        }
                        break;
                    }
                }

                if (hasCryst1) {
                    const auto& box = result.parameters->get_space_info().box_size;  // already internal nm
                    const bool hasBox = (box[0] > 0.0f || box[1] > 0.0f || box[2] > 0.0f);
                    if (hasBox) {
                        const std::array<float, 3> crystNm = {0.1f * crystA[0], 0.1f * crystA[1], 0.1f * crystA[2]};
                        auto ratioIn = [&](float inpNm, float pdbNm) -> float {
                            return (pdbNm > 0.0f) ? (inpNm / pdbNm) : 1.0f;
                        };
                        const float rx = ratioIn(box[0], crystNm[0]);
                        const float ry = ratioIn(box[1], crystNm[1]);
                        const float rz = ratioIn(box[2], crystNm[2]);

                        const auto isFactor10 = [&](float r) {
                            return ((r > 8.0f && r < 12.0f) || (r > 0.08f && r < 0.12f));
                        };
                        if (isFactor10(rx) || isFactor10(ry) || isFactor10(rz)) {
                            std::ostringstream oss;
                            oss << "PDB CRYST1=" << crystA[0] << " " << crystA[1] << " " << crystA[2] << " Å, "
                                << "INP box_size=" << box[0] << " " << box[1] << " " << box[2] << " nm "
                                << "(ratio≈" << rx << "," << ry << "," << rz << "); likely a 10x unit mismatch.";
                            pushWarning("UNIT_MISMATCH_PDB_CRYST1_VS_INP_BOX_SIZE", oss.str());
                        }
                    }
                }
            }
        } catch (const std::exception& e) {
            log("ERROR: Failed to load PDB: " + std::string(e.what()));
            throw;
        }
    }
    
    // Step 3: Load TOP topology if available
    std::shared_ptr<Topology> topology;
    if (config_.loadTopology && !fileInfo.topology_file.empty()) {
        std::string topPath = resolveFilePath(fileInfo.topology_file, baseDir);
        log("Loading topology: " + topPath);
        try {
            topology = loadTopology(topPath);
            if (topology) {
                result.topologyLoaded = true;
                log("Loaded topology with " + std::to_string(topology->get_num_atoms()) + " atoms");
            }
        } catch (const std::exception& e) {
            log("ERROR: Failed to load topology: " + std::string(e.what()));
            throw;
        }
    }
    
    // Step 4: Load PAR force field parameters if available
    // NOTE: gcmc_gpu uses GROMACS .itp for nonbonded; existing PRMParser targets CHARMM .prm/.str.
    if (config_.loadParameters && !fileInfo.par_files.empty()) {
        log("Loading force field parameters");
        try {
            std::vector<std::string> resolvedPrmFiles;
            std::vector<std::string> resolvedItpFiles;
            std::vector<std::string> resolvedDefaultsFiles;

            for (const auto& parFile : fileInfo.par_files) {
                std::string resolved = resolveFilePath(parFile, baseDir);
                std::filesystem::path p(resolved);
                std::string ext = p.extension().string();
                std::transform(ext.begin(), ext.end(), ext.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

                if (ext == ".itp") {
                    resolvedItpFiles.push_back(resolved);
                } else {
                    resolvedPrmFiles.push_back(resolved);
                }
            }

            for (const auto& defaultsFile : fileInfo.itp_defaults_files) {
                if (defaultsFile.empty()) continue;
                resolvedDefaultsFiles.push_back(resolveFilePath(defaultsFile, baseDir));
            }

            if (!resolvedPrmFiles.empty()) {
                result.forceField = loadParameters(resolvedPrmFiles);
                if (result.forceField) {
                    const bool anyLJ = result.forceField->get_num_lj_params() > 0;
                    const bool anyNBFix = result.forceField->get_num_nbfix() > 0;
                    result.parametersLoaded = anyLJ || anyNBFix;
                    if (result.parametersLoaded) {
                        log("Loaded force field parameters from PRM/STR files");
                    } else {
                        log("Warning: No usable LJ/NBFIX parameters found in PRM/STR files");
                    }
                }
            }

            if (!resolvedItpFiles.empty() || !resolvedDefaultsFiles.empty()) {
                std::vector<std::string> itpInputs = resolvedItpFiles;
                itpInputs.insert(itpInputs.begin(), resolvedDefaultsFiles.begin(), resolvedDefaultsFiles.end());
                itpNonbonded = io::ItpNonbondedParser::parse_files(itpInputs);
                itpNonbondedLoaded = !itpNonbonded.atomTypes.empty() ||
                                     !itpNonbonded.nbfixOverrides.empty() ||
                                     !itpNonbonded.pairtypesOverrides.empty() ||
                                     itpNonbonded.defaults.present;
                if (itpNonbondedLoaded) {
                    log("Loaded nonbonded parameters from GROMACS ITP files");
                } else {
                    log("Warning: No atomtypes/pair overrides found in GROMACS ITP files");
                }
                if (result.parameters) {
                    auto& basic = result.parameters->get_basic_info();
                    basic.gromacs_defaults_present = itpNonbonded.defaults.present;
                    basic.gromacs_nbfunc = itpNonbonded.defaults.nbfunc;
                    basic.gromacs_comb_rule = itpNonbonded.defaults.combRule;
                    basic.gromacs_gen_pairs = itpNonbonded.defaults.genPairs;
                    basic.gromacs_fudge_lj = itpNonbonded.defaults.fudgeLJ;
                    basic.gromacs_fudge_qq = itpNonbonded.defaults.fudgeQQ;
                    basic.gromacs_gen_pairs_present = itpNonbonded.defaults.genPairsPresent;
                    basic.gromacs_fudge_present = itpNonbonded.defaults.fudgePresent;

                    if (itpNonbonded.defaults.genPairsPresent) {
                        std::string gen = itpNonbonded.defaults.genPairs;
                        std::transform(gen.begin(), gen.end(), gen.begin(),
                                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                        if (gen != "yes" && gen != "y" && gen != "1") {
                            pushUnique(basic.inp_keys_ignored, "gromacs_gen_pairs");
                        }
                    }
                    if (itpNonbonded.defaults.fudgePresent) {
                        if (std::abs(itpNonbonded.defaults.fudgeLJ - 1.0) > 1e-6) {
                            pushUnique(basic.inp_keys_ignored, "gromacs_fudge_lj");
                        }
                        if (std::abs(itpNonbonded.defaults.fudgeQQ - 1.0) > 1e-6) {
                            pushUnique(basic.inp_keys_ignored, "gromacs_fudge_qq");
                        }
                    }
                }
            }
        } catch (const std::exception& e) {
            log("Warning: Failed to load parameters: " + std::string(e.what()));
        }
    }
    
    // Step 4b: Load fragment templates if available
    if (!fileInfo.fragment_top_files.empty()) {
        log("Loading fragment templates");
        std::vector<std::string> resolvedFragFiles;
        for (const auto& fragFile : fileInfo.fragment_top_files) {
            resolvedFragFiles.push_back(resolveFilePath(fragFile, baseDir));
        }
        result.fragmentTemplates = loadFragmentTemplates(resolvedFragFiles, 
                                                         result.parameters,
                                                         baseDir);
        log("Loaded " + std::to_string(result.fragmentTemplates.size()) + " fragment templates");
    }
    
    // Step 5: Combine molecular data if we have BOTH structure and topology with actual data
    if (result.structureLoaded && result.topologyLoaded) {
        log("Combining molecular data");
        // Use already loaded structure instead of re-loading
        
        // Only combine if we have actual atoms
        if (structure && structure->get_atoms().size() > 0 && 
            topology && topology->get_num_atoms() > 0) {
            try {
                result.molecular = combineMolecular(structure, topology, result.forceField);
            } catch (const std::exception& e) {
                // Compatibility: allow extra residues in the PDB (e.g., initial guests/solvent)
                // when they can be mapped to fragment templates (fragitp).
                if (result.fragmentTemplates.empty()) {
                    throw;
                }

                auto lower = [](std::string s) {
                    std::transform(s.begin(), s.end(), s.begin(),
                                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                    return s;
                };
                auto residueKey = [&](const std::string& name, int number) {
                    return lower(name) + ":" + std::to_string(number);
                };

                std::unordered_set<std::string> topoResidues;
                topoResidues.reserve(static_cast<size_t>(topology->get_num_residues()));
                for (int i = 0; i < topology->get_num_residues(); ++i) {
                    const auto& tr = topology->get_residue(i);
                    topoResidues.insert(residueKey(tr.name, tr.number));
                }

                std::unordered_set<std::string> templateNames;
                templateNames.reserve(result.fragmentTemplates.size());
                for (const auto& [name, tmpl] : result.fragmentTemplates) {
                    (void)tmpl;
                    templateNames.insert(lower(name));
                }

                auto filtered = std::make_shared<Structure>();
                filtered->set_box_dimensions(structure->get_box_dimensions());
                extraStructureResidues.clear();

                for (const auto& res : structure->get_residues()) {
                    if (!res) continue;
                    const std::string key = residueKey(res->get_resname(), res->get_ires());
                    if (topoResidues.find(key) != topoResidues.end()) {
                        filtered->add_residue(res);
                        for (const auto& atom : res->get_atoms()) {
                            if (atom) filtered->add_atom(atom);
                        }
                    } else {
                        extraStructureResidues.push_back(res);
                    }
                }

                // Only apply this fallback when the structure contains ALL topology residues
                // plus additional residues that can be mapped to fragment templates.
                if (filtered->get_residues().size() != static_cast<size_t>(topology->get_num_residues()) ||
                    extraStructureResidues.empty()) {
                    throw;
                }
                for (const auto& res : extraStructureResidues) {
                    const std::string nm = lower(res->get_resname());
                    if (templateNames.find(nm) == templateNames.end()) {
                        std::ostringstream oss;
                        oss << "Inconsistent structure/topology: PDB contains residue '" << res->get_resname()
                            << "' (ires=" << res->get_ires()
                            << ") which is not present in the topology and has no matching fragment template.";
                        throw std::runtime_error(oss.str());
                    }
                }

                log("Info: PDB contains residues beyond the TOP; loading framework from TOP and seeding "
                    "initial fragment residues from PDB.");
                result.molecular = combineMolecular(filtered, topology, result.forceField);
            }
        } else {
            log("Warning: Structure or topology is empty, skipping molecular combination");
            result.structureLoaded = false;  // Mark as not loaded if empty
            result.topologyLoaded = false;
        }
    }
    
    // Step 6: Initialize MC state
    log("Initializing MC state");
    result.mcState = std::make_shared<montecarlo::MCState>();
    
    // Use MCInitializer if we have complete molecular data with actual atoms.
    // NOTE: This should work even when force field parameters are provided as GROMACS ITP
    // (via ItpNonbondedParser) and no CHARMM ForceField object is available.
    bool useMolecularInit = false;
    if (result.molecular && result.structureLoaded && result.topologyLoaded &&
        result.molecular->get_num_atoms() > 0) {
        auto hasConsistentTopologyMapping = [&]() -> bool {
            const size_t nRes = result.molecular->get_num_residues();
            if (nRes == 0) return false;
            if (result.molecular->residues.size() < nRes) return false;
            if (result.molecular->topology_residues.size() < nRes) return false;

            const size_t nTopAtoms = result.molecular->topology_atoms.size();
            for (size_t ri = 0; ri < nRes; ++ri) {
                const auto& molRes = result.molecular->residues[ri];
                if (!molRes) return false;
                const auto& topRes = result.molecular->topology_residues[ri];
                const auto& molAtoms = molRes->get_atoms();
                if (topRes.atoms.size() != molAtoms.size()) return false;
                for (int idx : topRes.atoms) {
                    if (idx < 0 || static_cast<size_t>(idx) >= nTopAtoms) return false;
                }
            }
            return true;
        };

        useMolecularInit = hasConsistentTopologyMapping();
        if (!useMolecularInit) {
            log("Warning: Incomplete topology mapping for structure residues; "
                "skipping MCInitializer and using INP-only initialization");
        }
    }

    if (useMolecularInit) {
        try {
            // Pre-allocate MCState arrays for MCCore::addInitialResidues() capacity checks.
            const size_t initialAtoms = result.molecular->get_num_atoms();
            const size_t initialResidues = result.molecular->get_num_residues();

            result.mcState->info.maxAtoms = static_cast<int>(initialAtoms);
            result.mcState->info.maxResidues = static_cast<int>(initialResidues);
            result.mcState->atoms.resize(initialAtoms);
            result.mcState->residues.resize(initialResidues);

            // Use MCInitializer to populate from real molecular data
            system::montecarlo::MCInitializer initializer;
            initializer.initializeFromMolecular(*result.mcState, result.molecular);

            // Set beta/cutoff from INP (units already normalized by InpParserGCMC).
            const auto& mcInfo = result.parameters->get_mc_info();
            result.mcState->info.beta = mcInfo.beta;
            const auto& spaceInfo = result.parameters->get_space_info();
            result.mcState->info.cutoff = spaceInfo.cutoff;

            // Periodic box semantics (unit-mode aware):
            // - inp_units == nm: INP box_size is the intended periodic box; CRYST1 in minimal/legacy PDBs may be a stub.
            // - inp_units == gcmc_gpu: PDB CRYST1 is the periodic box; INP box_size commonly describes the *GCMC region*
            //   (active box / grid) rather than the periodic unit cell.
            const bool inpBoxProvided = (spaceInfo.box_size[0] > 0.0f ||
                                         spaceInfo.box_size[1] > 0.0f ||
                                         spaceInfo.box_size[2] > 0.0f);
            const bool pdbBoxProvided = (result.mcState->info.box[0] > 0.0f &&
                                         result.mcState->info.box[1] > 0.0f &&
                                         result.mcState->info.box[2] > 0.0f);

            auto toLower = [](std::string s) {
                std::transform(s.begin(), s.end(), s.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                return s;
            };

            bool preferPdbBox = false;
            if (result.parameters) {
                const std::string u = toLower(result.parameters->get_basic_info().inp_units);
                preferPdbBox = (u == "gcmc_gpu" || u == "charmm" ||
                                u == "angstrom" || u == "ang" || u == "a" ||
                                u == "a_kcal" || u == "a/kcal" || u == "akcal");
            }

            if (!preferPdbBox) {
                // nm/native: prefer INP box_size when present.
                if (inpBoxProvided) {
                    result.mcState->setBoxDimensions(spaceInfo.box_size[0],
                                                     spaceInfo.box_size[1],
                                                     spaceInfo.box_size[2]);
                    log("Set periodic box from INP box_size: " +
                        std::to_string(spaceInfo.box_size[0]) + " x " +
                        std::to_string(spaceInfo.box_size[1]) + " x " +
                        std::to_string(spaceInfo.box_size[2]) + " nm");
                } else if (pdbBoxProvided) {
                    result.mcState->setBoxDimensions(result.mcState->info.box[0],
                                                     result.mcState->info.box[1],
                                                     result.mcState->info.box[2]);
                } else {
                    // Ensure periodicBox exists even if box is unset.
                    result.mcState->setBoxDimensions(result.mcState->info.box[0],
                                                     result.mcState->info.box[1],
                                                     result.mcState->info.box[2]);
                }
            } else {
                // gcmc_gpu/opencl: prefer CRYST1 when present.
                if (pdbBoxProvided) {
                    result.mcState->setBoxDimensions(result.mcState->info.box[0],
                                                     result.mcState->info.box[1],
                                                     result.mcState->info.box[2]);
                } else if (inpBoxProvided) {
                    result.mcState->setBoxDimensions(spaceInfo.box_size[0],
                                                     spaceInfo.box_size[1],
                                                     spaceInfo.box_size[2]);
                    log("Set periodic box from INP box_size (no CRYST1): " +
                        std::to_string(spaceInfo.box_size[0]) + " x " +
                        std::to_string(spaceInfo.box_size[1]) + " x " +
                        std::to_string(spaceInfo.box_size[2]) + " nm");
                } else {
                    // Ensure periodicBox is initialized (some movement/energy code indexes it directly).
                    result.mcState->setBoxDimensions(result.mcState->info.box[0],
                                                     result.mcState->info.box[1],
                                                     result.mcState->info.box[2]);
                }
            }

            // Mark pre-existing residues as fixed (protein) unless they match INP fragment names.
            // This supports cavity-bias "exclude protein" and prevents protein from being moved/deleted.
            if (result.parameters) {
                std::unordered_set<std::string> fragmentNamesLower;
                for (const auto& n : result.parameters->get_file_info().fragment_names) {
                    std::string key = n;
                    std::transform(key.begin(), key.end(), key.begin(),
                                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                    fragmentNamesLower.insert(key);
                }

                size_t fixedCount = 0;
                size_t movableCount = 0;
                const int nRes = std::min(result.mcState->activeResidueCount,
                                          static_cast<int>(result.mcState->residues.size()));
                for (int i = 0; i < nRes; ++i) {
                    auto& res = result.mcState->residues[i];
                    std::string resName = res.resname;
                    if (resName.empty()) {
                        resName = result.mcState->residueTypes.getTypeName(res.type);
                    }
                    std::transform(resName.begin(), resName.end(), resName.begin(),
                                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                    const bool isFragment = fragmentNamesLower.find(resName) != fragmentNamesLower.end();
                    res.fixed = !isFragment;
                    if (res.fixed) {
                        ++fixedCount;
                    } else {
                        ++movableCount;
                    }
                }
                log("Residue fixed flags: fixed=" + std::to_string(fixedCount) +
                    " movable=" + std::to_string(movableCount));
            }

            if (result.forceField) {
                initializer.initializeForceField(*result.mcState, *result.forceField);
            }

            log("MC state initialized with real molecular data");
        } catch (const std::exception& e) {
            log(std::string("Warning: Failed to initialize MC state from molecular data: ") +
                e.what() + " (falling back to INP-only initialization)");
            useMolecularInit = false;
        }
    }

    if (!useMolecularInit) {
        // Fall back to empty state with box dimensions from INP
        const auto& spaceInfo = result.parameters->get_space_info();
        result.mcState->setBoxDimensions(spaceInfo.box_size[0], spaceInfo.box_size[1], spaceInfo.box_size[2]);
        
        // Set temperature
        const auto& mcInfo = result.parameters->get_mc_info();
        result.mcState->info.beta = mcInfo.beta;

        // Set cutoff from space info
        result.mcState->info.cutoff = spaceInfo.cutoff;

        // Also sync cutoff to energy parameters if parameters exist
        if (result.parameters) {
            auto& energyInfo = const_cast<model::param::EnergyInfo&>(result.parameters->get_energy_info());
            energyInfo.fragment_cutoff = spaceInfo.cutoff;
            energyInfo.protein_cutoff = spaceInfo.cutoff;
            energyInfo.fragment_cutoff_squared = spaceInfo.cutoff * spaceInfo.cutoff;
            energyInfo.protein_cutoff_squared = spaceInfo.cutoff * spaceInfo.cutoff;
        }

        log("MC state initialized with INP parameters only");
    }

    // Enforce minimum-image convention: cutoff must be strictly less than half the smallest box dimension.
    // This prevents "runs but silently wrong" truncated interactions under PBC.
    if (result.mcState) {
        const float cutoff = result.mcState->info.cutoff;
        const float bx = result.mcState->info.box[0];
        const float by = result.mcState->info.box[1];
        const float bz = result.mcState->info.box[2];
        const float minBox = std::min(bx, std::min(by, bz));
        if (cutoff > 0.0f && minBox > 0.0f && cutoff > 0.5f * minBox) {
            std::ostringstream oss;
            oss << "Invalid cutoff: cutoff=" << cutoff << " nm must be <= 0.5 * min(box) where box="
                << bx << " " << by << " " << bz << " nm.";
            throw std::runtime_error(oss.str());
        }
    }

    // Step 6b: Apply GROMACS ITP nonbonded parameters to MC state if available.
    // Only apply when the MCState force field is still empty to avoid overriding a fully
    // initialized CHARMM (PRM/STR) force field path.
    if (itpNonbondedLoaded && !itpNonbonded.atomTypes.empty() && result.mcState &&
        (result.mcState->forcefield.numTotalTypes == 0 || result.mcState->forcefield.ljSigma.empty())) {
        // Deterministic type indexing (std::map iteration order).
        for (const auto& [typeName, lj] : itpNonbonded.atomTypes) {
            (void)lj;
            result.mcState->getOrAddAtomType(typeName);
        }

        const int n = static_cast<int>(result.mcState->atomTypes.atomTypes.size());
        result.mcState->forcefield.numTotalTypes = n;
        result.mcState->forcefield.numMovementTypes = n;
        auto mixingRule = montecarlo::MCForceField::MixingRule::LorentzBerthelot;
        if (itpNonbonded.defaults.combRule == 1 || itpNonbonded.defaults.combRule == 3) {
            mixingRule = montecarlo::MCForceField::MixingRule::Geometric;
        }
        result.mcState->forcefield.mixingRule = mixingRule;
        result.mcState->forcefield.ljSigmaType.assign(static_cast<size_t>(n), 0.0f);
        result.mcState->forcefield.ljEpsType.assign(static_cast<size_t>(n), 0.0f);
        result.mcState->forcefield.nbfix.clear();
        result.mcState->forcefield.clearPairtypes14();

        for (const auto& [typeName, lj] : itpNonbonded.atomTypes) {
            auto it = result.mcState->atomTypes.atomTypeIndices.find(typeName);
            if (it == result.mcState->atomTypes.atomTypeIndices.end()) {
                continue;
            }
            const int idx = it->second;
            if (idx >= 0 && idx < n) {
                result.mcState->forcefield.ljSigmaType[static_cast<size_t>(idx)] =
                    static_cast<float>(lj.sigma_nm);
                result.mcState->forcefield.ljEpsType[static_cast<size_t>(idx)] =
                    static_cast<float>(lj.epsilon_kj);
            }
        }

        const auto& basicInfo = result.parameters->get_basic_info();
        std::string pairtypesMode = basicInfo.itp_pairtypes_mode;
        std::transform(pairtypesMode.begin(), pairtypesMode.end(), pairtypesMode.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        const bool strictPairtypes = (pairtypesMode == "strict");

        if (!strictPairtypes) {
            for (const auto& [pair, lj] : itpNonbonded.pairtypesOverrides) {
                auto it1 = result.mcState->atomTypes.atomTypeIndices.find(pair.first);
                auto it2 = result.mcState->atomTypes.atomTypeIndices.find(pair.second);
                if (it1 == result.mcState->atomTypes.atomTypeIndices.end() ||
                    it2 == result.mcState->atomTypes.atomTypeIndices.end()) {
                    continue;
                }
                result.mcState->forcefield.addNBFix(
                    it1->second,
                    it2->second,
                    static_cast<float>(lj.sigma_nm),
                    static_cast<float>(lj.epsilon_kj)
                );
            }
        }

        for (const auto& [pair, lj] : itpNonbonded.nbfixOverrides) {
            auto it1 = result.mcState->atomTypes.atomTypeIndices.find(pair.first);
            auto it2 = result.mcState->atomTypes.atomTypeIndices.find(pair.second);
            if (it1 == result.mcState->atomTypes.atomTypeIndices.end() ||
                it2 == result.mcState->atomTypes.atomTypeIndices.end()) {
                continue;
            }
            result.mcState->forcefield.addNBFix(
                it1->second,
                it2->second,
                static_cast<float>(lj.sigma_nm),
                static_cast<float>(lj.epsilon_kj)
            );
        }

        if (strictPairtypes) {
            for (const auto& [pair, lj] : itpNonbonded.pairtypesOverrides) {
                auto it1 = result.mcState->atomTypes.atomTypeIndices.find(pair.first);
                auto it2 = result.mcState->atomTypes.atomTypeIndices.find(pair.second);
                if (it1 == result.mcState->atomTypes.atomTypeIndices.end() ||
                    it2 == result.mcState->atomTypes.atomTypeIndices.end()) {
                    continue;
                }
                result.mcState->forcefield.setPairtype14(
                    it1->second,
                    it2->second,
                    static_cast<float>(lj.sigma_nm),
                    static_cast<float>(lj.epsilon_kj)
                );
            }
        }

        result.mcState->forcefield.rebuildLJMatrix();
        log("Applied GROMACS ITP nonbonded parameters to MCState (types=" + std::to_string(n) + ")");
    }

    // Step 6c: If the initial PDB contained extra residues (not present in TOP),
    // seed them into the MCState as fragment instances using the loaded fragment templates.
    if (!extraStructureResidues.empty() && result.mcState && !result.fragmentTemplates.empty()) {
        auto lower = [](std::string s) {
            std::transform(s.begin(), s.end(), s.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return s;
        };

        std::unordered_map<std::string, const platform::cpu::movement::FragmentTemplate*> tmplByLower;
        tmplByLower.reserve(result.fragmentTemplates.size());
        for (const auto& [name, tmpl] : result.fragmentTemplates) {
            tmplByLower.emplace(lower(name), &tmpl);
        }

        size_t seededFragments = 0;
        for (const auto& pdbRes : extraStructureResidues) {
            if (!pdbRes) continue;

            const auto itT = tmplByLower.find(lower(pdbRes->get_resname()));
            if (itT == tmplByLower.end() || !itT->second) {
                continue;
            }
            const auto& tmpl = *itT->second;
            if (tmpl.atoms.empty()) {
                throw std::runtime_error("Fragment template has no atoms for: " + pdbRes->get_resname());
            }

            const auto& pdbAtoms = pdbRes->get_atoms();
            if (pdbAtoms.empty()) {
                continue;
            }

            std::vector<size_t> pdbIndexForTemplate;
            pdbIndexForTemplate.resize(tmpl.atoms.size(), 0);
            bool mapped = false;

            if (pdbAtoms.size() == tmpl.atoms.size()) {
                for (size_t i = 0; i < tmpl.atoms.size(); ++i) {
                    pdbIndexForTemplate[i] = i;
                }
                mapped = true;
            } else {
                // Name-based sequential matching fallback.
                size_t p = 0;
                mapped = true;
                for (size_t i = 0; i < tmpl.atoms.size(); ++i) {
                    const std::string want = tmpl.atoms[i].name;
                    while (p < pdbAtoms.size() && pdbAtoms[p] && pdbAtoms[p]->get_type() != want) {
                        ++p;
                    }
                    if (p >= pdbAtoms.size()) {
                        mapped = false;
                        break;
                    }
                    pdbIndexForTemplate[i] = p;
                    ++p;
                }
            }

            if (!mapped) {
                std::ostringstream oss;
                oss << "Failed to map PDB residue '" << pdbRes->get_resname() << "' (ires=" << pdbRes->get_ires()
                    << ") atoms onto fragment template '" << tmpl.name << "' (template_atoms=" << tmpl.atoms.size()
                    << ", pdb_atoms=" << pdbAtoms.size() << ").";
                throw std::runtime_error(oss.str());
            }

            model::montecarlo::MCResidue mcRes;
            mcRes.active = true;
            mcRes.fixed = false;
            mcRes.resname = pdbRes->get_resname();
            mcRes.resid = pdbRes->get_ires();
            mcRes.type = result.mcState->residueTypes.getOrAddType(mcRes.resname);
            mcRes.radius = static_cast<float>(tmpl.radius);

            mcRes.atomStart = result.mcState->activeAtomCount;
            mcRes.atomCount = static_cast<int>(tmpl.atoms.size());
            mcRes.atoms.clear();
            mcRes.atoms.reserve(tmpl.atoms.size());

            double cx = 0.0, cy = 0.0, cz = 0.0;
            for (size_t i = 0; i < tmpl.atoms.size(); ++i) {
                const size_t pdbIdx = pdbIndexForTemplate[i];
                if (pdbIdx >= pdbAtoms.size() || !pdbAtoms[pdbIdx]) {
                    throw std::runtime_error("Invalid PDB atom mapping for residue: " + pdbRes->get_resname());
                }
                const auto& coorA = pdbAtoms[pdbIdx]->get_coor();

                model::montecarlo::MCAtom atom;
                atom.name = tmpl.atoms[i].name;
                atom.charge = tmpl.atoms[i].charge;
                atom.mass = tmpl.atoms[i].mass;

                // Map ITP type name -> MCState type index.
                if (!tmpl.atomTypeNames.empty() && tmpl.atomTypeNames.size() == tmpl.atoms.size()) {
                    const std::string& typeName = tmpl.atomTypeNames[i];
                    const auto it = result.mcState->atomTypes.atomTypeIndices.find(typeName);
                    if (it == result.mcState->atomTypes.atomTypeIndices.end()) {
                        std::ostringstream oss;
                        oss << "Missing atom type '" << typeName << "' in MCState while seeding initial residue '"
                            << pdbRes->get_resname() << "' (ires=" << pdbRes->get_ires() << ").";
                        throw std::runtime_error(oss.str());
                    }
                    atom.type = it->second;
                } else {
                    atom.type = tmpl.atoms[i].type;
                }

                atom.x = static_cast<float>(coorA[0] * 0.1);
                atom.y = static_cast<float>(coorA[1] * 0.1);
                atom.z = static_cast<float>(coorA[2] * 0.1);
                atom.updatePosition();

                mcRes.atoms.push_back(atom);
                if (result.mcState->activeAtomCount < static_cast<int>(result.mcState->atoms.size())) {
                    result.mcState->atoms[result.mcState->activeAtomCount] = atom;
                } else {
                    result.mcState->atoms.push_back(atom);
                }
                result.mcState->activeAtomCount++;

                cx += atom.x;
                cy += atom.y;
                cz += atom.z;
            }

            const double inv = 1.0 / static_cast<double>(tmpl.atoms.size());
            mcRes.center[0] = static_cast<float>(cx * inv);
            mcRes.center[1] = static_cast<float>(cy * inv);
            mcRes.center[2] = static_cast<float>(cz * inv);

            result.mcState->residues.push_back(std::move(mcRes));
            result.mcState->activeResidueCount = static_cast<int>(result.mcState->residues.size());
            seededFragments++;
        }

        log("Seeded " + std::to_string(seededFragments) + " initial fragment residues from PDB into MCState");
    }
    
    return result;
}

std::map<std::string, platform::cpu::movement::FragmentTemplate> SimulationInputBuilder::loadFragmentTemplates(
    const std::vector<std::string>& fragItpFiles,
    const std::shared_ptr<model::param::Param>& parameters,
    const std::filesystem::path& baseDir) {
    
    std::map<std::string, platform::cpu::movement::FragmentTemplate> templates;
    
    if (!parameters) {
        throw std::runtime_error("Parameters not available for fragment template loading");
    }
    
    // Get fragment info from parameters for matching
    const auto& fragNames = parameters->get_file_info().fragment_names;
    const auto& fragConcs = parameters->get_fragment_info().conc_list;
    const auto& fragMuexs = parameters->get_fragment_info().muex_list;
    const auto& monomerDir = parameters->get_file_info().monomer_dir;

    std::filesystem::path monomerDirResolved;
    if (!monomerDir.empty()) {
        monomerDirResolved = std::filesystem::path(resolveFilePath(monomerDir, baseDir));
    }
    
    for (size_t i = 0; i < fragItpFiles.size(); ++i) {
        const auto& itpFile = fragItpFiles[i];
        log("Loading fragment template from: " + itpFile);
        
        // Check if file exists
        std::ifstream file(itpFile);
        if (!file.good()) {
            log("ERROR: Fragment template not found: " + itpFile);
            throw std::runtime_error("Fragment template not found: " + itpFile);
        }
        file.close();
        
        // Resolve fragment name: prefer INP fragname list (aligned by index with fragitp),
        // fallback to ITP filename stem if missing.
        std::filesystem::path itpPath(itpFile);
        const std::string itpStem = itpPath.stem().string();
        std::string fragName = (i < fragNames.size() && !fragNames[i].empty()) ? fragNames[i] : itpStem;

        // Try to find fragment coordinate PDB from monomerdir/<frag>.pdb (gcmc_gpu behavior).
        // Many inputs use fragname that differs from the ITP filename stem (e.g., WAT vs sol.itp),
        // so we try multiple name candidates deterministically.
        auto lower = [](std::string s) {
            std::transform(s.begin(), s.end(), s.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return s;
        };
        const std::vector<std::string> pdbNameCandidates = {
            fragName,
            lower(fragName),
            itpStem,
            lower(itpStem),
        };
        std::filesystem::path coordPdb;
        for (const auto& n : pdbNameCandidates) {
            if (n.empty()) continue;
            std::filesystem::path candidate;
            if (!monomerDirResolved.empty()) {
                candidate = monomerDirResolved / (n + ".pdb");
            } else {
                candidate = itpPath.parent_path() / (n + ".pdb");
            }
            if (std::filesystem::exists(candidate)) {
                coordPdb = candidate;
                break;
            }
        }
        const std::string coordPdbPath = !coordPdb.empty() ? coordPdb.string() : "";
        
        // Use FragmentLibrary to load the ITP file
        io::topology::FragmentLibrary fragLib;
        if (!fragLib.loadFromITP(itpFile, fragName, templates.size(), coordPdbPath)) {
            throw std::runtime_error("Failed to parse fragment template: " + itpFile);
        }
        
        // Get the fragment data
        auto fragmentData = fragLib.get(fragName);
        if (!fragmentData) {
            throw std::runtime_error("Failed to retrieve fragment data for: " + fragName);
        }
        
        // Build complete FragmentTemplate
        platform::cpu::movement::FragmentTemplate tmpl;
        tmpl.name = fragmentData->name;
        tmpl.typeId = fragmentData->typeId;
        tmpl.atoms = fragmentData->atoms;
        tmpl.atomTypeNames = fragmentData->atomTypeNames;
        tmpl.molecularWeight = fragmentData->molecularWeight;
        tmpl.radius = fragmentData->radius;
        tmpl.bonds.clear();
        tmpl.bonds.reserve(fragmentData->bonds.size());
        for (const auto& b : fragmentData->bonds) {
            platform::cpu::movement::FragmentTemplate::Bond bond;
            bond.atom1 = b.atom1;
            bond.atom2 = b.atom2;
            bond.length = 0.0;
            tmpl.bonds.push_back(bond);
        }
        
        // Assign concentration and chemical potential by index (fragitp order).
        if (i < fragConcs.size()) {
            tmpl.concentration = fragConcs[i];
        }
        if (i < fragMuexs.size()) {
            tmpl.chemicalPotential = fragMuexs[i];
        }
        // Calculate activity: z = exp(β*μ) (acceptance uses its own activity model; this is for completeness)
        const double beta = parameters->get_mc_info().beta;
        tmpl.activity = std::exp(beta * tmpl.chemicalPotential);
        
        templates[tmpl.name] = tmpl;
        log("Loaded fragment " + tmpl.name + " with " + 
            std::to_string(tmpl.atoms.size()) + " atoms, " +
            "conc=" + std::to_string(tmpl.concentration) + " M, " +
            "μ=" + std::to_string(tmpl.chemicalPotential) + " kJ/mol");
    }
    
    if (templates.empty() && !fragItpFiles.empty()) {
        throw std::runtime_error("Failed to load any fragment templates");
    }
    
    return templates;
}

std::shared_ptr<param::Param> SimulationInputBuilder::parseINP(const std::string& filename) {
    auto params = std::make_shared<param::Param>();
    parameters::InpParserGCMC::parse_to_param(filename, *params);
    return params;
}

std::shared_ptr<model::Structure> SimulationInputBuilder::loadPDB(const std::string& filename) {
    // Check if file exists
    std::ifstream file(filename);
    if (!file.good()) {
        throw std::runtime_error("PDB file not found: " + filename);
    }
    file.close();
    
    auto structure = std::make_shared<model::Structure>();
    io::structure::PdbParserMain::parse_to_structure(filename, *structure);
    return structure;
}

std::shared_ptr<model::Topology> SimulationInputBuilder::loadTopology(const std::string& filename) {
    // Check if file exists
    std::ifstream file(filename);
    if (!file.good()) {
        throw std::runtime_error("Topology file not found: " + filename);
    }
    file.close();
    
    auto topology = std::make_shared<model::Topology>();

    std::string ext = std::filesystem::path(filename).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (ext == ".psf") {
        io::PSFParser parser;
        if (!parser.parse_to_topology(filename, *topology)) {
            throw std::runtime_error("Failed to parse PSF topology file: " + filename);
        }
    } else {
        io::TOPParser parser;
        if (!parser.parse_to_topology(filename, *topology)) {
            throw std::runtime_error("Failed to parse topology file: " + filename);
        }
    }
    return topology;
}

std::shared_ptr<model::ForceField> SimulationInputBuilder::loadParameters(
    const std::vector<std::string>& parFiles) {
    
    auto forceField = std::make_shared<model::ForceField>();
    
    for (const auto& parFile : parFiles) {
        try {
            log("Loading parameters from: " + parFile);
            io::PRMParser::parse_file_to_forcefield(parFile, *forceField);
        } catch (const std::exception& e) {
            log("Warning: Failed to load " + parFile + ": " + e.what());
        }
    }
    
    return forceField;
}

std::shared_ptr<model::Molecular> SimulationInputBuilder::combineMolecular(
    const std::shared_ptr<model::Structure>& structure,
    const std::shared_ptr<model::Topology>& topology,
    const std::shared_ptr<model::ForceField>& /*forceField*/) {
    
    system::molecular::MolecularCombiner combiner;
    
    // Create molecular system
    auto molecular = std::make_shared<model::Molecular>();
    
    if (structure && topology) {
        // Combine structure and topology
        molecular = combiner.combine(structure, topology);
    } else if (structure) {
        // Structure only - copy atoms and residues
        molecular->atoms = structure->get_atoms();
        molecular->residues = structure->get_residues();
    } else if (topology) {
        // Topology only - copy available topology data
        // Note: Topology class doesn't expose atoms/residues directly
        molecular->bonds = topology->get_bonds();
        molecular->angles = topology->get_angles();
        molecular->dihedrals = topology->get_dihedrals();
    }
    
    // Note: ForceField is separate and not stored in Molecular
    
    return molecular;
}

std::shared_ptr<model::MCState> SimulationInputBuilder::initializeMCState(
    const std::shared_ptr<model::Molecular>& molecular,
    const std::shared_ptr<model::ForceField>& forceField) {
    
    auto mcState = std::make_shared<model::MCState>();
    
    if (molecular) {
        system::montecarlo::MCInitializer initializer;
        initializer.initializeFromMolecular(*mcState, molecular);
        
        if (forceField) {
            initializer.initializeForceField(*mcState, *forceField);
        }
    }
    
    return mcState;
}

std::string SimulationInputBuilder::resolveFilePath(const std::string& path, 
                                                    const std::filesystem::path& baseDir) const {
    std::filesystem::path filePath(path);
    
    // If path is already absolute, return as-is
    if (filePath.is_absolute()) {
        return path;
    }
    
    // Otherwise resolve relative to base directory
    std::filesystem::path resolved = baseDir / filePath;
    return resolved.string();
}

void SimulationInputBuilder::log(const std::string& message) const {
    if (config_.verbose) {
        system::log::LogMain::info("[SimulationInputBuilder] ", message);
    }
}

} // namespace setup
} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc

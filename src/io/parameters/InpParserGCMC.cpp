#include "InpParserGCMC.hpp"
#include "../../model/param/ParamOperations.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace pygcmc {
namespace io {
namespace parameters {

using ::pygcmc::io::parameters::InpParserStructures;

void InpParserGCMC::parse_to_param(const std::string& filename, model::param::Param& param) {
    // First parse with the base parser
    InpParserMain::parse_to_param(filename, param);

    // Then parse GCMC-specific keys
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open input file: " + filename);
    }
    std::string line;
    while (std::getline(file, line)) {
        line = InpParserStructures::trim(line);
        if (line.empty() || line[0] == '#') continue;
        auto tokens = InpParserStructures::split(line, ':');
        if (tokens.size() != 2) continue;
        const std::string key = InpParserStructures::trim(tokens[0]);
        const std::string value = InpParserStructures::trim(tokens[1]);
        parse_line_ext(key, value, param);
    }
    enhance_param(param);
}

void InpParserGCMC::parse_line_ext(const std::string& key, const std::string& value, model::param::Param& param) {
    auto& mc_info = param.get_mc_info();
    auto& frag_info = param.get_fragment_info();
    auto& bias_info = param.get_bias_info();
    auto& space_info = param.get_space_info();
    auto& energy_info = param.get_energy_info();
    auto& file_info = param.get_file_info();
    auto& basic_info = param.get_basic_info();

    if (key == "mctime" || key == "mc_time") {
        // Support accumulation of multiple mctime lines
        auto times = InpParserStructures::parse_float_vector(value);
        for (float t : times) {
            mc_info.mc_time_list.push_back(t);
        }
    } else if (key == "energy_cutoff") {
        // Legacy gcmc_opencl key: treat as a shared cutoff for fragment/protein.
        float cutoff = std::stof(value);
        if (cutoff > 0.0f) {
            energy_info.fragment_cutoff = cutoff;
            energy_info.fragment_cutoff_squared = cutoff * cutoff;
            energy_info.protein_cutoff = cutoff;
            energy_info.protein_cutoff_squared = cutoff * cutoff;
            // Mirror into the global cutoff so that energy calculations stay consistent.
            space_info.cutoff = cutoff;
        }
    } else if (key == "energy_cutoff_frag" || key == "energy_cutoff_fragment") {
        float cutoff = std::stof(value);
        if (cutoff > 0.0f) {
            energy_info.fragment_cutoff = cutoff;
            energy_info.fragment_cutoff_squared = cutoff * cutoff;
            // Mirror fragment cutoff into the global cutoff so that energy calculations stay consistent
            space_info.cutoff = cutoff;
        }
    } else if (key == "energy_cutoff_prot" || key == "energy_cutoff_protein") {
        float cutoff = std::stof(value);
        if (cutoff > 0.0f) {
            energy_info.protein_cutoff = cutoff;
            energy_info.protein_cutoff_squared = cutoff * cutoff;
        }
    } else if (key == "fragradius") {
        // Raw value; normalized to internal nm in enhance_param.
        frag_info.radius_list = InpParserStructures::parse_float_vector(value);
    } else if (key == "fragconf" || key == "fragconfs") {
        // keep both for compatibility
        frag_info.conf_list = InpParserStructures::parse_int_vector(value);
        frag_info.fragconf_list = frag_info.conf_list;
    } else if (key == "num_conf_bias_trial" || key == "confbias_trials") {
        bias_info.num_conf_bias_trials = static_cast<unsigned int>(std::stoi(value));
    } else if (key == "cavity_grid_dx") {
        // Map to grid_dx if given (fallback)
        // Raw value; normalized to internal nm in enhance_param.
        space_info.grid_spacing = std::stof(value);
    } else if (key == "cavity_grid_dx_frag") {
        frag_info.cavity_grid_dx_list = InpParserStructures::parse_float_vector(value);
    } else if (key == "probe_radius") {
        // Map to sigma (approximate) if present
        // Raw value; normalized to internal nm in enhance_param.
        float r = std::stof(value);
        bias_info.sigma = r;
        bias_info.sigma_squared = r * r;
    } else if (key == "cavity_probe_radius_frag") {
        frag_info.cavity_probe_radius_list = InpParserStructures::parse_float_vector(value);
    } else if (key == "wdens") {
        // Water density output control
        mc_info.wdens = std::stof(value);
    } else if (key == "eps" || key == "epsilon") {
        // Dielectric constant
        frag_info.epsilon = std::stof(value);
    } else if (key == "target_numwaters" || key == "target_num_waters") {
        // Target number of water molecules - store in both places
        frag_info.target_num_waters = std::stoi(value);
    } else if (key == "cavity_mask_frag") {
        frag_info.cavity_mask_list = InpParserStructures::parse_int_vector(value);
    } else if (key == "gcmc_region") {
        // GCMC insertion region (sphere/box specification)
        space_info.gcmc_region = value;
    } else if (key == "exclude_protein_volume") {
        // Exclude protein volume from cavity bias
        space_info.exclude_protein_volume = (value == "yes" || value == "true" || value == "1");
    } else if (key == "use_vdw_radius_for_grid" || key == "use_vdw_radii_for_grid") {
        // Use VDW radii for grid generation
        space_info.use_vdw_radius_for_grid = (value == "yes" || value == "true" || value == "1");
    } else if (key == "exclude_hydrogens_from_grid") {
        // Exclude hydrogens from grid occupancy
        space_info.exclude_hydrogens_from_grid = (value == "yes" || value == "true" || value == "1");
    } else if (key == "use_switching") {
        // Enable switching function
        mc_info.use_switching = (value == "yes" || value == "true" || value == "1");
    } else if (key == "switch_r_on" || key == "switch_ron") {
        // Switching function r_on
        mc_info.switch_r_on = std::stof(value);  // Already in nm
    } else if (key == "switch_r_off" || key == "switch_roff") {
        // Switching function r_off
        mc_info.switch_r_off = std::stof(value);  // Already in nm
    } else if (key == "switch_dist_frag" || key == "switch_dist_fragment") {
        float dist = std::stof(value);
        energy_info.switch_dist_fragment = dist;
        energy_info.switch_dist_fragment_squared = dist * dist;
        energy_info.use_switching = true;
    } else if (key == "switch_dist_prot" || key == "switch_dist_protein") {
        float dist = std::stof(value);
        energy_info.switch_dist_protein = dist;
        energy_info.switch_dist_protein_squared = dist * dist;
        energy_info.use_switching = true;
    } else if (key == "pairlist_freq") {
        // Pairlist update frequency
        energy_info.pairlist_freq = static_cast<unsigned int>(std::stoi(value));
    } else if (key == "use_group_cutoff") {
        // Use group-based cutoff instead of atom-based
        energy_info.use_group_cutoff = (value == "yes" || value == "true" || value == "1");
    } else if (key == "pairlist_cutoff") {
        // Pairlist cutoff distance for fragments
        energy_info.pairlist_cutoff = std::stof(value);  // Already in nm
        energy_info.pairlist_cutoff_squared = energy_info.pairlist_cutoff * energy_info.pairlist_cutoff;
        energy_info.pair_list_cutoff_fragment = energy_info.pairlist_cutoff;
        energy_info.pair_list_cutoff_fragment_squared = energy_info.pairlist_cutoff_squared;
    } else if (key == "pairlist_cutoff_protein") {
        // Pairlist cutoff distance for protein
        float cutoff = std::stof(value);  // Already in nm
        energy_info.pair_list_cutoff_protein = cutoff;
        energy_info.pair_list_cutoff_protein_squared = cutoff * cutoff;
    } else if (key == "attempt_prob_ins") {
        // Per-fragment insertion attempt probabilities
        mc_info.attempt_prob_ins = InpParserStructures::parse_float_vector(value);
    } else if (key == "attempt_prob_del") {
        // Per-fragment deletion attempt probabilities
        mc_info.attempt_prob_del = InpParserStructures::parse_float_vector(value);
    } else if (key == "attempt_prob_trn") {
        // Per-fragment translation attempt probabilities
        mc_info.attempt_prob_trn = InpParserStructures::parse_float_vector(value);
    } else if (key == "attempt_prob_rot") {
        // Per-fragment rotation attempt probabilities
        mc_info.attempt_prob_rot = InpParserStructures::parse_float_vector(value);
    } else if (key == "mc_move_prob") {
        // Legacy format: four weights [insert, delete, translate, rotate]
        // Store in attempt_prob_* temporarily, will be broadcasted in enhance_param
        auto probs = InpParserStructures::parse_float_vector(value);

        if (probs.size() >= 4) {
            // Clear any existing values and store the four probabilities
            // These will be broadcasted to all fragments in enhance_param()
            mc_info.attempt_prob_ins = {probs[0]};
            mc_info.attempt_prob_del = {probs[1]};
            mc_info.attempt_prob_trn = {probs[2]};
            mc_info.attempt_prob_rot = {probs[3]};

            // Set a flag to indicate mc_move_prob was used (store in first element as negative to mark)
            // This is a hack to avoid adding new fields to MCInfo
            // We'll check this in enhance_param and broadcast to all fragments
            std::cout << "[INP] Parsed mc_move_prob: "
                      << probs[0] << " (ins), "
                      << probs[1] << " (del), "
                      << probs[2] << " (trn), "
                      << probs[3] << " (rot)"
                      << " - will be applied to all fragments" << std::endl;
        } else {
            std::cerr << "[WARNING] mc_move_prob requires 4 values, got "
                      << probs.size() << std::endl;
        }
    } else if (key == "const_water_nbar") {
        // Fixed target number of water molecules
        frag_info.use_const_water_nbar = true;
        frag_info.const_water_nbar = static_cast<int>(std::stoi(value));
    } else if (key == "number_water_nbar") {
        // Toggle number-based nbar (use current water count as target)
        frag_info.use_number_water_nbar = (value == "yes" || value == "true" || value == "1");
        if (frag_info.use_number_water_nbar) {
            frag_info.use_const_water_nbar = false;
        }
    } else if (key == "volume_water_nbar") {
        // Volume-based nbar via target concentration (M)
        // Use this to override water concentration for volume projection mode
        // (default mode already uses concentration, but we record explicitly)
        const float concM = std::stof(value);
        if (!frag_info.conc_list.empty()) {
            // Override first (water) entry if present
            frag_info.conc_list[0] = concM;
        } else {
            frag_info.conc_list.push_back(concM);
        }
    } else if (key == "insdel_frac" || key == "insdel_fraction") {
        float frac = std::stof(value);
        mc_info.insertion_deletion_frac = std::max(0.0f, std::min(1.0f, frac));
        mc_info.translation_rotation_frac = 1.0f - mc_info.insertion_deletion_frac;
    } else if (key == "attempt_prob_frag") {
        mc_info.fragment_prob = InpParserStructures::parse_float_vector(value);
    } else if (key == "attempt_prob_water") {
        mc_info.water_prob = InpParserStructures::parse_float_vector(value);
    } else if (key == "attempt_prob_atom") {
        mc_info.atom_prob = InpParserStructures::parse_float_vector(value);
    } else if (key == "test_energy") {
        energy_info.test_energy = (value == "yes" || value == "true" || value == "1");
    } else if (key == "test_sw_filters" || key == "test_SW_filters") {
        energy_info.test_sw_filters = (value == "yes" || value == "true" || value == "1");
    } else if (key == "apply_sw_filters" || key == "apply_SW_filters") {
        energy_info.apply_sw_filters = (value == "yes" || value == "true" || value == "1");
    } else if (key == "sw_reference" || key == "SW_reference") {
        // Legacy inputs provide kcal/mol – convert to kJ/mol for internal use
        energy_info.energy_sw_ref = std::stof(value) * 4.184f;
    } else if (key == "sw_scale" || key == "SW_scale") {
        energy_info.energy_sw_scale = std::stof(value) * 4.184f;
    } else if (key == "rotate_dihedral" || key == "rotate_dih_status") {
        // Legacy switch: enable/disable dihedral rotation
        if (value == "yes" || value == "true" || value == "1") {
            mc_info.rotate_dih_status = 1;
        } else if (value == "no" || value == "false" || value == "0") {
            mc_info.rotate_dih_status = 0;
        } else {
            // Allow numeric value passthrough
            try {
                mc_info.rotate_dih_status = std::stoi(value);
            } catch (const std::exception&) {
                mc_info.rotate_dih_status = 0;
            }
        }
    } else if (key == "remove_init") {
        frag_info.remove_init = InpParserStructures::parse_int_vector(value);
        frag_info.flag_remove_init = 1;
    } else if (key == "remove_excess") {
        frag_info.remove_excess = InpParserStructures::parse_int_vector(value);
        frag_info.flag_remove_excess = 1;
    } else if (key == "initial_fragments_cutoff") {
        const float cutoff = std::stof(value);
        frag_info.init_cutoff = cutoff;
        frag_info.init_cutoff_squared = cutoff * cutoff;
    } else if (key == "excess_fragments_threshold") {
        frag_info.excess_threshold = std::stof(value);
    } else if (key == "gcmc_cutoff") {
        frag_info.gcmc_cutoff = std::stof(value);
        frag_info.gcmc_cutoff_squared = frag_info.gcmc_cutoff * frag_info.gcmc_cutoff;
    } else if (key == "use_gcmc_cutoff") {
        frag_info.use_gcmc_cutoff = (value == "yes" || value == "true" || value == "1");
    } else if (key == "target_volume") {
        space_info.target_volume = std::stof(value);
    } else if (key == "use_const_water_nbar") {
        // gcmc_gpu compatibility: allow either yes/no or an integer value
        try {
            const int n = std::stoi(value);
            if (n > 0) {
                frag_info.use_const_water_nbar = true;
                frag_info.const_water_nbar = n;
            } else {
                frag_info.use_const_water_nbar = false;
            }
        } catch (const std::exception&) {
            frag_info.use_const_water_nbar = (value == "yes" || value == "true" || value == "1");
        }
    } else if (key == "use_number_water_nbar") {
        frag_info.use_number_water_nbar = (value == "yes" || value == "true" || value == "1");
        if (frag_info.use_number_water_nbar) {
            frag_info.use_const_water_nbar = false;
        }
    } else if (key == "fragmqtr") {
        // Legacy: per-fragment MQTR file(s)
        file_info.fragment_mqtr_files.push_back(value);
    } else if (key == "inp_units" || key == "units") {
        // Override unit system ("auto", "nm", "gcmc_gpu"/"angstrom"/"a")
        std::string v = value;
        std::transform(v.begin(), v.end(), v.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        basic_info.inp_units = v;
        basic_info.inp_units_explicit = true;
    }
}

void InpParserGCMC::enhance_param(model::param::Param& param) {
    auto& mc_info = param.get_mc_info();
    auto& file_info = param.get_file_info();
    auto& energy_info = param.get_energy_info();
    auto& bias_info = param.get_bias_info();
    auto& space_info = param.get_space_info();
    auto& fragment_info = param.get_fragment_info();
    auto& basic_info = param.get_basic_info();

    auto toLower = [](std::string s) {
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        return s;
    };

    // --- Unit conversion layer (gcmc_gpu compatibility) ---
    // Internal convention: nm + kJ/mol. Legacy gcmc_gpu INP commonly uses Å + kcal/mol.
    // We make conversion idempotent to tolerate repeated enhance_param calls.
    if (!basic_info.inp_units_converted) {
        enum class UnitMode { NmKj, AngstromKcal };

        UnitMode mode = UnitMode::NmKj;
        std::string u = toLower(basic_info.inp_units);

        // If the INP declares a legacy gcmc_gpu-style version (e.g., "gcmc_2.0"),
        // assume Å + kcal/mol unless the user explicitly overrides via inp_units.
        // Note: The native pygcmc_dev default version is "gcmc_v2.0" (with a 'v').
        const auto looksLikeLegacyGcmcGpuVersion = [&]() -> bool {
            const std::string v = toLower(basic_info.version);
            if (v == "gcmc_gpu") {
                return true;
            }
            if (v.rfind("gcmc_", 0) == 0) {
                // Treat "gcmc_2.0", "gcmc_3.1", ... as legacy (but not "gcmc_v2.0").
                if (v.size() > 5) {
                    const unsigned char c = static_cast<unsigned char>(v[5]);
                    return std::isdigit(c) != 0;
                }
            }
            return false;
        };

        if (!basic_info.inp_units_explicit && looksLikeLegacyGcmcGpuVersion()) {
            u = "gcmc_gpu";
            basic_info.inp_units = u;
        }

        const auto isNm = [&]() {
            return (u == "nm" || u == "nm_kj" || u == "nm/kj" || u == "native");
        };
        const auto isAngstrom = [&]() {
            return (u == "a" || u == "ang" || u == "angstrom" || u == "gcmc_gpu" || u == "a_kcal");
        };

        if (isNm()) {
            mode = UnitMode::NmKj;
        } else if (isAngstrom()) {
            mode = UnitMode::AngstromKcal;
        } else {
            // auto heuristic: legacy gcmc_gpu inputs are typically in Å and have "large" boxes/centers (~30-50).
            // Prefer not converting ambiguous small boxes (<= ~25) unless the user explicitly requests it.
            const float maxBox = std::max({std::abs(space_info.box_size[0]),
                                           std::abs(space_info.box_size[1]),
                                           std::abs(space_info.box_size[2])});
            const float maxGc = std::max({std::abs(space_info.gc_center[0]),
                                          std::abs(space_info.gc_center[1]),
                                          std::abs(space_info.gc_center[2])});
            const float maxSys = std::max({std::abs(space_info.sys_center[0]),
                                           std::abs(space_info.sys_center[1]),
                                           std::abs(space_info.sys_center[2])});
            const float maxCrystal = std::max({std::abs(space_info.crystal_dim[0]),
                                               std::abs(space_info.crystal_dim[1]),
                                               std::abs(space_info.crystal_dim[2])});
            const float maxLen = std::max({maxBox, maxGc, maxSys, maxCrystal});

            if (maxLen > 25.0f) {
                mode = UnitMode::AngstromKcal;
            } else if (maxLen == 0.0f && space_info.cutoff > 5.0f) {
                // No geometry provided but cutoff resembles Å-style defaults.
                mode = UnitMode::AngstromKcal;
            }
        }

        const float LEN = (mode == UnitMode::AngstromKcal) ? 0.1f : 1.0f;      // Å -> nm
        const float VOL = (mode == UnitMode::AngstromKcal) ? 0.001f : 1.0f;    // Å^3 -> nm^3
        const float ENE = (mode == UnitMode::AngstromKcal) ? 4.184f : 1.0f;    // kcal -> kJ

        auto scale3 = [&](std::array<float, 3>& a) {
            a[0] *= LEN;
            a[1] *= LEN;
            a[2] *= LEN;
        };
        auto scaleList = [&](std::vector<float>& v) {
            for (auto& x : v) x *= LEN;
        };

        // SpaceInfo
        space_info.grid_spacing *= LEN;
        scale3(space_info.gc_center);
        scale3(space_info.sys_center);
        scale3(space_info.crystal_dim);
        scale3(space_info.box_size);
        space_info.cutoff *= LEN;
        space_info.target_volume *= VOL;
        model::param::ParamOperations::updateVolume(space_info);

        // MC switching distances
        mc_info.switch_r_on *= LEN;
        mc_info.switch_r_off *= LEN;

        // Energy cutoffs and pairlists
        energy_info.fragment_cutoff *= LEN;
        energy_info.protein_cutoff *= LEN;
        energy_info.pairlist_cutoff *= LEN;
        energy_info.pair_list_cutoff_fragment *= LEN;
        energy_info.pair_list_cutoff_protein *= LEN;
        energy_info.switch_dist_fragment *= LEN;
        energy_info.switch_dist_protein *= LEN;

        // Bias radii
        bias_info.sigma *= LEN;

        // Fragment radii & cavity params
        scaleList(fragment_info.radius_list);
        scaleList(fragment_info.cavity_grid_dx_list);
        scaleList(fragment_info.cavity_probe_radius_list);
        fragment_info.init_cutoff *= LEN;
        fragment_info.gcmc_cutoff *= LEN;

        // Chemical potentials
        for (auto& mu : fragment_info.muex_list) mu *= ENE;

        // Region constraint string (sphere/box/cylinder) numeric values
        if (mode == UnitMode::AngstromKcal && !space_info.gcmc_region.empty()) {
            std::istringstream iss(space_info.gcmc_region);
            std::vector<std::string> toks;
            std::string t;
            while (iss >> t) toks.push_back(t);

            auto toFloat = [&](const std::string& s) -> float {
                return std::stof(s);
            };

            try {
                if (!toks.empty()) {
                    const std::string type = toks[0];
                    std::ostringstream oss;
                    oss << type;
                    if (type == "sphere" && toks.size() == 5) {
                        for (size_t i = 1; i < toks.size(); ++i) {
                            oss << " " << (toFloat(toks[i]) * LEN);
                        }
                        space_info.gcmc_region = oss.str();
                    } else if (type == "box" && toks.size() == 7) {
                        for (size_t i = 1; i < toks.size(); ++i) {
                            oss << " " << (toFloat(toks[i]) * LEN);
                        }
                        space_info.gcmc_region = oss.str();
                    } else if (type == "cylinder" && toks.size() == 7) {
                        // cylinder x y z r h axis
                        for (size_t i = 1; i < 6; ++i) {
                            oss << " " << (toFloat(toks[i]) * LEN);
                        }
                        oss << " " << toks[6];
                        space_info.gcmc_region = oss.str();
                    }
                }
            } catch (const std::exception&) {
                // Best-effort: keep original region string if parsing fails.
            }
        }

        // gcmc_gpu compatibility: if gcmc_region not provided, derive a box region from gc_center + box_size.
        if (space_info.gcmc_region.empty() &&
            (space_info.box_size[0] > 0.0f || space_info.box_size[1] > 0.0f || space_info.box_size[2] > 0.0f)) {
            const bool hasGcCenter = (space_info.gc_center[0] != 0.0f || space_info.gc_center[1] != 0.0f || space_info.gc_center[2] != 0.0f);
            const bool hasSysCenter = (space_info.sys_center[0] != 0.0f || space_info.sys_center[1] != 0.0f || space_info.sys_center[2] != 0.0f);
            const auto& c = hasGcCenter ? space_info.gc_center : space_info.sys_center;
            if (hasGcCenter || hasSysCenter) {
                const float hx = 0.5f * space_info.box_size[0];
                const float hy = 0.5f * space_info.box_size[1];
                const float hz = 0.5f * space_info.box_size[2];
                std::ostringstream oss;
                oss << "box "
                    << (c[0] - hx) << " " << (c[1] - hy) << " " << (c[2] - hz) << " "
                    << (c[0] + hx) << " " << (c[1] + hy) << " " << (c[2] + hz);
                space_info.gcmc_region = oss.str();
            }
        }

        // Update derived squared values after conversion
        model::param::ParamOperations::updateEnergySquaredValues(energy_info);
        model::param::ParamOperations::updateFragmentSquaredValues(fragment_info);
        model::param::ParamOperations::updateBiasSquaredValues(bias_info);

        basic_info.inp_units_converted = true;
    }

    // Clamp insertion/deletion ratio into [0,1] and mirror to translation/rotation ratio
    mc_info.insertion_deletion_frac = std::max(0.0f, std::min(1.0f, mc_info.insertion_deletion_frac));
    mc_info.translation_rotation_frac = std::max(0.0f, 1.0f - mc_info.insertion_deletion_frac);

    // Ensure switching flags are synchronized between MC and energy configs
    if (mc_info.use_switching) {
        energy_info.use_switching = true;
        if (energy_info.switch_dist_fragment == 0.0f) {
            energy_info.switch_dist_fragment = mc_info.switch_r_on;
        }
        if (energy_info.switch_dist_protein == 0.0f) {
            energy_info.switch_dist_protein = mc_info.switch_r_off;
        }
    } else if (energy_info.use_switching) {
        mc_info.use_switching = true;
        if (mc_info.switch_r_on <= 0.0f) {
            mc_info.switch_r_on = energy_info.switch_dist_fragment;
        }
        if (mc_info.switch_r_off <= 0.0f) {
            mc_info.switch_r_off = energy_info.switch_dist_protein > 0.0f ?
                energy_info.switch_dist_protein : energy_info.switch_dist_fragment;
        }
    }

    // Broadcast mc_move_prob to all fragments if it was used
    // If attempt_prob_* vectors have size 1, it means mc_move_prob was set
    // We need to broadcast to all fragments
    if (mc_info.attempt_prob_ins.size() == 1 &&
        mc_info.attempt_prob_del.size() == 1 &&
        mc_info.attempt_prob_trn.size() == 1 &&
        mc_info.attempt_prob_rot.size() == 1) {

        size_t fragment_count = file_info.fragment_names.size();
        if (fragment_count > 1) {
            // Broadcast the single values to all fragments
            float ins_val = mc_info.attempt_prob_ins[0];
            float del_val = mc_info.attempt_prob_del[0];
            float trn_val = mc_info.attempt_prob_trn[0];
            float rot_val = mc_info.attempt_prob_rot[0];

            mc_info.attempt_prob_ins.assign(fragment_count, ins_val);
            mc_info.attempt_prob_del.assign(fragment_count, del_val);
            mc_info.attempt_prob_trn.assign(fragment_count, trn_val);
            mc_info.attempt_prob_rot.assign(fragment_count, rot_val);

            std::cout << "[INP] Broadcasted mc_move_prob to " << fragment_count << " fragments" << std::endl;
        }
    }

    // If no per-fragment move probabilities were provided, derive them from insdel_frac
    if (mc_info.attempt_prob_ins.empty() &&
        mc_info.attempt_prob_del.empty() &&
        mc_info.attempt_prob_trn.empty() &&
        mc_info.attempt_prob_rot.empty()) {

        size_t fragment_count = file_info.fragment_names.empty() ? 1 : file_info.fragment_names.size();
        const float ins_prob = mc_info.insertion_deletion_frac * 0.5f;
        const float del_prob = mc_info.insertion_deletion_frac * 0.5f;
        const float trn_prob = mc_info.translation_rotation_frac * 0.5f;
        const float rot_prob = mc_info.translation_rotation_frac * 0.5f;

        mc_info.attempt_prob_ins.assign(fragment_count, ins_prob);
        mc_info.attempt_prob_del.assign(fragment_count, del_prob);
        mc_info.attempt_prob_trn.assign(fragment_count, trn_prob);
        mc_info.attempt_prob_rot.assign(fragment_count, rot_prob);
    }

    const size_t fragment_count = file_info.fragment_names.empty() ? 1 : file_info.fragment_names.size();
    auto ensureFloatList = [&](std::vector<float>& vec, float fallback) {
        if (vec.empty()) {
            vec.assign(fragment_count, fallback);
        } else if (vec.size() < fragment_count) {
            vec.resize(fragment_count, fallback);
        }
    };
    auto ensureIntList = [&](std::vector<int>& vec, int fallback) {
        if (vec.empty()) {
            vec.assign(fragment_count, fallback);
        } else if (vec.size() < fragment_count) {
            vec.resize(fragment_count, fallback);
        }
    };

    const float defaultGrid = space_info.grid_spacing > 0.0f ? space_info.grid_spacing : 0.2f;
    const float defaultProbe = bias_info.sigma > 0.0f ? bias_info.sigma : 0.14f;
    ensureFloatList(fragment_info.cavity_grid_dx_list, defaultGrid);
    ensureFloatList(fragment_info.cavity_probe_radius_list, defaultProbe);
    ensureIntList(fragment_info.cavity_mask_list, -1);

    // Ensure mc_time_cumulative if mc_time_list provided
    if (!mc_info.mc_time_list.empty()) {
        mc_info.mc_time_cumulative.clear();
        float s = 0.0f;
        for (float w : mc_info.mc_time_list) {
            s += w;
            mc_info.mc_time_cumulative.push_back(s);
        }
        if (s > 0.0f) {
            for (auto& v : mc_info.mc_time_cumulative) v /= s;
        }
    }

    // Update derived squared quantities (safe to repeat; units conversion is idempotent)
    model::param::ParamOperations::updateEnergySquaredValues(energy_info);
    model::param::ParamOperations::updateFragmentSquaredValues(fragment_info);
    model::param::ParamOperations::updateBiasSquaredValues(bias_info);
}

} // namespace parameters
} // namespace io
} // namespace pygcmc

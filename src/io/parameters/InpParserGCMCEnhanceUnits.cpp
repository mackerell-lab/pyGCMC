#include "InpParserGCMC.hpp"
#include "../../model/param/ParamOperations.hpp"
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace pygcmc {
namespace io {
namespace parameters {

bool InpParserGCMC::enhance_param_units(model::param::Param& param) {
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

    enum class UnitMode { NmKj, AngstromKcal };

    UnitMode mode = UnitMode::NmKj;
    std::string u = toLower(basic_info.inp_units);

    // Alias: "openmm" is treated as native nm + kJ/mol (internal unit system).
    if (u == "openmm") {
        u = "nm";
        basic_info.inp_units = u;
    }

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
        return (u == "a" || u == "ang" || u == "angstrom" ||
                u == "gcmc_gpu" || u == "charmm" ||
                u == "a_kcal" || u == "a/kcal" || u == "akcal");
    };

    if (isNm()) {
        mode = UnitMode::NmKj;
    } else if (isAngstrom()) {
        mode = UnitMode::AngstromKcal;
    } else if (u == "auto") {
        // Default to gcmc_gpu-style units (Å + kcal/mol). Native nm/kJ decks should explicitly set inp_units:nm.
        mode = UnitMode::AngstromKcal;
    } else {
        // Unknown unit tag: fall back to internal nm/kJ without conversion.
    }

    // Materialize the inferred unit mode for downstream consumers (e.g., --dump-params, periodic-box policy).
    // Preserve explicit concrete units; resolve "auto" to a concrete mode.
    if (u == "auto") {
        basic_info.inp_units = (mode == UnitMode::AngstromKcal) ? "gcmc_gpu" : "nm";
    }

    const float LEN = (mode == UnitMode::AngstromKcal) ? 0.1f : 1.0f;      // Å -> nm
    const float VOL = (mode == UnitMode::AngstromKcal) ? 0.001f : 1.0f;    // Å^3 -> nm^3
    const float ENE = (mode == UnitMode::AngstromKcal) ? 4.184f : 1.0f;    // kcal -> kJ

    auto pushWarning = [&](const std::string& code, const std::string& message) {
        auto& ws = basic_info.inp_warnings;
        const auto it = std::find_if(ws.begin(), ws.end(), [&](const auto& w) { return w.code == code; });
        if (it == ws.end()) {
            ws.push_back({code, message});
        }
    };

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
    // MC move step sizes
    mc_info.max_translation_dist *= LEN;

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
    // NOTE: Most legacy INP decks (gcmc_gpu/opencl) are written in Å (length) + kcal/mol (energy).
    // For inp_units:nm decks, we expect nm + kJ/mol, but keep a small heuristic here so that
    // legacy defaults (e.g., BiasInfo.sigma=2.4 which historically meant 2.4 Å) do not become
    // an unphysical 2.4 nm when the user enables cavity bias without an explicit probe radius.
    if (mode == UnitMode::NmKj && bias_info.sigma > 1.0f) {
        const float original = bias_info.sigma;
        const float converted = bias_info.sigma * 0.1f;
        std::ostringstream oss;
        oss << "probe_radius/sigma=" << original
            << " interpreted as nm (inp_units=" << basic_info.inp_units
            << "); this is unusually large, assuming legacy Å value and converting to " << converted << " nm.";
        pushWarning("UNIT_HEURISTIC_PROBE_RADIUS_ASSUMED_ANGSTROM", oss.str());
        bias_info.sigma *= 0.1f;
    }

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

    (void)file_info;
    return mode == UnitMode::NmKj;
}

} // namespace parameters
} // namespace io
} // namespace pygcmc


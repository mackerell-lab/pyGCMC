#include "InpParserGCMC.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace pygcmc {
namespace io {
namespace parameters {

void InpParserGCMC::enhance_param_warnings(model::param::Param& param, bool nm_mode) {
    auto& energy_info = param.get_energy_info();
    auto& space_info = param.get_space_info();
    auto& frag_info = param.get_fragment_info();
    auto& bias_info = param.get_bias_info();
    auto& file_info = param.get_file_info();
    auto& basic_info = param.get_basic_info();

    auto pushWarning = [&](const std::string& code, const std::string& message) {
        auto& ws = basic_info.inp_warnings;
        const auto it = std::find_if(ws.begin(), ws.end(), [&](const auto& w) { return w.code == code; });
        if (it == ws.end()) {
            ws.push_back({code, message});
        }
    };

    // --- Unit/geometry heuristics (diagnostic warnings, not fatal by default) ---
    // Goal: prevent "runs but silently wrong" unit mistakes when users opt into nm/openmm decks.
    if (nm_mode && basic_info.inp_units_explicit) {
        // In nm mode, common accidental legacy values are cutoff≈12 (Å) and box_size≈30 (Å).
        // These become 12 nm / 30 nm which are almost always unintended.
        if (space_info.cutoff > 4.0f) {
            std::ostringstream oss;
            oss << "cutoff=" << space_info.cutoff
                << " interpreted as nm (inp_units=" << basic_info.inp_units
                << "); this is unusually large and may indicate Å values were provided (e.g., 12Å -> 1.2nm).";
            pushWarning("UNIT_SUSPECT_CUTOFF_TOO_LARGE_FOR_NM", oss.str());
        }
        const float maxBox = std::max(space_info.box_size[0], std::max(space_info.box_size[1], space_info.box_size[2]));
        if (maxBox > 20.0f && (space_info.box_size[0] > 0.0f || space_info.box_size[1] > 0.0f || space_info.box_size[2] > 0.0f)) {
            std::ostringstream oss;
            oss << "box_size=" << space_info.box_size[0] << " " << space_info.box_size[1] << " " << space_info.box_size[2]
                << " interpreted as nm (inp_units=" << basic_info.inp_units
                << "); this is unusually large and may indicate Å values were provided (e.g., 30Å -> 3nm).";
            pushWarning("UNIT_SUSPECT_BOX_TOO_LARGE_FOR_NM", oss.str());
        }

        // Cavity grid spacing: typical ~0.05-0.5 nm. Values >= 1.0 often mean legacy 1.0 Å.
        if (space_info.grid_spacing > 0.8f) {
            std::ostringstream oss;
            oss << "grid_dx/grid_spacing=" << space_info.grid_spacing
                << " interpreted as nm (inp_units=" << basic_info.inp_units
                << "); this is unusually coarse and may indicate Å values were provided (e.g., 1Å -> 0.1nm).";
            pushWarning("UNIT_SUSPECT_GRID_SPACING_TOO_LARGE_FOR_NM", oss.str());
        }

        // Region constraints: prevent accidental 10x scale errors in nm decks.
        if (!space_info.gcmc_region.empty()) {
            std::istringstream iss(space_info.gcmc_region);
            std::vector<std::string> toks;
            std::string t;
            while (iss >> t) toks.push_back(t);

            auto parseFloat = [](const std::string& s, float& out) -> bool {
                try {
                    out = std::stof(s);
                    return true;
                } catch (const std::exception&) {
                    return false;
                }
            };

            bool suspect = false;
            if (!toks.empty()) {
                std::string type = toks[0];
                std::transform(type.begin(), type.end(), type.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

                if (type == "sphere" && toks.size() == 5) {
                    float r = 0.0f;
                    if (parseFloat(toks[4], r) && r > 10.0f) {
                        suspect = true;
                    }
                } else if (type == "box" && toks.size() == 7) {
                    float x1 = 0.0f, y1 = 0.0f, z1 = 0.0f, x2 = 0.0f, y2 = 0.0f, z2 = 0.0f;
                    if (parseFloat(toks[1], x1) && parseFloat(toks[2], y1) && parseFloat(toks[3], z1) &&
                        parseFloat(toks[4], x2) && parseFloat(toks[5], y2) && parseFloat(toks[6], z2)) {
                        const float dx = std::abs(x2 - x1);
                        const float dy = std::abs(y2 - y1);
                        const float dz = std::abs(z2 - z1);
                        if (dx > 20.0f || dy > 20.0f || dz > 20.0f) {
                            suspect = true;
                        }
                    }
                } else if (type == "cylinder" && toks.size() == 7) {
                    float r = 0.0f, h = 0.0f;
                    // cylinder x y z r h axis
                    if (parseFloat(toks[4], r) && parseFloat(toks[5], h)) {
                        if (r > 10.0f || h > 20.0f) {
                            suspect = true;
                        }
                    }
                }
            }

            if (suspect) {
                std::ostringstream oss;
                oss << "gcmc_region=\"" << space_info.gcmc_region
                    << "\" interpreted as nm (inp_units=" << basic_info.inp_units
                    << "); this looks unusually large and may indicate Å values were provided (10x scale error).";
                pushWarning("UNIT_SUSPECT_GCMC_REGION_TOO_LARGE_FOR_NM", oss.str());
            }
        }
    }

    // --- CBMC configuration (unit-agnostic) ---
    // Prevent "use_conf_bias enabled but trial-count source missing" from silently defaulting.
    // gcmc_gpu typically provides per-fragment `fragconf`, while gcmc_opencl uses a global
    // `num_conf_bias_trial`. We accept both, but warn when neither is present.
    if (bias_info.use_conf_bias) {
        const auto& handled = basic_info.inp_keys_handled;
        auto hasHandledKey = [&](const std::string& k) {
            return std::find(handled.begin(), handled.end(), k) != handled.end();
        };

        const bool hasFragconfKey = hasHandledKey("fragconf") || hasHandledKey("fragconfs");
        const bool hasGlobalTrialsKey = hasHandledKey("num_conf_bias_trial") || hasHandledKey("confbias_trials");
        const bool hasAnyTrialsKey = hasFragconfKey || hasGlobalTrialsKey;

        if (!hasAnyTrialsKey) {
            std::ostringstream oss;
            oss << "use_conf_bias=yes but no CBMC trial-count key provided "
                << "(fragconf/num_conf_bias_trial/confbias_trials); defaulting to "
                << "num_conf_bias_trials=" << bias_info.num_conf_bias_trials << " for all fragments.";
            pushWarning("CBMC_TRIAL_COUNT_DEFAULTED", oss.str());
        }

        // If per-fragment fragconf is provided, its length should match fragname count.
        // We still accept shorter lists by falling back to the global value, but warn.
        if (hasFragconfKey && !file_info.fragment_names.empty() &&
            !frag_info.fragconf_list.empty() &&
            frag_info.fragconf_list.size() != file_info.fragment_names.size()) {
            std::ostringstream oss;
            oss << "fragconf list length (" << frag_info.fragconf_list.size()
                << ") does not match fragname count (" << file_info.fragment_names.size()
                << "); missing entries will fall back to num_conf_bias_trials="
                << bias_info.num_conf_bias_trials << ".";
            pushWarning("CBMC_TRIAL_COUNT_LIST_LENGTH_MISMATCH", oss.str());
        }
    }

    // --- Cross-field consistency heuristics (unit-agnostic, still non-fatal by default) ---
    // target_volume is rarely intended to differ from the INP box_size by orders of magnitude.
    const float boxVol = space_info.box_size[0] * space_info.box_size[1] * space_info.box_size[2];
    if (space_info.target_volume > 0.0f && boxVol > 0.0f) {
        const float ratio = space_info.target_volume / boxVol;
        if (ratio > 10.0f || ratio < 0.1f) {
            std::ostringstream oss;
            oss << "target_volume=" << space_info.target_volume << " nm^3 vs box_size volume=" << boxVol
                << " nm^3 (ratio=" << ratio
                << "); this looks inconsistent and may indicate a unit mismatch (e.g., Å^3 provided in nm mode, or nm^3 provided in Å mode).";
            pushWarning("UNIT_SUSPECT_TARGET_VOLUME_BOX_MISMATCH", oss.str());
        }
    }

    // pairlist_cutoff should typically be >= cutoff (often cutoff + skin). Flag suspicious relations.
    if (energy_info.pairlist_cutoff > 0.0f && space_info.cutoff > 0.0f) {
        const float diff = energy_info.pairlist_cutoff - space_info.cutoff;
        if (diff < 0.0f) {
            std::ostringstream oss;
            oss << "pairlist_cutoff=" << energy_info.pairlist_cutoff << " nm is smaller than cutoff=" << space_info.cutoff
                << " nm; this is likely unintended and could cause missing interactions if pairlists are used.";
            pushWarning("PAIRLIST_CUTOFF_SMALLER_THAN_CUTOFF", oss.str());
        } else if (diff > 2.0f) {
            std::ostringstream oss;
            oss << "pairlist_cutoff=" << energy_info.pairlist_cutoff << " nm is much larger than cutoff=" << space_info.cutoff
                << " nm (diff=" << diff
                << "); this may waste compute and could indicate a unit mismatch.";
            pushWarning("PAIRLIST_CUTOFF_EXCESSIVELY_LARGE", oss.str());
        }
    }
}

} // namespace parameters
} // namespace io
} // namespace pygcmc

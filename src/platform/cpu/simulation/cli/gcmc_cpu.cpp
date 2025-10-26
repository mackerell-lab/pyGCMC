// Simple command-line GCMC driver (CPU)
// Usage example:
//   gcmc_cpu --inp run.inp --prefix out --seed 123 --verbose
//            --print-freq 1000 --traj-freq 10000 --checkpoint-freq 50000
//            --stats-interval 1000 --adaptive --store-probabilities

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdlib>

#include "platform/cpu/simulation/impl/GCMCSimulation.hpp"

using pygcmc::platform::cpu::simulation::GCMCSimulation;

namespace {

struct Options {
    std::string inp;
    std::string prefix = "gcmc";
    std::string resume_checkpoint;  // Checkpoint file to resume from
    std::string dump_accept;  // Acceptance log output file
    int print_freq = -1; // -1 means: use INP nprint
    int traj_freq = 10000;
    int checkpoint_freq = 0;  // Disabled by default
    int stats_interval = 1000;
    int seed = -1;
    bool verbose = false;
    bool enable_stats = true;
    bool adaptive = false;
    bool store_prob = false;
    int max_molecules_per_type = 10000;  // Default value, -1 to disable
};

void print_usage(const char* prog) {
    std::cout << "GCMC CPU Driver (pygcmc_dev)\n";
    std::cout << "Usage: " << prog << " --inp <file.inp> [options]\n\n";
    std::cout << "Required:\n";
    std::cout << "  --inp <path>               Input INP file (gcmc_gpu-compatible keys)\n\n";
    std::cout << "Options:\n";
    std::cout << "  --prefix <str>             Output prefix (default: gcmc)\n";
    std::cout << "  --resume <path>            Resume from checkpoint file\n";
    std::cout << "  --dump-accept <path>       Dump acceptance log (JSONL format)\n";
    std::cout << "  --seed <int>               RNG seed (-1 = auto)\n";
    std::cout << "  --verbose                  Verbose logging\n";
    std::cout << "  --no-stats                 Disable statistics collection\n";
    std::cout << "  --stats-interval <int>     Statistics sampling interval (default: 1000)\n";
    std::cout << "  --print-freq <int>         Print frequency (default: use INP nprint)\n";
    std::cout << "  --traj-freq <int>          Trajectory output frequency (default: 10000)\n";
    std::cout << "  --checkpoint-freq <int>    Checkpoint frequency (default: 0 = disabled)\n";
    std::cout << "  --max-molecules-per-type <int> Max molecules per fragment type (default: 10000, -1 = disabled)\n";
    std::cout << "  --adaptive                 Enable adaptive sampling of move probs\n";
    std::cout << "  --store-probabilities      Store acceptance probabilities\n";
}

bool parse_int(const char* arg, int& out) {
    char* end = nullptr;
    long v = std::strtol(arg, &end, 10);
    if (end == arg || *end != '\0') return false;
    out = static_cast<int>(v);
    return true;
}

bool parse_args(int argc, char** argv, Options& opt) {
    if (argc < 3) return false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for option '" << name << "'\n";
                return nullptr;
            }
            return argv[++i];
        };
        if (a == "--inp") {
            const char* v = need("--inp");
            if (!v) return false;
            opt.inp = v;
        } else if (a == "--prefix") {
            const char* v = need("--prefix");
            if (!v) return false;
            opt.prefix = v;
        } else if (a == "--resume") {
            const char* v = need("--resume");
            if (!v) return false;
            opt.resume_checkpoint = v;
        } else if (a == "--dump-accept") {
            const char* v = need("--dump-accept");
            if (!v) return false;
            opt.dump_accept = v;
        } else if (a == "--print-freq") {
            const char* v = need("--print-freq");
            if (!v || !parse_int(v, opt.print_freq)) return false;
        } else if (a == "--traj-freq") {
            const char* v = need("--traj-freq");
            if (!v || !parse_int(v, opt.traj_freq)) return false;
        } else if (a == "--checkpoint-freq") {
            const char* v = need("--checkpoint-freq");
            if (!v || !parse_int(v, opt.checkpoint_freq)) return false;
        } else if (a == "--stats-interval") {
            const char* v = need("--stats-interval");
            if (!v || !parse_int(v, opt.stats_interval)) return false;
        } else if (a == "--seed") {
            const char* v = need("--seed");
            if (!v || !parse_int(v, opt.seed)) return false;
        } else if (a == "--verbose") {
            opt.verbose = true;
        } else if (a == "--no-stats") {
            opt.enable_stats = false;
        } else if (a == "--adaptive") {
            opt.adaptive = true;
        } else if (a == "--store-probabilities") {
            opt.store_prob = true;
        } else if (a == "--max-molecules-per-type") {
            const char* v = need("--max-molecules-per-type");
            if (!v || !parse_int(v, opt.max_molecules_per_type)) return false;
        } else if (a == "-h" || a == "--help") {
            return false;
        } else {
            std::cerr << "Unknown option: " << a << "\n";
            return false;
        }
    }
    return !opt.inp.empty();
}

} // namespace

int main(int argc, char** argv) {
    Options opt;
    if (!parse_args(argc, argv, opt)) {
        print_usage(argv[0]);
        return 1;
    }

    // Configure simulation
    GCMCSimulation::Config cfg;
    cfg.inputFile = opt.inp;
    cfg.outputPrefix = opt.prefix;
    cfg.printFrequency = opt.print_freq;
    cfg.trajectoryFrequency = opt.traj_freq;
    cfg.checkpointFrequency = opt.checkpoint_freq;
    cfg.verbose = opt.verbose;
    cfg.randomSeed = opt.seed;
    cfg.enableStatistics = opt.enable_stats;
    cfg.statisticsInterval = opt.stats_interval;
    cfg.storeProbabilities = opt.store_prob;
    cfg.enableAdaptiveSampling = opt.adaptive;
    cfg.maxMoleculesPerType = opt.max_molecules_per_type;

    // Create and run
    GCMCSimulation sim(cfg);
    if (!sim.initialize()) {
        // Mirror error to stdout to satisfy tests that check stdout
        std::cout << "ERROR: Failed to initialize GCMC simulation." << std::endl;
        std::cerr << "Failed to initialize GCMC simulation." << std::endl;
        return 2;
    }

    // Enable diagnostics if acceptance log requested
    if (!opt.dump_accept.empty()) {
        sim.enableDiagnostics(65536);  // Buffer size for acceptance log
        std::cout << "Diagnostics enabled for acceptance log export" << std::endl;
    }

    // Load checkpoint if specified
    if (!opt.resume_checkpoint.empty()) {
        std::cout << "Loading checkpoint from: " << opt.resume_checkpoint << std::endl;
        if (!sim.loadCheckpoint(opt.resume_checkpoint)) {
            std::cout << "ERROR: Failed to load checkpoint file: " << opt.resume_checkpoint << std::endl;
            std::cerr << "Failed to load checkpoint file: " << opt.resume_checkpoint << std::endl;
            return 4;
        }
        std::cout << "Checkpoint loaded successfully. Resuming simulation..." << std::endl;
    }

    if (!sim.run()) {
        std::cout << "ERROR: GCMC simulation failed during execution." << std::endl;
        std::cerr << "GCMC simulation failed during execution." << std::endl;
        return 3;
    }
    sim.finalize();

    // Dump acceptance log if requested
    if (!opt.dump_accept.empty()) {
        std::cout << "Dumping acceptance log to: " << opt.dump_accept << std::endl;
        sim.dumpAcceptanceLog(opt.dump_accept);
        std::cout << "Acceptance log exported" << std::endl;
    }

    return 0;
}

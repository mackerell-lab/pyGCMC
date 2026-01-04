#include "GCMCEngine.hpp"
#include "GCMCAcceptance.hpp"
#include "GCMCConfig.hpp"
#include "../reservoir/MultiTypeReservoir.hpp"
#include "../../energy/EnergyModule.hpp"
#include "../../energy/common/EnergyDirectCore.hpp"
#include "../../energy/pme/PMEGlobal.hpp"
#include "../../energy/pme/PMEGridCharge.hpp"
#include "../../energy/pme/PMESetup.hpp"
#include "../../energy/pme/PMESpline.hpp"
#include "../../energy/pme/PMESystemCore.hpp"
#include "../../energy/pgp/PGPComposite.hpp"
#include "../../energy/pgp/PGPGlobal.hpp"
#include "../../energy/pgp/PGPInterpolation.hpp"
#include "../../energy/pgp/PGPSelf.hpp"
#include "../../energy/pgp/PGPPrecompute.hpp"
#include "../../energy/common/EnergyConstants.hpp"
#include "../../energy/drude/DrudeMain.hpp"
#include <cmath>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <unordered_set>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

// Type aliases for clarity
using MCState = model::montecarlo::MCState;
using MCResidue = model::montecarlo::MCResidue;
using MCAtom = model::montecarlo::MCAtom;

namespace {

struct MovementResiduesGuard {
    MCState& state;
    std::vector<model::MCMovementResidueInfo> saved;

    explicit MovementResiduesGuard(MCState& stateIn)
        : state(stateIn), saved(stateIn.movementResidues) {
    }

    void setSingleResidue(int residueIdx) {
        state.movementResidues.clear();
        model::MCMovementResidueInfo info;
        info.startIndex = residueIdx;
        info.activeCount = 1;
        state.movementResidues.push_back(info);
    }

    ~MovementResiduesGuard() {
        state.movementResidues = saved;
    }
};

struct ResidueActiveGuard {
    MCState& state;
    int residueIdx = -1;
    bool savedActive = false;

    ResidueActiveGuard(MCState& stateIn, int residueIdxIn)
        : state(stateIn), residueIdx(residueIdxIn) {
        if (residueIdx >= 0 && residueIdx < static_cast<int>(state.residues.size())) {
            savedActive = state.residues[residueIdx].active;
        }
    }

    void setInactive() {
        if (residueIdx >= 0 && residueIdx < static_cast<int>(state.residues.size())) {
            state.residues[residueIdx].active = false;
        }
    }

    void restore() {
        if (residueIdx >= 0 && residueIdx < static_cast<int>(state.residues.size())) {
            state.residues[residueIdx].active = savedActive;
        }
    }

    ~ResidueActiveGuard() {
        restore();
    }
};

bool hasAnyActiveFixedResidue(const MCState& state) {
    for (int i = 0; i < state.activeResidueCount; ++i) {
        const auto& res = state.residues[i];
        if (res.active && res.fixed) {
            return true;
        }
    }
    return false;
}

} // namespace

void GCMCEngine::setDrudeHostTopology(
    std::vector<::pygcmc::platform::cpu::DrudeParticle> particles,
    std::vector<::pygcmc::platform::cpu::ScreenedPair> screenedPairs) {
    drudeHostParticles_ = std::move(particles);
    drudeHostScreenedPairs_ = std::move(screenedPairs);
    drudeTopologyDirty_ = true;
    ::pygcmc::platform::cpu::DrudeComplete::clearHistory();
}

void GCMCEngine::setDrudeForceFieldModel(const ::pygcmc::model::forcefield::ForceField* forceField) {
    drudeForceFieldModel_ = forceField;
    drudeTopologyDirty_ = true;
}

// Constructor
GCMCEngine::GCMCEngine()
    : state_(nullptr),
      reservoir_(nullptr),
      cavityManager_(nullptr),
      configBias_(nullptr),
      acceptanceCalculator_(nullptr),
      energyCallback_(std::make_unique<GCMCEnergyCallback>()),
      temperature_(300.0),
      cutoff_(12.0),
      energyMethod_(EnergyMethod::DIRECT),
      rng_(0u),  // Deterministic initial value, will be set via setSeed()
      uniform_(0.0, 1.0),
      normal_(0.0, 1.0),
      totalMoves_(0),
      acceptedMoves_(0),
      maxTranslationStep_(2.0),
      maxRotationAngleRad_(0.5),
      useCavityBias_(true) {
    energyCache_.valid = false;
    // Configure default energy callback
    energyCallback_->setEnergyMethod(energyMethod_);
    energyCallback_->setParameters(true, true);  // Use cutoff and PBC by default
}

// Destructor
GCMCEngine::~GCMCEngine() {
}

// Initialize
void GCMCEngine::initialize(MCState* state, FragmentReservoir* reservoir) {
    state_ = state;
    reservoir_ = reservoir;
    energyCache_.invalidate();
    if (state_ && reservoir_) {
        reservoir_->reserveInstanceIds(std::max(0, state_->activeResidueCount));
    }
}

void GCMCEngine::setEnergyBackend(GCMCEnergyBackend backend) {
    energyBackend_ = backend;

    // Invalidate cached grids when switching backends.
    pgpInitialized_ = false;
    pgpHostGridReady_ = false;
    pgpFullGridReady_ = false;
    pgpFullGridExcludedResidue_ = -999;

    // Keep the legacy EnergyMethod in a sensible state for fallback paths.
    switch (backend) {
        case GCMCEnergyBackend::DirectCutoff:
            setEnergyMethod(EnergyMethod::DIRECT);
            break;
        case GCMCEnergyBackend::Ewald:
            setEnergyMethod(EnergyMethod::EWALD);
            break;
        case GCMCEnergyBackend::Pme:
            setEnergyMethod(EnergyMethod::PME);
            break;
        case GCMCEnergyBackend::PgpHost:
        case GCMCEnergyBackend::PgpFull:
        case GCMCEnergyBackend::PgpFullPme:
            // PGP is handled inside GCMCEngine; keep fallback set to DIRECT.
            setEnergyMethod(EnergyMethod::DIRECT);
            break;
    }
}

// Set seed - unified for all RNG components
void GCMCEngine::setSeed(unsigned int seed) {
    // Set engine's RNG seed
    rng_.seed(seed);
    lastSeed_ = seed;  // Store for auto-seeding acceptance
    
    // Also set acceptance calculator's RNG seed if present
    if (acceptanceCalculator_) {
        acceptanceCalculator_->setSeed(seed + 1);  // Use different but deterministic seed
    }
    
    // Set reservoir's RNG seed if it has one
    if (reservoir_) {
        // Note: Add setSeed to FragmentReservoir if it needs random operations
        // reservoir_->setSeed(seed + 2);
    }
    
    // Set cavity manager's RNG seed if it has one
    if (cavityManager_) {
        // Note: Add setSeed to CavityManager if it needs random operations
        // cavityManager_->setSeed(seed + 3);
    }
}

void GCMCEngine::ensurePgpInitialized() {
    if (pgpInitialized_) {
        return;
    }
    if (!state_) {
        throw std::runtime_error("PGP backend requires MCState to be initialized");
    }

    const double cutoff = static_cast<double>(state_->info.cutoff);
    if (!(cutoff > 0.0)) {
        throw std::runtime_error("PGP backend requires a positive cutoff (nm)");
    }

    const double box[3] = {
        static_cast<double>(state_->info.box[0]),
        static_cast<double>(state_->info.box[1]),
        static_cast<double>(state_->info.box[2]),
    };
    if (!(box[0] > 0.0 && box[1] > 0.0 && box[2] > 0.0)) {
        throw std::runtime_error("PGP backend requires a valid periodic box (nm)");
    }

    // Auto-tune PME parameters and reuse them for PGP setup.
    // This keeps PGP consistent with the existing PME error tolerance model.
    autoAdjustPMEParameters(pgpTolerance_, cutoff, box);
    // PGP's setup path relies on PME-global box state for charge spreading and for copying
    // PME parameters into the PGP parameter block. Without this, PGP can silently run with
    // the default unit box and produce wildly mis-scaled electrostatics.
    setPMEBox(box);
    const auto& pme = getPMEParams();
    int meshSize[3] = {pme.meshSize[0], pme.meshSize[1], pme.meshSize[2]};
    // Guard against pathological auto-tuning (e.g., selecting a 4×4×4 mesh for a multi-nm box),
    // which can severely under-resolve the reciprocal potential and silently mis-scale PGP energies.
    // Use a conservative minimum grid spacing for PGP potential interpolation.
    auto nextPow2 = [](int n) {
        int p = 1;
        while (p < n) p <<= 1;
        return p;
    };
    // Heuristic: require a finer reciprocal mesh when alpha is large to keep
    // PGP(interpolated potential) consistent with PME(system energy).
    // For typical alpha~5 1/nm, 0.5/alpha -> 0.1 nm grid spacing (64^3 for a 4 nm box).
    constexpr double kMinGridDxLowerNm = 0.10;
    constexpr double kMinGridDxUpperNm = 0.25;
    double minGridDxNm = kMinGridDxUpperNm;
    if (pme.alpha > 0.0) {
        minGridDxNm = std::min(kMinGridDxUpperNm, std::max(kMinGridDxLowerNm, 0.5 / pme.alpha));
    }
    for (int d = 0; d < 3; ++d) {
        const int minSize = nextPow2(static_cast<int>(std::ceil(box[d] / minGridDxNm)));
        if (meshSize[d] < minSize) {
            meshSize[d] = minSize;
        }
    }

    // Use the same mesh for potential grid size in the first implementation.
    PGPComposite::initialize(
        cutoff,
        box,
        pme.alpha,
        meshSize,
        cutoff,
        meshSize,
        pgpSplineOrder_,
        pgpTolerance_
    );

    pgpInitialized_ = true;
}

void GCMCEngine::buildPgpReciprocalKernel() {
    const auto& pme = getPMEParams();
    const int nx = pme.meshSize[0];
    const int ny = pme.meshSize[1];
    const int nz = pme.meshSize[2];
    if (!(nx > 0 && ny > 0 && nz > 0)) {
        pgpRecipKernel_.clear();
        pgpRecipKernelMesh_[0] = pgpRecipKernelMesh_[1] = pgpRecipKernelMesh_[2] = 0;
        pgpRecipKernelSplineOrder_ = 0;
        return;
    }

    const int totalGridSize = nx * ny * nz;
    if (!(totalGridSize > 0)) {
        pgpRecipKernel_.clear();
        pgpRecipKernelMesh_[0] = pgpRecipKernelMesh_[1] = pgpRecipKernelMesh_[2] = 0;
        pgpRecipKernelSplineOrder_ = 0;
        return;
    }

    auto& pmeMutable = getPMEParams();
    std::vector<std::complex<double>> pmeGridBackup = std::move(pmeMutable.pmeGrid);

    pmeMutable.pmeGrid.assign(static_cast<size_t>(totalGridSize), std::complex<double>(0.0, 0.0));
    pmeMutable.pmeGrid[0] = std::complex<double>(1.0, 0.0);

    performFFTForward();

    const double boxNm[3] = {pmeMutable.box[0], pmeMutable.box[1], pmeMutable.box[2]};
    double unusedRecipEnergy = 0.0;
    computeEnergyFromGrid(unusedRecipEnergy, boxNm);
    if (!pmeMutable.pmeGrid.empty()) {
        pmeMutable.pmeGrid[0] = std::complex<double>(0.0, 0.0);
    }

    performFFTBackward();

    // Undo inverse FFT normalization (1/(nx*ny*nz)) to match the physical scaling used by PGP
    // potential grids and by the quadratic-form mesh-self evaluation.
    const double fftScale = static_cast<double>(totalGridSize);
    for (int i = 0; i < totalGridSize; ++i) {
        pmeMutable.pmeGrid[static_cast<size_t>(i)] *= fftScale;
    }

    pgpRecipKernel_.assign(static_cast<size_t>(totalGridSize), 0.0);
    for (int i = 0; i < totalGridSize; ++i) {
        pgpRecipKernel_[static_cast<size_t>(i)] = pmeMutable.pmeGrid[static_cast<size_t>(i)].real();
    }

    pgpRecipKernelMesh_[0] = nx;
    pgpRecipKernelMesh_[1] = ny;
    pgpRecipKernelMesh_[2] = nz;
    pgpRecipKernelSplineOrder_ = pmeMutable.splineOrder;

    pmeMutable.pmeGrid = std::move(pmeGridBackup);
}

void GCMCEngine::ensurePgpHostGridReady() {
    ensurePgpInitialized();
    if (pgpHostGridReady_) {
        return;
    }
    if (!state_) {
        throw std::runtime_error("PGP backend requires MCState");
    }
    if (!hasAnyActiveFixedResidue(*state_)) {
        throw std::runtime_error(
            "energy_method=pgp_host requires at least one active fixed residue (host/framework)");
    }
    precomputeGridPotential(*state_, true);
    pgpHostGridReady_ = true;
}

void GCMCEngine::ensurePgpFullGridReadyExcluding(int excludedResidueIdx) {
    ensurePgpInitialized();
    if (!state_) {
        throw std::runtime_error("PGP backend requires MCState");
    }

    int exclude = excludedResidueIdx;
    if (exclude < 0 || exclude >= state_->activeResidueCount) {
        exclude = -1;
    }

    if (pgpFullGridReady_ && pgpFullGridExcludedResidue_ == exclude) {
        return;
    }

    // Fast path: if the background (all residues except the excluded one) is empty,
    // the reciprocal potential grid is identically zero. Skipping the full
    // precomputeGridPotential() avoids unnecessary FFT work and prevents any chance
    // of transient PME-state mutations affecting later diagnostics.
    bool hasAnyBackgroundCharge = false;
    for (int r = 0; r < state_->activeResidueCount; ++r) {
        if (r == exclude) continue;
        if (state_->residues[r].active) {
            hasAnyBackgroundCharge = true;
            break;
        }
    }
    if (!hasAnyBackgroundCharge) {
        auto& pgp = getPGPParams();
        std::fill(pgp.potentialGrid.begin(), pgp.potentialGrid.end(), std::complex<double>(0.0, 0.0));
        pgpFullGridReady_ = true;
        pgpFullGridExcludedResidue_ = exclude;
        return;
    }

    if (exclude >= 0) {
        ResidueActiveGuard guard(*state_, exclude);
        guard.setInactive();
        precomputeGridPotential(*state_, false);
    } else {
        precomputeGridPotential(*state_, false);
    }

    pgpFullGridReady_ = true;
    pgpFullGridExcludedResidue_ = exclude;
}

void GCMCEngine::computePgpFullGridExcludingNoCache(int excludedResidueIdx) {
    ensurePgpInitialized();
    if (!state_) {
        throw std::runtime_error("PGP backend requires MCState");
    }

    int exclude = excludedResidueIdx;
    if (exclude < 0 || exclude >= state_->activeResidueCount) {
        exclude = -1;
    }

    if (exclude >= 0) {
        ResidueActiveGuard guard(*state_, exclude);
        guard.setInactive();
        precomputeGridPotential(*state_, false);
    } else {
        precomputeGridPotential(*state_, false);
    }
}

double GCMCEngine::calculatePgpReciprocalEnergyFromGrid(int residueIdx) {
    if (!state_) {
        return 0.0;
    }
    MovementResiduesGuard movementGuard(*state_);
    movementGuard.setSingleResidue(residueIdx);
    double gridEnergy = 0.0;
    interpolateMoleculeEnergy(*state_, gridEnergy);
    return gridEnergy;
}

double GCMCEngine::calculatePgpReciprocalMeshSelfEnergy(int residueIdx) {
    if (!state_) {
        return 0.0;
    }
    const auto& pmeParams = getPMEParams();
    const int nx = pmeParams.meshSize[0];
    const int ny = pmeParams.meshSize[1];
    const int nz = pmeParams.meshSize[2];
    if (!(nx > 0 && ny > 0 && nz > 0)) {
        return 0.0;
    }
    const int totalGridSize = nx * ny * nz;
    if (!(totalGridSize > 0)) {
        return 0.0;
    }
    if (residueIdx < 0 || residueIdx >= static_cast<int>(state_->residues.size())) {
        return 0.0;
    }
    const auto& residue = state_->residues[residueIdx];
    if (!residue.active) {
        return 0.0;
    }

    const float* box = state_->info.box;
    if (!(box[0] > 0.0f && box[1] > 0.0f && box[2] > 0.0f)) {
        return 0.0;
    }

    // Mode E (energy_method=pgp_full_pme) needs the exact PME mesh-self term so that per-move ΔU
    // matches `energy_method=pme`. We compute it by running the same PME reciprocal pipeline as
    // `computeReciprocalPME`, but restricting charge spreading to this residue only.
    //
    // This is intentionally slower than the sparse-kernel form; it is only used in strict
    // PME-compat mode and avoids subtle charge-spreading/normalization mismatches.

    std::vector<char> residueActiveBackup;
    residueActiveBackup.reserve(static_cast<size_t>(state_->activeResidueCount));
    for (int r = 0; r < state_->activeResidueCount; ++r) {
        residueActiveBackup.push_back(static_cast<char>(state_->residues[r].active ? 1 : 0));
    }
    for (int r = 0; r < state_->activeResidueCount; ++r) {
        if (r == residueIdx) continue;
        state_->residues[r].active = false;
    }

    auto& pmeMutable = getPMEParams();
    std::vector<std::complex<double>> gridBackup = std::move(pmeMutable.pmeGrid);
    pmeMutable.pmeGrid.assign(static_cast<size_t>(totalGridSize), std::complex<double>(0.0, 0.0));
    spreadChargesOntoGrid(*state_, false);

    performFFTForward();

    const double boxNm[3] = {
        static_cast<double>(box[0]),
        static_cast<double>(box[1]),
        static_cast<double>(box[2]),
    };
    pmeMutable.setBox(boxNm);

    double recipEnergy = 0.0;
    computeEnergyFromGrid(recipEnergy, boxNm);

    pmeMutable.pmeGrid = std::move(gridBackup);

    for (int r = 0; r < state_->activeResidueCount && r < static_cast<int>(residueActiveBackup.size()); ++r) {
        state_->residues[r].active = (residueActiveBackup[static_cast<size_t>(r)] != 0);
    }
    return recipEnergy;
}

double GCMCEngine::calculatePgpSelfEnergyMovementResidue(int residueIdx) {
    if (!state_) {
        return 0.0;
    }
    MovementResiduesGuard movementGuard(*state_);
    movementGuard.setSingleResidue(residueIdx);
    return computeSelfEnergyPGPImpl(*state_, true);
}

double GCMCEngine::calculatePgpRealSpaceElectrostaticsFixedOnly(int residueIdx) {
    if (!state_) {
        return 0.0;
    }
    if (residueIdx < 0 || residueIdx >= static_cast<int>(state_->residues.size())) {
        return 0.0;
    }
    const auto& residues = state_->residues;
    const auto& atoms = state_->atoms;
    const auto& resI = residues[residueIdx];
    if (!resI.active) {
        return 0.0;
    }

    const double cutoff2 = static_cast<double>(state_->info.cutoff) * static_cast<double>(state_->info.cutoff);
    const double alpha = getPGPParams().alpha;
    const float* box = state_->info.box;

    double energy = 0.0;
    for (int atom_i = resI.atomStart; atom_i < resI.atomStart + resI.atomCount; ++atom_i) {
        if (atom_i < 0 || atom_i >= static_cast<int>(atoms.size())) continue;
        const auto& ai = atoms[atom_i];
        if (ai.name == "LP" || ai.name == "LPA") continue;
        const double qi = static_cast<double>(ai.charge);
        if (std::abs(qi) < 1e-12) continue;

        const double xi = static_cast<double>(ai.x);
        const double yi = static_cast<double>(ai.y);
        const double zi = static_cast<double>(ai.z);

        for (int r = 0; r < state_->activeResidueCount; ++r) {
            const auto& resJ = residues[r];
            if (!resJ.active || !resJ.fixed) continue;

            for (int atom_j = resJ.atomStart; atom_j < resJ.atomStart + resJ.atomCount; ++atom_j) {
                if (atom_j < 0 || atom_j >= static_cast<int>(atoms.size())) continue;
                const auto& aj = atoms[atom_j];
                if (aj.name == "LP" || aj.name == "LPA") continue;
                const double qj = static_cast<double>(aj.charge);
                if (std::abs(qj) < 1e-12) continue;

                double dx = static_cast<double>(aj.x) - xi;
                double dy = static_cast<double>(aj.y) - yi;
                double dz = static_cast<double>(aj.z) - zi;

                dx -= static_cast<double>(box[0]) * std::round(dx / static_cast<double>(box[0]));
                dy -= static_cast<double>(box[1]) * std::round(dy / static_cast<double>(box[1]));
                dz -= static_cast<double>(box[2]) * std::round(dz / static_cast<double>(box[2]));

                const double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > cutoff2 || r2 < 1e-12) continue;

                const double rDist = std::sqrt(r2);
                energy += COULOMB * qi * qj * std::erfc(alpha * rDist) / rDist;
            }
        }
    }
    return energy;
}

double GCMCEngine::calculatePgpRealSpaceElectrostaticsAllPartners(int residueIdx) {
    if (!state_) {
        return 0.0;
    }
    if (residueIdx < 0 || residueIdx >= static_cast<int>(state_->residues.size())) {
        return 0.0;
    }
    const auto& residues = state_->residues;
    const auto& atoms = state_->atoms;
    const auto& resI = residues[residueIdx];
    if (!resI.active) {
        return 0.0;
    }

    const double cutoff2 = static_cast<double>(state_->info.cutoff) * static_cast<double>(state_->info.cutoff);
    const double alpha = getPGPParams().alpha;
    const float* box = state_->info.box;

    double energy = 0.0;
    for (int atom_i = resI.atomStart; atom_i < resI.atomStart + resI.atomCount; ++atom_i) {
        if (atom_i < 0 || atom_i >= static_cast<int>(atoms.size())) continue;
        const auto& ai = atoms[atom_i];
        if (ai.name == "LP" || ai.name == "LPA") continue;
        const double qi = static_cast<double>(ai.charge);
        if (std::abs(qi) < 1e-12) continue;

        const double xi = static_cast<double>(ai.x);
        const double yi = static_cast<double>(ai.y);
        const double zi = static_cast<double>(ai.z);

        for (int r = 0; r < state_->activeResidueCount; ++r) {
            if (r == residueIdx) continue;
            const auto& resJ = residues[r];
            if (!resJ.active) continue;

            for (int atom_j = resJ.atomStart; atom_j < resJ.atomStart + resJ.atomCount; ++atom_j) {
                if (atom_j < 0 || atom_j >= static_cast<int>(atoms.size())) continue;
                const auto& aj = atoms[atom_j];
                if (aj.name == "LP" || aj.name == "LPA") continue;
                const double qj = static_cast<double>(aj.charge);
                if (std::abs(qj) < 1e-12) continue;

                double dx = static_cast<double>(aj.x) - xi;
                double dy = static_cast<double>(aj.y) - yi;
                double dz = static_cast<double>(aj.z) - zi;

                dx -= static_cast<double>(box[0]) * std::round(dx / static_cast<double>(box[0]));
                dy -= static_cast<double>(box[1]) * std::round(dy / static_cast<double>(box[1]));
                dz -= static_cast<double>(box[2]) * std::round(dz / static_cast<double>(box[2]));

                const double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > cutoff2 || r2 < 1e-12) continue;

                const double rDist = std::sqrt(r2);
                energy += COULOMB * qi * qj * std::erfc(alpha * rDist) / rDist;
            }
        }
    }
    return energy;
}

double GCMCEngine::calculateFragmentEnergyPgpHost(int residueIdx) {
    if (!state_) {
        return 0.0;
    }
    ensurePgpHostGridReady();

    if (residueIdx < 0 || residueIdx >= static_cast<int>(state_->residues.size())) {
        return 0.0;
    }
    auto& residue = state_->residues[residueIdx];
    if (!residue.active) {
        return 0.0;
    }

    // 1) Direct (cutoff) interactions against non-fixed residues (guest↔guest), incl. intramolecular 1-4 LJ once.
    computeResidueNonbondedEnergy(
        *state_,
        residueIdx,
        true,
        true,
        false,
        true,
        ResiduePartnerFilter::NonFixedOnly
    );
    const double vdwNonFixed = residue.energy_vdw;
    const double elecNonFixed = residue.energy_elec;

    // 2) LJ interactions against fixed residues (host↔guest) using the same LJ backend; skip intra 1-4 here.
    computeResidueNonbondedEnergy(
        *state_,
        residueIdx,
        true,
        true,
        true,
        false,
        ResiduePartnerFilter::FixedOnly
    );
    const double vdwFixed = residue.energy_vdw;

    // 3) Real-space erfc term against fixed residues only.
    const double elecFixedReal = calculatePgpRealSpaceElectrostaticsFixedOnly(residueIdx);

    // 4) Reciprocal term from the fixed-only potential grid.
    const double elecFixedRecip = calculatePgpReciprocalEnergyFromGrid(residueIdx);

    const double vdwTotal = vdwNonFixed + vdwFixed;
    const double elecTotal = elecNonFixed + elecFixedReal + elecFixedRecip;

    residue.energy_vdw = static_cast<float>(vdwTotal);
    residue.energy_elec = static_cast<float>(elecTotal);

    return vdwTotal + elecTotal;
}

double GCMCEngine::calculateFragmentEnergyPgpFullUsingCurrentGrid(int residueIdx) {
    if (!state_) {
        return 0.0;
    }
    if (residueIdx < 0 || residueIdx >= static_cast<int>(state_->residues.size())) {
        return 0.0;
    }
    auto& residue = state_->residues[residueIdx];
    if (!residue.active) {
        return 0.0;
    }

    // LJ: always DIRECT+cutoff against all partners (host+guests), incl. intramolecular 1-4 LJ.
    computeResidueNonbondedEnergy(
        *state_,
        residueIdx,
        true,
        true,
        true,
        true,
        ResiduePartnerFilter::All
    );
    const double vdwTotal = residue.energy_vdw;

    const double elecReal = calculatePgpRealSpaceElectrostaticsAllPartners(residueIdx);
    const double elecRecip = calculatePgpReciprocalEnergyFromGrid(residueIdx);
    const bool includeMeshSelf = (energyBackend_ == GCMCEnergyBackend::PgpFullPme);
    const double elecRecipSelf = includeMeshSelf ? calculatePgpReciprocalMeshSelfEnergy(residueIdx) : 0.0;
    const double elecSelf = calculatePgpSelfEnergyMovementResidue(residueIdx);

    const double elecTotal = elecReal + elecRecip + elecRecipSelf + elecSelf;
    residue.energy_elec = static_cast<float>(elecTotal);

    return vdwTotal + elecTotal;
}

// Attempt insertion
GCMCEngine::MoveResult GCMCEngine::attemptInsertion(int typeId) {
    MoveResult result;
    result.type = MoveResult::INSERT;
    result.fragmentType = typeId;

    // Diagnostics: CBMC trial energies for this move (only populated when CBMC is used).
    lastCbmcTrialEnergies_.clear();
    
    if (!state_ || !reservoir_) {
        result.accepted = false;
        return result;
    }
    
    // Get template
    FragmentTemplate* tmpl = reservoir_->getTemplate(typeId);
    if (!tmpl) {
        result.accepted = false;
        return result;
    }
    
    double baseVolume = getBoxVolume();
    if (baseVolume <= 0.0) {
        baseVolume = 1.0;
    }
    if (acceptanceCalculator_) {
        double configuredVolume = acceptanceCalculator_->getVolume();
        if (configuredVolume > 0.0) {
            baseVolume = configuredVolume;
        } else {
            acceptanceCalculator_->setVolume(baseVolume);
            baseVolume = acceptanceCalculator_->getVolume();
            if (baseVolume <= 0.0) {
                baseVolume = 1.0;
            }
        }
    }
    result.effectiveVolume = baseVolume;

    // CRITICAL FIX: Get N BEFORE insertion for correct acceptance calculation
    int N_before = reservoir_->getActiveCount(typeId);

    // Configurable capacity check to prevent runaway growth
    // Can be disabled by setting maxMoleculesPerType to -1
    double maxMolecules = getConfigValue("maxMoleculesPerType");
    if (maxMolecules <= 0) {
        // Default: large but reasonable limit for safety
        maxMolecules = 10000;
    }

    if (N_before >= static_cast<int>(maxMolecules)) {
        // At limit, reject insertion immediately
        MoveResult early;
        early.type = MoveResult::INSERT;
        early.fragmentType = typeId;
        early.accepted = false;
        early.deltaE = 0.0;
        early.energyBefore = 0.0;
        early.energyAfter = 0.0;
        early.bias = 0.0;
        early.acceptanceProbability = 0.0;
        totalMoves_++;
        return early;
    }

    // Determine number of CBMC trials for this fragment type
    int numTrials = 1;
    if (useConfBias_ && typeId < static_cast<int>(cbmcTrialsPerType_.size())) {
        numTrials = cbmcTrialsPerType_[typeId];
    }

    Vector3 position;
    Quaternion orientation;
    double cbmcRosen = 1.0;
    double cavityVolumeFraction = 1.0;
    int trialsUsed = 1;

    if (useConfBias_ && numTrials > 1) {
        // Use CBMC to select configuration
        TrialConfiguration selected = performCBMCInsertion(typeId, numTrials);
        if (!selected.valid) {
            MoveResult early;
            early.type = MoveResult::INSERT;
            early.fragmentType = typeId;
            early.accepted = false;
            early.deltaE = 0.0;
            early.energyBefore = 0.0;
            early.energyAfter = 0.0;
            early.bias = 0.0;
            early.cbmcTrialsUsed = numTrials;
            early.rosenbluthWeight = 0.0;
            early.acceptanceProbability = shouldStoreProbability() ? 0.0 : -1.0;
            totalMoves_++;
            return early;
        }
        position = selected.position;
        orientation = selected.orientation;
        cbmcRosen = selected.weight;
        result.cbmcSelectedEnergy = selected.energy;
        result.cbmcLogWOverK = selected.logWOverK;
        if (selected.trialsUsed > 0) {
            trialsUsed = selected.trialsUsed;
        } else {
            trialsUsed = numTrials;
        }
    } else {
        // Original single configuration generation
        // Try multiple times to find a placement fully inside region (if configured)
        const int maxTrials = 50;
        int tries = 0;
        do {
            position = (useCavityBias_ && cavityManager_) ?
                      generateCavityPosition() : generateRandomPosition();
            applyPeriodicBoundary(position);  // Ensure position is within PBC
            orientation = generateRandomOrientation();
            tries++;
        } while (regionConstraint_ && !isMoleculeWithinRegion(typeId, position, orientation) && tries < maxTrials);

        // If no valid placement found, reject early (hard region constraint)
        if (regionConstraint_ && !isMoleculeWithinRegion(typeId, position, orientation)) {
            MoveResult early;
            early.type = MoveResult::INSERT;
            early.fragmentType = typeId;
            early.accepted = false;
            early.position = position;
            early.deltaE = 0.0;
            early.energyBefore = 0.0;
            early.energyAfter = 0.0;
            early.bias = cbmcRosen;
            early.acceptanceProbability = shouldStoreProbability() ? 0.0 : -1.0;
            totalMoves_++;
            return early;
        }
    }

    if (useCavityBias_ && cavityManager_) {
        cavityVolumeFraction = std::max(1e-12, cavityManager_->getCavityVolumeFraction(*state_));
    }
    result.effectiveVolume = baseVolume * cavityVolumeFraction;

    result.position = position;
    
    // Precompute any PGP background grid before adding the new residue, so the grid does not
    // accidentally include the trial/inserted charges.
    if (energyBackend_ == GCMCEnergyBackend::PgpHost) {
        ensurePgpHostGridReady();
    } else if (energyBackend_ == GCMCEnergyBackend::PgpFull ||
               energyBackend_ == GCMCEnergyBackend::PgpFullPme) {
        ensurePgpFullGridReadyExcluding(-1);
    }

    const bool useSystemEnergyDelta =
        useDrude_ ||
        energyBackend_ == GCMCEnergyBackend::Ewald ||
        energyBackend_ == GCMCEnergyBackend::Pme;

    // Calculate energy before insertion (only needed for full-system ΔU paths).
    if (useSystemEnergyDelta) {
        result.energyBefore = calculateSystemEnergy();
    } else {
        result.energyBefore = 0.0;
    }
    
    // Create instance
    int instanceId = reservoir_->createInstance(typeId, position, orientation);
    if (instanceId < 0) {
        result.accepted = false;
        return result;
    }
    
    // Synchronize MCState with the new instance
    synchronizeStateWithReservoir(instanceId, true);
    
    // Calculate energy change using the selected backend.
    if (useSystemEnergyDelta) {
        result.energyAfter = calculateSystemEnergy();
        result.deltaE = result.energyAfter - result.energyBefore;
    } else if (energyBackend_ == GCMCEnergyBackend::DirectCutoff) {
        // Fast local ΔE calculation - only compute interaction of new residue with system.
        cpu::computeResidueEnergyCutoffPBC(*state_, instanceId);
        const auto& residue = state_->residues[instanceId];
        result.deltaE = residue.energy_vdw + residue.energy_elec;
        result.energyAfter = result.deltaE;  // For consistency
    } else if (energyBackend_ == GCMCEnergyBackend::PgpHost) {
        result.deltaE = calculateFragmentEnergyPgpHost(instanceId);
        result.energyAfter = result.deltaE;
    } else if (energyBackend_ == GCMCEnergyBackend::PgpFull ||
               energyBackend_ == GCMCEnergyBackend::PgpFullPme) {
        // Use the precomputed background grid (exclude=-1) prepared above / in CBMC.
        result.deltaE = calculateFragmentEnergyPgpFullUsingCurrentGrid(instanceId);
        result.energyAfter = result.deltaE;
    } else {
        // Should be unreachable: non-local backends are covered by useSystemEnergyDelta.
        result.energyAfter = calculateSystemEnergy();
        result.deltaE = result.energyAfter - result.energyBefore;
    }
    result.cbmcTrialsUsed = std::max(trialsUsed, 1);

    // Calculate scheduler/config bias components (exclude cavity handled separately)
    double schedulerBias = calculateInsertionBias(*tmpl, position, orientation);

    // Store individual components for detailed balance verification
    double cavityFraction = cavityVolumeFraction;
    result.rosenbluthWeight = cbmcRosen;
    result.cavityBiasComponent = cavityFraction;
    result.bias = schedulerBias * cbmcRosen;

    // Calculate acceptance probability using proper GCMC formula
    bool accept = false;
    double prob = 0.0;
    if (acceptanceCalculator_) {
        double proposalBias = getConfigValue("proposalBias");
        double proposalLogRatio = 0.0;
        if (proposalBias > 0.0) {
            proposalLogRatio = std::log(proposalBias);
        }
        GCMCAcceptance::GrandCanonicalInsertionTerms terms;
        terms.typeId = typeId;
        terms.countBefore = N_before;
        terms.deltaE = result.deltaE;
        terms.cavityFraction = cavityFraction;
        terms.lambdaNm = acceptanceCalculator_->getThermalLambda(typeId);
        terms.rosenbluthWeight = std::max(cbmcRosen, 1e-30);
        terms.cbmcTrials = result.cbmcTrialsUsed;
        terms.proposalLogRatio = proposalLogRatio;
        prob = acceptanceCalculator_->calculateInsertionProbabilityDetailed(terms);
        
        // Only store probability if configured
        result.acceptanceProbability = shouldStoreProbability() ? prob : -1.0;
        
        accept = acceptanceCalculator_->acceptMove(prob);
    } else {
        throw std::runtime_error(
            "GCMCEngine::attemptInsertion requires an acceptance calculator; "
            "log-space fallback acceptance has been removed.");
    }
    
    if (accept) {
        result.accepted = true;
        result.residueIndex = instanceId;
        acceptedMoves_++;
        energyCache_.invalidate();
        if (energyBackend_ == GCMCEnergyBackend::PgpFull ||
            energyBackend_ == GCMCEnergyBackend::PgpFullPme) {
            // Background changed; force a rebuild on the next evaluation.
            pgpFullGridReady_ = false;
            pgpFullGridExcludedResidue_ = -999;
        }
    } else {
        // Remove the instance and revert state
        reservoir_->deleteInstance(instanceId);
        synchronizeStateWithReservoir(instanceId, false);
        relaxDrudeIfEnabled();
        result.accepted = false;
    }
    
    totalMoves_++;
    
    // Sample statistics if configured using a lightweight countdown
    if (collectStats_) {
        if (--statsCountdown_ <= 0) {
            int particleCount = reservoir_ ? reservoir_->getActiveCount() : 0;
            double energy = calculateSystemEnergy();
            statistics_.addSample(totalMoves_, particleCount, energy,
                                  getAcceptanceRate(), temperature_, 100.0);
            statsCountdown_ = std::max(1, statsInterval_);
        }
    }
    
    return result;
}

// Attempt deletion
GCMCEngine::MoveResult GCMCEngine::attemptDeletion(int typeId) {
    MoveResult result;
    result.type = MoveResult::DELETE;
    result.fragmentType = typeId;

    // Diagnostics: CBMC trial energies for this move (only populated when CBMC is used).
    lastCbmcTrialEnergies_.clear();
    
    if (!state_ || !reservoir_) {
        result.accepted = false;
        return result;
    }
    
    // CRITICAL FIX: Get N BEFORE deletion for correct acceptance calculation
    int N_before = reservoir_->getActiveCount(typeId);
    
    // Check if any instances exist
    if (N_before == 0) {
        result.accepted = false;
        // Set probability semantics: 0.0 when storing, -1.0 when not
        result.acceptanceProbability = shouldStoreProbability() ? 0.0 : -1.0;
        return result;
    }
    
    // Select random instance of this type
    int instanceId = selectRandomInstance(typeId);
    if (instanceId < 0) {
        result.accepted = false;
        // Consistent probability storage: 0 when N=0, -1 when disabled
        result.acceptanceProbability = shouldStoreProbability() ? 0.0 : -1.0;
        return result;
    }
    
    result.residueIndex = instanceId;
    
    // Get instance info before deletion
    FragmentInstance* instance = reservoir_->getInstance(instanceId);
    if (!instance) {
        result.accepted = false;
        return result;
    }
    
    // Save position and orientation for potential restoration
    Vector3 savedPosition = instance->position;
    Quaternion savedOrientation = instance->orientation;
    result.position = savedPosition;

    double baseVolume = getBoxVolume();
    if (baseVolume <= 0.0) {
        baseVolume = 1.0;
    }
    if (acceptanceCalculator_) {
        double configuredVolume = acceptanceCalculator_->getVolume();
        if (configuredVolume > 0.0) {
            baseVolume = configuredVolume;
        } else {
            acceptanceCalculator_->setVolume(baseVolume);
            baseVolume = acceptanceCalculator_->getVolume();
            if (baseVolume <= 0.0) {
                baseVolume = 1.0;
            }
        }
    }
    result.effectiveVolume = baseVolume;

    auto residueEnergyForBackend = [&](int residueIdx) -> double {
        if (!state_) {
            return 0.0;
        }
        if (energyBackend_ == GCMCEnergyBackend::DirectCutoff) {
            cpu::computeResidueEnergyCutoffPBC(*state_, residueIdx);
            const auto& residue = state_->residues[residueIdx];
            return residue.energy_vdw + residue.energy_elec;
        }
        if (energyBackend_ == GCMCEnergyBackend::PgpHost) {
            return calculateFragmentEnergyPgpHost(residueIdx);
        }
        if (energyBackend_ == GCMCEnergyBackend::PgpFull ||
            energyBackend_ == GCMCEnergyBackend::PgpFullPme) {
            ensurePgpFullGridReadyExcluding(residueIdx);
            return calculateFragmentEnergyPgpFullUsingCurrentGrid(residueIdx);
        }
        return calculateFragmentEnergy(residueIdx);
    };

    // Calculate CBMC bias for deletion if enabled
    double cbmcRosen = 1.0;
    double cbmcSelectedEnergy = 0.0;
    double cbmcLogWOverK = 0.0;
    double cavityVolumeFraction = 1.0;
    int numTrials = 1;
    if (useConfBias_ && typeId < static_cast<int>(cbmcTrialsPerType_.size())) {
        numTrials = cbmcTrialsPerType_[typeId];
    }

    if (useConfBias_ && numTrials > 1) {
        if (energyBackend_ == GCMCEnergyBackend::Pme) {
            // Mode D (PME): compute CBMC retracing weights using full-system PME energy differences.
            // This is intentionally slow and serves as a validation/reference backend.
            std::vector<TrialConfiguration> trials;
            trials.reserve(numTrials);

            const double energyWithCurrent = calculateSystemEnergy();

            // Temporarily deactivate this residue in MCState so trial energies do not
            // include interactions with the molecule being deleted.
            bool restoredActive = false;
            int restoredAtomCount = 0;
            if (state_ && instanceId >= 0 && instanceId < static_cast<int>(state_->residues.size())) {
                auto& residue = state_->residues[instanceId];
                restoredActive = residue.active;
                restoredAtomCount = residue.atomCount;
                residue.active = false;
                residue.atomCount = 0;
            }

            const double baselineEnergy = calculateSystemEnergy();

            TrialConfiguration current;
            current.weight = 0.0;
            current.logWOverK = 0.0;
            current.position = savedPosition;
            current.orientation = savedOrientation;
            current.energy = energyWithCurrent - baselineEnergy;
            trials.push_back(current);
            cbmcSelectedEnergy = current.energy;

            FragmentTemplate* tmpl = reservoir_->getTemplate(typeId);
            if (tmpl) {
                for (int k = 1; k < numTrials; ++k) {
                    TrialConfiguration trial;
                    trial.weight = 0.0;
                    trial.logWOverK = 0.0;
                    trial.position = (useCavityBias_ && cavityManager_) ?
                                    generateCavityPosition() : generateRandomPosition();
                    applyPeriodicBoundary(trial.position);
                    trial.orientation = generateRandomOrientation();

                    int tempId = reservoir_->createInstance(typeId, trial.position, trial.orientation);
                    if (tempId >= 0) {
                        synchronizeStateWithReservoir(tempId, true);
                        const double energyWithTrial = calculateSystemEnergy();
                        trial.energy = energyWithTrial - baselineEnergy;

                        reservoir_->deleteInstance(tempId);
                        synchronizeStateWithReservoir(tempId, false);
                        trials.push_back(trial);
                    }
                }
            }

            // Restore residue activity for the actual deletion attempt.
            if (state_ && instanceId >= 0 && instanceId < static_cast<int>(state_->residues.size())) {
                auto& residue = state_->residues[instanceId];
                residue.active = restoredActive;
                residue.atomCount = restoredAtomCount;
            }

            if (trials.size() == static_cast<size_t>(numTrials)) {
                lastCbmcTrialEnergies_.clear();
                lastCbmcTrialEnergies_.reserve(trials.size());
                for (const auto& trial : trials) {
                    lastCbmcTrialEnergies_.push_back(trial.energy);
                }

                const double beta = 1.0 / (8.314e-3 * temperature_);
                double minEnergy = std::numeric_limits<double>::max();
                for (const auto& trial : trials) {
                    minEnergy = std::min(minEnergy, trial.energy);
                }

                double sumScaled = 0.0;
                for (const auto& trial : trials) {
                    sumScaled += std::exp(-beta * (trial.energy - minEnergy));
                }

                const double avgScaled = sumScaled / static_cast<double>(numTrials);
                const double safeAvgScaled = std::max(avgScaled, 1e-30);
                cbmcLogWOverK = std::log(safeAvgScaled) - beta * minEnergy;

                // rosen = (W/K)/exp(-β u_current) = exp(log(W/K) + β u_current)
                const double logLower = std::log(1e-30);
                const double logUpper = 700.0;
                const double logRosen = cbmcLogWOverK + beta * current.energy;
                const double logRosenClamped = std::min(std::max(logRosen, logLower), logUpper);
                cbmcRosen = std::exp(logRosenClamped);
            }
        } else {
        // For deletion, compute CBMC retracing weights in the post-deletion environment
        // (i.e., excluding the current molecule) to match the reverse insertion proposal.
        std::vector<TrialConfiguration> trials;
        trials.reserve(numTrials);

        TrialConfiguration current;
        current.weight = 0.0;
        current.logWOverK = 0.0;
        current.position = savedPosition;
        current.orientation = savedOrientation;
        current.energy = residueEnergyForBackend(instanceId);
        trials.push_back(current);
        cbmcSelectedEnergy = current.energy;

        // Temporarily deactivate this residue in MCState so trial energies do not
        // include interactions with the molecule being deleted.
        bool restoredActive = false;
        int restoredAtomCount = 0;
        if (state_ && instanceId >= 0 && instanceId < static_cast<int>(state_->residues.size())) {
            auto& residue = state_->residues[instanceId];
            restoredActive = residue.active;
            restoredAtomCount = residue.atomCount;
            residue.active = false;
            residue.atomCount = 0;
        }

        FragmentTemplate* tmpl = reservoir_->getTemplate(typeId);
        if (tmpl) {
            for (int k = 1; k < numTrials; ++k) {
                TrialConfiguration trial;
                trial.weight = 0.0;
                trial.logWOverK = 0.0;
                trial.position = (useCavityBias_ && cavityManager_) ?
                                generateCavityPosition() : generateRandomPosition();
                applyPeriodicBoundary(trial.position);
                trial.orientation = generateRandomOrientation();

                // Create temporary instance for energy calculation
                int tempId = reservoir_->createInstance(typeId, trial.position, trial.orientation);
                if (tempId >= 0) {
                    synchronizeStateWithReservoir(tempId, true);

                    if (energyBackend_ == GCMCEnergyBackend::DirectCutoff) {
                        cpu::computeResidueEnergyCutoffPBC(*state_, tempId);
                        const auto& residue = state_->residues[tempId];
                        trial.energy = residue.energy_vdw + residue.energy_elec;
                    } else if (energyBackend_ == GCMCEnergyBackend::PgpHost) {
                        trial.energy = calculateFragmentEnergyPgpHost(tempId);
                    } else if (energyBackend_ == GCMCEnergyBackend::PgpFull ||
                               energyBackend_ == GCMCEnergyBackend::PgpFullPme) {
                        // Reuse the background grid built excluding the molecule being deleted.
                        trial.energy = calculateFragmentEnergyPgpFullUsingCurrentGrid(tempId);
                    } else {
                        trial.energy = calculateFragmentEnergy(tempId);
                    }

                    reservoir_->deleteInstance(tempId);
                    synchronizeStateWithReservoir(tempId, false);
                    trials.push_back(trial);
                }
            }
        }

        // Restore residue activity for the actual deletion attempt.
        if (state_ && instanceId >= 0 && instanceId < static_cast<int>(state_->residues.size())) {
            auto& residue = state_->residues[instanceId];
            residue.active = restoredActive;
            residue.atomCount = restoredAtomCount;
        }

        if (trials.size() == static_cast<size_t>(numTrials)) {
            lastCbmcTrialEnergies_.clear();
            lastCbmcTrialEnergies_.reserve(trials.size());
            for (const auto& trial : trials) {
                lastCbmcTrialEnergies_.push_back(trial.energy);
            }

            const double beta = 1.0 / (8.314e-3 * temperature_);
            double minEnergy = std::numeric_limits<double>::max();
            for (const auto& trial : trials) {
                minEnergy = std::min(minEnergy, trial.energy);
            }

            double sumScaled = 0.0;
            for (const auto& trial : trials) {
                sumScaled += std::exp(-beta * (trial.energy - minEnergy));
            }

            const double avgScaled = sumScaled / static_cast<double>(numTrials);
            const double safeAvgScaled = std::max(avgScaled, 1e-30);
            cbmcLogWOverK = std::log(safeAvgScaled) - beta * minEnergy;

            // rosen = (W/K)/exp(-β u_current) = exp(log(W/K) + β u_current)
            const double logLower = std::log(1e-30);
            const double logUpper = 700.0;
            const double logRosen = cbmcLogWOverK + beta * current.energy;
            const double logRosenClamped = std::min(std::max(logRosen, logLower), logUpper);
            cbmcRosen = std::exp(logRosenClamped);
        }
        }
    }

    const bool cbmcEnabled = (useConfBias_ && numTrials > 1);
    const bool useSystemEnergyDelta =
        useDrude_ ||
        energyBackend_ == GCMCEnergyBackend::Ewald ||
        energyBackend_ == GCMCEnergyBackend::Pme;

    if (!useSystemEnergyDelta &&
        (energyBackend_ == GCMCEnergyBackend::DirectCutoff ||
         energyBackend_ == GCMCEnergyBackend::PgpHost ||
         energyBackend_ == GCMCEnergyBackend::PgpFull ||
         energyBackend_ == GCMCEnergyBackend::PgpFullPme)) {
        const double residueEnergy = cbmcEnabled ? cbmcSelectedEnergy : residueEnergyForBackend(instanceId);
        result.deltaE = -residueEnergy;  // Removing this energy from system
        result.energyBefore = residueEnergy;
        result.energyAfter = 0.0;
    } else {
        // Full system energy for Drude/EWALD/PME modes.
        result.energyBefore = calculateSystemEnergy();
    }
    
    // Temporarily delete (convert to ghost)
    reservoir_->deleteInstance(instanceId);
    synchronizeStateWithReservoir(instanceId, false);
    
    // Calculate energy after deletion for full-system modes.
    if (useSystemEnergyDelta) {
        result.energyAfter = calculateSystemEnergy();
        result.deltaE = result.energyAfter - result.energyBefore;
    }
    
    // Calculate cavity bias component separately for detailed balance tracking
    if (useCavityBias_ && cavityManager_) {
        cavityVolumeFraction = std::max(1e-12, cavityManager_->getCavityVolumeFraction(*state_));
    }
    result.effectiveVolume = baseVolume * cavityVolumeFraction;

    const int trialsUsed = std::max(numTrials, 1);
    result.cbmcTrialsUsed = trialsUsed;

    // CRITICAL FIX: Calculate deletion bias using saved position for robustness
    // Include CBMC bias in total bias (cavity handled separately)
    double schedulerBias = calculateDeletionBiasAtPosition(savedPosition);
    result.bias = schedulerBias * cbmcRosen;

    // Store individual components for detailed balance verification
    double cavityFraction = cavityVolumeFraction;
    result.rosenbluthWeight = cbmcRosen;
    result.cbmcSelectedEnergy = cbmcSelectedEnergy;
    result.cbmcLogWOverK = cbmcLogWOverK;
    result.cavityBiasComponent = cavityFraction;

    // Calculate acceptance probability using proper GCMC formula
    bool accept = false;
    double prob = 0.0;
    if (acceptanceCalculator_) {
        double proposalBias = getConfigValue("proposalBias");
        double proposalLogRatio = 0.0;
        if (proposalBias > 0.0) {
            proposalLogRatio = std::log(proposalBias);
        }
        GCMCAcceptance::GrandCanonicalDeletionTerms terms;
        terms.typeId = typeId;
        terms.countBefore = N_before;
        terms.deltaE = result.deltaE;
        terms.cavityFraction = cavityFraction;
        terms.lambdaNm = acceptanceCalculator_->getThermalLambda(typeId);
        terms.rosenbluthWeight = std::max(cbmcRosen, 1e-30);
        terms.cbmcTrials = trialsUsed;
        terms.proposalLogRatio = proposalLogRatio;
        prob = acceptanceCalculator_->calculateDeletionProbabilityDetailed(terms);
        
        // Only store probability if configured
        result.acceptanceProbability = shouldStoreProbability() ? prob : -1.0;
        
        accept = acceptanceCalculator_->acceptMove(prob);
    } else {
        throw std::runtime_error(
            "GCMCEngine::attemptDeletion requires an acceptance calculator; "
            "log-space fallback acceptance has been removed.");
    }
    
    if (accept) {
        result.accepted = true;
        acceptedMoves_++;
        energyCache_.invalidate();
        if (energyBackend_ == GCMCEnergyBackend::PgpFull ||
            energyBackend_ == GCMCEnergyBackend::PgpFullPme) {
            // Background changed; force a rebuild on the next insertion/background evaluation.
            pgpFullGridReady_ = false;
            pgpFullGridExcludedResidue_ = -999;
        }
    } else {
        // Restore the deleted instance at the same position (keeps same ID)
        bool restored = reservoir_->restoreInstance(instanceId, savedPosition, savedOrientation);
        if (restored) {
            // Successfully restored with same instance ID
            synchronizeStateWithReservoir(instanceId, true);
            result.residueIndex = instanceId;  // Keep original index
            relaxDrudeIfEnabled();
        } else {
            // Fallback: create a new instance if restoration failed
            // This can happen if the instance was already purged
            int restoredId = reservoir_->createInstance(typeId, savedPosition, savedOrientation);
            if (restoredId >= 0) {
                synchronizeStateWithReservoir(restoredId, true);
                result.residueIndex = restoredId;
                relaxDrudeIfEnabled();
            } else {
                std::cerr << "WARNING: Failed to restore deleted instance in GCMC deletion rejection" << std::endl;
            }
        }
        result.accepted = false;
    }
    
    totalMoves_++;
    
    // Sample statistics if configured using a lightweight countdown
    if (collectStats_) {
        if (--statsCountdown_ <= 0) {
            int particleCount = reservoir_ ? reservoir_->getActiveCount() : 0;
            double energy = calculateSystemEnergy();
            statistics_.addSample(totalMoves_, particleCount, energy,
                                  getAcceptanceRate(), temperature_, 100.0);
            statsCountdown_ = std::max(1, statsInterval_);
        }
    }
    
    return result;
}

// Attempt translation
GCMCEngine::MoveResult GCMCEngine::attemptTranslation(int residueIdx) {
    MoveResult result;
    result.type = MoveResult::TRANSLATE;
    result.residueIndex = residueIdx;

    // Diagnostics: non-CBMC move; prevent leaking previous trial energies.
    lastCbmcTrialEnergies_.clear();
    
    if (!state_ || !reservoir_) {
        result.accepted = false;
        return result;
    }
    
    FragmentInstance* instance = reservoir_->getInstance(residueIdx);
    if (!instance || !instance->isActive) {
        result.accepted = false;
        return result;
    }
    
    result.fragmentType = instance->templateId;
    
    // Store old position
    Vector3 oldPos = instance->position;
    result.position = oldPos;
    
    // Calculate energy before move
    if (useDrude_ ||
        energyBackend_ == GCMCEnergyBackend::Ewald ||
        energyBackend_ == GCMCEnergyBackend::Pme) {
        result.energyBefore = calculateSystemEnergy();
    } else {
        result.energyBefore = calculateFragmentEnergy(residueIdx);
    }
    
    // Generate translation using configured step size
    Vector3 displacement = generateTranslationVector(maxTranslationStep_);
    Vector3 newPos = oldPos + displacement;
    applyPeriodicBoundary(newPos);

    // Enforce region constraint: reject moves that leave region (atom-level)
    if (regionConstraint_) {
        // Use current orientation of the instance
        Quaternion orient = instance->orientation;
        if (!isMoleculeWithinRegion(result.fragmentType, newPos, orient)) {
            // Reject without changing position
            result.deltaE = 0.0;
            result.energyAfter = result.energyBefore;
            result.accepted = false;
            // Probability bookkeeping (only if configured)
            result.acceptanceProbability = shouldStoreProbability() ? 0.0 : -1.0;

            totalMoves_++;
            // Sample statistics if configured using a lightweight countdown
            if (collectStats_) {
                if (--statsCountdown_ <= 0) {
                    int particleCount = reservoir_ ? reservoir_->getActiveCount() : 0;
                    double energy = calculateSystemEnergy();
                    statistics_.addSample(totalMoves_, particleCount, energy,
                                          getAcceptanceRate(), temperature_, 100.0);
                    statsCountdown_ = std::max(1, statsInterval_);
                }
            }
            return result;
        }
    }

    // Update position
    updateFragmentPosition(residueIdx, newPos);
    
    // Calculate energy after move
    if (useDrude_ ||
        energyBackend_ == GCMCEnergyBackend::Ewald ||
        energyBackend_ == GCMCEnergyBackend::Pme) {
        result.energyAfter = calculateSystemEnergy();
    } else {
        result.energyAfter = calculateFragmentEnergy(residueIdx);
    }
    result.deltaE = result.energyAfter - result.energyBefore;
    
    // Accept or reject using unified acceptance calculator
    bool accept = false;
    double prob = 0.0;
    if (acceptanceCalculator_) {
        prob = acceptanceCalculator_->calculateTranslationProbability(result.deltaE, 1.0);
        
        // Only store probability if configured
        result.acceptanceProbability = shouldStoreProbability() ? prob : -1.0;
        
        accept = acceptanceCalculator_->acceptMove(prob);
    } else {
        // Metropolis criterion
        double beta = 1.0 / (8.314e-3 * temperature_);
        prob = std::min(1.0, std::exp(-beta * result.deltaE));
        // Apply same probability storage logic as main branch
        result.acceptanceProbability = shouldStoreProbability() ? prob : -1.0;
        accept = acceptMove(result.deltaE, 1.0, temperature_);
    }
    
    if (accept) {
        result.accepted = true;
        result.position = newPos;
        acceptedMoves_++;
        energyCache_.invalidate();
    } else {
        // Restore old position
        updateFragmentPosition(residueIdx, oldPos);
        relaxDrudeIfEnabled();
        result.accepted = false;
    }
    
    totalMoves_++;
    
    // Sample statistics if configured using a lightweight countdown
    if (collectStats_) {
        if (--statsCountdown_ <= 0) {
            int particleCount = reservoir_ ? reservoir_->getActiveCount() : 0;
            double energy = calculateSystemEnergy();
            statistics_.addSample(totalMoves_, particleCount, energy,
                                  getAcceptanceRate(), temperature_, 100.0);
            statsCountdown_ = std::max(1, statsInterval_);
        }
    }
    
    return result;
}

// Attempt rotation
GCMCEngine::MoveResult GCMCEngine::attemptRotation(int residueIdx) {
    MoveResult result;
    result.type = MoveResult::ROTATE;
    result.residueIndex = residueIdx;

    // Diagnostics: non-CBMC move; prevent leaking previous trial energies.
    lastCbmcTrialEnergies_.clear();
    
    if (!state_ || !reservoir_) {
        result.accepted = false;
        return result;
    }
    
    FragmentInstance* instance = reservoir_->getInstance(residueIdx);
    if (!instance || !instance->isActive) {
        result.accepted = false;
        return result;
    }
    
    result.fragmentType = instance->templateId;
    result.position = instance->position;
    
    // Store old orientation
    Quaternion oldOrient = instance->orientation;
    
    // Calculate energy before rotation
    if (useDrude_ ||
        energyBackend_ == GCMCEnergyBackend::Ewald ||
        energyBackend_ == GCMCEnergyBackend::Pme) {
        result.energyBefore = calculateSystemEnergy();
    } else {
        result.energyBefore = calculateFragmentEnergy(residueIdx);
    }
    
    // Generate rotation using configured angle
    Quaternion rotation = generateRotationQuaternion(maxRotationAngleRad_);
    // CRITICAL FIX: Actually apply the rotation by quaternion multiplication
    Quaternion newOrient = oldOrient * rotation;
    newOrient.normalize();
    
    // Enforce region constraint: reject rotations that push atoms outside region
    if (regionConstraint_) {
        if (!isMoleculeWithinRegion(result.fragmentType, result.position, newOrient)) {
            // Reject early without state updates
            result.deltaE = 0.0;
            result.energyAfter = result.energyBefore;
            result.accepted = false;
            result.acceptanceProbability = shouldStoreProbability() ? 0.0 : -1.0;
            totalMoves_++;
            if (collectStats_) {
                if (--statsCountdown_ <= 0) {
                    int particleCount = reservoir_ ? reservoir_->getActiveCount() : 0;
                    double energy = calculateSystemEnergy();
                    statistics_.addSample(totalMoves_, particleCount, energy,
                                          getAcceptanceRate(), temperature_, 100.0);
                    statsCountdown_ = std::max(1, statsInterval_);
                }
            }
            return result;
        }
    }

    // Update orientation
    updateFragmentOrientation(residueIdx, newOrient);
    
    // Calculate energy after rotation
    if (useDrude_ ||
        energyBackend_ == GCMCEnergyBackend::Ewald ||
        energyBackend_ == GCMCEnergyBackend::Pme) {
        result.energyAfter = calculateSystemEnergy();
    } else {
        result.energyAfter = calculateFragmentEnergy(residueIdx);
    }
    result.deltaE = result.energyAfter - result.energyBefore;
    
    // Use unified acceptance calculation for rotation
    bool accept = false;
    double prob = 0.0;
    if (acceptanceCalculator_) {
        // Rotation uses standard Metropolis criterion (no N dependence)
        prob = acceptanceCalculator_->calculateTranslationProbability(
            result.deltaE, 1.0);  // bias = 1.0 for rotation
        
        // Only store probability if configured
        result.acceptanceProbability = shouldStoreProbability() ? prob : -1.0;
        
        accept = acceptanceCalculator_->acceptMove(prob);
    } else {
        // Fallback to direct calculation
        double beta = 1.0 / (8.314e-3 * temperature_);
        prob = (result.deltaE <= 0) ? 1.0 : std::exp(-beta * result.deltaE);
        // Apply same probability storage logic as main branch
        result.acceptanceProbability = shouldStoreProbability() ? prob : -1.0;
        accept = acceptMove(result.deltaE, 1.0, temperature_);
    }
    
    if (accept) {
        result.accepted = true;
        acceptedMoves_++;
        energyCache_.invalidate();
    } else {
        // Restore old orientation
        updateFragmentOrientation(residueIdx, oldOrient);
        relaxDrudeIfEnabled();
        result.accepted = false;
    }
    
    totalMoves_++;
    
    // Sample statistics if configured using a lightweight countdown
    if (collectStats_) {
        if (--statsCountdown_ <= 0) {
            int particleCount = reservoir_ ? reservoir_->getActiveCount() : 0;
            double energy = calculateSystemEnergy();
            statistics_.addSample(totalMoves_, particleCount, energy,
                                  getAcceptanceRate(), temperature_, 100.0);
            statsCountdown_ = std::max(1, statsInterval_);
        }
    }
    
    return result;
}

// Attempt swap
GCMCEngine::MoveResult GCMCEngine::attemptSwap(int typeId1, int typeId2) {
    MoveResult result;
    result.type = MoveResult::SWAP;
    
    if (!state_ || !reservoir_) {
        result.accepted = false;
        return result;
    }
    
    // Select instances to swap
    int idx1 = selectRandomInstance(typeId1);
    int idx2 = selectRandomInstance(typeId2);
    
    if (idx1 < 0 || idx2 < 0) {
        result.accepted = false;
        return result;
    }
    
    // For simplicity, swap is delete type1 + insert type2
    // In practice, would swap identities directly
    
    result.accepted = false;  // Not fully implemented
    totalMoves_++;
    return result;
}

// Attempt regrowth
GCMCEngine::MoveResult GCMCEngine::attemptRegrowth(int residueIdx) {
    MoveResult result;
    result.type = MoveResult::REGROWTH;
    result.residueIndex = residueIdx;
    
    if (!state_ || !reservoir_) {
        result.accepted = false;
        return result;
    }
    
    // Regrowth = deletion + insertion at new position
    // Simplified implementation
    
    result.accepted = false;  // Not fully implemented
    totalMoves_++;
    return result;
}

// Attempt cluster move
GCMCEngine::MoveResult GCMCEngine::attemptClusterMove(int residueIdx, double cutoff) {
    MoveResult result;
    result.type = MoveResult::CLUSTER;
    result.residueIndex = residueIdx;
    
    if (!state_ || !reservoir_) {
        result.accepted = false;
        return result;
    }
    
    // Find cluster
    std::vector<int> cluster = selectCluster(residueIdx, cutoff);
    
    if (cluster.empty()) {
        result.accepted = false;
        return result;
    }
    
    // Move entire cluster together
    // Simplified - not fully implemented
    
    result.accepted = false;
    totalMoves_++;
    return result;
}

// Select random fragment type
int GCMCEngine::selectRandomFragment() {
    if (!reservoir_) return -1;
    
    int nTypes = reservoir_->getTemplateCount();
    if (nTypes == 0) return -1;
    
    std::uniform_int_distribution<int> dist(0, nTypes - 1);
    return dist(rng_);
}

// Select random instance
int GCMCEngine::selectRandomInstance(int typeId) {
    if (!reservoir_) return -1;

    std::vector<int> instances = reservoir_->getActiveInstances(typeId);
    if (instances.empty()) return -1;

    // If region constraint is set, filter instances to those within the region
    if (regionConstraint_) {
        std::vector<int> filteredInstances;
        for (int id : instances) {
            FragmentInstance* inst = reservoir_->getInstance(id);
            if (inst) {
                movement::Vector3 pos(inst->position.x, inst->position.y, inst->position.z);
                if (regionConstraint_->isInRegion(pos)) {
                    filteredInstances.push_back(id);
                }
            }
        }

        // Use filtered list if not empty; otherwise (optionally) fall back to all instances.
        // NOTE: Falling back violates strict region detailed balance when the target
        // distribution is defined only over the constrained region.
        if (!filteredInstances.empty()) {
            instances = filteredInstances;
        } else if (getConfigValue("strict_region_balance") > 0.5) {
            return -1;
        }
    }

    std::uniform_int_distribution<int> dist(0, instances.size() - 1);
    return instances[dist(rng_)];
}

// Select cluster
std::vector<int> GCMCEngine::selectCluster(int seedIdx, double cutoff) {
    std::vector<int> cluster;
    if (!reservoir_) return cluster;
    
    FragmentInstance* seed = reservoir_->getInstance(seedIdx);
    if (!seed) return cluster;
    
    cluster.push_back(seedIdx);
    
    // Find neighbors within cutoff
    std::vector<int> allInstances = reservoir_->getActiveInstances();
    for (int idx : allInstances) {
        if (idx == seedIdx) continue;
        
        FragmentInstance* instance = reservoir_->getInstance(idx);
        if (!instance) continue;
        
        double dist = minimumImageDistance(seed->position, instance->position);
        if (dist < cutoff) {
            cluster.push_back(idx);
        }
    }
    
    return cluster;
}

// Generate random position
Vector3 GCMCEngine::generateRandomPosition() {
    if (!state_) {
        return Vector3(uniform_(rng_) * 100,
                      uniform_(rng_) * 100,
                      uniform_(rng_) * 100);
    }
    
    // Get box dimensions with fallback
    double boxX, boxY, boxZ;
    if (state_->periodicBox.size() >= 3) {
        boxX = state_->periodicBox[0];
        boxY = state_->periodicBox[1];
        boxZ = state_->periodicBox[2];
    } else if (state_->info.box[0] > 0 && state_->info.box[1] > 0 && state_->info.box[2] > 0) {
        // Fallback to info.box if periodicBox not set
        boxX = state_->info.box[0];
        boxY = state_->info.box[1];
        boxZ = state_->info.box[2];
    } else {
        // No valid box dimensions, use default
        boxX = boxY = boxZ = 100.0;
    }
    
    // If region constraint is set, sample from constrained region
    if (regionConstraint_) {
        movement::Vector3 movPos = regionConstraint_->samplePosition();
        return Vector3(movPos.x, movPos.y, movPos.z);
    }

    // Use [0, L) coordinate system to match CavityManager
    return Vector3(uniform_(rng_) * boxX,
                  uniform_(rng_) * boxY,
                  uniform_(rng_) * boxZ);
}

// Generate cavity position
Vector3 GCMCEngine::generateCavityPosition() {
    if (!cavityManager_) {
        return generateRandomPosition();
    }

    // Try to find a cavity position within the region constraint
    const int maxAttempts = 100;
    for (int attempt = 0; attempt < maxAttempts; ++attempt) {
        // Get cavity from manager (returns movement::Vector3)
        movement::Vector3 movCavity = cavityManager_->selectCavity();

        // Convert to montecarlo::Vector3
        Vector3 cavity(movCavity.x, movCavity.y, movCavity.z);

        // Add small random displacement
        Vector3 displacement(normal_(rng_) * 0.5,
                            normal_(rng_) * 0.5,
                            normal_(rng_) * 0.5);

        Vector3 position = cavity + displacement;

        // Check if position is within region constraint
        if (!regionConstraint_ || regionConstraint_->isInRegion(movement::Vector3(position.x, position.y, position.z))) {
            return position;
        }
    }

    // If no valid cavity found in region, fall back to random position in region
    if (regionConstraint_) {
        movement::Vector3 movPos = regionConstraint_->samplePosition();
        return Vector3(movPos.x, movPos.y, movPos.z);
    }

    // Last resort: random position in box
    return generateRandomPosition();
}

// Generate random orientation
Quaternion GCMCEngine::generateRandomOrientation() {
    // Use the correct Shoemake algorithm for uniform quaternion distribution
    // This matches gcmc_gpu's create_random_quarternion() implementation
    double u = uniform_(rng_);
    double v = uniform_(rng_);
    double w = uniform_(rng_);
    
    double sqrt_1_minus_u = std::sqrt(1.0 - u);
    double sqrt_u = std::sqrt(u);
    double two_pi_v = 2.0 * M_PI * v;
    double two_pi_w = 2.0 * M_PI * w;
    
    Quaternion q(
        sqrt_u * std::cos(two_pi_w),           // w component
        sqrt_1_minus_u * std::sin(two_pi_v),   // x component
        sqrt_1_minus_u * std::cos(two_pi_v),   // y component
        sqrt_u * std::sin(two_pi_w)            // z component
    );
    q.normalize();
    
    return q;
}

// Generate translation vector
Vector3 GCMCEngine::generateTranslationVector(double maxDist) {
    // Random direction
    double theta = uniform_(rng_) * 2 * M_PI;
    double phi = std::acos(2 * uniform_(rng_) - 1);
    
    // Random magnitude
    double r = uniform_(rng_) * maxDist;
    
    return Vector3(r * std::sin(phi) * std::cos(theta),
                  r * std::sin(phi) * std::sin(theta),
                  r * std::cos(phi));
}

// Generate rotation quaternion
Quaternion GCMCEngine::generateRotationQuaternion(double maxAngle) {
    // Random axis
    Vector3 axis(normal_(rng_), normal_(rng_), normal_(rng_));
    double norm = axis.norm();
    if (norm > 0) {
        axis = axis * (1.0 / norm);
    }
    
    // Random angle
    double angle = uniform_(rng_) * maxAngle;
    
    // Create quaternion from axis-angle
    double halfAngle = angle / 2;
    double s = std::sin(halfAngle);
    
    Quaternion q(std::cos(halfAngle), s * axis.x, s * axis.y, s * axis.z);
    q.normalize();
    
    return q;
}

void GCMCEngine::relaxDrudeIfEnabled() {
    if (!useDrude_ || !state_) {
        return;
    }
    rebuildDrudeTopologyIfNeeded();
    ::pygcmc::platform::cpu::DrudeComplete::calculateEnergy(*state_);
}

void GCMCEngine::rebuildDrudeTopologyIfNeeded() {
    if (!useDrude_ || !state_) {
        return;
    }
    if (!drudeTopologyDirty_) {
        return;
    }

    ::pygcmc::platform::cpu::DrudeComplete::clear();
    ::pygcmc::platform::cpu::DrudeComplete::clearHistory();

    auto& core = ::pygcmc::platform::cpu::DrudeComplete::getDrudeCore();
    for (const auto& particle : drudeHostParticles_) {
        core.addParticle(particle);
    }
    for (const auto& pair : drudeHostScreenedPairs_) {
        ::pygcmc::platform::cpu::DrudeComplete::addScreenedPair(pair);
    }

    appendDrudeForActiveFragments();

    drudeTopologyDirty_ = false;
}

void GCMCEngine::appendDrudeForActiveFragments() {
    if (!reservoir_ || !state_) {
        return;
    }
    for (int instanceId : reservoir_->getActiveInstances()) {
        appendDrudeForFragmentInstance(instanceId);
    }
}

void GCMCEngine::appendDrudeForFragmentInstance(int instanceId) {
    if (!reservoir_ || !state_) {
        return;
    }
    if (instanceId < 0 || instanceId >= static_cast<int>(state_->residues.size())) {
        return;
    }

    const auto& residue = state_->residues[instanceId];
    if (!residue.active || residue.atomCount <= 0) {
        return;
    }

    const FragmentInstance* instance = reservoir_->getInstance(instanceId);
    if (!instance) {
        return;
    }
    const FragmentTemplate* tmpl = reservoir_->getTemplate(instance->templateId);
    if (!tmpl) {
        return;
    }
    const int nAtoms = static_cast<int>(tmpl->atoms.size());
    if (nAtoms <= 0 || residue.atomCount != nAtoms) {
        return;
    }
    if (tmpl->bonds.empty()) {
        return;
    }

    auto isDrudeAtom = [&](int localIdx) -> bool {
        if (localIdx < 0 || localIdx >= nAtoms) {
            return false;
        }
        const auto& a = tmpl->atoms[static_cast<size_t>(localIdx)];
        if (std::abs(static_cast<double>(a.mass) - ::pygcmc::platform::cpu::DrudeConstants::DRUDE_MASS) < 1e-3) {
            return true;
        }
        if (static_cast<size_t>(localIdx) < tmpl->atomTypeNames.size()) {
            std::string t = tmpl->atomTypeNames[static_cast<size_t>(localIdx)];
            for (auto& c : t) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            if (t == "DRUD" || t == "DRUDE") {
                return true;
            }
        }
        return false;
    };

    std::vector<std::vector<int>> adjacency(static_cast<size_t>(nAtoms));
    adjacency.reserve(static_cast<size_t>(nAtoms));
    for (const auto& b : tmpl->bonds) {
        if (b.atom1 < 0 || b.atom2 < 0) {
            continue;
        }
        if (b.atom1 >= nAtoms || b.atom2 >= nAtoms) {
            continue;
        }
        adjacency[static_cast<size_t>(b.atom1)].push_back(b.atom2);
        adjacency[static_cast<size_t>(b.atom2)].push_back(b.atom1);
    }

    std::vector<int> parentLocalToDipole(static_cast<size_t>(nAtoms), -1);
    std::vector<double> parentThole(static_cast<size_t>(nAtoms), 0.0);

    for (int drudeLocal = 0; drudeLocal < nAtoms; ++drudeLocal) {
        if (!isDrudeAtom(drudeLocal)) {
            continue;
        }

        int parentLocal = -1;
        for (int neighbor : adjacency[static_cast<size_t>(drudeLocal)]) {
            if (neighbor < 0 || neighbor >= nAtoms) {
                continue;
            }
            if (!isDrudeAtom(neighbor)) {
                parentLocal = neighbor;
                break;
            }
        }
        if (parentLocal < 0) {
            continue;
        }

        if (!drudeForceFieldModel_) {
            continue;
        }
        if (static_cast<size_t>(parentLocal) >= tmpl->atomTypeNames.size()) {
            continue;
        }

        const std::string& parentTypeName = tmpl->atomTypeNames[static_cast<size_t>(parentLocal)];
        if (!drudeForceFieldModel_->has_alpha_params(parentTypeName)) {
            continue;
        }
        const auto& alphaParams = drudeForceFieldModel_->get_alpha_params(parentTypeName);
        const double alphaNm3 = alphaParams.alpha * 1e-3;
        if (!(alphaNm3 > 0.0)) {
            continue;
        }

        const int parentIdx = residue.atomStart + parentLocal;
        const int drudeIdx = residue.atomStart + drudeLocal;
        if (parentIdx < 0 || drudeIdx < 0) {
            continue;
        }
        if (parentIdx >= state_->activeAtomCount || drudeIdx >= state_->activeAtomCount) {
            continue;
        }

        ::pygcmc::platform::cpu::DrudeParticle particle;
        particle.parentIndex = parentIdx;
        particle.drudeIndex = drudeIdx;
        particle.charge = static_cast<double>(state_->atoms[static_cast<size_t>(drudeIdx)].charge);
        particle.polarizability = alphaNm3;
        particle.computeSpringConstants();

        const int dipoleIdx = ::pygcmc::platform::cpu::DrudeComplete::getDrudeCore().addParticle(particle);
        parentLocalToDipole[static_cast<size_t>(parentLocal)] = dipoleIdx;
        parentThole[static_cast<size_t>(parentLocal)] = alphaParams.thole;
    }

    auto makeKey = [](int a, int b) -> uint64_t {
        const uint32_t x = static_cast<uint32_t>(std::min(a, b));
        const uint32_t y = static_cast<uint32_t>(std::max(a, b));
        return (static_cast<uint64_t>(x) << 32) | static_cast<uint64_t>(y);
    };

    std::unordered_set<uint64_t> screenedParentPairs;
    screenedParentPairs.reserve(static_cast<size_t>(nAtoms));

    // 1-2 (bonded) polarizable parents.
    for (const auto& b : tmpl->bonds) {
        if (b.atom1 < 0 || b.atom2 < 0) {
            continue;
        }
        if (b.atom1 >= nAtoms || b.atom2 >= nAtoms) {
            continue;
        }
        if (parentLocalToDipole[static_cast<size_t>(b.atom1)] >= 0 &&
            parentLocalToDipole[static_cast<size_t>(b.atom2)] >= 0) {
            screenedParentPairs.insert(makeKey(b.atom1, b.atom2));
        }
    }

    // 1-3 (topological distance 2) polarizable parents.
    for (size_t mid = 0; mid < adjacency.size(); ++mid) {
        const auto& neighbors = adjacency[mid];
        for (size_t ia = 0; ia < neighbors.size(); ++ia) {
            for (size_t ib = ia + 1; ib < neighbors.size(); ++ib) {
                const int a = neighbors[ia];
                const int b = neighbors[ib];
                if (a < 0 || b < 0) {
                    continue;
                }
                if (a >= nAtoms || b >= nAtoms) {
                    continue;
                }
                if (parentLocalToDipole[static_cast<size_t>(a)] >= 0 &&
                    parentLocalToDipole[static_cast<size_t>(b)] >= 0) {
                    screenedParentPairs.insert(makeKey(a, b));
                }
            }
        }
    }

    // Register screened pairs with per-pair NBTHOLE override when available.
    for (const uint64_t key : screenedParentPairs) {
        const int parentA = static_cast<int>(key >> 32);
        const int parentB = static_cast<int>(key & 0xffffffffu);
        if (parentA < 0 || parentB < 0 || parentA >= nAtoms || parentB >= nAtoms) {
            continue;
        }
        const int dipoleA = parentLocalToDipole[static_cast<size_t>(parentA)];
        const int dipoleB = parentLocalToDipole[static_cast<size_t>(parentB)];
        if (dipoleA < 0 || dipoleB < 0 || dipoleA == dipoleB) {
            continue;
        }

        double thole = std::abs(parentThole[static_cast<size_t>(parentA)]) +
                       std::abs(parentThole[static_cast<size_t>(parentB)]);

        if (drudeForceFieldModel_ &&
            static_cast<size_t>(parentA) < tmpl->atomTypeNames.size() &&
            static_cast<size_t>(parentB) < tmpl->atomTypeNames.size()) {
            const auto [nbthole, hasNbthole] = drudeForceFieldModel_->get_nbthole(
                tmpl->atomTypeNames[static_cast<size_t>(parentA)],
                tmpl->atomTypeNames[static_cast<size_t>(parentB)]);
            if (hasNbthole) {
                thole = std::abs(nbthole);
            }
        }

        ::pygcmc::platform::cpu::ScreenedPair pair;
        pair.dipole1 = dipoleA;
        pair.dipole2 = dipoleB;
        pair.thole = thole;
        ::pygcmc::platform::cpu::DrudeComplete::addScreenedPair(pair);
    }
}

// Calculate system energy
double GCMCEngine::calculateSystemEnergy() {
    if (!state_) return 0.0;

    // Full-system long-range backends (EWALD/PME) may evaluate multiple trial states
    // within a single MC move; caching is unsafe there unless we implement fine-grained
    // invalidation. Keep caching enabled for local backends only (DIRECT/PGP).
    if (energyCache_.valid &&
        energyBackend_ != GCMCEnergyBackend::Ewald &&
        energyBackend_ != GCMCEnergyBackend::Pme) {
        return energyCache_.totalEnergy;
    }

    // Relax Drude oscillators (Born–Oppenheimer surface) before evaluating nonbonded energy.
    // This mutates Drude particle coordinates in-place.
    double drudeEnergy = 0.0;
    if (useDrude_) {
        rebuildDrudeTopologyIfNeeded();
        drudeEnergy = ::pygcmc::platform::cpu::DrudeComplete::calculateEnergy(*state_);
    }

    double totalEnergy = 0.0;

    // Use energy callback if available
    if (energyCallback_) {
        totalEnergy = energyCallback_->calculateSystemEnergy(*state_);
    } else {
        // Fallback to direct energy module usage
        if (energyMethod_ == EnergyMethod::PME) {
            computeSystemEnergy(*state_, EnergyMethod::PME, true, true);
        } else if (energyMethod_ == EnergyMethod::EWALD) {
            computeSystemEnergy(*state_, EnergyMethod::EWALD, true, true);
        } else {
            // DIRECT with cutoff and PBC
            computeSystemEnergy(*state_, EnergyMethod::DIRECT, true, true);
        }
        
        // Get total energy from state
        totalEnergy = energy::getTotalEnergyUniquePairs(*state_, energyMethod_);
    }

    totalEnergy += drudeEnergy;
    
    energyCache_.totalEnergy = totalEnergy;
    energyCache_.valid = true;
    
    return totalEnergy;
}

// Calculate fragment energy
double GCMCEngine::calculateFragmentEnergy(int residueIdx) {
    if (!state_ || residueIdx < 0 || residueIdx >= static_cast<int>(state_->residues.size())) {
        return 0.0;
    }
    
    // Get residue from state
    auto& residue = state_->residues[residueIdx];
    if (!residue.active) return 0.0;

    // PGP backends bypass the generic EnergyMethod dispatch and compute per-residue energies
    // directly against the current background/grid configuration.
    if (energyBackend_ == GCMCEnergyBackend::PgpHost) {
        return calculateFragmentEnergyPgpHost(residueIdx);
    }
    if (energyBackend_ == GCMCEnergyBackend::PgpFull ||
        energyBackend_ == GCMCEnergyBackend::PgpFullPme) {
        // For translation/rotation/deletion, the background grid must exclude the moved residue.
        // (CBMC insertion/deletion retracing handle grid reuse explicitly and call the
        // calculateFragmentEnergyPgpFullUsingCurrentGrid helper directly.)
        ensurePgpFullGridReadyExcluding(residueIdx);
        return calculateFragmentEnergyPgpFullUsingCurrentGrid(residueIdx);
    }
    
    // Check cache
    if (reservoir_) {
        FragmentInstance* instance = reservoir_->getInstance(residueIdx);
        if (instance && instance->lastEnergyUpdate == totalMoves_) {
            return instance->energy_total;
        }
    }
    
    double energy = 0.0;
    
    // Use energy callback if available
    if (energyCallback_) {
        energy = energyCallback_->calculateResidueEnergy(*state_, residueIdx);
    } else {
        // Fallback to direct energy module usage
        // Mark residue as moved for movement energy calculation
        // residue.moved /* moved flag not in MCResidue */ = true;
        
        // Calculate movement energy for this residue
        if (energyMethod_ == EnergyMethod::PME) {
            computeMovementEnergy(*state_, EnergyMethod::PME, true, true);
        } else if (energyMethod_ == EnergyMethod::EWALD) {
            computeMovementEnergy(*state_, EnergyMethod::EWALD, true, true);
        } else {
            computeMovementEnergy(*state_, EnergyMethod::DIRECT, true, true);
        }
        
        // Reset moved flag
        // residue.moved /* moved flag not in MCResidue */ = false;
        
        energy = residue.energy_vdw + residue.energy_elec;
    }
    
    // Update cache
    if (reservoir_) {
        FragmentInstance* instance = reservoir_->getInstance(residueIdx);
        if (instance) {
            instance->energy_total = energy;
            instance->lastEnergyUpdate = totalMoves_;
        }
    }
    
    return energy;
}

// Calculate interaction energy
double GCMCEngine::calculateInteractionEnergy(int residueIdx) {
    // This function is now deprecated - use calculateFragmentEnergy instead
    // which properly uses the energy module
    return calculateFragmentEnergy(residueIdx);
}

// Calculate pair energy
double GCMCEngine::calculatePairEnergy(int idx1, int idx2) {
    // Suppress unused parameter warnings
    (void)idx1;
    (void)idx2;
    
    // This function is now deprecated - the energy module handles
    // all pair interactions properly with the correct force field parameters
    // Use calculateFragmentEnergy or calculateSystemEnergy instead
    return 0.0;
}

// Calculate insertion bias
double GCMCEngine::calculateInsertionBias(const FragmentTemplate& tmpl,
                                         const Vector3& position,
                                         const Quaternion& orientation) {
    // Suppress unused parameter warnings
    (void)tmpl;
    (void)position;
    (void)orientation;

    double bias = 1.0;

    if (configBias_) {
        // Config bias calculation would go here
        bias *= 1.0;
    }

    // Proposal bias for detailed balance when using target_numwaters
    // For insertion: multiply by p_delete/p_insert ratio
    double proposalBias = getConfigValue("proposalBias");
    if (proposalBias > 0) {
        bias *= proposalBias;  // This is p_delete/p_insert for insertion
    }

    return bias;
}

// Calculate deletion bias at specific position (more robust)
double GCMCEngine::calculateDeletionBiasAtPosition(const Vector3& position) {
    double bias = 1.0;

    // CRITICAL: For detailed balance, deletion bias must match insertion bias
    // at the same position
    (void)position;

    if (configBias_) {
        // Config bias calculation would go here (must match insertion)
        bias *= 1.0;
    }

	// Proposal bias for detailed balance when using target_numwaters
	// For deletion: multiply by p_insert/p_delete ratio
	double proposalBias = getConfigValue("proposalBias");
	if (proposalBias > 0) {
	    bias *= proposalBias;
	}

	return bias;
}

// Calculate deletion bias (legacy - depends on reservoir state)
double GCMCEngine::calculateDeletionBias(int residueIdx) {
    // Try to get position from reservoir
    FragmentInstance* instance = reservoir_->getInstance(residueIdx);
    if (instance) {
        Vector3 pos(instance->position.x, instance->position.y, instance->position.z);
        return calculateDeletionBiasAtPosition(pos);
    }
    // Fallback if instance not accessible
    return 1.0;
}

// Calculate regrowth bias
double GCMCEngine::calculateRegrowthBias(int residueIdx) {
    // Combination of deletion and insertion biases
    double deletionBias = calculateDeletionBias(residueIdx);
    
    // Would calculate insertion bias at new position
    double insertionBias = 1.0;
    
    return deletionBias * insertionBias;
}

double GCMCEngine::getBoxVolume() const {
    if (!state_) {
        return 0.0;
    }
    double ax = static_cast<double>(state_->info.box[0]);
    double by = static_cast<double>(state_->info.box[1]);
    double cz = static_cast<double>(state_->info.box[2]);
    if (ax > 0.0 && by > 0.0 && cz > 0.0) {
        return ax * by * cz;
    }
    double fallback = static_cast<double>(state_->info.volume);
    if (fallback > 0.0) {
        return fallback;
    }
    return 0.0;
}

// Accept move
bool GCMCEngine::acceptMove(double deltaE, double bias, double temperature) {
    if (deltaE <= 0) return true;
    
    double beta = 1.0 / (8.314e-3 * temperature);  // kJ/(mol*K)
    double probability = bias * std::exp(-beta * deltaE);
    
    return uniform_(rng_) < probability;
}

// Calculate acceptance probability
double GCMCEngine::calculateAcceptanceProbability(const MoveResult& result,
                                                 double temperature) {
    double beta = 1.0 / (8.314e-3 * temperature);
    return std::min(1.0, result.bias * std::exp(-beta * result.deltaE));
}

// Synchronize MCState with reservoir
void GCMCEngine::synchronizeStateWithReservoir(int instanceId, bool isInsertion) {
    if (!state_ || !reservoir_) return;
    
    FragmentInstance* instance = reservoir_->getInstance(instanceId);
    if (!instance) return;
    
    if (isInsertion) {
        // Add residue to MCState
        if (instanceId >= static_cast<int>(state_->residues.size())) {
            // Need to expand residues vector
            state_->residues.resize(instanceId + 1);
        }
        
        // Get template for atom information
        const FragmentTemplate* tmpl = reservoir_->getTemplate(instance->templateId);
        if (!tmpl) return;
        
        // Update residue in state
        auto& residue = state_->residues[instanceId];
        residue.active = true;
        residue.resid = instanceId;
        residue.resname = tmpl->name;
        // MCResidue doesn't have chainid field
        // residue.moved /* moved flag not in MCResidue */ = false;
        residue.energy_vdw = 0.0;
        residue.energy_elec = 0.0;
        
        // CRITICAL: Set atomStart and atomCount for energy calculations
        residue.atomStart = state_->activeAtomCount;
        residue.atomCount = static_cast<int>(tmpl->atoms.size());
        
        // Clear and add atoms
        residue.atoms.clear();
        residue.atoms.reserve(tmpl->atoms.size());
        
        // Transform template atoms by instance position and orientation
        for (const auto& tmplAtom : tmpl->atoms) {
            MCAtom atom;
            // MCAtom doesn't have active or atomId fields
            // atom.active = true;
            // atom.atomId = tmplAtom.id;  // tmplAtom may not have id field
            atom.type = tmplAtom.type;
            atom.charge = tmplAtom.charge;
            atom.mass = tmplAtom.mass;
            atom.name = tmplAtom.name;
            
            // Apply rotation and translation
            // Quaternion rotation of vector: q * v * q^-1
            // For unit quaternion, q^-1 = (w, -x, -y, -z)
            Vector3 v(tmplAtom.x, tmplAtom.y, tmplAtom.z);
            Quaternion q = instance->orientation;
            
            // Use the movement layer's Quaternion::rotate() method
            double vx = v.x, vy = v.y, vz = v.z;
            Vector3 templatePos(vx, vy, vz);
            Vector3 rotatedPos = q.rotate(templatePos);
            
            atom.x = instance->position.x + rotatedPos.x;
            atom.y = instance->position.y + rotatedPos.y;
            atom.z = instance->position.z + rotatedPos.z;
            
            // Sync position vector with x,y,z coordinates
            atom.updatePosition();
            
            residue.atoms.push_back(atom);
            
            // CRITICAL: Also add to global atoms array for energy calculations
            if (state_->activeAtomCount < static_cast<int>(state_->atoms.size())) {
                state_->atoms[state_->activeAtomCount] = atom;
            } else {
                state_->atoms.push_back(atom);
            }
            state_->activeAtomCount++;
        }
        
        // Update residue index in fragment instance
        instance->residueIndex = instanceId;
        
        // CRITICAL: Update activeResidueCount
        // Find the highest active residue index + 1
        state_->activeResidueCount = 0;
        for (int i = 0; i < static_cast<int>(state_->residues.size()); i++) {
            if (state_->residues[i].active) {
                state_->activeResidueCount = i + 1;
            }
        }
    } else {
        // Mark residue as inactive (deletion case)
        if (instanceId < static_cast<int>(state_->residues.size())) {
            auto& residue = state_->residues[instanceId];
            residue.active = false;
            
            // CRITICAL: Don't actually remove atoms from global array to avoid shifting indices
            // Just mark the residue as inactive so energy calculations skip it
            residue.atoms.clear();
            
            // CRITICAL: Update activeResidueCount
            // Find the highest active residue index + 1
            state_->activeResidueCount = 0;
            for (int i = 0; i < static_cast<int>(state_->residues.size()); i++) {
                if (state_->residues[i].active) {
                    state_->activeResidueCount = i + 1;
                }
            }
            
            // Note: We don't decrement activeAtomCount here to avoid index shifting
            // This is a simplification for now - a production system would compact arrays
        }
    }

    // State atom/residue arrays were mutated (insert/delete); cached energies are invalid.
    energyCache_.invalidate();
    if (useDrude_) {
        drudeTopologyDirty_ = true;
    }
}

// Update fragment position
void GCMCEngine::updateFragmentPosition(int residueIdx, const Vector3& newPos) {
    if (!reservoir_) return;
    
    reservoir_->updatePosition(residueIdx, newPos);
    if (auto* instance = reservoir_->getInstance(residueIdx)) {
        // Position change invalidates any per-instance cached energy for this MC move.
        instance->lastEnergyUpdate = -1.0;
    }
    energyCache_.invalidate();
    
    // Update atom coordinates without adding new atoms
    updateAtomCoordinates(residueIdx);
}

// Update fragment orientation
void GCMCEngine::updateFragmentOrientation(int residueIdx, const Quaternion& newOrient) {
    if (!reservoir_) return;
    
    reservoir_->updateOrientation(residueIdx, newOrient);
    if (auto* instance = reservoir_->getInstance(residueIdx)) {
        // Orientation change invalidates any per-instance cached energy for this MC move.
        instance->lastEnergyUpdate = -1.0;
    }
    energyCache_.invalidate();
    
    // Update atom coordinates without adding new atoms
    updateAtomCoordinates(residueIdx);
}

// Update atom coordinates for an existing residue without changing atom count
void GCMCEngine::updateAtomCoordinates(int residueIdx) {
    if (!state_ || !reservoir_) return;
    
    if (residueIdx >= static_cast<int>(state_->residues.size())) return;
    
    auto& residue = state_->residues[residueIdx];
    if (!residue.active) return;
    
    FragmentInstance* instance = reservoir_->getInstance(residueIdx);
    if (!instance) return;
    
    const FragmentTemplate* tmpl = reservoir_->getTemplate(instance->templateId);
    if (!tmpl) return;
    
    // Update atoms in both residue.atoms and state->atoms arrays
    int atomIdx = 0;
    for (const auto& tmplAtom : tmpl->atoms) {
        if (atomIdx >= residue.atomCount) break;
        
        // Apply rotation and translation
        Vector3 v(tmplAtom.x, tmplAtom.y, tmplAtom.z);
        Quaternion q = instance->orientation;
        
        // Use the verified Quaternion::rotate() method
        Vector3 rotatedPos = q.rotate(v);
        
        // Update coordinates in residue.atoms
        if (atomIdx < static_cast<int>(residue.atoms.size())) {
            residue.atoms[atomIdx].x = instance->position.x + rotatedPos.x;
            residue.atoms[atomIdx].y = instance->position.y + rotatedPos.y;
            residue.atoms[atomIdx].z = instance->position.z + rotatedPos.z;
            residue.atoms[atomIdx].updatePosition();
        }
        
        // Update coordinates in global atoms array
        int globalIdx = residue.atomStart + atomIdx;
        if (globalIdx < static_cast<int>(state_->atoms.size())) {
            state_->atoms[globalIdx].x = instance->position.x + rotatedPos.x;
            state_->atoms[globalIdx].y = instance->position.y + rotatedPos.y;
            state_->atoms[globalIdx].z = instance->position.z + rotatedPos.z;
            state_->atoms[globalIdx].updatePosition();
        }
        
        atomIdx++;
    }
}

// Apply periodic boundary conditions
void GCMCEngine::applyPeriodicBoundary(Vector3& position) {
    if (!state_) return;
    
    // Get box dimensions with fallback
    double boxX, boxY, boxZ;
    if (state_->periodicBox.size() >= 3) {
        boxX = state_->periodicBox[0];
        boxY = state_->periodicBox[1];
        boxZ = state_->periodicBox[2];
    } else if (state_->info.box[0] > 0 && state_->info.box[1] > 0 && state_->info.box[2] > 0) {
        // Fallback to info.box if periodicBox not set
        boxX = state_->info.box[0];
        boxY = state_->info.box[1];
        boxZ = state_->info.box[2];
    } else {
        // No valid box dimensions, skip PBC
        return;
    }
    
    // Use [0, L) coordinate system
    while (position.x < 0) position.x += boxX;
    while (position.x >= boxX) position.x -= boxX;
    while (position.y < 0) position.y += boxY;
    while (position.y >= boxY) position.y -= boxY;
    while (position.z < 0) position.z += boxZ;
    while (position.z >= boxZ) position.z -= boxZ;
}

// Check whether all atoms of a fragment (template) lie within the region after transform
bool GCMCEngine::isMoleculeWithinRegion(int typeId, const Vector3& position, const Quaternion& orientation) const {
    if (!regionConstraint_ || !reservoir_) return true;

    const FragmentTemplate* tmpl = reservoir_->getTemplate(typeId);
    if (!tmpl) return true;

    for (const auto& a : tmpl->atoms) {
        Vector3 v(a.x, a.y, a.z);
        Vector3 vr = orientation.rotate(v);
        movement::Vector3 world(position.x + vr.x, position.y + vr.y, position.z + vr.z);
        if (!regionConstraint_->isInRegion(world)) {
            return false;
        }
    }
    return true;
}

// Calculate minimum image distance
double GCMCEngine::minimumImageDistance(const Vector3& r1, const Vector3& r2) {
    if (!state_) {
        return (r1 - r2).norm();
    }
    
    Vector3 dr = r1 - r2;
    
    // Get box dimensions with fallback (same as applyPeriodicBoundary)
    double boxX, boxY, boxZ;
    if (state_->periodicBox.size() >= 3) {
        boxX = state_->periodicBox[0];
        boxY = state_->periodicBox[1];
        boxZ = state_->periodicBox[2];
    } else if (state_->info.box[0] > 0 && state_->info.box[1] > 0 && state_->info.box[2] > 0) {
        // Fallback to info.box if periodicBox not set
        boxX = state_->info.box[0];
        boxY = state_->info.box[1];
        boxZ = state_->info.box[2];
    } else {
        // No valid box dimensions, return direct distance
        return dr.norm();
    }
    
    // Apply minimum image convention for [0, L) coordinate system
    if (std::abs(dr.x) > boxX/2) {
        dr.x = dr.x - std::copysign(boxX, dr.x);
    }
    if (std::abs(dr.y) > boxY/2) {
        dr.y = dr.y - std::copysign(boxY, dr.y);
    }
    if (std::abs(dr.z) > boxZ/2) {
        dr.z = dr.z - std::copysign(boxZ, dr.z);
    }
    
    return dr.norm();
}

// Dynamic configuration implementation
void GCMCEngine::setConfigValue(const std::string& key, double value) {
    configMap_[key] = value;
    
    // Apply specific configuration changes
    if (key == "temperature") {
        temperature_ = value;
        if (acceptanceCalculator_) {
            acceptanceCalculator_->setTemperature(value);
        }
    } else if (key == "cutoff") {
        cutoff_ = value;
    } else if (key == "statsInterval") {
        statsInterval_ = static_cast<int>(value);
        statistics_.setSamplingInterval(statsInterval_);
    } else if (key == "autoAdjustStats") {
        statistics_.setAutoAdjust(value > 0.5);
    } else if (key == "collectStats") {
        collectStats_ = (value > 0.5);
    } else if (key == "maxTranslation") {
        maxTranslationStep_ = value;
    } else if (key == "maxRotation" || key == "maxRotationAngle") {
        maxRotationAngleRad_ = value;  // Support both key names
    } else if (key == "useCavityBias") {
        useCavityBias_ = (value > 0.5);
    } else if (key == "useConfBias") {
        useConfBias_ = (value > 0.5);
    } else if (key == "storeProbabilities") {
        // Clear cache when config changes
        storeProbabilityCached_ = false;
    } else if (key == "storeAcceptanceProbability") {
        // Clear cache when config changes
        storeProbabilityCached_ = false;
    }
}

double GCMCEngine::getConfigValue(const std::string& key) const {
    auto it = configMap_.find(key);
    if (it != configMap_.end()) {
        return it->second;
    }
    
    // Return current values for known keys
    if (key == "temperature") return temperature_;
    if (key == "cutoff") return cutoff_;
    if (key == "statsInterval") return static_cast<double>(statsInterval_);
    if (key == "collectStats") return collectStats_ ? 1.0 : 0.0;
    if (key == "maxTranslation") return maxTranslationStep_;
    if (key == "maxRotation" || key == "maxRotationAngle") return maxRotationAngleRad_;
    if (key == "useCavityBias") return useCavityBias_ ? 1.0 : 0.0;
    if (key == "useConfBias") return useConfBias_ ? 1.0 : 0.0;
    if (key == "storeProbabilities") return shouldStoreProbability() ? 1.0 : 0.0;
    
    return 0.0;  // Default for unknown keys
}

void GCMCEngine::setStatisticsInterval(int interval) {
    statsInterval_ = std::max(1, interval);
    statistics_.setSamplingInterval(statsInterval_);
    configMap_["statsInterval"] = static_cast<double>(statsInterval_);
    // Reset countdown to align with the new interval
    statsCountdown_ = statsInterval_;
}

// Check if should store probability - controlled by configuration
bool GCMCEngine::shouldStoreProbability() const {
    // Use cached value for performance
    if (!storeProbabilityCached_) {
        // First time: check environment and config directly
        // Avoid calling the global function to prevent cross-test contamination
        const char* env_value = std::getenv("GCMC_STORE_PROB");
        if (env_value != nullptr) {
            storeProbabilityValue_ = true;
        } else {
            // Check runtime config keys
            auto it = configMap_.find("storeProbabilities");
            if (it != configMap_.end()) {
                storeProbabilityValue_ = (it->second > 0.5);
            } else {
                auto it2 = configMap_.find("storeAcceptanceProbability");
                if (it2 != configMap_.end()) {
                    storeProbabilityValue_ = (it2->second > 0.5);
                } else {
                    // Default: don't store probabilities for performance
                    storeProbabilityValue_ = false;
                }
            }
        }
        storeProbabilityCached_ = true;
    }
    return storeProbabilityValue_;
}

// CBMC insertion - generate K trials and select based on Boltzmann weights
GCMCEngine::TrialConfiguration GCMCEngine::performCBMCInsertion(int typeId, int numTrials) {
    // Diagnostics: capture the trial energies used to compute log(W/K).
    lastCbmcTrialEnergies_.clear();

    std::vector<TrialConfiguration> trials;
    trials.reserve(numTrials);

    FragmentTemplate* tmpl = reservoir_->getTemplate(typeId);
    if (!tmpl) {
        // Return default configuration if template not found
        return TrialConfiguration{Vector3(0,0,0), Quaternion(1,0,0,0), 0.0, 1.0, 0.0, false, 0};
    }

    // Prepare any background grid once per CBMC move (never per trial).
    if (energyBackend_ == GCMCEnergyBackend::PgpHost) {
        ensurePgpHostGridReady();
    } else if (energyBackend_ == GCMCEnergyBackend::PgpFull ||
               energyBackend_ == GCMCEnergyBackend::PgpFullPme) {
        // Background is the current system (trial residues are temporary and must not enter the grid).
        ensurePgpFullGridReadyExcluding(-1);
    }

    const bool useSystemEnergyDelta = (energyBackend_ == GCMCEnergyBackend::Pme);
    const double baselineSystemEnergy = useSystemEnergyDelta ? calculateSystemEnergy() : 0.0;

    // Generate K trial configurations
    double minEnergy = std::numeric_limits<double>::max();
    for (int k = 0; k < numTrials; ++k) {
        TrialConfiguration trial;
        trial.weight = 0.0;
        trial.logWOverK = 0.0;
        trial.valid = false;
        trial.trialsUsed = 0;

        // Generate position and orientation
        trial.position = (useCavityBias_ && cavityManager_) ?
                        generateCavityPosition() : generateRandomPosition();
        applyPeriodicBoundary(trial.position);
        trial.orientation = generateRandomOrientation();

        if (regionConstraint_ && !isMoleculeWithinRegion(typeId, trial.position, trial.orientation)) {
            continue;
        }

        // Create temporary instance for energy calculation
        int tempId = reservoir_->createInstance(typeId, trial.position, trial.orientation);
        if (tempId < 0) continue;

        synchronizeStateWithReservoir(tempId, true);

        // Calculate energy for this configuration
        if (useSystemEnergyDelta) {
            const double energyWithTrial = calculateSystemEnergy();
            trial.energy = energyWithTrial - baselineSystemEnergy;
        } else if (energyBackend_ == GCMCEnergyBackend::DirectCutoff) {
            cpu::computeResidueEnergyCutoffPBC(*state_, tempId);
            const auto& residue = state_->residues[tempId];
            trial.energy = residue.energy_vdw + residue.energy_elec;
        } else if (energyBackend_ == GCMCEnergyBackend::PgpHost) {
            trial.energy = calculateFragmentEnergyPgpHost(tempId);
        } else if (energyBackend_ == GCMCEnergyBackend::PgpFull ||
                   energyBackend_ == GCMCEnergyBackend::PgpFullPme) {
            trial.energy = calculateFragmentEnergyPgpFullUsingCurrentGrid(tempId);
        } else {
            trial.energy = calculateFragmentEnergy(tempId);
        }

        // Clean up temporary instance
        reservoir_->deleteInstance(tempId);
        synchronizeStateWithReservoir(tempId, false);

        // Track minimum energy for numerical stability
        minEnergy = std::min(minEnergy, trial.energy);
        trial.valid = true;
        trials.push_back(trial);
    }

    // Check if we got any valid trials
    if (trials.empty()) {
        // No valid CBMC trials - return invalid to let caller reject
        return TrialConfiguration{
            generateRandomPosition(),
            generateRandomOrientation(),
            0.0,  // energy
            0.0,  // rosenbluthWeight (invalid)
            0.0,  // log(W/K)
            false,
            0,
        };
    }

    lastCbmcTrialEnergies_.reserve(trials.size());
    for (const auto& trial : trials) {
        lastCbmcTrialEnergies_.push_back(trial.energy);
    }

    // Calculate Boltzmann weights (subtract minEnergy for numerical stability)
    const double beta = 1.0 / (8.314e-3 * temperature_);
    double totalWeight = 0.0;
    for (auto& trial : trials) {
        trial.weight = std::exp(-beta * (trial.energy - minEnergy));
        totalWeight += trial.weight;
    }

    // Select configuration based on weights
    double r = uniform_(rng_) * totalWeight;
    double cumWeight = 0.0;
    TrialConfiguration selected = trials.back();
    for (const auto& trial : trials) {
        cumWeight += trial.weight;
        if (cumWeight >= r) {
            selected = trial;
            break;
        }
    }

    // Convert selection weight into a Rosenbluth factor compatible with exp(-βΔU)
    // to avoid double-counting the selected configuration's Boltzmann term:
    //   rosen = (W/K) / exp(-β u_selected)
    const int effectiveTrials = static_cast<int>(trials.size());
    const double avgScaled = totalWeight / static_cast<double>(effectiveTrials);
    const double safeAvgScaled = std::max(avgScaled, 1e-30);
    const double logWOverK = std::log(safeAvgScaled) - beta * minEnergy;

    // rosen = (W/K)/exp(-β u_selected) = exp(log(W/K) + β u_selected)
    const double logLower = std::log(1e-30);
    const double logUpper = 700.0;
    const double logRosen = logWOverK + beta * selected.energy;
    const double logRosenClamped = std::min(std::max(logRosen, logLower), logUpper);
    selected.weight = std::exp(logRosenClamped);
    selected.logWOverK = logWOverK;
    selected.valid = true;
    selected.trialsUsed = effectiveTrials;

    return selected;
}

// Calculate CBMC bias factor
double GCMCEngine::calculateCBMCBias(const std::vector<TrialConfiguration>& trials, int selectedIdx) {
    if (trials.empty() || selectedIdx < 0 || selectedIdx >= static_cast<int>(trials.size())) {
        return 1.0;
    }

    // Calculate average Boltzmann factor
    double beta = 1.0 / (8.314e-3 * temperature_);
    double minEnergy = std::numeric_limits<double>::max();
    for (const auto& trial : trials) {
        minEnergy = std::min(minEnergy, trial.energy);
    }

    double sumBoltzmann = 0.0;
    for (const auto& trial : trials) {
        sumBoltzmann += std::exp(-beta * (trial.energy - minEnergy));
    }

    // Return W_new / K for insertion
    return (sumBoltzmann / trials.size()) * std::exp(-beta * minEnergy);
}

// Get residue position
Vector3 GCMCEngine::getResiduePosition(int residueIdx) {
    if (!reservoir_) return Vector3(0, 0, 0);
    
    FragmentInstance* instance = reservoir_->getInstance(residueIdx);
    if (!instance) return Vector3(0, 0, 0);
    
    return Vector3(instance->position.x, instance->position.y, instance->position.z);
}

// Get residue orientation
Quaternion GCMCEngine::getResidueOrientation(int residueIdx) {
    if (!reservoir_) return Quaternion(1, 0, 0, 0);
    
    FragmentInstance* instance = reservoir_->getInstance(residueIdx);
    if (!instance) return Quaternion(1, 0, 0, 0);
    
    return instance->orientation;
}

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

/**
 * @file DrudeMain.cpp
 * @brief Main entry point implementation for Drude force calculations
 */

#include "DrudeMain.hpp"
#include "DrudeCore.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

double DrudeComplete::calculateEnergy(model::MCState& state) {
    return DrudeCore::getInstance().calculateEnergy(state);
}

double DrudeComplete::calculateEnergy(model::MCState& state, DrudeAlgorithm algorithm) {
    auto& core = DrudeCore::getInstance();
    core.setAlgorithm(algorithm);
    return core.calculateEnergy(state);
}

void DrudeComplete::setParameters(const DrudeSCFParams& params) {
    DrudeCore::getInstance().setParameters(params);
}

void DrudeComplete::addParticle(const DrudeParticle& particle) {
    DrudeCore::getInstance().addParticle(particle);
}

void DrudeComplete::addScreenedPair(const ScreenedPair& pair) {
    DrudeCore::getInstance().addScreenedPair(pair);
}

void DrudeComplete::clear() {
    DrudeCore::getInstance().clear();
}

size_t DrudeComplete::getNumParticles() {
    return DrudeCore::getInstance().getNumParticles();
}

void DrudeComplete::enableASPC(bool enable) {
    DrudeCore::getInstance().enableASPC(enable);
}

bool DrudeComplete::isASPCEnabled() {
    return DrudeCore::getInstance().isASPCEnabled();
}

void DrudeComplete::clearHistory() {
    DrudeCore::getInstance().clearHistory();
}

DrudeFastFBP* DrudeComplete::getFastFBPOptimizer() {
    return DrudeCore::getInstance().getFastFBPOptimizer();
}

DrudeCore& DrudeComplete::getDrudeCore() {
    return DrudeCore::getInstance();
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
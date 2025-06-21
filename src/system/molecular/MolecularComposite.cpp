#include "MolecularComposite.hpp"
#include "MolecularCombiner.hpp"

namespace pygcmc {
namespace system {
namespace molecular {

MolecularComposite::MolecularComposite() 
    : combiner_(std::make_unique<MolecularCombiner>()) {
}

std::shared_ptr<model::Molecular> MolecularComposite::buildMolecular(
    const std::shared_ptr<model::Structure>& structure,
    const std::shared_ptr<model::Topology>& topology) {
    
    current_molecular_ = combiner_->combine(structure, topology);
    return current_molecular_;
}

std::shared_ptr<model::Molecular> MolecularComposite::buildMolecularMultiple(
    const std::shared_ptr<model::Structure>& structure,
    const std::vector<std::shared_ptr<model::Topology>>& topologies) {
    
    current_molecular_ = combiner_->combineMultiple(structure, topologies);
    return current_molecular_;
}

} // namespace molecular
} // namespace system
} // namespace pygcmc 
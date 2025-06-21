#include "MolecularMain.hpp"
#include "MolecularComposite.hpp"

namespace pygcmc {
namespace system {
namespace molecular {

MolecularMain::MolecularMain() 
    : impl_(std::make_unique<MolecularComposite>()) {
}

MolecularMain::~MolecularMain() = default;

std::shared_ptr<model::Molecular> MolecularMain::combine(
    const std::shared_ptr<model::Structure>& structure,
    const std::shared_ptr<model::Topology>& topology) {
    
    return impl_->buildMolecular(structure, topology);
}

std::shared_ptr<model::Molecular> MolecularMain::combine_multiple(
    const std::shared_ptr<model::Structure>& structure,
    const std::vector<std::shared_ptr<model::Topology>>& topologies) {
    
    return impl_->buildMolecularMultiple(structure, topologies);
}

const std::shared_ptr<model::Molecular>& MolecularMain::get_molecular() const {
    return impl_->getCurrentMolecular();
}

} // namespace molecular
} // namespace system
} // namespace pygcmc 
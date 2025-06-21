#include "MolecularMain.hpp"
#include "MolecularComposite.hpp"

namespace pygcmc {
namespace system {
namespace molecular {

MolecularMain::MolecularMain() 
    : impl_(std::make_unique<MolecularComposite>()) {
}

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

// Backward compatibility: provide the original MolecularSystem class
namespace system {

/**
 * @brief Backward compatible MolecularSystem class
 * 
 * This class maintains the original API while using the new modular implementation.
 */
class MolecularSystem {
public:
    MolecularSystem() = default;
    ~MolecularSystem() = default;

    /**
     * @brief Combine Structure and Topology data into a Molecular object
     * @param structure Structure data
     * @param topology Topology data
     * @return Combined Molecular object
     */
    std::shared_ptr<model::Molecular> combine(
        const std::shared_ptr<model::Structure>& structure,
        const std::shared_ptr<model::Topology>& topology) {
        return impl_.combine(structure, topology);
    }

    /**
     * @brief Combine multiple Structure and Topology data into a Molecular object
     * @param structure Structure data
     * @param topologies Topology data
     * @return Combined Molecular object
     */
    std::shared_ptr<model::Molecular> combine_multiple(
        const std::shared_ptr<model::Structure>& structure,
        const std::vector<std::shared_ptr<model::Topology>>& topologies) {
        return impl_.combine_multiple(structure, topologies);
    }

    /**
     * @brief Get the current Molecular object
     * @return Current Molecular object
     */
    const std::shared_ptr<model::Molecular>& get_molecular() const { 
        return impl_.get_molecular(); 
    }

private:
    molecular::MolecularMain impl_;
};

} // namespace system
} // namespace pygcmc 
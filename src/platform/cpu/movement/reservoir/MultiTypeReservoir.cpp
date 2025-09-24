#include "MultiTypeReservoir.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

void MultiTypeReservoir::addType(const TypeInfo& info, const FragmentTemplate& tpl) {
    TypeInfo newInfo = info;
    newInfo.currentCount = 0;  // Initialize count to 0
    size_t idx = types_.size();
    types_.push_back(newInfo);
    idToIdx_[newInfo.typeId] = idx;  // Fill idToIdx_ mapping
    FragmentReservoir::addTemplate(tpl);
}

int MultiTypeReservoir::selectTypeForInsertion(std::mt19937& rng) const {
    if (types_.empty()) return -1;
    
    // Build list of available types (not at max capacity)
    std::vector<double> weights;
    std::vector<int> availableTypes;
    
    for (const auto& t : types_) {
        // Use actual active counts from base reservoir to enforce caps
        const int active = activeCount(t.typeId);
        if (active < t.maxCount) {  // Skip types at max capacity
            weights.push_back(t.probability);
            availableTypes.push_back(t.typeId);
        }
    }
    
    // All types at max capacity
    if (availableTypes.empty()) return -1;
    
    // Single available type
    if (availableTypes.size() == 1) return availableTypes[0];
    
    // Select based on probability weights (no chemical potential weighting)
    double sum = 0.0;
    for (double w : weights) sum += w;
    
    if (sum <= 0.0) return availableTypes[0];
    
    std::uniform_real_distribution<double> u(0.0, sum);
    double r = u(rng);
    double acc = 0.0;
    
    for (size_t i = 0; i < weights.size(); ++i) {
        acc += weights[i];
        if (r <= acc) return availableTypes[i];
    }
    
    return availableTypes.back();
}

int MultiTypeReservoir::selectInstanceForDeletion(std::mt19937& rng) const {
    auto all = getActiveInstances(-1);
    if (all.empty()) return -1;
    std::uniform_int_distribution<size_t> u(0, all.size() - 1);
    const int inst = all[u(rng)];
    auto* p = getInstance(inst);
    if (!p) return -1;
    return p->templateId;
}

int MultiTypeReservoir::activeCount(int typeId) const {
    return getActiveCount(typeId);
}

// Internal helper methods implementation
MultiTypeReservoir::TypeInfo* MultiTypeReservoir::getTypeById(int typeId) {
    auto it = idToIdx_.find(typeId);
    if (it == idToIdx_.end()) {
        return nullptr;
    }
    return &types_[it->second];
}

void MultiTypeReservoir::incCount(int typeId) {
    TypeInfo* info = getTypeById(typeId);
    if (info && info->currentCount < info->maxCount) {
        info->currentCount++;
    }
}

void MultiTypeReservoir::decCount(int typeId) {
    TypeInfo* info = getTypeById(typeId);
    if (info && info->currentCount > 0) {
        info->currentCount--;
    }
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

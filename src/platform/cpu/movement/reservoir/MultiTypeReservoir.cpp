#include "MultiTypeReservoir.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

void MultiTypeReservoir::addType(const TypeInfo& info, const FragmentTemplate& tpl) {
    types_.push_back(info);
    FragmentReservoir::addTemplate(tpl);
}

int MultiTypeReservoir::selectTypeForInsertion(std::mt19937& rng) const {
    if (types_.empty()) return -1;
    double sum = 0.0;
    for (auto& t : types_) sum += t.probability;
    if (sum <= 0.0) return types_.front().typeId;
    std::uniform_real_distribution<double> u(0.0, sum);
    double r = u(rng);
    double acc = 0.0;
    for (auto& t : types_) {
        acc += t.probability;
        if (r <= acc) return t.typeId;
    }
    return types_.back().typeId;
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

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc
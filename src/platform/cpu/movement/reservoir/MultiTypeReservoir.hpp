#pragma once
#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_RESERVOIR_MULTITYPERESERVOIR_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_RESERVOIR_MULTITYPERESERVOIR_HPP

#include "fragment_reservoir.hpp"
#include <random>
#include <unordered_map>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

class MultiTypeReservoir : public FragmentReservoir {
public:
    struct TypeInfo {
        int typeId = -1;
        std::string name;
        double chemicalPotential = 0.0;
        double activity = 1.0;
        double probability = 1.0;
        int maxCount = 0;
        int currentCount = 0;
        double radius = 0.0;
    };

    void addType(const TypeInfo& info, const FragmentTemplate& tpl);
    int selectTypeForInsertion(std::mt19937& rng) const;
    int selectInstanceForDeletion(std::mt19937& rng) const;

    const std::vector<TypeInfo>& types() const { return types_; }
    int activeCount(int typeId) const;

private:
    std::vector<TypeInfo> types_;
    std::unordered_map<int, size_t> idToIdx_;  // typeId to index in types_ vector
    
    // Internal helper methods (not exposed to Python)
    TypeInfo* getTypeById(int typeId);
    void incCount(int typeId);
    void decCount(int typeId);
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_RESERVOIR_MULTITYPERESERVOIR_HPP
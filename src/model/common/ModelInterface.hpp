#pragma once

#ifndef PYGCMC_MODEL_COMMON_INTERFACE_HPP
#define PYGCMC_MODEL_COMMON_INTERFACE_HPP

#include <memory>
#include <string>

namespace pygcmc {
namespace model {
namespace common {

/**
 * @brief Interface for objects that can be cloned
 */
template<typename T>
class ICloneable {
public:
    virtual ~ICloneable() = default;
    virtual std::unique_ptr<T> clone() const = 0;
};

/**
 * @brief Interface for objects that can be serialized
 */
class ISerializable {
public:
    virtual ~ISerializable() = default;
    virtual std::string serialize() const = 0;
    virtual void deserialize(const std::string& data) = 0;
};

/**
 * @brief Interface for objects that can be validated
 */
class IValidatable {
public:
    virtual ~IValidatable() = default;
    virtual bool is_valid() const = 0;
    virtual std::string get_validation_error() const { return ""; }
};

/**
 * @brief Interface for objects with numeric identifiers
 */
class IIdentifiable {
public:
    virtual ~IIdentifiable() = default;
    virtual int get_id() const = 0;
    virtual void set_id(int id) = 0;
};

} // namespace common
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_COMMON_INTERFACE_HPP 
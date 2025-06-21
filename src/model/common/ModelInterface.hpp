#pragma once

#ifndef PYGCMC_MODEL_INTERFACE_HPP
#define PYGCMC_MODEL_INTERFACE_HPP

#include <memory>
#include <string>

namespace pygcmc {
namespace model {

/**
 * @brief Base interface for cloneable objects
 */
template<typename T>
class ICloneable {
public:
    virtual ~ICloneable() = default;
    virtual std::unique_ptr<T> clone() const = 0;
};

/**
 * @brief Base interface for serializable objects  
 */
class ISerializable {
public:
    virtual ~ISerializable() = default;
    
    /**
     * @brief Serialize object to string
     */
    virtual std::string serialize() const = 0;
    
    /**
     * @brief Deserialize object from string
     * @param data Serialized data
     * @return true if successful
     */
    virtual bool deserialize(const std::string& data) = 0;
    
    /**
     * @brief Get object type name for serialization
     */
    virtual std::string get_type_name() const = 0;
};

/**
 * @brief Base interface for validatable objects
 */
class IValidatable {
public:
    virtual ~IValidatable() = default;
    
    /**
     * @brief Check if object is in valid state
     */
    virtual bool is_valid() const = 0;
    
    /**
     * @brief Get validation error message
     */
    virtual std::string get_validation_error() const { return ""; }
};

/**
 * @brief Base interface for objects with unique identifiers
 */
class IIdentifiable {
public:
    virtual ~IIdentifiable() = default;
    
    /**
     * @brief Get unique identifier
     */
    virtual std::string get_id() const = 0;
    
    /**
     * @brief Set unique identifier
     */
    virtual void set_id(const std::string& id) = 0;
};

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_INTERFACE_HPP 
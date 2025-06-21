/**
 * @file system.hpp
 * @brief Legacy compatibility header for System
 * 
 * This file provides backward compatibility by including the new modular
 * System implementation. All original functionality is preserved.
 */

#pragma once

// Include the new modular implementation
#include "common/SystemMain.hpp"

// Re-export the System class for backward compatibility
// The class is now implemented in common/SystemMain.hpp but the API remains the same

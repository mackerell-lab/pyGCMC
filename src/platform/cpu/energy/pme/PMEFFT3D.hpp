#pragma once

#include "PMEFFTCore.hpp"
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace CustomFFT {

/**
 * @brief Perform 1D FFT along specified dimension of 3D grid
 */
void fft_1d_batch(cmplx* data, int dimension, int nx, int ny, int nz, bool inverse);

/**
 * @brief 3D forward FFT - requires grid dimensions to be powers of 2
 */
void fft3D_forward(cmplx* data, int nx, int ny, int nz);

/**
 * @brief 3D inverse FFT - requires grid dimensions to be powers of 2
 */
void fft3D_backward(cmplx* data, int nx, int ny, int nz);

} // namespace CustomFFT
} // namespace cpu
} // namespace platform
} // namespace pygcmc

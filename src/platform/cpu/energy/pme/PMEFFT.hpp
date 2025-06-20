#pragma once

#include <complex>
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {

// Type alias for complex numbers
using cmplx = std::complex<double>;

/**
 * @brief Custom FFT implementation for PME calculations
 * 
 * This namespace provides custom FFT functions optimized for PME grid operations.
 * The implementation uses radix-2 FFT algorithms and is designed to work with
 * power-of-2 grid sizes commonly used in PME.
 */
namespace CustomFFT {

// Core FFT algorithm functions
int bit_reverse(int x, int n);
void apply_bit_reverse(cmplx* A, int k);
void tfft_genw(int i, int b, cmplx z, cmplx *w);
void tfft_init(int k, cmplx *w);
void tfft_fft(int k, cmplx *A, const cmplx *w);
void tfft_ifft(int k, cmplx *A, const cmplx *w);

// Utility functions
int log2_power_of_2(int n);
void padded_fft(cmplx* data, int actual_size, bool inverse);

// 1D batch FFT operations for 3D transformations
void fft_1d_batch(cmplx* data, int dimension, int nx, int ny, int nz, bool inverse);

// High-level 3D FFT interfaces
void fft3D_forward(cmplx* data, int nx, int ny, int nz);
void fft3D_backward(cmplx* data, int nx, int ny, int nz);

// Specialized convolution operations (optional)
void tfft_convolver(int k, cmplx *A, const cmplx *w);

} // namespace CustomFFT

// <agent-hook:pme_fft>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 
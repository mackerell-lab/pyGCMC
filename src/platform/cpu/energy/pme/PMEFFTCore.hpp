#pragma once

#include <complex>
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {

// Type alias for complex numbers
using cmplx = std::complex<double>;

namespace CustomFFT {

// Constants
extern const double PI2;

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

// Specialized convolution operations
void tfft_convolver(int k, cmplx *A, const cmplx *w);

} // namespace CustomFFT
} // namespace cpu
} // namespace platform
} // namespace pygcmc
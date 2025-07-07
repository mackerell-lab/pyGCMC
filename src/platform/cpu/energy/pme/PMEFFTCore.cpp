#include "PMEFFTCore.hpp"
#include "platform/platform.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace CustomFFT {

// Define PI2 constant
const double PI2 = 6.28318530717958647692;

// Static weights vector for FFT operations
static std::vector<cmplx> weights;

/**
 * @brief Calculate bit-reversed value of n-bit number x
 */
int bit_reverse(int x, int n) {
    int result = 0;
    for (int i = 0; i < n; i++) {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    return result;
}

/**
 * @brief Apply bit-reverse sorting to make FFT results compatible with FFTW
 */
void apply_bit_reverse(cmplx* A, int k) {
    const int m = 1 << k;
    std::vector<cmplx> temp(m);
    
    // Copy data to temporary array, reordering by bit-reversed order
    for (int i = 0; i < m; i++) {
        int j = bit_reverse(i, k);
        temp[i] = A[j];
    }
    
    // Copy back to original array
    for (int i = 0; i < m; i++) {
        A[i] = temp[i];
    }
}

/**
 * @brief Generate FFT weights recursively
 */
void tfft_genw(int i, int b, cmplx z, cmplx *w) {
    if(b == 0)
        w[i] = z;
    else {
        tfft_genw(i, b>>1, z, w);
        tfft_genw(i|b, b>>1, z*w[b], w);
    }
}

/**
 * @brief Initialize FFT weights
 */
void tfft_init(int k, cmplx *w) {
    int i, j;
    const int m = 1<<k;
    const double arg = -PI2/m;
    for(i=1, j=m/4; j; i<<=1, j>>=1) {
        w[i] = std::exp(std::complex<double>(0, arg * j));
    }
    tfft_genw(0, m/4, 1, w);
}

/**
 * @brief Perform forward FFT
 */
void tfft_fft(int k, cmplx *A, const cmplx *w) {
    const int m = 1 << k;
    int u = 1;
    int v = m/4;
    int i, j;
    if(k&1) {
        for(j=0; j<m/2; j++) {
            cmplx Ajv = A[j+(m/2)];
            A[j+(m/2)] = A[j] - Ajv;
            A[j] += Ajv;
        }
        u <<= 1;
        v >>= 1;
    }
    for(i=k&~1;i>0;i-=2) {
        int jh;
        for(jh=0;jh<u;jh++) {
            cmplx wj = w[jh<<1];
            cmplx wj2 = w[jh];
            cmplx wj3 = wj2 * wj;
            int je;
            for(j = jh << i, je = j+v;j<je; j++) {
                cmplx tmp0 = A[j];
                cmplx tmp1 = wj * A[j+v];
                cmplx tmp2 = wj2 * A[j+2*v];
                cmplx tmp3 = wj3 * A[j+3*v];

                cmplx ttmp0 = tmp0 + tmp2;
                cmplx ttmp2 = tmp0 - tmp2;
                cmplx ttmp1 = tmp1 + tmp3;
                cmplx ttmp3 = -std::complex<double>(0, 1) * (tmp1 - tmp3);

                A[j] = ttmp0 + ttmp1;
                A[j+v] = ttmp0 - ttmp1;
                A[j+2*v] = ttmp2 + ttmp3;
                A[j+3*v] = ttmp2 - ttmp3;
            }
        }
        u <<= 2;
        v >>= 2;
    }
    
    // Add bit-reverse sorting to make results compatible with FFTW
    apply_bit_reverse(A, k);
}

/**
 * @brief Perform inverse FFT using conjugate method
 */
void tfft_ifft(int k, cmplx *A, const cmplx *w) {
    const int m = 1 << k;
    
    // Step 1: Conjugate input data
    for (int i = 0; i < m; i++) {
        A[i] = std::conj(A[i]);
    }
    
    // Step 2: Use forward FFT transform
    tfft_fft(k, A, w);
    
    // Step 3: Conjugate result again and normalize
    for (int i = 0; i < m; i++) {
        A[i] = std::conj(A[i]) / static_cast<double>(m);
    }
}

/**
 * @brief Calculate log2(n), n must be a power of 2
 */
int log2_power_of_2(int n) {
    int k = 0;
    while (n > 1) {
        n >>= 1;
        k++;
    }
    return k;
}

/**
 * @brief Perform FFT with automatic padding for power-of-2 sizes
 */
void padded_fft(cmplx* data, int actual_size, bool inverse) {
    // Check if it's a power of two
    bool is_power_of_two = (actual_size & (actual_size - 1)) == 0;
    
    if (!is_power_of_two) {
        throw std::runtime_error("FFT size must be a power of 2");
    }
    
    // Calculate k where 2^k = actual_size
    int k = log2_power_of_2(actual_size);
    
    // Ensure weights vector is large enough
    if (weights.size() < static_cast<size_t>(actual_size)) {
        weights.resize(actual_size);
        tfft_init(k, weights.data());
    }
    
    // Execute FFT
    if (inverse) {
        tfft_ifft(k, data, weights.data());
    } else {
        tfft_fft(k, data, weights.data());
    }
}

/**
 * @brief Convolver function for specialized operations
 */
void tfft_convolver(int k, cmplx *A, const cmplx *w) {
    int i, y;
    const int m = 1 << k;

    tfft_fft(k, A, w);
    A[0] = 4 * std::real(A[0]) * std::imag(A[0]) * std::complex<double>(0, 1);
    A[1] = 4 * std::real(A[1]) * std::imag(A[1]) * std::complex<double>(0, 1);
    i = 2;
    for(y = 2; y < m; y <<= 1) {
        for(; i < 2*y; i+=2) {
            int j = i^(y-1);
            A[i] = (A[i] + std::conj(A[j]))*(A[i] - std::conj(A[j]));
            A[j] = -std::conj(A[i]);
        }
    }

    for(i = 0; i < m; i+=2) {
        double scale = 4.0 * m;
        A[i/2] = (-(A[i]+A[i^1])*std::complex<double>(0, 1) + 
                  (A[i]-A[i^1])*std::conj(w[i/2]))/scale;
    }

    tfft_ifft(k-1, A, w);
}

/**
 * @brief Clear static FFT weights to force regeneration on next use
 */
void clearFFTWeights() {
    weights.clear();
    weights.shrink_to_fit();
}

// <agent-hook:fft_implementation>

} // namespace CustomFFT
} // namespace cpu
} // namespace platform
} // namespace pygcmc

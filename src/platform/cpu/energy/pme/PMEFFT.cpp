#include "PMEFFT.hpp"
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
 * @brief Perform 1D FFT along specified dimension of 3D grid
 */
void fft_1d_batch(cmplx* data, int dimension, int nx, int ny, int nz, bool inverse) {
    // Choose transform size and processing method based on dimension
    int transform_size;
    int num_transforms;
    
    switch (dimension) {
        case 0: // X direction
            transform_size = nx;
            num_transforms = ny * nz;
            break;
        case 1: // Y direction
            transform_size = ny;
            num_transforms = nx * nz;
            break;
        case 2: // Z direction
            transform_size = nz;
            num_transforms = nx * ny;
            break;
        default:
            throw std::invalid_argument("Invalid dimension for 3D FFT");
    }
    
    // Check if transform size is a power of 2
    if ((transform_size & (transform_size - 1)) != 0) {
        throw std::runtime_error("FFT size must be a power of 2");
    }
    
    // Process each 1D transform
    for (int t = 0; t < num_transforms; t++) {
        int y, z, x;
        
        // Temporary buffer for single transform
        std::vector<cmplx> buffer(transform_size);
        
        // Calculate corresponding 2D index based on dimension
        switch (dimension) {
            case 0: // X direction
                y = t % ny;
                z = t / ny;
                
                // Read data into buffer
                for (int x = 0; x < nx; x++) {
                    buffer[x] = data[x + y*nx + z*nx*ny];
                }
                
                // Perform FFT/IFFT
                padded_fft(buffer.data(), transform_size, inverse);
                
                // Write result back
                for (int x = 0; x < nx; x++) {
                    data[x + y*nx + z*nx*ny] = buffer[x];
                }
                break;
                
            case 1: // Y direction
                x = t % nx;
                z = t / nx;
                
                // Read data into buffer
                for (int y = 0; y < ny; y++) {
                    buffer[y] = data[x + y*nx + z*nx*ny];
                }
                
                // Perform FFT/IFFT
                padded_fft(buffer.data(), transform_size, inverse);
                
                // Write result back
                for (int y = 0; y < ny; y++) {
                    data[x + y*nx + z*nx*ny] = buffer[y];
                }
                break;
                
            case 2: // Z direction
                x = t % nx;
                y = t / nx;
                
                // Read data into buffer
                for (int z = 0; z < nz; z++) {
                    buffer[z] = data[x + y*nx + z*nx*ny];
                }
                
                // Perform FFT/IFFT
                padded_fft(buffer.data(), transform_size, inverse);
                
                // Write result back
                for (int z = 0; z < nz; z++) {
                    data[x + y*nx + z*nx*ny] = buffer[z];
                }
                break;
        }
    }
}

/**
 * @brief 3D forward FFT - requires grid dimensions to be powers of 2
 */
void fft3D_forward(cmplx* data, int nx, int ny, int nz) {
    // Execute three 1D FFTs in X, Y, Z order
    fft_1d_batch(data, 0, nx, ny, nz, false); // X direction
    fft_1d_batch(data, 1, nx, ny, nz, false); // Y direction
    fft_1d_batch(data, 2, nx, ny, nz, false); // Z direction
}

/**
 * @brief 3D inverse FFT - requires grid dimensions to be powers of 2
 */
void fft3D_backward(cmplx* data, int nx, int ny, int nz) {
    // Execute three 1D IFFTs in Z, Y, X order (reverse of forward order)
    fft_1d_batch(data, 2, nx, ny, nz, true); // Z direction
    fft_1d_batch(data, 1, nx, ny, nz, true); // Y direction
    fft_1d_batch(data, 0, nx, ny, nz, true); // X direction
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

// <agent-hook:fft_implementation>

} // namespace CustomFFT
} // namespace cpu
} // namespace platform
} // namespace pygcmc 
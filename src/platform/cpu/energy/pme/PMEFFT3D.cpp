#include "PMEFFT3D.hpp"
#include "PMEFFTCore.hpp"
#include "platform/platform.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace CustomFFT {

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

} // namespace CustomFFT
} // namespace cpu
} // namespace platform
} // namespace pygcmc

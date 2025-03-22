#include "energyPME.hpp"
#include "platform/platform.hpp"
#include "model/residue.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <complex>
#include <chrono>
#include <set>
#include <tuple>

// Optional: Include library for FFT if needed
// #include <fftw3.h>

// Define complex number type alias
using cmplx = std::complex<double>;
// Coulomb constant (kcal·mol⁻¹·e⁻²)
static const double COULOMB = 332.0716;

namespace {
    // FFTW Plan related constants
    static const unsigned FFTW_MEASURE = 0;
    static const int FFTW_FORWARD = -1;
    static const int FFTW_BACKWARD = 1;
    
    // Timer for diagnostics
    auto t_start = std::chrono::high_resolution_clock::now();
}

namespace {
// Define PI2
const double PI2 = 6.28318530717958647692;

// In C++, we use std::complex<double> instead of C's complex type
using cmplx = std::complex<double>;

namespace CustomFFT {

// Add static weights vector
static std::vector<cmplx> weights;

// Calculate bit-reversed value of n-bit number x
int bit_reverse(int x, int n) {
    int result = 0;
    for (int i = 0; i < n; i++) {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    return result;
}

// Apply bit-reverse sorting to make FFT results compatible with FFTW
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

// Generate FFT weights
void tfft_genw(int i, int b, cmplx z, cmplx *w) {
    if(b == 0)
        w[i] = z;
    else {
        tfft_genw(i, b>>1, z, w);
        tfft_genw(i|b, b>>1, z*w[b], w);
    }
}

// Initialize FFT weights
void tfft_init(int k, cmplx *w) {
    int i, j;
    const int m = 1<<k;
    const double arg = -PI2/m;
    for(i=1, j=m/4; j; i<<=1, j>>=1) {
        w[i] = std::exp(std::complex<double>(0, arg * j));
    }
    tfft_genw(0, m/4, 1, w);
}

// Perform forward FFT
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

// Perform inverse FFT (using conjugate method to ensure FFTW compatibility)
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

// Restore log2_power_of_2 function
// Calculate log2(n), n must be a power of 2
int log2_power_of_2(int n) {
    int k = 0;
    while (n > 1) {
        n >>= 1;
        k++;
    }
    return k;
}

// Fix padded_fft function, using log2_power_of_2 to calculate k value
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
        tfft_init(k, weights.data());  // Use k instead of actual_size
    }
    
    // Execute FFT
    if (inverse) {
        tfft_ifft(k, data, weights.data());  // Use k instead of actual_size
    } else {
        tfft_fft(k, data, weights.data());  // Use k instead of actual_size
    }
}

// Improved fft_1d_batch function, only supporting power-of-2 sized grids
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
        
        // Calculate corresponding 2D index based on dimension (t -> x,y,z, two of them)
        switch (dimension) {
            case 0: // X direction: t = y + z*ny, transform for each y,z pair, iterate over all x
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
                
            case 1: // Y direction: t = x + z*nx, transform for each x,z pair, iterate over all y
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
                
            case 2: // Z direction: t = x + y*nx, transform for each x,y pair, iterate over all z
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

// 3D forward FFT - requires grid dimensions to be powers of 2
void fft3D_forward(cmplx* data, int nx, int ny, int nz) {
    // Execute three 1D FFTs in X, Y, Z order
    fft_1d_batch(data, 0, nx, ny, nz, false); // X direction
    fft_1d_batch(data, 1, nx, ny, nz, false); // Y direction
    fft_1d_batch(data, 2, nx, ny, nz, false); // Z direction
}

// 3D inverse FFT - requires grid dimensions to be powers of 2
void fft3D_backward(cmplx* data, int nx, int ny, int nz) {
    // Execute three 1D IFFTs in Z, Y, X order (reverse of forward order)
    fft_1d_batch(data, 2, nx, ny, nz, true); // Z direction
    fft_1d_batch(data, 1, nx, ny, nz, true); // Y direction
    fft_1d_batch(data, 0, nx, ny, nz, true); // X direction
}

// Convolver function
[[maybe_unused]]
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
        A[i/2] = (-(A[i]+A[i^1])*std::complex<double>(0, 1) + (A[i]-A[i^1])*std::conj(w[i/2]))/scale;
    }

    tfft_ifft(k-1, A, w);
}

} // namespace CustomFFT

} // anonymous namespace

namespace pygcmc {
namespace platform {
namespace cpu {

// Global parameters instance
PMEParams pme_params;

// Add FFT related member variables
std::vector<std::complex<double>> fft_weights;

// Global variable to hold FFT'd grid data
std::vector<std::complex<double>> fftGridBackup;

/**
 * @brief Initialize lookup tables for erfc and scaling functions
 * 
 * @param cutoff Cutoff distance
 */
void PMEParams::initializeTables(double cutoff) {
    this->cutoff = cutoff;
    initialized = true;
    
    // Initialize erfc table for optimization
    erfcTable.resize(NUM_TABLE_POINTS);
    ewaldScaleTable.resize(NUM_TABLE_POINTS);
    
    double tableRange = cutoff;
    ewaldDX = tableRange / (NUM_TABLE_POINTS - 1);
    ewaldDXInv = (NUM_TABLE_POINTS - 1) / tableRange;
    erfcDXInv = ewaldDXInv;
    
    // Populate lookup tables
    for (int i = 0; i < NUM_TABLE_POINTS; i++) {
        double r = i * ewaldDX;
        double alphaR = alpha * r;
        erfcTable[i] = std::erfc(alphaR);
        
        // Store the scale factor for electrostatics
        if (r > 1e-6) {
            ewaldScaleTable[i] = erfcTable[i] / r;
        } else {
            ewaldScaleTable[i] = 2.0 * alpha / std::sqrt(M_PI);
        }
    }
    
    platform::log(LogLevel::INFO, "PME tables initialized with cutoff = ", cutoff,
                 ", alpha = ", alpha, ", mesh size = [", 
                 meshSize[0], ",", meshSize[1], ",", meshSize[2], "]");
}

/**
 * @brief Set box dimensions for PME calculations
 * 
 * @param newBox New box dimensions
 */
void PMEParams::setBox(const double newBox[3]) {
    for (int i = 0; i < 3; i++) {
        box[i] = newBox[i];
    }
    platform::log(LogLevel::INFO, "PME box dimensions set to [", box[0], ", ", box[1], ", ", box[2], "]");
}

/**
 * @brief Initialize B-splines for PME - exactly following pme_calculate_bsplines_moduli implementation in pme.cpp
 */
void PMEParams::initializeBsplines() {
    // Exactly replicate pme_calculate_bsplines_moduli implementation in pme.cpp
    platform::log(LogLevel::INFO, "Initializing B-splines with order = ", splineOrder, 
                 " and mesh size = [", meshSize[0], ",", meshSize[1], ",", meshSize[2], "]");
    
    // Ensure spline order is at least 2
    if (splineOrder < 2) {
        platform::log(LogLevel::WARNING, "B-spline order must be at least 2, setting to 2");
        splineOrder = 2;
    }
    
    // Calculate volume and boxfactor - consistent with orthorhombic box in pme.cpp
    double boxVolume = box[0] * box[1] * box[2];
    
    // If volume near zero (box not set), use unit volume
    if (boxVolume < 1e-10) {
        platform::log(LogLevel::WARNING, "Box volume near zero, using unit volume for boxfactor");
        boxVolume = 1.0;
    }
    
    platform::log(LogLevel::INFO, "Initializing B-splines with box volume = ", boxVolume);
    
    // Find maximum grid size
    int nmax = 0;
    for (int dim = 0; dim < 3; dim++) {
        nmax = (meshSize[dim] > nmax) ? meshSize[dim] : nmax;
        bsplineModuli[dim].resize(meshSize[dim]);
    }
    
    // Exactly replicate initialization logic from pme.cpp
    std::vector<double> data(splineOrder, 0.0);
    std::vector<double> ddata(splineOrder, 0.0);
    std::vector<double> bsplines_data(nmax, 0.0);
    
    // Exactly follow initialization order from pme.cpp
    data[splineOrder-1] = 0.0; // Explicitly set tail to 0
    data[1] = 0.0; // Explicitly set data[1] to 0, same as pme.cpp
    data[0] = 1.0; // Initial condition
    
    // Calculate B-spline coefficients - exactly replicate pme.cpp
    for (int k = 3; k < splineOrder; k++) {
        double div = 1.0/(k-1.0);
        data[k-1] = 0.0;
        for (int l = 1; l < (k-1); l++) {
            data[k-l-1] = div*(l*data[k-l-2] + (k-l)*data[k-l-1]);
        }
        data[0] = div*data[0];
    }
    
    // Calculate derivative - exactly replicate pme.cpp
    ddata[0] = -data[0];
    for (int k = 1; k < splineOrder; k++) {
        ddata[k] = data[k-1] - data[k];
    }
    
    // Calculate final coefficients - exactly replicate pme.cpp
    double div = 1.0/(splineOrder-1.0); // Ensure floating-point division
    data[splineOrder-1] = 0.0;
    for (int l = 1; l < (splineOrder-1); l++) {
        data[splineOrder-l-1] = div*(l*data[splineOrder-l-2] + (splineOrder-l)*data[splineOrder-l-1]);
    }
    data[0] = div*data[0];
    
    // Initialize bsplines_data - exactly replicate pme.cpp
    for (int i = 0; i < nmax; i++) {
        bsplines_data[i] = 0.0;
    }
    for (int i = 1; i <= splineOrder; i++) {
        bsplines_data[i] = data[i-1];
    }
    
    // Calculate B-spline moduli for each dimension - exactly replicate pme.cpp
    for (int dim = 0; dim < 3; dim++) {
        int ndata = meshSize[dim];
        for (int i = 0; i < ndata; i++) {
            double sc = 0.0, ss = 0.0;
            for (int j = 0; j < ndata; j++) {
                double arg = (2.0*M_PI*i*j)/ndata;
                sc += bsplines_data[j]*cos(arg);
                ss += bsplines_data[j]*sin(arg);
            }
            bsplineModuli[dim][i] = sc*sc + ss*ss;
        }
        
        // Improve numerical stability - exactly replicate pme.cpp
        for (int i = 0; i < ndata; i++) {
            if (bsplineModuli[dim][i] < 1.0e-7) {
                bsplineModuli[dim][i] = (bsplineModuli[dim][(i-1+ndata)%ndata] + 
                                      bsplineModuli[dim][(i+1)%ndata])/2.0;
            }
        }
    }
    
    // Output B-spline moduli and statistics
    double maxModuli[3] = {0.0, 0.0, 0.0};
    double minModuli[3] = {std::numeric_limits<double>::max(), 
                          std::numeric_limits<double>::max(), 
                          std::numeric_limits<double>::max()};
    
    for (int dim = 0; dim < 3; dim++) {
        for (size_t i = 0; i < bsplineModuli[dim].size(); i++) {
            maxModuli[dim] = std::max(maxModuli[dim], bsplineModuli[dim][i]);
            minModuli[dim] = std::min(minModuli[dim], bsplineModuli[dim][i]);
        }
        platform::log(LogLevel::DEBUG, "Dimension " + std::to_string(dim) + " B-spline moduli range: [" 
                     + std::to_string(minModuli[dim]) + ", " + std::to_string(maxModuli[dim]) + "]");
    }
    
    // Key debug output: display B-spline moduli values at [5,5,5] point
    if (meshSize[0] > 5 && meshSize[1] > 5 && meshSize[2] > 5) {
        platform::log(LogLevel::DEBUG, "B-spline moduli at [5,5,5]: ["
                     + std::to_string(bsplineModuli[0][5]) + ", "
                     + std::to_string(bsplineModuli[1][5]) + ", " 
                     + std::to_string(bsplineModuli[2][5]) + "]");
    }
    
    // Allocate PME grid
    int totalGridPoints = meshSize[0] * meshSize[1] * meshSize[2];
    pmeGrid.resize(totalGridPoints);
    
    // Clear grid
    std::fill(pmeGrid.begin(), pmeGrid.end(), std::complex<double>(0.0, 0.0));
    
    platform::log(LogLevel::INFO, "PME B-splines initialized with order = ", splineOrder);
}

/**
 * @brief Approximate erfc function using lookup table
 * 
 * @param r Distance
 * @return double erfc(alpha*r)
 */
double PMEParams::erfcApprox(double r) const {
    if (r >= cutoff) return 0.0;
    
    double x = r * erfcDXInv;
    int index = static_cast<int>(x);
    double fraction = x - index;
    
    // Linear interpolation
    return erfcTable[index] + fraction * (erfcTable[index+1] - erfcTable[index]);
}

/**
 * @brief Approximate electrostatic scaling function using lookup table
 * 
 * @param r Distance
 * @return double erfc(alpha*r)/r
 */
double PMEParams::ewaldScaleApprox(double r) const {
    if (r >= cutoff) return 0.0;
    
    double x = r * ewaldDXInv;
    int index = static_cast<int>(x);
    double fraction = x - index;
    
    // Linear interpolation
    return ewaldScaleTable[index] + fraction * (ewaldScaleTable[index+1] - ewaldScaleTable[index]);
}

/**
 * @brief Auto-adjust PME parameters to achieve desired accuracy
 * 
 * @param error_tolerance Desired error tolerance
 * @param cutoff_distance Cutoff distance
 * @param box Box dimensions
 */
void autoAdjustPMEParameters(double error_tolerance, double cutoff_distance, const double box[3]) {
    // Based on OpenMM's autoAdjustParameters implementation
    
    platform::log(LogLevel::INFO, "Auto-adjusting PME parameters for error tolerance ", 
                 error_tolerance, " and cutoff ", cutoff_distance);
    
    // Find the optimal alpha parameter and mesh dimensions

    // First determine alpha to achieve the desired real-space error
    double alpha = 0.0;
    
    // Determine alpha by targeting the real-space error
    // Real-space error: erfc(alpha*cutoff)
    double realSpaceError = error_tolerance / 2.0; // Split error between real and reciprocal space
    
    // Solve for alpha: erfc(alpha*cutoff) = realSpaceError
    // Use approximate inverse of erfc function
    if (realSpaceError < 1e-6) {
        alpha = 3.5 / cutoff_distance;  // Minimum reasonable alpha
    } else {
        double x = -std::log(realSpaceError);
        alpha = std::sqrt(x) / cutoff_distance;
        
        // Refine the value with a few Newton iterations
        for (int i = 0; i < 3; i++) {
            double currentError = std::erfc(alpha * cutoff_distance);
            double derivative = -2.0 / std::sqrt(M_PI) * std::exp(-alpha*alpha*cutoff_distance*cutoff_distance) * cutoff_distance;
            alpha -= (currentError - realSpaceError) / derivative;
        }
    }
    
    // New: limit the range of alpha parameter change to reduce wild fluctuations in self-energy with error tolerance
    // Use a reasonable range of alpha (1.0 - 3.0)/cutoff
    double minAlpha = 1.0 / cutoff_distance;
    double maxAlpha = 3.0 / cutoff_distance;
    
    if (alpha < minAlpha) {
        platform::log(LogLevel::WARNING, "PME alpha parameter too small, adjusting from ", 
                    alpha, " to ", minAlpha);
        alpha = minAlpha;
    } else if (alpha > maxAlpha) {
        platform::log(LogLevel::WARNING, "PME alpha parameter too large, adjusting from ", 
                    alpha, " to ", maxAlpha);
        alpha = maxAlpha;
    }
    
    // Now determine the mesh dimensions
    int meshSize[3];
    double minBoxSize = std::min(box[0], std::min(box[1], box[2]));
    
    // Base mesh size on the desired reciprocal-space error and size of the system
    double reciprocalSpaceError = error_tolerance / 2.0;
    
    // Using OpenMM's formula for estimating mesh size
    double logTerm = -std::log(reciprocalSpaceError);
    double xterm = alpha * minBoxSize / std::sqrt(logTerm);
    double meshSizeFactor = std::pow(M_PI, 1.0/6.0) * std::pow(6.0 * logTerm, 1.0/3.0) / xterm;
    
    // Determine the mesh dimensions based on the box dimensions
    for (int i = 0; i < 3; i++) {
        double meshScale = box[i] / minBoxSize;
        int size = (int) std::ceil(meshSizeFactor * meshScale);
        
        // Ensure the mesh size is even - following OpenMM's approach
        if (size % 2 != 0)
            size++;
        
        // Make sure the mesh size is at least 4
        if (size < 4)
            size = 4;
        
        // PME mesh size should ideally be a power of 2 for efficient FFT
        int power = 1;
        while (power < size)
            power *= 2;
        
        meshSize[i] = power;
    }
    
    // Log the selected parameters
    platform::log(LogLevel::INFO, "Auto-selected PME parameters: alpha = ", alpha,
                 ", mesh size = [", meshSize[0], ",", meshSize[1], ",", meshSize[2], "]");
    
    // Set parameters - setPMEParameters will initialize B-splines
    setPMEParameters(alpha, meshSize);
    
    // Initialize lookup tables
    pme_params.initializeTables(cutoff_distance);
    
    // Mark PME initialization complete
    pme_params.initialized = true;
}

/**
 * @brief Set PME parameters explicitly - aligned with pme_init in pme.cpp
 * 
 * @param alpha Ewald separation parameter
 * @param meshSize Grid dimensions for PME
 * @param splineOrder B-spline order
 * @param tolerance Precision parameter
 */
void setPMEParameters(double alpha, const int meshSize[3], int splineOrder, double tolerance) {
    pme_params.alpha = alpha;
    
    // Ensure grid size is a power of 2 - crucial for FFT implementation
    for (int i = 0; i < 3; i++) {
        if ((meshSize[i] & (meshSize[i] - 1)) != 0) {
            // If not a power of 2, take the next higher power of 2
            int power = 1;
            while (power < meshSize[i]) {
                power *= 2;
            }
            pme_params.meshSize[i] = power;
            platform::log(LogLevel::WARNING, "PME mesh size must be a power of 2. Adjusting dimension ", 
                         i, " from ", meshSize[i], " to ", pme_params.meshSize[i]);
        } else {
            pme_params.meshSize[i] = meshSize[i];
        }
    }
    
    // Ensure spline order is within reasonable range (usually 3-6)
    if (splineOrder < 3) {
        platform::log(LogLevel::WARNING, "B-spline order less than 3 may lead to poor accuracy. Setting to 3.");
        pme_params.splineOrder = 3;
    } else if (splineOrder > 6) {
        platform::log(LogLevel::WARNING, "B-spline orders > 6 may be computationally expensive. Consider using order 4-6 for optimal performance.");
        pme_params.splineOrder = std::min(splineOrder, 10);  // Limit upper bound to 10
    } else {
        pme_params.splineOrder = splineOrder;
    }
    
    pme_params.tolerance = tolerance;
    
    // Set dielectric constant - usually 1.0
    pme_params.epsilon_r = 1.0;
    
    // Detailed log output
    platform::log(LogLevel::INFO, "PME parameters set: alpha = ", alpha,
                 ", mesh size = [", pme_params.meshSize[0], ",", pme_params.meshSize[1], ",", pme_params.meshSize[2], "]",
                 ", spline order = ", pme_params.splineOrder,
                 ", tolerance = ", tolerance);
    
    // Initialize B-splines - same as pme_init function in pme.cpp
    pme_params.initializeBsplines();
}

/**
 * @brief Calculate pair energy using PME for real space
 */
std::pair<double, double> calcPairEnergyPME(
    double r2, double sigma, double eps, double q1, double q2, 
    const model::MCInfo& info,
    bool is_excluded)
{
    double r = std::sqrt(r2);
    double lj_energy = 0.0;
    double elec_energy = 0.0;
    
    // LJ energy calculation - same as Ewald
    if (r2 < info.cutoff * info.cutoff) {
        double inv_r2 = 1.0 / r2;
        double inv_r6 = inv_r2 * inv_r2 * inv_r2;
        double inv_r12 = inv_r6 * inv_r6;
        double sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
        double sigma12 = sigma6 * sigma6;
        lj_energy = 4.0 * eps * (sigma12 * inv_r12 - sigma6 * inv_r6);
    }
    
    // Electrostatic energy - use PME approximation
    if (!is_excluded && r < pme_params.cutoff) {
        // Real-space contribution for PME
        elec_energy = q1 * q2 * pme_params.ewaldScaleApprox(r);
    }
    
    return {lj_energy, elec_energy};
}

/**
 * @brief Compute B-spline coefficients exactly like pme.cpp
 * 
 * @param fractional Fractional position (0-1)
 * @param order Spline order
 * @param coefficients Output coefficients
 */
void computeBSplineCoefficients(double fractional, int order, std::vector<double>& coefficients) {
    // Ensure coefficients vector is of correct size
    coefficients.resize(order);
    
    // Zero out all coefficients
    for (int i = 0; i < order; i++) {
        coefficients[i] = 0.0;
    }
    
    // Get fractional part
    double dr = fractional;
    
    // Initialize second-order B-spline basis coefficients
    coefficients[0] = 1.0 - dr;
    coefficients[1] = dr;
    
    // Recursively compute B-spline coefficients from third order to order-1 (excluding last step)
    for (int k = 3; k < order; k++) {
        double div = 1.0 / (k - 1.0);
        coefficients[k-1] = div * dr * coefficients[k-2];
        
        for (int i = 1; i < (k-1); i++) {
            coefficients[k-i-1] = div * ((dr+i) * coefficients[k-i-2] + 
                                         (k-i-dr) * coefficients[k-i-1]);
        }
        
        coefficients[0] = div * (1.0-dr) * coefficients[0];
    }
    
    // Last step: special handling for k=order case
    double div = 1.0 / (order - 1);
    coefficients[order-1] = div * dr * coefficients[order-2];
    
    for (int i = 1; i < (order-1); i++) {
        coefficients[order-i-1] = div * ((dr+i) * coefficients[order-i-2] + 
                                        (order-i-dr) * coefficients[order-i-1]);
    }
    coefficients[0] = div * (1.0-dr) * coefficients[0];
    
    // Verify coefficient sum
    double sum = 0.0;
    for (int i = 0; i < order; i++) {
        sum += coefficients[i];
    }
    
    // Warn only if deviation is significant
    if (std::abs(sum - 1.0) > 1e-5) {
        platform::log(LogLevel::WARNING, "Warning: B-spline coefficient sum (" + std::to_string(sum) + ") deviates significantly from 1");
    }
}

/**
 * @brief Spread charges onto the PME grid
 * 
 * @param state MC state
 * @param movement_only Whether to process only moving atoms
 */
void spreadChargesOntoGrid(model::MCState& state, [[maybe_unused]] bool movement_only) {
    // Log the start of processing
    platform::log(LogLevel::DEBUG, "Spreading charges onto PME grid");
    
    // 只在debug_mode启用时执行以下代码
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Spreading charges onto PME grid");
        
        // Calculate total system charge
        double totalCharge = 0.0;
        for (int i = 0; i < state.activeAtomCount; i++) {
            totalCharge += state.atoms[i].charge;
        }
        platform::log(LogLevel::DEBUG, "Total system charge: " + std::to_string(totalCharge));
        
        // Output initial grid values
        platform::log(LogLevel::DEBUG, "Initial values of the first 10 grid points:");
        for (int i = 0; i < 10 && i < static_cast<int>(pme_params.pmeGrid.size()); i++) {
            platform::log(LogLevel::DEBUG, "  Grid point[" + std::to_string(i) + "] = " + std::to_string(pme_params.pmeGrid[i].real()));
        }
        
        // Output charge values
        platform::log(LogLevel::DEBUG, "Charge values of the first 10 atoms:");
        for (int i = 0; i < 10 && i < state.activeAtomCount; i++) {
            platform::log(LogLevel::DEBUG, "  Atom[" + std::to_string(i) + "] charge = " + std::to_string(state.atoms[i].charge));
        }
    }
    
    // Directly access member variables instead of using getter methods
    const auto& atoms = state.atoms;
    // For positions, use atoms to directly access coordinates
    const auto& info = state.info;
    // Change type from double to float to match info.box type
    const float* box = info.box;
    
    // Calculate total system charge
    double totalCharge = 0.0;
    for (int i = 0; i < state.activeAtomCount; ++i) {
        totalCharge += atoms[i].charge;
    }
    
    // Output total system charge - keep this key information
    platform::log(LogLevel::INFO, "Total system charge: ", totalCharge);
    platform::log(LogLevel::DEBUG, "Total system charge: " + std::to_string(totalCharge));
    
    // Check the initial values of the first 10 grid points
    platform::log(LogLevel::DEBUG, "Initial values of the first 10 grid points:");
    for (int i = 0; i < 10 && i < static_cast<int>(pme_params.pmeGrid.size()); i++) {
        platform::log(LogLevel::DEBUG, "  Grid point[" + std::to_string(i) + "] = " + std::to_string(pme_params.pmeGrid[i].real()));
    }
    
    // Print some atom charge values to verify if there are non-zero charges
    platform::log(LogLevel::DEBUG, "Charge values of the first 10 atoms:");
    for (int i = 0; i < 10 && i < state.activeAtomCount; i++) {
        platform::log(LogLevel::DEBUG, "  Atom[" + std::to_string(i) + "] charge = " + std::to_string(atoms[i].charge));
    }
    
    // Calculate reciprocal lattice vectors - ensure consistent with pme.cpp
    double recipBoxVectors[3][3] = {{0}};
    
    // Handle periodic box vectors - simplify for diagonal boxes
    double periodicBoxVectors[3][3] = {
        {box[0], 0.0, 0.0},
        {0.0, box[1], 0.0},
        {0.0, 0.0, box[2]}
    };
    
    // Check if it's a diagonal box - no need to output
    bool isDiagonalBox = true;  // Always true since we force it to be diagonal
    
    if (isDiagonalBox) {
        // Simple calculations for diagonal boxes
        recipBoxVectors[0][0] = 2.0 * M_PI / box[0]; // 2π/a
        recipBoxVectors[1][1] = 2.0 * M_PI / box[1]; // 2π/b 
        recipBoxVectors[2][2] = 2.0 * M_PI / box[2]; // 2π/c
    } else {
        // Non-diagonal boxes require full calculation of reciprocal lattice vectors
        double det = periodicBoxVectors[0][0] * (periodicBoxVectors[1][1] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][1]) -
                     periodicBoxVectors[0][1] * (periodicBoxVectors[1][0] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][0]) +
                     periodicBoxVectors[0][2] * (periodicBoxVectors[1][0] * periodicBoxVectors[2][1] - periodicBoxVectors[1][1] * periodicBoxVectors[2][0]);
        
        // Calculate cross products and reciprocal lattice vectors
        recipBoxVectors[0][0] = 2.0 * M_PI * (periodicBoxVectors[1][1] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][1]) / det;
        recipBoxVectors[0][1] = 2.0 * M_PI * (periodicBoxVectors[0][2] * periodicBoxVectors[2][1] - periodicBoxVectors[0][1] * periodicBoxVectors[2][2]) / det;
        recipBoxVectors[0][2] = 2.0 * M_PI * (periodicBoxVectors[0][1] * periodicBoxVectors[1][2] - periodicBoxVectors[0][2] * periodicBoxVectors[1][1]) / det;
        
        recipBoxVectors[1][0] = 2.0 * M_PI * (periodicBoxVectors[1][2] * periodicBoxVectors[2][0] - periodicBoxVectors[1][0] * periodicBoxVectors[2][2]) / det;
        recipBoxVectors[1][1] = 2.0 * M_PI * (periodicBoxVectors[0][0] * periodicBoxVectors[2][2] - periodicBoxVectors[0][2] * periodicBoxVectors[2][0]) / det;
        recipBoxVectors[1][2] = 2.0 * M_PI * (periodicBoxVectors[0][2] * periodicBoxVectors[1][0] - periodicBoxVectors[0][0] * periodicBoxVectors[1][2]) / det;
        
        recipBoxVectors[2][0] = 2.0 * M_PI * (periodicBoxVectors[1][0] * periodicBoxVectors[2][1] - periodicBoxVectors[1][1] * periodicBoxVectors[2][0]) / det;
        recipBoxVectors[2][1] = 2.0 * M_PI * (periodicBoxVectors[0][1] * periodicBoxVectors[2][0] - periodicBoxVectors[0][0] * periodicBoxVectors[2][1]) / det;
        recipBoxVectors[2][2] = 2.0 * M_PI * (periodicBoxVectors[0][0] * periodicBoxVectors[1][1] - periodicBoxVectors[0][1] * periodicBoxVectors[1][0]) / det;
    }
    
    // Reset grid - ensure all points are initialized to 0
    std::fill(pme_params.pmeGrid.begin(), pme_params.pmeGrid.end(), std::complex<double>(0.0, 0.0));
    
    // Get grid dimensions
    int nx = pme_params.meshSize[0];
    int ny = pme_params.meshSize[1];
    int nz = pme_params.meshSize[2];
    int order = pme_params.splineOrder;
    
    // Create temporary arrays - no resizeTempSplineArrays function
    // Create arrays for B-spline coefficients
    std::vector<std::vector<double>> bsplines_theta(3);
    for (int d = 0; d < 3; d++) {
        bsplines_theta[d].resize(order * state.activeAtomCount, 0.0);
    }
    
    // Create arrays for grid indices and fractional parts
    std::vector<std::vector<int>> gridIndices(state.activeAtomCount, std::vector<int>(3, 0));
    std::vector<std::vector<double>> gridFractions(state.activeAtomCount, std::vector<double>(3, 0.0));
    
    // Print grid index information for every 100 atoms - if there are enough atoms
    for (int i = 0; i < std::min(500, state.activeAtomCount); i += 100) {
        if (i < static_cast<int>(gridIndices.size())) {
            platform::log(LogLevel::DEBUG, "Atom " + std::to_string(i) + " grid index: [" 
                    + std::to_string(gridIndices[i][0]) + ", " 
                    + std::to_string(gridIndices[i][1]) + ", " 
                    + std::to_string(gridIndices[i][2]) + "]");
        }
    }
    
    // Calculate grid indices and fractional offsets for all atoms
    for (int atomIdx = 0; atomIdx < state.activeAtomCount; atomIdx++) {
        // Get position from atoms
        float pos[3] = {atoms[atomIdx].x, atoms[atomIdx].y, atoms[atomIdx].z};
        
        // Convert position to fractional coordinates
        double fractional[3];
        for (int d = 0; d < 3; d++) {
            // Calculate fractional coordinate - fix: ensure consistent with pme.cpp
            // Use reciprocal lattice vectors to calculate fractional coordinates, not simple division
            fractional[d] = 0.0;
            for (int j = 0; j < 3; j++) {
                fractional[d] += pos[j] * recipBoxVectors[j][d] / (2.0 * M_PI);
            }
            
            // Ensure in [0,1) range, handle periodic boundary conditions
            fractional[d] -= floor(fractional[d]);
            // Scale fractional coordinates to grid
            fractional[d] *= pme_params.meshSize[d];
        }
        
        // Calculate grid indices and fractional parts - fix: remove incorrect offset
        for (int d = 0; d < 3; d++) {
            gridFractions[atomIdx][d] = fractional[d] - floor(fractional[d]);
            // Fix: remove incorrect -order/2 offset, consistent with pme.cpp
            gridIndices[atomIdx][d] = static_cast<int>(floor(fractional[d]));
            // Ensure grid indices are within correct range
            if (gridIndices[atomIdx][d] < 0) 
                gridIndices[atomIdx][d] += pme_params.meshSize[d];
        }
    }
    
    // Calculate B-spline coefficients for all atoms
    for (int atomIdx = 0; atomIdx < state.activeAtomCount; atomIdx++) {
        double* thetax = &bsplines_theta[0][atomIdx * order];
        double* thetay = &bsplines_theta[1][atomIdx * order];
        double* thetaz = &bsplines_theta[2][atomIdx * order];
        
        // Calculate B-spline coefficients for each dimension
        std::vector<double> coefficients(order);
        
        // X dimension B-spline
        computeBSplineCoefficients(gridFractions[atomIdx][0], order, coefficients);
        for (int i = 0; i < order; i++) {
            thetax[i] = coefficients[i];
        }
        
        // Y dimension B-spline
        computeBSplineCoefficients(gridFractions[atomIdx][1], order, coefficients);
        for (int i = 0; i < order; i++) {
            thetay[i] = coefficients[i];
        }
        
        // Z dimension B-spline
        computeBSplineCoefficients(gridFractions[atomIdx][2], order, coefficients);
        for (int i = 0; i < order; i++) {
            thetaz[i] = coefficients[i];
        }
    }
    
    // Distribute charges to grid
    double totalGridCharge = 0.0;
    int nonZeroPoints = 0;
    int updatedPoints = 0;
    
    for (int atomIdx = 0; atomIdx < state.activeAtomCount; atomIdx++) {
        double charge = atoms[atomIdx].charge;
        
        // Get grid indices and B-spline coefficients
        int x0index = gridIndices[atomIdx][0];
        int y0index = gridIndices[atomIdx][1];
        int z0index = gridIndices[atomIdx][2];
        
        double* thetax = &bsplines_theta[0][atomIdx * order];
        double* thetay = &bsplines_theta[1][atomIdx * order];
        double* thetaz = &bsplines_theta[2][atomIdx * order];
        
        // Distribute charge to grid - exactly following pme_grid_spread_charge
        for (int ix = 0; ix < order; ix++) {
            int xindex = (x0index + ix) % nx;
            
            for (int iy = 0; iy < order; iy++) {
                int yindex = (y0index + iy) % ny;
                
                for (int iz = 0; iz < order; iz++) {
                    int zindex = (z0index + iz) % nz;
                    
                    // Calculate grid index - ensure exact match with pme.cpp
                    // Original is xindex * ny * nz + yindex * nz + zindex, which is correct
                    int index = xindex * ny * nz + yindex * nz + zindex;
                    
                    // Ensure index doesn't go out of bounds
                    if (index >= 0 && static_cast<size_t>(index) < pme_params.pmeGrid.size()) {
                        // Calculate B-spline weight (product of three directions)
                        double weight = thetax[ix] * thetay[iy] * thetaz[iz];
                        
                        // Distribute charge to grid point - use same method as pme.cpp
                        double chargeContribution = charge * weight;
                        
                        // Key fix: add contribution directly to grid point, consistent with pme.cpp
                        // In pme.cpp: pme->grid[index] += chargeContribution;
                        // This only affects the real part, as chargeContribution is real
                        pme_params.pmeGrid[index] += chargeContribution;
                        
                        // Update statistics
                        totalGridCharge += chargeContribution;
                        if (std::abs(chargeContribution) > 1e-10) {
                            nonZeroPoints++;
                            
                            // Track first 10 updated grid points - modified to use same format as pme.cpp
                            if (updatedPoints < 10) {
                                platform::log(LogLevel::DEBUG, "Updated grid point[" + std::to_string(index) + "]: charge=" + std::to_string(charge) 
                                          + ", weight=" + std::to_string(weight)
                                          + ", contribution=" + std::to_string(chargeContribution));
                                updatedPoints++;
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Output processing progress
    int atomIdx = state.activeAtomCount - 1; // Use index of last processed atom
    if ((atomIdx + 1) % 1000000 == 0) {
        platform::log(LogLevel::DEBUG, "Processed ", atomIdx + 1, " atoms");
    }
    
    // Output charge distribution completion info
    platform::log(LogLevel::INFO, "Charge spreading complete: ", state.activeAtomCount, 
                 " atoms processed, total grid charge = ", totalGridCharge,
                 ", non-zero grid points = ", nonZeroPoints);
    
    platform::log(LogLevel::DEBUG, "Charge spreading complete: " + std::to_string(state.activeAtomCount) 
              + " atoms processed, total grid charge = " + std::to_string(totalGridCharge)
              + ", non-zero grid points = " + std::to_string(nonZeroPoints));
    
    // Analyze grid information
    platform::log(LogLevel::DEBUG, "Grid size: " + std::to_string(pme_params.pmeGrid.size()));
    
    // Add lookup and display for grid points with maximum values
    platform::log(LogLevel::DEBUG, "\n===== Maximum value grid points after charge distribution =====");
    std::vector<std::pair<size_t, double>> topValues;
    for (size_t i = 0; i < pme_params.pmeGrid.size(); i++) {
        double realVal = std::abs(pme_params.pmeGrid[i].real());
        if (realVal > 1e-8) {  // Use larger threshold to find obvious non-zero values
            topValues.push_back({i, realVal});
        }
    }
    
    // Sort by value size
    std::sort(topValues.begin(), topValues.end(), 
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    // Output first 10 maximum value points
    int maxValueCount = 0;
    for (const auto& [idx, val] : topValues) {
        if (maxValueCount >= 10) break;
        // Calculate 3D indices
        int x = (idx / (ny * nz));
        int y = (idx - x * ny * nz) / nz;
        int z = idx - x * ny * nz - y * nz;
        
        platform::log(LogLevel::DEBUG, "Top " + std::to_string(maxValueCount + 1) + ": grid point[" + std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z) 
                  + "] (index=" + std::to_string(idx) + "): " + std::to_string(pme_params.pmeGrid[idx].real()) 
                  + " + " + std::to_string(pme_params.pmeGrid[idx].imag()) + "i, |val| = " + std::to_string(val));
        maxValueCount++;
    }
    
    // Count and output non-zero grid points
    int displayCount = 0;
    for (size_t i = 0; i < pme_params.pmeGrid.size(); i++) {
        if (std::abs(pme_params.pmeGrid[i].real()) > 1e-10) {
            if (displayCount < 5) {
                platform::log(LogLevel::DEBUG, "Non-zero grid point " + std::to_string(displayCount) + ": index=" + std::to_string(i) 
                          + ", value=" + std::to_string(pme_params.pmeGrid[i].real()));
            }
            displayCount++;
        }
    }
    
    // Analyze charge distribution for atom 0
    if (state.activeAtomCount > 0) {
        int atomIdx = 0;
        int x0 = gridIndices[atomIdx][0];
        int y0 = gridIndices[atomIdx][1];
        int z0 = gridIndices[atomIdx][2];
        
        platform::log(LogLevel::DEBUG, "Atom 0: charge=" + std::to_string(atoms[atomIdx].charge) 
                  + ", grid index=[" + std::to_string(x0) + "," + std::to_string(y0) + "," + std::to_string(z0) + "]");
        
        // Analyze grid points around atom 0
        for (int ix = 0; ix < order; ix++) {
            int xindex = (x0 + ix) % nx;
            
            for (int iy = 0; iy < order; iy++) {
                int yindex = (y0 + iy) % ny;
                
                for (int iz = 0; iz < order; iz++) {
                    int zindex = (z0 + iz) % nz;
                    int index = xindex * ny * nz + yindex * nz + zindex;
                    
                    if (index >= 0 && static_cast<size_t>(index) < pme_params.pmeGrid.size()) {
                        platform::log(LogLevel::DEBUG, "  Grid point[" + std::to_string(xindex) + "," + std::to_string(yindex) + "," + std::to_string(zindex) 
                                  + "] (index " + std::to_string(index) + "): " 
                                  + std::to_string(pme_params.pmeGrid[index].real()));
                    }
                }
            }
        }
    }
    
    // Add standard grid point output - for comparison with pme.cpp
    platform::log(LogLevel::DEBUG, "\n===== Standard grid point values comparison after charge distribution =====");
    const int keyIndices[] = {0, 1, nx, ny, nz, nx*ny, nx*nz, ny*nz};
    platform::log(LogLevel::DEBUG, "Total grid points: " + std::to_string(pme_params.pmeGrid.size()));
    for (int i : keyIndices) {
        if (i < static_cast<int>(pme_params.pmeGrid.size())) {
            platform::log(LogLevel::DEBUG, "Grid point[" + std::to_string(i) + "]: " + std::to_string(pme_params.pmeGrid[i].real()) 
                      + " + " + std::to_string(pme_params.pmeGrid[i].imag()) + "i");
        }
    }
    
    // Specific 3D coordinates
    const int keyCoords[][3] = {{0,0,1}, {0,1,0}, {1,0,0}, {1,1,1}, {2,2,2}};
    for (const auto& coord : keyCoords) {
        int idx = ((coord[0] % nx) * ny * nz) + ((coord[1] % ny) * nz) + (coord[2] % nz);
        if (idx < static_cast<int>(pme_params.pmeGrid.size())) {
            platform::log(LogLevel::DEBUG, "Grid point[" + std::to_string(coord[0]) + "," + std::to_string(coord[1]) + "," + std::to_string(coord[2]) 
                      + "] (index=" + std::to_string(idx) + "): " + std::to_string(pme_params.pmeGrid[idx].real()) 
                      + " + " + std::to_string(pme_params.pmeGrid[idx].imag()) + "i");
        }
    }
    
    // After the original processing, add detailed analysis of the first atom's charge distribution
    // Use correct member variable names
    if (state.activeAtomCount > 0) {
        platform::log(LogLevel::DEBUG, "\n===== [energyPME] Detailed analysis of first atom's charge distribution =====");
        
        // Assume the index of the first atom is 0
        int atomIndex = 0;
        
        // Directly get charge and position from atoms
        double atomCharge = state.atoms[atomIndex].charge;
        float posX = state.atoms[atomIndex].x;
        float posY = state.atoms[atomIndex].y;
        float posZ = state.atoms[atomIndex].z;
        
        platform::log(LogLevel::DEBUG, "Atom index: " + std::to_string(atomIndex) + ", charge: " + std::to_string(atomCharge) 
                  + ", position: [" + std::to_string(posX) + "," + std::to_string(posY) + "," + std::to_string(posZ) + "]");
        
        // Calculate the atom's position on the PME grid
        // First calculate fractional coordinates - in [0,1) range
        double fractionPosX = posX / box[0];
        double fractionPosY = posY / box[1];
        double fractionPosZ = posZ / box[2];
        
        // Ensure in [0,1) range
        fractionPosX -= floor(fractionPosX);
        fractionPosY -= floor(fractionPosY);
        fractionPosZ -= floor(fractionPosZ);
        
        // Then calculate grid coordinates
        double gridX = fractionPosX * nx;
        double gridY = fractionPosY * ny;
        double gridZ = fractionPosZ * nz;
        
        // Calculate grid indices and fractional parts
        int gridIX = static_cast<int>(floor(gridX));
        int gridIY = static_cast<int>(floor(gridY));
        int gridIZ = static_cast<int>(floor(gridZ));
        
        double fractionX = gridX - gridIX;
        double fractionY = gridY - gridIY; 
        double fractionZ = gridZ - gridIZ;
        
        // Calculate the starting index for B-spline, considering the order
        int startIX = gridIX - order/2;
        if (startIX < 0) startIX += nx;
        
        int startIY = gridIY - order/2;
        if (startIY < 0) startIY += ny;
        
        int startIZ = gridIZ - order/2;
        if (startIZ < 0) startIZ += nz;
        
        platform::log(LogLevel::DEBUG, "Grid coordinates: [" + std::to_string(gridX) + "," + std::to_string(gridY) + "," + std::to_string(gridZ) + "]");
        platform::log(LogLevel::DEBUG, "Grid integer indices: [" + std::to_string(gridIX) + "," + std::to_string(gridIY) + "," + std::to_string(gridIZ) + "]");
        platform::log(LogLevel::DEBUG, "Grid fractional parts: [" + std::to_string(fractionX) + "," + std::to_string(fractionY) + "," + std::to_string(fractionZ) + "]");
        
        // Calculate and display B-spline coefficients
        std::vector<double> bsCoeffsX(order), bsCoeffsY(order), bsCoeffsZ(order);
        
        // Use existing function to calculate B-spline coefficients
        computeBSplineCoefficients(fractionX, order, bsCoeffsX);
        computeBSplineCoefficients(fractionY, order, bsCoeffsY);
        computeBSplineCoefficients(fractionZ, order, bsCoeffsZ);
        
        platform::log(LogLevel::DEBUG, "X direction B-spline coefficients: " + std::to_string(bsCoeffsX[0]) + " " + std::to_string(bsCoeffsX[1]) + " " + std::to_string(bsCoeffsX[2]) + " " + std::to_string(bsCoeffsX[3]));
        platform::log(LogLevel::DEBUG, "Y direction B-spline coefficients: " + std::to_string(bsCoeffsY[0]) + " " + std::to_string(bsCoeffsY[1]) + " " + std::to_string(bsCoeffsY[2]) + " " + std::to_string(bsCoeffsY[3]));
        platform::log(LogLevel::DEBUG, "Z direction B-spline coefficients: " + std::to_string(bsCoeffsZ[0]) + " " + std::to_string(bsCoeffsZ[1]) + " " + std::to_string(bsCoeffsZ[2]) + " " + std::to_string(bsCoeffsZ[3]));
        
        // Add charge distribution output consistent with pme.cpp
        platform::log(LogLevel::DEBUG, "\nGrid points and their values where charge is distributed:");
        
        // Use different variable names to avoid conflicts
        int gridIndexX = gridIX;
        int gridIndexY = gridIY;
        int gridIndexZ = gridIZ;
        
        // Output charge distribution around the atom
        for (int ix = 0; ix < order; ix++) {
            int xindex = (gridIndexX + ix) % nx;
            
            for (int iy = 0; iy < order; iy++) {
                int yindex = (gridIndexY + iy) % ny;
                
                for (int iz = 0; iz < order; iz++) {
                    int zindex = (gridIndexZ + iz) % nz;
                    int index = xindex * ny * nz + yindex * nz + zindex;
                    
                    if (index >= 0 && static_cast<size_t>(index) < pme_params.pmeGrid.size()) {
                        double weight = bsCoeffsX[ix] * bsCoeffsY[iy] * bsCoeffsZ[iz];
                        double chargeContribution = atomCharge * weight;
                        
                        platform::log(LogLevel::DEBUG, "Grid point[" + std::to_string(xindex) + "," + std::to_string(yindex) + "," + std::to_string(zindex)
                                  + "], index: " + std::to_string(index) 
                                  + ", received charge: " + std::to_string(chargeContribution)
                                  + ", charge coefficient: " + std::to_string(weight));
                    }
                }
            }
        }
    }
    
    // Output grid index and B-spline coefficients for some atoms
    if (platform::verbose_ && platform::log_level_ <= LogLevel::DEBUG && platform::is_debug_mode()) {
        // Output info for a few atoms at DEBUG log level
        for (int i = 0; i < std::min(3, state.activeAtomCount); i++) {
            platform::log(LogLevel::DEBUG, "Atom ", i, " grid index: [", 
                         gridIndices[i][0], ", ", gridIndices[i][1], ", ", gridIndices[i][2], "]");
        }
    }
}

/**
 * @brief Perform forward FFT on the grid
 * 
 * Uses custom FFT implementation
 */
void performFFTForward() {
    // Console output - same as pme.cpp
    platform::log(LogLevel::DEBUG, "Performing forward FFT on PME grid");
    
    // Only execute test code in debug mode
    if (platform::is_debug_mode()) {
        // Calculate non-zero points before FFT
        int nonZeroBeforeFFT = 0;
        for (size_t i = 0; i < pme_params.pmeGrid.size(); i++) {
            if (std::abs(pme_params.pmeGrid[i]) > 1e-10) nonZeroBeforeFFT++;
        }
        platform::log(LogLevel::DEBUG, "Grid before FFT: non-zero points = " + std::to_string(nonZeroBeforeFFT));
    }
    
    // Get grid dimensions
    int nx = pme_params.meshSize[0];
    int ny = pme_params.meshSize[1];
    int nz = pme_params.meshSize[2];
    
    // Use custom FFT implementation
    CustomFFT::fft3D_forward(pme_params.pmeGrid.data(), nx, ny, nz);
    
    // Calculate non-zero points after FFT (moved outside debug check to ensure it's always calculated)
    int nonZeroAfterFFT = 0;
    for (size_t i = 0; i < pme_params.pmeGrid.size(); i++) {
        if (std::abs(pme_params.pmeGrid[i]) > 1e-10) nonZeroAfterFFT++;
    }
    platform::log(LogLevel::INFO, "Grid after FFT: non-zero points = ", nonZeroAfterFFT);
    
    // Only execute detailed test output in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Grid after FFT: non-zero points = " + std::to_string(nonZeroAfterFFT));
        
        // Standard grid point value comparison
        platform::log(LogLevel::DEBUG, "\n===== Standard grid point values comparison after FFT =====");
        platform::log(LogLevel::DEBUG, "Total grid points: " + std::to_string(pme_params.pmeGrid.size()));
        
        // Output first few grid points
        for (int i = 0; i < std::min(5, static_cast<int>(pme_params.pmeGrid.size())); i++) {
            platform::log(LogLevel::DEBUG, "Grid point[" + std::to_string(i) + "]: " + std::to_string(pme_params.pmeGrid[i].real()) 
                        + " + " + std::to_string(pme_params.pmeGrid[i].imag()) + "i, |val|² = " 
                        + std::to_string(std::norm(pme_params.pmeGrid[i])));
        }
        
        // Output a few standard grid coordinates
        const int numCoords = 4;
        const int standardCoords[numCoords][3] = {{0,0,0}, {5,5,5}, {10,10,10}, {20,20,20}};
        for (int i = 0; i < numCoords; i++) {
            const int* coord = standardCoords[i];
            int index = coord[0] * ny * nz + coord[1] * nz + coord[2];
            if (index >= 0 && static_cast<size_t>(index) < pme_params.pmeGrid.size()) {
                platform::log(LogLevel::DEBUG, "Grid point[" + std::to_string(coord[0]) + "," + std::to_string(coord[1]) + "," + std::to_string(coord[2]) 
                            + "] (index=" + std::to_string(index) + "): " + std::to_string(pme_params.pmeGrid[index].real()) 
                            + " + " + std::to_string(pme_params.pmeGrid[index].imag()) + "i, |val|² = " 
                            + std::to_string(std::norm(pme_params.pmeGrid[index])));
            }
        }
        
        // Find max value grid points
        platform::log(LogLevel::DEBUG, "\n===== Maximum value grid points after FFT =====");
        
        // Find top 10 grid values by magnitude
        std::vector<std::pair<int, double>> topValues;
        for (size_t i = 0; i < pme_params.pmeGrid.size(); i++) {
            double val = std::norm(pme_params.pmeGrid[i]);
            if (val > 1e-10) {
                topValues.push_back({static_cast<int>(i), val});
            }
        }
        
        // Sort by magnitude
        std::sort(topValues.begin(), topValues.end(), 
                 [](const auto& a, const auto& b) { return a.second > b.second; });
        
        // Output top values
        int count = 0;
        for (const auto& [idx, val] : topValues) {
            if (count >= 10) break;
            // Calculate 3D indices
            int x = (idx / (ny * nz));
            int y = (idx - x * ny * nz) / nz;
            int z = idx - x * ny * nz - y * nz;
            
            platform::log(LogLevel::DEBUG, "Top " + std::to_string(count + 1) + ": grid point[" + std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z) 
                        + "] (index=" + std::to_string(idx) + "): " + std::to_string(pme_params.pmeGrid[idx].real()) 
                        + " + " + std::to_string(pme_params.pmeGrid[idx].imag()) + "i, |val|² = " + std::to_string(val));
            count++;
        }
    }
    
    // Output grid index and B-spline coefficients for some atoms
    if (platform::verbose_ && platform::log_level_ <= LogLevel::DEBUG && platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "First few grid points after FFT:");
        for (int i = 0; i < 3 && i < static_cast<int>(pme_params.pmeGrid.size()); i++) {
            platform::log(LogLevel::DEBUG, "  Grid point[", i, "] = ", 
                         pme_params.pmeGrid[i].real(), " + ", pme_params.pmeGrid[i].imag(), "i");
        }
    }
}

/**
 * @brief Compute energy from the PME grid after FFT
 * 
 * @param energy Output energy
 */
void computeEnergyFromGrid(double& energy, const double box[3]) {
    // Log beginning of processing
    platform::log(LogLevel::DEBUG, "Computing energy from PME grid");
    
    // Use exactly the same calculation method as pme.cpp
    double volume = box[0] * box[1] * box[2];
    // Use exactly the same constants as pme.cpp
    double one_4pi_eps = 138.935456/pme_params.epsilon_r; // Ensure using the same Coulomb constant as pme.cpp
    double factor = M_PI*M_PI/(pme_params.alpha*pme_params.alpha);
    // Calculate boxfactor: exactly like pme.cpp
    double boxfactor = M_PI * volume;
    
    platform::log(LogLevel::DEBUG, "Computing energy from grid with box = [" + std::to_string(box[0]) + "," + 
                 std::to_string(box[1]) + "," + std::to_string(box[2]) + "], alpha = " + 
                 std::to_string(pme_params.alpha) + ", volume = " + std::to_string(volume));
    
    platform::log(LogLevel::DEBUG, "Energy parameters: one_4pi_eps = " + std::to_string(one_4pi_eps) + 
                 ", factor = " + std::to_string(factor) + ", boxfactor = " + 
                 std::to_string(boxfactor));
                 
    platform::log(LogLevel::DEBUG, "Updating grid data before energy calculation");
    
    // Get grid dimensions
    int nx = pme_params.meshSize[0];
    int ny = pme_params.meshSize[1];
    int nz = pme_params.meshSize[2];
    
    // Calculate reciprocal lattice vectors
    double recipBoxVectors[3][3] = {{0}};
    
    // Process box vectors correctly - maintain complete consistency with pme.cpp
    double periodicBoxVectors[3][3] = {
        {box[0], 0.0, 0.0},
        {0.0, box[1], 0.0},
        {0.0, 0.0, box[2]}
    };
    
    // Diagonal box check
    bool isDiagonalBox = true;  // Force diagonal box
    
    // Calculate reciprocal vectors exactly as in pme.cpp
    if (isDiagonalBox) {
        recipBoxVectors[0][0] = 1.0 / box[0]; 
        recipBoxVectors[1][1] = 1.0 / box[1]; 
        recipBoxVectors[2][2] = 1.0 / box[2]; 
    } else {
        // Non-diagonal box calculation also needs adjustment
        double det = periodicBoxVectors[0][0] * (periodicBoxVectors[1][1] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][1]) -
                     periodicBoxVectors[0][1] * (periodicBoxVectors[1][0] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][0]) +
                     periodicBoxVectors[0][2] * (periodicBoxVectors[1][0] * periodicBoxVectors[2][1] - periodicBoxVectors[1][1] * periodicBoxVectors[2][0]);
        
        // Calculate cross product and reciprocal lattice vectors
        recipBoxVectors[0][0] = (periodicBoxVectors[1][1] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][1]) / det;
        recipBoxVectors[0][1] = (periodicBoxVectors[0][2] * periodicBoxVectors[2][1] - periodicBoxVectors[0][1] * periodicBoxVectors[2][2]) / det;
        recipBoxVectors[0][2] = (periodicBoxVectors[0][1] * periodicBoxVectors[1][2] - periodicBoxVectors[0][2] * periodicBoxVectors[1][1]) / det;
        
        recipBoxVectors[1][0] = (periodicBoxVectors[1][2] * periodicBoxVectors[2][0] - periodicBoxVectors[1][0] * periodicBoxVectors[2][2]) / det;
        recipBoxVectors[1][1] = (periodicBoxVectors[0][0] * periodicBoxVectors[2][2] - periodicBoxVectors[0][2] * periodicBoxVectors[2][0]) / det;
        recipBoxVectors[1][2] = (periodicBoxVectors[0][2] * periodicBoxVectors[1][0] - periodicBoxVectors[0][0] * periodicBoxVectors[1][2]) / det;
        
        recipBoxVectors[2][0] = (periodicBoxVectors[1][0] * periodicBoxVectors[2][1] - periodicBoxVectors[1][1] * periodicBoxVectors[2][0]) / det;
        recipBoxVectors[2][1] = (periodicBoxVectors[0][1] * periodicBoxVectors[2][0] - periodicBoxVectors[0][0] * periodicBoxVectors[2][1]) / det;
        recipBoxVectors[2][2] = (periodicBoxVectors[0][0] * periodicBoxVectors[1][1] - periodicBoxVectors[0][1] * periodicBoxVectors[1][0]) / det;
    }
    
    // Output reciprocal lattice vectors
    platform::log(LogLevel::DEBUG, "Reciprocal lattice vectors:");
    platform::log(LogLevel::DEBUG, "  b1 = [" + std::to_string(recipBoxVectors[0][0]) + ", " + 
                 std::to_string(recipBoxVectors[0][1]) + ", " + std::to_string(recipBoxVectors[0][2]) + "]");
    platform::log(LogLevel::DEBUG, "  b2 = [" + std::to_string(recipBoxVectors[1][0]) + ", " + 
                 std::to_string(recipBoxVectors[1][1]) + ", " + std::to_string(recipBoxVectors[1][2]) + "]");
    platform::log(LogLevel::DEBUG, "  b3 = [" + std::to_string(recipBoxVectors[2][0]) + ", " + 
                 std::to_string(recipBoxVectors[2][1]) + ", " + std::to_string(recipBoxVectors[2][2]) + "]");
    
    // Count original grid data
    int nonZeroGridBefore = 0;
    for (size_t i = 0; i < pme_params.pmeGrid.size(); i++) {
        if (std::norm(pme_params.pmeGrid[i]) > 1e-10) {
            nonZeroGridBefore++;
        }
    }
    platform::log(LogLevel::DEBUG, "Grid before energy calculation: non-zero points = " + std::to_string(nonZeroGridBefore));
    
    // 只在debug_mode启用时执行以下代码
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Computing energy from grid with box = [" + std::to_string(box[0]) + "," +
                    std::to_string(box[1]) + "," + std::to_string(box[2]) + "]");
        
        platform::log(LogLevel::DEBUG, "Energy parameters: one_4pi_eps = " + std::to_string(one_4pi_eps) +
                    ", factor = " + std::to_string(factor));
        
        platform::log(LogLevel::DEBUG, "Updating grid data before energy calculation");
    }
    
    // Initialize energy and point counters
    energy = 0.0;
    int pointsProcessed = 0;
    int significantPoints = 0;
    
    int maxkx = (nx+1)/2;
    int maxky = (ny+1)/2;
    int maxkz = (nz+1)/2;
    
    // 定义监控点数据结构
    struct GridPointData {
        int kx, ky, kz;
        double mx, my, mz;
        double mhx, mhy, mhz;
        double m2;
        double bx, by, bz;
        double denom;
        double eterm;
        std::complex<double> originalValue;
        std::complex<double> updatedValue;
        double struct2;
        double energyContrib;
        bool isSignificant;
    };
    
    // 监控点数据和收集变量，仅在debug模式下使用
    std::vector<GridPointData> monitoredPoints;
    std::vector<GridPointData> significantEnergyPoints;
    std::set<std::tuple<int,int,int>> monitorIndices;
    
    // 仅在debug模式下初始化监控点数据
    if (platform::is_debug_mode()) {
        // Monitor specific points
        monitorIndices = {
            {0,0,1}, {0,1,0}, {1,0,0}, {1,1,1}, {2,2,2}, {5,5,5}, {10,10,10}
        };
    }
    
    // Calculate energy exactly as in pme.cpp
    for (int kx = 0; kx < nx; kx++) {
        // Calculate frequency
        double mx = (kx < maxkx) ? kx : (kx-nx);
        double mhx = mx * recipBoxVectors[0][0];
        
        // Key modification: ensure completely consistent calculation with pme.cpp
        // Exactly reproduce the B-spline moduli application method from pme.cpp
        double bx = boxfactor * pme_params.bsplineModuli[0][kx];
        
        for (int ky = 0; ky < ny; ky++) {
            double my = (ky < maxky) ? ky : (ky-ny);
            // Consistent with pme.cpp, considering non-diagonal terms
            double mhy = mx*recipBoxVectors[1][0] + my*recipBoxVectors[1][1];
            
            // Note: Don't apply boxfactor to by, consistent with pme.cpp
            double by = pme_params.bsplineModuli[1][ky];
            
            for (int kz = 0; kz < nz; kz++) {
                // Skip zero frequency term - for neutral systems
                if (kx == 0 && ky == 0 && kz == 0) {
                    continue;
                }
                
                double mz = (kz < maxkz) ? kz : (kz-nz);
                // Consistent with pme.cpp, considering non-diagonal terms
                double mhz = mx*recipBoxVectors[2][0] + my*recipBoxVectors[2][1] + mz*recipBoxVectors[2][2];
                
                // Get grid data
                int index = kx * ny * nz + ky * nz + kz;
                double d1 = pme_params.pmeGrid[index].real();
                double d2 = pme_params.pmeGrid[index].imag();
                
                // Calculate convolution
                double m2 = mhx * mhx + mhy * mhy + mhz * mhz;
                
                // Don't apply boxfactor to bz, consistent with pme.cpp
                double bz = pme_params.bsplineModuli[2][kz];
                
                // Calculate denom exactly as in pme.cpp
                double denom = m2 * bx * by * bz;
                
                // Improve numerical stability
                if (denom < 1e-10) {
                    denom = 1e-10;
                }
                
                // Ensure energy calculation formula is exactly the same as in pme.cpp
                double eterm = one_4pi_eps * exp(-factor * m2) / denom;
                double struct2 = d1*d1 + d2*d2;
                double energyContrib = eterm * struct2;
                
                // Build monitoring point data - 仅在debug模式下执行
                if (platform::is_debug_mode()) {
                    GridPointData pointData;
                    pointData.kx = kx;
                    pointData.ky = ky;
                    pointData.kz = kz;
                    pointData.mx = mx;
                    pointData.my = my;
                    pointData.mz = mz;
                    pointData.mhx = mhx;
                    pointData.mhy = mhy;
                    pointData.mhz = mhz;
                    pointData.m2 = m2;
                    pointData.bx = bx;
                    pointData.by = by;
                    pointData.bz = bz;
                    pointData.denom = denom;
                    pointData.eterm = eterm;
                    pointData.originalValue = std::complex<double>(d1, d2);
                    pointData.struct2 = struct2;
                    pointData.energyContrib = energyContrib;
                    pointData.isSignificant = (energyContrib > 1e-4);
                    
                    // 检查是否为监控点
                    bool isMonitorPoint = (monitorIndices.find(std::make_tuple(kx, ky, kz)) != monitorIndices.end());
                    if (isMonitorPoint) {
                        monitoredPoints.push_back(pointData);
                    }
                    
                    // 收集能量贡献显著的点
                    if (pointData.isSignificant) {
                        significantEnergyPoints.push_back(pointData);
                    }
                    
                    // 特殊输出点[5,5,5] - 保持与pme.cpp一致的格式
                    if (kx == 5 && ky == 5 && kz == 5) {
                        platform::log(LogLevel::DEBUG, "Grid point[" + std::to_string(kx) + "," + std::to_string(ky) + "," + std::to_string(kz) + "] before processing:");
                        platform::log(LogLevel::DEBUG, "  Grid index = " + std::to_string(index));
                        platform::log(LogLevel::DEBUG, "  mx,my,mz = [" + std::to_string(mx) + "," + std::to_string(my) + "," + std::to_string(mz) + "]");
                        platform::log(LogLevel::DEBUG, "  mhx,mhy,mhz = [" + std::to_string(mhx) + "," + std::to_string(mhy) + "," + std::to_string(mhz) + "]");
                        platform::log(LogLevel::DEBUG, "  m2 = " + std::to_string(m2));
                        platform::log(LogLevel::DEBUG, "  bx,by,bz = [" + std::to_string(bx) + "," + std::to_string(by) + "," + std::to_string(bz) + "]");
                        platform::log(LogLevel::DEBUG, "  boxfactor = " + std::to_string(boxfactor));
                        platform::log(LogLevel::DEBUG, "  B-spline moduli = [" + std::to_string(pme_params.bsplineModuli[0][kx]) + "," +
                                   std::to_string(pme_params.bsplineModuli[1][ky]) + "," + std::to_string(pme_params.bsplineModuli[2][kz]) + "]");
                        platform::log(LogLevel::DEBUG, "  denom = " + std::to_string(denom));
                        platform::log(LogLevel::DEBUG, "  eterm = " + std::to_string(eterm));
                        platform::log(LogLevel::DEBUG, "  one_4pi_eps = " + std::to_string(one_4pi_eps)); 
                        platform::log(LogLevel::DEBUG, "  exp(-factor*m2) = " + std::to_string(exp(-factor*m2)));
                        platform::log(LogLevel::DEBUG, "  Original grid value = " + std::to_string(d1) + " + " + std::to_string(d2) + "i");
                    }
                }
                
                // Update grid value - exactly reproduce pme.cpp method
                std::complex<double> updatedValue(d1 * eterm, d2 * eterm);
                pme_params.pmeGrid[index] = updatedValue;
                
                // 更新的网格值输出 - 仅在debug模式下执行
                if (platform::is_debug_mode() && kx == 5 && ky == 5 && kz == 5) {
                    platform::log(LogLevel::DEBUG, "  Updated grid value = " + std::to_string(updatedValue.real()) + " + " + std::to_string(updatedValue.imag()) + "i");
                    platform::log(LogLevel::DEBUG, "  struct2 = " + std::to_string(struct2));
                    platform::log(LogLevel::DEBUG, "  Energy contribution = " + std::to_string(energyContrib));
                    platform::log(LogLevel::DEBUG, "  Accumulated energy = " + std::to_string(energy));
                    platform::log(LogLevel::DEBUG, "  Significant energy? " + std::string(energyContrib > 1e-8 ? "Yes" : "No"));
                    platform::log(LogLevel::DEBUG, "");
                }
                
                // 保存更新的值到监控点数据 - 仅在debug模式下执行
                if (platform::is_debug_mode() && (monitorIndices.find(std::make_tuple(kx, ky, kz)) != monitorIndices.end())) {
                    for (auto& point : monitoredPoints) {
                        if (point.kx == kx && point.ky == ky && point.kz == kz) {
                            point.updatedValue = updatedValue;
                            break;
                        }
                    }
                }
                
                // Accumulate energy
                energy += energyContrib;
                pointsProcessed++;
                
                if (energyContrib > 1e-8) {
                    significantPoints++;
                }
                
                // Output after point update - only for point [5,5,5]
                if (platform::is_debug_mode() && kx == 5 && ky == 5 && kz == 5) {
                    platform::log(LogLevel::DEBUG, "  Updated grid value = " + std::to_string(updatedValue.real()) + " + " + std::to_string(updatedValue.imag()) + "i");
                    platform::log(LogLevel::DEBUG, "  struct2 = " + std::to_string(struct2)); 
                    platform::log(LogLevel::DEBUG, "  Energy contribution = " + std::to_string(energyContrib));
                    platform::log(LogLevel::DEBUG, "  Accumulated energy = " + std::to_string(energy));
                    platform::log(LogLevel::DEBUG, "  Significant energy? " + std::string(energyContrib > 1e-8 ? "Yes" : "No"));
                    platform::log(LogLevel::DEBUG, "");
                }
                
                // Special debug output - similar to monitoring points in pme.cpp
                if (platform::is_debug_mode() && ((kx == 0 && ky == 0 && kz == 1) || 
                    (kx == 0 && ky == 1 && kz == 0) || 
                    (kx == 1 && ky == 0 && kz == 0) ||
                    (kx == 1 && ky == 1 && kz == 1) ||
                    (kx == 2 && ky == 2 && kz == 2))) {
                    platform::log(LogLevel::DEBUG, "Grid point[" + std::to_string(kx) + "," + std::to_string(ky) + "," + std::to_string(kz) + "] processing:");
                    platform::log(LogLevel::DEBUG, "  mx,my,mz = [" + std::to_string(mx) + "," + std::to_string(my) + "," + std::to_string(mz) + "]");
                    platform::log(LogLevel::DEBUG, "  mhx,mhy,mhz = [" + std::to_string(mhx) + "," + std::to_string(mhy) + "," + std::to_string(mhz) + "]");
                    platform::log(LogLevel::DEBUG, "  m2 = " + std::to_string(m2) + ", eterm = " + std::to_string(eterm));
                    platform::log(LogLevel::DEBUG, "  grid = [" + std::to_string(d1) + "," + std::to_string(d2) + "]");
                    platform::log(LogLevel::DEBUG, "  energy contrib = " + std::to_string(energyContrib));
                }
            }
        }
    }
    
    // Consistent with pme.cpp: multiply by 0.5
    double rawEnergy = energy;
    energy *= 0.5;
    
    // Output important statistics
    platform::log(LogLevel::INFO, "PME reciprocal energy = ", energy);
    
    // Count non-zero updated points
    int nonZeroUpdated = 0;
    for (size_t i = 0; i < pme_params.pmeGrid.size(); i++) {
        if (std::norm(pme_params.pmeGrid[i]) > 1e-10) {
            nonZeroUpdated++;
        }
    }
    
    // 只在debug_mode模式下输出详细统计信息
    if (platform::is_debug_mode()) {
        // Add reciprocal space energy results title and format consistent with pme.cpp
        platform::log(LogLevel::DEBUG, "\n===== [energyPME.cpp] Reciprocal Space Energy Calculation Results =====");
        platform::log(LogLevel::DEBUG, "Total points processed: " + std::to_string(pointsProcessed));
        platform::log(LogLevel::DEBUG, "Significant energy points: " + std::to_string(significantPoints));
        platform::log(LogLevel::DEBUG, "Raw energy sum: " + std::to_string(rawEnergy));
        platform::log(LogLevel::DEBUG, "Final reciprocal space energy: " + std::to_string(energy) + " (multiplied by 0.5)");
        platform::log(LogLevel::DEBUG, "Non-zero grid points: before calculation=" + std::to_string(nonZeroGridBefore) + ", after calculation=" + std::to_string(nonZeroUpdated));
        
        // Output detailed information for all monitoring points
        if (!monitoredPoints.empty()) {
            platform::log(LogLevel::DEBUG, "\nMonitoring points energy contributions:");
            for (const auto& point : monitoredPoints) {
                platform::log(LogLevel::DEBUG, "  [" + std::to_string(point.kx) + "," 
                            + std::to_string(point.ky) + "," + std::to_string(point.kz) + "] = " 
                            + std::to_string(point.energyContrib) + " (significant: " 
                            + std::string(point.isSignificant ? "yes" : "no") + ")");
            }
        }
    }
}

/**
 * @brief Compute reciprocal space energy using PME
 */
double computeReciprocalPME(model::MCState& state, bool movement_only) {
    platform::log(LogLevel::INFO, "Computing PME reciprocal space energy");
    
    // Simplified system information output
    platform::log(LogLevel::INFO, "Box: [", state.info.box[0], ", ", 
                 state.info.box[1], ", ", state.info.box[2], "], Alpha: ", pme_params.alpha);
    
    // Simplified system charge check
    double totalCharge = 0.0;
    for (int i = 0; i < state.activeAtomCount; i++) {
        totalCharge += state.atoms[i].charge;
    }
    
    if (std::abs(totalCharge) > 1e-6) {
        platform::log(LogLevel::WARNING, "System is not neutral! Total charge = ", totalCharge);
    }
    
    // PME grid check
    if (pme_params.pmeGrid.empty()) {
        platform::log(LogLevel::ERROR, "PME grid not initialized!");
        return 0.0;
    }
    
    // Initialize B-spline functions
    if (pme_params.bsplineModuli[0].empty() || 
        pme_params.bsplineModuli[1].empty() || 
        pme_params.bsplineModuli[2].empty()) {
        platform::log(LogLevel::INFO, "Initializing B-splines for PME calculation...");
        pme_params.initializeBsplines();
    }
    
    // Reset grid
    std::fill(pme_params.pmeGrid.begin(), pme_params.pmeGrid.end(), std::complex<double>(0.0, 0.0));
    
    // Execute PME calculation steps
    spreadChargesOntoGrid(state, movement_only);
    performFFTForward();
    
    // Convert box to double array
    double box[3];
    for (int i = 0; i < 3; i++) {
        box[i] = static_cast<double>(state.info.box[i]);
    }
    
    // Update box size in PME parameters, ensure B-splines and energy calculation use the same volume
    pme_params.setBox(box);
    
    // Calculate energy
    double energy = 0.0;
    computeEnergyFromGrid(energy, box);
    
    // Energy already includes Coulomb constant, no need to multiply again
    double reciprocal_energy = energy;
    
    platform::log(LogLevel::INFO, "Reciprocal space energy = ", reciprocal_energy);
    // Add direct output to console
    platform::log(LogLevel::DEBUG, "Reciprocal space energy = " + std::to_string(reciprocal_energy));
    
    return reciprocal_energy;
}

/**
 * @brief Compute self-energy term for PME
 * 
 * @param state MC state
 * @param movement_only Whether to compute only for moving atoms
 * @return double Self-energy
 */
double computeSelfEnergyPME(model::MCState& state, bool movement_only) {
    platform::log(LogLevel::INFO, "Computing self energy with alpha = ", pme_params.alpha);
    
    // Self-energy calculation same as Ewald
    double self_energy = 0.0;
    
    // Track sum of squared charges
    double sum_q2 = 0.0;
    int count = 0;
    
    if(movement_only) {
        for(const auto& movementInfo : state.movementResidues) {
            for(int i = movementInfo.startIndex; 
                i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                if(!state.residues[i].active) continue;
                
                for(int j = state.residues[i].atomStart;
                    j < state.residues[i].atomStart + state.residues[i].atomCount; j++) {
                    if (j < state.activeAtomCount) { // Ensure index is valid
                        double charge = state.atoms[j].charge;
                        double q2 = charge * charge;
                        sum_q2 += q2;
                        
                        if (count < 5) {
                            platform::log(LogLevel::INFO, "Atom ", j, 
                                        " charge = ", charge, 
                                        ", q² = ", q2);
                            count++;
                        }
                    }
                }
            }
        }
    }
    else {
        for(int i = 0; i < state.activeAtomCount; i++) {
            // No longer check active member
            double charge = state.atoms[i].charge;
            double q2 = charge * charge;
            sum_q2 += q2;
            
            if (i < 5) {
                platform::log(LogLevel::INFO, "Atom ", i, 
                            " charge = ", charge, 
                            ", q² = ", q2);
            }
        }
    }
    
    platform::log(LogLevel::INFO, "Sum of q² = ", sum_q2);
    
    // Self-energy formula from pme.cpp: -ONE_4PI_EPS0 * alpha / sqrt(M_PI) * sum_q2
    // Ensure we use exactly the same formula, including the COULOMB constant
    double prefactor = -COULOMB * pme_params.alpha / sqrt(M_PI);
    self_energy = prefactor * sum_q2;
    
    platform::log(LogLevel::INFO, "Self energy prefactor = ", prefactor, 
                 ", resulting self energy = ", self_energy);
    
    // Self energy already includes COULOMB constant, no need to multiply elsewhere
    return self_energy;
}

/**
 * @brief Calculate real-space part of PME
 * 
 * @param state MC state
 * @param movement_only Whether to calculate only for moving residues
 * @param store_in_residues Whether to store energy in residues
 */
void computeRealSpacePME(model::MCState& state, bool movement_only, bool store_in_residues) {
    const auto& box = state.info.box;
    auto& atoms = state.atoms;
    auto& residues = state.residues;
    const float cutoff2 = pme_params.cutoff * pme_params.cutoff;

    // Reset electrostatic energy
    for(auto& residue : residues) {
        if(residue.active) {
            residue.energy_elec = 0.0f;
        }
    }
    
    // Real-space total energy
    double real_space_total = 0.0;
    
    // Add debug information
    int debug_count = 0;
    const int max_debug_pairs = 5;

    // Loop over all residue pairs - maintain existing residue loop structure
    for(int r1 = 0; r1 < state.activeResidueCount; r1++) {
        if(!residues[r1].active) continue;
        if(movement_only) {
            bool in_movement = false;
            for(const auto& movementInfo : state.movementResidues) {
                if(r1 >= movementInfo.startIndex && 
                   r1 < movementInfo.startIndex + movementInfo.activeCount) {
                    in_movement = true;
                    break;
                }
            }
            if(!in_movement) continue;
        }
        
        for(int r2 = r1 + 1; r2 < state.activeResidueCount; r2++) {
            if(!residues[r2].active) continue;
            
            // Loop over atoms in each residue
            for(int i = residues[r1].atomStart; 
                i < residues[r1].atomStart + residues[r1].atomCount; i++) {
                // Ensure atom index is valid
                if(i >= state.activeAtomCount) continue;
                
                for(int j = residues[r2].atomStart; 
                    j < residues[r2].atomStart + residues[r2].atomCount; j++) {
                    // Ensure atom index is valid
                    if(j >= state.activeAtomCount) continue;
                    
                    float dx = atoms[i].x - atoms[j].x;
                    float dy = atoms[i].y - atoms[j].y;
                    float dz = atoms[i].z - atoms[j].z;
                    
                    // Apply periodic boundary conditions
                    dx -= box[0] * round(dx / box[0]);
                    dy -= box[1] * round(dy / box[1]);
                    dz -= box[2] * round(dz / box[2]);
                    
                    float r2 = dx*dx + dy*dy + dz*dz;
                    
                    // Skip pairs beyond cutoff
                    if(r2 > cutoff2) continue;
                    
                    // Compute energy
                    float r = sqrt(r2);
                    float qi = atoms[i].charge;
                    float qj = atoms[j].charge;
                    
                    // Skip neutral atoms
                    if(std::abs(qi) < 1e-6 || std::abs(qj) < 1e-6) continue;
                    
                    // Calculate real space contribution for PME - only erfc part
                    double term = pme_params.erfcApprox(r);
                    double pair_energy = qi * qj * term / r;
                    
                    // Print debug information
                    if (debug_count < max_debug_pairs) {
                        platform::log(LogLevel::INFO, 
                            "Debug energyPME: Atom pair (", i, ",", j, "): ",
                            "r = ", r, " nm, ",
                            "q1*q2 = ", qi * qj, ", ",
                            "erfc term = ", term, ", ",
                            "energy = ", pair_energy,
                            ", with COULOMB = ", COULOMB * pair_energy, " kJ/mol");
                        debug_count++;
                    }

                    // Accumulate to total energy
                    real_space_total += pair_energy;
                    
                    // Decide how to store energy based on parameters
                    if (store_in_residues) {
                        // Each residue gets half of the pair interaction energy
                        residues[r1].energy_elec += pair_energy / 2.0f;
                        residues[r2].energy_elec += pair_energy / 2.0f;
                    }
                }
            }
        }
    }
    
    // Store total real-space energy (not yet multiplied by COULOMB)
    state.ewald_energy.real_space = real_space_total;
}

/**
 * @brief Calculate system energy using PME method
 */
void computeSystemEnergyPME(model::MCState& state) {
    if (!pme_params.initialized) {
        throw std::runtime_error("PME parameters not initialized. Call initializePMEParameters() first.");
    }
    
    // Calculate energy for the entire system
    // 1. First calculate real space part - needs to be multiplied by COULOMB factor
    computeRealSpacePME(state, false, true);
    
    // 2. Then calculate reciprocal space part - computeReciprocalPME already includes COULOMB factor
    state.ewald_energy.reciprocal = computeReciprocalPME(state, false);
    
    // 3. Finally calculate self energy part - computeSelfEnergyPME already includes COULOMB factor
    state.ewald_energy.self = computeSelfEnergyPME(state, false);
    
    // Only multiply real space energy by COULOMB coefficient
    state.ewald_energy.real_space *= COULOMB;
    
    // Calculate total energy
    state.ewald_energy.total = state.ewald_energy.real_space + 
                             state.ewald_energy.reciprocal + 
                             state.ewald_energy.self;
    
    platform::log(LogLevel::INFO, "PME system energy components: real_space=", state.ewald_energy.real_space,
                 " reciprocal=", state.ewald_energy.reciprocal,
                 " self=", state.ewald_energy.self,
                 " total=", state.ewald_energy.total);
}

/**
 * @brief Calculate movement energy using PME method
 */
void computeMovementEnergyPME(model::MCState& state) {
    if (!pme_params.initialized) {
        throw std::runtime_error("PME parameters not initialized. Call initializePMEParameters() first.");
    }
    
    // Only calculate for residues that moved
    computeRealSpacePME(state, true, true);
    state.ewald_energy.reciprocal = computeReciprocalPME(state, true);
    state.ewald_energy.self = computeSelfEnergyPME(state, true);
    
    // Apply Coulomb factor only to real-space component
    // Note: computeReciprocalPME and computeSelfEnergyPME already include COULOMB factor
    state.ewald_energy.real_space *= COULOMB;
    
    // Total energy is the sum of all components
    state.ewald_energy.total = state.ewald_energy.real_space + 
                             state.ewald_energy.reciprocal + 
                             state.ewald_energy.self;
    
    platform::log(LogLevel::INFO, "PME movement energy components: real_space=", state.ewald_energy.real_space,
                 " reciprocal=", state.ewald_energy.reciprocal,
                 " self=", state.ewald_energy.self,
                 " total=", state.ewald_energy.total);
}

/**
 * @brief Perform backward FFT on the grid
 * 
 * Uses custom FFT implementation, matches pme.cpp implementation
 */
void performFFTBackward() {
    // Get grid dimensions
    int nx = pme_params.meshSize[0];
    int ny = pme_params.meshSize[1];
    int nz = pme_params.meshSize[2];
    
    // Only execute test code in debug mode
    if (platform::is_debug_mode()) {
        // Calculate non-zero points before backward FFT
        int nonZeroBeforeFFT = 0;
        for (size_t i = 0; i < pme_params.pmeGrid.size(); i++) {
            if (std::abs(pme_params.pmeGrid[i]) > 1e-10) nonZeroBeforeFFT++;
        }
        platform::log(LogLevel::DEBUG, "Grid before backward FFT: non-zero points = " + std::to_string(nonZeroBeforeFFT));
    }
    
    // Use custom FFT implementation
    CustomFFT::fft3D_backward(pme_params.pmeGrid.data(), nx, ny, nz);
    
    // Only execute test code in debug mode
    if (platform::is_debug_mode()) {
        // Calculate non-zero points after backward FFT
        int nonZeroAfterFFT = 0;
        for (size_t i = 0; i < pme_params.pmeGrid.size(); i++) {
            if (std::abs(pme_params.pmeGrid[i]) > 1e-10) nonZeroAfterFFT++;
        }
        platform::log(LogLevel::DEBUG, "Grid after backward FFT: non-zero points = " + std::to_string(nonZeroAfterFFT));
        
        // Output grid values after backward FFT for comparison
        platform::log(LogLevel::DEBUG, "\nFirst few grid points after backward FFT:");
        for (int i = 0; i < 5 && i < static_cast<int>(pme_params.pmeGrid.size()); i++) {
            platform::log(LogLevel::DEBUG, "  Grid point[" + std::to_string(i) + "] = " 
                        + std::to_string(pme_params.pmeGrid[i].real()) + " + " 
                        + std::to_string(pme_params.pmeGrid[i].imag()) + "i");
        }
    }
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc





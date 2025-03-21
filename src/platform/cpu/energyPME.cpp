#include "energyPME.hpp"
#include "platform/platform.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <complex>
#include <chrono>

// Optional: Include library for FFT if needed
// #include <fftw3.h>

namespace {
// 定义PI2
const double PI2 = 6.28318530717958647692;

// 在C++中，我们使用std::complex<double>替代C的complex类型
using cmplx = std::complex<double>;

namespace CustomFFT {

// 添加静态权重向量
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

// 恢复log2_power_of_2函数
// Calculate log2(n), n must be a power of 2
int log2_power_of_2(int n) {
    int k = 0;
    while (n > 1) {
        n >>= 1;
        k++;
    }
    return k;
}

// 修复padded_fft函数，使用log2_power_of_2计算k值
void padded_fft(cmplx* data, int actual_size, bool inverse) {
    // Check if it's a power of two
    bool is_power_of_two = (actual_size & (actual_size - 1)) == 0;
    
    if (!is_power_of_two) {
        throw std::runtime_error("FFT size must be a power of 2");
    }
    
    // Calculate k where 2^k = actual_size
    int k = log2_power_of_2(actual_size);
    
    // 确保权重向量足够大
    if (weights.size() < static_cast<size_t>(actual_size)) {
        weights.resize(actual_size);
        tfft_init(k, weights.data());  // 使用k而不是actual_size
    }
    
    // 执行FFT
    if (inverse) {
        tfft_ifft(k, data, weights.data());  // 使用k而不是actual_size
    } else {
        tfft_fft(k, data, weights.data());  // 使用k而不是actual_size
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

// 添加FFT相关成员变量
std::vector<std::complex<double>> fft_weights;

// 全局变量用于保存FFT后的网格数据
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
 * @brief Initialize B-splines for PME - 完全按照pme.cpp中的pme_calculate_bsplines_moduli实现
 */
void PMEParams::initializeBsplines() {
    // 初始化B-spline模数 - 完全基于pme.cpp中的pme_calculate_bsplines_moduli实现
    platform::log(LogLevel::INFO, "Initializing B-splines with order = ", splineOrder, 
                 " and mesh size = [", meshSize[0], ",", meshSize[1], ",", meshSize[2], "]");
    
    // 确保样条阶数至少为2
    if (splineOrder < 2) {
        platform::log(LogLevel::WARNING, "B-spline order must be at least 2, setting to 2");
        splineOrder = 2;
    }
    
    // 计算体积和boxfactor - 与pme.cpp中的正交盒子实现保持一致
    double boxVolume = box[0] * box[1] * box[2];
    
    // 如果体积为零（未设置盒子），使用单位体积
    if (boxVolume < 1e-10) {
        platform::log(LogLevel::WARNING, "Box volume near zero, using unit volume for boxfactor");
        boxVolume = 1.0;
    }
    
    double boxfactor = M_PI * boxVolume;
    platform::log(LogLevel::INFO, "Initializing B-splines with boxfactor = ", boxfactor, 
                 " (box volume = ", boxVolume, ")");
    
    // 初始化 bsplineModuli 数组 - 精确对齐pme.cpp
    for (int dim = 0; dim < 3; dim++) {
        int size = meshSize[dim];
        bsplineModuli[dim].resize(size);
        
        // 对于每个网格点，计算B样条模数
        for (int i = 0; i < size; i++) {
            // 计算m值 - 与pme.cpp一致
            int m = (i < size/2) ? i : (size - i);
            
            // 零频率点特殊处理 - 与pme.cpp一致
            if (m == 0) {
                bsplineModuli[dim][i] = 1.0;
                // 不再对第一个维度预先应用boxfactor - 将在能量计算时应用，与pme.cpp保持一致
                // if (dim == 0) {
                //     bsplineModuli[dim][i] *= boxfactor;
                // }
                continue;
            }
            
            // 计算B样条模数 (sin(πm/N)/(πm/N))^p - 与pme.cpp完全一致
            double arg = M_PI * m / size;
            // 对于非常小的角度，使用泰勒展开近似
            double sinc = (arg < 1e-7) ? 1.0 - arg*arg/6.0 : sin(arg) / arg;
            
            // 计算sinc函数的p次方 - 与pme.cpp一致
            bsplineModuli[dim][i] = pow(sinc, splineOrder);
            
            // 检查是否需要近似0值 - 保持数值稳定性
            if (bsplineModuli[dim][i] < 1e-10) {
                bsplineModuli[dim][i] = 1e-10;  // 避免除零
            }
            
            // 计算B样条模值的倒数平方 - 与pme.cpp一致
            bsplineModuli[dim][i] = 1.0 / (bsplineModuli[dim][i] * bsplineModuli[dim][i]);
            
            // 不再对第一个维度预先应用boxfactor - 将在能量计算时应用，与pme.cpp保持一致
            // if (dim == 0) {
            //     bsplineModuli[dim][i] *= boxfactor;
            // }
        }
    }
    
    // 输出B样条模数和统计信息
    double maxModuli[3] = {0.0, 0.0, 0.0};
    double minModuli[3] = {std::numeric_limits<double>::max(), 
                          std::numeric_limits<double>::max(), 
                          std::numeric_limits<double>::max()};
    
    for (int dim = 0; dim < 3; dim++) {
        for (size_t i = 0; i < bsplineModuli[dim].size(); i++) {
            maxModuli[dim] = std::max(maxModuli[dim], bsplineModuli[dim][i]);
            minModuli[dim] = std::min(minModuli[dim], bsplineModuli[dim][i]);
        }
        std::cout << "Dimension " << dim << " B-spline moduli range: [" 
                 << minModuli[dim] << ", " << maxModuli[dim] << "]" << std::endl;
    }
    
    // 分配PME网格
    int totalGridPoints = meshSize[0] * meshSize[1] * meshSize[2];
    pmeGrid.resize(totalGridPoints);
    
    // 清零网格
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
    
    // 新增：限制alpha参数的变化范围，减少自能量随误差容忍度的剧烈波动
    // 根据经验，使用适度的alpha范围 (1.0 - 3.0)/cutoff
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
    
    // 设置参数 - setPMEParameters会初始化B样条
    setPMEParameters(alpha, meshSize);
    
    // 初始化查找表
    pme_params.initializeTables(cutoff_distance);
    
    // 标记PME初始化完成
    pme_params.initialized = true;
}

/**
 * @brief Set PME parameters explicitly - 与pme.cpp中的pme_init对齐
 * 
 * @param alpha Ewald separation parameter
 * @param meshSize Grid dimensions for PME
 * @param splineOrder B-spline order
 * @param tolerance Precision parameter
 */
void setPMEParameters(double alpha, const int meshSize[3], int splineOrder, double tolerance) {
    pme_params.alpha = alpha;
    
    // 确保网格尺寸是2的幂 - 这对FFT实现至关重要
    for (int i = 0; i < 3; i++) {
        if ((meshSize[i] & (meshSize[i] - 1)) != 0) {
            // 如果不是2的幂，取下一个更大的2的幂
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
    
    // 确保样条阶数在合理范围内(通常3-6)
    if (splineOrder < 3) {
        platform::log(LogLevel::WARNING, "B-spline order less than 3 may lead to poor accuracy. Setting to 3.");
        pme_params.splineOrder = 3;
    } else if (splineOrder > 6) {
        platform::log(LogLevel::WARNING, "B-spline orders > 6 may be computationally expensive. Consider using order 4-6 for optimal performance.");
        pme_params.splineOrder = std::min(splineOrder, 10);  // 限制上限为10
    } else {
        pme_params.splineOrder = splineOrder;
    }
    
    pme_params.tolerance = tolerance;
    
    // 设置介电常数 - 通常为1.0
    pme_params.epsilon_r = 1.0;
    
    // 详细日志输出
    platform::log(LogLevel::INFO, "PME parameters set: alpha = ", alpha,
                 ", mesh size = [", pme_params.meshSize[0], ",", pme_params.meshSize[1], ",", pme_params.meshSize[2], "]",
                 ", spline order = ", pme_params.splineOrder,
                 ", tolerance = ", tolerance);
    
    // 初始化B样条 - 与pme.cpp中的pme_init函数一致
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
    // 确保coefficients大小正确
    coefficients.resize(order);
    
    // 清零所有系数
    for (int i = 0; i < order; i++) {
        coefficients[i] = 0.0;
    }
    
    // 获取小数部分
    double dr = fractional;
    
    // 初始化二阶B样条基本系数
    coefficients[0] = 1.0 - dr;
    coefficients[1] = dr;
    
    // 递归计算B样条系数，从三阶到order阶(不包括最后一步)
    for (int k = 3; k < order; k++) {
        double div = 1.0 / (k - 1.0);
        coefficients[k-1] = div * dr * coefficients[k-2];
        
        for (int i = 1; i < (k-1); i++) {
            coefficients[k-i-1] = div * ((dr+i) * coefficients[k-i-2] + 
                                         (k-i-dr) * coefficients[k-i-1]);
        }
        
        coefficients[0] = div * (1.0-dr) * coefficients[0];
    }
    
    // 最后一步：特殊处理k=order的情况
    double div = 1.0 / (order - 1);
    coefficients[order-1] = div * dr * coefficients[order-2];
    
    for (int i = 1; i < (order-1); i++) {
        coefficients[order-i-1] = div * ((dr+i) * coefficients[order-i-2] + 
                                        (order-i-dr) * coefficients[order-i-1]);
    }
    coefficients[0] = div * (1.0-dr) * coefficients[0];
    
    // 验证系数总和
    double sum = 0.0;
    for (int i = 0; i < order; i++) {
        sum += coefficients[i];
    }
    
    // 只在偏差大时警告
    if (std::abs(sum - 1.0) > 1e-5) {
        std::cerr << "警告: B样条系数总和 (" << sum << ") 与1相差较大" << std::endl;
    }
}

/**
 * @brief Spread charges onto the PME grid
 * 
 * @param state MC state
 * @param movement_only Whether to process only moving atoms
 */
void spreadChargesOntoGrid(model::MCState& state, [[maybe_unused]] bool movement_only) {
    // 记录开始处理
    platform::log(LogLevel::DEBUG, "Spreading charges onto PME grid");
    
    // 添加控制台输出，与pme.cpp保持一致
    std::cout << "Spreading charges onto PME grid" << std::endl;
    
    // 直接访问成员变量而不是使用getter方法
    const auto& atoms = state.atoms;
    // 对于positions，使用atoms直接访问坐标
    const auto& info = state.info;
    // 修改类型从double到float以匹配info.box的类型
    const float* box = info.box;
    
    // 计算总系统电荷
    double totalCharge = 0.0;
    for (int i = 0; i < state.activeAtomCount; ++i) {
        totalCharge += atoms[i].charge;
    }
    
    // 输出总系统电荷 - 保留这个关键信息
    platform::log(LogLevel::INFO, "Total system charge: ", totalCharge);
    std::cout << "Total system charge: " << totalCharge << std::endl;
    
    // 检查前10个网格点的初始值
    std::cout << "初始10个网格点的值:" << std::endl;
    for (int i = 0; i < 10 && i < static_cast<int>(pme_params.pmeGrid.size()); i++) {
        std::cout << "  网格点[" << i << "] = " << pme_params.pmeGrid[i].real() << std::endl;
    }
    
    // 打印一些原子的电荷值，验证是否有非零电荷
    std::cout << "前10个原子的电荷值:" << std::endl;
    for (int i = 0; i < 10 && i < state.activeAtomCount; i++) {
        std::cout << "  原子[" << i << "] 电荷 = " << atoms[i].charge << std::endl;
    }
    
    // 计算倒易晶格矢量 - 确保与pme.cpp一致
    double recipBoxVectors[3][3] = {{0}};
    
    // 正确处理盒子向量 - 对角盒子简化处理
    double periodicBoxVectors[3][3] = {
        {box[0], 0.0, 0.0},
        {0.0, box[1], 0.0},
        {0.0, 0.0, box[2]}
    };
    
    // 检查是否是对角盒子 - 不需要输出
    bool isDiagonalBox = true;  // 由于我们强制设置为对角盒子，所以始终为true
    
    if (isDiagonalBox) {
        // 对角盒子的简单计算
        recipBoxVectors[0][0] = 2.0 * M_PI / box[0]; // 2π/a
        recipBoxVectors[1][1] = 2.0 * M_PI / box[1]; // 2π/b 
        recipBoxVectors[2][2] = 2.0 * M_PI / box[2]; // 2π/c
    } else {
        // 非对角盒子需要完全计算倒格矢向量
        double det = periodicBoxVectors[0][0] * (periodicBoxVectors[1][1] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][1]) -
                     periodicBoxVectors[0][1] * (periodicBoxVectors[1][0] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][0]) +
                     periodicBoxVectors[0][2] * (periodicBoxVectors[1][0] * periodicBoxVectors[2][1] - periodicBoxVectors[1][1] * periodicBoxVectors[2][0]);
        
        // 计算叉积和倒易晶格向量
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
    
    // 重置网格 - 确保所有点初始化为0
    std::fill(pme_params.pmeGrid.begin(), pme_params.pmeGrid.end(), std::complex<double>(0.0, 0.0));
    
    // 获取网格维度
    int nx = pme_params.meshSize[0];
    int ny = pme_params.meshSize[1];
    int nz = pme_params.meshSize[2];
    int order = pme_params.splineOrder;
    
    // 创建临时数组 - 因为没有resizeTempSplineArrays函数
    // 创建用于B样条系数的临时数组
    std::vector<std::vector<double>> bsplines_theta(3);
    for (int d = 0; d < 3; d++) {
        bsplines_theta[d].resize(order * state.activeAtomCount, 0.0);
    }
    
    // 创建用于网格索引和分数部分的临时数组
    std::vector<std::vector<int>> gridIndices(state.activeAtomCount, std::vector<int>(3, 0));
    std::vector<std::vector<double>> gridFractions(state.activeAtomCount, std::vector<double>(3, 0.0));
    
    // 打印每100个原子的网格索引信息 - 如果有足够多的原子
    for (int i = 0; i < std::min(500, state.activeAtomCount); i += 100) {
        if (i < static_cast<int>(gridIndices.size())) {
            std::cout << "原子 " << i << " 的网格索引: [" 
                    << gridIndices[i][0] << ", " 
                    << gridIndices[i][1] << ", " 
                    << gridIndices[i][2] << "]" << std::endl;
        }
    }
    
    // 计算所有原子的网格索引和分数偏移
    for (int atomIdx = 0; atomIdx < state.activeAtomCount; atomIdx++) {
        // 从atoms直接获取位置
        float pos[3] = {atoms[atomIdx].x, atoms[atomIdx].y, atoms[atomIdx].z};
        
        // 将位置转换为分数坐标
        double fractional[3];
        for (int d = 0; d < 3; d++) {
            // 计算分数坐标 - 修复：确保与pme.cpp一致
            // 使用倒易格矢量来计算分数坐标，而不是简单的除法
            fractional[d] = 0.0;
            for (int j = 0; j < 3; j++) {
                fractional[d] += pos[j] * recipBoxVectors[j][d] / (2.0 * M_PI);
            }
            
            // 确保在[0,1)范围内，处理周期性边界条件
            fractional[d] -= floor(fractional[d]);
            // 将分数坐标缩放到网格上
            fractional[d] *= pme_params.meshSize[d];
        }
        
        // 计算网格索引和分数部分 - 修复：移除错误的offset
        for (int d = 0; d < 3; d++) {
            gridFractions[atomIdx][d] = fractional[d] - floor(fractional[d]);
            // 修复：移除错误的-order/2偏移，与pme.cpp保持一致
            gridIndices[atomIdx][d] = static_cast<int>(floor(fractional[d]));
            // 确保网格索引在正确范围内
            if (gridIndices[atomIdx][d] < 0) 
                gridIndices[atomIdx][d] += pme_params.meshSize[d];
        }
    }
    
    // 计算所有原子的B样条系数
    for (int atomIdx = 0; atomIdx < state.activeAtomCount; atomIdx++) {
        double* thetax = &bsplines_theta[0][atomIdx * order];
        double* thetay = &bsplines_theta[1][atomIdx * order];
        double* thetaz = &bsplines_theta[2][atomIdx * order];
        
        // 计算三个维度的B样条系数
        std::vector<double> coefficients(order);
        
        // X维度B样条
        computeBSplineCoefficients(gridFractions[atomIdx][0], order, coefficients);
        for (int i = 0; i < order; i++) {
            thetax[i] = coefficients[i];
        }
        
        // Y维度B样条
        computeBSplineCoefficients(gridFractions[atomIdx][1], order, coefficients);
        for (int i = 0; i < order; i++) {
            thetay[i] = coefficients[i];
        }
        
        // Z维度B样条
        computeBSplineCoefficients(gridFractions[atomIdx][2], order, coefficients);
        for (int i = 0; i < order; i++) {
            thetaz[i] = coefficients[i];
        }
    }
    
    // 将电荷分布到网格上
    double totalGridCharge = 0.0;
    int nonZeroPoints = 0;
    int updatedPoints = 0;
    
    for (int atomIdx = 0; atomIdx < state.activeAtomCount; atomIdx++) {
        double charge = atoms[atomIdx].charge;
        
        // 获取网格索引和B样条系数
        int x0index = gridIndices[atomIdx][0];
        int y0index = gridIndices[atomIdx][1];
        int z0index = gridIndices[atomIdx][2];
        
        double* thetax = &bsplines_theta[0][atomIdx * order];
        double* thetay = &bsplines_theta[1][atomIdx * order];
        double* thetaz = &bsplines_theta[2][atomIdx * order];
        
        // 将电荷分布到网格上 - 完全按照pme_grid_spread_charge
        for (int ix = 0; ix < order; ix++) {
            int xindex = (x0index + ix) % nx;
            
            for (int iy = 0; iy < order; iy++) {
                int yindex = (y0index + iy) % ny;
                
                for (int iz = 0; iz < order; iz++) {
                    int zindex = (z0index + iz) % nz;
                    
                    // 计算网格索引 - 确保与pme.cpp完全一致
                    // 原始是xindex * ny * nz + yindex * nz + zindex，这个公式是正确的
                    int index = xindex * ny * nz + yindex * nz + zindex;
                    
                    // 确保索引不越界
                    if (index >= 0 && static_cast<size_t>(index) < pme_params.pmeGrid.size()) {
                        // 计算B样条权重（三个方向的乘积）
                        double weight = thetax[ix] * thetay[iy] * thetaz[iz];
                        
                        // 将电荷分布到网格点 - 使用与pme.cpp相同的方式
                        double chargeContribution = charge * weight;
                        
                        // 关键修复：直接将贡献添加到网格点，与pme.cpp一致
                        // 在pme.cpp中使用的是: pme->grid[index] += chargeContribution;
                        // 这只影响实部，因为chargeContribution是real
                        pme_params.pmeGrid[index] += chargeContribution;
                        
                        // 更新统计信息
                        totalGridCharge += chargeContribution;
                        if (std::abs(chargeContribution) > 1e-10) {
                            nonZeroPoints++;
                            
                            // 追踪前10个更新的网格点 - 修改为使用与pme.cpp相同的格式
                            if (updatedPoints < 10) {
                                std::cout << "更新网格点[" << index << "]: 电荷=" << charge 
                                          << ", 权重=" << weight
                                          << ", 贡献=" << chargeContribution << std::endl;
                                updatedPoints++;
                            }
                        }
                    }
                }
            }
        }
    }
    
    // 输出处理进度
    int atomIdx = state.activeAtomCount - 1; // 使用最后一个处理的原子索引
    if ((atomIdx + 1) % 1000000 == 0) {
        platform::log(LogLevel::DEBUG, "Processed ", atomIdx + 1, " atoms");
    }
    
    // 输出电荷分布完成信息
    platform::log(LogLevel::INFO, "Charge spreading complete: ", state.activeAtomCount, 
                 " atoms processed, total grid charge = ", totalGridCharge,
                 ", non-zero grid points = ", nonZeroPoints);
    
    std::cout << "Charge spreading complete: " << state.activeAtomCount 
              << " atoms processed, total grid charge = " << totalGridCharge
              << ", non-zero grid points = " << nonZeroPoints << std::endl;
    
    // 分析网格信息
    std::cout << "网格大小: " << pme_params.pmeGrid.size() << std::endl;
    
    // 统计并输出非零网格点
    int displayCount = 0;
    for (size_t i = 0; i < pme_params.pmeGrid.size(); i++) {
        if (std::abs(pme_params.pmeGrid[i].real()) > 1e-10) {
            if (displayCount < 5) {
                std::cout << "非零网格点 " << displayCount << ": 索引=" << i 
                          << ", 值=" << pme_params.pmeGrid[i].real() << std::endl;
            }
            displayCount++;
        }
    }
    
    // 分析原子0的电荷分布
    if (state.activeAtomCount > 0) {
        int atomIdx = 0;
        int x0 = gridIndices[atomIdx][0];
        int y0 = gridIndices[atomIdx][1];
        int z0 = gridIndices[atomIdx][2];
        
        std::cout << "原子0: 电荷=" << atoms[atomIdx].charge 
                  << ", 网格索引=[" << x0 << "," << y0 << "," << z0 << "]" << std::endl;
        
        // 分析原子0周围的网格点
        for (int ix = 0; ix < order; ix++) {
            int xindex = (x0 + ix) % nx;
            
            for (int iy = 0; iy < order; iy++) {
                int yindex = (y0 + iy) % ny;
                
                for (int iz = 0; iz < order; iz++) {
                    int zindex = (z0 + iz) % nz;
                    int index = xindex * ny * nz + yindex * nz + zindex;
                    
                    if (index >= 0 && static_cast<size_t>(index) < pme_params.pmeGrid.size()) {
                        std::cout << "  网格点[" << xindex << "," << yindex << "," << zindex 
                                  << "] (索引 " << index << "): " 
                                  << pme_params.pmeGrid[index].real() << std::endl;
                    }
                }
            }
        }
    }
    
    // 添加标准位置格点输出 - 用于与pme.cpp比较
    std::cout << "\n===== 电荷分布后的标准格点值比较 =====\n";
    const int keyIndices[] = {0, 1, nx, ny, nz, nx*ny, nx*nz, ny*nz};
    std::cout << "格点总数: " << pme_params.pmeGrid.size() << std::endl;
    for (int i : keyIndices) {
        if (i < static_cast<int>(pme_params.pmeGrid.size())) {
            std::cout << "格点[" << i << "]: " << pme_params.pmeGrid[i].real() 
                      << " + " << pme_params.pmeGrid[i].imag() << "i" << std::endl;
        }
    }
    
    // 特定的三维坐标
    const int keyCoords[][3] = {{0,0,1}, {0,1,0}, {1,0,0}, {1,1,1}, {2,2,2}};
    for (const auto& coord : keyCoords) {
        int idx = ((coord[0] % nx) * ny * nz) + ((coord[1] % ny) * nz) + (coord[2] % nz);
        if (idx < static_cast<int>(pme_params.pmeGrid.size())) {
            std::cout << "格点[" << coord[0] << "," << coord[1] << "," << coord[2] 
                      << "] (索引=" << idx << "): " << pme_params.pmeGrid[idx].real() 
                      << " + " << pme_params.pmeGrid[idx].imag() << "i" << std::endl;
        }
    }
    
    // 在原本的处理之后，添加详细的第一个原子电荷分布分析
    // 使用正确的成员变量名称
    if (state.activeAtomCount > 0) {
        std::cout << "\n===== [energyPME] 第一个原子电荷分布分析 =====\n";
        
        // 假设第一个原子的索引是0
        int atomIndex = 0;
        
        // 直接从atoms中获取电荷和位置
        double atomCharge = state.atoms[atomIndex].charge;
        float posX = state.atoms[atomIndex].x;
        float posY = state.atoms[atomIndex].y;
        float posZ = state.atoms[atomIndex].z;
        
        std::cout << "原子索引: " << atomIndex << ", 电荷: " << atomCharge 
                  << ", 位置: [" << posX << "," << posY << "," << posZ << "]\n";
        
        // 正确计算原子在PME网格中的位置
        // 首先计算分数坐标 - 在[0,1)范围内
        double fractionPosX = posX / box[0];
        double fractionPosY = posY / box[1];
        double fractionPosZ = posZ / box[2];
        
        // 确保在[0,1)范围内
        fractionPosX -= floor(fractionPosX);
        fractionPosY -= floor(fractionPosY);
        fractionPosZ -= floor(fractionPosZ);
        
        // 然后计算网格坐标
        double gridX = fractionPosX * nx;
        double gridY = fractionPosY * ny;
        double gridZ = fractionPosZ * nz;
        
        // 计算网格索引和小数部分
        int gridIX = static_cast<int>(floor(gridX));
        int gridIY = static_cast<int>(floor(gridY));
        int gridIZ = static_cast<int>(floor(gridZ));
        
        double fractionX = gridX - gridIX;
        double fractionY = gridY - gridIY; 
        double fractionZ = gridZ - gridIZ;
        
        // 计算B样条的起始索引，考虑B样条的顺序
        int startIX = gridIX - order/2;
        if (startIX < 0) startIX += nx;
        
        int startIY = gridIY - order/2;
        if (startIY < 0) startIY += ny;
        
        int startIZ = gridIZ - order/2;
        if (startIZ < 0) startIZ += nz;
        
        std::cout << "网格坐标位置: [" << gridX << "," << gridY << "," << gridZ << "]\n";
        std::cout << "网格整数索引: [" << gridIX << "," << gridIY << "," << gridIZ << "]\n";
        std::cout << "网格小数部分: [" << fractionX << "," << fractionY << "," << fractionZ << "]\n";
        
        // 计算并显示B样条系数
        std::vector<double> bsCoeffsX(order), bsCoeffsY(order), bsCoeffsZ(order);
        
        // 使用已有的函数计算B样条系数
        computeBSplineCoefficients(fractionX, order, bsCoeffsX);
        computeBSplineCoefficients(fractionY, order, bsCoeffsY);
        computeBSplineCoefficients(fractionZ, order, bsCoeffsZ);
        
        std::cout << "X方向B样条系数: ";
        for (int i = 0; i < order; i++) {
            std::cout << bsCoeffsX[i] << " ";
        }
        std::cout << "\n";
        
        std::cout << "Y方向B样条系数: ";
        for (int i = 0; i < order; i++) {
            std::cout << bsCoeffsY[i] << " ";
        }
        std::cout << "\n";
        
        std::cout << "Z方向B样条系数: ";
        for (int i = 0; i < order; i++) {
            std::cout << bsCoeffsZ[i] << " ";
        }
        std::cout << "\n";
        
        // 添加与pme.cpp一致的电荷分布输出
        std::cout << "\n电荷分布到的网格点及其值:\n";
        
        // 使用不同的变量名，避免冲突
        int gridIndexX = gridIX;
        int gridIndexY = gridIY;
        int gridIndexZ = gridIZ;
        
        // 输出原子周围的电荷分布
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
                        
                        std::cout << "网格点[" << xindex << "," << yindex << "," << zindex
                                  << "], 索引: " << index 
                                  << ", 接收电荷: " << chargeContribution
                                  << ", 电荷系数: " << weight << std::endl;
                    }
                }
            }
        }
    }
    
    // 对某些原子的网格索引和B样条系数输出日志
    if (platform::verbose_ && platform::log_level_ <= LogLevel::DEBUG) {
        // 仅在DEBUG日志级别时输出少数原子的信息
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
    platform::log(LogLevel::INFO, "Performing forward FFT on PME grid");
    std::cout << "Performing forward FFT on PME grid" << std::endl;
    
    // 确定网格尺寸
    int nx = pme_params.meshSize[0];
    int ny = pme_params.meshSize[1];
    int nz = pme_params.meshSize[2];
    
    // 验证网格尺寸是2的幂
    if ((nx & (nx - 1)) != 0 || (ny & (ny - 1)) != 0 || (nz & (nz - 1)) != 0) {
        throw std::runtime_error("PME grid size must be a power of 2 for FFT");
    }
    
    // 保存FFT前的网格统计数据 - 使用与pme.cpp一致的方法：只检查实部
    int nonZeroBeforeFFT = 0;
    for (const auto& val : pme_params.pmeGrid) {
        if (std::abs(val.real()) > 1e-10) {
            nonZeroBeforeFFT++;
        }
    }
    
    platform::log(LogLevel::INFO, "Grid before FFT: non-zero points = ", nonZeroBeforeFFT);
    // 添加直接输出到控制台
    std::cout << "Grid before FFT: non-zero points = " << nonZeroBeforeFFT << std::endl;
    
    // 创建FFT前的备份
    fftGridBackup = pme_params.pmeGrid;
    
    // 执行3D FFT - 完全按照pme.cpp中的实现
    // 在pme.cpp中，这是通过fftw_execute调用CustomFFT::fft3D_forward实现的
    CustomFFT::fft3D_forward(pme_params.pmeGrid.data(), nx, ny, nz);
    
    // 简化FFT后的统计信息
    int nonZeroAfterFFT = 0;
    for (const auto& val : pme_params.pmeGrid) {
        if (std::norm(val) > 1e-10) {
            nonZeroAfterFFT++;
        }
    }
    
    // 只保留日志输出，移除控制台输出
    platform::log(LogLevel::INFO, "Grid after FFT: non-zero points = ", nonZeroAfterFFT);
    // 添加控制台输出，确保这个信息一定会显示
    std::cout << "Grid after FFT: non-zero points = " << nonZeroAfterFFT << std::endl;
    
    // 添加详细的FFT后关键网格点值输出 - 与pme.cpp一致
    std::cout << "\n===== FFT后的标准格点值比较 =====\n";
    std::cout << "格点总数: " << pme_params.pmeGrid.size() << std::endl;
    
    // 输出关键索引点
    const int keyIndices[] = {0, 1, 32, nx, ny, nz, nx*ny, nx*nz, ny*nz};
    for (int i : keyIndices) {
        if (i < static_cast<int>(pme_params.pmeGrid.size())) {
            std::cout << "格点[" << i << "]: " << pme_params.pmeGrid[i].real() 
                      << " + " << pme_params.pmeGrid[i].imag() << "i" << std::endl;
        }
    }
    
    // 特定的三维坐标
    const int keyCoords[][3] = {{0,0,1}, {0,1,0}, {1,0,0}, {1,1,1}, {2,2,2}};
    for (const auto& coord : keyCoords) {
        int idx = ((coord[0] % nx) * ny * nz) + ((coord[1] % ny) * nz) + (coord[2] % nz);
        if (idx < static_cast<int>(pme_params.pmeGrid.size())) {
            std::cout << "格点[" << coord[0] << "," << coord[1] << "," << coord[2] 
                      << "] (索引=" << idx << "): " << pme_params.pmeGrid[idx].real() 
                      << " + " << pme_params.pmeGrid[idx].imag() << "i" << std::endl;
        }
    }
    
    // 移除大部分详细的网格点值输出，仅在DEBUG级别保留少量信息
    if (platform::verbose_ && platform::log_level_ <= LogLevel::DEBUG) {
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
    // 记录开始处理
    platform::log(LogLevel::DEBUG, "Computing energy from PME grid");
    
    // 计算盒子体积并输出关键参数 - 与pme.cpp一致
    double volume = box[0] * box[1] * box[2];
    const double boxfactor = M_PI * volume;  // 用于调试输出和B样条调制因子
    
    std::cout << "Computing energy from grid with box = [" << box[0] << "," << box[1] << "," 
              << box[2] << "], alpha = " << pme_params.alpha << ", volume = " << volume << std::endl;
    
    // 初始化计数器和能量
    // 这些变量已经被替换为pointsProcessed和significantPoints
    // int nonZeroPoints = 0;
    // int significantEnergyPoints = 0;
    energy = 0.0;
    
    // 获取网格尺寸
    int nx = pme_params.meshSize[0];
    int ny = pme_params.meshSize[1];
    int nz = pme_params.meshSize[2];
    
    // 获取并计算常数
    double alpha = pme_params.alpha;
    
    // 计算倒易晶格矢量
    double recipBoxVectors[3][3] = {{0}};
    
    // 正确处理盒子向量
    double periodicBoxVectors[3][3] = {
        {box[0], 0.0, 0.0},
        {0.0, box[1], 0.0},
        {0.0, 0.0, box[2]}
    };
    
    // 对角盒子检查
    bool isDiagonalBox = true;  // 强制为对角盒子
    
    // 重要修复 - 使用与pme.cpp完全一致的倒格矢计算
    if (isDiagonalBox) {
        // 简化的对角盒子计算 - 与pme.cpp保持一致
        recipBoxVectors[0][0] = 2.0 * M_PI / box[0]; 
        recipBoxVectors[1][1] = 2.0 * M_PI / box[1]; 
        recipBoxVectors[2][2] = 2.0 * M_PI / box[2]; 
    } else {
        // 非对角盒子需要完全计算倒格矢向量
        double det = periodicBoxVectors[0][0] * (periodicBoxVectors[1][1] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][1]) -
                     periodicBoxVectors[0][1] * (periodicBoxVectors[1][0] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][0]) +
                     periodicBoxVectors[0][2] * (periodicBoxVectors[1][0] * periodicBoxVectors[2][1] - periodicBoxVectors[1][1] * periodicBoxVectors[2][0]);
        
        // 计算叉积和倒易晶格向量
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
    
    // 添加倒格矢输出 - 与pme.cpp保持一致
    std::cout << "Reciprocal lattice vectors:" << std::endl;
    for(int i=0; i<3; i++) {
        std::cout << "  b" << (i+1) << " = [" 
                  << recipBoxVectors[i][0] << ", "
                  << recipBoxVectors[i][1] << ", " 
                  << recipBoxVectors[i][2] << "]" << std::endl;
    }
    
    // 计算必要的因子 - 与pme.cpp保持一致
    const double one_4pi_eps = COULOMB / pme_params.epsilon_r;
    const double factor = M_PI * M_PI / (alpha * alpha);  // 修改为与pme.cpp一致的公式
    
    // 输出关键能量计算参数 - 与pme.cpp一致
    std::cout << "Energy parameters: one_4pi_eps = " << one_4pi_eps 
              << ", factor = " << factor << ", boxfactor = " << boxfactor << std::endl;
    
    std::cout << "Updating grid data before energy calculation" << std::endl;
    
    // 计算能量 - 精确匹配pme.cpp中的实现
    double energySum = 0.0;
    int pointsProcessed = 0;
    int significantPoints = 0;
    
    // 添加能量贡献记录，用于调试比较
    double contributions[3][3][3] = {{{0}}};
    
    // 确定最大k值，用于周期性处理 - 与pme.cpp保持一致
    int maxkx = (nx + 1) / 2;
    int maxky = (ny + 1) / 2;
    int maxkz = (nz + 1) / 2;
    
    // 因为我们将修改网格数据，所以不再需要副本
    // std::vector<std::complex<double>> gridCopy = pme_params.pmeGrid;
    
    // 初始化粘度相关变量 - 与pme.cpp一致
    double virxx = 0, virxy = 0, virxz = 0;
    double viryy = 0, viryz = 0, virzz = 0;
    
    for (int kx = 0; kx < nx; kx++) {
        // Calculate frequency - 使用与pme.cpp相同的变量名
        double mx = (kx < maxkx) ? kx : (kx - nx);
        double mhx = mx * recipBoxVectors[0][0];
        // 关键修复：在这里应用boxfactor，与pme.cpp保持一致
        double bx = boxfactor * pme_params.bsplineModuli[0][kx];
        
        for (int ky = 0; ky < ny; ky++) {
            double my = (ky < maxky) ? ky : (ky - ny);
            // 关键修复：mhy的计算方式需要与pme.cpp完全一致
            double mhy = mx * recipBoxVectors[1][0] + my * recipBoxVectors[1][1];
            double by = pme_params.bsplineModuli[1][ky];
            
            for (int kz = 0; kz < nz; kz++) {
                // Skip zero frequency term for neutral systems
                if (kx == 0 && ky == 0 && kz == 0) {
                    continue;
                }
                
                double mz = (kz < maxkz) ? kz : (kz - nz);
                // 关键修复：mhz的计算方式需要与pme.cpp完全一致
                double mhz = mx * recipBoxVectors[2][0] + my * recipBoxVectors[2][1] + mz * recipBoxVectors[2][2];
                
                // 网格索引
                int gridIndex = kx * ny * nz + ky * nz + kz;
                
                // 获取当前网格值
                double gridReal = pme_params.pmeGrid[gridIndex].real();
                double gridImag = pme_params.pmeGrid[gridIndex].imag();
                
                // 计算倒格点的长度平方 m²
                double m2 = mhx * mhx + mhy * mhy + mhz * mhz;
                
                // 获取第三维B样条调制因子
                double bz = pme_params.bsplineModuli[2][kz];
                
                // 计算分母 - 注意现在bx已经包含了boxfactor
                double denom = m2 * bx * by * bz;
                
                // 数值稳定性 - 与pme.cpp一致
                if (denom < 1e-10) {
                    denom = 1e-10;
                }
                
                // 关键修复 - 完全按照pme.cpp的方式计算eterm
                double eterm = one_4pi_eps * std::exp(-factor * m2) / denom;
                
                // 修改网格值 - 与pme.cpp保持一致
                // 这是关键的修复点: pme.cpp在这里更新了网格数据
                pme_params.pmeGrid[gridIndex].real(gridReal * eterm);
                pme_params.pmeGrid[gridIndex].imag(gridImag * eterm);
                
                // 计算能量贡献 - 修改以匹配pme.cpp
                double struct2 = gridReal * gridReal + gridImag * gridImag;
                double ets2 = eterm * struct2;
                energySum += ets2;
                
                // 统计处理点和有显著贡献的点
                pointsProcessed++;
                if (ets2 > 1e-8) {
                    significantPoints++;
                }
                
                // 记录关键点的能量贡献
                if (kx < 3 && ky < 3 && kz < 3) {
                    contributions[kx][ky][kz] = ets2;
                }
                
                // 调试输出 - 关键观察点，与pme.cpp一致
                if ((kx == 0 && ky == 0 && kz == 1) || 
                    (kx == 0 && ky == 1 && kz == 0) || 
                    (kx == 1 && ky == 0 && kz == 0)) {
                    std::cout << "比较点[" << kx << "," << ky << "," << kz << "]:" << std::endl;
                    std::cout << "  mx,my,mz = [" << mx << "," << my << "," << mz << "]" << std::endl;
                    std::cout << "  mhx,mhy,mhz = [" << mhx << "," << mhy << "," << mhz << "]" << std::endl;
                    std::cout << "  m2 = " << m2 << std::endl;
                    std::cout << "  bx,by,bz = [" << bx << "," << by << "," << bz << "]" << std::endl;
                    std::cout << "  denom = " << denom << std::endl;
                    std::cout << "  eterm = " << eterm << std::endl;
                    std::cout << "  原始网格值 = " << gridReal << " + " << gridImag << "i" << std::endl;
                    std::cout << "  更新后网格值 = " << pme_params.pmeGrid[gridIndex].real() << " + " << pme_params.pmeGrid[gridIndex].imag() << "i" << std::endl;
                    std::cout << "  struct2 = " << struct2 << std::endl;
                    std::cout << "  能量贡献 = " << ets2 << std::endl;
                    std::cout << "  累计能量 = " << energySum << std::endl;
                }
                
                // 计算粘度贡献 - 与pme.cpp一致
                double vfactor = 2.0 * (1.0 - factor * m2);
                virxx += vfactor * mhx * mhx * ets2;
                virxy += vfactor * mhx * mhy * ets2;
                virxz += vfactor * mhx * mhz * ets2;
                viryy += vfactor * mhy * mhy * ets2;
                viryz += vfactor * mhy * mhz * ets2;
                virzz += vfactor * mhz * mhz * ets2;
            }
        }
    }
    
    // 输出处理的点数和有显著贡献的点数 - 与pme.cpp一致
    std::cout << "Processed " << pointsProcessed << " grid points, "
              << significantPoints << " have significant energy contributions" << std::endl;
    
    // 输出关键点的能量贡献 - 与pme.cpp保持一致
    std::cout << "Energy contributions from key grid points:" << std::endl;
    for (int kx = 0; kx < 3; kx++) {
        for (int ky = 0; ky < 3; ky++) {
            for (int kz = 0; kz < 3; kz++) {
                if (kx == 0 && ky == 0 && kz == 0) continue; // 跳过零频率点
                std::cout << "  (" << kx << "," << ky << "," << kz << ") = " 
                          << contributions[kx][ky][kz] << std::endl;
            }
        }
    }
    
    // 计算并输出网格更新后的非零点数
    int nonZeroUpdated = 0;
    for (int i = 0; i < nx * ny * nz; i++) {
        if (std::norm(pme_params.pmeGrid[i]) > 1e-10) {
            nonZeroUpdated++;
        }
    }
    std::cout << "Grid after update: non-zero points = " << nonZeroUpdated << std::endl;
    
    // 计算最终能量 - 乘以0.5因子，与pme.cpp一致
    energy = 0.5 * energySum;
    
    // 输出能量总结 - 与pme.cpp一致
    std::cout << "Raw energy sum = " << energySum << std::endl;
    std::cout << "Final energy from grid (with scaling) = " << energy << std::endl;
    
    // 输出粘度项 - 与pme.cpp一致
    std::cout << "Virial terms: xx=" << virxx << ", xy=" << virxy << ", xz=" << virxz
              << ", yy=" << viryy << ", yz=" << viryz << ", zz=" << virzz << std::endl;
    
    // 添加详细的倒空间能量计算结果输出 - 与pme.cpp一致
    std::cout << "\n===== [energyPME] 倒空间能量计算结果 =====\n";
    std::cout << "总计算点数: " << pointsProcessed << std::endl;
    std::cout << "有意义能量点数: " << significantPoints << std::endl;
    std::cout << "原始能量和: " << energySum << std::endl;
    std::cout << "最终倒空间能量: " << energy << " (已乘以0.5)" << std::endl;
    
    // 添加能量计算后的标准格点输出 - 与pme.cpp保持一致
    std::cout << "\n===== [energyPME] 能量计算后的标准格点值比较 =====\n";
    const int keyIndices[] = {0, 1, nx, ny, nz, nx*ny, nx*nz, ny*nz};
    for (int i : keyIndices) {
        if (i < nx * ny * nz) {
            std::cout << "格点[" << i << "]: " << pme_params.pmeGrid[i].real() 
                      << " + " << pme_params.pmeGrid[i].imag() << "i" << std::endl;
        }
    }
    
    // 特定的三维坐标 - 与pme.cpp保持一致
    const int keyCoords[][3] = {{0,0,1}, {0,1,0}, {1,0,0}, {1,1,1}, {2,2,2}};
    for (const auto& coord : keyCoords) {
        int idx = ((coord[0] % nx) * ny * nz) + ((coord[1] % ny) * nz) + (coord[2] % nz);
        if (idx < nx * ny * nz) {
            std::cout << "格点[" << coord[0] << "," << coord[1] << "," << coord[2] 
                      << "] (索引=" << idx << "): " << pme_params.pmeGrid[idx].real() 
                      << " + " << pme_params.pmeGrid[idx].imag() << "i" << std::endl;
        }
    }
}

/**
 * @brief Compute reciprocal space energy using PME
 */
double computeReciprocalPME(model::MCState& state, bool movement_only) {
    platform::log(LogLevel::INFO, "Computing PME reciprocal space energy");
    
    // 精简系统信息输出
    platform::log(LogLevel::INFO, "Box: [", state.info.box[0], ", ", 
                 state.info.box[1], ", ", state.info.box[2], "], Alpha: ", pme_params.alpha);
    
    // 简化系统电荷检查
    double totalCharge = 0.0;
    for (int i = 0; i < state.activeAtomCount; i++) {
        totalCharge += state.atoms[i].charge;
    }
    
    if (std::abs(totalCharge) > 1e-6) {
        platform::log(LogLevel::WARNING, "System is not neutral! Total charge = ", totalCharge);
    }
    
    // PME网格检查
    if (pme_params.pmeGrid.empty()) {
        platform::log(LogLevel::ERROR, "PME grid not initialized!");
        return 0.0;
    }
    
    // 初始化B样条函数
    if (pme_params.bsplineModuli[0].empty() || 
        pme_params.bsplineModuli[1].empty() || 
        pme_params.bsplineModuli[2].empty()) {
        platform::log(LogLevel::INFO, "Initializing B-splines for PME calculation...");
        pme_params.initializeBsplines();
    }
    
    // 重置网格
    std::fill(pme_params.pmeGrid.begin(), pme_params.pmeGrid.end(), std::complex<double>(0.0, 0.0));
    
    // 执行PME计算步骤
    spreadChargesOntoGrid(state, movement_only);
    performFFTForward();
    
    // 将box转换为double数组
    double box[3];
    for (int i = 0; i < 3; i++) {
        box[i] = static_cast<double>(state.info.box[i]);
    }
    
    // 更新PME参数中的盒子尺寸，确保B样条和能量计算使用相同的体积
    pme_params.setBox(box);
    
    // 计算能量
    double energy = 0.0;
    computeEnergyFromGrid(energy, box);
    
    // 能量已经包含Coulomb常数，不需要再次乘以
    double reciprocal_energy = energy;
    
    platform::log(LogLevel::INFO, "Reciprocal space energy = ", reciprocal_energy);
    // 添加直接输出到控制台
    std::cout << "Reciprocal space energy = " << reciprocal_energy << std::endl;
    
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
    
    // 跟踪电荷平方和
    double sum_q2 = 0.0;
    int count = 0;
    
    if(movement_only) {
        for(const auto& movementInfo : state.movementResidues) {
            for(int i = movementInfo.startIndex; 
                i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                if(!state.residues[i].active) continue;
                
                for(int j = state.residues[i].atomStart;
                    j < state.residues[i].atomStart + state.residues[i].atomCount; j++) {
                    if (j < state.activeAtomCount) { // 确保索引有效
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
            // 不再检查active成员
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
    // 确保我们使用完全相同的公式，包括COULOMB常数
    double prefactor = -COULOMB * pme_params.alpha / sqrt(M_PI);
    self_energy = prefactor * sum_q2;
    
    platform::log(LogLevel::INFO, "Self energy prefactor = ", prefactor, 
                 ", resulting self energy = ", self_energy);
    
    // 自能量已经乘以COULOMB常数，不需要在别处再乘
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
                // 确保原子索引有效
                if(i >= state.activeAtomCount) continue;
                
                for(int j = residues[r2].atomStart; 
                    j < residues[r2].atomStart + residues[r2].atomCount; j++) {
                    // 确保原子索引有效
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
    
    // 对整个系统计算能量
    // 1. 首先计算实空间部分 - 需要乘以COULOMB因子
    computeRealSpacePME(state, false, true);
    
    // 2. 然后计算倒空间部分 - computeReciprocalPME已经包含COULOMB因子
    state.ewald_energy.reciprocal = computeReciprocalPME(state, false);
    
    // 3. 最后计算自能部分 - computeSelfEnergyPME已经包含COULOMB因子
    state.ewald_energy.self = computeSelfEnergyPME(state, false);
    
    // 仅对实空间能量乘以COULOMB系数
    state.ewald_energy.real_space *= COULOMB;
    
    // 计算总能量
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
    // 注意: computeReciprocalPME和computeSelfEnergyPME现在已经包含COULOMB因子
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
    platform::log(LogLevel::INFO, "Performing backward FFT on PME grid");
    
    // 确定网格尺寸
    int nx = pme_params.meshSize[0];
    int ny = pme_params.meshSize[1];
    int nz = pme_params.meshSize[2];
    
    // 验证网格尺寸是2的幂
    if ((nx & (nx - 1)) != 0 || (ny & (ny - 1)) != 0 || (nz & (nz - 1)) != 0) {
        throw std::runtime_error("PME grid size must be a power of 2 for FFT");
    }
    
    // 记录FFT前的网格统计
    int nonZeroBeforeFFT = 0;
    for (const auto& val : pme_params.pmeGrid) {
        if (std::norm(val) > 1e-10) {
            nonZeroBeforeFFT++;
        }
    }
    
    platform::log(LogLevel::INFO, "Grid before backward FFT: non-zero points = ", nonZeroBeforeFFT);
    
    // 执行3D反向FFT - 完全按照pme.cpp中的实现
    // 在pme.cpp中，这是通过fftw_execute调用CustomFFT::fft3D_backward实现的
    CustomFFT::fft3D_backward(pme_params.pmeGrid.data(), nx, ny, nz);
    
    // FFT后的简要统计
    int nonZeroAfterFFT = 0;
    for (const auto& val : pme_params.pmeGrid) {
        if (std::abs(val.real()) > 1e-10) {
            nonZeroAfterFFT++;
        }
    }
    
    platform::log(LogLevel::INFO, "Grid after backward FFT: non-zero points = ", nonZeroAfterFFT);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc





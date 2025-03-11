#include "energyPME.hpp"
#include "platform/platform.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <complex>

// Optional: Include library for FFT if needed
// #include <fftw3.h>

namespace {
// 定义PI2
const double PI2 = 6.28318530717958647692;

// 在C++中，我们使用std::complex<double>替代C的complex类型
using cmplx = std::complex<double>;


// 计算n位二进制数x的位反转值
int bit_reverse(int x, int n) {
    int result = 0;
    for (int i = 0; i < n; i++) {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    return result;
}

// 应用位反转排序使FFT结果与FFTW兼容
void apply_bit_reverse(cmplx* A, int k) {
    const int m = 1 << k;
    std::vector<cmplx> temp(m);
    
    // 将数据复制到临时数组，按位反转顺序重排
    for (int i = 0; i < m; i++) {
        int j = bit_reverse(i, k);
        temp[i] = A[j];
    }
    
    // 复制回原数组
    for (int i = 0; i < m; i++) {
        A[i] = temp[i];
    }
}

// 生成FFT权重
void tfft_genw(int i, int b, cmplx z, cmplx *w) {
    if(b == 0)
        w[i] = z;
    else {
        tfft_genw(i, b>>1, z, w);
        tfft_genw(i|b, b>>1, z*w[b], w);
    }
}

// 初始化FFT权重
void tfft_init(int k, cmplx *w) {
    int i, j;
    const int m = 1<<k;
    const double arg = -PI2/m;
    for(i=1, j=m/4; j; i<<=1, j>>=1) {
        w[i] = std::exp(std::complex<double>(0, arg * j));
    }
    tfft_genw(0, m/4, 1, w);
}

// 执行前向FFT
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
    
    // 添加位反转排序，使结果与FFTW兼容
    apply_bit_reverse(A, k);
}

// 执行反向FFT (使用共轭方法确保与FFTW一致)
void tfft_ifft(int k, cmplx *A, const cmplx *w) {
    const int m = 1 << k;
    
    // 步骤1: 对输入数据取共轭
    for (int i = 0; i < m; i++) {
        A[i] = std::conj(A[i]);
    }
    
    // 步骤2: 使用前向FFT变换
    tfft_fft(k, A, w);
    
    // 步骤3: 对结果再次取共轭并归一化
    for (int i = 0; i < m; i++) {
        A[i] = std::conj(A[i]) / static_cast<double>(m);
    }
}

// 卷积器
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

// 修正后的fft_1d_batch函数，正确处理3D索引
void fft_1d_batch(cmplx* data, int dimension, int nx, int ny, int nz, bool inverse) {
    // 根据维度选择变换大小和处理方式
    int transform_size;
    int num_transforms;
    
    switch (dimension) {
        case 0: // X方向
            transform_size = nx;
            num_transforms = ny * nz;
            break;
        case 1: // Y方向
            transform_size = ny;
            num_transforms = nx * nz;
            break;
        case 2: // Z方向
            transform_size = nz;
            num_transforms = nx * ny;
            break;
        default:
            throw std::invalid_argument("Invalid dimension for 3D FFT");
    }
    
    // 检查是否为2的幂
    if ((transform_size & (transform_size - 1)) != 0) {
        throw std::runtime_error("FFT size must be a power of 2");
    }
    
    // 计算log2(transform_size)
    int k = 0;
    int temp = transform_size;
    while (temp >>= 1) ++k;
    
    // 初始化FFT权重
    std::vector<cmplx> weights(transform_size);
    tfft_init(k, weights.data());
    
    // 临时缓冲区用于单个变换
    std::vector<cmplx> buffer(transform_size);
    
    // 处理每个一维变换
    for (int t = 0; t < num_transforms; t++) {
        int y, z, x;
        
        // 根据维度计算对应的二维索引 (t -> x,y,z 中的两个)
        switch (dimension) {
            case 0: // X方向: t = y + z*ny, 对每个y,z执行变换，遍历所有x
                y = t % ny;
                z = t / ny;
                
                // 将数据读入缓冲区
                for (int x = 0; x < nx; x++) {
                    buffer[x] = data[x + y*nx + z*nx*ny];
                }
                
                // 执行FFT/IFFT
                if (inverse) {
                    tfft_ifft(k, buffer.data(), weights.data());
                } else {
                    tfft_fft(k, buffer.data(), weights.data());
                }
                
                // 将结果写回
                for (int x = 0; x < nx; x++) {
                    data[x + y*nx + z*nx*ny] = buffer[x];
                }
                break;
                
            case 1: // Y方向: t = x + z*nx, 对每个x,z执行变换，遍历所有y
                x = t % nx;
                z = t / nx;
                
                // 将数据读入缓冲区
                for (int y = 0; y < ny; y++) {
                    buffer[y] = data[x + y*nx + z*nx*ny];
                }
                
                // 执行FFT/IFFT
                if (inverse) {
                    tfft_ifft(k, buffer.data(), weights.data());
                } else {
                    tfft_fft(k, buffer.data(), weights.data());
                }
                
                // 将结果写回
                for (int y = 0; y < ny; y++) {
                    data[x + y*nx + z*nx*ny] = buffer[y];
                }
                break;
                
            case 2: // Z方向: t = x + y*nx, 对每个x,y执行变换，遍历所有z
                x = t % nx;
                y = t / nx;
                
                // 将数据读入缓冲区
                for (int z = 0; z < nz; z++) {
                    buffer[z] = data[x + y*nx + z*nx*ny];
                }
                
                // 执行FFT/IFFT
                if (inverse) {
                    tfft_ifft(k, buffer.data(), weights.data());
                } else {
                    tfft_fft(k, buffer.data(), weights.data());
                }
                
                // 将结果写回
                for (int z = 0; z < nz; z++) {
                    data[x + y*nx + z*nx*ny] = buffer[z];
                }
                break;
        }
    }
}

// 使用修正的fft_1d_batch函数重写3D前向FFT
[[maybe_unused]]
void fft3D_forward(cmplx* data, int nx, int ny, int nz) {
    // 验证所有尺寸都是2的幂
    if ((nx & (nx - 1)) != 0 || (ny & (ny - 1)) != 0 || (nz & (nz - 1)) != 0) {
        throw std::runtime_error("All 3D FFT dimensions must be powers of 2");
    }
    
    // 按X、Y、Z顺序执行三次1D FFT
    fft_1d_batch(data, 0, nx, ny, nz, false); // X方向
    fft_1d_batch(data, 1, nx, ny, nz, false); // Y方向
    fft_1d_batch(data, 2, nx, ny, nz, false); // Z方向
}

// 使用修正的fft_1d_batch函数重写3D逆FFT
[[maybe_unused]]
void fft3D_backward(cmplx* data, int nx, int ny, int nz) {
    // 验证所有尺寸都是2的幂
    if ((nx & (nx - 1)) != 0 || (ny & (ny - 1)) != 0 || (nz & (nz - 1)) != 0) {
        throw std::runtime_error("All 3D FFT dimensions must be powers of 2");
    }
    
    // 按Z、Y、X顺序执行三次1D IFFT (与前向顺序相反)
    fft_1d_batch(data, 2, nx, ny, nz, true); // Z方向
    fft_1d_batch(data, 1, nx, ny, nz, true); // Y方向
    fft_1d_batch(data, 0, nx, ny, nz, true); // X方向
}

} // anonymous namespace

namespace pygcmc {
namespace platform {
namespace cpu {

// Global parameters instance
PMEParams pme_params;

// 添加FFT相关成员变量
std::vector<std::complex<double>> fft_weights;

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
 * @brief Initialize B-splines for PME
 */
void PMEParams::initializeBsplines() {
    // 初始化B-spline模数
    platform::log(LogLevel::INFO, "Initializing B-splines with order = ", splineOrder, 
                 " and mesh size = [", meshSize[0], ",", meshSize[1], ",", meshSize[2], "]");
    
    // 确保样条阶数至少为2
    if (splineOrder < 2) {
        platform::log(LogLevel::WARNING, "B-spline order must be at least 2, setting to 2");
        splineOrder = 2;
    }
    
    // 确保样条阶数在合理范围内
    if (splineOrder > 6) {
        platform::log(LogLevel::WARNING, "B-spline orders > 6 may require special handling, consider using order 4-6 for optimal performance");
    }
    
    // 初始化 bsplineModuli 数组 - 精确匹配OpenMM的实现
    for (int dim = 0; dim < 3; dim++) {
        int size = meshSize[dim];
        bsplineModuli[dim].resize(size);
        
        // 初始化为0
        for (int i = 0; i < size; i++)
            bsplineModuli[dim][i] = 0.0;
            
        // 设置DC分量(k=0)
        bsplineModuli[dim][0] = 1.0;
        
        // 计算B样条模数 - 精确匹配OpenMM的计算方法
        double factor = M_PI / size;
        
        // 对每个k向量计算B样条模数
        for (int i = 1; i < size; i++) {
            int m = (i < size/2) ? i : (size - i);
            
            if (m == 0) {
                bsplineModuli[dim][i] = 1.0;
                continue;
            }
            
            // 计算B样条模数
            double numerator = 0.0;
            // 处理小角度情况
            double w = m * factor;
            if (w < 1e-7) {
                // 小角度近似
                numerator = 1.0;
            } else {
                // 标准计算 - 精确匹配OpenMM
                numerator = std::sin(splineOrder * w) / (splineOrder * std::sin(w));
            }
            
            // 计算B样条函数的傅里叶变换
            double bspline;
            // 对于任意阶数的B样条，使用相应的幂
            bspline = std::pow(numerator, splineOrder);
            
            // 修正系数 - 关键是这里的算法
            if (splineOrder > 4) {
                // 对于高阶样条的特殊处理
                double eps = 1.0e-7;
                if (bspline < eps && m <= splineOrder) {
                    // 高阶低频处理 - 这是OpenMM使用的方法
                    double sum = 0.0;
                    for (int j = 1; j <= splineOrder; j++) {
                        double term = std::sin(w*j) / (w*j);
                        sum += term*term;
                    }
                    // 避免除以零
                    if (sum < eps)
                        bspline = 0.0;
                    else
                        bspline = 1.0 / (sum * size * size);
                }
            } else {
                // 对于4阶及以下的处理
                // 低阶样条的额外检查
                if (bspline < 1e-10 && m <= splineOrder) {
                    double sum = 0.0;
                    for (int j = 1; j <= splineOrder; j++) {
                        double term = std::sin(w*j) / (w*j);
                        sum += term*term;
                    }
                    if (sum > 1e-10)
                        bspline = 1.0 / (sum * size * size);
                }
            }
            
            // 存储结果
            bsplineModuli[dim][i] = bspline;
            
            // 对于重要频率，确保有合理的非零值
            if (m <= splineOrder && bspline < 1e-10) {
                platform::log(LogLevel::WARNING, "Very small B-spline modulus detected for m=", m, 
                             ", setting to minimum value");
                bsplineModuli[dim][i] = 1e-10;
            }
            
            // 记录调试信息
            if (i <= 10 || i >= size-10 || m <= splineOrder) {
                platform::log(LogLevel::DEBUG, "B-spline[", dim, "][", i, "] = ", bsplineModuli[dim][i],
                            " (m=", m, ")");
            }
        }
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
    
    // Set the parameters
    setPMEParameters(alpha, meshSize);
    
    // 关键修复：初始化查找表并设置initialized标志
    pme_params.initializeTables(cutoff_distance);
    
    // 确保B样条也被初始化
    pme_params.initializeBsplines();
}

/**
 * @brief Set PME parameters explicitly
 * 
 * @param alpha Ewald separation parameter
 * @param meshSize Grid dimensions for PME
 * @param splineOrder B-spline order
 * @param tolerance Precision parameter
 */
void setPMEParameters(double alpha, const int meshSize[3], int splineOrder, double tolerance) {
    pme_params.alpha = alpha;
    
    // 确保网格尺寸是2的幂
    for (int i = 0; i < 3; i++) {
        if ((meshSize[i] & (meshSize[i] - 1)) != 0) {
            // 如果不是2的幂，找到最近的2的幂
            int log2_size = 0;
            while ((1 << log2_size) < meshSize[i]) log2_size++;
            pme_params.meshSize[i] = 1 << log2_size;
            platform::log(LogLevel::WARNING, "PME mesh size must be a power of 2. Adjusting dimension ", 
                         i, " from ", meshSize[i], " to ", pme_params.meshSize[i]);
        } else {
            pme_params.meshSize[i] = meshSize[i];
        }
    }
    
    pme_params.splineOrder = splineOrder;
    pme_params.tolerance = tolerance;
    
    platform::log(LogLevel::INFO, "PME parameters set: alpha = ", alpha,
                 ", mesh size = [", pme_params.meshSize[0], ",", pme_params.meshSize[1], ",", pme_params.meshSize[2], "]",
                 ", spline order = ", splineOrder);
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
 * @brief Compute B-spline coefficients
 * 
 * @param fractional Fractional position (0-1)
 * @param order Spline order
 * @param coefficients Output coefficients
 */
void computeBSplineCoefficients(double fractional, int order, std::vector<double>& coefficients) {
    coefficients.resize(order + 1);
    
    // For order 1, it's just linear interpolation
    if (order == 1) {
        coefficients[0] = 1.0 - fractional;
        coefficients[1] = fractional;
        return;
    }
    
    // B-spline recursive formula for higher orders
    coefficients[0] = 1.0;
    
    for (int k = 1; k <= order; k++) {
        double div = 1.0 / k;
        double term = fractional * coefficients[0];
        coefficients[0] *= (1.0 - fractional) * div;
        
        for (int i = 1; i < k; i++) {
            double saved = term;
            term = fractional * coefficients[i];
            coefficients[i] = ((1.0 - fractional) * coefficients[i-1] + saved) * div;
        }
        
        coefficients[k] = term * div;
    }
}

/**
 * @brief Spread charges onto the PME grid
 * 
 * @param state MC state
 * @param movement_only Whether to process only moving atoms
 */
void spreadChargesOntoGrid(model::MCState& state, bool movement_only) {
    const auto& box = state.info.box;
    const auto& atoms = state.atoms;
    
    platform::log(LogLevel::INFO, "Spreading charges onto PME grid for ", 
                 movement_only ? "moving atoms" : "all atoms");
    platform::log(LogLevel::INFO, "Grid dimensions: [",
                 pme_params.meshSize[0], ",", pme_params.meshSize[1], ",", pme_params.meshSize[2], "]");
    
    // Reset the grid
    std::fill(pme_params.pmeGrid.begin(), pme_params.pmeGrid.end(), std::complex<double>(0.0, 0.0));
    
    // 计数器 - 了解有多少原子的电荷被扩散
    int processedAtoms = 0;
    double totalCharge = 0.0;
    
    // Calculate spline coefficients for each atom
    const int order = pme_params.splineOrder;
    std::vector<double> splineCoefficients(order+1);
    
    // For each atom, spread its charge on the grid
    for (int n = 0; n < state.activeAtomCount; n++) {
        // Skip if not in movement group for movement-only calculation
        if (movement_only) {
            bool in_movement = false;
            for (const auto& movementInfo : state.movementResidues) {
                if (n >= movementInfo.startIndex && 
                    n < movementInfo.startIndex + movementInfo.activeCount) {
                    in_movement = true;
                    break;
                }
            }
            if (!in_movement) continue;
        }
        
        double charge = atoms[n].charge;
        totalCharge += charge;
        
        if (std::abs(charge) < 1e-10) continue; // Skip atoms with zero charge
        
        // Convert atom coordinates to fractional coordinates (0-1 range)
        double posInBox[3];
        posInBox[0] = atoms[n].x;
        posInBox[1] = atoms[n].y;
        posInBox[2] = atoms[n].z;
        
        // Apply periodic boundary conditions
        // Wrap coordinates to primary box (0 to box size)
        while (posInBox[0] < 0) posInBox[0] += box[0];
        while (posInBox[0] >= box[0]) posInBox[0] -= box[0];
        while (posInBox[1] < 0) posInBox[1] += box[1];
        while (posInBox[1] >= box[1]) posInBox[1] -= box[1];
        while (posInBox[2] < 0) posInBox[2] += box[2];
        while (posInBox[2] >= box[2]) posInBox[2] -= box[2];
        
        // Convert to fractional coordinates (0-1)
        double scale[3] = { 1.0 / box[0], 1.0 / box[1], 1.0 / box[2] };
        double t[3] = { posInBox[0] * scale[0], posInBox[1] * scale[1], posInBox[2] * scale[2] };
        
        // Convert to grid coordinates and compute B-spline coefficients
        int gridIndex[3];
        double dr[3];
        std::vector<std::vector<double>> thetai(3, std::vector<double>(order));
        
        for (int dim = 0; dim < 3; dim++) {
            // Grid coordinates
            t[dim] = t[dim] * pme_params.meshSize[dim];
            
            // Integer and fractional parts
            gridIndex[dim] = (int) std::floor(t[dim]);
            dr[dim] = t[dim] - gridIndex[dim];
            
            // Ensure grid index is within bounds
            gridIndex[dim] = gridIndex[dim] % pme_params.meshSize[dim];
            if (gridIndex[dim] < 0) gridIndex[dim] += pme_params.meshSize[dim];
            
            // Compute B-spline coefficients - optimized version based on OpenMM
            if (order == 4) {
                // Replace incorrect optimization with correct 4th-order B-spline calculation
                double w = dr[dim];
                double w2 = w * w;
                double w3 = w2 * w;
                double oneSixth = 1.0 / 6.0;
                thetai[dim][0] = oneSixth * (1.0 - w) * (1.0 - w) * (1.0 - w);
                thetai[dim][1] = oneSixth * (4.0 - 6.0 * w2 + 3.0 * w3);
                thetai[dim][2] = oneSixth * (1.0 + 3.0 * w + 3.0 * w2 - 3.0 * w3);
                thetai[dim][3] = oneSixth * w3;
            } else {
                // General B-spline calculation for other orders
                computeBSplineCoefficients(dr[dim], order, splineCoefficients);
                
                // Copy coefficients to thetai array
                for (int i = 0; i < order; i++) {
                    thetai[dim][i] = splineCoefficients[i];
                }
            }
        }
        
        // Spread charge to nearby grid points using B-spline weights
        // Optimized for better cache performance and numerical stability
        for (int ix = 0; ix < order; ix++) {
            int xindex = (gridIndex[0] + ix) % pme_params.meshSize[0];
            double xterm = charge * thetai[0][ix];
            
            for (int iy = 0; iy < order; iy++) {
                int yindex = (gridIndex[1] + iy) % pme_params.meshSize[1];
                double xyterm = xterm * thetai[1][iy];
                
                for (int iz = 0; iz < order; iz++) {
                    int zindex = (gridIndex[2] + iz) % pme_params.meshSize[2];
                    double weight = xyterm * thetai[2][iz];
                    
                    // Calculate linear index
                    int index = xindex * pme_params.meshSize[1] * pme_params.meshSize[2] + 
                                yindex * pme_params.meshSize[2] + zindex;
                    
                    // Atomic increment of grid value to ensure thread safety
                    pme_params.pmeGrid[index] += std::complex<double>(weight, 0.0);
                }
            }
        }
        
        processedAtoms++;
    }
    
    platform::log(LogLevel::INFO, "Processed ", processedAtoms, " atoms for PME. Total charge = ", totalCharge);
}

/**
 * @brief Perform forward FFT on the grid
 * 
 * Uses custom FFT implementation
 */
void performFFTForward() {
    platform::log(LogLevel::DEBUG, "Performing 3D forward FFT on PME grid");
    
    // 确定FFT的大小 (确保是2的幂)
    int nx = pme_params.meshSize[0];
    int ny = pme_params.meshSize[1];
    int nz = pme_params.meshSize[2];
    
    // 验证网格尺寸是2的幂
    if ((nx & (nx - 1)) != 0 || (ny & (ny - 1)) != 0 || (nz & (nz - 1)) != 0) {
        throw std::runtime_error("PME grid size must be a power of 2 for this FFT implementation");
    }
    
    // 直接使用fft3D_forward函数进行3D FFT
    fft3D_forward(pme_params.pmeGrid.data(), nx, ny, nz);
}

/**
 * @brief Perform backward FFT on the grid
 * 
 * Uses custom FFT implementation
 */
void performFFTBackward() {
    platform::log(LogLevel::DEBUG, "Performing 3D backward FFT on PME grid");
    
    // 确定FFT的大小 (确保是2的幂)
    int nx = pme_params.meshSize[0];
    int ny = pme_params.meshSize[1];
    int nz = pme_params.meshSize[2];
    
    // 验证网格尺寸是2的幂
    if ((nx & (nx - 1)) != 0 || (ny & (ny - 1)) != 0 || (nz & (nz - 1)) != 0) {
        throw std::runtime_error("PME grid size must be a power of 2 for this FFT implementation");
    }
    
    // 直接使用fft3D_backward函数进行3D逆FFT
    fft3D_backward(pme_params.pmeGrid.data(), nx, ny, nz);
    
    // 应用额外的缩放因子 - fft3D_backward中有些缩放，但我们额外需要全局缩放
    double scale = static_cast<double>(nx * ny * nz);
    for (size_t i = 0; i < pme_params.pmeGrid.size(); i++) {
        pme_params.pmeGrid[i] /= scale;
    }
}

/**
 * @brief Compute energy from the PME grid after FFT
 * 
 * @param energy Output energy
 */
void computeEnergyFromGrid(double& energy, const double box[3]) {
    double volume = box[0] * box[1] * box[2];
    double scaleFactor = COULOMB * 4.0 * M_PI / volume;
    
    energy = 0.0;
    
    platform::log(LogLevel::INFO, "Computing energy from grid with box = [", 
                 box[0], ",", box[1], ",", box[2], "], volume = ", volume);
    platform::log(LogLevel::INFO, "Mesh size = [", 
                 pme_params.meshSize[0], ",", pme_params.meshSize[1], ",", pme_params.meshSize[2], "]");
    
    // 计算倒空间能量贡献总和
    double totalContribution = 0.0;
    
    // 检查 B 样条模数是否正确初始化
    bool bsplinesInitialized = true;
    for (int dim = 0; dim < 3; dim++) {
        if (pme_params.bsplineModuli[dim].empty()) {
            bsplinesInitialized = false;
            platform::log(LogLevel::ERROR, "B-spline moduli for dimension ", dim, " not initialized!");
        }
    }
    
    if (!bsplinesInitialized) {
        platform::log(LogLevel::WARNING, "Initializing B-splines before energy calculation");
        pme_params.initializeBsplines();
    }
    
    // Compute reciprocal space energy from the grid
    for (int ix = 0; ix < pme_params.meshSize[0]; ix++) {
        int kx = ix;
        if (kx > pme_params.meshSize[0]/2) kx -= pme_params.meshSize[0];
        
        for (int iy = 0; iy < pme_params.meshSize[1]; iy++) {
            int ky = iy;
            if (ky > pme_params.meshSize[1]/2) ky -= pme_params.meshSize[1];
            
            for (int iz = 0; iz < pme_params.meshSize[2]; iz++) {
                int kz = iz;
                if (kz > pme_params.meshSize[2]/2) kz -= pme_params.meshSize[2];
                
                // Skip k = 0 case
                if (kx == 0 && ky == 0 && kz == 0) continue;
                
                int gridIndex = ix * pme_params.meshSize[1] * pme_params.meshSize[2] +
                                iy * pme_params.meshSize[2] + iz;
                
                double kx_sq = kx * kx;
                double ky_sq = ky * ky;
                double kz_sq = kz * kz;
                
                // Compute squared magnitude of k-vector
                double msq = (4.0 * M_PI * M_PI) * (
                           (kx_sq / (box[0] * box[0])) + 
                           (ky_sq / (box[1] * box[1])) + 
                           (kz_sq / (box[2] * box[2])));
                
                // B-spline influence function
                double bx = pme_params.bsplineModuli[0][ix];
                double by = pme_params.bsplineModuli[1][iy];
                double bz = pme_params.bsplineModuli[2][iz];
                
                // 确保 B 样条权重不为零 (关键修复)
                if (msq > 0 && (bx == 0.0 || by == 0.0 || bz == 0.0)) {
                    // 为非零频率重新计算 B 样条权重
                    double factor_x = 2.0 * M_PI / pme_params.meshSize[0];
                    double factor_y = 2.0 * M_PI / pme_params.meshSize[1];
                    double factor_z = 2.0 * M_PI / pme_params.meshSize[2];
                    
                    int m_x = ix > pme_params.meshSize[0]/2 ? pme_params.meshSize[0] - ix : ix;
                    int m_y = iy > pme_params.meshSize[1]/2 ? pme_params.meshSize[1] - iy : iy;
                    int m_z = iz > pme_params.meshSize[2]/2 ? pme_params.meshSize[2] - iz : iz;
                    
                    if (m_x > 0 && bx == 0.0) {
                        double w = factor_x * m_x;
                        bx = 1.0;
                        for (int j = 2; j <= pme_params.splineOrder; j++) {
                            double sin_term = std::sin(j * w / 2.0) / w;
                            bx = 4.0 * sin_term * sin_term * bx;
                        }
                    }
                    
                    if (m_y > 0 && by == 0.0) {
                        double w = factor_y * m_y;
                        by = 1.0;
                        for (int j = 2; j <= pme_params.splineOrder; j++) {
                            double sin_term = std::sin(j * w / 2.0) / w;
                            by = 4.0 * sin_term * sin_term * by;
                        }
                    }
                    
                    if (m_z > 0 && bz == 0.0) {
                        double w = factor_z * m_z;
                        bz = 1.0;
                        for (int j = 2; j <= pme_params.splineOrder; j++) {
                            double sin_term = std::sin(j * w / 2.0) / w;
                            bz = 4.0 * sin_term * sin_term * bz;
                        }
                    }
                }
                
                // 添加(2π)²因子以修正波矢计算
                double m2 = (4.0 * M_PI * M_PI) * msq;
                
                // B样条修正：B样条模数应该在分母而不是乘在分母的m2上
                // 只需m2用于exponent term，B样条模数单独作为分母
                double bsplineProduct = bx * by * bz;
                // 处理极小的分母值
                if (bsplineProduct < 1e-10) {
                    // Skip this term if it would cause instability
                    if (std::abs(bsplineProduct) < 1e-12)
                        continue;
                    
                    // For very small but non-zero bspline product, use a minimum value
                    bsplineProduct = 1e-10;
                }
                
                // 由于m2现已包含(2π)²因子，应相应调整指数项
                // 从 exp(-π²·m2/(α²)) 更改为 exp(-m2/(4·α²))
                double m2_term = m2 != 0.0 ? std::exp(-m2 / (4.0 * pme_params.alpha * pme_params.alpha)) / m2 : 0.0;
                
                // Get squared magnitude of complex grid value
                double gridMagnitudeSq = std::norm(pme_params.pmeGrid[gridIndex]);
                
                // 增加0.5系数并额外除以一次bsplineProduct，以确保正确的B样条模数幂次
                // 原始公式: double energyTerm = scaleFactor * m2_term * gridMagnitudeSq / bsplineProduct;
                // 修正为使用与splineOrder匹配的幂次
                double bsplinePower = std::pow(bsplineProduct, pme_params.splineOrder/2.0);
                if (bsplinePower < 1e-12) bsplinePower = 1e-12;
                double energyTerm = 0.5 * scaleFactor * m2_term * gridMagnitudeSq / bsplinePower;
                
                // 记录一些能量贡献值用于调试
                if ((ix <= 2 && iy <= 2 && iz <= 2) || gridMagnitudeSq > 1e-6) {
                    platform::log(LogLevel::DEBUG, "Grid[", ix, ",", iy, ",", iz, "] = ", 
                                 gridMagnitudeSq, ", bx*by*bz = ", bx*by*bz, 
                                 ", m2_term = ", m2_term, ", term = ", energyTerm);
                }
                
                energy += energyTerm;
                totalContribution += std::abs(energyTerm);
            }
        }
    }
    
    platform::log(LogLevel::INFO, "Total reciprocal energy = ", energy, 
                 ", total contribution = ", totalContribution);
}

/**
 * @brief Compute reciprocal space energy using PME
 * 
 * @param state MC state
 * @param movement_only Whether to compute only for moving atoms
 * @return double Reciprocal space energy
 */
double computeReciprocalPME(model::MCState& state, bool movement_only) {
    // 使用state.info.box，这是正确的盒子尺寸
    const auto& box = state.info.box;
    // 将float盒子尺寸转换为double类型
    double box_double[3] = {static_cast<double>(box[0]), 
                           static_cast<double>(box[1]), 
                           static_cast<double>(box[2])};
    const auto& atoms = state.atoms;
    
    // Check system neutrality
    double totalCharge = 0.0;
    for(const auto& atom : atoms) {
        totalCharge += static_cast<double>(atom.charge);
    }
    if (std::abs(totalCharge) > 1e-10) {
        throw std::runtime_error("System must be charge neutral for PME calculation");
    }
    
    platform::log(LogLevel::INFO, "Computing PME reciprocal energy with alpha = ", 
                 pme_params.alpha, ", mesh size = [", 
                 pme_params.meshSize[0], ",", pme_params.meshSize[1], ",", pme_params.meshSize[2], "]",
                 ", spline order = ", pme_params.splineOrder);
    
    // Spread charges onto grid
    spreadChargesOntoGrid(state, movement_only);
    
    // 调试 - 检查网格上的电荷分布
    double gridChargeSum = 0.0;
    for (const auto& val : pme_params.pmeGrid) {
        gridChargeSum += val.real();
    }
    platform::log(LogLevel::DEBUG, "Total charge on grid before FFT: ", gridChargeSum);
    
    // Perform forward FFT
    performFFTForward();
    
    // Compute reciprocal space energy
    double recipEnergy = 0.0;
    double volume = box_double[0] * box_double[1] * box_double[2];
    
    if (volume < 1e-10)
        throw std::runtime_error("Box volume too small for PME calculation");
    
    // Scale factor - 直接从OpenMM复制的公式
    double scaleFactor = (1.0/(2.0*M_PI*volume)) * 4.0 * M_PI * COULOMB;
    
    // 跟踪最大能量贡献，用于调试
    double maxEnergyTerm = 0.0;
    int maxTermIndices[3] = {0, 0, 0};
    
    // Calculate reciprocal space energy - 精确匹配OpenMM的计算方法
    for (int mx = 0; mx < pme_params.meshSize[0]; mx++) {
        int kx = (mx < pme_params.meshSize[0]/2) ? mx : mx - pme_params.meshSize[0];
        double mhx = kx / box_double[0];
        
        for (int my = 0; my < pme_params.meshSize[1]; my++) {
            int ky = (my < pme_params.meshSize[1]/2) ? my : my - pme_params.meshSize[1];
            double mhy = ky / box_double[1];
            
            for (int mz = 0; mz < pme_params.meshSize[2]; mz++) {
                int kz = (mz < pme_params.meshSize[2]/2) ? mz : mz - pme_params.meshSize[2];
                double mhz = kz / box_double[2];
                
                // Skip k=0 (DC component)
                if (kx == 0 && ky == 0 && kz == 0)
                    continue;
                
                // Calculate squared reciprocal vector length
                double msq = mhx*mhx + mhy*mhy + mhz*mhz;
                
                // Get grid index
                int gridIndex = mx*pme_params.meshSize[1]*pme_params.meshSize[2] + 
                               my*pme_params.meshSize[2] + mz;
                
                // Get B-spline values - 注意这里我们使用mx/my/mz而不是kx/ky/kz
                double bx = pme_params.bsplineModuli[0][mx];
                double by = pme_params.bsplineModuli[1][my];
                double bz = pme_params.bsplineModuli[2][mz];
                
                // Calculate influence function
                // 添加(2π)²因子以修正波矢计算
                double m2 = (4.0 * M_PI * M_PI) * msq;
                // B样条修正：B样条模数应该在分母而不是乘在分母的m2上
                // 只需m2用于exponent term，B样条模数单独作为分母
                double bsplineProduct = bx * by * bz;
                // 处理极小的分母值
                if (bsplineProduct < 1e-10) {
                    // Skip this term if it would cause instability
                    if (std::abs(bsplineProduct) < 1e-12)
                        continue;
                    
                    // For very small but non-zero bspline product, use a minimum value
                    bsplineProduct = 1e-10;
                }
                
                // 由于m2现已包含(2π)²因子，应相应调整指数项
                // 从 exp(-π²·m2/(α²)) 更改为 exp(-m2/(4·α²))
                double m2_term = m2 != 0.0 ? std::exp(-m2 / (4.0 * pme_params.alpha * pme_params.alpha)) / m2 : 0.0;
                
                // Get squared magnitude of complex grid value
                double gridMagnitudeSq = std::norm(pme_params.pmeGrid[gridIndex]);
                
                // 增加0.5系数并额外除以一次bsplineProduct，以确保正确的B样条模数幂次
                // 原始公式: double energyTerm = scaleFactor * m2_term * gridMagnitudeSq / bsplineProduct;
                // 修正为使用与splineOrder匹配的幂次
                double bsplinePower = std::pow(bsplineProduct, pme_params.splineOrder/2.0);
                if (bsplinePower < 1e-12) bsplinePower = 1e-12;
                double energyTerm = 0.5 * scaleFactor * m2_term * gridMagnitudeSq / bsplinePower;
                
                // 跟踪最大能量贡献
                if (std::abs(energyTerm) > std::abs(maxEnergyTerm)) {
                    maxEnergyTerm = energyTerm;
                    maxTermIndices[0] = mx;
                    maxTermIndices[1] = my;
                    maxTermIndices[2] = mz;
                }
                
                // Accumulate energy
                recipEnergy += energyTerm;
            }
        }
    }
    
    // 记录最大能量贡献，用于诊断问题
    platform::log(LogLevel::DEBUG, "Max energy term = ", maxEnergyTerm,
                 " at indices [", maxTermIndices[0], ",", 
                 maxTermIndices[1], ",", maxTermIndices[2], "]");
    
    platform::log(LogLevel::INFO, "PME reciprocal energy = ", recipEnergy);
    return recipEnergy;
}

/**
 * @brief Compute self-energy term for PME
 * 
 * @param state MC state
 * @param movement_only Whether to compute only for moving atoms
 * @return double Self-energy
 */
double computeSelfEnergyPME(model::MCState& state, bool movement_only) {
    // Self-energy calculation same as Ewald
    double self_energy = 0.0;
    
    if(movement_only) {
        for(const auto& movementInfo : state.movementResidues) {
            for(int i = movementInfo.startIndex; 
                i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                if(!state.residues[i].active) continue;
                
                for(int j = state.residues[i].atomStart;
                    j < state.residues[i].atomStart + state.residues[i].atomCount; j++) {
                    double charge = state.atoms[j].charge;
                    self_energy += charge * charge;
                }
            }
        }
    } else {
        for(int i = 0; i < state.activeAtomCount; i++) {
            double charge = state.atoms[i].charge;
            self_energy += charge * charge;
        }
    }
    
    // Self-energy formula same as Ewald
    self_energy = -COULOMB * pme_params.alpha / std::sqrt(M_PI) * self_energy;
    
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
                for(int j = residues[r2].atomStart; 
                    j < residues[r2].atomStart + residues[r2].atomCount; j++) {
                    
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
    // Output Coulomb constant value for debugging
    platform::log(LogLevel::INFO, "COULOMB constant in energyPME.cpp = ", COULOMB);

    if (!pme_params.initialized) {
        throw std::runtime_error("PME parameters not initialized. Call initializePMEParameters() first.");
    }
    
    // Calculate all components of PME energy
    computeRealSpacePME(state, false, true);
    state.ewald_energy.reciprocal = computeReciprocalPME(state, false);
    state.ewald_energy.self = computeSelfEnergyPME(state, false);
    
    // Apply Coulomb factor to all components
    state.ewald_energy.real_space *= COULOMB;
    
    // Total energy is the sum of all components
    state.ewald_energy.total = state.ewald_energy.real_space + 
                             state.ewald_energy.reciprocal + 
                             state.ewald_energy.self;
    
    platform::log(LogLevel::INFO, "PME energy components: real_space=", state.ewald_energy.real_space,
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
    
    // Apply Coulomb factor to all components
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

} // namespace cpu
} // namespace platform
} // namespace pygcmc




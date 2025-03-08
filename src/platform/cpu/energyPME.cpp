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
}

// 执行反向FFT
void tfft_ifft(int k, cmplx *A, const cmplx *w) {
    const int m = 1 << k;
    int u = m/4;
    int v = 1;
    int i, j;
    for(i=2;i<=k;i+=2) {
        int jh;
        for(jh=0;jh<u;jh++) {
            cmplx wj = std::conj(w[jh<<1]);
            cmplx wj2 = std::conj(w[jh]);
            cmplx wj3 = wj2 * wj;
            int je;
            for(j = jh << i, je = j+v;j<je; j++) {
                cmplx tmp0 = A[j];
                cmplx tmp1 = A[j+v];
                cmplx tmp2 = A[j+2*v];
                cmplx tmp3 = A[j+3*v];

                cmplx ttmp0 = tmp0 + tmp1;
                cmplx ttmp1 = tmp0 - tmp1;
                cmplx ttmp2 = tmp2 + tmp3;
                cmplx ttmp3 = std::complex<double>(0, 1) * (tmp2 - tmp3);

                A[j] = ttmp0 + ttmp2;
                A[j+v] = wj * (ttmp1 + ttmp3);
                A[j+2*v] = wj2 * (ttmp0 - ttmp2);
                A[j+3*v] = wj3 * (ttmp1 - ttmp3);
            }
        }
        u >>= 2;
        v <<= 2;
    }
    if(k&1) {
        for(j = 0;j<m/2; j++) {
            cmplx Ajv = A[j+(m/2)];
            A[j+(m/2)] = A[j] - Ajv;
            A[j] += Ajv;
        }
    }
}

// 卷积器
[[maybe_unused]] void tfft_convolver(int k, cmplx *A, const cmplx *w) {
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
    for (int dim = 0; dim < 3; dim++) {
        bsplineModuli[dim].resize(meshSize[dim]);
        
        // 计算B样条模数
        double factor = 2.0 * M_PI / meshSize[dim];
        for (int i = 0; i < meshSize[dim]; i++) {
            int m = i;
            // 折叠到Nyquist频率
            if (i > meshSize[dim]/2)
                m = meshSize[dim] - i;
                
            double tmp = 0.0;
            if (m != 0) {
                double w = factor * m;
                // 对4阶B样条优化
                tmp = 1.0;
                for (int j = 2; j <= splineOrder; j++) {
                    double sin_term = std::sin(j * w / 2.0) / w;
                    tmp = 4.0 * sin_term * sin_term * tmp;
                }
            }
            bsplineModuli[dim][i] = tmp;
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
 * @brief Auto-adjust PME parameters to achieve desired error tolerance
 * 
 * @param error_tolerance Error tolerance
 * @param cutoff_distance Cutoff distance
 * @param box Box dimensions
 */
void autoAdjustPMEParameters(double error_tolerance, double cutoff_distance, const double box[3]) {
    // Start with reasonable defaults
    double alpha = 0.0;
    int meshSize[3] = {0, 0, 0};
    
    // Approximate optimal alpha based on cutoff and tolerance
    alpha = std::sqrt(-std::log(2.0 * error_tolerance)) / cutoff_distance;
    
    // Calculate mesh size based on alpha and box size
    for (int i = 0; i < 3; i++) {
        // First approximation for mesh size
        meshSize[i] = static_cast<int>(std::ceil(2.0 * alpha * box[i] / M_PI));
        
        // Make sure it's a power of 2 for FFT efficiency
        int log2_size = 0;
        while ((1 << log2_size) < meshSize[i]) log2_size++;
        meshSize[i] = 1 << log2_size;
        
        // Ensure minimum size
        meshSize[i] = std::max(meshSize[i], 16);
    }
    
    platform::log(LogLevel::INFO, "Auto-adjusted PME parameters: alpha = ", alpha,
                 ", mesh size = [", meshSize[0], ",", meshSize[1], ",", meshSize[2], "]");
    
    // Set calculated parameters
    setPMEParameters(alpha, meshSize);
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
    
    // Reset the grid
    std::fill(pme_params.pmeGrid.begin(), pme_params.pmeGrid.end(), std::complex<double>(0.0, 0.0));
    
    // For each atom, spread its charge on the grid
    for (int n = 0; n < state.activeAtomCount; n++) {
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
        double x = atoms[n].x;
        double y = atoms[n].y;
        double z = atoms[n].z;
        
        // Convert to fractional coordinates
        double fracX = x / box[0];
        double fracY = y / box[1];
        double fracZ = z / box[2];
        
        // Get grid indices and weights
        double gridX = fracX * pme_params.meshSize[0];
        double gridY = fracY * pme_params.meshSize[1];
        double gridZ = fracZ * pme_params.meshSize[2];
        
        int ix = static_cast<int>(std::floor(gridX));
        int iy = static_cast<int>(std::floor(gridY));
        int iz = static_cast<int>(std::floor(gridZ));
        
        // Compute B-spline coefficients
        std::vector<double> splineX, splineY, splineZ;
        computeBSplineCoefficients(gridX - ix, pme_params.splineOrder, splineX);
        computeBSplineCoefficients(gridY - iy, pme_params.splineOrder, splineY);
        computeBSplineCoefficients(gridZ - iz, pme_params.splineOrder, splineZ);
        
        // Distribute charge to grid points
        for (int i = 0; i <= pme_params.splineOrder; i++) {
            int gridIndexX = (ix + i) % pme_params.meshSize[0];
            if (gridIndexX < 0) gridIndexX += pme_params.meshSize[0];
            
            for (int j = 0; j <= pme_params.splineOrder; j++) {
                int gridIndexY = (iy + j) % pme_params.meshSize[1];
                if (gridIndexY < 0) gridIndexY += pme_params.meshSize[1];
                
                for (int k = 0; k <= pme_params.splineOrder; k++) {
                    int gridIndexZ = (iz + k) % pme_params.meshSize[2];
                    if (gridIndexZ < 0) gridIndexZ += pme_params.meshSize[2];
                    
                    int gridIndex = gridIndexX * pme_params.meshSize[1] * pme_params.meshSize[2] +
                                    gridIndexY * pme_params.meshSize[2] + gridIndexZ;
                    
                    double weight = splineX[i] * splineY[j] * splineZ[k];
                    pme_params.pmeGrid[gridIndex] += std::complex<double>(charge * weight, 0.0);
                }
            }
        }
    }
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
    
    int log2_nx = 0;
    while ((1 << log2_nx) < nx) log2_nx++;
    
    int log2_ny = 0;
    while ((1 << log2_ny) < ny) log2_ny++;
    
    int log2_nz = 0;
    while ((1 << log2_nz) < nz) log2_nz++;
    
    // 初始化FFT权重
    fft_weights.resize(nx > ny ? (nx > nz ? nx : nz) : (ny > nz ? ny : nz));
    
    // 执行X方向FFT
    std::vector<cmplx> row_data(nx);
    tfft_init(log2_nx, fft_weights.data());
    for (int y = 0; y < ny; y++) {
        for (int z = 0; z < nz; z++) {
            // 提取一行数据
            for (int x = 0; x < nx; x++) {
                int index = x * ny * nz + y * nz + z;
                row_data[x] = pme_params.pmeGrid[index];
            }
            
            // 执行1D FFT
            tfft_fft(log2_nx, row_data.data(), fft_weights.data());
            
            // 写回结果
            for (int x = 0; x < nx; x++) {
                int index = x * ny * nz + y * nz + z;
                pme_params.pmeGrid[index] = row_data[x];
            }
        }
    }
    
    // 执行Y方向FFT
    row_data.resize(ny);
    tfft_init(log2_ny, fft_weights.data());
    for (int x = 0; x < nx; x++) {
        for (int z = 0; z < nz; z++) {
            // 提取一列数据
            for (int y = 0; y < ny; y++) {
                int index = x * ny * nz + y * nz + z;
                row_data[y] = pme_params.pmeGrid[index];
            }
            
            // 执行1D FFT
            tfft_fft(log2_ny, row_data.data(), fft_weights.data());
            
            // 写回结果
            for (int y = 0; y < ny; y++) {
                int index = x * ny * nz + y * nz + z;
                pme_params.pmeGrid[index] = row_data[y];
            }
        }
    }
    
    // 执行Z方向FFT
    row_data.resize(nz);
    tfft_init(log2_nz, fft_weights.data());
    for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny; y++) {
            // 提取一行数据
            for (int z = 0; z < nz; z++) {
                int index = x * ny * nz + y * nz + z;
                row_data[z] = pme_params.pmeGrid[index];
            }
            
            // 执行1D FFT
            tfft_fft(log2_nz, row_data.data(), fft_weights.data());
            
            // 写回结果
            for (int z = 0; z < nz; z++) {
                int index = x * ny * nz + y * nz + z;
                pme_params.pmeGrid[index] = row_data[z];
            }
        }
    }
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
    
    int log2_nx = 0;
    while ((1 << log2_nx) < nx) log2_nx++;
    
    int log2_ny = 0;
    while ((1 << log2_ny) < ny) log2_ny++;
    
    int log2_nz = 0;
    while ((1 << log2_nz) < nz) log2_nz++;
    
    // 初始化FFT权重
    fft_weights.resize(nx > ny ? (nx > nz ? nx : nz) : (ny > nz ? ny : nz));
    
    // 执行Z方向IFFT
    std::vector<cmplx> row_data(nz);
    tfft_init(log2_nz, fft_weights.data());
    for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny; y++) {
            // 提取一行数据
            for (int z = 0; z < nz; z++) {
                int index = x * ny * nz + y * nz + z;
                row_data[z] = pme_params.pmeGrid[index];
            }
            
            // 执行1D IFFT
            tfft_ifft(log2_nz, row_data.data(), fft_weights.data());
            
            // 写回结果
            for (int z = 0; z < nz; z++) {
                int index = x * ny * nz + y * nz + z;
                pme_params.pmeGrid[index] = row_data[z];
            }
        }
    }
    
    // 执行Y方向IFFT
    row_data.resize(ny);
    tfft_init(log2_ny, fft_weights.data());
    for (int x = 0; x < nx; x++) {
        for (int z = 0; z < nz; z++) {
            // 提取一列数据
            for (int y = 0; y < ny; y++) {
                int index = x * ny * nz + y * nz + z;
                row_data[y] = pme_params.pmeGrid[index];
            }
            
            // 执行1D IFFT
            tfft_ifft(log2_ny, row_data.data(), fft_weights.data());
            
            // 写回结果
            for (int y = 0; y < ny; y++) {
                int index = x * ny * nz + y * nz + z;
                pme_params.pmeGrid[index] = row_data[y];
            }
        }
    }
    
    // 执行X方向IFFT
    row_data.resize(nx);
    tfft_init(log2_nx, fft_weights.data());
    for (int y = 0; y < ny; y++) {
        for (int z = 0; z < nz; z++) {
            // 提取一行数据
            for (int x = 0; x < nx; x++) {
                int index = x * ny * nz + y * nz + z;
                row_data[x] = pme_params.pmeGrid[index];
            }
            
            // 执行1D IFFT
            tfft_ifft(log2_nx, row_data.data(), fft_weights.data());
            
            // 写回结果
            for (int x = 0; x < nx; x++) {
                int index = x * ny * nz + y * nz + z;
                double scale = static_cast<double>(nx * ny * nz);
                pme_params.pmeGrid[index] = row_data[x] / scale; // 需要重新缩放
            }
        }
    }
}

/**
 * @brief Compute energy from the PME grid after FFT
 * 
 * @param energy Output energy
 */
void computeEnergyFromGrid(double& energy) {
    const auto& box = pme_params.cutoff > 0.0 ? 
                    std::array<double, 3>{pme_params.cutoff, pme_params.cutoff, pme_params.cutoff} : 
                    std::array<double, 3>{1.0, 1.0, 1.0};
    double volume = box[0] * box[1] * box[2];
    double recipCoeff = COULOMB * 4.0 * M_PI / volume;
    
    energy = 0.0;
    
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
                double m2 = (kx_sq / (box[0] * box[0])) + 
                           (ky_sq / (box[1] * box[1])) + 
                           (kz_sq / (box[2] * box[2]));
                
                // B-spline influence function
                double bx = pme_params.bsplineModuli[0][ix];
                double by = pme_params.bsplineModuli[1][iy];
                double bz = pme_params.bsplineModuli[2][iz];
                double m2_term = m2 != 0.0 ? std::exp(-M_PI * M_PI * m2 / (pme_params.alpha * pme_params.alpha)) / m2 : 0.0;
                
                // Complex conjugate of structure factor - not actually needed in this calculation
                // std::complex<double> conjStructureFactor = std::conj(pme_params.pmeGrid[gridIndex]);
                
                // Energy contribution for this k-vector
                double term = recipCoeff * m2_term * bx * by * bz * 
                             std::norm(pme_params.pmeGrid[gridIndex]) * 0.5;
                
                energy += term;
            }
        }
    }
}

/**
 * @brief Compute reciprocal space energy using PME
 * 
 * @param state MC state
 * @param movement_only Whether to compute only for moving atoms
 * @return double Reciprocal space energy
 */
double computeReciprocalPME(model::MCState& state, bool movement_only) {
    // 不再需要box变量，直接使用state.info.box时需要时再引用
    const auto& atoms = state.atoms;
    
    // Check system neutrality
    double totalCharge = 0.0;
    for(const auto& atom : atoms) {
        totalCharge += static_cast<double>(atom.charge);
    }
    if (std::abs(totalCharge) > 1e-10) {
        throw std::runtime_error("System must be charge neutral for PME calculation");
    }
    
    // Spread charges onto grid
    spreadChargesOntoGrid(state, movement_only);
    
    // Perform forward FFT
    performFFTForward();
    
    // Compute energy from the grid
    double recipEnergy = 0.0;
    computeEnergyFromGrid(recipEnergy);
    
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
        
        for(int r2 = 0; r2 < state.activeResidueCount; r2++) {
            if(!residues[r2].active) continue;
            
            // Skip self-interactions in real space
            if(r1 == r2) continue;
            
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

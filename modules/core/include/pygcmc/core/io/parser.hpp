// modules/core/include/pygcmc/core/io/parser.hpp

#ifndef PYGCMC_CORE_IO_PARSER_HPP
#define PYGCMC_CORE_IO_PARSER_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <utility>

namespace pygcmc {
namespace core {
namespace io {

// ---------------------------
// 1. PDB Parser
// ---------------------------

/**
 * @brief 代表 PDB 文件中的原子信息
 */
struct PDBAtom {
    int serial;              ///< 原子序列号
    std::string name;        ///< 原子名称
    std::string residue;     ///< 残基名称
    int sequence;            ///< 残基序号
    double x, y, z;          ///< 原子坐标
    double charge;           ///< 原子电荷（可选）
    int type;                ///< 原子类型（可选）

    PDBAtom(int serial_, const std::string& name_, const std::string& residue_,
            int sequence_, double x_, double y_, double z_,
            double charge_ = 0.0, int type_ = 0)
        : serial(serial_), name(name_), residue(residue_),
          sequence(sequence_), x(x_), y(y_), z(z_),
          charge(charge_), type(type_) {}
};

/**
 * @brief PDB 文件解析器类
 * 
 * 负责解析 PDB 文件，提取晶胞信息和原子列表。
 * 实现位于 `pdb_parser.cpp`
 */
class PDBParser {
public:
    /**
     * @brief 解析 PDB 文件
     * @param filename PDB 文件路径
     * @return std::pair<std::vector<double>, std::vector<PDBAtom>>
     *         第一部分为晶胞信息 [a, b, c]，
     *         第二部分为原子列表
     */
    static std::pair<std::vector<double>, std::vector<PDBAtom>> parse(const std::string& filename);
};

// ---------------------------
// 2. PSF Parser
// ---------------------------

/**
 * @brief 代表 PSF 文件中的键信息
 */
struct PSFBond {
    int atom1; ///< 键的第一个原子序号
    int atom2; ///< 键的第二个原子序号
};

/**
 * @brief PSF 文件解析器类
 * 
 * 负责解析 PSF 文件，提取键等拓扑信息。
 * 实现位于 `psf_parser.cpp`
 */
struct PSFTopology {
    std::vector<PSFBond> bonds; ///< 键列表
    // 可以添加更多拓扑信息，如角度、二面角等
};

class PSFParser {
public:
    /**
     * @brief 解析 PSF 文件
     * @param filename PSF 文件路径
     * @return PSFTopology 拓扑信息
     */
    static PSFTopology parse(const std::string& filename);
};

// ---------------------------
// 3. TOP Parser
// ---------------------------

/**
 * @brief 代表 TOP 文件中的原子类型信息
 */
struct TopAtomType {
    std::string name; ///< 原子类型名称
    int type;         ///< 原子类型编号
    double charge;    ///< 原子电荷
    double mass;      ///< 原子质量
};

/**
 * @brief TOP 文件解析器类
 * 
 * 负责解析 TOP 文件，提取原子类型等拓扑信息。
 * 实现位于 `top_parser.cpp`
 */
struct Topology {
    std::vector<TopAtomType> atom_types; ///< 原子类型列表
    // 可以添加更多拓扑信息，如键、角度等
};

class TopParser {
public:
    /**
     * @brief 解析 TOP 文件
     * @param filename TOP 文件路径
     * @return Topology 拓扑信息
     */
    static Topology parse(const std::string& filename);
};

// ---------------------------
// 4. ITP Parser
// ---------------------------

/**
 * @brief 代表 ITP 文件中的原子信息
 */
struct ITPAtom {
    std::string name;    ///< 原子名称
    std::string type;    ///< 原子类型
    int resid;           ///< 残基序号
    std::string resname; ///< 残基名称
    double charge;       ///< 原子电荷
};

/**
 * @brief ITP 文件解析器类
 * 
 * 负责解析 ITP 文件，提取分子信息等。
 * 实现位于 `itp_parser.cpp`
 */
class ITPParser {
public:
    /**
     * @brief 解析 ITP 文件
     * @param filename ITP 文件路径
     * @return std::vector<ITPAtom> 原子列表
     */
    static std::vector<ITPAtom> parse(const std::string& filename);
};

// ---------------------------
// 5. Force Field Parser
// ---------------------------

/**
 * @brief 代表势能参数中的非键参数
 */
struct ForceFieldPair {
    double sigma;  ///< Lennard-Jones 势的 sigma 参数
    double epsilon;///< Lennard-Jones 势的 epsilon 参数
};

/**
 * @brief 势能文件解析器类
 * 
 * 负责解析势能文件，提取非键参数和修正参数。
 * 实现位于 `ff_parser.cpp`
 */
class FFParser {
public:
    /**
     * @brief 解析势能文件
     * @param filename 势能文件路径
     * @return std::pair<std::unordered_map<std::string, ForceFieldPair>, 
     *                 std::unordered_map<std::pair<std::string, std::string>, ForceFieldPair>>
     *         第一部分为非键参数映射，
     *         第二部分为修正参数映射（键为原子类型对）
     */
    static std::pair<
        std::unordered_map<std::string, ForceFieldPair>, 
        std::unordered_map<std::pair<std::string, std::string>, ForceFieldPair>
    > parse(const std::string& filename);
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_PARSER_HPP

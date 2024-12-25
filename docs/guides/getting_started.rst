Getting Started
===============

欢迎使用 `pyGCMC`！本指南将帮助您快速构建、测试和运行最基础的功能。

安装依赖
---------

确保您已经安装了以下依赖：

- **CMake**：版本 3.15 或更高
- **编译器**：支持C++17的编译器（如 `g++` 或 `clang++`）
- **Git**：用于获取项目代码和依赖
- **Catch2**：单元测试框架，通过CMake的FetchContent自动获取

构建项目
--------

1. **克隆项目仓库**

    ```bash
    git clone https://github.com/您的用户名/pyGCMC.git
    cd pyGCMC
    ```

2. **创建构建目录并配置项目**

    ```bash
    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    ```

    - 您可以将 `-DCMAKE_BUILD_TYPE` 设置为 `Debug` 以启用调试信息。
    - CMake 会自动下载并配置Catch2作为测试框架。

3. **编译项目**

    ```bash
    make -j$(nproc)
    ```

    - `-j$(nproc)` 选项会利用所有可用的CPU核心加速编译过程。

运行单元测试
------------

在构建目录下，运行以下命令执行单元测试：

```bash
ctest

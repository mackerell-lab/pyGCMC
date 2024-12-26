# modules/core/CMakeLists.txt

# 添加核心静态库
add_library(pygcmc_core STATIC
    src/system.cpp
    src/box.cpp
    src/force_field.cpp
    src/particles.cpp
)

# 启用位置无关代码
set_target_properties(pygcmc_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

# 包含目录
target_include_directories(pygcmc_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# 创建 Python 绑定模块
pybind11_add_module(pyGCMC MODULE src/bindings.cpp)

# 链接核心库到绑定模块
target_link_libraries(pyGCMC PRIVATE pygcmc_core)

# 设置模块输出路径为构建目录的 `python` 子目录
set_target_properties(pyGCMC PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
)

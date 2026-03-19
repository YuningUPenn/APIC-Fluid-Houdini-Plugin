cmake_minimum_required(VERSION 3.14)
project(apic_fluid VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -------------------------------------------------------
# Houdini path configuration
#
# Priority (high to low):
#   1. cmake -DHOUDINI_INSTALL_DIR="C:/Program Files/Side Effects Software/Houdini 21.0.xxx" ..
#   2. environment variable HFS (auto-set by Houdini installer)
#   3. 下面的 HOUDINI_FALLBACK_PATHS 列表（常见默认路径）
#
# No need to modify this file once path is found.
# -------------------------------------------------------

# Can be passed from CMakeSettings.json or command line
set(HOUDINI_INSTALL_DIR "" CACHE PATH
    "Houdini 安装根目录，例如 C:/Program Files/Side Effects Software/Houdini 21.0.741")

# Add your Houdini version path here if not listed below
set(HOUDINI_FALLBACK_PATHS
    # === ADD YOUR VERSION HERE ===
    # "C:/Program Files/Side Effects Software/Houdini 21.0.800"

    # === Known versions (no need to edit) ===
    "C:/Program Files/Side Effects Software/Houdini 21.0.741"
    "C:/Program Files/Side Effects Software/Houdini 21.0.718"
    "C:/Program Files/Side Effects Software/Houdini 21.0.506"
    "C:/Program Files/Side Effects Software/Houdini 20.5.584"
    "C:/Program Files/Side Effects Software/Houdini 20.5.550"
    "C:/Program Files/Side Effects Software/Houdini 20.0.688"
)

# Resolve final HFS path
if(HOUDINI_INSTALL_DIR)
    set(_HFS "${HOUDINI_INSTALL_DIR}")
elseif(DEFINED ENV{HFS} AND NOT "$ENV{HFS}" STREQUAL "")
    set(_HFS "$ENV{HFS}")
else()
    foreach(_path IN LISTS HOUDINI_FALLBACK_PATHS)
        if(EXISTS "${_path}/toolkit/cmake")
            set(_HFS "${_path}")
            break()
        endif()
    endforeach()
endif()

if(NOT _HFS OR NOT EXISTS "${_HFS}/toolkit/cmake")
    message(FATAL_ERROR
        "Houdini HDK not found。请用以下方式之一指定路径：\n"
        "  方式1: cmake -DHOUDINI_INSTALL_DIR=\"C:/Program Files/Side Effects Software/Houdini 21.0.xxx\" ..\n"
        "  方式2: 在 CMakeSettings.json 的 cmakeCommandArgs 里加入 -DHOUDINI_INSTALL_DIR=...\n"
        "  方式3: 设置环境变量 HFS 指向 Houdini 安装目录\n"
        "Current _HFS = '${_HFS}'"
    )
endif()

message(STATUS "Houdini HDK path: ${_HFS}")
list(APPEND CMAKE_PREFIX_PATH "${_HFS}/toolkit/cmake")
find_package(Houdini REQUIRED)

# -------------------------------------------------------
# Eigen (header-only, auto-detected)
#
# Priority:
#   1. vcpkg 集成（VS 打开 folder 时自动生效，推荐）
#   2. cmake -DEIGEN3_ROOT=C:/path/to/eigen ..
#   3. third_party/eigen-3.4.0/ in project (auto-detected after extraction)
#   4. system install (Linux apt / macOS brew)
# -------------------------------------------------------
set(EIGEN3_ROOT "" CACHE PATH "Eigen3 根目录（可选，也可用 vcpkg 或放 third_party/）")

if(EIGEN3_ROOT)
    list(APPEND CMAKE_PREFIX_PATH "${EIGEN3_ROOT}")
endif()

# Check for Eigen in third_party/
set(_EIGEN_LOCAL "${CMAKE_SOURCE_DIR}/third_party/eigen-3.4.0")
if(EXISTS "${_EIGEN_LOCAL}/Eigen/Core")
    message(STATUS "Using bundled Eigen: ${_EIGEN_LOCAL}")
    set(EIGEN3_INCLUDE_DIR "${_EIGEN_LOCAL}")
    # Create import target for consistent interface
    add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
    target_include_directories(Eigen3::Eigen INTERFACE "${_EIGEN_LOCAL}")
else()
    find_package(Eigen3 3.4 REQUIRED NO_MODULE)
    message(STATUS "Using system/vcpkg Eigen: ${EIGEN3_INCLUDE_DIRS}")
endif()

# -------------------------------------------------------
# Core APIC library (pure C++, no Houdini dependency)
# -------------------------------------------------------
add_library(apic_core STATIC
    src/ApicParticleSet.cpp
    src/ApicGrid.cpp
    src/ApicSolver.cpp
    src/PressureSolver.cpp
    src/CollisionHandler.cpp
)

target_include_directories(apic_core PUBLIC include)
target_link_libraries(apic_core PUBLIC Eigen3::Eigen)

# MSVC-specific settings
if(MSVC)
    target_compile_options(apic_core PRIVATE /W3 /wd4244 /wd4267 /wd4305)
    # Multi-processor compilation
    target_compile_options(apic_core PRIVATE /MP)
    # Disable MSVC min/max macros (conflict with Eigen)
    target_compile_definitions(apic_core PUBLIC NOMINMAX WIN32_LEAN_AND_MEAN)
endif()

# -------------------------------------------------------
# Houdini SOP DSO
# -------------------------------------------------------
add_library(SOP_ApicFluid SHARED
    src/SOP_ApicFluidSolver.cpp
)

target_link_libraries(SOP_ApicFluid PRIVATE apic_core Houdini)
target_include_directories(SOP_ApicFluid PRIVATE include)
houdini_configure_target(SOP_ApicFluid)

if(MSVC)
    target_compile_options(SOP_ApicFluid PRIVATE /MP /wd4244 /wd4267)
    target_compile_definitions(SOP_ApicFluid PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
endif()

# -------------------------------------------------------
# Install: copy .dll/.so to Houdini dso directory
# cmake --build . --target install
# -------------------------------------------------------
set(HOUDINI_DSO_DIR "$ENV{USERPROFILE}/Documents/houdini21.0/dso"
    CACHE PATH "Houdini 用户 dso 目录（Windows 默认）")

install(TARGETS SOP_ApicFluid
    RUNTIME DESTINATION "${HOUDINI_DSO_DIR}"  # .dll
    LIBRARY DESTINATION "${HOUDINI_DSO_DIR}"  # .so
)

# -------------------------------------------------------
# Standalone unit tests (no Houdini required)
# Runnable directly from VS Test Explorer
# -------------------------------------------------------
enable_testing()

foreach(_test p2g_conservation angular_momentum pressure_solver)
    add_executable(test_${_test} tests/test_${_test}.cpp)
    target_link_libraries(test_${_test} PRIVATE apic_core)
    if(MSVC)
        target_compile_options(test_${_test} PRIVATE /W3 /wd4244 /wd4267)
        target_compile_definitions(test_${_test} PRIVATE NOMINMAX)
    endif()
    add_test(NAME ${_test} COMMAND test_${_test})
endforeach()

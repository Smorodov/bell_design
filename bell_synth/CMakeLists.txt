# =============================================================================
# CMakeLists.txt - СБОРКА ПРОЕКТА С FFTW
# =============================================================================

cmake_minimum_required(VERSION 3.14)
project(BellSynthesizer VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---- ВСЕ БИНАРНИКИ В КОРЕНЬ СБОРКИ ----
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# ---- ПОДКЛЮЧЕНИЕ ЗАВИСИМОСТЕЙ ЧЕРЕЗ INTERNET ----
include(FetchContent)

# nlohmann/json (если ещё не подключен)
if(NOT TARGET nlohmann_json::nlohmann_json)
    FetchContent_Declare(
        json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.3
    )
    FetchContent_MakeAvailable(json)
endif()

# fmt (если ещё не подключен)
if(NOT TARGET fmt::fmt)
    FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG 10.2.1
    )
    FetchContent_MakeAvailable(fmt)
endif()

# ---- ОПЦИИ ----
option(USE_OPENMP "Use OpenMP parallelization" ON)

# ---- ТИП СБОРКИ ----
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# ---- ОПТИМИЗАЦИИ ----
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -mtune=native -ffast-math")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -Wall -Wextra -pedantic")
endif()

# ---- ПОИСК EIGEN3 ----
find_package(Eigen3 REQUIRED)

# ---- ПОИСК FFTW (БЕЗ ЖЕСТКИХ ПУТЕЙ) ----
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED fftw3)

if(NOT FFTW3_FOUND)
    find_path(FFTW3_INCLUDE_DIRS
        NAMES fftw3.h
        HINTS /usr/include /usr/local/include
    )
    find_library(FFTW3_LIBRARIES
        NAMES fftw3 libfftw3 fftw3-3
        HINTS /usr/lib /usr/local/lib
    )
endif()

if(FFTW3_INCLUDE_DIRS AND FFTW3_LIBRARIES)
    set(FFTW3_FOUND TRUE)
    message(STATUS "FFTW3 found:")
    message(STATUS "  Include: ${FFTW3_INCLUDE_DIRS}")
    message(STATUS "  Library: ${FFTW3_LIBRARIES}")
else()
    message(FATAL_ERROR "FFTW3 not found! Install: sudo apt install libfftw3-dev")
endif()

# ---- OPENMP ----
if(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        set(USE_OPENMP_FLAGS ${OpenMP_CXX_FLAGS})
        set(USE_OPENMP_LIBS ${OpenMP_CXX_LIBRARIES})
        message(STATUS "OpenMP found, enabling parallelization")
    else()
        message(STATUS "OpenMP not found")
        set(USE_OPENMP_FLAGS "")
        set(USE_OPENMP_LIBS "")
    endif()
endif()

# =============================================================================
# ИСПОЛНЯЕМЫЙ ФАЙЛ
# =============================================================================

add_executable(bell_synth
    src/main.cpp
    src/bell_synthesis.cpp
)

# ---- EIGEN3 ----
if(TARGET Eigen3::Eigen)
    target_link_libraries(bell_synth PRIVATE Eigen3::Eigen)
else()
    target_include_directories(bell_synth PRIVATE ${EIGEN3_INCLUDE_DIRS})
endif()

# ---- FFTW ----
target_include_directories(bell_synth PRIVATE ${FFTW3_INCLUDE_DIRS})
target_link_libraries(bell_synth PRIVATE ${FFTW3_LIBRARIES})
target_compile_definitions(bell_synth PRIVATE USE_FFTW)

# ---- JSON & FMT ----
target_link_libraries(bell_synth PRIVATE
    nlohmann_json::nlohmann_json
    fmt::fmt
)

target_link_libraries(bell_synth PRIVATE ${USE_OPENMP_LIBS})
target_compile_options(bell_synth PRIVATE ${USE_OPENMP_FLAGS})

if(WIN32)
    target_compile_definitions(bell_synth PRIVATE _USE_MATH_DEFINES)
endif()

target_include_directories(bell_synth PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# ---- ИНФОРМАЦИЯ ----
message(STATUS "")
message(STATUS "========================================")
message(STATUS "BellSynthesizer Configuration:")
message(STATUS "  Build Type:     ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ Standard:   ${CMAKE_CXX_STANDARD}")
message(STATUS "  Eigen3 Include: ${EIGEN3_INCLUDE_DIRS}")
message(STATUS "  OpenMP:         ${USE_OPENMP_FLAGS}")
message(STATUS "  FFTW:           Enabled")
if(TARGET nlohmann_json::nlohmann_json)
    message(STATUS "  JSON:           FOUND")
else()
    message(STATUS "  JSON:           NOT FOUND")
endif()
if(TARGET fmt::fmt)
    message(STATUS "  FMT:            FOUND")
else()
    message(STATUS "  FMT:            NOT FOUND")
endif()
message(STATUS "========================================")
message(STATUS "")
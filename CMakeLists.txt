cmake_minimum_required(VERSION 3.10)
project(BellOptimizer VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- ВСЕ БИНАРНИКИ В КОРЕНЬ СБОРКИ ----
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# ---- ПОДКЛЮЧЕНИЕ ЗАВИСИМОСТЕЙ ЧЕРЕЗ INTERNET ----
include(FetchContent)

# nlohmann/json
FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
)
FetchContent_MakeAvailable(json)

# fmt
FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 10.2.1
)
FetchContent_MakeAvailable(fmt)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -Wall -Wextra -Wno-unused-parameter")
endif()

# ============================================================
# 1. ПОИСК GSL
# ============================================================
find_package(GSL REQUIRED)

# ============================================================
# 2. ПОИСК NLopt
# ============================================================
find_package(NLopt REQUIRED)

# ============================================================
# 3. ИСХОДНЫЙ КОД
# ============================================================
set(SOURCES 
    src/main.cpp 
    src/bell_generator.cpp
    src/bell_optimizer.cpp
)

# ============================================================
# 4. СОЗДАНИЕ ИСПОЛНЯЕМОГО ФАЙЛА
# ============================================================
add_executable(bell_opt ${SOURCES})

target_compile_definitions(bell_opt PRIVATE NLOPT_FOUND)

target_link_libraries(bell_opt PRIVATE
    ${GSL_LIBRARIES}
    NLopt::nlopt
    nlohmann_json::nlohmann_json
    fmt::fmt
)

target_include_directories(bell_opt PRIVATE
    ${GSL_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}
)

# ============================================================
# 5. КОПИРОВАНИЕ DLL ДЛЯ WINDOWS
# ============================================================
if(WIN32)
    find_file(GSL_DLL libgsl.dll)
    find_file(GSL_CBLAS_DLL libgslcblas.dll)
    find_file(NLOPT_DLL libnlopt.dll)
    
    if(GSL_DLL AND GSL_CBLAS_DLL)
        add_custom_command(TARGET bell_opt POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${GSL_DLL} ${GSL_CBLAS_DLL} ${CMAKE_BINARY_DIR}
            COMMENT "Copying GSL DLLs"
        )
    endif()
    
    if(NLOPT_DLL)
        add_custom_command(TARGET bell_opt POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${NLOPT_DLL} ${CMAKE_BINARY_DIR}
            COMMENT "Copying NLopt DLL"
        )
    endif()
endif()

# ============================================================
# 6. ИНФОРМАЦИЯ
# ============================================================
message(STATUS "")
message(STATUS "========================================")
message(STATUS "BellOptimizer Configuration:")
message(STATUS "  Build Type:     ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ Standard:   ${CMAKE_CXX_STANDARD}")
message(STATUS "  GSL:            ${GSL_INCLUDE_DIRS}")
if(TARGET NLopt::nlopt)
    message(STATUS "  NLopt:          FOUND")
else()
    message(STATUS "  NLopt:          NOT FOUND")
endif()
if(TARGET nlohmann_json::nlohmann_json)
    message(STATUS "  JSON:           FOUND (FetchContent)")
else()
    message(STATUS "  JSON:           NOT FOUND")
endif()
if(TARGET fmt::fmt)
    message(STATUS "  FMT:            FOUND (FetchContent)")
else()
    message(STATUS "  FMT:            NOT FOUND")
endif()
message(STATUS "========================================")
message(STATUS "")

add_subdirectory(bell_synth)
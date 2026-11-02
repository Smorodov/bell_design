cmake_minimum_required(VERSION 3.16)

# Detached research module. It is intentionally excluded from the normal build.
# It consumes the validated BellOptimizer core from the parent project but does
# not add Chinese-bell code to the production library.

add_library(bell_chinese_research STATIC
    src/excitation/chinese_bell_strike.cpp
    src/excitation/chinese_bell_objective.cpp
    src/chinese_modal_evaluator.cpp
)
target_include_directories(bell_chinese_research PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/core
)
target_link_libraries(bell_chinese_research PUBLIC bell_core)

if(MSVC)
    target_compile_options(bell_chinese_research PRIVATE /W4)
else()
    target_compile_options(bell_chinese_research PRIVATE -Wall -Wextra -Wpedantic)
endif()

function(add_chinese_research_test name source)
    add_executable(${name} tests/${source})
    target_link_libraries(${name} PRIVATE bell_chinese_research)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

if(NOT BUILD_TESTING)
    enable_testing()
endif()

add_chinese_research_test(chinese_bell_strike_surface_normal_contract chinese_bell_strike_surface_normal_contract.cpp)
add_chinese_research_test(chinese_bell_strike_mesh_convergence chinese_bell_strike_mesh_convergence.cpp)
add_chinese_research_test(chinese_bell_strike_angular_convergence chinese_bell_strike_angular_convergence.cpp)
add_chinese_research_test(chinese_bell_contact_surface_metric_audit chinese_bell_contact_surface_metric_audit.cpp)
add_chinese_research_test(chinese_bell_temporal_impulse_audit chinese_bell_temporal_impulse_audit.cpp)
add_chinese_research_test(chinese_modal_family_contract chinese_modal_family_contract.cpp)
add_chinese_research_test(chinese_bell_objective_contract chinese_bell_objective_contract.cpp)

cmake_minimum_required(VERSION 3.16)
project(BellOptimizer VERSION 1.2.96 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # GCC 14 can emit a false-positive -Wfree-nonheap-object from the
    # inlined header-only PCG allocator path at -O2/-O3. Keep the remaining
    # warning classes enabled while suppressing only this known diagnostic.
    add_compile_options(-Wno-free-nonheap-object)
endif()

option(BUILD_SMOKE_TESTS "Build validation and smoke tests" ON)
set(CCX_EXECUTABLE "ccx" CACHE STRING "CalculiX executable used by the reference modal solver")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(BELL_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(BELL_CORE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/core)

add_library(bell_core STATIC
    src/core/bell_generator.cpp
    src/core/bell_optimizer.cpp
    src/core/modal_tracking.cpp
    src/models/bell_model_profile.cpp
    src/math/coordinate_morphing.cpp
    src/math/modal_coordinate_morphing.cpp
    src/math/cylindrical_strain.cpp
    src/fem/cylindrical_q4_element.cpp
    src/fem/cylindrical_global_assembly.cpp
    src/fem/cylindrical_global_modes.cpp
    src/fem/cylindrical_bell_model.cpp
    src/fem/cylindrical_modal_evaluator.cpp
    src/fem/cylindrical_modal_tracking.cpp
    src/fem/internal_modal_solver.cpp
    src/excitation/localized_strike.cpp
    src/excitation/chinese_bell_strike.cpp
    src/excitation/meridional_influence.cpp
    src/ccx/ccx_solver.cpp
    src/ccx/ccx_input_writer.cpp
)
target_include_directories(bell_core PUBLIC ${BELL_SRC} ${BELL_CORE_SRC})
target_compile_definitions(bell_core PRIVATE BELL_CCX_EXECUTABLE="${CCX_EXECUTABLE}")

if(MSVC)
    target_compile_options(bell_core PRIVATE /W4)
else()
    target_compile_options(bell_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(bell
    src/main.cpp
    src/commands/optimize.cpp
    src/commands/synthesize.cpp
    src/commands/synthesize_geometry.cpp
    src/synthesis/bell_synthesis.cpp
)
target_include_directories(bell PRIVATE ${BELL_SRC} ${BELL_CORE_SRC} ${BELL_SRC}/synthesis)
target_link_libraries(bell PRIVATE bell_core)


if(BUILD_SMOKE_TESTS)
    enable_testing()

    # Experimental diagnostics are deliberately kept out of the production library.
    add_library(bell_experimental STATIC
        src/core/modal_state_tracker.cpp
        src/fem/cylindrical_modal_field.cpp
        src/excitation/chinese_bell_objective.cpp
    )
    add_library(bell_validation STATIC
        src/validation/nonaxisymmetric_reference.cpp
    )
    target_include_directories(bell_validation PUBLIC ${BELL_SRC} ${BELL_CORE_SRC})
    target_link_libraries(bell_validation PUBLIC bell_core)
    if(MSVC)
        target_compile_options(bell_validation PRIVATE /W4)
    else()
        target_compile_options(bell_validation PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    target_include_directories(bell_experimental PUBLIC ${BELL_SRC} ${BELL_CORE_SRC})
    target_link_libraries(bell_experimental PUBLIC bell_core)
    if(MSVC)
        target_compile_options(bell_experimental PRIVATE /W4)
    else()
        target_compile_options(bell_experimental PRIVATE -Wall -Wextra -Wpedantic)
    endif()

    add_executable(direct_elliptical_fem_evaluation tests/direct_elliptical_fem_evaluation.cpp)
    target_link_libraries(direct_elliptical_fem_evaluation PRIVATE bell_core)
    add_executable(rcm_ordering_audit tests/rcm_ordering_audit.cpp)
    target_link_libraries(rcm_ordering_audit PRIVATE bell_core)

    add_executable(sparse_matrix_contract tests/sparse_matrix_contract.cpp)
    target_include_directories(sparse_matrix_contract PRIVATE ${BELL_SRC})

    add_executable(sparse_eigensolver_contract tests/sparse_eigensolver_contract.cpp)
    add_executable(rcm_initial_vector_contract tests/rcm_initial_vector_contract.cpp)
    target_include_directories(sparse_eigensolver_contract PRIVATE ${BELL_SRC})
    target_include_directories(rcm_initial_vector_contract PRIVATE ${BELL_SRC})

    add_executable(matrix_assembly_audit tests/matrix_assembly_audit.cpp)
    target_link_libraries(matrix_assembly_audit PRIVATE bell_core)

    add_executable(sparse_fem_solver_audit tests/sparse_fem_solver_audit.cpp)
    target_link_libraries(sparse_fem_solver_audit PRIVATE bell_core)

    add_executable(sparse_dense_mode_comparison tests/sparse_dense_mode_comparison.cpp)
    target_link_libraries(sparse_dense_mode_comparison PRIVATE bell_core)

    add_executable(sparse_solver_robustness tests/sparse_solver_robustness.cpp)
    target_link_libraries(sparse_solver_robustness PRIVATE bell_core)

    add_executable(sparse_solver_production_contract tests/sparse_solver_production_contract.cpp)
    target_include_directories(sparse_solver_production_contract PRIVATE ${BELL_SRC})

    add_executable(sparse_solver_scaling tests/sparse_solver_scaling.cpp)
    add_executable(sparse_solver_scaling_contract tests/sparse_solver_scaling_contract.cpp)
    target_link_libraries(sparse_solver_scaling_contract PRIVATE bell_core)
    target_include_directories(sparse_solver_scaling PRIVATE ${BELL_SRC})

    add_executable(sparse_conditioning_audit tests/sparse_conditioning_audit.cpp)
    add_executable(rcm_solver_benchmark tests/rcm_solver_benchmark.cpp)
    target_link_libraries(rcm_solver_benchmark PRIVATE bell_core)
    add_executable(sparse_pruning_audit tests/sparse_pruning_audit.cpp)
    target_include_directories(sparse_pruning_audit PRIVATE ${BELL_SRC})
    target_include_directories(sparse_conditioning_audit PRIVATE ${BELL_SRC})

    add_executable(eigenpair_integrity tests/eigenpair_integrity.cpp)
    target_link_libraries(eigenpair_integrity PRIVATE bell_core)

    add_executable(low_frequency_mode_policy tests/low_frequency_mode_policy.cpp)
    target_link_libraries(low_frequency_mode_policy PRIVATE bell_core)

    add_executable(rotational_reference_frame_audit tests/rotational_reference_frame_audit.cpp)
    add_executable(meridional_influence_audit tests/meridional_influence_audit.cpp)

    add_executable(modal_physics_validation tests/modal_physics_validation.cpp)
    target_link_libraries(modal_physics_validation PRIVATE bell_core)

    add_executable(modal_objective_convergence tests/modal_objective_convergence.cpp)
    target_link_libraries(modal_objective_convergence PRIVATE bell_core)

    add_executable(modal_family_tracking tests/modal_family_tracking.cpp)
    target_link_libraries(modal_family_tracking PRIVATE bell_core)

    add_executable(optimizer_shape_scale_audit tests/optimizer_shape_scale_audit.cpp)
    target_link_libraries(optimizer_shape_scale_audit PRIVATE bell_core)

    add_executable(profile_spline_convergence tests/profile_spline_convergence.cpp)
    target_link_libraries(profile_spline_convergence PRIVATE bell_core)

    add_executable(castability_geometry_contract tests/castability_geometry_contract.cpp)
    target_link_libraries(castability_geometry_contract PRIVATE bell_core)

    add_executable(acoustic_modal_validation tests/acoustic_modal_validation.cpp)
    target_link_libraries(acoustic_modal_validation PRIVATE bell_core)

    add_executable(acoustic_radiation_audit tests/acoustic_radiation_audit.cpp)
    target_link_libraries(acoustic_radiation_audit PRIVATE bell_core)

    add_executable(acoustic_multipole_audit tests/acoustic_multipole_audit.cpp)
    target_link_libraries(acoustic_multipole_audit PRIVATE bell_core)

    add_executable(acoustic_pressure_contract tests/acoustic_pressure_contract.cpp)
    target_include_directories(acoustic_pressure_contract PRIVATE ${BELL_SRC} ${BELL_CORE_SRC})

    add_executable(acoustic_pressure_scaling_audit tests/acoustic_pressure_scaling_audit.cpp)
    target_include_directories(acoustic_pressure_scaling_audit PRIVATE ${BELL_SRC} ${BELL_CORE_SRC})

    add_executable(acoustic_pressure_phase_contract tests/acoustic_pressure_phase_contract.cpp)
    add_executable(acoustic_amplitude_contract tests/acoustic_amplitude_contract.cpp)
    target_include_directories(acoustic_amplitude_contract PRIVATE ${BELL_SRC} ${BELL_CORE_SRC})
    target_include_directories(acoustic_pressure_phase_contract PRIVATE ${BELL_SRC} ${BELL_CORE_SRC})

    add_executable(rotational_excitation_contract tests/rotational_excitation_contract.cpp)
    add_executable(excitation_model_audit tests/excitation_model_audit.cpp)
    add_executable(localized_strike_spectrum_audit tests/localized_strike_spectrum_audit.cpp)
    target_link_libraries(localized_strike_spectrum_audit PRIVATE bell_core)
    add_executable(localized_strike_direction_contract tests/localized_strike_direction_contract.cpp)
    target_link_libraries(localized_strike_direction_contract PRIVATE bell_core)

    add_executable(localized_strike_mesh_convergence tests/localized_strike_mesh_convergence.cpp)
    target_link_libraries(localized_strike_mesh_convergence PRIVATE bell_core)

    add_executable(inner_normal_offset_contract tests/inner_normal_offset_contract.cpp)
    target_link_libraries(inner_normal_offset_contract PRIVATE bell_core)

    add_executable(internal_fem_compare tests/internal_fem_compare.cpp)
    target_link_libraries(internal_fem_compare PRIVATE bell_core)
    target_compile_definitions(internal_fem_compare PRIVATE BELL_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

    add_executable(ccx_current_geometry tests/ccx_current_geometry.cpp)
    target_link_libraries(ccx_current_geometry PRIVATE bell_core)

    add_executable(mesh_convergence tests/mesh_convergence.cpp)
    target_link_libraries(mesh_convergence PRIVATE bell_core)

    add_executable(thickness_mesh_convergence tests/thickness_mesh_convergence.cpp)
    target_link_libraries(thickness_mesh_convergence PRIVATE bell_core)

    add_executable(input_safety_contract tests/input_safety_contract.cpp)
    target_link_libraries(input_safety_contract PRIVATE bell_core)

    add_executable(tradition_assessment tests/tradition_assessment.cpp)
    target_link_libraries(tradition_assessment PRIVATE bell_core)

    add_executable(model_layer_contract tests/model_layer_contract.cpp)
    target_link_libraries(model_layer_contract PRIVATE bell_core)

    add_executable(modal_influence_contract tests/modal_influence_contract.cpp)
    target_link_libraries(modal_influence_contract PRIVATE bell_core)

    add_executable(coordinate_morphing_contract tests/coordinate_morphing_contract.cpp)
    target_link_libraries(coordinate_morphing_contract PRIVATE bell_core)

    add_executable(mechanical_coordinate_morphing_contract tests/mechanical_coordinate_morphing_contract.cpp)
    target_link_libraries(mechanical_coordinate_morphing_contract PRIVATE bell_core)

    add_executable(coordinate_morphing_ellipticity_sweep tests/coordinate_morphing_ellipticity_sweep.cpp)
    target_link_libraries(coordinate_morphing_ellipticity_sweep PRIVATE bell_core)



    add_executable(cylindrical_strain_contract tests/cylindrical_strain_contract.cpp)
    target_link_libraries(cylindrical_strain_contract PRIVATE bell_core)

    add_executable(cylindrical_q4_element_contract tests/cylindrical_q4_element_contract.cpp)
    target_link_libraries(cylindrical_q4_element_contract PRIVATE bell_core)
    add_executable(cylindrical_global_assembly_contract tests/cylindrical_global_assembly_contract.cpp)
    target_link_libraries(cylindrical_global_assembly_contract PRIVATE bell_core)

    add_executable(cylindrical_bell_geometry_path_contract tests/cylindrical_bell_geometry_path_contract.cpp)
    target_link_libraries(cylindrical_bell_geometry_path_contract PRIVATE bell_core)

    add_executable(cylindrical_global_modes_contract tests/cylindrical_global_modes_contract.cpp)
    target_link_libraries(cylindrical_global_modes_contract PRIVATE bell_core)
    add_executable(cylindrical_internal_solver_compare tests/cylindrical_internal_solver_compare.cpp)
    target_link_libraries(cylindrical_internal_solver_compare PRIVATE bell_core)
    add_executable(cylindrical_modal_family_tracking tests/cylindrical_modal_family_tracking.cpp)
    target_link_libraries(cylindrical_modal_family_tracking PRIVATE bell_core)
    add_executable(modal_subspace_similarity_contract tests/modal_subspace_similarity_contract.cpp)
    target_link_libraries(modal_subspace_similarity_contract PRIVATE bell_core)
    add_executable(modal_subspace_conditioning_contract tests/modal_subspace_conditioning_contract.cpp)
    target_link_libraries(modal_subspace_conditioning_contract PRIVATE bell_core)
    add_executable(cylindrical_bell_geometry_modal_tracking tests/cylindrical_bell_geometry_modal_tracking.cpp)
    target_link_libraries(cylindrical_bell_geometry_modal_tracking PRIVATE bell_core)
    add_executable(cylindrical_bell_geometry_convergence tests/cylindrical_bell_geometry_convergence.cpp)
    target_link_libraries(cylindrical_bell_geometry_convergence PRIVATE bell_core)
    add_executable(cylindrical_harmonic_family_continuation tests/cylindrical_harmonic_family_continuation.cpp)
    target_link_libraries(cylindrical_harmonic_family_continuation PRIVATE bell_core)
    add_executable(cylindrical_ellipticity_convergence tests/cylindrical_ellipticity_convergence.cpp)
    target_link_libraries(cylindrical_ellipticity_convergence PRIVATE bell_core)
    add_executable(cylindrical_unified_ellipticity_audit tests/cylindrical_unified_ellipticity_audit.cpp)
    target_link_libraries(cylindrical_unified_ellipticity_audit PRIVATE bell_core)


    add_executable(ellipticity_search_policy_contract tests/ellipticity_search_policy_contract.cpp)
    target_link_libraries(ellipticity_search_policy_contract PRIVATE bell_core)

    add_executable(nonaxisymmetric_reference_mesh_contract tests/nonaxisymmetric_reference_mesh_contract.cpp)
    target_link_libraries(nonaxisymmetric_reference_mesh_contract PRIVATE bell_validation)

    add_executable(bell_core_smoke tests/bell_core_smoke.cpp)
    target_link_libraries(bell_core_smoke PRIVATE bell_core)
    target_compile_definitions(bell_core_smoke PRIVATE BELL_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

    add_test(NAME bell_core_smoke COMMAND bell_core_smoke)
    add_test(NAME rcm_ordering_audit COMMAND rcm_ordering_audit)
    add_test(NAME sparse_matrix_contract COMMAND sparse_matrix_contract)
    add_test(NAME sparse_eigensolver_contract COMMAND sparse_eigensolver_contract)
    add_test(NAME rcm_initial_vector_contract COMMAND rcm_initial_vector_contract)
    add_test(NAME matrix_assembly_audit COMMAND matrix_assembly_audit)
    add_test(NAME sparse_fem_solver_audit COMMAND sparse_fem_solver_audit)
    add_test(NAME sparse_dense_mode_comparison COMMAND sparse_dense_mode_comparison)
    add_test(NAME sparse_solver_robustness COMMAND sparse_solver_robustness)
    add_test(NAME sparse_solver_production_contract COMMAND sparse_solver_production_contract)
    add_test(NAME sparse_solver_scaling COMMAND sparse_solver_scaling)
    add_test(NAME sparse_solver_scaling_contract COMMAND sparse_solver_scaling_contract)
    add_test(NAME sparse_conditioning_audit COMMAND sparse_conditioning_audit)
    add_test(NAME rcm_solver_benchmark COMMAND rcm_solver_benchmark)
    add_test(NAME sparse_pruning_audit COMMAND sparse_pruning_audit)
    add_test(NAME eigenpair_integrity COMMAND eigenpair_integrity)
    add_test(NAME low_frequency_mode_policy COMMAND low_frequency_mode_policy)
    add_test(NAME rotational_reference_frame_audit COMMAND rotational_reference_frame_audit)
    add_test(NAME modal_physics_validation COMMAND modal_physics_validation)
    add_test(NAME modal_objective_convergence COMMAND modal_objective_convergence)
    add_test(NAME modal_family_tracking COMMAND modal_family_tracking)
    add_test(NAME optimizer_shape_scale_audit COMMAND optimizer_shape_scale_audit)
    add_test(NAME profile_spline_convergence COMMAND profile_spline_convergence)
    add_test(NAME castability_geometry_contract COMMAND castability_geometry_contract)
    add_test(NAME acoustic_modal_validation COMMAND acoustic_modal_validation)
    add_test(NAME acoustic_pressure_contract COMMAND acoustic_pressure_contract)
    add_test(NAME acoustic_pressure_scaling_audit COMMAND acoustic_pressure_scaling_audit)
    add_test(NAME acoustic_pressure_phase_contract COMMAND acoustic_pressure_phase_contract)
    add_test(NAME acoustic_amplitude_contract COMMAND acoustic_amplitude_contract)
    add_test(NAME rotational_excitation_contract COMMAND rotational_excitation_contract)
    add_test(NAME excitation_model_audit COMMAND excitation_model_audit)
    add_test(NAME localized_strike_spectrum_audit COMMAND localized_strike_spectrum_audit)
    add_test(NAME localized_strike_direction_contract COMMAND localized_strike_direction_contract)
    add_test(NAME localized_strike_mesh_convergence COMMAND localized_strike_mesh_convergence)
    add_test(NAME inner_normal_offset_contract COMMAND inner_normal_offset_contract)
    add_test(NAME acoustic_multipole_audit COMMAND acoustic_multipole_audit)
    add_test(NAME acoustic_radiation_audit COMMAND acoustic_radiation_audit)
    add_test(NAME internal_fem_sanity COMMAND internal_fem_compare)
    add_test(NAME ccx_current_geometry_usage COMMAND ccx_current_geometry)
    add_test(NAME mesh_convergence COMMAND mesh_convergence "${CMAKE_BINARY_DIR}/mesh_convergence.csv")
    add_test(NAME thickness_mesh_convergence COMMAND thickness_mesh_convergence "${CMAKE_BINARY_DIR}/thickness_mesh_convergence.csv")
    add_test(NAME input_safety_contract COMMAND input_safety_contract)
    add_test(NAME tradition_assessment COMMAND tradition_assessment)
    add_test(NAME model_layer_contract COMMAND model_layer_contract)
    add_test(NAME modal_influence_contract COMMAND modal_influence_contract)
    add_test(NAME coordinate_morphing_contract COMMAND coordinate_morphing_contract)
    add_test(NAME mechanical_coordinate_morphing_contract COMMAND mechanical_coordinate_morphing_contract)
    add_test(NAME coordinate_morphing_ellipticity_sweep COMMAND coordinate_morphing_ellipticity_sweep)
    add_test(NAME nonaxisymmetric_reference_mesh_contract COMMAND nonaxisymmetric_reference_mesh_contract)
    add_test(NAME ellipticity_search_policy_contract COMMAND ellipticity_search_policy_contract)
    add_test(NAME cylindrical_strain_contract COMMAND cylindrical_strain_contract)
    add_test(NAME cylindrical_q4_element_contract COMMAND cylindrical_q4_element_contract)
    add_test(NAME cylindrical_global_assembly_contract COMMAND cylindrical_global_assembly_contract)
    add_test(NAME cylindrical_bell_geometry_path_contract COMMAND cylindrical_bell_geometry_path_contract)
    add_test(NAME cylindrical_global_modes_contract COMMAND cylindrical_global_modes_contract)
    add_test(NAME cylindrical_internal_solver_compare COMMAND cylindrical_internal_solver_compare)
    add_test(NAME cylindrical_ellipticity_convergence COMMAND cylindrical_ellipticity_convergence)
    add_test(NAME cylindrical_modal_family_tracking COMMAND cylindrical_modal_family_tracking)
    add_test(NAME modal_subspace_similarity_contract COMMAND modal_subspace_similarity_contract)
    add_test(NAME modal_subspace_conditioning_contract COMMAND modal_subspace_conditioning_contract)
    add_test(NAME cylindrical_bell_geometry_modal_tracking COMMAND cylindrical_bell_geometry_modal_tracking)
    add_test(NAME cylindrical_harmonic_family_continuation COMMAND cylindrical_harmonic_family_continuation)
    add_test(NAME cylindrical_unified_ellipticity_audit COMMAND cylindrical_unified_ellipticity_audit)
    set_tests_properties(ccx_current_geometry_usage PROPERTIES WILL_FAIL TRUE)

    add_executable(chinese_bell_strike_selectivity tests/chinese_bell_strike_selectivity.cpp)
    target_link_libraries(chinese_bell_strike_selectivity PRIVATE bell_core)
    add_test(NAME chinese_bell_strike_selectivity COMMAND chinese_bell_strike_selectivity)
    add_executable(chinese_bell_strike_surface_normal_contract tests/chinese_bell_strike_surface_normal_contract.cpp)
    target_link_libraries(chinese_bell_strike_surface_normal_contract PRIVATE bell_core)
    add_test(NAME chinese_bell_strike_surface_normal_contract COMMAND chinese_bell_strike_surface_normal_contract)
    add_executable(chinese_bell_strike_position_search tests/chinese_bell_strike_position_search.cpp)
    target_link_libraries(chinese_bell_strike_position_search PRIVATE bell_core)
    add_test(NAME chinese_bell_strike_position_search COMMAND chinese_bell_strike_position_search)
    add_executable(chinese_bell_strike_mesh_convergence tests/chinese_bell_strike_mesh_convergence.cpp)
    target_link_libraries(chinese_bell_strike_mesh_convergence PRIVATE bell_core)
    add_test(NAME chinese_bell_strike_mesh_convergence COMMAND chinese_bell_strike_mesh_convergence)

    add_executable(chinese_bell_strike_width_sensitivity tests/chinese_bell_strike_width_sensitivity.cpp)
    target_link_libraries(chinese_bell_strike_width_sensitivity PRIVATE bell_core)
    add_test(NAME chinese_bell_strike_width_sensitivity COMMAND chinese_bell_strike_width_sensitivity)
    add_executable(chinese_bell_objective_contract tests/chinese_bell_objective_contract.cpp)
    target_link_libraries(chinese_bell_objective_contract PRIVATE bell_experimental)
    add_test(NAME chinese_bell_objective_contract COMMAND chinese_bell_objective_contract)
    add_executable(modal_state_tracker_contract tests/modal_state_tracker_contract.cpp)
    target_link_libraries(modal_state_tracker_contract PRIVATE bell_experimental)
    add_test(NAME modal_state_tracker_contract COMMAND modal_state_tracker_contract)

    add_executable(modal_tracker_fem_continuation tests/modal_tracker_fem_continuation.cpp)
    target_link_libraries(modal_tracker_fem_continuation PRIVATE bell_experimental)
    add_test(NAME modal_tracker_fem_continuation COMMAND modal_tracker_fem_continuation)

    add_executable(cylindrical_modal_field_cross_mesh tests/cylindrical_modal_field_cross_mesh.cpp)
    target_link_libraries(cylindrical_modal_field_cross_mesh PRIVATE bell_experimental)
    add_test(NAME cylindrical_modal_field_cross_mesh COMMAND cylindrical_modal_field_cross_mesh)

    add_executable(cylindrical_modal_field_fem_matrix tests/cylindrical_modal_field_fem_matrix.cpp)
    target_link_libraries(cylindrical_modal_field_fem_matrix PRIVATE bell_experimental)
    add_test(NAME cylindrical_modal_field_fem_matrix COMMAND cylindrical_modal_field_fem_matrix)

    add_executable(chinese_bell_readiness_sweep tests/chinese_bell_readiness_sweep.cpp)
    target_link_libraries(chinese_bell_readiness_sweep PRIVATE bell_core)
    add_test(NAME chinese_bell_readiness_sweep COMMAND chinese_bell_readiness_sweep)

endif() # BUILD_SMOKE_TESTS

message(STATUS "")
message(STATUS "BellOptimizer ${PROJECT_VERSION}")
message(STATUS "  C++:             ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build type:      ${CMAKE_BUILD_TYPE}")
message(STATUS "  CLI:             bell (optimize, synth, synth-geometry)")
message(STATUS "  Smoke tests:     ${BUILD_SMOKE_TESTS}")
message(STATUS "  CCX executable:  ${CCX_EXECUTABLE}")
message(STATUS "  External C++ libs: none")
message(STATUS "")


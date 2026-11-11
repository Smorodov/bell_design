cmake_minimum_required(VERSION 3.16)
project(BellOptimizer VERSION 1.2.184 LANGUAGES CXX)

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
option(BUILD_CHINESE_RESEARCH "Build detached Chinese-bell research module" OFF)
set(CCX_EXECUTABLE "ccx" CACHE STRING "CalculiX executable used by the reference modal solver")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(BELL_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(BELL_CORE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/core)

add_library(bell_core STATIC
    src/core/bell_generator.cpp
    src/core/bell_material.cpp
    src/core/bell_optimizer.cpp
    src/core/modal_tracking.cpp
    src/core/multifidelity_pipeline.cpp
    src/core/physical_modal_selection.cpp
    research/round_bell_audit/modal_family_adapter.cpp
    research/round_bell_audit/round_bell_candidate_analysis.cpp
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
    src/excitation/modal_impulse_response.cpp
    src/excitation/modal_decay.cpp
    src/excitation/hammer_contact.cpp
    src/excitation/internal_clapper_geometry.cpp
    src/excitation/clapper_dynamics.cpp
    src/excitation/clapper_modal_evaluator.cpp
    src/excitation/clapper_modal_state_provider.cpp
    src/excitation/clapper_acoustic_objective.cpp
    src/excitation/meridional_influence.cpp
    src/acoustics/acoustic_directional_field.cpp
    src/acoustics/directional_radiation.cpp
    src/acoustics/acoustic_field_cache.cpp
    src/acoustics/bell_modal_profile.cpp
    src/audio/wav_writer.cpp
    src/audio/modal_pressure_synthesis.cpp
    src/acoustics/listener_geometry.cpp
    src/acoustics/finite_distance_radiation.cpp
    src/ccx/ccx_solver.cpp
    src/ccx/ccx_input_writer.cpp
)
target_include_directories(bell_core PUBLIC ${BELL_SRC} ${BELL_CORE_SRC})
target_compile_definitions(bell_core PRIVATE BELL_CCX_EXECUTABLE="${CCX_EXECUTABLE}")

# Standalone synthesis SDK. This target deliberately excludes optimizer, FEM,
# mesh generation, CalculiX integration and command-line application code.
add_library(bell_synthesis
    src/excitation/internal_clapper_geometry.cpp
    src/excitation/clapper_dynamics.cpp
    src/excitation/clapper_modal_evaluator.cpp
    src/acoustics/acoustic_directional_field.cpp
    src/acoustics/directional_radiation.cpp
    src/acoustics/acoustic_field_cache.cpp
    src/acoustics/bell_modal_profile.cpp
    src/acoustics/listener_geometry.cpp
    src/acoustics/finite_distance_radiation.cpp
    src/audio/modal_pressure_synthesis.cpp
    src/audio/wav_writer.cpp
)
add_library(BellSynthesis::bell_synthesis ALIAS bell_synthesis)
target_include_directories(bell_synthesis PUBLIC
    $<BUILD_INTERFACE:${BELL_SRC}>
    $<INSTALL_INTERFACE:include/bell_synthesis>
)
target_compile_features(bell_synthesis PUBLIC cxx_std_17)
set_target_properties(bell_synthesis PROPERTIES POSITION_INDEPENDENT_CODE ON WINDOWS_EXPORT_ALL_SYMBOLS ON)
target_compile_definitions(bell_synthesis PUBLIC BELL_SYNTHESIS_VERSION="${PROJECT_VERSION}")
if(MSVC)
    target_compile_options(bell_synthesis PRIVATE /W4)
else()
    target_compile_options(bell_synthesis PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
install(TARGETS bell_synthesis EXPORT BellSynthesisTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(FILES
    src/core/bell_types.h
    src/core/modal_types.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/bell_synthesis/core
)
install(DIRECTORY src/acoustics/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/bell_synthesis/acoustics
    FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY src/excitation/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/bell_synthesis/excitation
    FILES_MATCHING PATTERN "clapper_dynamics.hpp" PATTERN "clapper_modal_evaluator.hpp" PATTERN "internal_clapper_geometry.hpp")
install(DIRECTORY src/audio/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/bell_synthesis/audio
    FILES_MATCHING PATTERN "*.hpp")
install(EXPORT BellSynthesisTargets NAMESPACE BellSynthesis:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BellSynthesis)
configure_package_config_file(cmake/BellSynthesisConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/BellSynthesisConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BellSynthesis)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/BellSynthesisConfigVersion.cmake
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/BellSynthesisConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/BellSynthesisConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BellSynthesis)

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
    src/commands/synthesize_contact.cpp
    src/synthesis/bell_synthesis.cpp
)
target_include_directories(bell PRIVATE ${BELL_SRC} ${BELL_CORE_SRC} ${BELL_SRC}/synthesis)
target_link_libraries(bell PRIVATE bell_core)
target_compile_definitions(bell PRIVATE BELL_VERSION="${PROJECT_VERSION}")
if(MSVC)
    target_compile_options(bell PRIVATE /W4)
else()
    target_compile_options(bell PRIVATE -Wall -Wextra -Wpedantic)
endif()


if(BUILD_SMOKE_TESTS)
    enable_testing()

    # This contract links ONLY the distributable synthesis SDK target.
    add_executable(synthesis_sdk_link_contract tests/synthesis_sdk_link_contract.cpp)
    target_link_libraries(synthesis_sdk_link_contract PRIVATE BellSynthesis::bell_synthesis)
    add_test(NAME synthesis_sdk_link_contract COMMAND synthesis_sdk_link_contract)

    # Experimental diagnostics are deliberately kept out of the production library.
    add_library(bell_experimental STATIC
        src/core/modal_state_tracker.cpp
        src/fem/cylindrical_modal_field.cpp
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

    add_executable(wav_writer_contract tests/wav_writer_contract.cpp)
    target_link_libraries(wav_writer_contract PRIVATE bell_core)

    add_executable(directional_radiation_equivalence tests/directional_radiation_equivalence.cpp)
    target_link_libraries(directional_radiation_equivalence PRIVATE bell_core)

    add_executable(cylindrical_axis_regularity_contract tests/cylindrical_axis_regularity_contract.cpp)
    target_link_libraries(cylindrical_axis_regularity_contract PRIVATE bell_core)

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

    add_executable(multifidelity_calibration_probe tests/multifidelity_calibration_probe.cpp)
    target_link_libraries(multifidelity_calibration_probe PRIVATE bell_core)

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
    add_executable(acoustic_directional_field_rotation_contract tests/acoustic_directional_field_rotation_contract.cpp)
    add_executable(stereo_listener_geometry_contract tests/stereo_listener_geometry_contract.cpp)
    add_executable(stereo_coordinate_covariance_contract tests/stereo_coordinate_covariance_contract.cpp)
    add_executable(clapper_stereo_frame_covariance_contract tests/clapper_stereo_frame_covariance_contract.cpp)
    add_executable(stereo_polar_angle_approximation_contract tests/stereo_polar_angle_approximation_contract.cpp)
    add_executable(polar_directional_field_interpolation_contract tests/polar_directional_field_interpolation_contract.cpp)
    add_executable(finite_distance_radiation_contract tests/finite_distance_radiation_contract.cpp)
    add_executable(finite_distance_real_fem_contract tests/finite_distance_real_fem_contract.cpp)
    add_executable(real_fem_polar_interpolation_contract tests/real_fem_polar_interpolation_contract.cpp)
    add_executable(acoustic_field_cache_contract tests/acoustic_field_cache_contract.cpp)
    add_executable(bell_modal_profile_contract tests/bell_modal_profile_contract.cpp)
    target_link_libraries(bell_modal_profile_contract PRIVATE bell_core)
    add_executable(stereo_phase_delay_contract tests/stereo_phase_delay_contract.cpp)
    target_link_libraries(stereo_phase_delay_contract PRIVATE bell_core)
    target_link_libraries(acoustic_field_cache_contract PRIVATE bell_core)
    target_link_libraries(finite_distance_real_fem_contract PRIVATE bell_core)
    target_link_libraries(real_fem_polar_interpolation_contract PRIVATE bell_core)
    target_link_libraries(finite_distance_radiation_contract PRIVATE bell_core)
    target_link_libraries(stereo_listener_geometry_contract PRIVATE bell_core)
    target_include_directories(stereo_coordinate_covariance_contract PRIVATE ${BELL_SRC} ${BELL_CORE_SRC})
    target_link_libraries(stereo_coordinate_covariance_contract PRIVATE bell_core)
    target_link_libraries(clapper_stereo_frame_covariance_contract PRIVATE bell_core)
    target_link_libraries(stereo_polar_angle_approximation_contract PRIVATE bell_core)
    target_link_libraries(polar_directional_field_interpolation_contract PRIVATE bell_core)
    target_include_directories(acoustic_directional_field_rotation_contract PRIVATE ${BELL_SRC} ${BELL_CORE_SRC})
    target_link_libraries(acoustic_directional_field_rotation_contract PRIVATE bell_core)
    add_executable(excitation_model_audit tests/excitation_model_audit.cpp)
    add_executable(localized_strike_spectrum_audit tests/localized_strike_spectrum_audit.cpp)
    target_link_libraries(localized_strike_spectrum_audit PRIVATE bell_core)
    add_executable(localized_strike_direction_contract tests/localized_strike_direction_contract.cpp)
    target_link_libraries(localized_strike_direction_contract PRIVATE bell_core)
    add_executable(circular_localized_strike_projection tests/circular_localized_strike_projection.cpp)
    target_link_libraries(circular_localized_strike_projection PRIVATE bell_core)

    add_executable(modal_impulse_response_contract tests/modal_impulse_response_contract.cpp)
    target_link_libraries(modal_impulse_response_contract PRIVATE bell_core)
    add_test(NAME modal_impulse_response_contract COMMAND modal_impulse_response_contract)

    add_executable(modal_decay_contract tests/modal_decay_contract.cpp)
    target_link_libraries(modal_decay_contract PRIVATE bell_core)
    add_test(NAME modal_decay_contract COMMAND modal_decay_contract)

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
    add_executable(material_preset_contract tests/material_preset_contract.cpp)
    target_link_libraries(material_preset_contract PRIVATE bell_core)
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

    add_executable(target_profile_loader_contract tests/target_profile_loader_contract.cpp)
    target_include_directories(target_profile_loader_contract PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/core)
    target_compile_definitions(target_profile_loader_contract PRIVATE BELL_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(NAME target_profile_loader_contract COMMAND target_profile_loader_contract)

    add_executable(multifidelity_pipeline_contract tests/multifidelity_pipeline_contract.cpp)
    target_link_libraries(multifidelity_pipeline_contract PRIVATE bell_core)
    add_test(NAME multifidelity_pipeline_contract COMMAND multifidelity_pipeline_contract)

    add_executable(physical_modal_selection_contract tests/physical_modal_selection_contract.cpp)
    target_link_libraries(physical_modal_selection_contract PRIVATE bell_core)
    add_test(NAME physical_modal_selection_contract COMMAND physical_modal_selection_contract)

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
    add_test(NAME acoustic_directional_field_rotation_contract COMMAND acoustic_directional_field_rotation_contract)
    add_test(NAME directional_radiation_equivalence COMMAND directional_radiation_equivalence)
    add_test(NAME wav_writer_contract COMMAND wav_writer_contract)
    add_test(NAME stereo_listener_geometry_contract COMMAND stereo_listener_geometry_contract)
    add_test(NAME stereo_coordinate_covariance_contract COMMAND stereo_coordinate_covariance_contract)
    add_test(NAME clapper_stereo_frame_covariance_contract COMMAND clapper_stereo_frame_covariance_contract)
    add_test(NAME stereo_polar_angle_approximation_contract COMMAND stereo_polar_angle_approximation_contract)
    add_test(NAME polar_directional_field_interpolation_contract COMMAND polar_directional_field_interpolation_contract)
    add_test(NAME finite_distance_radiation_contract COMMAND finite_distance_radiation_contract)
    add_test(NAME finite_distance_real_fem_contract COMMAND finite_distance_real_fem_contract)
    add_test(NAME real_fem_polar_interpolation_contract COMMAND real_fem_polar_interpolation_contract)
    add_test(NAME acoustic_field_cache_contract COMMAND acoustic_field_cache_contract)
    add_test(NAME bell_modal_profile_contract COMMAND bell_modal_profile_contract)
    add_test(NAME stereo_phase_delay_contract COMMAND stereo_phase_delay_contract)
    add_test(NAME excitation_model_audit COMMAND excitation_model_audit)
    add_test(NAME localized_strike_spectrum_audit COMMAND localized_strike_spectrum_audit)
    add_test(NAME localized_strike_direction_contract COMMAND localized_strike_direction_contract)
    add_test(NAME circular_localized_strike_projection COMMAND circular_localized_strike_projection)
    add_test(NAME cylindrical_axis_regularity_contract COMMAND cylindrical_axis_regularity_contract)
    add_test(NAME multifidelity_calibration_probe COMMAND multifidelity_calibration_probe)
    set_tests_properties(multifidelity_calibration_probe PROPERTIES TIMEOUT 300)
    add_test(NAME direct_elliptical_fem_evaluation COMMAND direct_elliptical_fem_evaluation)
    add_test(NAME localized_strike_mesh_convergence COMMAND localized_strike_mesh_convergence)
    add_test(NAME inner_normal_offset_contract COMMAND inner_normal_offset_contract)
    add_test(NAME acoustic_multipole_audit COMMAND acoustic_multipole_audit)
    add_test(NAME acoustic_radiation_audit COMMAND acoustic_radiation_audit)
    add_test(NAME internal_fem_sanity COMMAND internal_fem_compare)
    add_test(NAME ccx_current_geometry_usage COMMAND ccx_current_geometry)
    add_test(NAME mesh_convergence COMMAND mesh_convergence "${CMAKE_BINARY_DIR}/mesh_convergence.csv")
    add_test(NAME thickness_mesh_convergence COMMAND thickness_mesh_convergence "${CMAKE_BINARY_DIR}/thickness_mesh_convergence.csv")
    add_test(NAME input_safety_contract COMMAND input_safety_contract)
    add_test(NAME material_preset_contract COMMAND material_preset_contract)
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
    add_executable(hammer_contact_contract tests/hammer_contact_contract.cpp)
    target_link_libraries(hammer_contact_contract PRIVATE bell_core)
    add_test(NAME hammer_contact_contract COMMAND hammer_contact_contract)

    add_executable(internal_clapper_geometry_contract tests/internal_clapper_geometry_contract.cpp)
    target_link_libraries(internal_clapper_geometry_contract PRIVATE bell_core)
    add_test(NAME internal_clapper_geometry_contract COMMAND internal_clapper_geometry_contract)

    add_executable(inner_surface_modal_export_contract tests/inner_surface_modal_export_contract.cpp)
    target_link_libraries(inner_surface_modal_export_contract PRIVATE bell_core)
    add_test(NAME inner_surface_modal_export_contract COMMAND inner_surface_modal_export_contract)

    add_executable(clapper_dynamics_contract tests/clapper_dynamics_contract.cpp)
    target_link_libraries(clapper_dynamics_contract PRIVATE bell_core)
    add_test(NAME clapper_dynamics_contract COMMAND clapper_dynamics_contract)

    add_executable(clapper_modal_evaluator_contract tests/clapper_modal_evaluator_contract.cpp)
    target_link_libraries(clapper_modal_evaluator_contract PRIVATE bell_core)
    add_test(NAME clapper_modal_evaluator_contract COMMAND clapper_modal_evaluator_contract)

    add_executable(clapper_candidate_data_reuse_contract tests/clapper_candidate_data_reuse_contract.cpp)
    target_link_libraries(clapper_candidate_data_reuse_contract PRIVATE bell_core)
    add_test(NAME clapper_candidate_data_reuse_contract COMMAND clapper_candidate_data_reuse_contract)
    add_executable(clapper_acoustic_objective_contract tests/clapper_acoustic_objective_contract.cpp)
    target_link_libraries(clapper_acoustic_objective_contract PRIVATE bell_core)
    add_executable(clapper_scaling_contract tests/clapper_scaling_contract.cpp)
    target_link_libraries(clapper_scaling_contract PRIVATE bell_core)
    add_test(NAME clapper_acoustic_objective_contract COMMAND clapper_acoustic_objective_contract)
    add_test(NAME clapper_scaling_contract COMMAND clapper_scaling_contract)

    add_executable(clapper_scale_energy_audit tests/clapper_scale_energy_audit.cpp)
    target_link_libraries(clapper_scale_energy_audit PRIVATE bell_core)
    add_test(NAME clapper_scale_energy_audit COMMAND clapper_scale_energy_audit)

    add_executable(clapper_modal_spectrum_family_audit tests/clapper_modal_spectrum_family_audit.cpp)
    target_link_libraries(clapper_modal_spectrum_family_audit PRIVATE bell_core)
    add_test(NAME clapper_modal_spectrum_family_audit COMMAND clapper_modal_spectrum_family_audit)
    set_tests_properties(clapper_modal_spectrum_family_audit PROPERTIES TIMEOUT 180)
    add_test(NAME contact_directional_audio_smoke COMMAND $<TARGET_FILE:bell> synth-contact
        "${CMAKE_CURRENT_SOURCE_DIR}/bells/best_geometry.json"
        --output "${CMAKE_BINARY_DIR}/contact_directional_smoke.wav"
        --mesh-nodes 15 --duration 0.25 --contact-duration 0.02
        --strike-azimuth-deg 37 --listener-polar-deg 25 --listener-azimuth-deg 80)
    set_tests_properties(contact_directional_audio_smoke PROPERTIES TIMEOUT 120)
    add_test(NAME internal_clapper_audio_smoke COMMAND $<TARGET_FILE:bell> synth-clapper
        "${CMAKE_CURRENT_SOURCE_DIR}/bells/best_geometry.json"
        --output "${CMAKE_BINARY_DIR}/internal_clapper_smoke.wav"
        --mesh-nodes 15 --duration 0.25 --impact-horizon 0.25 --contact-step 0.000002)
    set_tests_properties(internal_clapper_audio_smoke PROPERTIES TIMEOUT 180)

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

    

endif() # BUILD_SMOKE_TESTS

if(BUILD_CHINESE_RESEARCH)
    add_subdirectory(research/chinese_bell)
endif()

message(STATUS "")
message(STATUS "BellOptimizer ${PROJECT_VERSION}")
message(STATUS "  C++:             ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build type:      ${CMAKE_BUILD_TYPE}")
message(STATUS "  CLI:             bell (optimize, synth, synth-geometry, synth-contact, synth-clapper)")
message(STATUS "  Smoke tests:     ${BUILD_SMOKE_TESTS}")
message(STATUS "  Chinese research: ${BUILD_CHINESE_RESEARCH}")
message(STATUS "  CCX executable:  ${CCX_EXECUTABLE}")
message(STATUS "  External C++ libs: none")
message(STATUS "")






if(BUILD_SMOKE_TESTS)
    add_executable(round_geometry_dataset research/round_bell_audit/round_geometry_dataset.cpp)
    target_link_libraries(round_geometry_dataset PRIVATE bell_core)
endif()

if(BUILD_SMOKE_TESTS)
    add_executable(round_geometry_modal_convergence research/round_bell_audit/round_geometry_modal_convergence.cpp)
    target_link_libraries(round_geometry_modal_convergence PRIVATE bell_core)
endif()

if(BUILD_SMOKE_TESTS)
    add_executable(round_geometry_multifidelity_probe research/round_bell_audit/round_geometry_multifidelity_probe.cpp)
    target_link_libraries(round_geometry_multifidelity_probe PRIVATE bell_core)
endif()

# Differential diagnostic: one-to-one full-field matching with the full consistent FEM mass matrix.
# This remains research-only and does not alter production optimizer behavior.
if(BUILD_SMOKE_TESTS)
    add_executable(round_solver_full_volume_probe research/round_bell_audit/round_solver_full_volume_probe.cpp)
    target_link_libraries(round_solver_full_volume_probe PRIVATE bell_core)
endif()

# Full consistent-mass rank-2 subspace audit on a circular bell.
if(BUILD_SMOKE_TESTS)
    add_executable(round_bell_consistent_mass_subspace_probe research/round_bell_audit/round_bell_consistent_mass_subspace_probe.cpp)
    target_link_libraries(round_bell_consistent_mass_subspace_probe PRIVATE bell_core)
    add_test(NAME round_bell_consistent_mass_subspace_probe COMMAND round_bell_consistent_mass_subspace_probe)
endif()

if(BUILD_SMOKE_TESTS)
    add_executable(round_geometry_fullfield_match research/round_bell_audit/round_geometry_fullfield_match.cpp)
    target_link_libraries(round_geometry_fullfield_match PRIVATE bell_core bell_experimental)
endif()

if(BUILD_SMOKE_TESTS)
    add_executable(round_geometry_multiharmonic_match research/round_bell_audit/round_geometry_multiharmonic_match.cpp)
    target_link_libraries(round_geometry_multiharmonic_match PRIVATE bell_core bell_experimental)
endif()

if(BUILD_SMOKE_TESTS)
    add_executable(round_geometry_thickness_sweep research/round_bell_audit/round_geometry_thickness_sweep.cpp)
    target_link_libraries(round_geometry_thickness_sweep PRIVATE bell_core)
endif()

if(BUILD_SMOKE_TESTS)
    add_executable(round_geometry_thickness_mass_sweep research/round_bell_audit/round_geometry_thickness_mass_sweep.cpp)
    target_link_libraries(round_geometry_thickness_mass_sweep PRIVATE bell_core bell_experimental)
endif()

# Round-bell optimizer integration diagnostic. This is intentionally a research
# executable: it validates the direct cylindrical FEM + localized circular
# strike hand-off without changing the production optimizer/objective.
if(BUILD_SMOKE_TESTS)
    add_executable(round_optimizer_integration_probe
        research/round_bell_audit/round_optimizer_integration_probe.cpp
        research/round_bell_audit/modal_family_adapter.cpp
        research/round_bell_audit/round_bell_candidate_analysis.cpp)
    target_link_libraries(round_optimizer_integration_probe PRIVATE bell_core)

    add_executable(modal_family_adapter_contract
        tests/modal_family_adapter_contract.cpp
        research/round_bell_audit/modal_family_adapter.cpp)
    target_link_libraries(modal_family_adapter_contract PRIVATE bell_core)
    add_test(NAME modal_family_adapter_contract COMMAND modal_family_adapter_contract)
endif()

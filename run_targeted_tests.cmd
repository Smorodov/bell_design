@echo off
setlocal
chcp 65001 >nul

if not exist build (
    echo ERROR: build directory not found. Run build.cmd first.
    exit /b 1
)

pushd build

rem Focused regression set for the current cylindrical modal-continuation audit.
rem Keep output limited to the affected path and its nearest contracts.
ctest -R "^(rcm_initial_vector_contract|modal_subspace_similarity_contract|modal_subspace_conditioning_contract|cylindrical_ellipticity_convergence|cylindrical_modal_family_tracking|cylindrical_bell_geometry_modal_tracking|cylindrical_harmonic_family_continuation|cylindrical_unified_ellipticity_audit|modal_tracker_fem_continuation|cylindrical_modal_field_cross_mesh|cylindrical_modal_field_fem_matrix)$" --output-on-failure
set ERR=%ERRORLEVEL%

popd
endlocal & exit /b %ERR%

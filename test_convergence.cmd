@echo off
setlocal EnableExtensions
chcp 65001 >nul
set "ROOT=%~dp0"
set "MODE=%~1"
if not defined MODE set "MODE=all"
if /I not "%MODE%"=="all" if /I not "%MODE%"=="mesh" if /I not "%MODE%"=="thickness" (
  echo Usage: test_convergence.cmd [all^|mesh^|thickness]
  exit /b 2
)
if not defined CMAKE_GENERATOR set "CMAKE_GENERATOR=MinGW Makefiles"
if not defined BUILD_JOBS set "BUILD_JOBS=4"

cmake -S "%ROOT%" -B "%ROOT%build" -G "%CMAKE_GENERATOR%" -DBUILD_SMOKE_TESTS=ON
if errorlevel 1 exit /b 1

if /I "%MODE%"=="all" (
  cmake --build "%ROOT%build" --target mesh_convergence thickness_mesh_convergence -- -j%BUILD_JOBS%
  if errorlevel 1 exit /b 1
  "%ROOT%build\mesh_convergence.exe" "%ROOT%build\mesh_convergence.csv"
  if errorlevel 1 exit /b 1
  "%ROOT%build\thickness_mesh_convergence.exe" "%ROOT%build\thickness_mesh_convergence.csv"
  if errorlevel 1 exit /b 1
) else if /I "%MODE%"=="mesh" (
  cmake --build "%ROOT%build" --target mesh_convergence -- -j%BUILD_JOBS%
  if errorlevel 1 exit /b 1
  "%ROOT%build\mesh_convergence.exe" "%ROOT%build\mesh_convergence.csv"
  if errorlevel 1 exit /b 1
) else (
  cmake --build "%ROOT%build" --target thickness_mesh_convergence -- -j%BUILD_JOBS%
  if errorlevel 1 exit /b 1
  "%ROOT%build\thickness_mesh_convergence.exe" "%ROOT%build\thickness_mesh_convergence.csv"
  if errorlevel 1 exit /b 1
)

echo.
echo Convergence reports are in build\
endlocal

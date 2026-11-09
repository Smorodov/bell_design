@echo off
setlocal
rem Run this script from the MSYS2 UCRT64 shell using: cmd.exe /c test_windows_ucrt64.cmd
rem Alternatively run the equivalent commands in UCRT64 bash shown in docs/WINDOWS_UCRT64_AUDIT.md.
set ROOT=%~dp0
cd /d "%ROOT%"
if not defined MSYSTEM (
  echo ERROR: MSYSTEM is not set. Start the MSYS2 UCRT64 environment first.
  exit /b 2
)
if not "%MSYSTEM%"=="UCRT64" (
  echo ERROR: Expected MSYSTEM=UCRT64, got %MSYSTEM%.
  exit /b 2
)
if not defined CXX set CXX=g++
if not defined CC set CC=gcc
if not defined CMAKE_GENERATOR set CMAKE_GENERATOR=MinGW Makefiles
if not exist audit_results mkdir audit_results
(echo DATE=%DATE% %TIME% & echo MSYSTEM=%MSYSTEM% & echo MINGW_PREFIX=%MINGW_PREFIX% & echo CXX=%CXX% & echo CC=%CC% & cmake --version & %CXX% --version & ninja --version 2>nul) > audit_results/windows_ucrt64_environment.txt 2>&1
cmake -S . -B build_ucrt64 -G "%CMAKE_GENERATOR%" -DCMAKE_BUILD_TYPE=Release -DBUILD_SMOKE_TESTS=ON -DBUILD_CHINESE_RESEARCH=OFF > audit_results/windows_ucrt64_configure.log 2>&1
if errorlevel 1 exit /b %errorlevel%
cmake --build build_ucrt64 --parallel 4 > audit_results/windows_ucrt64_build.log 2>&1
if errorlevel 1 exit /b %errorlevel%
ctest --test-dir build_ucrt64 --output-on-failure -j2 > audit_results/windows_ucrt64_ctest.log 2>&1
set TEST_RC=%errorlevel%
certutil -hashfile CMakeLists.txt SHA256 > audit_results/windows_ucrt64_sha256.txt
certutil -hashfile src\core\modal_tracking.cpp SHA256 >> audit_results/windows_ucrt64_sha256.txt
certutil -hashfile src\core\bell_optimizer.cpp SHA256 >> audit_results/windows_ucrt64_sha256.txt
exit /b %TEST_RC%

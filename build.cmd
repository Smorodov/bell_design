@echo off
setlocal

if not exist build mkdir build

cmake -S . -B build -DBUILD_SMOKE_TESTS=ON -DCMAKE_BUILD_TYPE=Release
if errorlevel 1 exit /b %errorlevel%

cmake --build build --parallel
if errorlevel 1 exit /b %errorlevel%

endlocal

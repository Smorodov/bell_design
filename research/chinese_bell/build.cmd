@echo off
setlocal

set "ROOT=%~dp0..\..\.."
set "BUILD=%ROOT%\build-chinese"

cmake -S "%ROOT%" -B "%BUILD%" -DCMAKE_BUILD_TYPE=Release -DBUILD_SMOKE_TESTS=ON -DBUILD_CHINESE_RESEARCH=ON
if errorlevel 1 exit /b 1

cmake --build "%BUILD%" --target bell_chinese_research --parallel 2
if errorlevel 1 exit /b 1

ctest --test-dir "%BUILD%" -R "chinese_" --output-on-failure
exit /b %errorlevel%

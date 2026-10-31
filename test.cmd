@echo off
setlocal
ctest --test-dir "%~dp0build" --output-on-failure
endlocal

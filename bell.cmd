@echo off
setlocal
"%~dp0build\bell.exe" %*
set RC=%ERRORLEVEL%
endlocal & exit /b %RC%

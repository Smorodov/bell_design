@echo off
setlocal
for %%C in (run_western_220.cmd run_western_440.cmd run_western_880.cmd run_orthodox_220.cmd run_orthodox_440.cmd run_orthodox_880.cmd) do (
  echo.
  echo ===== %%C =====
  call "%~dp0%%C"
  if errorlevel 1 goto failed
)
echo.
echo All six demonstration optimizations completed.
endlocal & exit /b 0

:failed
set "RC=%ERRORLEVEL%"
endlocal & exit /b %RC%

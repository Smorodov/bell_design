@echo off
setlocal EnableExtensions
chcp 65001 >nul
set "ROOT=%~dp0"
if not defined CCX_EXE set "CCX_EXE=%ROOT%tools\ccx\ccx.exe"
set "TEST_EXE=%ROOT%build\ccx_current_geometry.exe"
set "WORK=%ROOT%build\ccx_current_geometry_run"

if not exist "%TEST_EXE%" (
  echo [1/4] Building BellOptimizer test module...
  call "%ROOT%build.cmd"
  if errorlevel 1 exit /b 1
)
if not exist "%TEST_EXE%" (
  echo ERROR: test executable was not built: "%TEST_EXE%"
  exit /b 2
)
if not exist "%CCX_EXE%" (
  echo ERROR: CalculiX executable not found: "%CCX_EXE%"
  echo Copy the Windows CalculiX distribution to "%ROOT%tools\ccx\"
  echo The executable must be named ccx.exe. Keep its required DLLs beside it.
  echo Or set CCX_EXE to the full path of ccx.exe before running this script.
  exit /b 3
)

if not exist "%WORK%" mkdir "%WORK%"
if errorlevel 1 exit /b 4

echo [2/4] Preparing the corrected geometry and internal FEM reference...
"%TEST_EXE%" prepare "%WORK%"
if errorlevel 1 exit /b 5
if not exist "%WORK%\bell_current_geometry.inp" (
  echo ERROR: CalculiX input file was not created.
  exit /b 6
)

echo [3/4] Running CalculiX. This can take a while...
pushd "%WORK%"
if exist bell_current_geometry.dat del /q bell_current_geometry.dat
if exist bell_current_geometry.frd del /q bell_current_geometry.frd
if exist bell_current_geometry.sta del /q bell_current_geometry.sta
if exist bell_current_geometry.msg del /q bell_current_geometry.msg
"%CCX_EXE%" bell_current_geometry > ccx_console.log 2>&1
set "CCX_RC=%ERRORLEVEL%"
popd
if not "%CCX_RC%"=="0" (
  echo ERROR: CalculiX returned code %CCX_RC%.
  echo Inspect "%WORK%\ccx_console.log" and the .sta/.msg files.
  exit /b 7
)
if not exist "%WORK%\bell_current_geometry.dat" (
  echo ERROR: CalculiX did not create the expected .dat file.
  echo Inspect "%WORK%\ccx_console.log" and the .sta/.msg files.
  exit /b 8
)

echo [4/4] Comparing modal frequencies, effective masses, and mode shapes...
"%TEST_EXE%" compare "%WORK%"
if errorlevel 1 (
  echo ERROR: comparison failed. Inspect .dat and ccx_console.log in:
  echo "%WORK%"
  exit /b 9
)
echo.
echo Finished. Please send these files back for review:
echo   "%WORK%\comparison.csv"
echo   "%WORK%\mode_shape_comparison.csv"
echo   "%WORK%\internal_mode_shapes.csv"
echo   "%WORK%\comparison_summary.txt"
echo   "%WORK%\mesh_info.txt"
echo   "%WORK%\ccx_console.log"
echo   "%WORK%\bell_current_geometry.dat"
echo   "%WORK%\bell_current_geometry.frd"
endlocal

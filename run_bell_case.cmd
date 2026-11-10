@echo off
setlocal EnableExtensions
chcp 65001 >nul

rem Usage: run_bell_case.cmd western|orthodox 220|440|880
set "ROOT=%~dp0"
set "PROFILE=%~1"
set "FREQ=%~2"
if not defined PROFILE goto usage
if not defined FREQ goto usage

if /I "%PROFILE%"=="western" (
  set "EVOLUTION=%ROOT%bells\evolution_european.json"
) else if /I "%PROFILE%"=="orthodox" (
  set "EVOLUTION=%ROOT%bells\evolution_orthodox.json"
) else goto usage

if not "%FREQ%"=="220" if not "%FREQ%"=="440" if not "%FREQ%"=="880" goto usage
set "TARGETS=%ROOT%bells\targets_%PROFILE%_%FREQ%.json"
if not exist "%TARGETS%" (
  echo ERROR: Target preset not found: "%TARGETS%"
  exit /b 2
)
if not exist "%ROOT%build\bell.exe" (
  echo ERROR: build\bell.exe not found. Run build.cmd first.
  exit /b 2
)

set "RUNDIR=%ROOT%runs\%PROFILE%_%FREQ%Hz"
if not exist "%RUNDIR%" mkdir "%RUNDIR%"
copy /Y "%EVOLUTION%" "%RUNDIR%\evolution.json" >nul
if errorlevel 1 exit /b 1
copy /Y "%TARGETS%" "%RUNDIR%\targets.json" >nul
if errorlevel 1 exit /b 1

pushd "%RUNDIR%"
echo.
echo ===== OPTIMIZATION: %PROFILE% / %FREQ% Hz =====
"%ROOT%build\bell.exe" optimize --solver=internal --evolution evolution.json --targets targets.json
set "RC=%ERRORLEVEL%"
if not "%RC%"=="0" goto failed

if not exist best_geometry.json (
  echo ERROR: optimizer did not produce best_geometry.json
  set "RC=3"
  goto failed
)
if not exist audio mkdir audio

echo.
echo ===== MODAL SYNTHESIS PREVIEW =====
"%ROOT%build\bell.exe" synth-geometry best_geometry.json --output "audio\%PROFILE%_%FREQ%Hz_modal.wav"
set "RC=%ERRORLEVEL%"
if not "%RC%"=="0" goto failed
if not exist "audio\%PROFILE%_%FREQ%Hz_modal.wav" (
  echo ERROR: modal WAV was not created.
  set "RC=4"
  goto failed
)

echo.
echo ===== LOCALIZED STRIKE SYNTHESIS =====
"%ROOT%build\bell.exe" synth-contact best_geometry.json --output "audio\%PROFILE%_%FREQ%Hz_contact.wav" --mesh-nodes 31 --strike-height 0.72 --strike-azimuth-deg 0 --listener-polar-deg 30 --listener-azimuth-deg 0 --striker-mass 0.045 --striker-velocity 1.2 --hertz-stiffness 2e7 --contact-step 2e-6 --contact-duration 0.025 --t60 5 --duration 6 --sample-rate 48000
set "RC=%ERRORLEVEL%"
if not "%RC%"=="0" goto failed
if not exist "audio\%PROFILE%_%FREQ%Hz_contact.wav" (
  echo ERROR: contact WAV was not created.
  set "RC=4"
  goto failed
)

echo.
echo ===== COMPLETE =====
echo Geometry: "%RUNDIR%\best_geometry.json"
echo Audio files:
dir /b audio\*.wav
popd
endlocal & exit /b 0

:failed
echo.
echo Run failed with exit code %RC%.
popd
endlocal & exit /b %RC%

:usage
echo Usage: run_bell_case.cmd western^|orthodox 220^|440^|880
echo Example: run_bell_case.cmd orthodox 440
exit /b 2

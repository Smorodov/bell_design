@echo off
setlocal
chcp 65001 >nul
pushd "%~dp0"

if not exist build\bell.exe (
    echo ERROR: build\bell.exe not found. Run build.cmd first.
    popd
    exit /b 1
)
if not exist audio mkdir audio

build\bell.exe synth-contact bells\best_geometry.json --output audio\bell_contact.wav --mesh-nodes 71 --strike-height 0.72 --strike-azimuth-deg 0 --listener-polar-deg 30 --listener-azimuth-deg 0 --striker-mass 0.045 --striker-velocity 1.2 --hertz-stiffness 2e7 --contact-step 2e-6 --contact-duration 0.025 --t60 5 --duration 6 --sample-rate 48000
set ERR=%ERRORLEVEL%
if not "%ERR%"=="0" (
    popd
    endlocal & exit /b %ERR%
)

echo.
echo Generated: audio\bell_contact.wav
echo Metadata:  audio\bell_contact_metadata.txt
popd
endlocal

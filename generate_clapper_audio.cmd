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

build\bell.exe synth-clapper bells\best_geometry.json --output audio\bell_clapper.wav --mesh-nodes 31 --clapper-mass 5.0 --com-length 0.46 --tip-length 0.52 --tip-radius 0.035 --initial-angle-deg 0 --initial-angular-velocity 3.2 --strike-azimuth-deg 0 --listener-polar-deg 30 --listener-azimuth-deg 0 --contact-step 2e-6 --impact-horizon 0.25 --t60 5 --duration 6 --sample-rate 48000
set ERR=%ERRORLEVEL%
if not "%ERR%"=="0" (
    popd
    endlocal & exit /b %ERR%
)

echo.
echo Generated: audio\bell_clapper.wav (31-node listening preview)
echo Metadata:  audio\bell_clapper_metadata.txt
echo Modal data: audio\bell_clapper_modes.csv
popd
endlocal

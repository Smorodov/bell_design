@echo off
setlocal
chcp 65001 >nul

if not exist build\bell.exe (
    echo === Build ===
    call build.cmd
    if errorlevel 1 exit /b 1
)

if not exist audio mkdir audio

pushd bells
..\build\bell.exe synth best_bell.dat
if errorlevel 1 (
    popd
    exit /b 1
)

copy /Y bell_*.wav ..\audio\ >nul
copy /Y spectrum_*.csv ..\audio\ >nul
copy /Y envelopes_all.csv ..\audio\ >nul
popd

echo.
echo Audio files are in audio\
dir /b audio\*.wav
endlocal

@echo off
setlocal
if exist build rmdir /s /q build
if exist bells\bell.inp del /q bells\bell.inp
if exist bells\bell.sta del /q bells\bell.sta
if exist bells\bell.frd del /q bells\bell.frd
if exist bells\bell.12d del /q bells\bell.12d
if exist bells\bell.cvg del /q bells\bell.cvg
if exist bells\bell.eig del /q bells\bell.eig
if exist bells\spectrum_*.csv del /q bells\spectrum_*.csv
if exist bells\envelopes_all.csv del /q bells\envelopes_all.csv
endlocal

@ECHO OFF
REM create input files for GP calculation
copy %1.dat calculate.bat
REM erase %1.dat
copy %1-1.dat %1.CPS_001
REM erase %1-1.dat
copy %1-2.dat %1.GOM
REM erase %1-2.dat
copy %1-3.dat %1.OPD
REM erase %1-3.dat
REM create GP folder and copy input files, executable and dll's into it
md ..\%1.GP
copy calculate.bat "..\%1.GP\calculate.bat"
copy %1.CPS_001 "..\%1.GP\%1.CPS_001"
copy %1.GOM "..\%1.GP\%1.GOM"
REM copy %1.OPD "..\%1.GP\%1.OPD"
copy %3\exec\ADSIM_2025.exe "..\%1.GP\ADSIM_2025.exe"
xcopy %3\dll\*.dll "..\%1.GP\*.dll" /sy
@ECHO OFF
MKDIR "backup"

del "*.~*"
del "*.obj"                                                         

@ECHO OFF
setlocal enabledelayedexpansion
set first=1
set one = "_"
for /f "tokens=5" %%a in (constants.h) do (
      if !first! == 1 (
           rem set first=0
           rem ECHO !first!
           rem ECHO %%a
           set one=%%a
           SET ONE=!one:~1,7!
           GOTO :RAR
      )
      set /a first=first+1
)
rem Echo !one!
:RAR
rar.exe a -u -mc10:20t+ -Mm5 -ag-YYYYMMDD-NN -xrar.exe -x*.mzML -xscopy.exe -x*.lib backup/HDXconsole_%one%.rar  *.*
@ECHO OFF
setlocal enabledelayedexpansion
set first=1
set one = "_"
for /f "tokens=5" %%a in (constants.h) do (
      if !first! == 1 (
           rem set first=0
           rem ECHO !first!
           ECHO %%a
           set one=%%a
           SET ONE=!one:~1,7!
           rem GOTO :EOF
      )
      set /a first=first+1
)
Echo !one!
@echo off
setlocal
cd /d "%~dp0"
if /I "%~1"=="--yes" goto AUTO
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0tools\setup_windows.ps1"
set "RC=%ERRORLEVEL%"
echo.
pause
exit /b %RC%
:AUTO
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0tools\setup_windows.ps1" -Yes
exit /b %ERRORLEVEL%

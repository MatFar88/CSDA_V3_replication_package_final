@echo off
setlocal
cd /d "%~dp0"
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0tools\setup_windows.ps1"
if errorlevel 1 goto FAILED
"%~dp0.venv\Scripts\python.exe" "%~dp0run_replication.py" REPORT
if errorlevel 1 goto FAILED
echo.
echo REPORT completed. Open REPORT.html at the path printed above.
pause
exit /b 0
:FAILED
echo.
echo Setup or REPORT failed. No new fit was started.
pause
exit /b 1

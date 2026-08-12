@echo off
if "%~2"=="" (
  echo Usage: Import-Counts.cmd COUNT_FILE PROJECT_ID
  pause
  exit /b 2
)
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0Import-Counts.ps1" "%~1" "%~2"

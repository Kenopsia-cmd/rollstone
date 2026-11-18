@echo off
setlocal
cd /d "%~dp0"
if not exist build\rollingstone.exe call build.cmd
if errorlevel 1 (
  pause
  exit /b 1
)
build\rollingstone.exe %*
if errorlevel 1 pause

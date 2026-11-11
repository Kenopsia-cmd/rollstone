@echo off
setlocal
cd /d "%~dp0"
set "VSROOT="
for /f "usebackq tokens=*" %%i in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VSROOT=%%i"
if not defined VSROOT (
  echo Visual Studio C++ tools were not found.
  exit /b 1
)
call "%VSROOT%\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 exit /b 1
if not exist build mkdir build
cl /nologo /EHsc /std:c++17 /utf-8 /W4 /O2 /D_UNICODE /DUNICODE /DNOMINMAX /I"%VSROOT%\VC\Auxiliary\VS\include" rollingstone.cpp /Fe:build\rollingstone.exe /Fo:build\rollingstone.obj /link /LIBPATH:"%VSROOT%\VC\Auxiliary\VS\lib\x64" EasyXw.lib user32.lib gdi32.lib shell32.lib
exit /b %errorlevel%

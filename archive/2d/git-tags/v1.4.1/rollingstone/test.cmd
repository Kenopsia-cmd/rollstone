@echo off
setlocal
cd /d "%~dp0"
set "VSROOT="
for /f "usebackq tokens=*" %%i in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VSROOT=%%i"
if not defined VSROOT exit /b 1
call "%VSROOT%\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 exit /b 1
if not exist build mkdir build
cl /nologo /EHsc /std:c++17 /utf-8 /O2 batch_tests.cpp /Fe:build\batch_tests.exe /Fo:build\batch_tests.obj
if errorlevel 1 exit /b 1
cl /nologo /EHsc /std:c++17 /utf-8 /O2 physics_tests.cpp /Fe:build\physics_tests.exe /Fo:build\physics_tests.obj
if errorlevel 1 exit /b 1
cl /nologo /EHsc /std:c++17 /utf-8 /O2 studies_tests.cpp /Fe:build\studies_tests.exe /Fo:build\studies_tests.obj
if errorlevel 1 exit /b 1
pushd build
batch_tests.exe
if errorlevel 1 goto fail
physics_tests.exe
if errorlevel 1 goto fail
studies_tests.exe
if errorlevel 1 goto fail
popd
echo All C++ checks passed.
exit /b 0
:fail
popd
exit /b 1

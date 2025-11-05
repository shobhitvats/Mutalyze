@echo off
REM Mutalyze Executable Builder for Windows
REM This script builds a standalone .exe with full OpenMM support

echo ========================================
echo  MUTALYZE EXECUTABLE BUILDER
echo ========================================
echo.

REM Check if conda is available
where conda >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Conda not found! Please install Miniconda first.
    echo Download from: https://docs.conda.io/en/latest/miniconda.html
    pause
    exit /b 1
)

REM Check if environment exists
conda env list | find "mutalyze_build" >nul
if %ERRORLEVEL% NEQ 0 (
    echo Creating conda environment...
    conda create -n mutalyze_build python=3.10 -y
)

echo.
echo Activating environment...
call conda activate mutalyze_build

echo.
echo Installing OpenMM...
conda install -c conda-forge openmm -y

echo.
echo Installing Python dependencies...
pip install -r requirements_exe.txt

echo.
echo Verifying installations...
python -c "import openmm; print(f'OpenMM {openmm.version.version} OK')"
python -c "import PyInstaller; print(f'PyInstaller {PyInstaller.__version__} OK')"

echo.
echo ========================================
echo  BUILDING EXECUTABLE
echo ========================================
echo.
echo This may take 5-10 minutes...
echo.

REM Build using spec file
pyinstaller build_exe_advanced.spec

echo.
if exist "dist\Mutalyze.exe" (
    echo ========================================
    echo  BUILD SUCCESSFUL!
    echo ========================================
    echo.
    echo Executable location: dist\Mutalyze.exe
    echo.
    echo File size:
    dir "dist\Mutalyze.exe" | find "Mutalyze.exe"
    echo.
    echo To test, run: dist\Mutalyze.exe
    echo.
) else (
    echo ========================================
    echo  BUILD FAILED!
    echo ========================================
    echo.
    echo Check the output above for errors.
    echo Common issues:
    echo   - OpenMM not installed
    echo   - Missing dependencies
    echo   - File path too long
    echo.
)

pause

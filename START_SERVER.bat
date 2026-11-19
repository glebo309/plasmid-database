@echo off
REM Plasmid Library Server Startup Script for Windows

echo ============================================================
echo 🧬 Starting Plasmid Library Server...
echo ============================================================

REM Check if Python 3 is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo ❌ Python is not installed!
    echo Please install Python 3 from https://www.python.org/downloads/
    pause
    exit /b 1
)

REM Check if Flask is installed
python -c "import flask" >nul 2>&1
if errorlevel 1 (
    echo 📦 Flask not found. Installing...
    pip install -r requirements.txt
)

REM Start the server
python plasmid_server.py
pause

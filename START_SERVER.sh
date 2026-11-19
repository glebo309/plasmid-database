#!/bin/bash
# Plasmid Library Server Startup Script

echo "============================================================"
echo "🧬 Starting Plasmid Library Server..."
echo "============================================================"

# Check if Python 3 is installed
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is not installed!"
    echo "Please install Python 3 from https://www.python.org/downloads/"
    exit 1
fi

# Check if Flask is installed
if ! python3 -c "import flask" &> /dev/null; then
    echo "📦 Flask not found. Installing..."
    pip3 install -r requirements.txt
fi

# Start the server
python3 plasmid_server.py

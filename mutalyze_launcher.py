"""
Mutalyze Desktop Launcher
Standalone launcher that starts Streamlit server and opens browser
Works as standalone executable with full OpenMM support
"""

import streamlit.web.cli as stcli
import sys
import os
from pathlib import Path
import socket
import time
import webbrowser
import threading

def find_free_port():
    """Find an available port for Streamlit"""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('', 0))
        s.listen(1)
        port = s.getsockname()[1]
    return port

def open_browser(port, delay=2):
    """Open browser after delay to let server start"""
    time.sleep(delay)
    webbrowser.open(f'http://localhost:{port}')

def main():
    """Main launcher function"""
    
    # Get the directory where the executable is located
    if getattr(sys, 'frozen', False):
        # Running as compiled executable
        application_path = Path(sys.executable).parent
    else:
        # Running as script
        application_path = Path(__file__).parent
    
    # Change to application directory
    os.chdir(application_path)
    
    # Find the app.py file
    app_file = application_path / 'app.py'
    if not app_file.exists():
        print(f"❌ ERROR: app.py not found at {app_file}")
        input("Press Enter to exit...")
        sys.exit(1)
    
    # Find free port
    port = find_free_port()
    
    print("="*70)
    print("🧬 MUTALYZE - Protein Mutation Stability Analysis Platform")
    print("="*70)
    print(f"\n✅ Starting server on port {port}...")
    print(f"📂 Working directory: {application_path}")
    print(f"🔬 OpenMM support: Enabled (full accuracy)")
    print(f"\n🌐 Browser will open at: http://localhost:{port}")
    print("\n💡 To stop: Close this window or press Ctrl+C")
    print("="*70 + "\n")
    
    # Open browser in separate thread
    browser_thread = threading.Thread(target=open_browser, args=(port,))
    browser_thread.daemon = True
    browser_thread.start()
    
    # Start Streamlit
    sys.argv = [
        "streamlit",
        "run",
        str(app_file),
        f"--server.port={port}",
        "--server.headless=true",
        "--browser.gatherUsageStats=false",
        "--server.fileWatcherType=none",  # Disable file watcher in exe
    ]
    
    sys.exit(stcli.main())

if __name__ == '__main__':
    main()

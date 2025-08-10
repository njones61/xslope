#!/usr/bin/env python3
"""
Script to run mkdocs serve with warnings suppressed.
This avoids all the type annotation warnings from griffe.
"""

import warnings
import subprocess
import sys

# Suppress all warnings
warnings.filterwarnings("ignore")

# Suppress specific griffe warnings
import logging
logging.getLogger("griffe").setLevel(logging.ERROR)

# Run mkdocs serve
try:
    result = subprocess.run([sys.executable, "-m", "mkdocs", "serve"], 
                          check=False, capture_output=False)
    sys.exit(result.returncode)
except KeyboardInterrupt:
    print("\nMkDocs stopped by user")
    sys.exit(0)
except Exception as e:
    print(f"Error running mkdocs: {e}")
    sys.exit(1) 
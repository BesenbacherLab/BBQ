"""
Entry-point module, in case you use `python -m betterbasequals`.
"""

import sys

from betterbasequals.cli import main

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

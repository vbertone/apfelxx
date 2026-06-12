#!/usr/bin/env python3
"""
Non-regression test: run a test executable and compare its filtered output
against a stored benchmark file.

Usage: regression.py <executable> <benchmark>
Exit 0 = pass, Exit 1 = fail or error
"""

import sys
import subprocess
import re
import os
import difflib


def filter_output(text):
    """Remove machine-dependent noise from test output."""
    # Strip ANSI escape codes (colour warnings, etc.)
    text = re.sub(r'\x1b\[[0-9;]*[A-Za-z]', '', text)
    # Strip timing suffixes: "... Time elapsed: 0.123456 seconds"
    text = re.sub(r'\s*Time elapsed:\s*[\d.]+\s*seconds', '', text)
    # Normalise memory addresses (e.g. "Grid: 0x16dd8ee40"), which change run to run
    text = re.sub(r'\b0x[0-9a-fA-F]+\b', '0x*', text)
    return text


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <executable> <benchmark>", file=sys.stderr)
        sys.exit(1)

    executable, benchmark = sys.argv[1], sys.argv[2]
    name = os.path.basename(executable)

    if not os.path.exists(executable):
        print(f"SKIP: {name} — executable not found: {executable}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(benchmark):
        print(f"ERROR: {name} — benchmark not found: {benchmark}", file=sys.stderr)
        print("Run  python3 tests/generate_benchmarks.py <build_tests_dir> tests/benchmarks",
              file=sys.stderr)
        sys.exit(1)

    result = subprocess.run([executable], capture_output=False,
                            stdout=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"FAIL: {name} — exited with code {result.returncode}", file=sys.stderr)
        sys.exit(1)

    actual = filter_output(result.stdout)

    with open(benchmark) as f:
        expected = f.read()

    if actual == expected:
        print(f"PASS: {name}")
        sys.exit(0)

    print(f"FAIL: {name}")
    diff = difflib.unified_diff(
        expected.splitlines(keepends=True),
        actual.splitlines(keepends=True),
        fromfile='benchmark',
        tofile='actual',
        n=3,
    )
    sys.stdout.writelines(diff)
    sys.exit(1)


if __name__ == '__main__':
    main()

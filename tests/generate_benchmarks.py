#!/usr/bin/env python3
"""
Generate (or refresh) benchmark reference files for the non-regression suite.

Runs every test executable found in <build_tests_dir>, strips machine-dependent
noise (ANSI codes, timing lines), and writes the result to
<benchmarks_dir>/<test_name>.out.

Usage:
    python3 tests/generate_benchmarks.py <build_tests_dir> tests/benchmarks

Typical one-liner from the repo root after building:
    python3 tests/generate_benchmarks.py build/tests tests/benchmarks
"""

import sys
import subprocess
import re
import os


def filter_output(text):
    """Remove machine-dependent noise from test output."""
    text = re.sub(r'\x1b\[[0-9;]*[A-Za-z]', '', text)
    text = re.sub(r'\s*Time elapsed:\s*[\d.]+\s*seconds', '', text)
    # Normalise memory addresses (e.g. "Grid: 0x16dd8ee40"), which change run to run
    text = re.sub(r'\b0x[0-9a-fA-F]+\b', '0x*', text)
    return text


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <build_tests_dir> <benchmarks_dir>")
        sys.exit(1)

    build_dir, benchmarks_dir = sys.argv[1], sys.argv[2]

    if not os.path.isdir(build_dir):
        print(f"Error: build directory not found: {build_dir}")
        sys.exit(1)

    os.makedirs(benchmarks_dir, exist_ok=True)

    # Test names are derived from the .cc sources next to this script,
    # mirroring the glob in tests/CMakeLists.txt
    tests_dir = os.path.dirname(os.path.abspath(__file__))
    names = sorted(
        os.path.splitext(f)[0] for f in os.listdir(tests_dir) if f.endswith('.cc')
    )

    executables = []
    for name in names:
        path = os.path.join(build_dir, name)
        if os.path.isfile(path) and os.access(path, os.X_OK):
            executables.append(name)
        else:
            print(f"Warning: executable not found (not built?): {path}")

    if not executables:
        print(f"No test executables found in {build_dir}")
        sys.exit(1)

    print(f"Generating benchmarks for {len(executables)} tests in {benchmarks_dir}/\n")

    failed = []
    for name in executables:
        path = os.path.join(build_dir, name)
        print(f"  {name} ...", end='', flush=True)
        result = subprocess.run([path], stdout=subprocess.PIPE, text=True)
        filtered = filter_output(result.stdout)
        out_path = os.path.join(benchmarks_dir, name + '.out')
        with open(out_path, 'w') as f:
            f.write(filtered)
        lines = filtered.count('\n')
        status = f" {lines} lines"
        if result.returncode != 0:
            status += f"  [exit {result.returncode}]"
            failed.append(name)
        print(status)

    print(f"\nDone. Benchmarks written to: {benchmarks_dir}")
    if failed:
        print(f"Warning: {len(failed)} test(s) exited non-zero: {', '.join(failed)}")


if __name__ == '__main__':
    main()

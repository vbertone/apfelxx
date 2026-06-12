#!/usr/bin/env python3
"""
Non-regression test: run a test executable and compare its filtered output
against a stored benchmark file.

The comparison is numeric-aware: floating-point numbers appearing in the
output are compared with a relative and an absolute tolerance, so that
platform-dependent rounding noise (different architectures, libm, FMA usage)
does not produce spurious failures. Non-numeric text must match exactly.

Usage: regression.py <executable> <benchmark>
Exit 0 = pass, Exit 1 = fail or error
"""

import sys
import subprocess
import re
import os
import math

# Tolerances for numeric comparison. The absolute tolerance absorbs
# "zero plus rounding noise" values (e.g. 1.7e-15 vs 1.1e-15 for sum rules
# that vanish identically); the relative tolerance bounds genuine values.
REL_TOL = 1e-5
ABS_TOL = 1e-10

# Matches a floating-point or integer literal, possibly with exponent
NUMBER = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')


def filter_output(text):
    """Remove machine-dependent noise from test output."""
    # Strip ANSI escape codes (colour warnings, etc.)
    text = re.sub(r'\x1b\[[0-9;]*[A-Za-z]', '', text)
    # Strip timing suffixes: "... Time elapsed: 0.123456 seconds"
    text = re.sub(r'\s*Time elapsed:\s*[\d.]+\s*seconds', '', text)
    # Normalise memory addresses (e.g. "Grid: 0x16dd8ee40"), which change run to run
    text = re.sub(r'\b0x[0-9a-fA-F]+\b', '0x*', text)
    return text


def compare_lines(expected, actual):
    """
    Compare two lines treating numbers numerically and the rest textually.
    Return None if they match, otherwise a string describing the mismatch.
    """
    if expected == actual:
        return None

    exp_numbers = NUMBER.findall(expected)
    act_numbers = NUMBER.findall(actual)

    # The non-numeric skeleton of the line must be identical
    if NUMBER.sub('#', expected) != NUMBER.sub('#', actual):
        return "text differs"

    if len(exp_numbers) != len(act_numbers):
        return f"number count differs ({len(exp_numbers)} vs {len(act_numbers)})"

    for e, a in zip(exp_numbers, act_numbers):
        if not math.isclose(float(e), float(a), rel_tol=REL_TOL, abs_tol=ABS_TOL):
            return f"value {a} not within tolerance of benchmark {e}"

    return None


def compare_outputs(expected, actual):
    """
    Compare two whole outputs line by line. Return a list of mismatch
    descriptions (empty list means the outputs are equivalent).
    """
    exp_lines = expected.splitlines()
    act_lines = actual.splitlines()

    mismatches = []
    if len(exp_lines) != len(act_lines):
        mismatches.append(f"line count differs: benchmark has {len(exp_lines)}, "
                          f"actual has {len(act_lines)}")

    for i, (e, a) in enumerate(zip(exp_lines, act_lines), start=1):
        problem = compare_lines(e, a)
        if problem is not None:
            mismatches.append(f"line {i}: {problem}\n"
                              f"  benchmark: {e}\n"
                              f"  actual:    {a}")

    return mismatches


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

    mismatches = compare_outputs(expected, actual)

    if not mismatches:
        print(f"PASS: {name}")
        sys.exit(0)

    print(f"FAIL: {name} — {len(mismatches)} mismatch(es)")
    # Cap the report so a badly broken test does not flood the log
    for m in mismatches[:50]:
        print(m)
    if len(mismatches) > 50:
        print(f"... and {len(mismatches) - 50} more")
    sys.exit(1)


if __name__ == '__main__':
    main()

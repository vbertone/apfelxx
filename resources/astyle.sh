#!/bin/bash

astyle --verbose --recursive --suffix=none --style=gnu --indent=spaces=2 --indent-namespaces --max-continuation-indent=120 "../src/*.cc"
astyle --verbose --recursive --suffix=none --style=gnu --indent=spaces=2 --indent-namespaces --max-continuation-indent=120 "../tests/*.cc"
astyle --verbose --recursive --suffix=none --style=gnu --indent=spaces=2 --indent-namespaces --max-continuation-indent=120 "../tests/ACOT_tests/*.cc"
astyle --verbose --recursive --suffix=none --style=gnu --indent=spaces=2 --indent-namespaces --max-continuation-indent=120 "../pywrap/*.cc"
astyle --verbose --recursive --suffix=none --style=gnu --indent=spaces=2 --indent-namespaces --max-continuation-indent=120 --keep-one-line-blocks "../inc/apfel/*.h"

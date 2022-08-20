#!/bin/bash

astyle --recursive --suffix=none --style=gnu --indent=spaces=2 --indent-namespaces --max-continuation-indent=120 "../src/*.cc"
astyle --recursive --suffix=none --style=gnu --indent=spaces=2 --indent-namespaces --max-continuation-indent=120 "../tests/*.cc"
astyle --recursive --suffix=none --style=gnu --indent=spaces=2 --indent-namespaces --max-continuation-indent=120 "../pywrap/*.cc"
astyle --recursive --suffix=none --style=gnu --indent=spaces=2 --indent-namespaces --max-continuation-indent=120 --keep-one-line-blocks "../inc/apfel/*.h"

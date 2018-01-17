#!/bin/sh
cmake -E make_directory build/Release
cd build/Release
cmake -DCMAKE_BUILD_TYPE=Release -G "CodeBlocks - Unix Makefiles" ../../
cmake --build . --target DMRG++ -- -j 4


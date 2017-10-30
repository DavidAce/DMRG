#!/bin/sh
cmake -E make_directory build/Release
cd build/Release
cmake -j4  -Bbuild/Release --build build -config Release ../../

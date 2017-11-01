#!/bin/sh
cmake -E make_directory build/Release
cd build/Release
cmake -Bbuild/Release --build build -config Release ../../
make


#!/bin/bash

cmake --preset=clang-static \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=build/install \
      -DDMRG_PACKAGE_MANAGER=conan \
       .

#cmake --build . --parallel 4
#ctest -C Debug --output-on-failure
#!/bin/bash
  
# build on Mac, from cmake directory run this as ./build.mac.bash

cmake -B ../cmake_build -S . -DCMAKE_BUILD_TYPE=Sundials_Debug
cmake --build ../cmake_build --target all


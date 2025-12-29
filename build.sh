#!/usr/bin/env bash
if [ ! -d build ]; then
  mkdir build
  cd build
  cmake ..
  cd ..
fi
make -C build

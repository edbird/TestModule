#!/bin/bash
cmake -DCMAKE_PREFIX_PATH="/home/ecb/Falaise-local;$(brew --prefix);$(brew --prefix qt5-base)" $@
#cmake -DCMAKE_PREFIX_PATH="$(brew --prefix);$(brew --prefix qt5-base)" $@

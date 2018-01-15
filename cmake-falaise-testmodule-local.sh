#!/bin/bash
cmake -DCMAKE_PREFIX_PATH="/home/ecb/Falaise-local;$(brew --prefix);$(brew --prefix qt5-base)" $@

#!/bin/bash
(cd ./data && ./run-flsimulate.sh -n 1000 && ./run-flreconstruct.sh)
./run-test-module.sh

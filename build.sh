#!/bin/sh
g++ src/collision.cpp src/matrix.cpp src/physics_simulation.cpp src/potential.cpp src/timestep.cpp -lfftw3_threads -lfftw3 -lm -O3 -o helios

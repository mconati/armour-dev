g++ -o test PZ_tests.cpp Dynamics.cpp Trajectory.cpp PZsparse.cpp -fopenmp -std=c++14 -O2 -lstdc++ -ldl -lm -lc -lgomp

nvcc -o collision_checker collision_checker.cu FastCollisionChecking.cu -Xcompiler -fopenmp -std=c++14 -O3 -lquadmath -lstdc++ -ldl -lm -lc -lgomp

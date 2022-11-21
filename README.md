# armour-dev
Can’t Touch This: Real-Time, Safe Motion Planning and Control for Manipulators Under Uncertainty

## Start
The repo has been verified on Matlab R>=2021b and Ubuntu >= 20.04

Run 
 - initialize.m
 - kinova_src/initialize.m
in Matlab before you run any other scripts!

### Dependency
This repo depends on the following repos:
 - [CORA](https://tumcps.github.io/CORA/)

This repo assumes that you have installed the following libraries:
 - libboost-dev (1.71)
 - libeigen3-dev (3.3.7) in eigen3/
 - libipopt
 - libcoinhsl
 
 libipopt and libcoinhsl could be very annoying to install and to work with Matlab. 
 Suppose libipopt and libcoinhsl are both installed in /usr/local/lib.
 You need to add that path to both user's environmental variable 'LD_LIBRARY_PATH' and Matlab's environment variable 'LD_LIBRARY_PATH'
 Check [here](https://www.mathworks.com/help/matlab/matlab_external/set-run-time-library-path-on-linux-systems.html) and [here](https://stackoverflow.com/questions/13428910/how-to-set-the-environmental-variable-ld-library-path-in-linux) for more information.

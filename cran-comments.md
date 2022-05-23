2022-05-18
Valgrind detected some errors in the previous version. There's a 'vertex' type 
in my code to which one has to associate some indices. I associated indices 
starting from 1. My guess is that this is the cause of the errors and the 
indices must start from 0. I did the correction.

2022-05-20
I found another to enter the coordinates of the vertices and their index. I 
hope this will solve the issues spotted by Valgrind.

2022-05-23
There was a declared but uninitialized (and unused) vector in the C++ code. 
Perhaps it was the cause of the Valgrind issues. There also was a missing 
template argument in a function.

**Note:** The latest version of the 'rgl' package has serious issues.


## Testing environments

- Local R-4.1.2 on Windows 10
- win-builder devel
- mac-builder
- Ubuntu 20 via Github action


## R CMD check results

OK

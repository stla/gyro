The previous submission has been archived because of an issue spotted by Valgrind 
(triggered by a call to the 'hdelaunay' function). 
Actually, my C++ code is correct: this issue occurs because of the CGAL library, 
not because of my code. So, I have no way to fix it, and then I removed the 'hdelaunay' 
function from the package. Now the package does not use CGAL. 

## Testing environments

- Local R-4.1.2 on Windows 10
- win-builder devel
- mac-builder
- Ubuntu 20 via Github action


## R CMD check results

OK

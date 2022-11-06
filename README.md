# ode-tokimasu
A ode solver repository in C, C++, MATLAB and Python. Courtesy of "Numerical Methods For Solution of Differential Equations" by Tobias Ritschel.

# Quickstart
1) Install 
```
git clone https://github.com/smallpondtom/ode-tokimasu.git
```
2) Install the Eigen 
```
git submodule init
git submodule update
```


## ToDo
1) Implement DOP853 algorithm from Hairer.
2) Implement ESDIRK23 algorithm (implicit method) for stiff ODEs.
3) Understand and implement dense output for each algorithm to reduce computation cost.
4) Write ESDIRK23 method function for python used for stiff-functions.
5) Create code for C.
6) Create code for C++.
7) Create code for Julia.

# MonteCarloGeometryProcessing

## Introduction
Course project for [CSC2520 Geometry Processing](https://github.com/alecjacobson/geometry-processing-csc2520).

An implementation based on the paper by Rohan Sawhney and Keenan Crane \
["Monte Carlo Geometry Processing: A Grid-Free Approach to PDE-Based Methods on Volumetric Domains"](https://www.cs.cmu.edu/~kmcrane/Projects/MonteCarloGeometryProcessing/paper.pdf)

## Compilation
This project makes use of open-source libraries `libigl` and `Eigen` which will be included by cloning `libigl` repository.

To clone the repository, use 
```
git clone https://github.com/hgeorge21/MonteCarloGeometryProcessing.git
cd MonteCarloGeometryProcessing
git clone https://github.com/libigl/libigl.git
```

To build the project, go into the project directory and do the following
```
mkdir build
cd build
cmake --CMAKE_TYPE=Release ..
```

## Project Struction
TODO

## TODO Item List
1. Figure out how to define the boundary for input mesh... By user or using boundary loop?
2. What PDE to support? (Prob. list them out and let user choose)
	- $\delta u = 0, u = ||x|| on \del\omega$ (Laplace with Dirichlet BD)
    - $\delta u = g, u = 0 on \del\omega$ (Poisson, need to think about $g$)
    - Maybe biharmonic? (can show the parameterization using this biharmonic tehnique LOL)
    - let's have a couple gradient and hessian in there as well :)
3. Add the data folder, find some good meshes
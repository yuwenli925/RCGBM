# RCGBM
This repo is a collection of MATLAB/Octave codes for the paper "A reduced conjugate gradient basis method for fractional diffusion" by Yuwen Li, Ludmil Zikatanov, and Cheng Zuo.

The mesh refinement functions uniformrefine.m, uniformrefine3.m and auxiliary structure functions auxmesh3.m, gradbasis3.m are borrowed from the iFEM package developed Prof. Long Chen at UC Irvine.

OGA.m is to construct rational approximation for the fractional power function z^{-s} by the orthogonal greedy algorithm.

FracRBM3.m is the main file about the example of the RCGBM for solving fractional Laplacian on the unit cube.

FracSurface.m is the main file about the example of the RCGBM for solving fractional Laplacian on the unit sphere

FracGraphLaplacian.m is the main file about example of the RCGBM for solving fractional Laplacian on a random graph.



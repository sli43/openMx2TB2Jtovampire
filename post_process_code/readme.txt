This code written by Shaozhi Li, 09/17/2022
This code is used to compute dynamical magnetic structure
Please use makefile to compile this code. Please change the path of fftw library in the makefile


One needs to change the read-in path on line 140 and 216
Line 140: define the path for the spin configuration
Line 216: define the path for the atomic coordinate

This code does not compute the dynamical magnetic structure in the whole momentum space.
Instead, this code only computes the dynamical magnetic structure along a defined momentum path
The momentum path is defined from Line 162 to 202
The total momentum number is given by Nkpoint, which is defined on line 15

Output: spectral.dat

Plot: use plot_spectra.py to plot the dynamical magnetic structure
% python plot_spectra.py

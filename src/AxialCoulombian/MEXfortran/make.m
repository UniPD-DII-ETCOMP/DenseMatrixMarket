clear 
close all
clc
mex -O -largeArrayDims -output funLphiphi3_for_st L_funLphiphi3_mex.f90 funLphiphi3.f90 

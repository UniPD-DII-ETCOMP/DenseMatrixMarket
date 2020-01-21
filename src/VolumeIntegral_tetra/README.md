# VolumeIntegral code 

This directory contains a Volume Integral code for electromagnetic problems,


In this code, simply connected conductive domains can be considered. 
Tetrahedral elements are considered for the mesh discretization.
External fields can be used for the excitation.

-------------------------------------------------------------------

# Description
 
MAIN_dense.m is the main file you must run to start the code. 
                      The final system of equations is generated and then solved by Matlab backslash

-------------------------------------------------------------------

MAIN_to_low_rank.m is the main file you must run to start the code coupled with the low rank compression library https://github.com/numpi/hm-toolbox based on HODLR and HSS. HSS and HODLR methods can be used to compress the matrix.
Note that the compression performances are highly problem dependent and the compression may be very poor without a proper reordering of the DoFs.
In this regard, "fun_reo_DMM.m" is provided and it can be used to reorder the DoFs 
In order to test the application of HODLR/HSS to the VolumeIntegral  formulation follow these steps:

1. Download the Matlab hm-toolbox https://github.com/numpi/hm-toolbox
2. Extract it so that the "hm-toolbox-master" directory is at the same level of "test_cases", "fun", and "MEXfortran" 		 
3. Execute "MAIN_to_low_rank.m"

Note that this is only to show how to combine HODLR/HSS with the VolumeIntegral problems in a matrix-free form, i.e. the system matrix is never fully assembled/stored. 
The actual compression ratio depends on the specific problem features, in particular the ordering of the unknowns.

All user-settable quantities, e.g. frequency and resistivity, are contained in the block identified by the 
BEGIN USER SETTINGS / END USER SETTINGS comments.
-------------------------------------------------------------------

Available test cases
--------------------
Several test cases are contained in separate directories under "test_cases". 
Set the "test_case_dir" variable in "MAIN_dense.m" or "MAIN_to_low_rank.m " to the appropiate directory.

User-defined test cases
-----------------------
Follow the instuctions given in "README.txt" inside the "test_cases" directory.

Contacts
-----------------------
Riccardo Torchio (riccardo.torchio@unipd.it)

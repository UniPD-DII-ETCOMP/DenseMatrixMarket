# AxialCoulombian code by Riccardo Torchio (riccardo.torchio@studenti.unipd.it)

This directory contains an Axialsymmetric code for electromagnetic problems based on the Coulombian formulation presented in

* [R. Torchio et al., "Volume Integral Equation Methods for Axisymmetric Problems with Conductive and Magnetic Media," in IEEE Transactions on Magnetics. in press]()

A more general, optimized, parallel (OpenMP), fortran90 version of this code has been used in the above paper


In this code, conductive and magnetic (homogeneous) domains can be considered. 
Quadrilateral elements are considered for the mesh discretization.
External fields and current driven coils can be used for the excitation.

See the above references for more details and consider citing it.

-------------------------------------------------------------------

% Description
 
MAIN_dense.m is the main file you must run to start the code. 
                      The final system of equations is generated and then solved by Matlab backslash


MAIN_to_low_rank_cond_only.m is the main file you must run to start the code coupled with the low rank compression library https://github.com/numpi/hm-toolbox based on HODLR and HSS. 
                   		 HSS and HODLR methods can be used o compress the matrix.
                    		 Note that the compression performances are higly problem dependent and the compression may be very poor without a proper reordering of the DoFs.
		   		 In order to test the application of HSS to the Axysimmetric Coulombian formulation follow these steps:
	            			1. Download the Matlab hm-toolbox https://github.com/numpi/hm-toolbox
	            			2. Extract it so that the "hm-toolbox-master" directory is at the same level of "test_cases", "fun" 		 
		    			3. Execute "MAIN_to_low_rank_cond_only.m"
		    		 Note that this is only to show how to combine HODLR/HSS with the Axisymmetric problems in a matrix-free form, 
                   		 i.e. the system matrix is never fully assembled/stored. 
                   		 The actual compression ratio depends on the specific problem features, in particular the ordering of the unknowns.
                   		 WARNING: this version of the code only support conductive media (for the moment)

All user-settable quantities, e.g. frequency and resistivity, are contained in the block identified by the 
BEGIN USER SETTINGS / END USER SETTINGS comments.
-------------------------------------------------------------------

Available test cases
--------------------
Several test cases are contained in separate directories under "test_cases". 
Set the "test_case_dir" variable in "MAIN_dense.m" or "MAIN_to_low_rank_cond_only.m " to the appropiate directory.

User-defined test cases
-----------------------
Follow the instuctions given in "README.txt" inside the "test_cases" directory.

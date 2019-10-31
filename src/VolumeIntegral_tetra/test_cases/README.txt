Dense Matrix Market: VolumeIntegral code by Riccardo Torchio (riccardo.torchio@studenti.unipd.it)
-------------------------------------------------------------------------------------------------------
How to create a new user-defined test-case:

1. User must provide file "geo.mat" which contains:
	a. Matrix_P0 = Npoints x 3 (the coordianate of the Npoints nodes of the tetrahedral mesh) 
	b. VP = 4 x Ntetra incidence matrix (tetrahedra to nodes of the mesh) 

2. duplicate the directory "test_from_user" and rename it (e.g. "new_user_dir")
3. copy the user's file "geo.mat" inside "new_user_dir"
4. run "MAIN_geometry_runme_from_geo"

File "data.mat" is now created and you can select "new_user_dir" in "MAIN_dense.m" 
or "MAIN_to_low_rank.m" to select the user test case
-------------------------------------------------------------------------------------------------------

Subdirectories contain some examples of data generation:

"SphereShell1866", "SphereShell8940", "test1", "test2", "testTK" and "testDMM": geometrical data are taken from a COMSOL mesh file

"test_from_user": contains a sample user defined "geo.mat" file 

"test_spherical_shell_from_user": generate a user defined test case of a spherical shell



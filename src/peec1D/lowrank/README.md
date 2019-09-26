In order to test the application of HSS to the PEEC "stick" formulation follow these steps:
1. [Download the Matlab HM-toolbox](https://github.com/numpi/hm-toolbox)
2. Extract it so that the "hm-toolbox-master" directory is at the same level of "test_cases", "fun" and "lowrank"
3. Replace the file "hss.m" in "hm-toolbox-master/@hss" with the one provided in the "lowrank" directory
4. Execute "MAIN_PEEC_sticks_MF_HSS.m"

Note that this is only to show how to combine  HSS with the PEEC "stick" formulation in a matrix-free form, i.e. the system matrix is never fully assembled/stored. The actual compression ratio depends on the specific problem features, in particular the ordering of the unknowns.

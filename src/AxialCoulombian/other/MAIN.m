%% Main_AXIAL_EH
clear
close all
clear global
clc
%% load data
disp('-------------------------------------------------------------------')
main_load_object_data
disp('-------------------------------------------------------------------')
%% SIMULATION DATA
disp('-------------------------------------------------------------------')
disp('Simulation data set... ')
f=50; % set frequency
disp('... done!')
disp('-------------------------------------------------------------------')
%% Physical parameter
disp('-------------------------------------------------------------------')
disp('Physical parameters... ')
main_evalualte_physical_parameter;
disp('-------------------------------------------------------------------')
%% External field
disp('-------------------------------------------------------------------')
ext_Field_ext=0; %0=external field off ||| 1=extarnal field on
main_ext_field; % set the external field and the extarna coil current
if ext_Field_ext==1
    disp('External Field ON!')
else
    disp('External Field OFF!')    
end
disp('-------------------------------------------------------------------')
%% Quick change of simulation data (if they are not loaded yet)
main_set_material; % set the material of condcutive and magnetic domains
%% How many threads?
N.thread = 22; 
%% data_handling
disp('-------------------------------------------------------------------')
disp('DATA HANDLING...')
main_data_handling; % handling of geometry data 
disp('... done!')
disp('-------------------------------------------------------------------')
%% disp
disp('-------------------------------------------------------------------')
disp('SIMULATION DATA... ')
disp(['...frequency = ' num2str(f) ' [Hz]'])
disp(['...N_face_con = ',num2str(N.face_con),'...'])
disp(['...N_edge_con = ',num2str(N.edge_con),'...'])
disp(['...N_node_con = ',num2str(N.node_con),'...'])
disp(['...N_face_mag = ',num2str(N.face_mag),'...'])
disp(['...N_edge_mag = ',num2str(N.edge_mag),'...'])
disp(['...N_node_mag = ',num2str(N.node_mag),'...'])
disp('-------------------------------------------------------------------')
%% set material vector
disp('-------------------------------------------------------------------')
disp('Set material...')
main_set_material_vec; % set the value of rho_c and rho_m for each mesh element
disp('... done!')
disp('-------------------------------------------------------------------')
%% Boundary plot
disp('-------------------------------------------------------------------')
disp('geo-boundary plot...')
main_plot_boundary;
pause(0.5)
disp('... done!')
disp('-------------------------------------------------------------------')
%% Computing matrices RLP
disp('-------------------------------------------------------------------')
disp('computing matrices RLP...')
main_compute_matrices2; % compute matrices R L P N 
disp('...RLP done!')
disp('-------------------------------------------------------------------')
%% Computing matrices incidances
disp('-------------------------------------------------------------------')
disp('computing matrices GCD...')
main_computing_incidance_matrices; % compute incidence matrices
disp('... done!')
disp('-------------------------------------------------------------------')
%% size storage system
disp('-------------------------------------------------------------------')
size_SISTEM=N.face_con+N.edge_mag;
disp(['Size system= ',num2str(size_SISTEM)])
disp('-------------------------------------------------------------------')
%% Assebling System
disp('-------------------------------------------------------------------')
disp('Assembling system...')
main_assembling_system; % assembling the system of equations
disp('... done!')
disp('-------------------------------------------------------------------')
%% size storage system
disp('-------------------------------------------------------------------')
size_SISTEM=size(SYSTEM,1);
disp(['Size system = ',num2str(size_SISTEM)])
disp('-------------------------------------------------------------------')
%% Rhs-term External field
disp('-------------------------------------------------------------------')
disp('External coils...')
main_compute_ext_fiel_coils; % rhs: contribution of external coils
disp('... done!')
disp('-------------------------------------------------------------------')
%% Rhs-term External field
disp('-------------------------------------------------------------------')
disp('External fields...')
main_compute_rhs_ext_field; % rhs: contribution of external fields
disp('... done!')
disp('-------------------------------------------------------------------')
%% Change of variables
disp('-------------------------------------------------------------------')
disp('Change of variables (jc,m)-->(jc,phim)...')
main_change_of_variable; % project the system into a new set of equations
disp('... done!')
disp('-------------------------------------------------------------------')
%% size  system
disp('-------------------------------------------------------------------')
size_SISTEM2=size(SYSTEM2,1);
disp(['Size system= ',num2str(size_SISTEM2)])
disp('-------------------------------------------------------------------')
%% SOLVING
disp('-------------------------------------------------------------------')
disp('Solving ...')
tic
main_solving; % solve the system of equation
toc
disp('... done!')
disp('-------------------------------------------------------------------')
%% Extracting solution
disp('-------------------------------------------------------------------')
disp('Extracting solution ...')
tic
main_extracting_solution; % extract the solution jc jm
toc
disp('... done!')
disp('-------------------------------------------------------------------')
%%
disp('-------------------------------------------------------------------')
disp('Post-processing M ...')
if 1
main_post_M; % compute the Magnetization vector 
else
disp('... off ...')  
disp('... type "main_post_M" to evaluate M...')
end
disp('... done!')
disp('-------------------------------------------------------------------')
%%
disp('-------------------------------------------------------------------')
disp('Plot JM ...')
main_plot_JM; % plot the results
disp('... done!')
disp('-------------------------------------------------------------------')
%%
return

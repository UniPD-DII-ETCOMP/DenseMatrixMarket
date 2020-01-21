Dense Matrix Market: PEEC 1D (data) by Riccardo Torchio (riccardo.torchio@unipd.it)

How to create a new user-defined test-case:

1. create a new subdirectory
2. the new subdirectory must contain two files:
   a) "description.txt": contains a description of the test case
   b) "data.mat": contains all data stuctures (see below)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Description of "data.mat"
(nStick:number of stick elements, nNodes:number of geometrical points, )

"data.mat" MUST contains:
G = 2 x nStick incidence matrix ( G(1,j)=end node of j-th stick element, G(2,j)=start node of j-th stick element (current flows from start to end))
NN = 3 x nNodes matrix containing coordinates of geometrical points
Voltage_Source = array of struct for voltage generators and appended elements
    Voltage_Source.node_start=[]; % index of starting node 
    Voltage_Source.node_end=[];   % index of ending node
    Voltage_Source.R=[];          % value of series resistance [Ohm]
    Voltage_Source.L=[];          % value of series self inductance [H]
    Voltage_Source.Cinv=[];       % value of series inverse of self capacitance [F^-1]
    Voltage_Source.value=[];      % value of voltage excitation [V]
Current_Source = array of struct for injected currents 
    Current_Source.node=[]; % index of the node where the current is injected
    Current_Source.value=[];% value of the injected current [A]

N.B.: All this variables MUST be contained in "data.mat". For instance,
      even if no Current Sources  are involved in the simulations, the varibale 
      "Current_Source" must be present with empty fields 
      (i.e. Current_Source.node=[];  Current_Source.value=[];)

INFO 
Voltage_Source.node_start
Voltage_Source.node_end
Current_Source.node
can be positive numbers in the range 1:nNodes (nNodes is the number of nodes of the mesh)
when the element is conneted to the device
can be negative numbers if the component is connected to some "appended" circuit node
can be qual to zero if the component is connected to the infinity node

example:
    % FIRST COMPONENT
    Voltage_Source(1).node_start=45; % connected to node 45 (i.e. NN(1:3,45) 
    Voltage_Source(1).node_end=-1; % connected to the extra appended node 
    Voltage_Source(1).R=50;         % resistance value [Ohm]
    Voltage_Source(1).L=0;         % self inductance value [H]
    Voltage_Source(1).Cinv=0;      % inverse self capacitance value [F^-1]
    Voltage_Source(1).value=3;     % voltage excitation value [V]
    % SECOND COMPONENT (short element)
    Voltage_Source(2).node_start=0; % connected to the infinity node  
    Voltage_Source(2).node_end=-1;  % connected to the extra appended node
    Voltage_Source(2).R=0;          
    Voltage_Source(2).L=0;
    Voltage_Source(2).Cinv=0;
    Voltage_Source(2).value=0;      % no external excitacion 
% similar for "Current_Source" variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subdirectories contain some examples of data generation:

"grid": geometrical data are taken from a Comsol mesh file

"other directories": the data are generated in Matlab or uploaded from 
existing .mat files

"line1","line2","line3": can be considered as a easy starting point for
the generation of new user-defined data




 
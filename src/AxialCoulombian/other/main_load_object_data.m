%% main load data
%% set 
%%  CONDUCTIVE
ex.vol_con=true; % true/false
if ex.vol_con % 
load data_okok_volu_con.mat    
end
%%  MAGNETIC
ex.vol_mag=true; % true/false 
if ex.vol_mag % 
load data_okok_volu_mag.mat    
end
%%  EXT COIL
ex.vol_ext=true; %  true/false
if ex.vol_ext % 
load data_okok_volu_ext.mat    
end
%% disp
disp('...Materials involved in the simulation ...')
disp(ex)
%%
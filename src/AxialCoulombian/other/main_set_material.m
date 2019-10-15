%% MATERIAL DATA
%% VOLUME CONDUCTIVE
rrho_con=1/2e3; % electric resistivity of conductive media [Ohm/m^2]
%
skin_depth=sqrt(2*rrho_con/(w*mu_0)); %skin_depth=sqrt(2*rrho_con/(w*abs(mmu_r)*mu_0))
disp(['SKIN_DEPTH = ',num2str(skin_depth),' [m]'])
%% VOLUME MAGNETIC 
mmu_r=5;% permeability of magnetic  media
rrho_mag=1/(1j*w(1)*mu_0*(mmu_r-1)); % magnetic resistivity of magnetic volume media
%%
% If different values of rrho_con and  mmu_r want be set go to 
% "main_set_material_vec" 
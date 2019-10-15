%%
sol.jphi=x(1:N.face_con); % j DoFs
sol.Jphi=sol.jphi./Area_c; % J
%%
sol.jm=x(N.face_con+1:N.face_con+N.edge_mag); % jm DoFs
sol.m=(1/(1j*w(1)*mu_0))*sol.jm; % m DoFs
sol.rho_m=mu_0*Matrix_C_m*sol.m; % rhom DoFs
sol.sig_m=mu_0*Matrix_Cb_m*sol.m; % sigm DoFs



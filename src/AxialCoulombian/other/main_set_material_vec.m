%% COND
rho_c=zeros(N.face_con,1);
for ii = 1:N.face_con
    rho_c(ii)=rrho_con;
end

%% MAG 
rho_m=zeros(N.face_mag,1);
mu_r_vec=zeros(N.face_mag,1);
for ii = 1:N.face_mag
    rho_m(ii)=rrho_mag;
    mu_r_vec(ii)=mmu_r;
end

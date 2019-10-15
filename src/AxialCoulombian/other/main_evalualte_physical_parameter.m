w = 2*pi*f; % omega
eps_0 = 8.85418781762e-12; % vacuum permittivity
mu_0 = 4*pi*1e-7; % vacuum permeability 
c_light=299792458;
lambda = c_light/f; % wave length
k0=w*sqrt(eps_0*mu_0);
Z_cost=sqrt(mu_0/eps_0);
beta = w/c_light; % propagation constant
%% disp
disp(['...lambda = ', num2str(lambda),'[m]']) 
disp(['...k0 = ', num2str(k0),'[1/m]']) 
function [I_obj,I_app,Q,U_obj,rho_charge,Jrenorm,Jimnorm,Jre,Jim] = fun_PEEC_1D_postprocessing(nSticks,x,N,...
                            Aobj,A_eapp,Is,omega,ret,NN,G,radiusvec,np_gauss,...
                            ne_cap_max,Cap_Elem,NN_xx1,NN_xx2,PPg_xx,ll_xx,...
                            ll_tot_xx,gauss_W,nNodes,beta,lgt,radius_cap)
%% extract solution
disp('------------------------------------')
disp('extract solution...')
I_obj = x(1:nSticks,1); %currents on stick elements [A]
I_app = x(nSticks+1:nSticks+N.vol_sou,1); %currents on appended branches [A]
Q = (Aobj*I_obj+A_eapp*I_app+Is)/(1j*omega); %charge on object nodes [C]
% warning: P matrix is here generated (U_obj = P*Q;)
if ret == 1
U_obj =(P_stick_assemble_matlab_vec_MF(NN,G,radiusvec,np_gauss,ne_cap_max,Cap_Elem,...
               NN_xx1,NN_xx2,PPg_xx,ll_xx,ll_tot_xx,...
               gauss_W,1:nNodes,1:nNodes)...
               .*exp(-1j*beta*dist_node_node_MF(NN,1:nNodes,1:nNodes)))*Q; % electric potential on object nodes [V]
else
U_obj =(P_stick_assemble_matlab_vec_MF(NN,G,radiusvec,np_gauss,ne_cap_max,Cap_Elem,...
               NN_xx1,NN_xx2,PPg_xx,ll_xx,ll_tot_xx,...
               gauss_W,1:nNodes,1:nNodes))*Q; % electric potential on object nodes [V]
    
end
U_app = x(nSticks+N.vol_sou+1:end,1); % potentials on appended nodes [V]
disp('done')
disp('------------------------------------')
%% plot of result
figure
subplot(2,1,1)
plot(real(U_obj),'r o-')
hold on
plot(imag(U_obj),'b o-')
title('Potential [V]')
subplot(2,1,2)
plot(real(I_obj),'r o-')
hold on
plot(imag(I_obj),'b o-')
title('Current [A]')
drawnow
%% post processing
disp('------------------------------------')
disp('post-processing')
%current density vector
disp('...computing J...')
tic
Jre = zeros(3,nSticks); % [Am^-2]
Jim = zeros(3,nSticks); % [Am^-2]
ut=zeros(3,nSticks);
Jrenorm=zeros(nSticks,1);
Jimnorm=zeros(nSticks,1);
for k=1:nSticks
    ut(:,k)= (NN(:,G(2,k))-NN(:,G(1,k)))/lgt(k);
    Jre(:,k) = real(I_obj(k))*ut(:,k)./(pi*radiusvec(k).^2);%./lgt(k);
    Jim(:,k) = imag(I_obj(k))*ut(:,k)./(pi*radiusvec(k).^2);%/lgt(k);
    Jrenorm(k)=norm(Jre(:,k));
    Jimnorm(k)=norm(Jim(:,k));
end
toc
disp('...computing charge_density...') 
rho_charge=Q./(pi*(ll_tot_xx.*(radius_cap.^2)));
disp('...done!')
disp('------------------------------------')  
end


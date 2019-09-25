%% ****************************** data_comsol *****************************
clear
close all
clc
%%
C = textread('data_com.mphtxt', '%s','delimiter', '\n');
%% read device geo data from comsol txt
[NN,G,nNodes,nSticks,Matrix_G] = fun_extract_from_comsol_1D(C);%main_extract;
%%
bar_stick = 0.5*(NN(:,G(1,:))+NN(:,G(2,:)));
%% initialize generators
count_Vgen=0;
Voltage_Source=[];
Voltage_Source.node_start=[];
Voltage_Source.node_end=[];
Voltage_Source.R=[];
Voltage_Source.L=[];
Voltage_Source.Cinv=[];
Voltage_Source.value=[];
count_Igen=0;
Current_Source=[];
Current_Source.node=[];
Current_Source.value=[];
%% *************************  VOLTAGE generators and Lumped branches ******
disp('setting generators...')
tic
figure(1)
%connected to OBJ nodes
if 1
for k=1:nNodes
    %gen pos   
    if NN(1,k)<0.01 && NN(2,k)<0.01
        count_Vgen=count_Vgen+1;
        Voltage_Source(count_Vgen).node_start=0;
        Voltage_Source(count_Vgen).node_end=k;
        Voltage_Source(count_Vgen).R=0;
        Voltage_Source(count_Vgen).L=0;
        Voltage_Source(count_Vgen).Cinv=0;
        Voltage_Source(count_Vgen).value=1;
        plot3(NN(1,k),NN(2,k),NN(3,k),'*r')
    elseif NN(1,k)>(1-0.01) && NN(2,k)>(1-0.01)
        count_Vgen=count_Vgen+1;
        Voltage_Source(count_Vgen).node_start=0;
        Voltage_Source(count_Vgen).node_end=k;
        Voltage_Source(count_Vgen).R=0;
        Voltage_Source(count_Vgen).L=0;
        Voltage_Source(count_Vgen).Cinv=0;
        Voltage_Source(count_Vgen).value=-1;
        plot3(NN(1,k),NN(2,k),NN(3,k),'*r')        
    end
end
end
n_Vsource=count_Vgen
%% ***********************  CURRENT generators  **********************
if 0
for k=1:nNodes
    if NN(1,k)<0.01 && NN(2,k)<0.01
        count_Igen=count_Igen+1;
        Current_Source(count_Igen).node=k;
        Current_Source(count_Igen).value=1;
        plot3(NN(1,k),NN(2,k),NN(3,k),'*g')   
    end
    if NN(1,k)>(1-0.01) && NN(2,k)>(1-0.01)
        count_Igen=count_Igen+1;
        Current_Source(count_Igen).node=k;
        Current_Source(count_Igen).value=-1;
        plot3(NN(1,k),NN(2,k),NN(3,k),'*g')   
    end    
end
end
n_Isource=count_Igen
%%
legend('Voltage Source','Current Source')
toc
%%
app=[[Voltage_Source.node_start],[Voltage_Source.node_end],[Current_Source.node]];
app=unique(app);
app=find(app<0);
N.node_app=length(app);% number of appended nodes
%%
disp('done')
disp('---------------------------------')
disp(['SYSTEM SIZE: ' num2str(nSticks+n_Vsource+N.node_app) ' x ' num2str(nSticks+n_Vsource+N.node_app)])
disp('---------------------------------')
%% 
save data.mat G NN Voltage_Source Current_Source
%%
disp('data saved')


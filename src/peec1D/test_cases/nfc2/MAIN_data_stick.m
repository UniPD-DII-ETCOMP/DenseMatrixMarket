clear
close all
clc
%% load NFC antenna data
load data_NFC.mat
%%
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
%%
k1=1;
k2=176;
figure
plot3(0,0,0,'r')
hold on
plot3(0,0,0,'g')
plot3(0,0,0,'b')
hold on
axis equal 
view(3)
plot3([NN(1,G(1,:));NN(1,G(2,:))],[NN(2,G(1,:));NN(2,G(2,:))],[NN(3,G(1,:));NN(3,G(2,:))],'.-','color','k')
%% refine
if 1
Nref=3;
disp('refine...')
for ii = 1:Nref
disp(['...',num2str(ii),'...'])    
[G,NN] = fun_refine(NN,G);
end
disp('...done!')
close all
figure
plot3(0,0,0,'r')
hold on
plot3(0,0,0,'g')
plot3(0,0,0,'b')
hold on
axis equal 
view(3)
plot3([NN(1,G(1,:));NN(1,G(2,:))],[NN(2,G(1,:));NN(2,G(2,:))],[NN(3,G(1,:));NN(3,G(2,:))],'.-','color','k')
end
%% *************************  VOLTAGE generators **********************
disp('setting generators...')
tic
figure(1)
if 1
count_Vgen=1;
Voltage_Source(1).node_start=k1;
Voltage_Source(1).node_end=k2;
Voltage_Source(1).R=0;
Voltage_Source(1).L=0;
Voltage_Source(1).Cinv=0;
Voltage_Source(1).value=1;
plot3(NN(1,k1),NN(2,k1),NN(3,k1),'*r')
plot3(NN(1,k2),NN(2,k2),NN(3,k2),'*r')
end
%% ***********************  CURRENT generators  **********************
%connected to OBJ nodes 
if 0
count_Igen=1;
Current_Source(count_Igen).node=k1;
Current_Source(count_Igen).value=1-1j;
count_Igen=2;
Current_Source(count_Igen).node=k2;
Current_Source(count_Igen).value=-(1-1j);
plot3(NN(1,k1),NN(2,k1),NN(3,k1),'*g')   
plot3(NN(1,k2),NN(2,k2),NN(3,k2),'*g') 
end
%%
disp('...done!')
legend('Voltage Source','Current Source')
toc
%% 
save data.mat G NN Voltage_Source Current_Source
%%


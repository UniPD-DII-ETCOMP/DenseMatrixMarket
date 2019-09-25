clear
close all
clc
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
%% two wire line data
long=1; % length of the two_wire_line [m]
dist=0.1; % distance between the two lines [m]
%%
nsp=100; % number of points
NN1=zeros(3,nsp);
NN1(1,:)=linspace(0,long,nsp);
NN1(3,:)=dist*0.5;
NN2=zeros(3,nsp);
NN2(1,:)=linspace(0,long,nsp);
NN2(3,:)=-dist*0.5;
NN=[NN1,NN2];
G=[1:nsp-1;2:nsp];
G=[G,G+nsp];
nNodes=size(NN,2);
nSticks=size(G,2);
k1=1;
k2=nsp+1;
k3=nsp;
k4=nsp*2;
%%
bar_stick = 0.5*(NN(:,G(1,:))+NN(:,G(2,:)));
%%
figure(1)
plot3(0,0,0,'r')
hold on
plot3(0,0,0,'g')
plot3(0,0,0,'b')
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
plot3([NN(1,G(1,:));NN(1,G(2,:))],[NN(2,G(1,:));NN(2,G(2,:))],[NN(3,G(1,:));NN(3,G(2,:))],'.-','color','k')
%% *************************  VOLTAGE generators and Lumped branches ******
disp('setting generators...')
tic
figure(1)
if 1
    count_Vgen=count_Vgen+1;
    Voltage_Source(count_Vgen).node_start=k3; 
    Voltage_Source(count_Vgen).node_end=k4; 
    Voltage_Source(count_Vgen).R=0;
    Voltage_Source(count_Vgen).L=0;
    Voltage_Source(count_Vgen).Cinv=0;
    Voltage_Source(count_Vgen).value=0;
    plot3(NN(1,k3),NN(2,k3),NN(3,k3),'*r')  
    plot3(NN(1,k4),NN(2,k4),NN(3,k4),'*r')  
end
n_Vsource=count_Vgen
%% ***********************  CURRENT generators  **********************
%connected to OBJ nodes 
count_Igen=count_Igen+1;
Current_Source(count_Igen).node=k1;
Current_Source(count_Igen).value=1-1j;
plot3(NN(1,k1),NN(2,k1),NN(3,k1),'*g')  
count_Igen=count_Igen+1;
Current_Source(count_Igen).node=k2;
Current_Source(count_Igen).value=-(1-1j);
plot3(NN(1,k2),NN(2,k2),NN(3,k2),'*g')  
%connected to APPENDED nodes and INF  (!! appended node index is negative)
if 0
I_napp=[1 -1];
for k=1:length(I_napp)
    count_Igen=count_Igen+1;
    Current_Source(count_Igen).node=-k;
    Current_Source(count_Igen).value=I_napp(k);
end
end
n_Isource=count_Igen
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
save data.mat G NN Voltage_Source  Current_Source
%%
disp('data saved')


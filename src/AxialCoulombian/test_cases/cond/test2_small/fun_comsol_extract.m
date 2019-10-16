function [Matrix_P0,VP] = fun_comsol_extract(C)
for ii = 1:size(C,1)
    CC = char(C(ii));
   if size(CC) == size('# Mesh point coordinates')
   if CC ==  '# Mesh point coordinates'
       ind_start_point = ii+1;
   end
   end
   if size(CC) == size('4 # number of nodes per element')
   if CC ==  '4 # number of nodes per element'
       ind_start_q = ii+3;
   end
   end
end
Npoint= C(19)
Ntetra = C(ind_start_q-2)
CC=char(C(19));
ind=find(CC=='#');
Np=str2double(CC(1:ind-1))
CC=char(C(ind_start_q-2));
ind=find(CC=='#');
Nv=str2double(CC(1:ind-1))
a = str2num(char(C(ind_start_point:ind_start_point+Np)));
Matrix_P0(:,1)=a(:,1);
Matrix_P0(:,3)=a(:,2);
VP = str2num(char(C(ind_start_q:ind_start_q+Nv))); 
VP = VP+1;
check = max(max(VP)) == Np
VP = VP.';
reorder = 1;
permute_Mat=[2 1 3 4];
if reorder == 1 
    for ii=1:size(VP,2)
        VP(:,ii)=VP(permute_Mat,ii);
    end
end
end
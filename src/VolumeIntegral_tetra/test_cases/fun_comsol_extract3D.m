function [VP,Matrix_P0] = fun_comsol_extract3D(C)
for ii = 1:size(C,1)
    CC = char(C(ii));
   if size(CC) == size('# Mesh point coordinates')
   if CC ==  '# Mesh point coordinates'
       ind_start_point = ii+1;
   end
   end
   if size(CC) == size('4 # number of nodes per element')
   if CC ==  '4 # number of nodes per element'
       ind_start_tetra = ii+3;
   end
   end
end
Npoint= C(18);
Ntetra = C(ind_start_tetra-2);
CC=char(C(18));
ind=find(CC=='#');
Np=str2double(CC(1:ind-1));
CC=char(C(ind_start_tetra-2));
ind=find(CC=='#');
Nv=str2double(CC(1:ind-1))
Matrix_P0 = str2num(char(C(ind_start_point:ind_start_point+Np)));
VP = str2num(char(C(ind_start_tetra:ind_start_tetra+Nv))); 
VP = VP+1;
check = max(max(VP)) == Np;
VP = VP.';
reorder = 1;
%% rerodering
if reorder == 1 
a = zeros(size(VP,2),1);
for ii = 1:size(VP,2)
    e12 = Matrix_P0(VP(2,ii),:)-Matrix_P0(VP(1,ii),:);
    e13 = Matrix_P0(VP(3,ii),:)-Matrix_P0(VP(1,ii),:);
    e14 = Matrix_P0(VP(4,ii),:)-Matrix_P0(VP(1,ii),:);
    vec1=cross(e12,e13);
    a(ii,1) = dot(vec1,e14);
    if a(ii,1) < 0
        VP(:,ii) = VP([1,2,4,3],ii);
    end
end
end

end
 
 
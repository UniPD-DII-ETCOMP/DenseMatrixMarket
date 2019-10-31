function [Vertices Faces FaceID] = cubedsphere(n, prjtype)
% function [Vertices Faces FaceID] = cubedsphere(n, prjtype)
%
% Generate a cartesian coordinates of the cubed-sphere of radius 1
% the grid is gnomonic equiangular/equidistance central projection
%
% INPUT:
% n: number of linear subdivision of the face
% prjtype: optional, 'equiangular' (default) or 'equidistance'
% OUTPUT:
% Vertices: (nv x 3) array, where nv = n^2*6+2, 3D vertices coordinates
% Faces: (nf x 4) array, where nf = n^2*6, vertices indices of patches
% FaceID: (nf x 1) array, number where the faces belong (to the cube
% topology, 1-6)
% +---+
% | 6 |
% +---+---+---+---+
% | 4 | 1 | 2 | 3 |
% +---+---+---+---+
% | 5 |
% +---+
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
% 14-Aug-2009 original
% 15-Aug-2009, change the patch order (ndgrid -> meshgrid, code #74)
n = round(n);
if n<1
 error('cubedsphere: n must be stricly positive number');
end
if nargin<2 || isempty(prjtype)
 prjtype = 'equiangular';
end
n = n+1;
% Discretize basic face of a cube
switch lower(prjtype)
 case 'equiangular'
 % equidistance projection
 x = 1;
 theta = linspace(-pi/4, pi/4, n);
 y = tan(theta);
 z = y;
 [X Y Z] = ndgrid(x,y,z);
 case 'equidistance'
 % equidistance projection
 x = 1;
 y = linspace(-1,1,n);
 z = y;
 [X Y Z] = ndgrid(x,y,z);
 otherwise
 error('cubedsphere: prjtype must be ''equiangular'' or ''equidistance''');
end
% Project on sphere S2
S = sqrt(1./(X.^2+Y.^2+Z.^2));
X = X.*S;
Y = Y.*S;
Z = Z.*S;
% Generate six faces
C1 = [X(:) Y(:) Z(:)].';
% Other faces are rotations of the first
M = makehgtform('zrotate',pi/2);
C2 = M(1:3,1:3)*C1;
M = makehgtform('zrotate',pi);
C3 = M(1:3,1:3)*C1;
M = makehgtform('zrotate',-pi/2);
C4 = M(1:3,1:3)*C1;
M = makehgtform('yrotate',pi/2);
C5 = M(1:3,1:3)*C1;
M = makehgtform('yrotate',-pi/2);
C6 = M(1:3,1:3)*C1;
% Group the faces totheger
C = [C1 C2 C3 C4 C5 C6].';
clear C1 C2 C3 C4 C5 C6; % clean up
% Patch topology of the first face
[i j] = meshgrid(1:n-1,1:n-1); % 090815: Change to meshgrid from ndgrid
i = i(:); j=j(:);
P = [sub2ind([n n], i, j) ...
 sub2ind([n n], i+1, j) ...
 sub2ind([n n], i+1, j+1) ...
 sub2ind([n n], i, j+1)];
% Topology for all faces
P = reshape(P,[size(P,1) 1 size(P,2)]);
offset = (n^2)*(0:5);
P = bsxfun(@plus, offset, P);
P = reshape(P, [], 4);
FaceID = repmat(1:6,(n-1)^2,1);
FaceID = FaceID(:);
% Detect common vertices (between two neighboring faces)
tol = 1e-3/(n-1);
Cround = round(C/tol);
[trash I J] = unique(Cround,'rows'); %#ok mlint complains on trash
% Remove vertices and update face topology
Vertices = C(I,:);
Faces = J(P);
end % cubedsphere

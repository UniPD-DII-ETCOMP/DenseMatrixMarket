function [ col2 ] = tetra_mesh2(VP,P0,alpha,col)


patch('Faces',VP(1:3,:).','Vertices',P0,'Facecolor',col,'FaceAlpha',alpha)
hold on
patch('Faces',VP([1,2,4],:).','Vertices',P0,'Facecolor',col,'FaceAlpha',alpha)
patch('Faces',VP([1,3,4],:).','Vertices',P0,'Facecolor',col,'FaceAlpha',alpha)
patch('Faces',VP([2,3,4],:).','Vertices',P0,'Facecolor',col,'FaceAlpha',alpha)

col2=col;

% patch('Faces',5:8,'Vertices',Hexa_point,'Facecolor',col,'FaceAlpha',alpha)
% patch('Faces',[1,2,6,5],'Vertices',Hexa_point,'Facecolor',col,'FaceAlpha',alpha)
% patch('Faces',[4,3,7,8],'Vertices',Hexa_point,'Facecolor',col,'FaceAlpha',alpha)
% patch('Faces',[1,4,8,5],'Vertices',Hexa_point,'Facecolor',col,'FaceAlpha',alpha)
% patch('Faces',[2,3,7,6],'Vertices',Hexa_point,'Facecolor',col,'FaceAlpha',alpha)


end


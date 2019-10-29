function [map,invmap] = fun_reo_DMM(bared,N,ind,dims,base,lev,...
    xmin,xmax,ymin,ymax,zmin,zmax)
bar=bared(:,ind.cotree);
clust.n = N.cotree; 
clust.d=dims; 
clust.x=bar(1:clust.d,:).'; fflag=zeros(N.cotree,1).';
[map,invmap,leaf,clust]=fun_reo_RT(clust,lev,1,0,...
    base,fflag,length(unique(fflag))); map=map+1; 
figure
subplot(1,2,1)
scatter3(bar(1,:),bar(2,:),bar(3,:),5,'filled','CData',1:N.cotree)
axis equal
colormap jet
title('DoF index before reordering')
axis([xmin,xmax,ymin,ymax,zmin,zmax])
subplot(1,2,2)
scatter3(bar(1,map),bar(2,map),bar(3,map),5,'filled','CData',1:N.cotree)
axis equal
title('DoF index after reordering')
axis([xmin,xmax,ymin,ymax,zmin,zmax])
drawnow
end


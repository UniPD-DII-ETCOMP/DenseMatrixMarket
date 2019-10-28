function [N,Q,ind] = fun_tree_cotree(ind,G1,Matrix_C,N)
tic
%% DIE
reo=[ind.edge_ext,ind.edge_int];
tree=quicktree(G1(:,reo).');
tree=reo(tree);
N.dom=N.node-length(tree);
%%
ind.tree=[ind.edge_ext,tree];
ind.cotree=setdiff(1:N.edge,ind.tree);
Q=Matrix_C(:,ind.cotree);    
N.tree=length(ind.tree);
N.cotree=length(ind.cotree);    
%%
toc
end
%%

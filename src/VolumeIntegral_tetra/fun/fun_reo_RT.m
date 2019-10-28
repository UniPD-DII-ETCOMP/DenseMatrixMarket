function [index,inv_index,leaf,cluster]=fun_reo_RT(clust_fact,n_reo_lev,split_dir,id_start,base,flag,N_flag)
index=linspace(0,clust_fact.n-1,clust_fact.n);
%%
N_flags=zeros(N_flag,1);
for ii = 1:N_flag
   N_flags(ii)=sum((flag==(ii-1)));
end
%% type of reordering
[index,cluster]=reg_partition_mod(index,n_reo_lev,split_dir,clust_fact.n,id_start,clust_fact,base,flag,N_flag,N_flags);
%% inverse mapping (from original to reordered mesh)
[~,inv_index]=sort(index);
%% leaf clusters
n_leafs=N_flag*base^(n_reo_lev);
leaf=zeros(2,n_leafs);
st=size(cluster,2)-N_flag*base^(n_reo_lev-1);
for ii = 1:n_leafs/base
    for kk = 1:base
        leaf(1,(ii-1)*base+N_flag+kk)=[cluster(st+ii).start(kk)];
        leaf(2,(ii-1)*base+N_flag+kk)=[cluster(st+ii).size(kk)];
    end
end
end
function [index,cluster]=reg_partition_mod(index,n_lev,split_dir,sz,start,clust_fact,base,flag,N_flag,N_flags)
%initiate tree from root cluster
[index,cluster]=ind_bisection_flag_part(index,sz,start,clust_fact,N_flag,flag,N_flags);
count=1;
%splitting sons at level p
daddy=1;
for p=1:n_lev
%     daddy=base^(p-1):base^p-1;
    if p == 1
        base_loc=N_flag;
    else
        base_loc=base;
    end
    for k=1:length(daddy)
        for ss=1:base_loc
            [index,sson]=ind_bisection(index,cluster(daddy(k)).size(ss),cluster(daddy(k)).start(ss),clust_fact,split_dir,base);
            count=count+1;
            cluster(count)=sson;
        end
    end
    daddy=daddy(end)+[1:N_flag*base^(p-1)];
    split_dir=mod(split_dir,clust_fact.d)+1;
end
end
%%
function [index,cluster]=ind_bisection_flag_part(index,n,start,clust_fact,base,flag,N_flags)
d=clust_fact.d;
x=clust_fact.x;
%% build the bouding box
%preliminary reordering of incides
for ii = 1:base
    sons(ii).son=index(flag==(ii-1));
end
for ii = 1:base
    index(start+sum(N_flags(1:ii-1))+1:start+sum(N_flags(1:ii)))=sons(ii).son;
end
for ii = 1:base
    size_sons(ii).size=length(sons(ii).son);
    length_sons(1,ii)=length(sons(ii).son);
end
cluster.ns=base; % numero di figli?
cluster.size=length_sons;%[size_l,size_r];  % dimensione dei figli 
cluster.start(1)=start+0;
for ii = 2:base
    cluster.start(ii)=start+sum(length_sons(1,1:ii-1)); % indice di partenza dei figli (da 0) 
end
end
%%
function [index,cluster]=ind_bisection(index,n,start,clust_fact,split_dir,base)
d=clust_fact.d;
x=clust_fact.x;
%% build the bouding box
%preliminary reordering of incides
[~,iloc]=sort(x(index(start+1:start+n)+1,split_dir));
mid=floor(n/base);
for ii = 1:base-1
    sons(ii).son=index(start+iloc(mid*(ii-1)+1:ii*mid));
end
ii=base;
sons(ii).son=index(start+iloc(mid*(ii-1)+1:end));
for ii = 1:base-1
    index(start+(mid*(ii-1)+1):start+(mid*ii))=sons(ii).son;
end
ii=base;
index(start+(mid*(ii-1)+1):start+n)=sons(ii).son;
for ii = 1:base
    size_sons(ii).size=length(sons(ii).son);
    length_sons(1,ii)=length(sons(ii).son);
end
cluster.ns=base; % 
cluster.size=length_sons;
cluster.start(1)=start+0;
for ii = 2:base
    cluster.start(ii)=start+sum(length_sons(1,1:ii-1));
end
end
%%
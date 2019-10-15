function [g,c]=gcd_quad(fp,nf)
nn=max(max(fp(1:4,1:nf)));
esiz=max(nn*10,10000);
ehtpsiz=round((log2(esiz))+1);
elocsiz=2^ehtpsiz;
g=zeros(2,elocsiz);
c=zeros(4,nf);
%%
locedg_face =[1,2,3,4;...
              2,3,4,1].';
%G matrix
ne=0;
for i=1:nf  
  for j=1:4 
     nnd2(:)=fp(locedg_face(j,:),i); 
     nnd2ori(:)=nnd2(:);
     nnd2=sort(nnd2);
     ind=htsrcedge(nnd2,g);
     if ind==0
       ne=ne+1;
       if(ne>elocsiz)
         disp(['ne > locsiz in gcd, aborting'])
       end
       g(:,ne)=nnd2;   
     end
  end
end
g=g(:,1:ne);
%C matrix
c=zeros(4,nf);
for i=1:nf
  for j=1:4
     nnd2(1:2)=fp(locedg_face(j,1:2),i);
     nnd2ori(:)=nnd2(:);
     nnd2=ordsw2(nnd2(1),nnd2(2));
     ind=htsrcedge(nnd2,g);
     ied=ind;
     if(ind==0)
       disp('Immpossible error XXX in gcd')
       return
     end
     if(nnd2ori(1)==nnd2(1) && nnd2ori(2)==nnd2(2))
       c(j,i)=ied;
     else
       c(j,i)=-ied;
     end
   end
end

end 
%%
function [res]=ordsw2(aa,bb)
  if(aa>=bb)
    kk1=aa;
    aa=bb;
    bb=kk1;
  end
  res=[aa,bb];
end 
%%
function ind=htsrcedge ( ed, g)
who=find(g(1,:)==ed(1));
who2= g(2,who)==ed(2);
ind=who(who2);
if isempty(ind)
    ind=0;
end
end
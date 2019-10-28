function [gshort,c,d,f] = gcd_matlab(vp)
locedg=[1,2;2,3;3,1;4,1;4,2;4,3];
locfacedg=[2,3;3,1;1,2];
locfac=[2,1,1,1;3,4,2,3;4,3,4,2].';
nn=max(max(vp));
nv=size(vp,2);
esiz=max(nn*10,10000);
ehtpsiz=ceil(log2(esiz))+1;
elocsiz=2^ehtpsiz;
eht(1:elocsiz)=0;
elinks(1:elocsiz)=0;
g(1:2,1:elocsiz)=0;
%G matrix
ne=0;
for i=1:nv  
  for j=1:6
     nnd2(:)=vp(locedg(j,:),i);
     nnd2ori(:)=nnd2(:);
     nnd2 =ordsw2(nnd2(1),nnd2(2));
     ind=htsrcedge(nnd2,g);
     if ind==0
       ne=ne+1;
       if(ne>elocsiz)
         error(['ne > locsiz in gcd, aborting'])
       end
       g(:,ne)=nnd2;
     end
 end
end
%F and D matrices
fsiz=max(nv*3,10000);
fhtpsiz=ceil(log2(fsiz))+1;
flocsiz=2^fhtpsiz;
f(1:5,1:flocsiz)=0;
ford(1:3,1:flocsiz)=0;
d(1:4,1:nv)=0;
nf=0;
for i=1:nv
  for j=1:4
    nnd3(:)=vp(locfac(j,:),i);
    nnd3rot =nnd3;
    nnd3 = sort(nnd3);
    [nnd3rot(1),nnd3rot(2),nnd3rot(3)]=rotate3(nnd3rot(1),nnd3rot(2),nnd3rot(3));
    ind=htsrcface(nnd3,ford);
    if ind==0
      nf=nf+1;
      if nf > flocsiz
        error(['nf > flocsi in gcd, aborting'])
      end
      f(1:3,nf)=nnd3rot(1:3);
      ford(1:3,nf)=nnd3(1:3);
	  f(4,nf)=vp(j,i);
	  f(5,nf)=0;
      %regions(1,nf)=reg(i)
      d(j,i)=nf;
    else
      %regions(2,ind)=reg(i)
      d(j,i)=-ind;
	  f(5,ind)=vp(j,i);
    end
  end
end

%C matrix
c(1:3,1:nf)=0;
for i=1:nf
  for j=1:3
     nnd2(1:2)=f(locfacedg(j,1:2),i);
     nnd2ori(:)=nnd2(:);
     nnd2=sort(nnd2);
     ind=htsrcedge(nnd2,g);
     ied=ind;
     if ind == 0 
        error('impossible error')
     end
     if nnd2ori(1)>=nnd2(1) && nnd2ori(2)>=nnd2(2)
       c(j,i)=ied;
     else
       c(j,i)=-ied;
     end
  end
end
gshort(1:2,1:ne)=g(1:2,1:ne);
f=f(:,1:nf);
c=c(:,1:nf);
end 

function [aa,bb,cc] = rotate3(aa,bb,cc)
imin=min([aa,bb,cc]);
if(aa == imin) 
    return
end
if(bb == imin)
  kk1=aa;
  kk2=bb;
  kk3=cc;
  aa=kk2;
  bb=kk3;
  cc=kk1;
  return
end
  kk1=aa;
  kk2=bb;
  kk3=cc;
  aa=kk3;
  bb=kk1;
  cc=kk2;
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
who1=find(g(1,:)==ed(1));
who2= g(2,who1)==ed(2);
ind=who1(who2);
if isempty(ind)
    ind=0;
end
end

function ind=htsrcface(abc,f)
who1=find(f(1,:)==abc(1));
who2= f(2,who1)==abc(2);    
who3= f(3,who1)==abc(3); 
who2=find(who2.*who3);
ind=who1(who2);
if isempty(ind)
    ind=0;
% else
%    ah=1; 
end
end
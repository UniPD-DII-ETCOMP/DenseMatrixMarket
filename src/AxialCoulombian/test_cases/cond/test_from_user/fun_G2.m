function G2=fun_G2(F1,C1,N,G1)
%%
loc_face_edge(1 ,1:2)=[1,2];
loc_face_edge(2 ,1:2)=[2,3];
loc_face_edge(3 ,1:2)=[3,4];
loc_face_edge(4 ,1:2)=[4,1];
%%
G2=zeros(8,N.edge); % 
cont=-ones(1,N.edge); 
for ii = 1:N.face
    p_face=F1(1:4,ii); % 
    e_face=C1(:,ii); % 
    e_face_a=abs(e_face);
    cont(e_face_a)=cont(e_face_a)+1; % 
    %
    p_e1=G1(:,e_face_a(1));
    p_e2=G1(:,e_face_a(2));
    p_e3=G1(:,e_face_a(3));
    p_e4=G1(:,e_face_a(4));  
    % ed 1
    id_p=setdiff(p_face,p_e1);% 
    if sort(id_p)==sort(p_e2) %
     ed_o=2; % 
     p_ed_o=p_e2; %
    end
    if sort(id_p)==sort(p_e3) % 
     ed_o=3; % 
     p_ed_o=p_e3; %
    end  
    if sort(id_p)==sort(p_e4) % 
     ed_o=4; % 
     p_ed_o=p_e4; % 
    end     
    if cont(e_face_a(1))==0 % 
        G2(1:2,e_face_a(1))=p_e1; % 
        if sign(e_face(ed_o))*sign(e_face(1))==1 % 
            G2(3:4,e_face_a(1))=p_ed_o;
        elseif sign(e_face(ed_o))*sign(e_face(1))==-1 % 
            G2(3:4,e_face_a(1))=flip(p_ed_o);
        end
    elseif cont(e_face_a(1))==1 % 
        G2(5:6,e_face_a(1))=flip(p_e1); %
        if sign(e_face(ed_o))*sign(e_face(1))==1 %
            G2(7:8,e_face_a(1))=flip(p_ed_o);
        elseif sign(e_face(ed_o))*sign(e_face(1))==-1 % 
            G2(7:8,e_face_a(1))=p_ed_o;
        end
    end 
% ed 2
    id_p=setdiff(p_face,p_e2);% 
    if sort(id_p)==sort(p_e1) % 
     ed_o=1; % 
     p_ed_o=p_e1; % 
    end
    if sort(id_p)==sort(p_e3) % 
     ed_o=3; % 
     p_ed_o=p_e3; % 
    end  
    if sort(id_p)==sort(p_e4) %
     ed_o=4; % 
     p_ed_o=p_e4; % 
    end     
    if cont(e_face_a(2))==0 % 
        G2(1:2,e_face_a(2))=p_e2; % 
        if sign(e_face(ed_o))*sign(e_face(2))==1 % 
            G2(3:4,e_face_a(2))=p_ed_o;
        elseif sign(e_face(ed_o))*sign(e_face(2))==-1 %
            G2(3:4,e_face_a(2))=flip(p_ed_o);
        end
    elseif cont(e_face_a(2))==1 %
        G2(5:6,e_face_a(2))=flip(p_e2); % 
        if sign(e_face(ed_o))*sign(e_face(2))==1 % 
            G2(7:8,e_face_a(2))=flip(p_ed_o);
        elseif sign(e_face(ed_o))*sign(e_face(2))==-1 % 
            G2(7:8,e_face_a(2))=p_ed_o;
        end
    end       
% ed 3
    id_p=setdiff(p_face,p_e3);% 
    if sort(id_p)==sort(p_e2) % 
     ed_o=2; %
     p_ed_o=p_e2; %
    end
    if sort(id_p)==sort(p_e1) % 
     ed_o=1; % 
     p_ed_o=p_e1; % 
    end  
    if sort(id_p)==sort(p_e4) %
     ed_o=4; %
     p_ed_o=p_e4; %
    end     
    if cont(e_face_a(3))==0 % 
        G2(1:2,e_face_a(3))=p_e3; % 
        if sign(e_face(ed_o))*sign(e_face(3))==1 % 
            G2(3:4,e_face_a(3))=p_ed_o;
        elseif sign(e_face(ed_o))*sign(e_face(3))==-1 % 
            G2(3:4,e_face_a(3))=flip(p_ed_o);
        end
    elseif cont(e_face_a(3))==1 % 
        G2(5:6,e_face_a(3))=flip(p_e3); % 
        if sign(e_face(ed_o))*sign(e_face(3))==1 % 
            G2(7:8,e_face_a(3))=flip(p_ed_o);
        elseif sign(e_face(ed_o))*sign(e_face(3))==-1 % 
            G2(7:8,e_face_a(3))=p_ed_o;
        end
    end      
% lato 4
    id_p=setdiff(p_face,p_e4);% 
    if sort(id_p)==sort(p_e2) % 
     ed_o=2; % 
     p_ed_o=p_e2; %
    end
    if sort(id_p)==sort(p_e3) % 
     ed_o=3; %
     p_ed_o=p_e3; % 
    end  
    if sort(id_p)==sort(p_e1) % 
     ed_o=1; %
     p_ed_o=p_e1; % 
    end     
    if cont(e_face_a(4))==0 % 
        G2(1:2,e_face_a(4))=p_e4; % 
        if sign(e_face(ed_o))*sign(e_face(4))==1 % 
            G2(3:4,e_face_a(4))=p_ed_o;
        elseif sign(e_face(ed_o))*sign(e_face(4))==-1 % 
            G2(3:4,e_face_a(4))=flip(p_ed_o);
        end
    elseif cont(e_face_a(4))==1 % 
        G2(5:6,e_face_a(4))=flip(p_e4); % 
        if sign(e_face(ed_o))*sign(e_face(4))==1 % 
            G2(7:8,e_face_a(4))=flip(p_ed_o);
        elseif sign(e_face(ed_o))*sign(e_face(4))==-1 % 
            G2(7:8,e_face_a(4))=p_ed_o;
        end
    end    
end
for ii = 1:N.face
    p_face=F1(1:4,ii); % 
    e_face=C1(:,ii); %
    e_face_a=abs(e_face);
    for jj = 1:4 % 
        a=setdiff(p_face,G2(1:4,e_face_a(jj)));
        a=size(a,1)*size(a,2);
        if sign(e_face(jj))<0 && a>0 % 
            fl_1_4=G2(1:4,e_face_a(jj));
            fl_5_8=G2(5:8,e_face_a(jj));
            G2(1:4,e_face_a(jj))=fl_5_8;
            G2(5:8,e_face_a(jj))=fl_1_4;
        end     
    end
end  
Face_e_loc_mag=zeros(4,N.face);% Face_e_loc_mag
for ii = 1:N.face % 
    pf=F1(1:4,ii);
    e1=pf(loc_face_edge(1,:)); 
    e2=pf(loc_face_edge(2,:)); 
    e3=pf(loc_face_edge(3,:)); 
    e4=pf(loc_face_edge(4,:)); 
    for jj = 1:4 % 
       edge_a=abs(C1(jj,ii)); % 
       pejj=G1(:,edge_a);% 
       if sort(pejj)==sort(e1)
          Face_e_loc_mag(jj,ii)=1; 
       elseif sort(pejj)==sort(e2)
          Face_e_loc_mag(jj,ii)=2; 
       elseif sort(pejj)==sort(e3)
          Face_e_loc_mag(jj,ii)=3; 
       elseif sort(pejj)==sort(e4)
          Face_e_loc_mag(jj,ii)=4;           
       end
    end
end
%
end
function [L] = fun_L_interface(N_node_con,N_face_con,Area_c,Matrix_P0_c,F1_c,...
                              N_thread,ii,jj)

if size(ii,1)~=1
    ii=ii.';
end
if size(jj,1)~=1
    jj=jj.';
end
L=funLphiphi3_for_st(N_node_con,N_face_con,Area_c.',...
            1,2,Matrix_P0_c,F1_c,N_thread,length(ii),length(jj),ii,jj);
end


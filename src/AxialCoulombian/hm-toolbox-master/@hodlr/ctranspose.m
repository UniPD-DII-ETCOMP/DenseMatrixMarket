function Ht = ctranspose(H)
%TRANSPOSE Conjugate transpose of the H-matrix H.

Ht = H;

if is_leafnode(H)
    Ht.F = Ht.F';
    Ht.sz = H.sz(2:-1:1);
else
    Ht.A11 = H.A11';
    Ht.A22 = H.A22';
    Ht.U12 = H.V21;
    Ht.V12 = H.U21;
    Ht.U21 = H.V12;
    Ht.V21 = H.U12;
    Ht.sz = H.sz(2:-1:1);
end

end

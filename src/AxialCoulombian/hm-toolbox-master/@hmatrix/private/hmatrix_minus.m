function H = hmatrix_minus(H1, H2)

if any(H1.sz ~= H2.sz)
    error('HMATRIX_MINUS:: non conformal partitioning')
end

H = hmatrix;
H.sz = H1.sz;

if is_leafnode(H1)
    if is_leafnode(H2)
        if H1.admissible && H2.admissible
            H.U = [H1.U, -H2.U];
            H.V = [H1.V, H2.V];
            H.admissible = true;
        else
            H.F = full(H1) - full(H2);
        end
    else
        if H1.admissible
            H = hmatrix_rank_update(-H2, H1.U, H1.V, []);
        else
            H.F = full(H1) - full(H2);
        end
    end
elseif is_leafnode(H2)
    H = -hmatrix_minus(H2, H1);
else
    H.A11 = hmatrix_minus(H1.A11, H2.A11);
    H.A12 = hmatrix_minus(H1.A12, H2.A12);
    H.A21 = hmatrix_minus(H1.A21, H2.A21);
    H.A22 = hmatrix_minus(H1.A22, H2.A22);
end

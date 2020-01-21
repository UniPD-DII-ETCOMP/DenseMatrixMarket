function H = solve_lower_triangular(H1, H2)
if isa(H2, 'hmatrix')
    H = hmatrix;
    H.sz = H2.sz;
    if is_leafnode(H1)
        if H1.admissible
            error('SOLVE_LOWER_TRIANGULAR:: admissible block on the diagonal ---> dangerous inversion?');
        else
            if H2.admissible
                H.U = H1.F \ H2.U;
                H.V = H2.V;
                H.admissible = true;
            else
                H.F = H1.F\full(H2);
            end
        end
    else
        if is_leafnode(H2)
            if H2.admissible
                H.U = H1\ H2.U;
                H.V = H2.V;
                H.admissible = true;
            else
                H.F = H1 \ H2.F;
            end
        else
            H.A11 = H1.A11 \ H2.A11;
            H.A12 = H1.A11 \ H2.A12;
            H.A21 = H1.A22 \ (H2.A21 -  H1.A21 * H.A11);
            H.A22 = H1.A22 \ (H2.A22 -  H1.A21 * H.A12);
        end
    end
else
    if is_leafnode(H1)
        if H1.admissible
            error('SOLVE_LOWER_TRIANGULAR:: admissible block on the diagonal ---> dangerous inversion?');
        else
            H = H1.F \ H2;
        end
    else
        n1 = H1.A11.sz(2);
        X1 = H1.A11 \ H2(1:n1, :);
        X2 = H1.A22 \ (H2(n1 + 1:end, :) -  H1.A21 * X1);
        H = [X1; X2];
    end
end
end

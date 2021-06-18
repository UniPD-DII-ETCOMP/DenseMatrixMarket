function H = halr_from_low_rank(H, U, V)

m = size(U, 1);
n = size(V, 1);

if isempty(H)
	H = halr; H.admissible = true; H.sz = [m n];
end

H = halr_low_rank_ric(H, U, V);

end

function H = halr_low_rank_ric(H, U, V)

if is_leafnode(H)
    if H.admissible
        H.U = U;
        H.V = V;
    else
        H.F = U * V';
    end
else
    [m1, n1] = size(H.A11);
    
    H.A11 = halr_low_rank_ric(H.A11, U(1:m1, :), V(1:n1, :));
    H.A12 = halr_low_rank_ric(H.A12, U(1:m1, :), V(n1+1:end,:));
    H.A21 = halr_low_rank_ric(H.A21, U(m1+1:end,:), V(1:n1,:));
    H.A22 = halr_low_rank_ric(H.A22, U(m1+1:end,:), V(n1+1:end,:));
end

end

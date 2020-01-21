function H = solve_lower_triangular(H1, H2)
if isa(H2,'hodlr') %case of hierarchical right-hand side
    if is_leafnode(H1)
        H = H2;
        H.F = H1.F \ H2.F;
    else
        H = H2;
        
        if isempty(H1.A11.U12) && ~is_leafnode(H1.A11)
            H.U12 = solve_lower_triangular(H1.A11, H2.U12);
            H.A11 = solve_lower_triangular(H1.A11, H2.A11);
        else
            H.U12 = H1.A11\H2.U12;
            H.A11 = H1.A11\H2.A11;
        end
        
        H.U21 = [H2.U21, H1.U21];
        H.V21 = [H2.V21, -dense_mtimes_hodlr(H1.V21', H.A11)' ];
        
        if isempty(H1.A22.U12) && ~is_leafnode(H1.A22) % Check this!
            H.U21 = solve_lower_triangular(H1.A22, H.U21);
            H.A22 = solve_lower_triangular(H1.A22, ...
                hodlr_rank_update(H2.A22, -H1.U21 * (H1.V21' * H.U12), H2.V12));
        else
            H.U21 = H1.A22\H.U21;
            H.A22 = H1.A22\hodlr_rank_update(H2.A22, -H1.U21 * (H1.V21' * H.U12), H2.V12);
        end
    end
    
else % case of dense right-hand side
    if is_leafnode(H1)
        H = H1.F\H2;
    else
        mp = H1.A11.sz(2);
        
        if isempty(H1.A11.U12)
            x2 = solve_lower_triangular(H1.A11, H2(1:mp, :));
        else
            x2 = H1.A11 \ H2(1:mp,:);
        end
        
        if isempty(H1.A22.U12)
            x1 = solve_lower_triangular(H1.A22, H2(mp+1:end,:) - H1.U21 * (H1.V21' * x2));
        else
            x1 = H1.A22 \ (H2(mp+1:end,:) - H1.U21 * (H1.V21' * x2));
        end
        
        H = [x2; x1];
    end
end

function H = times(H1, H2)
%TIMES Scalar multiplication

if (isfloat(H1) && isscalar(H1)) || (isfloat(H2) && isscalar(H2)) % if one of the two is a scalar
    	H = H1 * H2;
elseif isa(H1, 'halr') && isa(H2, 'halr')
    	H = halr_hadamard_mul(H1, H2);
else
    error('A .* B: Unsupported operation');
end

end


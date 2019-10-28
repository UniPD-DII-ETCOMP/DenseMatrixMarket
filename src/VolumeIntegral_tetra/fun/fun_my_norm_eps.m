function [B] = fun_my_norm_eps(A,neps)
B = sqrt(A(1).^2 +A(2).^2 + A(3).^2+neps^2);
end


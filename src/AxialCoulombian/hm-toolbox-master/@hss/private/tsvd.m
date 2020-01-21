function [U,S,V] = tsvd(A,tol)
if min(size(A)) == 0
    U = zeros(size(A, 1), 0); S = []; V = zeros(size(A, 2), 0); return;
end
[U,S,V] = svd(A);

t = diag(S);

switch hssoption('norm')
    case 2
        r = sum(t > tol);
    case 'fro'
        tt = sqrt(cumsum(t(end:-1:1).^2));
        r = sum(tt > tol);
end

U = U(:,1:r);
V = V(:,1:r);
S = S(1:r,1:r);

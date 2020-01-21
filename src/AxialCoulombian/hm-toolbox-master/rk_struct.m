function [S, ST] = rk_struct(A, cholesky)
%BUILD_RK_STRUCT Build a struct for building rational Krylov subspaces.
%

S = struct(...
    'solve', @(nu, mu, x) (nu * A - mu * eye(size(A), 'like',A)) \ x, ...
    'multiply', @(rho, eta, x) rho * A * x - eta * x, ...
    'isreal', isreal(A), ...
    'nrm', normest(A, 1e-2));

ST = struct(...
    'solve', @(nu, mu, x) (nu * A' - mu * eye(size(A), 'like', A)) \ x, ...
    'multiply', @(rho, eta, x) rho * A' * x - eta * x, ...
    'isreal', S.isreal, ...
    'nrm', S.nrm);

end


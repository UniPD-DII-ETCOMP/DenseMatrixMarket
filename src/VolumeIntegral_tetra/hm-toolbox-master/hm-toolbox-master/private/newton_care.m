function [X, converged] = newton_care(A, B, C, X0, tol, maxit, debug)
% X = NEWTON_CARE(A,B,C,X0) solves the CARE C + XA + A'X - XBX = 0
% by means of Newton's method with linesearch
%    A, B, C: matrix coefficients
%    X0: initial approximation
%    X : solution of the CARE
if ~exist('tol','var')
	tol = 1e-13;
end
if ~exist('maxit', 'var')
	maxit = 30;
end
if ~exist('debug', 'var')
	debug = false;
end
X = X0;
err = 1; err_old = 1;
k = 1;
RX = C + X0 * A + A' * X0 - X0 * B * X0;
RX = .5 * (RX + RX');

converged = true;
linesearch = true;
while err > tol && k < maxit && err_old >= err 
    H = lyap(A' - X * B, RX); 
    if linesearch
    	V = H*B*H;
    	a = trace(RX*RX);
    	b = trace(RX*V);
    	c = trace(V*V);
    	tk = fminbnd(@(t) a*(1-t)^2-2*b*(1-t)*t^2+c*t^4,0,2);
    else
	tk = 1;
    end
    X = X + tk *H;
    RX = C + X * A + A' * X - X * B * X;
    RX = .5 * (RX + RX');
    err = norm(RX, 'fro') / norm(X, 'fro');
    if mod(k, 2) == 1
	err_old = err;
    end
    k = k + 1;
end

if k == maxit && debug
    fprintf('NEWTON_CARE::Warning: reached the maximum number of iterations: %d, Res = %e\n', maxit, err)
end
if err > 1e-10
	converged = false;
end


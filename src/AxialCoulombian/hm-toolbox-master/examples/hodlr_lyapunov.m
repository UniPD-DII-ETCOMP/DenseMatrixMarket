%% Solving a Lyapunov equation

%% Problem setting
%
% We consider the problem $AX + XA^T = C$, with $A$ symmetric positive
% definite, and $C$ with a HODLR structure. This can be solved by using the
% LYAP function in the HM format. If you are interested in the solver in
% the HSS format, check this <hss_lyapunov.html page>.
%
% Assume that $A = {\mathrm{diag}}_n(-1, 2, -1) \in {C}^{n \times n}$,
% and $C$ is the matrix with the sampling of $g(x,y) = \log(1 + |x - y|)$,
% scaled by a factor $h^2$, with $h = \frac{1}{n + 1}$.
%
% Then, the solution $X$ is a sampling of the function $u(x,y)$ that solves
% the PDE with zero Dirichlet boundary conditions
%
% $$
%   \frac{\partial^2}{\partial x^2} u(x,y) +
%    \frac{\partial^2}{\partial y^2} u(x,y) = \log(1 + |x - y|),
% $$
%
% where $(x,y) \in [0, 1]^2$, and $u(x,y) \equiv 0$ on $\partial [0, 1]^2$.
%
% For this test, we select $n = 2^{11}$.

n = 2^11;
h = 1 / (n + 1);

%% Building the matrices
% We start by constructing the matrices. For $A$, we can use the banded
% constructor:

A = hodlr('banded', spdiags(ones(n, 1) * [ -1, 2, -1 ], -1:1, n, n), 1);

%%
% To build $C$, instead, we rely on adaptive cross approximation

x = linspace(0, 1, n + 2);
C = hodlr('handle', @(i,j) log(1 + abs(x(j+1) - x(i+1)')), n, n);

%%
% Let us briefly check that all these matrices have a low-rank structure in
% their off-diagonal blocks:
fprintf('hodlrrank(A) = %d, hodlrrank(C) = %d\n', hodlrrank(A), hodlrrank(C));

%% Solving the equation
% We can now use LYAP to solve the Lyapunov equation. In this example, we
% use the sign iteration. This is also the default method used by the LYAP
% function.

tic; X = lyap(A, -h^2 * C, 'method', 'sign'); toc

%%
% Let us verify that the solution satifies the differential equation. Here,
% we compute the $L^2$ norm of the residual, relative to the $L^2$ norm of
% the solution $X$.
fprintf('Relative residual: %e\n', ...
    norm(A * X + X * A - h^2 * C, 'fro') / norm(X, 'fro'));

%%
% Let's have a look at the solution; we reduce the number of points since
% we do not need to zoom into the details.
XX = full(X);
x = linspace(0, 1, n);
x = x(1:10:end);
mesh(x, x, XX(1:10:end,1:10:end));

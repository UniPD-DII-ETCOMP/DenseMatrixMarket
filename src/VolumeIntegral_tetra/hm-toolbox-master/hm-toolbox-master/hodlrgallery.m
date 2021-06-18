function varargout = hodlrgallery(name, varargin)
%HODLRGALLERY Generate examples of H-matrices.
%
% H = HODLRGALLERY('laplacian', n) generates a n x n discretization of the
% Laplacian operator on [-1, 1] using a grid of n+1 points, assuming
% Dirichlet boundary conditions.
%
% H = HODLRGALLERY('haber', n) generates a 6n x 6n discretization for the
% Heat diffusion problem described in [1].
%
% H = HODLRGALLERY('rand', n, k) generates an n x n random HODLR matrix
% with HODLR rank equal to K.
%
% [1] Haber, Aleksandar, and Michel Verhaegen. "Sparse solution of the
%     Lyapunov equation for large-scale interconnected systems."
%     Automatica 73 (2016): 256-268.

switch name
    case 'laplacian'
        if isempty(varargin)
            n = 2048;
        else
            n = varargin{1};
        end
        
        varargout{1} = hodlr('tridiagonal', ...
            spdiags(ones(n,1) * [ -1 2 -1 ], -1:1, n, n));
        
    case 'haber'
        if isempty('varargin')
            n = 2048;
        else
            n = varargin{1};
        end
        
        m = 6 * n;
        
        % We first generate the matrices A and P as sparse matrices.
        A = sparse(m, m);
        P = sparse(m, m);
        
        a = -1.36;
        e = 0.34;
        
        % We build the diagonal blocks
        for i = 1 : n
            range = ((i-1)*6+1) : i*6;
            A(range, range) = spdiags(ones(6,1) * [ e, a, e ], -1:1, 6, 6);
            P(range, range) = -.2 * ones(6) - .8 * eye(6);
        end
        
        % And now super and sub-diagonals
        for i = 1 : n - 1
            range = ((i-1)*6+1) : i*6;
            
            A(range, range + 6) = e * speye(6);
            A(range + 6, range) = e * speye(6);
            
            P(range, range + 6) = -.1 * ones(6);
            P(range + 6, range) = -.1 * ones(6);
        end
        
        varargout{1} = hodlr('banded', A, 6);
        varargout{2} = hodlr('banded', P, 11);
    case 'rand'
        if isempty(varargin)
            n = 2048; k = 10;
        elseif length(varargin) == 1
            n = varargin{1}; k = 10;
        else
            n = varargin{1}; k = varargin{2};
        end
        H = hodlr('diagonal', ones(n, 1));
        H = hodlr_rand_ric(H, k);
        varargout{1} = H;
    otherwise
        error('Unsupported example name');
end
end

function H = hodlr_rand_ric(H, k)
if ~isempty(H.F)
    H.F = randn(size(H.F));
else
    H.U12 = randn(size(H.U12, 1), k);
    H.V12 = randn(size(H.V12, 1), k);
    H.U21 = randn(size(H.U21, 1), k);
    H.V21 = randn(size(H.V21, 1), k);
    H.A11 =	hodlr_rand_ric(H.A11, k);
    H.A22 = hodlr_rand_ric(H.A22, k);
end
end


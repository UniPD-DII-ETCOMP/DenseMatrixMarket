classdef hss
%HSS HSS matrices
%
% H = HSS(A) constructs an HSS representation of the matrix A, using the
%     algorithm described in [2]. This procedure has a cost O(n^2), where
%     n is the size of A, provided that the off-diagonal rank is negligible
%     with respect to n. 
%
%     If A is sparse, then the random sampling constructor described in
%     HSS('handle', ...) below is used.  
%
% H = HSS('banded', A) constructs an HSS representation of a banded matrix
%     A. The matrix A can be either sparse or dense. 
%
% H = HSS('banded', A, B) can be used to specify the symmetric bandwidth B 
%     of the matrix A.
%
% H = HSS('banded', A, BL, BU) specifies different lower and upper 
%     bandwidth BL and BU, respectively.
%
% H = HSS('chebfun2', F, XDOM, YDOM, M, N) constructs the M x N matrix
%     containing the samplings of the bivariate function F over a uniform
%     grid of the square XDOM x YDOM. The procedure relies on separable
%     approximation of F(X,Y) as provided by the Chebfun package.
%
% H = HSS('diagonal', D) constructs the diagonal matrix with the entries of
%     the vector D on the main diagonal. 
%
% H = HSS('eye', N) constructs an HSS representation of the N x N identity 
%     matrix. 
%
% H = HSS('handle', AFUN, AFUNT, AEVAL, M, N) constructs an HSS matrix
%     using the random sampling based algorithm in [1]. It requires the
%     handle function AFUN and AFUNT which perform the matrix-vector
%     products A*v and A'*v, respectively, and AEVAL which, given two
%     integer vectors I, J returns the submatrix A(I, J). M and N are the
%     number of rows and columns of A. 
%
% H = HSS('low-rank', U, V) construct an HSS representation of the low-rank
%     matrix U*V'.
%
% H = HSS('ones', M, N) constructs an HSS representation of the rank-1
%     M x N matrix of all ones. 
%
% H = HSS('toeplitz', C, R) constructs the Toeplitz matrix with C as first
%     column and R as first row. The representation is constructed using
%     the 'handle' constructor, and fast Toeplitz-vector multiplication. 
%
% H = HSS('zeros', M, N) constructs the HSS representation of the M x N
%     zero matrix. 
%
% All the constructors support an additional 'cluster' keyword that allows
% to specify custom row and column clusters. These are described as a
% vector of indices J = [J(1), ..., J(2^P)], such that the partitioning at
% the lowest level P is 
%
%        (1, J(1))    (J(1)+1, J(2))   ...   (J(2^(P-1)+1), J(2^P)), 
%
% J(2^P) = N. If J(I) = J(I+1) the corresponding leafnode is assumed to be
% missing from the tree. The cluster can be specified with the syntax
%
%   H = HSS(..., 'cluster', rowcluster, colcluster). 
%
% If colcluster is omitted then it is assumed that rowcluster == colcluster.
%
% The partitioning of an HSS matrix can be retrieved calling CLUSTER(H).  
%
%[1] Martinsson, P. G. (2011). A fast randomized algorithm for computing a
%    hierarchically semiseparable representation of a matrix. SIAM Journal
%    on Matrix Analysis and Applications, 32(4), 1251-1274.
%
%[2] Xia, J., Chandrasekaran, S., Gu, M., & Li, X. S. (2010). Fast 
%    algorithms for hierarchically semiseparable matrices. Numerical 
%    Linear Algebra with Applications, 17(6), 953-976.
    
    properties
        level
        blocktype
        
        % top and bottom blocks of the matrix.
        B12
        B21
        
        % Factorization of the upper triangular block as U12 * V12'
        U
        V
        
        % Factorization of the lower triangular block as U21 * V21'
        Rl
        Rr
        Wl
        Wr
        
        % Size of the matrix
        ml
        nl
        mr
        nr
        
        % Dense version of the matrix, if the size is smaller than the
        % minimum allowed block size.
        D
        
        topnode
        leafnode
        
        A11
        A22
        
    end
    
    methods
        
        function obj = hss(varargin)
            %HSS Create a new Hierarchical matrix.
            if nargin == 0
                return;
            end
            
            rowcluster = [];
            colcluster = [];
            
            % Find the first string parameter after varargin{1}
            charpos = 2;
            while charpos <= nargin && ~ischar(varargin{charpos})
                charpos = charpos + 1;
            end
            
            if charpos <= nargin && strcmp(varargin{charpos}, 'cluster')
                rowcluster = varargin{charpos + 1};
                if nargin >= charpos + 2
                    colcluster = varargin{charpos + 2};
                else
                    colcluster = rowcluster;
                end
            end
                      
            
            if ~ischar(varargin{1})
                A = varargin{1};
                
                if issparse(A)
                    obj = hss('handle', ...
                        @(v) A * v, @(v) A' * v, @(i,j) full(A(i,j)), ...
                        size(A, 1), size(A, 2), 'cluster', ...
                        rowcluster, colcluster);
                else
                    obj = hss_build_hss_tree(size(A, 1), size(A, 2), ...
                                hssoption('block-size'), rowcluster, ...
                                colcluster);
                    obj = hss_from_full(obj, A);
                end
                
                return;
            end
            
            if nargin > 1
                switch varargin{1} % added by Piergiorgio Alotto
                    case 'function'
                        if charpos < 5
                            error('Unsufficient parameters for the handle constructor');
                        end
                        

                        obj = hss_build_hss_tree(varargin{3}, varargin{4}, ...
                                hssoption('block-size'), rowcluster, ...
                                colcluster);                            

                        obj = hss_from_fun(obj, varargin{2});
                        
                    case 'banded'
                        obj = hss_build_hss_tree(size(varargin{2}, 1), ...
                            size(varargin{2}, 2), hssoption('block-size'), ...
                            rowcluster, colcluster);
                        obj = hss_from_banded(obj, varargin{2:charpos-1});

                    case 'cauchy'
                        %obj = hm2hss(hm('cauchy', varargin{2:end}));
                        obj = hss_from_cauchy(varargin{2:end});

                    case 'chebfun2'
                        obj = hm2hss(hm('chebfun2', varargin{2:end}));

                    case 'diagonal'
                        obj = hss_build_hss_tree(length(varargin{2}), ...
                            length(varargin{2}), hssoption('block-size'), ...
                            rowcluster, colcluster);
                        obj = hss_build_diagonal(obj, varargin{2:charpos-1});

                    case 'eye'
                        n = varargin{2};

                        if ~check_cluster_equality(rowcluster, colcluster)
                            error('row and column cluster must match for the identity matrix');
                        end

                        obj = hss('diagonal', ones(n, 1), 'cluster', rowcluster);

                    case 'handle'
                        if charpos < 7
                            error('Unsufficient parameters for the handle constructor');
                        end
                        
                        obj = hss_build_hss_tree(varargin{5}, varargin{6}, ...
                                hssoption('block-size'), rowcluster, ...
                                colcluster);

                        obj = hss_from_random_sampling(obj, varargin{2:charpos-1});

                    case 'low-rank'
                        obj = hss_build_hss_tree(size(varargin{2}, 1), ...
                            size(varargin{3}, 1), hssoption('block-size'), ...
                            rowcluster, colcluster);
                        obj = hss_build_low_rank(obj, varargin{2:charpos-1});

                    case 'ones'
                        m = varargin{2};
                        if charpos > 3
                            n = varargin{3};
                        else
                            n = m;
                        end

                        obj = hss('low-rank', ones(m, 1), ones(n, 1), ...
                            'cluster', rowcluster, colcluster);
                        
                    case 'toeplitz'
                        if charpos == 4
                            m = length(varargin{2});
                            n = length(varargin{3});
                        elseif charpos == 5
                            m = varargin{4};
                            n = m;
                        else
                            m = varargin{4};
                            n = varargin{5};
                        end

                        obj = hss_from_symbol(varargin{2:3}, m, n, ...
                            rowcluster, colcluster);

                    case 'zeros'
                        m = varargin{2};
                        if charpos > 3
                            n = varargin{3};
                        else
                            n = m;
                        end
                        obj = hss_build_hss_tree(m, n, ...
                            hssoption('block-size'), rowcluster, colcluster);
                    otherwise
                        error('Unsupported constructor mode');
                end
            end
        end
        
    end
    
    %
    % Start of the private methods used to instantiate the HSS objects
    %
    methods (Access = private)
        
    end
end

%%
function H = hss_from_fun(obj, fun)
%HSS_FROM_FUN Build an HSS representation of a dense matrix from a function
%providing its entries
%
% This functions is equivalent to calling H = HSS(A) with no other options,
% which is the preferred syntax.
%
% The implemenetation and the algorithm is based on the one presented in
% the paper
%
%     Xia, Jianlin, et al. "Fast algorithms for hierarchically
%     semiseparable matrices." Numerical Linear Algebra with Applications
%     17.6 (2010): 953-976.
%
% The complexity of the algorithm, assuming a low HSS rank, is quadratic in
% the dimension.


mA = size(obj, 1);
nA = size(obj, 2);

if mA ~= nA
    % error('Rectangular HSS matrices are not (yet) supported');
end

tol = hssoption('threshold');

% Prepare the tree for the HSS structure -- leaving all the blocks empty
% H = hss_build_hss_tree(m, n, block_size);
H = obj;

% Create the stack
rs = createStack();
cs = createStack();

row_offsets = [];
col_offsets = [];

H = BuildHSS_iter_fun(H, fun, mA, nA, tol, rs, cs, ...
    row_offsets, col_offsets, 0, 0);


end

function [H, rs, cs, row_offsets, col_offsets, rs1, cs1] = BuildHSS_iter_fun(...
    H, fun, mA, nA, tol, rs, cs, row_offsets, col_offsets, mh, nh)

m = size(H, 1);
n = size(H, 2);

if H.leafnode
    % You might want to uncomment this for debugging purposes
    % fprintf('> Leafnode, mh = %d, nh = %d, m = %d, n = %d\n', mh, nh, m, n);
    
    % We need to extract all the blocks in the stack, and put them
    % together to assemble the matrix.
    
    % Let's do it for the HSS block row first.
    %%%B = A(mh+1:mh+m, nh+n+1:end);
    B=getblock(fun,mh+1,mh+m, nh+n+1,nA);
    
    for j = length(col_offsets) - 1 : -1 : 0
        Z = hss_elementAt(cs, j);
        
        % FIXME: We should preallocate B
        B = [ Z(col_offsets(j+1)+1:col_offsets(j+1)+m,:) , B ];
    end
    
    [H.U, Z] = compress_hss_block(B, tol);
    
    % Push Z in the stack, and updates the elements in there as well
    counter = 0;
    for j = 0 : length(col_offsets) - 1
        W = hss_elementAt(cs, j);
        
        cs = hss_setElementAt(cs, [ W(1:col_offsets(j+1), :) ; ...
            Z(counter+1:counter+size(W, 2), :)' ; ...
            W(col_offsets(j+1)+m+1:end,:) ], j);
        
        counter = counter + size(W, 2);
        col_offsets(j+1) = col_offsets(j+1) + size(Z, 2);
        % col_offsets(j+1) = size(Z, 2);
    end
    
    % Last but not least, push the new compressed stuff into the stack,
    % and add an offset entry for it
    rs = hss_push(rs, Z(counter+1:end, :));
    rs1 = size(Z, 2);
    
    % And now do the columns
    %%%B = A(mh+m+1:end, nh+1:nh+n)';
    B=getblock(fun,mh+m+1,mA, nh+1,nh+n)';
    
    for j = length(row_offsets) - 1 : -1 : 0
        Z = hss_elementAt(rs, j);
        
        % FIXME: We should preallocate B
        B2 = [ Z(row_offsets(j+1)+1:row_offsets(j+1)+n,:) , B ];
        B = B2;
    end
    
    [H.V, Z] = compress_hss_block(B, tol);
    
    % Push Z in the stack, and updates the elements in the as well
    counter = 0;
    for j = 0 : length(row_offsets) - 1
        W = hss_elementAt(rs, j);
        
        rs = hss_setElementAt(rs, [ W(1:row_offsets(j+1), :) ; ...
            Z(counter+1:counter+size(W, 2), :)' ; ...
            W(row_offsets(j+1)+n+1:end,:) ], j);
        
        counter = counter + size(W, 2);
        row_offsets(j+1) = row_offsets(j+1) + size(Z, 2);
    end
    
    % Last but not least, push the new compressed stuff into the stack,
    % and add an offset entry for it
    cs = hss_push(cs, Z(counter+1:end, :));
    cs1 = size(Z, 2);
    
    col_offsets = [ col_offsets , 0 ];
    row_offsets = [ row_offsets , 0 ];
    
    %%%H.D = A(mh+1:mh+m, nh+1:nh+n);
    H.D = getblock(fun,mh+1,mh+m, nh+1,nh+n);
    
else
    % fprintf('> Non-leafnode, mh = %d, nh = %d, m = %d, n = %d\n', mh, nh, m, n);
    % Call the constructor recursively on the left and right childs
    [H.A11, rs, cs, row_offsets, col_offsets] = BuildHSS_iter_fun(...
        H.A11, fun, mA, nA, tol, rs, cs, ...
        row_offsets, col_offsets, mh, nh);
    
    [H.A22, rs, cs, row_offsets, col_offsets] = BuildHSS_iter_fun(...
        H.A22, fun, mA, nA, tol, rs, cs, row_offsets, col_offsets, ...
        mh + size(H.A11, 1), nh + size(H.A11, 2));
    
    % Extract Bl and Bu from the stacks, and merge the children
    [row1, rs] = hss_pop(rs); [row2, rs] = hss_pop(rs);
    [col1, cs] = hss_pop(cs); [col2, cs] = hss_pop(cs);
    
    rs1 = size(row2, 1) - size(row1, 1);
    cs1 = size(col2, 1) - size(col1, 1);
    
    [rU, rV, rlU, rlV, rrU, rrV] = generatorsRank(H);
    
    assert(rs1 == rrV)
    assert(cs1 == rrU)
    
    H.B12 = row2(1:rrV, 1:rlU)';
    H.B21 = col2(1:rrU, 1:rlV);
    
    if H.topnode
        return;
    end
    
    if length(col_offsets) > 2
        col_offsets = col_offsets(1:end-2) - col_offsets(end-2);
    else
        col_offsets = [];
    end
    
    if length(row_offsets) > 2
        row_offsets = row_offsets(1:end-2) - row_offsets(end-2);
    else
        row_offsets = [];
    end
    
    % Merge the rows and cols
    B = [ row2(rrV+1:end, :), row1 ]';
    
    for j = length(col_offsets) - 1 : -1 : 0
        Z = hss_elementAt(cs, j);
        
        % FIXME: We should preallocate B
        B = [ Z(col_offsets(j+1)+1:col_offsets(j+1)+rU,:) , B ];
    end
    
    [U, Z] = compress_hss_block(B, tol);
    
    H.Rl = U(1:rlU, :);
    H.Rr = U(rlU+1:end, :);
    
    % Push Z in the stack, and updates the elements in there as well
    counter = 0;
    for j = 0 : length(col_offsets) - 1
        W = hss_elementAt(cs, j);
        
        cs = hss_setElementAt(cs, [ W(1:col_offsets(j+1), :) ; ...
            Z(counter+1:counter+size(W, 2), :)' ; ...
            W(col_offsets(j+1)+rU+1:end,:) ], j);
        
        counter = counter + size(W, 2);
        col_offsets(j+1) = col_offsets(j+1) + size(Z, 2);
    end
    
    rs = hss_push(rs, Z(counter+1:end, :));
    
    B = [ col2(rrU+1:end, :), col1 ]';
    
    for j = length(row_offsets) - 1 : -1 : 0
        Z = hss_elementAt(rs, j);
        
        % FIXME: We should preallocate B
        B = [ Z(row_offsets(j+1)+1:row_offsets(j+1)+rV,:) , B ];
    end
    
    [U, Z] = compress_hss_block(B, tol);
    
    H.Wl = U(1:rlV, :);
    H.Wr = U(rlV+1:end, :);
    
    % Push Z in the stack, and updates the elements in the as well
    counter = 0;
    for j = 0 : length(row_offsets) - 1
        W = hss_elementAt(rs, j);
        
        rs = hss_setElementAt(rs, [ W(1:row_offsets(j+1), :) ; ...
            Z(counter+1:counter+size(W, 2), :)' ; ...
            W(row_offsets(j+1)+rV+1:end,:) ], j);
        
        counter = counter + size(W, 2);
        row_offsets(j+1) = row_offsets(j+1) + size(Z, 2);
    end
    
    cs = hss_push(cs, Z(counter+1:end, :));
    
    row_offsets = [ row_offsets, 0 ];
    col_offsets = [ col_offsets, 0 ];
end
end

function [rU, rV, rlU, rlV, rrU, rrV] = generatorsRank(H)
if H.A11.leafnode
    rlU = size(H.A11.U, 2);
    rlV = size(H.A11.V, 2);
else
    rlU = size(H.A11.Rl, 2);
    rlV = size(H.A11.Wl, 2);
end

if H.A22.leafnode
    rrU = size(H.A22.U, 2);
    rrV = size(H.A22.V, 2);
else
    rrU = size(H.A22.Rl, 2);
    rrV = size(H.A22.Wl, 2);
end

rU = rlU + rrU;
rV = rlV + rrV;
end



function [U, Z] = compress_hss_block(B, tol)

if isempty(B)
    U = zeros(size(B, 1), 0);
    Z = zeros(size(B, 2), 0);
    return
end

switch hssoption('compression')
    case 'qr'
        [Q, R, P] = qr(B, 0);
        
        % Select only the important columns out of Q
        rk = sum(abs(diag(R)) > abs(R(1,1)) * tol);
        IP = zeros(1, length(P)); IP(P) = 1 : length(P);
        U = Q(:, 1:rk);
        Z = R(1:rk, IP)';
        
    case 'svd'
        [U, B, Z] = svd(B, 'econ');
        
        if tol == 0
            rk = size(B, 1);
        else
            rk = sum(diag(B) > tol * B(1,1));
        end
        
        U = U(:, 1:rk);
        B = B(1:rk, 1:rk);
        
        Z = Z(:, 1:rk) * B;
end
end

% Simple home-made implementation of a stack -- if we really want to use
% this is still open for discussion. It does not matter much
% performance-wise in this context anyway.

function s = createStack()
s = {};
end

function s = hss_push(s, el)
s{end+1} = el;
end

function [el, s] = hss_pop(s)
el = s{end};
s = { s{1:end-1} };
end

function el = hss_elementAt(s, j)
el = s{j+1};
end

function s = hss_setElementAt(s, el, j)
s{j+1} = el;
end

function [M]=getblock(fun,a,b,c,d)

  M=zeros(b-a+1,d-c+1);
  if (d-c)>(b-a)
      for ii=a:b
            M(ii-a+1,1:d-c+1)=fun(ii,c:d); 
      end
  else
          for jj=c:d
              M(1:b-a+1,jj-c+1)=fun(a:b,jj); 
          end
  end
end
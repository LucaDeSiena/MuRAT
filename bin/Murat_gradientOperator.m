function D  =   Murat_gradientOperator(coords)
% function D  =   Murat_gradientOperator(coords)
%
% BUILDS the first-difference (gradient) operator of a regular 3-D grid.
%
% Each row of D takes the difference between two nodes that are adjacent
% along x, y or z, so D*m = 0 for a uniform model and norm(D*m) measures
% model roughness. Used as smoothing regulariser in the Tikhonov
% inversions of Qc and Q.
%
% Input parameters:
%    coords:    N-by-3 node coordinates of a regular grid (any order,
%               e.g. Murat.input.modv(:,1:3))
%
% Output parameters:
%    D:         sparse M-by-N operator, one row per adjacent node pair
%
% Example:
%    D      =   Murat_gradientOperator(Murat.input.modv(:,1:3));
%    Dk     =   D(:,rayCrosses_k);
%    Dk     =   Dk(sum(abs(Dk),2) == 2,:);   % pairs inside the subset

N           =   size(coords,1);
sub         =   zeros(N,3);
dims        =   zeros(1,3);
for j = 1:3
    c           =   coords(:,j);
    tol         =   1e-6 * max(1, max(abs(c)));
    [u,~,idx]   =   uniquetol(c, tol, 'DataScale', 1);
    sub(:,j)    =   idx;
    dims(j)     =   numel(u);
end

node        =   zeros(prod(dims),1);
node(sub2ind(dims, sub(:,1), sub(:,2), sub(:,3)))   =   1:N;

I           =   [];
J           =   [];
V           =   [];
nRows       =   0;
for j = 1:3
    shift       =   zeros(1,3);
    shift(j)    =   1;
    nb          =   sub + shift;
    inside      =   all(nb <= dims, 2);
    a           =   find(inside);
    b           =   node(sub2ind(dims, nb(inside,1), nb(inside,2),...
                    nb(inside,3)));
    valid       =   b > 0;
    a           =   a(valid);
    b           =   b(valid);
    m           =   numel(a);
    rows        =   nRows + (1:m)';
    I           =   [I; rows; rows];                    %#ok<AGROW>
    J           =   [J; a; b];                          %#ok<AGROW>
    V           =   [V; -ones(m,1); ones(m,1)];         %#ok<AGROW>
    nRows       =   nRows + m;
end

D           =   sparse(I, J, V, nRows, N);

end

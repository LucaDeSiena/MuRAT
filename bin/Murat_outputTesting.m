function m          =   Murat_outputTesting(G,d,reguParam,inversionMethod)
% function m        =   Murat_outputTesting(G,d,reguParam,inversionMethod)
%
% PRODUCES the outputs of checkerboard and spike tests
%
% Input parameters:
%    G:                 inversion matrix
%    d:                 data vector
%    reguParam:         regularization parameter
%    inversionMethod:   choice of inversion technique
%    
% Output parameters:
%    m:                 inversion solution

if isequal(inversionMethod,'Tikhonov')
    [U,S,V]         =   svd(G);
    m               =   tikhonov(U,diag(S),V,d,reguParam);
    
elseif isequal(inversionMethod,'Iterative')
    optionsI        =   IRset('MaxIter', 500,'RegMatrix','Identity',...
        'IterBar','off','RegParam', reguParam,'verbosity','off');
    
    [m,~]           =   IRcgls(G,d,optionsI);
    
elseif isequal(inversionMethod,'Hybrid')
    options         =   IRset('RegParam', reguParam,'IterBar','off',...
    'verbosity','off');
    
    [m,~]           =   IRhybrid_lsqr(G, d, options);
    
end
    
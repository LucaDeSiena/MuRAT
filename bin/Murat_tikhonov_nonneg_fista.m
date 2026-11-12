function x  =   Murat_tikhonov_nonneg_fista(Afun, Atfun, Lfun, Ltfun,...
    b, n, lambda, opts)
% Afun:  function handle y = Afun(x)
% Atfun: function handle x = Atfun(y)
% Lfun, Ltfun: handles for L and L'
% b:      data vector
% n:      length of x
% lambda: regularization parameter
% opts:  struct with fields maxit (200), tol (1e-6), verbose (false),
%        x0 (starting model, default zeros; negative values are set to 0)

if ~isfield(opts,'maxit'), opts.maxit = 200; end
if ~isfield(opts,'tol'),   opts.tol   = 1e-6; end
if ~isfield(opts,'verbose'), opts.verbose = false; end
if ~isfield(opts,'x0') || isempty(opts.x0), opts.x0 = zeros(n,1); end

% Gradient of quadratic: grad(x) = H*x - A'*b, where H = A'*A + lambda^2 L'*L
Ab  =   Atfun(b);

% Estimate Lipschitz constant Ls via power method on H (matrix-free)
x0  =   randn(n,1);
x0  =   x0 / norm(x0);
num_power   =   20;
for k = 1:num_power
    y       =   Atfun(Afun(x0)) + (lambda^2) * (Ltfun(Lfun(x0)));
    x0      =   y / norm(y);
end
Ls  = max(1e-6, (x0' * (Atfun(Afun(x0)) + (lambda^2)*(Ltfun(Lfun(x0)))))); 
t   = 1 / Ls;

% FISTA initialization
x   =   max(0, opts.x0(:));       % start at the projected x0 (feasible)
y_k =   x;
tk  =   1;

for k = 1:opts.maxit
    % gradient at y_k
    grad    =   Atfun(Afun(y_k)) + (lambda^2)*(Ltfun(Lfun(y_k))) - Ab;
    x_new   =   y_k - t * grad;
    % projection onto x >= 0
    x_new   =   max(0, x_new);
    tk1     =   (1 + sqrt(1 + 4*tk^2)) / 2;
    y_k     =   x_new + ((tk - 1)/tk1)*(x_new - x);
    % check convergence (relative step size)
    if norm(x_new - x) <= opts.tol * max(norm(x_new), eps)
        x   =   x_new;
        break;
    end
    x       =   x_new;
    tk      =   tk1;
end

if k == opts.maxit
    warning('MuRAT:fistaMaxIterations', ['Tikhonov inversion stopped at ' ...
        'the maximum of %d iterations without converging; increase ' ...
        'MaximumIterations or the damping.'], opts.maxit);
end

end
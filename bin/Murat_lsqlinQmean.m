function [d0, Q0, mDest]   =   Murat_lsqlinQmean(cf_k,l,time0_k,rapsp_k,...
    Qc_k, tCm, te)
% [d0, Q0, mDest]   =   Murat_lsqlinQmean(cf_k,l,time0_k,rapsp_k)
%
% INVERTS with minimum least squares to obtain average Q
%
% Input parameters:
%    cf_f:          central frequency
%    l:             total ray length in km
%    time0_k:       travel time per frequency
%    rapsp_k:       energy ratio per frequency
%
% Output parameters:
%    d0:            data for the average inversion
%    Q0:            average Q, uncertainties
%    mDest:         mean estimated diffusion constant
%%
% Data creation for the true inversion, removing the parameters
% pre-calculated using the diffusion model

data    =   log(rapsp_k.*l.^2)/2/pi/cf_k ;

G       =   ones(length(time0_k),1);
G(:,2)  =   -time0_k;

% Storing inverted parameters
Q0      =   lsqlin(G,data);

x       =   Q0; 
res     =   G*x - data;
m       =   size(G,1); 
n       =   size(G,2);
dof     =   max(m - n, 1);      % degrees of freedom (avoid div by zero)
sigma2  =   (res'*res) / dof;   % estimated residual variance

% If unconstrained or you treat all params as free:
cov_x   =   sigma2 * pinv(G'*G);  % use pinv(G'*G) if ill-conditioned

Q0(:,2) =   diag(cov_x);

d0      =   data-Q0(1);

cEst    =   Q0(1,1);

% Linear grid
Dmin    =   1e-3;
Dmax    =   1e3;
Ngrid   =   2000;
Ds      =   linspace(Dmin, Dmax, Ngrid);   % [1 x Ngrid]
Ds_col  =   Ds(:);                     % [Ngrid x 1]

% Ensure row vectors for broadcasting
Qc_row  =   reshape(Qc_k, 1, []);
l_row   =   reshape(l, 1, []);
tCm_row =   reshape(tCm, 1, []);
te_row  =   reshape(te, 1, []);


% Precompute log(4*pi*D) for Ds_col
log4piD =   log(4*pi*Ds_col);           % OK because Dmin > 0

% Compute A(D) -> [Ngrid x N]
A       =   (tCm_row + te_row).^(-1.5) .* ...
                exp( -l_row.^2 ./ (4 .* Ds_col .* te_row)...
                - 2*pi .* cf_k .* Qc_row .* te_row ) ...
                - tCm_row.^(-1.5) .* ...
                exp( -l_row.^2 ./ (4 .* Ds_col .* tCm_row)...
                - 2*pi .* cf_k .* Qc_row .* tCm_row );

% Valid mask and full RHS computation (cf_k scalar)
valid   =   A > 0;                % [Ngrid x N]
RHS_full=   (1.5 / (2*pi*cf_k)) .* log4piD ...
            - (1   / (2*pi*cf_k)) .* log(A); 
RHS_full(~valid) = NaN;       % avoid invalid log entries

% Residual r = cEst - RHS over grid
R       =   cEst - RHS_full;          % [Ngrid x N]
R_abs   =   abs(R);
R_abs(~isfinite(R_abs)) =   Inf;

% Pick Ds that minimize |residual| for each parameter column
[~, idxMin] = min(R_abs, [], 1);
Dest    =   Ds(idxMin);       % no transpose
mDest   =   mean(Dest);      % scalar mean

end
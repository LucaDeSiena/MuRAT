function g = lcfun(lambda,s,beta,xi,fifth)

% Auxiliary routine for l_corner; computes the NEGATIVE of the curvature.
% Note: lambda may be a vector.  PCH, DTU Compute, Jan. 31, 2015.

% Initialization.
phi = zeros(size(lambda)); dphi = phi; psi = phi; dpsi = phi;
eta = phi; rho = phi;
if length(beta) > length(s)  % A possible least squares residual.
    LS = true;
    rhoLS2 = beta(end)^2;
    beta = beta(1:end-1);
else
    LS = false;
end

% Compute some intermediate quantities.
for i = 1:length(lambda)
  if (nargin==4)
    f  = (s.^2)./(s.^2 + lambda(i)^2);
  else
    f  = s./(s + lambda(i));
  end
  cf = 1 - f;
  eta(i) = norm(f.*xi);
  rho(i) = norm(cf.*beta);
  f1 = -2*f.*cf/lambda(i);
  f2 = -f1.*(3-4*f)/lambda(i);
  phi(i)  = sum(f.*f1.*abs(xi).^2);
  psi(i)  = sum(cf.*f1.*abs(beta).^2);
  dphi(i) = sum((f1.^2 + f.*f2).*abs(xi).^2);
  dpsi(i) = sum((-f1.^2 + cf.*f2).*abs(beta).^2);
end
if LS  % Take care of a possible least squares residual.
    rho = sqrt(rho.^2 + rhoLS2);
end

% Now compute the first and second derivatives of eta and rho
% with respect to lambda;
deta  =  phi./eta;
drho  = -psi./rho;
ddeta =  dphi./eta - deta.*(deta./eta);
ddrho = -dpsi./rho - drho.*(drho./rho);

% Convert to derivatives of log(eta) and log(rho).
dlogeta  = deta./eta;
dlogrho  = drho./rho;
ddlogeta = ddeta./eta - (dlogeta).^2;
ddlogrho = ddrho./rho - (dlogrho).^2;

% Let g = curvature.
g = - (dlogrho.*ddlogeta - ddlogrho.*dlogeta)./...
      (dlogrho.^2 + dlogeta.^2).^(1.5);
function [storeD,residualD,outputSolverD,exitFlagSolverD]   =...
    Murat_diffusion(cf,sp_i,t,tW,r)
% function [storeD,residualD,outputSolverD,exitFlagSolverD]   =...
%   Murat_diffusion(cf,sp_i,t,tW,r)
%
% CALCULATES the diffusion constant
%
% Input parameters:
%    cf:            central frequency
%    sp_i:          envelopes at different frequencies
%    t:             lapse time
%    tW:            coda window
%    r:             ray length
%    
% Outputs:
%   storeD:         best-fit D per frequency
%   residualD:      objective value (squared error) per frequency
%   outputSolverD:  struct with solver info per frequency
%   exitFlagSolverD:struct with exit flag (0/1 info) per frequency
 

% E(x,t1) = W./(4*pi*D*t1)^1.5*exp(-r^2/D/t1)
% E(x,t2) = W./(4*pi*D*t2)^1.5*exp(-r^2/D/t2)
% E(x,t1)/E(x,t2) = (4*pi*D*t2)^1.5/(4*pi*D*t1)^1.5*exp(-r^2/D*[1/t1-1/t2])

nF              = numel(cf);
storeD          =   zeros(1,nF);
residualD       =   zeros(1,nF);
outputSolverD   =   struct();
exitFlagSolverD =   struct();
    
ratioE              =   sp_i(1,:)./sp_i(end,:);
tcW                 =   t + tW;
r                   =   r/1000;
for k = 1:nF
    ratioE_k        =   ratioE(k);
    D               =   optimvar('D',1,LowerBound=10^-7,UpperBound=10^3);
    obj             =   norm(ratioE_k-(4*pi*D*tcW)^1.5/(4*pi*D*t)^1.5*...
        exp(-r^2/D*(1/t-1/tcW)));
    prob            =   optimproblem('Objective',obj);
    x0.D            =   10^-2;
    fld             =   sprintf('Hz%g', cf(k));
    fld             =   strrep(fld, '.', '_');
    
    options         = optimoptions('fmincon','Algorithm','sqp');
    [sol,fval,exitflag,output]  = solve(prob,x0, 'Options',options);

    storeD(k)                   = sol.D;
    residualD(k)                = fval;
    outputSolverD.(fld)         = output;
    exitFlagSolverD.(fld)       = exitflag;

end
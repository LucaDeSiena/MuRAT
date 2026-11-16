function [sol,fval,eflag,output,dampValue] = Murat_inversionQc(Ac_k,...
    Qm_k,inversionMethod,iter,iterStall,RZZ_k,coordPrior,dampValue,...
    smoothValue,PlotI,pathFolder,D)
% function [sol,fval,eflag,output,dampValue] = Murat_inversionQc(Ac_k,...
%   Qm_k,inversionMethod,iter,iterStall,RZZ_k,coordPrior,dampValue,...
%   smoothValue,PlotI,pathFolder,D)
%
% INVERTS spatial values of Qc using pre-compiled kernels
%
% Input parameters:
%    Ac_k:              differential sensitivity kernel
%    Qm_k:              measured coda Q
%    inversionMethod:   choose between three optimizations
%    iter:              iterations for the non-linear optimization
%    iterStall:         max iterations for the non-linear optimization
%    RZZ_k:             std to weight inversion matrix
%    coordPrior:        coords for the starting solution
%    dampValue:         damping value, relative to the largest singular
%                       value of the (weighted) inversion matrix: the
%                       penalty is (dampValue*norm(A))^2*norm(m)^2, so 1
%                       weighs model norm and data misfit equally
%    smoothValue:       smoothing value. For 'Tikhonov' it is the weight
%                       of the first-difference regulariser D (empty or 0:
%                       no smoothing), also relative to norm(A); for the
%                       other methods it is the correlation length used
%                       by Murat_smoothing
%    PlotI:             plotting iterations during inversion
%    pathFolder:        folder to save inversion figure for Tikhonov
%    D:                 (optional) first-difference operator for the
%                       columns of Ac_k (see Murat_gradientOperator)
%
% Output parameters:
%    [sol,fval,exitflag,output]: standard optimization outputs

%%
if nargin < 12, D = []; end
useTikhonov =   strcmp(inversionMethod,'Tikhonov');
smoothTikh  =   useTikhonov && ~isempty(smoothValue) && smoothValue ~= 0 ...
                && ~isempty(D);

%Parameter model number
MQc         =   size(Ac_k,2);

%Covariance matrix for the data - weight by correlation coefficient squared
w           =   RZZ_k.^2;
%w           =   (w - min(w(:))) ./ (max(w(:)) - min(w(:)));
W           =   diag(w);


%Starting model from regionalization
Ac_boxes    =   sum(Ac_k,1);
m0Qc        =   (sum(Ac_k.*Qm_k,1)')./Ac_boxes';

if ~useTikhonov && ~isempty(smoothValue) && smoothValue~=0
    sigmaM  =   mean(m0Qc/10)^2;
    CovM    =   Murat_smoothing(coordPrior,smoothValue,sigmaM);
else
    CovM    =   0;
end

if ~isempty(dampValue)
    dampValue   =   abs(dampValue);
elseif smoothTikh
    %Smoothing alone regularises the problem
    dampValue   =   0;
else
    %Default relative damping
    dampValue   =   0.1;
end

%Damping and smoothing are relative to the largest singular value of the
%weighted matrix, so the same values act alike for Qc and Q
normA       =   norm(W*Ac_k);
lambdaDamp  =   dampValue*normA;


% --- optimization variables and objective ---
Qc          =   optimvar('Qc',MQc,'LowerBound',min(Qm_k));
x0.Qc       =   m0Qc;

obj0        =   norm(W*(Qm_k - Ac_k*m0Qc))^2 + ...
                lambdaDamp^2*norm(m0Qc)^2 +...
                norm(CovM*m0Qc);

obj         =   norm(W*(Qm_k - Ac_k*Qc))^2 + ...
                lambdaDamp^2*norm(Qc)^2 +...
                norm(CovM*Qc);

prob        =   optimproblem('Objective',obj);
maxStall    =   'MaxStallIterations';
maxIt       =   'MaxIterations';
switch inversionMethod

    case 'Tikhonov'
        [sol.Qc,fval,~] =...
            Murat_tikhonovQc(PlotI,W*Ac_k,W*Qm_k,lambdaDamp,x0.Qc,iter,...
            smoothValue*normA,D);
        eflag           =   'Tikhonov';
        output          =   [];
        saveFigureAsImage(pathFolder);
        savefig(gcf, [pathFolder '.fig']);
        close(gcf);
        fval            =   fval*obj0;

    case 'Particle'
        if PlotI == 1
            options     =   optimoptions(@particleswarm,maxStall,...
                iterStall,maxIt,iter,'PlotFcn','pswplotbestf');

        else
            options     =   optimoptions(@particleswarm,maxStall,...
                iterStall,maxIt,iter);

        end
        [sol,fval,eflag,output]  =  solve(prob,x0,'Solver',...
            'particleswarm','Options',options);


    case 'Annealing'
        if PlotI == 1
            options     =   optimoptions(@simulannealbnd,maxStall,...
                iterStall,maxIt,iter,'PlotFcns',...
                {@saplotbestx,@saplotbestf,@saplotx,@saplotf});
        else
            options     =   optimoptions(@simulannealbnd,maxStall,...
                iterStall,maxIt,iter);

        end
        [sol,fval,eflag,output] =   solve(prob,x0,"Solver",...
            "simulannealbnd",'Options',options);

    case 'Global'
        if PlotI == 1
            options     =   optimoptions(@ga,'PlotFcn','gaplotbestf');
            [sol,fval,eflag,output] = solve(prob, "Solver",...
                "ga",'Options',options);
        else
            [sol,fval,eflag,output] =   solve(prob, "Solver","ga");

        end

    otherwise
        error('Unknown inversion method.')
end

fval                            =   fval/obj0;
end
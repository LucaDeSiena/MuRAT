function [sol,fval,exitflag,output,dampValue]=  Murat_inversionQ(A_k,d0,...
    rapsp_k,inversionMethod,Q0,iter,iterStall,coordPrior,dampValue,...
    smoothValue,PlotI,pathFolder,D)
% function [sol,fval,exitflag,output,dampValue]=Murat_inversionQ(A_k,d0,...
%   rapsp_k,inversionMethod,Q0,iter,iterStall,coordPrior,dampValue,...
%   smoothValue,PlotI,pathFolder,D)
%
% INVERTS spatial values of Q using ray approximation
%
% Input parameters:
%    A_k:               linearised kernel from ray approximation
%    d0:                data out of average Q computation
%    rapsp_k:           energy ratios (no longer used; kept so existing
%                       calls keep working)
%    inversionMethod:   choose between three optimizations
%    Q0:                average Q with variance
%    iter:              iterations for the non-linear optimization
%    iterStall:         max iterations for the non-linear optimization
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
%                       columns of A_k (see Murat_gradientOperator)
%
% Output parameters:
%    [sol,fval,exitflag,output]: standard optimization outputs

%%
if nargin < 13, D = []; end
useTikhonov =   strcmp(inversionMethod,'Tikhonov');
smoothTikh  =   useTikhonov && ~isempty(smoothValue) && smoothValue ~= 0 ...
                && ~isempty(D);

%Parameter model number
MQ          =   size(A_k,2);

%Starting model from regionalization
d           =   d0./sum(A_k,2);
Ac_boxes    =   abs(A_k./sum(A_k,1));
Qmean       =   Q0(2,1);
m0Q         =   zeros(MQ,1);
for i = 1:MQ
    n       =   Ac_boxes(:,i);
    m       =   n(n>0);
    d1      =   d(n>0);
    m0Q(i)  =   sum(m.*d1);
end

if ~useTikhonov && ~isempty(smoothValue) && smoothValue~=0
    %Variance for models to square of 1/10 of the mean inverse Qc value
    sigmaM  =   Qmean;
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
%matrix, so the same values act alike for Qc and Q
normA       =   norm(A_k);
lambdaDamp  =   dampValue*normA;

Q       =   optimvar('Q',MQ,LowerBound=10^-7);
x0Q.Q   =   Qmean*ones(MQ,1);
obj0    =   norm(d0 - A_k*m0Q)^2 + ...
            lambdaDamp^2*norm(m0Q)^2 +...
            norm(CovM*m0Q);

objQ    =   norm(d0 - A_k*Q)^2 +...
            lambdaDamp^2*norm(Q)^2 + ...
            norm(CovM*Q);

probQ   =   optimproblem('Objective',objQ);
maxStall=   'MaxStallIterations';
maxIt   =   'MaxIterations';

switch inversionMethod

    case 'Tikhonov'
        [sol.Q,fval,~]  =   Murat_tikhonovQ(PlotI,A_k,d0,lambdaDamp,m0Q,iter,...
            smoothValue*normA,D);
        sol.Q           =   sol.Q;
        exitflag        =   'Tikhonov';
        output          =   [];
        saveFigureAsImage(pathFolder);
        savefig(gcf, [pathFolder '.fig']);
        close(gcf);
    
    case 'Particle'
        if PlotI == 1
        options     =   optimoptions(@particleswarm,maxStall,...
            iterStall,maxIt,iter,'PlotFcn','pswplotbestf');
        else
        options     =   optimoptions(@particleswarm,maxStall,...
            iterStall,maxIt,iter);
        
        end
        [sol,fval,exitflag,output]  =   solve(probQ,x0Q,'Solver',...
            'particleswarm','Options',options);
        fval                        =   fval/obj0;
        
    case 'Annealing'
        if PlotI == 1
        options     =   optimoptions(@simulannealbnd,maxStall,...
            iterStall,maxIt,iter,'PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf});
        else
        options     =   optimoptions(@simulannealbnd,maxStall,...
            iterStall,maxIt,iter);
        
        end
        [sol,fval,exitflag,output]  =   solve(probQ,x0Q,"Solver",...
            "simulannealbnd",'Options',options);
        fval                        =   fval/obj0;

    case 'Global'
        if PlotI == 1
        options = optimoptions(@ga,'PlotFcn','gaplotbestf');
        [sol,fval,exitflag,output]  =   solve(probQ, "Solver",...
            "ga",'Options',options);
        else
        [sol,fval,exitflag,output]  =   solve(probQ, "Solver","ga");
        
        end
        fval                        =   fval/obj0;

    otherwise
        error('Unknown inversion method.')
end

end
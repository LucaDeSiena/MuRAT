function [sol,fval,exitflag,output,dampValue]=  Murat_inversionQ(A_k,d0,...
    rapsp_k,inversionMethod,Q0,iter,iterStall,coordPrior,dampValue,...
    smoothValue,PlotI,pathFolder)
% function [sol,fval,exitflag,output,dampValue]=Murat_inversionQ(A_k,d0,...
%   inversionMethod,Q0,iter,iterStall,coordPrior,dampValue,...
%   smoothValue,PlotI,pathFolder)
%
% INVERTS spatial values of Q using ray approximation
%
% Input parameters:
%    A_k:               linearised kernel from ray approximation
%    d0:                data out of average Q computation
%    inversionMethod:   choose between three optimizations
%    Q0:                average Q with variance
%    iter:              iterations for the non-linear optimization
%    iterStall:         max iterations for the non-linear optimization
%    coordPrior:        coords for the starting solution
%    dampValue:         damping value
%    smoothValue:       smoothing value
%    PlotI:             plotting iterations during inversion
%    pathFolder:        folder to save inversion figure for Tikhonov
%
% Output parameters:
%    [sol,fval,exitflag,output]: standard optimization outputs

%%
%Parameter model number
MQ      =   size(A_k,2);

%Starting model from regionalization
d = d0./sum(A_k,2);
Ac_boxes    =   abs(A_k./sum(A_k,1));
m0Q = zeros(MQ,1);
for i = 1:MQ
    n = Ac_boxes(:,i);
    m = n(n>0);
    d1 = d(n>0);
    m0Q(i) = sum(m.*d1);
end

if smoothValue~=0
    %Variance for models to square of 1/10 of the mean inverse Qc value
    sigmaM      =   Q0(1);
    CovM        =   Murat_smoothing(coordPrior,smoothValue,sigmaM);
else
    CovM        =   0;
end


if ~isempty(dampValue)
    dampValue   =   abs(dampValue);
else
    %Damping is the mean energy ratio
    dampValue   =   0.1*mean(rapsp_k);
end

Q       =   optimvar('Q',MQ,LowerBound=10^-7);
x0Q.Q   =   Q0(1)*ones(MQ,1);
obj0    =   norm(d0 - A_k*m0Q)^2 + ...
            dampValue*norm(m0Q)^2 +...
            norm(CovM*m0Q);

objQ    =   norm(d0 - A_k*Q)^2 +...
            dampValue*norm(Q)^2 + ...
            norm(CovM*Q);

probQ   =   optimproblem('Objective',objQ);
maxStall=   'MaxStallIterations';
maxIt   =   'MaxIterations';

switch inversionMethod

    case 'Tikhonov'
        [sol.Q,fval,~]  =   Murat_tikhonovQ(PlotI,A_k,d0,dampValue,m0Q);
        sol.Q           =   sol.Q;
        exitflag        =   'Tikhonov';
        output          =   [];
        saveFigureAsImage(pathFolder);
        savefig(gcf, [pathFolder '.fig']);
        close(gcf);
        fval            =   fval*obj0;
    
    case 'Particle'
        if PlotI == 1
        options = optimoptions(@particleswarm,maxStall,...
            iterStall,maxIt,iter,'PlotFcn','pswplotbestf');
        else
        options = optimoptions(@particleswarm,maxStall,...
            iterStall,maxIt,iter);
        
        end
        [sol,fval,exitflag,output]  = solve(probQ,x0Q,'Solver',...
            'particleswarm','Options',options);

    case 'Annealing'
        if PlotI == 1
        options = optimoptions(@simulannealbnd,maxStall,...
            iterStall,maxIt,iter,'PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf});
        else
        options = optimoptions(@simulannealbnd,maxStall,...
            iterStall,maxIt,iter);
        
        end
        [sol,fval,exitflag,output] = solve(probQ,x0Q,"Solver",...
            "simulannealbnd",'Options',options);

    case 'Global'
        if PlotI == 1
        options = optimoptions(@ga,'PlotFcn','gaplotbestf');
        [sol,fval,exitflag,output] = solve(probQ, "Solver",...
            "ga",'Options',options);
        else
        [sol,fval,exitflag,output] = solve(probQ, "Solver","ga");
        
        end
        
    otherwise
        error('Unknown inversion method.')
end

fval                            =   fval/obj0;
end
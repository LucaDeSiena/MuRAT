function [sol,fval,eflag,output,dampValue] = Murat_inversionQc(Ac_k,...
    Qm_k,inversionMethod,iter,iterStall,RZZ_k,coordPrior,dampValue,...
    smoothValue,PlotI,pathFolder)
% function [sol,fval,eflag,output]=Murat_inversionQc(Ac_k,Qm_k,inversionMethod)
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
%    dampValue:         damping value
%    smoothValue:       smoothing value
%    PlotI:             plotting iterations during inversion
%    pathFolder:        folder to save inversion figure for Tikhonov
%
% Output parameters:
%    [sol,fval,exitflag,output]: standard optimization outputs

%%
%Parameter model number
MQc         =   size(Ac_k,2);

%Covariance matrix for the data
w           =   1 ./ RZZ_k;
%w           =   (w - min(w(:))) ./ (max(w(:)) - min(w(:))); 
W           =   diag(w.^2);


%Starting model from regionalization
Ac_boxes    =   sum(Ac_k,1);
m0Qc        =   (sum(Ac_k.*Qm_k,1)')./Ac_boxes';

if smoothValue~=0
    sigmaM  =   mean(m0Qc/10)^2;
    CovM    =   Murat_smoothing(coordPrior,smoothValue,sigmaM);
else
    CovM    =   0;
end

if ~isempty(dampValue)
    dampValue   =   abs(dampValue);
else
    %Damping is also 1/10 of the mean inverse Qc value
    dampValue   =   mean(m0Qc/10);
end


% --- optimization variables and objective ---
Qc          =   optimvar('Qc',MQc,'LowerBound',min(Qm_k));
x0.Qc       =   m0Qc;

obj0        =   norm(W*(Qm_k - Ac_k*m0Qc))^2 + ...
                dampValue*norm(m0Qc)^2 +...
                norm(CovM*m0Qc);

obj         =   norm(W*(Qm_k - Ac_k*Qc))^2 + ...
                dampValue*norm(Qc)^2 +...
                norm(CovM*Qc);

prob        =   optimproblem('Objective',obj);
maxStall    =   'MaxStallIterations';
maxIt       =   'MaxIterations';
switch inversionMethod

    case 'Tikhonov'
        [sol.Qc,fval,~] =...
            Murat_tikhonovQc(PlotI,W*Ac_k,W*Qm_k,dampValue,x0.Qc,iter);
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
            [sol,fval,~,output] =   solve(prob, "Solver","ga");

        end

    otherwise
        error('Unknown inversion method.')
end

fval                            =   fval/obj0;
end
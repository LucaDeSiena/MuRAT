function [image, minimizeVector,infoVector,minimizeValue]...
                                =   Murat_minimise(outputLCurve,G,bQ,name)
% function [image, minimizeVector,infoVector,minimizeValue]...
%                               =   Murat_minimise(G,Q,W,name)
% MINIMIZES both L2 and I cost functions with conjugate gradients
% CALCULATES the damping value and PLOTS the costs.
%
% Input parameters:
%    G:                 inversion matrix
%    Q:                 data vector
%    W:                 weighting matrix
%    name:              name of the figure
%    
% Output parameters:
%    image:             L-curves and cost functions
%    minimizeVector:    best solution minimizing both costs
%    infoVector:        informations related to minimization
%    minimizeValue:     regularization parameter

M                               =   size(G,2);
epsilon                         =   zeros(65,1);

for j=1:59
    if j<=6
        epsilon(j,1)            =   10^-9*10^j;
    elseif j>6 && j<=15
        epsilon(j,1)            =   10^-3 + 10^-3*(j-6);
    elseif j>15 && j<=24
        epsilon(j,1)            =   10^-2 + 10^-2*(j-15);
    elseif j>24 && j<=33
        epsilon(j,1)            =   10^-1 + 10^-1*(j-24);
    elseif j>33 && j<=42
        epsilon(j,1)            =   1 + j - 33;
    elseif j>42 && j<=51
        epsilon(j,1)            =   10 + (j-42)*10;    
    elseif j>51 && j<=60
        epsilon(j,1)            =   100 + (j-51)*100;   
    end
    
end
epsilon(60:65)                  =   10.^(3:8);
Nit                             =   length(epsilon);

epsilonRND                      =   zeros(Nit,7);
mestI                           =   zeros(M,Nit);
mestL2                          =   zeros(M,Nit);

epsilonRND(:,1)                 =   epsilon;
for i=1:Nit
    epsilon_i                   =   epsilon(i,1);
    
    optionsI                    =...
        IRset('MaxIter', 500,'RegMatrix','Identity','IterBar','off',...
        'verbosity','off','RegParam', epsilon_i);
        
    [mestI_i,infoI]             =   IRcgls(G,bQ,optionsI);
    epsilonRND(i,2)             =	infoI.StopReg.Rnrm;
    epsilonRND(i,3)             =   infoI.StopReg.Xnrm;
    epsilonRND(i,4)             =   infoI.StopReg.NE_Rnrm;
    
    optionsL2                   =...
        IRset('MaxIter', 500,'RegMatrix','Laplacian1D','IterBar','off',...
        'verbosity','off','RegParam', epsilon_i);
    
    [mestL2_i,infoL2]           =   IRcgls(G,bQ,optionsL2);
    
    epsilonRND(i,5)             =	infoL2.StopReg.Rnrm;
    epsilonRND(i,6)             =   infoL2.StopReg.Xnrm;
    epsilonRND(i,7)             =   infoL2.StopReg.NE_Rnrm;
    
    mestI(:,i)                  =   mestI_i;
    mestL2(:,i)                 =   mestL2_i;
end

costFunctionI                   =   0.5*sum(epsilonRND(:,2:4).^2,2);
[minCostI,indexCostI]           =   min(costFunctionI);

disp(['Minimum of I cost function is ',num2str(minCostI)...
    ' for epsilonI = ' num2str(epsilon(indexCostI))]);

costFunctionL2                  =   0.5*sum(epsilonRND(:,5:7).^2,2);

[minCostL2,indexCostL2]         =   min(costFunctionL2);

disp(['Minimum of L2 cost function is ',num2str(minCostL2)...
    ' for epsilonL2 = ' num2str(epsilon(indexCostL2))]);

image                           =   figure('Name',name,...
    'NumberTitle','off','Position',[20,400,1200,1000]);

subplot(2,2,1)
loglog(epsilonRND(:,2),epsilonRND(:,3),'k','LineWidth',2);
title('I regularization')
xlabel('Model perturbation');
ylabel('Data fit');
for k = 1:2:Nit
    text(epsilonRND(k,2),epsilonRND(k,3),num2str(epsilon(k)));
end
SetFDefaults()

subplot(2,2,2)
loglog(epsilonRND(:,5),epsilonRND(:,6),'k','LineWidth',2);
title('L2 regularization')
xlabel('Model perturbation');
ylabel('Data fit');
for k = 1:2:Nit
    text(epsilonRND(k,5),epsilonRND(k,6),num2str(epsilon(k)));
end
SetFDefaults()

subplot(2,2,3)
loglog(epsilon,costFunctionI,'k','LineWidth',2);
title('I Cost')
xlabel('Regularization Parameter');
ylabel('Data fit');
SetFDefaults()

subplot(2,2,4)
loglog(epsilon,costFunctionL2,'k','LineWidth',2);
title('L2 Cost')
xlabel('Regularization Parameter');
ylabel('Data fit');
SetFDefaults()

if outputLCurve == 1
    minimizeValue                               =...
        input('Your damping parameter for coda: ');
else
    minimizeValue                               =...
        mean(indexCostL2,indexCostI);
end


options                         =...
    IRset('MaxIter', 500,'RegMatrix','Identity','IterBar','off',...
    'RegParam', minimizeValue);

[minimizeVector,infoVector]     =   IRcgls(G,bQ,options);

end
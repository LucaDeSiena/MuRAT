function CN_analysis        =   Murat_imageCheckCN(rtQk,rcQk,cfk,...
        outputLCurve,energyRatio_k,residualQ_k,const_Qc_k,Edirect_k,...
        luntot,time0,A_i,lCurveQ,CN_title)
% function CN_analysis        =   Murat_imageCheckCN(rtQk,rcQk,cfk,...
%         outputLCurve,energyRatio_k,residualQ_k,const_Qc_k,Edirect_k,...
%         luntot,time0,A_i,lCurveQ,CN_title)
% 
% PLOTS the Qc checks
%
% Input parameters:
%    rtQk:          retained Q data
%    rcQk:          retained Q parameter models
%    cfk:           central frequency
%    outputLCurve:  decide if you want to see the L curve
%    energyRatio_k: energy ratios at all frequencies
%    residualQ_k:   residuals of Q inversions
%    const_Qc_k:    the constants coming from Qc analysis
%    Edirect_k:     direct energy
%    luntot:        total length of ray
%    time0:         travel time
%    A_i:         	CN inversion matrix
%    lcurveQ:       damping parameter
%    CN_title:      title of the figure
%
% Output parameters:
%    CN_analysis:   figure for coda normalization check

CN_analysis                 =   figure('Name',...
    CN_title,'NumberTitle','off','Position',[20,400,1200,1000]);

A                           =   A_i(rtQk,rcQk);
[~,~,~,~,d1k,spreadAverageQ]=   Murat_tikhonovQ(cfk,rtQk,outputLCurve,...
    energyRatio_k,const_Qc_k,luntot,time0,A,lCurveQ,0);

equationQ                   =   -log(const_Qc_k)...
    -spreadAverageQ(1,1)*log(luntot(rtQk))...
    -2*pi*cfk*time0(rtQk)*spreadAverageQ(2,1);

dRatio                      =	log(energyRatio_k);
luntotQ                     =   luntot(rtQk)/1000;
[U,S,~]                     =   svd(A);

%%
%Plot of the left and right hand sides of the CN equation.
subplot(3,1,1)
plot(time0(rtQk),dRatio,'o','MarkerSize',6,'MarkerEdgeColor',[0 0 0])
hold on
plot(time0(rtQk),equationQ,'r*','MarkerSize',6)
hold off
xlabel('Travel time (s)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Logarithmic direct-to-coda energy ratio','FontSize',12,...
    'FontWeight','bold','Color','k')
title(['Dependence of logarithmic energy ratios on travel time. '...
    'Geometrical spreading is: ' num2str(spreadAverageQ(1,1))...
    '+/- ' num2str(spreadAverageQ(1,2))]);
legend({'Q^{-1}',cat(2,'<Q> = ',num2str(1/spreadAverageQ(2,1)))},...
    'Location','northeast')
SetFDefaults()
%%
%Plot of the direct energy versus distance.
subplot(3,1,2)
plot(log(luntotQ),log(Edirect_k),'o','MarkerSize',6,...
    'MarkerEdgeColor',[0 0 0])
xlabel('Hypocentral distance (km)','FontSize',12,...
    'FontWeight','bold','Color','k')
ylabel('Logarithmic direct energy (J/m^2)','FontSize',12,...
    'FontWeight','bold','Color','k')
title('Dependence of logarithmic direct energy on hypocentral distance.');

xti                         =   xticks;
xt                          =   cell(length(xti),1);
for i = 1:length(xti)
    xt(i,1)                 =   {exp(xti(i))};
end
xticklabels(xt)
SetFDefaults()
%%
% Then it checks the inversion with a Picard
subplot(3,1,3)
picard(U,diag(S),d1k);
title(['Picard condition. Residual is: ' num2str(residualQ_k)]);
SetFDefaults()
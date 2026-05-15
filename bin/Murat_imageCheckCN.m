function CN_analysis    =   Murat_imageCheckCN(time0_k,Q0,d0,Ed_k,CN_title)
% function CN_analysis  =   Murat_imageCheckCN(time0_k,Q0,d0,Ed_k,CN_title)
% 
% PLOTS the Qc checks
%
% Input parameters:
%    time0_k:       travel time
%    Q0:            average Q with uncertainty
%    d0:            data for average Q inversion
%    Ed_k:     direct energy
%    CN_title:      title of the figure
%
% Output parameters:
%    CN_analysis:   figure for coda normalization check

CN_analysis             =   myfig(CN_title);

equationQ               =   -time0_k*Q0(2,1);
%%
%Plot of the left and right hand sides of the CN equation.
subplot(2,1,1)
plot(time0_k,d0,'o','MarkerSize',6,'MarkerEdgeColor',[0 0 0])
hold on
plot(time0_k,equationQ,'r-')
hold off
xlabel('Travel time (s)','FontSize',10,'FontWeight','bold','Color','k')
ylabel('Corrected log-energy ratio','FontSize',10,...
    'FontWeight','bold','Color','k')
title('Corrected log-energy ratio vs travel time');
legend({'Q^{-1}',cat(2,'<Q> = ',num2str(Q0(2,1)),'+- ',...
    num2str(Q0(2,2)))},'Location','northeast')

SetFDefaults()
%%
%Plot of the direct energy versus distance.
subplot(2,1,2)
plot(log(time0_k),log(Ed_k),'o','MarkerSize',6,...
    'MarkerEdgeColor',[0 0 0])
xlabel('Travel time (s)','FontSize',10,...
    'FontWeight','bold','Color','k')
ylabel('Log-energy (J/m^2)','FontSize',10,...
    'FontWeight','bold','Color','k')
title('Log. direct energy vs travel time.');

xti                     =   xticks;
xt                      =   cell(length(xti),1);
for i = 1:length(xti)
    xt(i,1)             =   {exp(xti(i))};
end
xticklabels(xt)
SetFDefaults()
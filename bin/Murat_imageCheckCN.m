function CN_analysis    =   Murat_imageCheckCN(equationQ,residualQ_k,d1,...
    spreadAverageQ,luntot_k,time0_k,energyRatio_k,A_k,Edirect_k,CN_title)
% function CN_analysis  =   Murat_imageCheckCN(equationQ,d1,...
%     spreadAverageQ,luntot_k,time0_k,energyRatio_k,A_k,Edirect_k,CN_title)
% 
% PLOTS the Qc checks
%
% Input parameters:
%    equationQ:     sum of all three terms of the CN equation for average Q
%    d1:            data for average Q inversion
%    spreadAverageQ:values of average spreading and Q for labels
%    luntot_k:      total length of ray
%    time0_k:       travel time
%    energyRatio_k: energy ratios
%    A_k:         	CN inversion matrix
%    Edirect_k:     direct energy
%    CN_title:      title of the figure
%
% Output parameters:
%    CN_analysis:   figure for coda normalization check

CN_analysis             =   figure('Name',CN_title,'NumberTitle','off',...
    'Position',[20,400,1200,1000],'visible','off');

dRatio                  =	log(energyRatio_k);
luntotQ                 =   luntot_k/1000;
[U,S,~]                 =   svd(A_k);

%%
%Plot of the left and right hand sides of the CN equation.
subplot(3,1,1)
plot(time0_k,dRatio,'o','MarkerSize',6,'MarkerEdgeColor',[0 0 0])
hold on
plot(time0_k,equationQ,'r*','MarkerSize',6)
hold off
xlabel('Travel time (s)','FontSize',10,'FontWeight','bold','Color','k')
ylabel('Log-energy ratio','FontSize',10,'FontWeight','bold','Color','k')
title(['Log-energy ratios vs travel time - \gamma = '...
    num2str(spreadAverageQ(1,1)) '+/- ' ...
    num2str(spreadAverageQ(1,2))]);
legend({'Q^{-1}',cat(2,'<Q> = ',num2str(1/spreadAverageQ(2,1)))},...
    'Location','northeast')

SetFDefaults()
%%
%Plot of the direct energy versus distance.
subplot(3,1,2)
plot(log(luntotQ),log(Edirect_k),'o','MarkerSize',6,...
    'MarkerEdgeColor',[0 0 0])
xlabel('Hypocentral distance (km)','FontSize',10,...
    'FontWeight','bold','Color','k')
ylabel('Log-energy (J/m^2)','FontSize',10,...
    'FontWeight','bold','Color','k')
title('Log. direct energy vs hypocentral distance.');

xti                     =   xticks;
xt                      =   cell(length(xti),1);
for i = 1:length(xti)
    xt(i,1)             =   {exp(xti(i))};
end
xticklabels(xt)
SetFDefaults()
%%
% Then it checks the inversion with a Picard
subplot(3,1,3)
picard(U,diag(S),d1);
title(['Picard condition - Residual = ' num2str(residualQ_k)]);
SetFDefaults()
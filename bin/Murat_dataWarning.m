function  [problempd,problemQc,problemRZZ,problemQ,yes_pd]...
                                        =...
            Murat_dataWarning(listaSac,nonlinear,tresholdnoise,...
            maPD,miPD,fT,peakd,Qm,RZZ,rapspcn)
% Warns about problems with the data and locates indices where this happens
lcf                                     =   size(Qm);
problempd                               =   cell(1,lcf(2));                                
problemQc                               =   cell(1,lcf(2));                                
problemRZZ                              =   cell(1,lcf(2));                                
problemQ                                =   cell(1,lcf(2));                                

yes_pd                                  =   peakd > miPD & peakd < maPD;
no_pd                                   =   peakd > maPD | peakd < miPD;
no_Qc                                   =   Qm == 0;
no_Q                                    =   rapspcn<tresholdnoise;

switch nonlinear
    case 0
        no_RZZ                          =   RZZ <= fT;
    case 1
        no_RZZ                          =   RZZ >= fT;
end

for i = 1:lcf(2)
    problempd{i}                        =   listaSac(no_pd(:,i));
    problemQc{i}                        =   listaSac(no_Qc(:,i));
    problemRZZ{i}                       =   listaSac(no_RZZ(:,i));
    problemQ{i}                         =   listaSac(no_Q(:,i));
end

%Conditions for the linear and non-linear inversions
displayNoQc                             =   sum(no_Qc)/lcf(1)*100;
disp(['In your frequency range, [',num2str(displayNoQc),...
    '] % of your Qc are = 0']);

displayNoRZZ                            =   sum(no_RZZ)/lcf(1)*100;
disp(['In your frequency range, [',num2str(displayNoRZZ),...
    ']% of your correlation coefficients are below accuracy treshold']);

displayNoPD                             =   sum(no_pd)/lcf(1)*100;
disp(['In your frequency range, [',num2str(displayNoPD),...
    '] % of your peak delays are outside peak delay limits']);

displayNoQ                              =   sum(no_Q)/lcf(1)*100;
disp(['In your frequency range, [',num2str(displayNoQ),...
    ']% of your coda-to-noise ratios are below noise treshold']);
end
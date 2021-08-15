function  [problempd,problemQc,problemRZZ,problemQ,yes_pd,compMissing,...
    flag]   =   Murat_dataWarning(listaSac,tresholdnoise,...
    maPD,miPD,fT,peakd,Qm,RZZ,rapspcn,comp,flag)
% function  [problempd,problemQc,problemRZZ,problemQ,yes_pd,compMissing,...
%     flag]   =   Murat_dataWarning(listaSac,tresholdnoise,...
%     maPD,miPD,fT,peakd,Qm,RZZ,rapspcn,comp,flag)
%
% WARNS about problems with the data and locates indices where this happens
%
% Input parameters:
%    listaSac:          list of SAC files
%    tresholdnoise:     allowed noise level for coda
%    maPD:              maximum peak delay allowed
%    miPD:              minimum peak delay allowed
%    fT:                treshold for uncertainty on Qc
%    peakd:             peak delay values
%    Qm:                coda attenuation values
%    RZZ:               uncertainty on Qc
%    rapspcn:           coda to noise ratio
%    comp:           	components
%    flag:           	flag to change between one or more components
%
% Output parameters:
%    problempd:         stores the data not used because of peak delays
%    problemQc:         stores the data not used because of Qc
%    problemRZZ:        stores the data not used because of RZZ
%    problemQ:          stores the data not used because of Q
%    yes_pd:            good data because of peak delay
%    compMissing:       missing components after processing
%    flag:              flag to discriminate one-component from three

lcf                                     =   size(Qm);
problempd                               =   cell(1,lcf(2));
problemQc                               =   cell(1,lcf(2));
problemRZZ                              =   cell(1,lcf(2));
problemQ                                =   cell(1,lcf(2));

yes_pd                                  =   peakd > miPD & peakd < maPD;
no_pd                                   =   peakd > maPD | peakd < miPD;
no_Qc                                   =   Qm == 0;
no_Q                                    =   rapspcn < tresholdnoise;

compMissing(:,:,1)                      =   no_pd;
compMissing(:,:,2)                      =   no_Qc;
compMissing(:,:,3)                      =   no_Q;

no_RZZ                                  =   RZZ <= fT;

for i = 1:lcf(2)
    problempd{i}                        =   listaSac(no_pd(:,i));
    problemQc{i}                        =   listaSac(no_Qc(:,i));
    problemRZZ{i}                       =   listaSac(no_RZZ(:,i));
    problemQ{i}                         =   listaSac(no_Q(:,i));
end

% Displays different messages in case of more than 1 component
if comp == 1
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
    
else
    switch flag
        case 0
            
            displayNoQc                     =   sum(no_Qc)/lcf(1)*100;
            disp(['In your frequency range, [',num2str(displayNoQc),...
                '] % of your original Qc are = 0']);
            
            displayNoRZZ                    =   sum(no_RZZ)/lcf(1)*100;
            disp(['In your frequency range, [',num2str(displayNoRZZ),...
                ']% of your original correlation coefficients are below treshold']);
            
            displayNoPD                     =   sum(no_pd)/lcf(1)*100;
            disp(['In your frequency range, [',num2str(displayNoPD),...
                '] % of your original peak delays are outside peak delay limits']);
            
            displayNoQ                      =   sum(no_Q)/lcf(1)*100;
            disp(['In your frequency range, [',num2str(displayNoQ),...
                ']% of your original coda-to-noise ratios are below treshold']);
            
            disp(['Processing to see how many data you have when considering '...
                str2double(comp) ' components'])
            
            flag                            =   1;
        case 1
            
            displayNoQc                     =   sum(no_Qc)/lcf(1)*100;
            disp(['In your frequency range, [',num2str(displayNoQc),...
                '] % of your Qc are = 0']);
            
            displayNoRZZ                    =   sum(no_RZZ)/lcf(1)*100;
            disp(['In your frequency range, [',num2str(displayNoRZZ),...
                ']% of your correlation coefficients are below accuracy treshold']);
            
            displayNoPD                     =   sum(no_pd)/lcf(1)*100;
            disp(['In your frequency range, [',num2str(displayNoPD),...
                '] % of your peak delays are outside peak delay limits']);
            
            displayNoQ                      =   sum(no_Q)/lcf(1)*100;
            disp(['In your frequency range, [',num2str(displayNoQ),...
                ']% of your coda-to-noise ratios are below noise treshold']);
    end
end
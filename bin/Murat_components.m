function [peakd,Qm,RZZ,rapsp,rapspcn]	=...
    Murat_components(comp,peakd1,Qm1,RZZ1,rapsp1,rapspcn1,...
    compMissing)
% function [peakd,Qm,RZZ,rapsp,rapspcn]	=...
%     Murat_components(comp,peakd1,Qm1,RZZ1,rapsp1,rapspcn1,...
%     compMissing)
%
% CREATES data vector for two- and three-component recording
%
% Input parameters:
%    comp:          origin times in seconds
%    peakd1:        peak delay before averaging
%    Qm1:           coda attenuation before averaging
%    RZZ1:       	uncertainty on Qc before averaging
%    rapsp1:        energy ratio before averaging
%    rapspcn1:      coda to noise ratio before averaging
%    compMissing:   takes into account which component values are missing
%
% Output parameters:
%    peakd:         peak delay after averaging
%    Qm:            coda attenuation after averaging
%    RZZ:       	uncertainty on Qc after averaging
%    rapsp:         energy ratio after averaging
%    rapspcn:       coda to noise ratio after averaging

[dataL,lcf]             =   size(peakd1);
nGroups                 =   dataL / comp;
idxOut                  =   0;

% Preallocate outputs
peakd                   =   nan(nGroups, lcf);
Qm                      =   nan(nGroups, lcf);
RZZ                     =   nan(nGroups, lcf);
rapsp                   =   nan(nGroups, lcf);
rapspcn                 =   nan(nGroups, lcf);

for i = 1:comp:dataL
    idxOut              =   idxOut+1;
    rows                =   i:(i+comp-1);
    
    % masks: valid (not missing) for each metric, size comp x lcf
    validPD             =   ~squeeze(compMissing(rows,:,1)); % comp x lcf
    validQc             =   ~squeeze(compMissing(rows,:,2));
    validRap            =   ~squeeze(compMissing(rows,:,3));
    
    % Extract data blocks comp x lcf
    blockPD             =   peakd1(rows, :);
    blockQm             =   Qm1(rows, :);
    blockRZZ            =   RZZ1(rows, :);
    blockRap            =   rapsp1(rows, :);
    blockRapC           =   rapspcn1(rows, :);
    
    % For each column, average available values; if none available -> NaN
    % peakd
    sumPD               =   sum(blockPD .* validPD, 1);
    nPD                 =   sum(validPD, 1);
    peakd(idxOut, :)    =   sumPD ./ nPD;
    peakd(idxOut, nPD == 0) =   NaN;

    % Qm and RZZ (independent averaging)
    sumQ                =   sum(blockQm .* validQc, 1);
    nQ                  =   sum(validQc, 1);
    Qm(idxOut, :)       =   sumQ ./ nQ;
    Qm(idxOut, nQ == 0) =   NaN;
    
    sumR                =   sum(blockRZZ .* validQc, 1); 
    RZZ(idxOut, :)      =   sumR ./ nQ;
    RZZ(idxOut, nQ == 0)=   NaN;

    % rapsp and rapspcn
    sumRp               =   sum(blockRap .* validRap, 1);
    nRp                 =   sum(validRap, 1);
    rapsp(idxOut, :)    =   sumRp ./ nRp;
    rapsp(idxOut, nRp == 0) =   NaN;
    
    sumRpc              =   sum(blockRapC .* validRap, 1);
    rapspcn(idxOut, :)  =   sumRpc ./ nRp;
    rapspcn(idxOut, nRp == 0) = NaN;

end


function  [problempd,problemQc,problemRZZ,problemQ,yes_pd,compMissing,...
    flag]   =   Murat_dataWarning(listaSac,tresholdnoise,...
    maPD,miPD,fT,peakd,Qm,RZZ,rapspcn,rapsp,comp,flag,QcM)
% function  [problempd,problemQc,problemRZZ,problemQ,yes_pd,compMissing,...
%     flag]   =   Murat_dataWarning(listaSac,tresholdnoise,...
%     maPD,miPD,fT,peakd,Qm,RZZ,rapspcn,rapsp,comp,flag,QcM)
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
%    rapsp:             P/S to coda ratio
%    comp:           	components
%    flag:           	flag to change between one or more components
%    QcM:           	chosen method to measure Qc between Lin e NonLin
%
% Output parameters:
%    problempd:         stores the data not used because of peak delays
%    problemQc:         stores the data not used because of Qc
%    problemRZZ:        stores the data not used because of RZZ
%    problemQ:          stores the data not used because of Q
%    yes_pd:            good data because of peak delay
%    compMissing:       missing components after processing
%    flag:              flag to discriminate one-component from three

% sizes
[nRows, nCols]      = size(Qm);

% logical masks
yes_pd              =   (peakd > miPD) & (peakd < maPD);
no_pd               =   ~yes_pd;
no_Qc               =   (Qm == 0);

if isequal(QcM,'Linearized')
    no_RZZ              =   (RZZ <= fT);

elseif isequal(QcM,'NonLinear')
    no_RZZ              =   (RZZ >= fT);

end



no_Q                =   (rapspcn < tresholdnoise)|(rapsp < tresholdnoise);

compMissing         =   false(nRows, nCols, 3);
compMissing(:,:,1)  =   no_pd;
compMissing(:,:,2)  =   no_Qc;
compMissing(:,:,3)  =   no_Q;

problempd           =   cell(1,nCols);
problemQc           =   cell(1,nCols);
problemRZZ          =   cell(1,nCols);
problemQ            =   cell(1,nCols);

for c = 1:nCols
    idx_pd          =   no_pd(:,c);
    idx_Qc          =   no_Qc(:,c);
    idx_RZZ         =   no_RZZ(:,c);
    idx_Q           =   no_Q(:,c);

    problempd{c}    =   listaSac(idx_pd);
    problemQc{c}    =   listaSac(idx_Qc);
    problemRZZ{c}   =   listaSac(idx_RZZ);
    problemQ{c}     =   listaSac(idx_Q);
    
end

% Percentages (per frequency / column)
pctNoQc             =   sum(no_Qc)  / nRows * 100;
pctNoRZZ            =   sum(no_RZZ) / nRows * 100;
pctNoPD             =   sum(no_pd)  / nRows * 100;
pctNoQ              =   sum(no_Q)   / nRows * 100;


% Displays different messages in case of more than 1 component
if comp == 1 && flag ~= 2
    % When single component processing and flag not equal 2, show overall percentages
    fprintf('[%s] %% of Qc are below threshold\n', num2str(pctNoQc));
    fprintf('[%s] %% of Qc uncertainties are below threshold\n', num2str(pctNoRZZ));
    fprintf('[%s] %% of peak delays are outside limits\n', num2str(pctNoPD));
    fprintf('[%s] %% of ratios are below threshold\n', num2str(pctNoQ));
    flag = 2;

else
    % compute means once
    meanQc  = mean(pctNoQc);
    meanRZZ = mean(pctNoRZZ);
    meanPD  = mean(pctNoPD);
    meanQ   = mean(pctNoQ);

    switch flag
        case 0
            fprintf('[%g] %% of Qc are below threshold\n', meanQc);
            fprintf('[%g] %% of Qc uncertainties are below threshold\n', meanRZZ);
            fprintf('[%g] %% of peak delays are outside limits\n', meanPD);
            fprintf('[%g] %% of ratios are below threshold\n', meanQ);

            if isnumeric(comp)
                ncompStr = num2str(comp);
            else
                ncompStr = comp;
            end
            fprintf('Processing to see how many data you have when considering %s components\n', ncompStr);
            flag = 1;

        case 1
            fprintf('[%g] %% of Qc are below threshold\n', meanQc);
            fprintf('[%g] %% of Qc uncertainties are below threshold\n', meanRZZ);
            fprintf('[%g] %% of peak delays are outside limits\n', meanPD);
            fprintf('[%g] %% of ratios are below threshold\n', meanQ);
    end
end
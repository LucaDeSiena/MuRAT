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

[dataL,lcf]                                     =   size(peakd1);
index                                           =   0;

for i = 1:comp:dataL
    index                                       =   index+1;
    
    for j = 1:lcf
        if comp ==  2
            
            noPD                                =   compMissing(i,j,1);
            noPD1                               =   compMissing(i+1,j,1);
            noQc                                =   compMissing(i,j,2);
            noQc1                               =   compMissing(i+1,j,2);
            noQ                                 =   compMissing(i,j,3);
            noQ1                                =   compMissing(i+1,j,3);
            
            if noPD && noPD1
                peakd(index,j)                  =   NaN; %#ok<*AGROW>
                
            elseif noPD && ~noPD1
                peakd(index,j)                  =   peakd1(i+1,j);
                
            elseif ~noPD && noPD1
                peakd(index,j)                  =   peakd1(i,j);
                
            else
                peakd(index,j)                  =...
                    (peakd1(i,j) + peakd1(i+1,j))/2;
                
            end
            
            if noQc && noQc1
                Qm(index,j)                     =   NaN;
                RZZ(index,j)                    =   NaN;
                
            elseif noQc && ~noQc1
                Qm(index,j)                     =   Qm1(i+1,j);
                RZZ(index,j)                    =   RZZ1(i+1,j);
                
            elseif ~noQc && noQc1
                Qm(index,j)                     =   Qm1(i,j);
                RZZ(index,j)                    =   RZZ1(i,j);
                
            else
                Qm(index,j)                     =...
                    (Qm1(i,j) + Qm1(i+1,j))/2;
                RZZ(index,j)                    =...
                    (RZZ1(i,j) + RZZ1(i+1,j))/2;
                
            end
            
            if noQ && noQ1
                rapsp(index,j)                  =   NaN;
                rapspcn(index,j)                =   NaN;
                
            elseif noQc && ~noQc1
                rapsp(index,j)                  =   rapsp1(i+1,j);
                rapspcn(index,j)                =   rapspcn1(i+1,j);
                
            elseif ~noQc && noQc1
                rapsp(index,j)                  =   rapsp1(i,j);
                rapspcn(index,j)                =   rapspcn1(i,j);
                
            else
                rapsp(index,j)                  =...
                    (rapsp1(i,j) + rapsp1(i+1,j))/2;
                rapspcn(index,j)                =...
                    (rapspcn1(i,j) + rapspcn1(i+1,j))/2;
                
            end
            
            
        elseif comp == 3
            
            noPD                                =   compMissing(i,j,1);
            noPD1                               =   compMissing(i+1,j,1);
            noPD2                               =   compMissing(i+2,j,1);
            noQc                                =   compMissing(i,j,2);
            noQc1                               =   compMissing(i+1,j,2);
            noQc2                               =   compMissing(i+2,j,2);
            noQ                                 =   compMissing(i,j,3);
            noQ1                                =   compMissing(i+1,j,3);
            noQ2                                =   compMissing(i+2,j,3);
            
            if noPD && noPD1 && noPD2
                peakd(index,j)                  =   NaN;
                
            elseif noPD && noPD1 && ~noPD2
                peakd(index,j)                  =   peakd1(i+2,j);
                
            elseif noPD && ~noPD1 && noPD2
                peakd(index,j)                  =   peakd1(i+1,j);
                
            elseif ~noPD && noPD1 && noPD2
                peakd(index,j)                  =   peakd1(i,j);
                
            elseif ~noPD && ~noPD1 && noPD2
                peakd(index,j)                  =...
                    (peakd1(i,j) + peakd1(i+1,j))/2;
            elseif ~noPD && noPD1 && ~noPD2
                peakd(index,j)                  =...
                    (peakd1(i,j) + peakd1(i+2,j))/2;
            elseif noPD && ~noPD1 && ~noPD2
                peakd(index,j)                  =...
                    (peakd1(i+1,j) + peakd1(i+2,j))/2;
            else
                peakd(index,j)                  =...
                    (peakd1(i,j) + peakd1(i+1,j) + peakd1(i+2,j))/3;
                
            end
            
            if noQc && noQc1 && noQc2
                Qm(index,j)                     =   NaN;
                RZZ(index,j)                    =   NaN;
                
            elseif noQc && noQc1 && ~noQc2
                Qm(index,j)                     =   Qm1(i+2,j);
                RZZ(index,j)                    =   RZZ1(i+2,j);
                
            elseif noQc && ~noQc1 && noQc2
                Qm(index,j)                     =   Qm1(i+1,j);
                RZZ(index,j)                    =   RZZ1(i+1,j);
                
            elseif ~noQc && noQc1 && noQc2
                Qm(index,j)                     =   Qm1(i,j);
                RZZ(index,j)                    =   RZZ1(i,j);
                
            elseif ~noQc && ~noQc1 && noQc2
                Qm(index,j)                     =...
                    (Qm1(i,j) + Qm1(i+1,j))/2;
                RZZ(index,j)                    =...
                    (RZZ1(i,j) + RZZ1(i+1,j))/2;
            elseif ~noQc && noQc1 && ~noQc2
                Qm(index,j)                     =...
                    (Qm1(i,j) + Qm1(i+2,j))/2;
                RZZ(index,j)                    =...
                    (RZZ1(i,j) + RZZ1(i+2,j))/2;
            elseif noQc && ~noQc1 && ~noQc2
                Qm(index,j)                     =...
                    (Qm1(i+1,j) + Qm1(i+2,j))/2;
                RZZ(index,j)                    =...
                    (RZZ1(i+1,j) + RZZ1(i+2,j))/2;
            else
                Qm(index,j)                     =...
                    (Qm1(i,j) + Qm1(i+1,j) + Qm1(i+2,j))/3;
                RZZ(index,j)                    =...
                    (RZZ1(i,j) + RZZ1(i+1,j) + RZZ1(i+2,j))/3;
                
            end
            
            if noQ && noQ1 && noQ2
                rapsp(index,j)                  =   NaN;
                rapspcn(index,j)                =   NaN;
                
            elseif noQ && noQ1 && ~noQ2
                rapsp(index,j)                  =   rapsp1(i+2,j);
                rapspcn(index,j)                =   rapspcn1(i+2,j);
                
            elseif noQ && ~noQ1 && noQ2
                rapsp(index,j)                  =   rapsp1(i+1,j);
                RZZ(index,j)                    =   rapspcn1(i+1,j);
                
            elseif ~noQ && noQ1 && noQ2
                rapsp(index,j)                  =   rapsp1(i,j);
                rapspcn(index,j)                =   rapspcn1(i,j);
                
            elseif ~noQ && ~noQ1 && noQ2
                rapsp(index,j)                  =...
                    (rapsp1(i,j) + rapsp1(i+1,j))/2;
                rapspcn(index,j)                =...
                    (rapspcn1(i,j) + rapspcn1(i+1,j))/2;
            elseif ~noQ && noQ1 && ~noQ2
                rapsp(index,j)                  =...
                    (rapsp1(i,j) + rapsp1(i+2,j))/2;
                rapspcn(index,j)                =...
                    (rapspcn1(i,j) + rapspcn1(i+2,j))/2;
            elseif noQ && ~noQ1 && ~noQ2
                rapsp(index,j)                  =...
                    (rapsp1(i+1,j) + rapsp1(i+2,j))/2;
                rapspcn(index,j)                =...
                    (rapspcn1(i+1,j) + rapspcn1(i+2,j))/2;
            else
                rapsp(index,j)                  =...
                    (rapsp1(i,j) + rapsp1(i+1,j) + rapsp1(i+2,j))/3;
                rapspcn(index,j)                =...
                    (rapspcn1(i,j) + rapspcn1(i+1,j) +...
                    rapspcn1(i+2,j))/3;
                
            end
            
        end
        
    end
    
end
end


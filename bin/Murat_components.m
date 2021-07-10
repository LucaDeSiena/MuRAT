function [peakd,Qm,RZZ,rapsp,rapspcn]	=...
    Murat_components(comp,peakd1,Qm1,RZZ1,rapsp1,rapspcn1,...
    compMissing)
%% Creates data vector for two- and three-component recording

[dataL,lcf]                                     =   size(peakd1);
index                                           =   0;
for i = 1:comp:dataL
    index                                       =   index+1;
    % 2 components (typically [WE, SN]) - average data that fullfill criteria
    
    for j = 1:lcf
        if comp ==  2
            
            noPD                                =   compMissing(i,j,1);
            noPD1                               =   compMissing(i+1,j,1);
            noQc                                =   compMissing(i,j,2);
            noQc1                               =   compMissing(i+1,j,2);
            noQ                                 =   compMissing(i,j,3);
            noQ1                                =   compMissing(i+1,j,3);
            
            if noPD && noPD1
                peakd(index,j)                  =   0; %#ok<*AGROW>
                
            elseif noPD && ~noPD1
                peakd(index,j)                  =   peakd1(i+1,j);
                
            elseif ~noPD && noPD1
                peakd(index,j)                  =   peakd1(i+1,j);
                
            else
                peakd(index,j)                  =...
                    (peakd1(i,j) + peakd1(i+1,j))/2;
                
            end
            
            if noQc && noQc1
                Qm(index,j)                     =   0;
                RZZ(index,j)                    =   0;
                
            elseif noQc && ~noQc1
                Qm(index,j)                     =   Qm1(i+1,j);
                RZZ(index,j)                    =   RZZ1(i+1,j);
                
            elseif ~noQc && noQc1
                Qm(index,j)                     =   Qm1(i+1,j);
                RZZ(index,j)                    =   RZZ1(i+1,j);
                
            else
                Qm(index,j)                     =...
                    (Qm1(i,j) + Qm1(i+1,j))/2;
                RZZ(index,j)                    =...
                    (RZZ1(i,j) + RZZ1(i+1,j))/2;
                
            end
            
            if noQ && noQ1
                rapsp(index,j)                  =   0;
                rapspcn(index,j)                =   0;
                
            elseif noQc && ~noQc1
                rapsp(index,j)                  =   rapsp1(i+1,j);
                rapspcn(index,j)                =   rapspcn1(i+1,j);
                
            elseif ~noQc && noQc1
                rapsp(index,j)                  =   rapsp1(i+1,j);
                rapspcn(index,j)                =   rapspcn1(i+1,j);
                
            else
                rapsp(index,j)                  =...
                    (rapsp1(i,j) + rapsp1(i+1,j))/2;
                rapspcn(index,j)                =...
                    (rapspcn1(i,j) + rapspcn1(i+1,j))/2;
                
            end
            
            
        elseif comp == 3
            
            %Averaging the 3 components
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
                peakd(index,j)                  =   0;
                
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
                    (peakd1(i+1,j) + peakd1(i+1,j) + peakd1(i+2,j))/3;
            
            end
            
            if noQc && noQc1 && noQc2
                Qm(index,j)                     =   0;
                RZZ(index,j)                    =   0;
                
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
                    (Qm1(i+1,j) + Qm1(i+1,j) + Qm1(i+2,j))/3;
                RZZ(index,j)                    =...
                    (RZZ1(i+1,j) + RZZ1(i+1,j) + RZZ1(i+2,j))/3;
                
            end
            
            if noQ && noQ1 && noQ2
                rapsp(index,j)                  =   0;
                rapspcn(index,j)                =   0;
                
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
                    (rapsp1(i+1,j) + rapsp1(i+1,j) + rapsp1(i+2,j))/3;
                rapspcn(index,j)                =...
                    (rapspcn1(i+1,j) + rapspcn1(i+1,j) +...
                    rapspcn1(i+2,j))/3;
                
            end
            
        end
        
    end
   
end



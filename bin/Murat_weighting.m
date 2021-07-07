function Wc             =   Murat_weighting(nonlinear,RZZ_k)
%% Creates weighting for the coda method
% Weighting depends on the values of RZZ - this is done only for the
% linear case.
    if nonlinear == 0
        
        W1              =   RZZ_k<0.3;
        W2              =   RZZ_k<0.5;
        W3              =   RZZ_k<0.7;
        W4              =   W3+W2+W1;
        Wc              =   diag(1./(W4+1));
        
    elseif nonlinear == 1
        
        Wc              =   1;
        
    end
end
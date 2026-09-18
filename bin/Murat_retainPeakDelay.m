function [pab, lpdelta_i, retain_pd_i, ray_crosses_pd_i] = ...
         Murat_retainPeakDelay(t_phase, l10pd_i, yesPD_i, Apd_i, stDevPD)  
% function  [pab,lpdelta_i,retain_pd_i,ray_crosses_pd_i] = ...  
%             Murat_retainPeakDelay(t_phase,l10pd_i,yesPD_i,Apd_i, stDevPD)  
%  
% CREATES all constraints for peak delay inversion  
%  
% Input parameters:  
%    t_phase:           time relative to envelope  
%    l10pd_i:           logarithm of peak delay  
%    yesPD_i:           sets waveform where this must be computed  
%    Apd_i:             peak delay matrix  
%    stDevPD:           standard deviation around trending line
%  
% Output parameters:  
%    pab:               coefficients of the linear relationship  
%    lpdelta_i:         residuals  
%    retain_pd_i:       keeps tab on which waveforms are kept for imaging  
%    ray_crosses_pd_i:  keeps tab on which rays are kept for imaging  

% Perform robust linear fitting on the data  
fitr_i  =   fit(t_phase(yesPD_i), l10pd_i(yesPD_i),'poly1','Robust','on');  
pab     =   [fitr_i.p1 fitr_i.p2];  

% Calculate fitted values and residuals  
lpdelta_i   =   l10pd_i - polyval(pab, t_phase); 

% Standard Deviation (3-sigma) method for outlier detection  
mu          =   mean(lpdelta_i(yesPD_i), 'omitnan');  
sigma       =   std(lpdelta_i(yesPD_i), 'omitnan');  

lB          =   mu - stDevPD * sigma;  
uB          =   mu + stDevPD * sigma;  

outliersS   =   (lpdelta_i(yesPD_i) < lB) | (lpdelta_i(yesPD_i) > uB);

% Determine which waveforms to retain  
retain_pd_i         =   false(size(yesPD_i));  
retain_pd_i(yesPD_i)=   ~outliersS;

% rays crossing retained waveforms: sum over retained rows of Apd_i
s_pd                =   sum(Apd_i(retain_pd_i, :), 1);
ray_crosses_pd_i    =   (s_pd ~= 0);

end
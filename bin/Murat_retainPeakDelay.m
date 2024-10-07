% function  [pab,lpdelta_i,retain_pd_i,ray_crosses_pd_i]...
%                                 =...
%             Murat_retainPeakDelay(t_phase,l10pd_i,yesPD_i,Apd_i)
% % function  [pab,lpdelta_i,retain_pd_i,ray_crosses_pd_i]...
% %                                 =...
% %             Murat_retainPeakDelay(t_phase,l10pd_i,yesPD_i,Apd_i)
% %
% % CREATES all constraints for peak delay inversion
% %
% % Input parameters:
% %    t_phase:           time relative to envelope
% %    l10pd_i:           logarithm of peak delay
% %    yesPD_i:           sets waveform where this must be computed
% %    Apd_i:             peak delay matrix
% %
% % Output parameters:
% %    pab:               coefficients of the linear relationship
% %    lpdelta_i:         residuals
% %    retain_pd_i:       keeps tab on which waveforms are kept for imaging
% %    ray_crosses_pd_i:  keeps tab on which rays are kept for imaging

% dataL                           =   length(t_phase);
% outlierspd                      =   false(dataL,1);
% fitrobust_i                     =...
%     fit(t_phase(yesPD_i),l10pd_i(yesPD_i),'poly1','Robust','on');
% pab                             =   [fitrobust_i.p1 fitrobust_i.p2];

% l10pdt                          =   polyval(pab,t_phase);
% lpdelta_i                       =   l10pd_i-l10pdt;
% I                               =   abs(lpdelta_i) >= 2*std(lpdelta_i,'omitnan');
% outliers                        =...
%     excludedata(t_phase(yesPD_i),lpdelta_i(yesPD_i),'indices',I(yesPD_i));
% outlierspd(yesPD_i)             =   outliers;
% retain_pd_i                     =   yesPD_i & ~outlierspd;
% Apd_retain_pd_i                 =   Apd_i(retain_pd_i,:);
% s_pd                            =   sum(Apd_retain_pd_i);
% ray_crosses_pd_i                =   s_pd~=0;

% end


function [pab, lpdelta_i, retain_pd_i, ray_crosses_pd_i] ...
                                = ...
         Murat_retainPeakDelay(t_phase, l10pd_i, yesPD_i, Apd_i)  
% function  [pab,lpdelta_i,retain_pd_i,ray_crosses_pd_i] = ...  
%             Murat_retainPeakDelay(t_phase,l10pd_i,yesPD_i,Apd_i)  
%  
% CREATES all constraints for peak delay inversion  
%  
% Input parameters:  
%    t_phase:           time relative to envelope  
%    l10pd_i:           logarithm of peak delay  
%    yesPD_i:           sets waveform where this must be computed  
%    Apd_i:             peak delay matrix  
%  
% Output parameters:  
%    pab:               coefficients of the linear relationship  
%    lpdelta_i:         residuals  
%    retain_pd_i:       keeps tab on which waveforms are kept for imaging  
%    ray_crosses_pd_i:  keeps tab on which rays are kept for imaging  

% Get the length of the input data  
dataL = length(t_phase);  

% Initialize outlier flag array  
outlierspd = false(dataL, 1);  

% Perform robust linear fitting on the data  
fitrobust_i = fit(t_phase(yesPD_i), l10pd_i(yesPD_i), 'poly1', 'Robust', 'on');  
pab = [fitrobust_i.p1 fitrobust_i.p2];  

% Calculate fitted values and residuals  
l10pdt = polyval(pab, t_phase);  
lpdelta_i = l10pd_i - l10pdt;  

% Standard Deviation (3-sigma) method for outlier detection  
mu = mean(lpdelta_i(yesPD_i), 'omitnan');  
sigma = std(lpdelta_i(yesPD_i), 'omitnan');  
lower_bound = mu - 1 * sigma;  
upper_bound = mu + 1 * sigma;  

% Flag outliers  
I = (lpdelta_i < lower_bound) | (lpdelta_i > upper_bound);  

% Update outlier flags  
outlierspd(yesPD_i) = I(yesPD_i);  

% Determine which waveforms to retain  
retain_pd_i = yesPD_i & ~outlierspd;  

% Process retained waveform data  
Apd_retain_pd_i = Apd_i(retain_pd_i, :);  
s_pd = sum(Apd_retain_pd_i);  
ray_crosses_pd_i = s_pd ~= 0;  

end
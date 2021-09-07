function  [pab,lpdelta_i,retain_pd_i,ray_crosses_pd_i]...
                                =...
            Murat_retainPeakDelay(t_phase,l10pd_i,yesPD_i,Apd_i)
% function  [pab,lpdelta_i,retain_pd_i,ray_crosses_pd_i]...
%                                 =...
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

dataL                           =   length(t_phase);
outlierspd                      =   false(dataL,1);
fitrobust_i                     =...
    fit(t_phase(yesPD_i),l10pd_i(yesPD_i),'poly1','Robust','on');
pab                             =   [fitrobust_i.p1 fitrobust_i.p2];

l10pdt                          =   polyval(pab,t_phase);
lpdelta_i                       =   l10pd_i-l10pdt;
I                               =   abs(lpdelta_i) >= 2*std(lpdelta_i);
outliers                        =...
    excludedata(t_phase(yesPD_i),lpdelta_i(yesPD_i),'indices',I(yesPD_i));
outlierspd(yesPD_i)             =   outliers;
retain_pd_i                     =   yesPD_i & ~outlierspd;
Apd_retain_pd_i                 =   Apd_i(retain_pd_i,:);
s_pd                            =   sum(Apd_retain_pd_i);
ray_crosses_pd_i                =   s_pd~=0;

end
function  [pab,lpdelta_i,retain_pd_i,ray_crosses_pd_i]...
                                =...
            Murat_retainPeakDealay(t_phase,l10pd_i,yesPD_i,Apd_i)
% Creates all constraints for peak delay inversion

fitrobust_i                     =...
    fit(t_phase,l10pd_i,'poly1','Robust','on');
pab                             =   [fitrobust_i.p1 fitrobust_i.p2];

l10pdt                          =   polyval(pab,t_phase);
lpdelta_i                       =   l10pd_i-l10pdt;

I                               =   abs(lpdelta_i) >= 2*std(lpdelta_i);
outlierspd                      =...
    excludedata(t_phase,lpdelta_i,'indices',I);

retain_pd_i                     =   yesPD_i & ~outlierspd;
Apd_retain_pd_i                 =   Apd_i(retain_pd_i,:);
s_pd                            =   sum(Apd_retain_pd_i);
ray_crosses_pd_i                =   s_pd~=0;

end
function peakDelay_i            =...
    Murat_peakDelay(sp_i,cursorPick_i,cursorPeakDelay_i,srate_i)
%CALCULATES peak delay time
pdsp                            =...
    sp_i(cursorPick_i:cursorPeakDelay_i);

[~,timeMaxAmplitude]            =	max(pdsp);

peakDelay_i                     =   timeMaxAmplitude/srate_i;


% % NEW FEATURE Takahashi 2009
% f                   = [3 6 12 18];
% 
% pdsp_ref            =...
%     Murat.temp.sp_ref(Murat.temp.cursor1:Murat.temp.pdcursor);
% 
% [~,tspm_ref]        = max(pdsp_ref);
% 
% % Reference peak delay at 6Hz (Takahashi et al. 2007)
% log_t_ref           = log(tspm_ref/tspm_ref(2));
% 
% [Af,Bf]             = polyfit(f,log_t_ref);
% Murat.temp.Af_i     = Af;
% Murat.temp.Bf_i     = Bf;

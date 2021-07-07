function peakDelay_i            =...
    Murat_peakDelay(sp_i,cursorPick_i,cursorPeakDelay_i,srate_i)

%CALCULATES peak delay time
lsp                             =   size(sp_i,2);
peakDelay_i                     =   zeros(lsp,1);
for i = 1:lsp
    pdsp                            =...
        sp_i(cursorPick_i:cursorPeakDelay_i,i);

    [~,timeMaxAmplitude]            =	max(pdsp);

    peakDelay_i(i)                     =   timeMaxAmplitude/srate_i;
end


% % NEW FEATURE Takahashi 2009
% cf                   = Murat.data.centralFrequency;
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

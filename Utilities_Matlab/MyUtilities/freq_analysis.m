function freq_analysis(path_s, path_w)
% 
% FREQ_ANALYSIS analyzes the waveform data in sac format returning
%               plot of the frequency content and arrival time of P and S.
%               It uses STA/LTA to exclude the noise before the analysis.
%               The same method is adopted to estimate the coda start which
%               is plotted as histogram.
%
% Input parameters:
%    path_s:         path of the folder for storing the output images
%    path_w:         path of the folder containing the waveforms
%
% Output parameters:
%    Frequency_distribution.png:  boxplot of the frequency content for all
%                                 the waveforms
%    Auto_Coda_Start.png:         histogram plot of the coda start 
%                                 (trigger-off values from STA/LTA analysis)
%    P_arrival_times.png:         histogram plot of the P-waves arrival
%                                 times (difference from the t0)
%    S_arrival_times:             histogram plot of the S-waves arrival
%                                 times (difference from the t0)
%    
%    Check_Plots folder:          random plot of waveforms
%
disp('Analyzing the frequency content with STA/LTA')

% List of the sac data
sacpath = fullfile(path_w, '*.sac');
sacfiles = dir(sacpath);

% STA/LTA parameters
sta_len = 1.0;
lta_len = 10.0;
thresh_on = 3.0;
thresh_off = 1.2;

% Spectra parameters
max_freq = 50;
freq_bins = 0:1:max_freq;
n_bins = length(freq_bins) - 1;
to_boxplot = zeros(length(sacfiles), n_bins);
p_arr = zeros(length(sacfiles),1);
s_arr = zeros(length(sacfiles),1);
t_trigger_on = NaN(length(sacfiles),1);
t_coda_start = NaN(length(sacfiles),1);

% Random plot settings
n_check_plots = 10;
indices_to_plot = randperm(length(sacfiles), min(n_check_plots, length(sacfiles)));
check_dir = fullfile(path_s, 'Check_Plots');
if ~exist(check_dir, 'dir')
    mkdir(check_dir);
end

% Iteration
for j = 1:length(sacfiles)
    % reading data
    sacfile = strcat(path_w,sacfiles(j).name);
    [~, data, hdr] = fread_sac(sacfile);
    fs = 1 / hdr.delta;
    p_arr(j) = hdr.a - hdr.o;
    s_arr(j) = hdr.t(1) - hdr.o;
    data = detrend(data);

    % STA/LTA
    cf = data.^2;
    n_sta = round(sta_len * fs);
    n_lta = round(lta_len * fs);
    if length(cf) < n_lta
        warning(['Waveform too short for LTA: ' sacfiles(j).name]);
        continue;
    end
    b_sta = ones(1, n_sta) / n_sta;
    b_lta = ones(1, n_lta) / n_lta;
    sta = filter(b_sta, 1, cf);
    lta = filter(b_lta, 1, cf);
    lta(lta == 0) = eps; 
    ratio = sta ./ lta;
    ratio(1:n_lta) = 0;
    idx_on = find(ratio > thresh_on, 1, 'first');
    if isempty(idx_on)
        warning(['No event detected in ', sacfiles(j).name]);
        continue;
    end
    idx_future = find(ratio(idx_on:end) < thresh_off, 1, 'first');
    if isempty(idx_future)
        idx_off = length(data);
    else
        idx_off = idx_on + idx_future - 1;
    end
    t_on_sec = (idx_on - 1) * hdr.delta;
    t_off_sec = (idx_off - 1) * hdr.delta;
    t_trigger_on(j) = t_on_sec;
    t_coda_start(j) = t_off_sec;

    % random plot
    if ismember(j, indices_to_plot)
        h_fig = figure('Visible', 'off');
        t_axis = (0:length(data)-1) * hdr.delta;
        subplot(2,1,1);
        plot(t_axis, data, 'k'); hold on;
        title(['Waveform: ' sacfiles(j).name], 'Interpreter', 'none');
        ylabel('Amplitude');
        grid on;
        t_start = 0;
        xline(hdr.o - t_start, 'b', 'T0', 'LineWidth', 1.5, 'LabelVerticalAlignment','bottom');
        xline(hdr.a - t_start, 'g', 'P (Head)', 'LineWidth', 1.5, 'LabelVerticalAlignment','bottom');
        xline(hdr.t(1) - t_start, 'r', 'S (Head)', 'LineWidth', 1.5, 'LabelVerticalAlignment','bottom');
        xline(t_on_sec, 'b--', 'Trig ON', 'LineWidth', 1.5);
        xline(t_off_sec, 'm--', 'Trig OFF', 'LineWidth', 1.5);
        legend('Data','T0','P','S','Trig ON','Trig OFF','Location','best');
        subplot(2,1,2);
        plot(t_axis, ratio, 'b'); hold on;
        yline(thresh_on, 'g--', 'Thresh ON');
        yline(thresh_off, 'r:', 'Thresh OFF');
        title('STA/LTA Ratio');
        xlabel('Time (s)');
        grid on;
        saveas(h_fig, fullfile(check_dir, ['Check_' sacfiles(j).name '.png']));
        close(h_fig);
    end

    % power spectra
    pad = round(1.0 * fs);
    i1 = max(1, idx_on - pad);
    i2 = min(length(data), idx_off + pad);
    data_cut = data(i1:i2);
    if length(data_cut) > fs
        [pw, fq] = pspectrum(data_cut,fs);
        for i = 1:n_bins
            idx_in_bin = fq >= freq_bins(i) & fq < freq_bins(i+1);
            if any(idx_in_bin)
                to_boxplot(j, i) = mean(pw(idx_in_bin));
            else
                to_boxplot(j, i) = NaN;
            end
        end
    end
end

disp('Saving results')
% Plot the frequency content
figure;
freq_centers = freq_bins(1:end-1) + 0.5;
to_boxplot = to_boxplot(~all(isnan(to_boxplot), 2), :);
boxplot(to_boxplot, 'Labels', string(freq_centers), 'Symbol','r.','OutlierSize',2);
xlabel('Frequency (Hz)');
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',6)
ylabel('Amplitude of the power spectra');
title('Box plot (1 Hz bin)');
yscale log
saveas(gcf, [path_s '/Frequency_distribution.png']);

% Plot the trigger off (coda onset)
t_coda_start = t_coda_start(~isnan(t_coda_start));
if ~isempty(t_coda_start)
    figure;
    histogram(t_coda_start, 20);
    xlabel('Time (s from start of file)');
    ylabel('Number of events');
    title(['Automatic Coda Start (LTA=' num2str(lta_len) 's, Thr=' num2str(thresh_off) ')']);
    grid on;
    saveas(gcf, fullfile(path_s, 'Auto_Coda_Start.png'));
end

% Plot the P arrival time
figure;
histogram(p_arr)
xlabel('Arrival time (s from t0)');
ylabel('Number of events');
title('Arrival time of P wave');
grid on;
saveas(gcf, [path_s '/P_arrival_times.png']);

% Plot the S arrival time
figure;
histogram(s_arr)
xlabel('Arrival time (s from t0)');
ylabel('Number of events');
title('Arrival time of S wave');
grid on;
saveas(gcf, [path_s '/S_arrival_times.png']);
end
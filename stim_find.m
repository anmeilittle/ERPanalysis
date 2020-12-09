function stim_times_cell = stim_find(dmat, Fs)
% function stim_find(dmat, Fs)
%
% Find the times at which stimulation occured.
%
% input: dmat - (N-channels x N-samples) - matrix of channel timeseries
%        Fs - sampling frequency
% David Huberdeau, ntb lab.

time = (1:length(dmat))/Fs;
% dmat_ = [zeros(size(dmat,1),1), diff(dmat, 1, 2)./(1/Fs)];
dmat_ = dmat;

% N_PC = 10;
WIN_TIME_LEN = .9;
WIN_SAMP_LEN = floor(WIN_TIME_LEN*Fs);
N_WIN = floor(size(dmat,2)/WIN_SAMP_LEN);
sig_inds = 1:size(dmat,2);
sig_inds = sig_inds(1:(end - mod(length(time), N_WIN)));

Z_TH_FREQ = 1;
Z_TH_TIME = 4.8;
% pmat = nan(N_WIN, N_PC, size(dmat,1)); %PC data
% fmat = nan(N_WIN, 2, size(dmat,1)); %factors about each window (time of onset and stim boolean)
% ps_diff_z = nan(N_WIN,size(dmat,1));
stim_times_cell = cell(1,size(dmat,1));
for i_ch = 1:size(dmat,1)
    dsignal_ = dmat_(i_ch,sig_inds);
    time_temp = time(sig_inds);

%     temp_winds = reshape(temp, length(sig_inds)/N_WIN, N_WIN)';
    signal_winds_ = reshape(dsignal_, length(sig_inds)/N_WIN, N_WIN)';
    signal_winds__ = signal_winds_ - repmat(mean(signal_winds_,2), 1, size(signal_winds_,2));
    time_winds = reshape(time_temp, length(sig_inds)/N_WIN, N_WIN)';
%     otrig_winds = reshape(odd_trig_temp, length(sig_inds)/N_WIN, N_WIN)';

%     [~, scr] = pca([temp_winds, temp_winds_]);
    [f_signal, p_signal] = simple_psd(signal_winds_(1,:), Fs);
    ps_winds = nan(N_WIN, length(p_signal));
    for i_wind = 1:N_WIN
        [f_signal, p_signal] = simple_psd(signal_winds__(i_wind,:), Fs);
        ps_winds(i_wind, :) = p_signal;
    end
    fs = f_signal;

    mps = mean(ps_winds,1);
    dps = ps_winds - repmat(mps, N_WIN, 1);
    ps_diff = sum(dps(:, fs > 150),2);
    ps_z = zscore(ps_diff);

    % do pca decomp for cleaner data:
    dps_ = dps(:, fs > 3);
    [coefs, scr] = pca(dps_);
    recon = scr(:, 1:2)*coefs(:, 1:2)';
    z_feat_recon = zscore(sum(recon(:, 50:end),2));

    % to be called stim, both must have ID'ed it. but this isn't enough..
    z_feat = (ps_z > Z_TH_FREQ) & (z_feat_recon > Z_TH_FREQ);
    inds = 1:length(z_feat);
    inds_th = inds(z_feat);

    % check within-window dynamics, and narrow down stim onset.
    stim_times = nan(length(inds_th), 1);
    inds_sig = 1:size(signal_winds__,2);
    for i_stim = 1:sum(z_feat)
        temp_sig = signal_winds__(inds_th(i_stim), :);
        time_sig = time_winds(inds_th(i_stim), :);
        temp_z = abs(zscore(temp_sig));
        if sum(temp_z > Z_TH_TIME) > 0
            % a true outlier is present in this window.. find it.
            inds_temp = inds_sig(temp_z > Z_TH_TIME);
            stim_times(i_stim) = time_sig(inds_temp(1));
        end
    end
    stim_times_cell{i_ch} = stim_times(~isnan(stim_times));
end

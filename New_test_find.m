hinfo = h5info('DC_stim_1_512hz.h5');
hdata = h5read('DC_stim_1_512hz.h5', '/data');
d2 = squeeze(hdata);
d2 = d2(1:193, :); 
% DC 193
% CH 137
% DR 233
% RW 81

load('labels_DC.mat')
chan = length(d2(:,1));

Fs = 512;
stim_wind = 10;

tt = (1/Fs):(1/Fs):(size(d2,2)/Fs);
inds = 1:length(tt);

%% finding times
times = stim_find(d2, 512); % stim_find func outputs the time points at which the stimulations were applied
time_sec = (1:length(d2))/512; % this array converts the frequency points into actual times

% combine and reduce all times
alltimes = times{1,1}.';
for i = 2:chan
    alltimes = sort(unique([alltimes times{1, i}.']));
end

% remove duplicates (within 4ms)
stim_times = 0;
%loop through all values
for i = 1:length(alltimes)-1
    if abs(alltimes(i+1) - alltimes(i)) <= .004
        continue
    else
        stim_times = [stim_times alltimes(i)];
    end
end
%remove the first 0
stim_times(1) = [];
stim_times(1) = [];
stim_times2 = stim_times';

% find the differences between stim points
time_dif = zeros(1, length(stim_times)-1);
for i = 1:length(stim_times)-1
    time_dif(i) = stim_times(i+1) - stim_times(i);
    if time_dif(i) > 0.999 && time_dif(i) < 1.001
        time_dif(i) = 1;
    else
    end    
end

% set threshold
n = 10; 
threshold = ones(1, n);
Hzstim = strfind(time_dif, threshold);

% extract each block start and end time
eachstart = Hzstim(1);
for i = 1:length(Hzstim)-1
    if Hzstim(i+1) - Hzstim(i) == 1
        continue
    else
        eachstart = [eachstart Hzstim(i)+10 Hzstim(i+1)];
    end
end

%if odd #, add last point
if mod(length(eachstart),2) == 1 
    eachstart = [eachstart Hzstim(length(Hzstim))+10];
else
end
%eachstart is an array of alternating start and end times

range = 0.7; %set range (no delay)
ERPmatrix3D = zeros(chan, round(range*512), length(eachstart)/2);

for i = 1:2:length(eachstart)
    [ERPmat] = ERP_func(range, chan, eachstart(i), eachstart(i+1), stim_times, d2);
    ERPmatrix3D(:,:,(i+1)/2) = ERPmat;
end

%save('ERPmatResults1_DC.mat','ERPmat');
%save('Eachstart_DC', 'eachstart');
 
%% find stim channels:
stim_channel_list = cell(2,length(stim_times2));
stim_channel_inds = nan(2,length(stim_times2));
stim_channel_mag = nan(2, length(stim_times2));
Everything = zeros(length(stim_times2), 3);

for i_stim = 1:length(stim_times2)
    [~, k_stim] = min((tt - stim_times2(i_stim)).^2);
    mag_list = nan(size(d2, 1), 1);
    for i_chan = 1:size(d2,1)
        sig_wind = d2(i_chan, k_stim + (-stim_wind:stim_wind));
        mag_list(i_chan) = max(abs(sig_wind));
    end

    [m_sort, k_sort] = sort(mag_list, 'descend');
    stim_channel_inds(:, i_stim) = k_sort(1:2);
    stim_channel_list{1, i_stim} = data.label{k_sort(1)};
    stim_channel_list{2, i_stim} = data.label{k_sort(2)};
    stim_channel_mag(:, i_stim) = m_sort(1:2);

    disp([data.label{k_sort(1)}, ' ', data.label{k_sort(2)}]);
    Everything(i_stim,:) = [stim_times2(i_stim) k_sort(1) k_sort(2)];
end

%save('StimLoc_DC, 'Everything');
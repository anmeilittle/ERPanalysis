% ERP script

d = h5read('/gpfs/milgram/project/turk-browne/projects/stimulation_behavior/intermediate_data/RW_stim_comb', '/data');
d2 = squeeze(d);
d2 = d2(1:81, :); 
% DC 193
% CH 137
% DR 233
% RW 81

load('RW_labels.mat')
chan = length(d2(:,1));

Fs = 512;
stim_wind = 10;

tt = (1/Fs):(1/Fs):(size(d2,2)/Fs);
inds = 1:length(tt);

%% finding times
times = stim_find(d2, 512); % stim_find func outputs the time points at which the stimulations were applied
time_sec = (1:length(d2))/512; % this array converts the frequency points into actual times

%find channel with most stimuations points
alltimes = zeros(1, chan);  % alltimes = times{1,1}.';
for i = 1:chan
    alltimes(i) = length(times{1, i});
end

%use that maximum # stim points as a dimension in stim_times matrix
maxlength = max(alltimes);
stim_times = zeros(chan, maxlength);
for i = 1:chan
    r = length(times{1,i});
    stim_times(i, 1:r) = times{1, i}.';
end

%find the differences between stim points
time_dif = zeros(chan, maxlength-1);
for i = 1:chan
    for j = 1:maxlength-1
        time_dif(i, j) = stim_times(i, j+1) - stim_times(i, j);
        if time_dif(i,j) > 0.999 && time_dif(i,j) < 1.001
            time_dif(i,j) = 1;
        else
        end
    end
end

% set threshold
n = 10; 
threshold = ones(1, n);
Hzstim = zeros(chan, maxlength);
% find where there are sequences of n 1s intervals
for i = 1:chan
    %get original stim_times not just the indices
    Hzstim(i, 1:length(strfind(time_dif(i,:), threshold))) = stim_times(i, strfind(time_dif(i,:), threshold));
end
Hzstim = Hzstim(:, any(Hzstim)); %remove empty columms

% we now have a matrix with all the stim times detected from each of the
% channels
% combine and reduce this matrix

alltimes = Hzstim(1,:);
for i = 2:chan
    alltimes = sort(unique([alltimes Hzstim(i,:)]));
end
if alltimes(1) == 0
    alltimes(1) = [];
else
end

%remove duplicates
alltimes2 = 0;
%loop through all values
for i = 1:length(alltimes)-1
    if abs(alltimes(i+1) - alltimes(i)) <= .5 %.004 %idk if this is valid
        continue
    else
        alltimes2 = [alltimes2 alltimes(i)];
    end
end
%remove the first 0
alltimes2(1) = [];

% get start and stop times for each block
eachstart = 1;
for i = 2:length(alltimes2)-1
    if (alltimes2(i+1) - alltimes2(i)) < 1.001 && (alltimes2(i) - alltimes2(i-1)) < 1.001 %within block
        continue
    elseif (alltimes2(i+1) - alltimes2(i)) > 1.001 && (alltimes2(i) - alltimes2(i-1)) > 1.001 %isolated block with length n
        eachstart = [eachstart i i];
    elseif (alltimes2(i) - alltimes2(i-1)) > 1.001 && (alltimes2(i+1) - alltimes2(i)) < 1.001 %start of block
        eachstart = [eachstart i]; 
    elseif (alltimes2(i) - alltimes2(i-1)) < 1.001 && (alltimes2(i+1) - alltimes2(i)) > 1.001 %end of block
        eachstart = [eachstart i]; 
    end
end
%if odd #, add last point
if mod(length(eachstart),2) == 1 
    eachstart = [eachstart length(alltimes2)];
else
end
%eachstart is an array of alternating start and end times

save('eachstart4.mat', 'eachstart');

range = 0.7; %set range (no delay)
ERPmatrix3D = zeros(chan, round(range*512), length(eachstart)/2);

for i = 1:2:length(eachstart)-1
    [ERPmat] = ERP_func(range, chan, eachstart(i), eachstart(i+1), alltimes2, d2);
    ERPmatrix3D(:,:,(i+1)/2) = ERPmat;
end

save('ERPmatResults4.mat','ERPmatrix3D');

%% find stim channels:
stim_channel_list = cell(2,length(alltimes2));
stim_channel_inds = nan(2,length(alltimes2));
stim_channel_mag = nan(2, length(alltimes2));
Everything = zeros(length(alltimes2), 3);

for i_stim = 1:length(alltimes2)
    [~, k_stim] = min((tt - alltimes2(i_stim)).^2);
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
    Everything(i_stim,:) = [alltimes2(i_stim) k_sort(1) k_sort(2)];
end

save('StimLoc4', 'Everything');

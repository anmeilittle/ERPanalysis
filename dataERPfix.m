% ERP script

% % load data
d = h5read('/gpfs/milgram/project/turk-browne/projects/stimulation_behavior/intermediate_data/DC_stim_comb', '/data');
d2 = squeeze(d);
chan = length(d2(:,1));

% % finding times
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
n = 20; 
threshold = ones(1, n);
Hzstim = zeros(chan, maxlength);
% find where there are sequences of 1s intervals
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
alltimes2 = 0
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

%so there are still stim times that seem to be overlapping
eachstart = alltimes2(1);
for i = 1:length(alltimes2)-1
    if alltimes2(i+1) - eachstart(end) < n
        continue
    else
        eachstart = [eachstart alltimes2(i+1)];
    end
end

range = 0.7; %set range (no delay)
ERPmatrix3D = zeros(chan, round(range*512), length(eachstart));


for i = 1:length(eachstart)
    [ERPmat] = ERP_func(range, chan, eachstart(i), eachstart(i)+n, stim_times, d2);
    ERPmatrix3D(:,:,i) = ERPmat;
end

save('ERPmatResults1.mat','ERPmatrix3D');

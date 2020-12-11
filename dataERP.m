
% ERP script

% % load data
d = h5read('/gpfs/milgram/project/turk-browne/projects/stimulation_behavior/intermediate_data/CH_stim_comb', '/data');
d2 = squeeze(d);
chan = length(d2(:,1));

% % finding times
times = stim_find(d2, 512); % stim_find func outputs the time points at which the stimulations were applied
time_sec = (1:length(d2))/512; % this array converts the frequency points into actual times

% combine and reduce all times
alltimes = times{1,1}.';
for i = 2:chan
    alltimes = sort(unique([alltimes times{1, i}.']));
end

%remove duplicates (within 4ms)
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

%find the differences between stim points
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
    ERPmatrix3D(:,:,i) = ERPmat;
end

save('ERPmatResults2','ERPmatrix3D');

% %plot stim times with a channel to check rounding
% startt = 606;
% stopt = 626;
% 
% timePlot = stim_times(startt)*512:stim_times(stopt)*512;
% figure;
% tiledlayout(2, 2);
% for i = [37,38,99,100]
%      nexttile;
%      plot(timePlot, d2(i, timePlot));
%      for j = startt:stopt
%         xline(stim_times(j)*512,'--r');
%      end
%      title(['Channel ' num2str(i)]);  
% end

 
% %single channel
% figure;
% plot(timePlot, d2(96, timePlot));
%      for j = startt:stopt
%         xline(stim_times(j)*512,'--r');
%      end
% 



function [ERPmat] = ERP_func(range, chan, stim_start, stim_end, stim_times, d2)
% this function creates a matrix displaying the average responses of each
% channel to a repeated stimulation in a given block

% points in stim_times that are 1s apart in one block
block = stim_start:stim_end;
for i = 1:10
    block = [block block(end)+1];
end

% for each channel, we want a .9sec (range) ERP
ERPmat = zeros(chan, round(range*512));

for i = 1:chan %loop over each channel
    
    %temporary matrix for all the stimulations within the block (rows) and
    %all the responses from a given channel following each stimulation
    temp_mat = zeros(length(block), round(range*512));
    
    % iterate through each row (stimulation) within the block
    for j = 1:length(block)
        % average the responses --> ERP
        ERPtime = [round(block(j)*512):round((block(j)+range)*512)-1];
        temp_mat(j, :) = d2(i, ERPtime);  
    end
 
    %get average of each column (time) and insert into ERPmatrix
    ERPmat(i, :) = mean(temp_mat);
   
end

% 
% figure;
% tiledlayout(6, 6);
% for i = 1:36
%      nexttile;
%      plot(1:round(range*512), ERPmat(i,:));
%      title(['Channel ' num2str(i)]);
%      ylim([-300 300]);  
% end

% figure;
% tiledlayout(6, 6);
% for i = 37:72
%      nexttile;
%      plot(1:round(range*512), ERPmat(i,:));
%      title(['Channel ' num2str(i)]);
%      ylim([-350 350]);  
% end
% 
% figure;
% tiledlayout(6, 6);
% for i = 72:107
%      nexttile;
%      plot(1:round(range*512), ERPmat(i,:));
%      title(['Channel ' num2str(i)]);
%      ylim([-350 350]);  
% end
% 
% figure;
% tiledlayout(6, 6);
% for i = 108:143
%      nexttile;
%      plot(1:round(range*512), ERPmat(i,:));
%      title(['Channel ' num2str(i)]);
%      ylim([-350 350]);  
% end
% 
% figure;
% tiledlayout(6, 6);
% for i = 144:179
%      nexttile;
%      plot(1:round(range*512), ERPmat(i,:));
%      title(['Channel ' num2str(i)]);
%      ylim([-350 350]);  
% end
% 
% figure;
% tiledlayout(6, 6);
% for i = 180:203
%      nexttile;
%      plot(1:round(range*512), ERPmat(i,:));
%      title(['Channel ' num2str(i)]);
%      ylim([-350 350]);  
% end
%     

end


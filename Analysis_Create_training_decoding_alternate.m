%% script from Eric and modified by jisoo
% to alternate training/decoding session.

function [training_ts] =Analysis_Create_training_decoding_alternate(ca_time,len)
data = ca_time';

%% get indicies in alternating blocks of a specific size.


true_arr = ones(1,len);
false_arr = true_arr*0;

keep_idx = logical(repmat([true_arr, false_arr], 1,ceil(size(data,2)/length([true_arr, false_arr]))));
keep_idx = keep_idx(1:size(data,2)); % cut off any outliers.
training_ts=keep_idx;

%% grab a random set of data points.
% [~, idx] = datasample(data, floor(size(data, 2)/2));
% keep_idx2 = ones(size(data));
% 
% keep_idx2(idx) = 0;
% keep_idx2 = logical(keep_idx2);

% % plot to check
% figure(101)
% hold on
% plot(1:length(data), data, 'k')
% plot(1:length(data), keep_idx+2, 'r')
% plot(1:length(data), keep_idx2-2, 'b')
% legend({'data', 'block idx', 'rand sample'})

end
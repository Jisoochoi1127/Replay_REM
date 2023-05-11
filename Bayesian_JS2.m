
 

function [decoding]= Bayesian_JS2(PARAMS, ms, behav, all_binary_pre_REM, all_binary_pre_SW, all_binary_post_SW, all_binary_post_REM)
%% Loading input

% clear;
% close all;

%Loading track data and convert to each variable
% load 'ms'

% load 'behav'
 ca_time= ms.time/1000;
 ca_data=ms.RawTraces ;
 behav_time=behav.time/1000;
 behav_vec=behav.position(:,1);
% 
% cell_used = logical(ones(size(ca_data,2),1)); % Use every cell

% % load 'behav'
ca_time= ms.time/1000;
ca_data=ms.RawTraces ;
behav_time=behav.time/1000;
behav_vec=behav.position(:,1);

 cell_used = logical(ones(size(ca_data,2),1)); % Use every cell

%Loading sleep data
% load 'all_binary_post_REM'
% load 'all_binary_post_SW'
% load 'all_binary_pre_REM'
% load 'all_binary_pre_SW'

%%  Binarizing calcium data during Track
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise

binarized_data = zeros(size(ca_data));
for cell_i = 1:size(ca_data,2)
    binarized_data(:,cell_i) = extract_binary(ca_data(:,cell_i), sampling_frequency, z_threshold);
end

%% Interpolate behav and calcium

[unique_behavtime, IAbehav, ICbehav]=unique(behav_time);
[unique_mstime, IAms, ICms]=unique(ca_time);
interp_behav_vec = interp1(IAbehav,behav_vec(IAbehav,1),IAms);
interp_behav_vec(end) = interp_behav_vec(end-1);

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[velocity] = extract_velocity(interp_behav_vec, ca_time);

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
min_speed_threshold = 5; % 2 cm.s-1
running_ts = velocity > min_speed_threshold;

%% Compute occupancy and joint probabilities
%bin_size = 3;
% Make sure that your binning vector includes every data point of
% interp_behav_vec using the min/max function:
bin_vector = min(interp_behav_vec):bin_size:max(interp_behav_vec)+bin_size; % start : bin_size : end
bin_centers_vector = bin_vector + bin_size/2;
bin_centers_vector(end) = [];

%% Decoding
% First let's binarize traces from all cells
binarized_data = zeros(size(ca_data));
for cell_i = 1:size(ca_data,2)
    binarized_data(:,cell_i) = extract_binary(ca_data(:,cell_i), sampling_frequency, z_threshold);
end

 %% Select random frames to train decoder on 50% of wake data
    training_ts = create_training_set(ca_time,'random',0.5);
training_ts(running_ts == 0) = 0;    
    %% Create tuning curves for every cell
    for cell_i = 1:size(binarized_data,2)
        [~, ~, occupancy_vector, prob_being_active(cell_i), tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, bin_vector, training_ts);
    end
    
    Occupancy_bin=occupancy_vector;
    
    %% Plot the tunning curves
[~,max_index] = max(tuning_curve_data,[],1);
[~,sorted_index] = sort(max_index);
sorted_tuning_curve_data = tuning_curve_data(:,sorted_index);

figure
imagesc(bin_centers_vector,1:size(ca_data,2),sorted_tuning_curve_data')
daspect([1 1 1])
caxis([0 0.3])
title 'Neuronal tuning curves'
xlabel 'Position on the track (cm)'
ylabel 'Cell ID'

    %% Decode position
    
    % First, let us establish the timestamps used for decoding.
    decoding_ts = ~training_ts; % Training timestamps are excluded
    decoding_ts(running_ts == 0) = 0; % Periods of immobility are excluded
    
    % Minimal a priori (use to remove experimental a priori)
    occupancy_vector = occupancy_vector./occupancy_vector*(1/length(occupancy_vector));
    
   
    
%% Extract  decoded position during Track

[WAKE_decoded_probabilities] = bayesian_decode1D(binarized_data, occupancy_vector, prob_being_active, tuning_curve_data, cell_used);

[max_decoded_prob, decoded_bin] = max(WAKE_decoded_probabilities,[],1);
decoded_position = bin_centers_vector(decoded_bin);

%added by jisoo
decoded_bin(isnan(max_decoded_prob)) = nan;
decoded_position(isnan(decoded_bin)) = nan;

whole_Wake_decoded_position=decoded_position;
whole_Wake_decoded_probabilities=WAKE_decoded_probabilities;


% Before looking at the error rate, we must first bin the actual data using the same bin vector used by
% the decoder
actual_bin = nan*interp_behav_vec;
actual_position = nan*interp_behav_vec;
for bin_i = 1:length(bin_vector)-1
    position_idx = find(interp_behav_vec>bin_vector(bin_i) & interp_behav_vec < bin_vector(bin_i+1));
    actual_bin(position_idx) = bin_i;
    actual_position(position_idx) = bin_centers_vector(bin_i);
end

%added by jisoo_calculate velcodity for each bin
for i=1:length(bin_vector)-1
    A=find(actual_bin==i);
    Velocity_bin(i)=mean(velocity(A));
end
Velocity_bin=Velocity_bin';

%% plotting 
figure
subplot(3,1,1)
imagesc(ca_time,bin_centers_vector,WAKE_decoded_probabilities)
title 'Posterior probabilities'
xlabel 'Time (s)'
ylabel 'Position on the track (cm)'
ax1 = gca;
colormap(hot)
% colorbar
caxis([0.02,0.08]);
%ax1.CLim = [0 0.1];
ax1.XLim = [515 521];
ax1.YDir = 'normal';

subplot(3,1,2)
scatter(ca_time,actual_position,'o')
hold on
scatter(ca_time, decoded_position,'*')
title 'Actual versus decoded position'
xlabel 'Time (s)'
ylabel 'Location on the track (cm)'
ax2 = gca;
ax2.XLim = [515 521]; % Let's plot a single trajectory
linkaxes([ax1 ax2], 'x')
legend('actual','decoded');

subplot(3,1,3)
plot(ca_time,velocity,'LineWidth',1)
title 'Speed'
xlabel 'Time (s)'
ylabel 'Velocity'
ax2 = gca;
ax2.XLim = [515 521]; % Let's plot a single trajectory
linkaxes([ax1 ax2], 'x')

% subplot(4,1,4)
% plot(ca_time,sum(binarized_data,2))
% title 'Sum_binarized'
% xlabel 'Time (s)'
% ylabel 'Total binary'
% ax2 = gca;
% ax2.XLim = [447 452]; % Let's plot a single trajectory
% linkaxes([ax1 ax2], 'x')

%% Remove training timestamps to assess decoding error rate
decoded_bin(~decoding_ts) = nan;

decoded_position(~decoding_ts) = nan;
%decoded_probabilities(:,~decoding_ts) = nan; % this is wrong.
WAKE_decoded_probabilities(:,~decoding_ts) = nan;

actual_bin(~decoding_ts) = nan;
actual_position(~decoding_ts) = nan;
actual_bin = actual_bin';
actual_position =  actual_position';




%% Compute decoding agreement
decoding_agreement_vector = double(decoded_bin == actual_bin);
decoding_agreement_vector(isnan(decoded_bin)) = nan;
decoding_agreement_vector(isnan(actual_bin)) = nan;
decoding_agreement_vector(isnan(decoding_agreement_vector)) = [];
decoding_agreement = sum(decoding_agreement_vector)./length(decoding_agreement_vector);

%% Compute decoding error
decoding_error = actual_position - decoded_position;
abs_decoding_error=abs(decoding_error);
mean_decoding_error = mean(abs(decoding_error), 'omitnan');

%% Compute confusion matrix
confusion_matrix = zeros(length(bin_centers_vector),length(bin_centers_vector));

for actual_i = 1:length(bin_centers_vector)
   for decoded_i = 1:length(bin_centers_vector)
       confusion_matrix(actual_i,decoded_i) = sum(decoded_bin == decoded_i & actual_bin == actual_i)./sum(actual_bin == actual_i);
   end
end

% Plot the confusion matrix
figure
imagesc(bin_centers_vector, bin_centers_vector, confusion_matrix)
set(gca,'YDir','normal') 
colormap hot
colorbar
caxis([0 1])
title 'Confusion matrix'
xlabel 'Actual position (cm)'
ylabel 'Decoded position (cm)'

Confusion_max=max(confusion_matrix)';
figure
plot(Confusion_max);




%% Extract decoded position during post REM sleep

    [REM_decoded_probabilities] = bayesian_decode1D(all_binary_post_REM, occupancy_vector, prob_being_active, tuning_curve_data, cell_used);
    
    % Maximum a posteriori
    [REM_max_prob, REM_decoded_bin] = max(REM_decoded_probabilities,[],1);
    REM_decoded_position = bin_centers_vector(REM_decoded_bin);
    
    REM_decoded_bin(isnan(REM_max_prob)) = nan;
    % added by jisoo _ nan value for position as well

    REM_decoded_position(isnan(REM_decoded_bin)) = nan;


%% Extract decoded position during post SWS sleep

[SWS_decoded_probabilities] = bayesian_decode1D(all_binary_post_SW, occupancy_vector, prob_being_active, tuning_curve_data, cell_used);

% Maximum a posteriori
[SWS_max_prob, SWS_decoded_bin] = max(SWS_decoded_probabilities,[],1);
SWS_decoded_position = bin_centers_vector(SWS_decoded_bin);

SWS_decoded_bin(isnan(SWS_max_prob)) = nan;
%added by jisoo

SWS_decoded_position(isnan(SWS_decoded_bin)) = nan;

  
%% Extract decoded position during pre REM sleep

[pre_REM_decoded_probabilities] = bayesian_decode1D(all_binary_pre_REM, occupancy_vector, prob_being_active, tuning_curve_data, cell_used);

% Maximum a posteriori
[pre_REM_max_prob, pre_REM_decoded_bin] = max(pre_REM_decoded_probabilities,[],1);
pre_REM_decoded_position = bin_centers_vector(pre_REM_decoded_bin);

pre_REM_decoded_bin(isnan(pre_REM_max_prob)) = nan;
%added by jisoo
pre_REM_decoded_position(isnan(pre_REM_decoded_bin)) = nan;



%% Extract decoded position during pre SWS sleep

[pre_SWS_decoded_probabilities] = bayesian_decode1D(all_binary_pre_SW, occupancy_vector, prob_being_active, tuning_curve_data, cell_used);

% Maximum a posteriori
[pre_SWS_max_prob, pre_SWS_decoded_bin] = max(pre_SWS_decoded_probabilities,[],1);
pre_SWS_decoded_position = bin_centers_vector(pre_SWS_decoded_bin);

pre_SWS_decoded_bin(isnan(pre_SWS_max_prob)) = nan;
%added by jisoo

pre_SWS_decoded_position(isnan(pre_SWS_decoded_bin)) = nan;

 
  %% Plotting decoded position during REM
    
figure;
subplot(2,1,1)
x1=1:size(all_binary_post_REM,1);

imagesc(x1,bin_centers_vector,REM_decoded_probabilities)
title 'Posterior probabilities during REM'
xlabel 'Time (s)'
ylabel 'Position on the track (cm)'
ax1 = gca;
ax1.CLim = [0 0.1];
%ax1.XLim = [447 452];
%ax1.XLim = [480 484]; % Modified by Jisoo
subplot(2,1,2)
plot(x1,sum(all_binary_post_REM,2))
title 'Sum_binarized'
xlabel 'Time (s)'
ylabel 'Total binary'
ax2 = gca;
linkaxes([ax1 ax2], 'x')
    
%% Plotting decoded position during SWS
    
figure;
subplot(2,1,1)
x2=1:size(all_binary_post_SW,1);

imagesc(x2,bin_centers_vector,SWS_decoded_probabilities)
title 'Posterior probabilities during SWS'
xlabel 'Time (s)'
ylabel 'Position on the track (cm)'
ax1 = gca;
ax1.CLim = [0 0.1];

subplot(2,1,2)
plot(x2,sum(all_binary_post_SW,2))
title 'Sum_binarized'
xlabel 'Time (s)'
ylabel 'Total binary'
ax2 = gca;
linkaxes([ax1 ax2], 'x')

%% Plotting decoded position during pre REM
    
figure;
subplot(2,1,1)
x3=1:size(all_binary_pre_REM,1);
imagesc(x3,bin_centers_vector,pre_REM_decoded_probabilities)
title 'Posterior probabilities during pre REM'
xlabel 'Time (s)'
ylabel 'Position on the track (cm)'
ax1 = gca;
ax1.CLim = [0 0.1];

subplot(2,1,2)
plot(x3,sum(all_binary_pre_REM,2))
title 'Sum_binarized'
xlabel 'Time (s)'
ylabel 'Total binary'
ax2 = gca;
linkaxes([ax1 ax2], 'x')

%% Plotting decoded position during pre SWS
    
figure;
subplot(2,1,1)
x4=1:size(all_binary_pre_SW,1);
imagesc(x4,bin_centers_vector,pre_SWS_decoded_probabilities)
title 'Posterior probabilities during pre SWS'
xlabel 'Time (s)'
ylabel 'Position on the track (cm)'
ax1 = gca;
ax1.CLim = [0 0.1];

subplot(2,1,2)
plot(x4,sum(all_binary_pre_SW,2))
title 'Sum_binarized'
xlabel 'Time (s)'
ylabel 'Total binary'
ax2 = gca;
linkaxes([ax1 ax2], 'x')

%% Save the output of decoded results
decoding.bin_centers_vector=bin_centers_vector;
decoding.WAKE_decoded_probabilities =whole_Wake_decoded_probabilities;
decoding.WAKE_decoded_position=whole_Wake_decoded_position;

decoding.Testing_wake_decoded_probabilities=WAKE_decoded_probabilities;
decoding.Testing_wake_decoded_position=decoded_position;


decoding.REM_decoded_probabilities=REM_decoded_probabilities;
decoding.REM_max_prob=REM_max_prob;
decoding.REM_decoded_position=REM_decoded_position;

decoding.SWS_decoded_probabilities=SWS_decoded_probabilities;
decoding.SWS_max_prob=SWS_max_prob;
decoding.SWS_decoded_position=SWS_decoded_position;

decoding.pre_REM_decoded_probabilities=pre_REM_decoded_probabilities;
decoding.pre_REM_max_prob=pre_REM_max_prob;
decoding.pre_REM_decoded_position=pre_REM_decoded_position;

decoding.pre_SWS_decoded_probabilities=pre_SWS_decoded_probabilities;
decoding.pre_SWS_max_prob=pre_SWS_max_prob;
decoding.pre_SWS_decoded_position=pre_SWS_decoded_position;
decoding.decoding_error=decoding_error;
decoding.abs_decoding_error=abs_decoding_error;
    

%save('decoding.mat','decoding')
%save('decoding_training90.mat','decoding')

 end



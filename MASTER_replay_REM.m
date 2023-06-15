%%% MASTER_REM_replay

%% Jisoo to do list

% for master scrip
% DONE_ collect data first as ms is needed for params
% DONE_ Set up parameters -- May 11th 
% DONE_need to add line(line 45) for passing some sessions doens't have sleep data
% DONE_ add Plotting part from Bayesian_JS script
% DONE_ add replay part
% DONE_Run data session for whole session 5/26
% Add Place cell script - pv1254_HATD1


% for Bayesian_JS2 script
% DONE_ change Bayesian JS2 variable using params -- May 11th 
% DONE_ add blocking method to get training_set 
% DONE_ Separate plotting section from Bayesian JS2
% DONE_ Remove all_binary_pre_SW,all_binary_post_SW in the analysis
% DONE_ add shuffling error
% DONE_ save decoding as id, data
% DONE_ check shuffling error with GE
% DONE_ check shuffling error line by line 5/26

% DONE_ Check [Shuffled_decoded_probabilities] = bayesian_decode1D(Shuffled_binarized_data, occupancy_vector_S, prob_being_active_S, tuning_curve_data_shuffled, PARAMS.decoding.cell_used);
% DONE_ add half left, right decoding error 5/26
% DONE_ Remove plotting part from Bayesian_JS2-confusion matrix. etc
% DONE_add decoding.tuningcurve
%Check rng seed - for bayesianjs2
% add decoding_parameters..? different length., cell number.
% half left, right error with matching frames _ HAT..

% for Replay script
% DONE_step size, window size, threshold etc put PARAMS.
% DONE_Check script with decoding data
% DONE_Replay- info is not loaded. 
% separate the selection part from replay - from line 227
% save all the variables from replay script for plotting figures
% sampling_percentage_time ; less than threshold... - that cause error
% jumpiness - parameter.
% remove plotting from replay script.

% etc
% Add place cell script before decoding
% Add theta amplitude/frequency during replay script
% Add Bias index script 
% Check the selction of replay _by slope/ jumpiness
%% something to discuss



%% set path for data
data_dir = '/Users/jisoo/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/Data';
inter_dir = '/Users/jisoo/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/inter/inter_data';
mkdir(inter_dir);
decoding_dir = '/Users/jisoo/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/inter/decoding';
mkdir(decoding_dir);
replay_dir = '/Users/jisoo/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/inter/replay';
mkdir(replay_dir);


%% collect data and generate intermediate files.

cd(data_dir)
sub_list = dir('pv*');

for iSub = length(sub_list):-1:1
    cd([data_dir filesep sub_list(iSub).name]);
    
    
    LT_list = dir('*LTD*');
    HAT_list= dir('*HATD*');
    sess_list=[LT_list; HAT_list];
    
    
    for iS  = length(sess_list):-1:1
        
        cd([data_dir filesep sub_list(iSub).name filesep sess_list(iS).name]);
        
       
        warning off
        load('ms.mat', 'ms')
        load('behav.mat')
        load('all_binary_pre_REM.mat')
        load('all_binary_post_REM.mat')
        
        
        info = [];
        info.subject = sub_list(iSub).name;
        info.session = sess_list(iS).name;
        info.nCells = ms.numNeurons;
        info.dur_preREM = size(all_binary_pre_REM,1)/30/60;%min
        info.dur_postREM = size(all_binary_post_REM,1)/30/60;%min
        
        
        fprintf('Saving %s   %s ...\n', info.subject, info.session)
        save([inter_dir filesep info.subject '_' info.session '_data.mat'], 'ms', 'behav', 'all_binary_pre_REM', 'all_binary_post_REM', 'info')
        warning on
    end
end

%% PARAMS

% set RNG for consistency
PARAMS.rng = 1234;

% set the parameters for decoding
PARAMS.decoding.bin_size = 3;
PARAMS.decoding.cell_used = 'Whole_cell'; %Whole_cell, Place_cell, Non_Place_cell
PARAMS.decoding.sampling_frequency = 30; % This data set has been sampled at 30 images per second
PARAMS.decoding.z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
PARAMS.decoding.min_speed_threshold = 5; % 2 cm.s-1
PARAMS.decoding.training_set_creation_method = 'Block'; 
%'quarter_portion'%'Non_overlapped';%divided ; quarter%'random'; % 'odd', odd timestamps; 'first_portion', first portion of the recording; 3, 'random' random frames
PARAMS.decoding.training_set_portion = 1; % Portion of the recording used to train the decoder for method 2 and 3
PARAMS.decoding.len = 10; % length of block frames(alternating training/decoding for this length of block)
PARAMS.decoding.numshuffles = 100; % calculate average decoding error with 1000 shuffled data

% set the parameters for replay
PARAMS.replay.step_size=10; %step_size ; overlapped portion. 
PARAMS.replay.windowsize=29; %windowsize=duration of time window.
PARAMS.replay.numshuffles = 100; 
PARAMS.replay.sampling_threshold=0.7;%how many data points for each time window
%PARAMS.replay.slope
PARAMS.replay.jumpiness=50;

% style parameters for figures

%% Decoding
cd(inter_dir);

fnames = dir('*_data.mat');

for iF = 1:length(fnames)
    warning off
    load(fnames(iF).name)
    warning on
    fprintf('Generating TCs for: %s   %s....', info.subject, info.session)
    tic
    
    % set the parameters for data
    PARAMS.data.ca_time= ms.time/1000;
    PARAMS.data.ca_data=ms.RawTraces ;
    PARAMS.data.behav_time=behav.time/1000;
    PARAMS.data.behav_vec=behav.position(:,1);
    
    [decoding] = Bayesian_JS2(decoding_dir, info, PARAMS, ms, behav, all_binary_pre_REM, all_binary_post_REM);
    toc
    
    % where and how do we save this decoding?
    save([decoding_dir filesep info.subject '_' info.session '_decoding.mat'], 'decoding')
    
    
    clear decoding ms behav info all_binary_pre_REM all_binary_pre_SW all_binary_post_REM all_binary_post_SW
    
    
    fprintf('done\n')
end


% function for generating TCs goes here. Give it the Inter_dir as an input
% along with the paramters and let it go.


%% Replay

cd(decoding_dir);
fnames = dir('*_decoding.mat');
for iF = 1:length(fnames)
    warning off
    load(fnames(iF).name)
    warning on
    fprintf('Generating TCs for: %s   %s....', decoding.info.subject, decoding.info.session)
    tic
   
    [out] = Replay_Fidelity_linear_regression(replay_dir, PARAMS,decoding);
    toc
    
    save([replay_dir filesep decoding.info.subject '_' decoding.info.session '_replay.mat'], 'out')
        
    
    fprintf('done\n')
end




%% Figures

% Plot the tunning curves
[~,max_index] = max(decoding.tuning_curve_data,[],1);
[~,sorted_index] = sort(max_index);
sorted_tuning_curve_data = decoding.tuning_curve_data(:,sorted_index);

figure
imagesc(bin_centers_vector,1:size(PARAMS.data.ca_data,2),sorted_tuning_curve_data')
daspect([1 1 1])
caxis([0 0.3])
title 'Neuronal tuning curves'
xlabel 'Position on the track (cm)'
ylabel 'Cell ID'

% Plot the confusion matrix
figure
imagesc(bin_centers_vector, bin_centers_vector, decoding.confusion_matrix)
set(gca,'YDir','normal') 
colormap hot
colorbar
caxis([0 1])
title 'Confusion matrix'
xlabel 'Actual position (cm)'
ylabel 'Decoded position (cm)'

Confusion_max=max(decoding.confusion_matrix)';
figure
plot(Confusion_max);

% plotting running (actual position and decoded position)
figure
subplot(3,1,1)
imagesc(PARAMS.data.ca_time,bin_centers_vector,WAKE_decoded_probabilities)
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
scatter(PARAMS.data.ca_time,actual_position,'o')
hold on
scatter(PARAMS.data.ca_time, decoded_position,'*')
title 'Actual versus decoded position'
xlabel 'Time (s)'
ylabel 'Location on the track (cm)'
ax2 = gca;
ax2.XLim = [515 521]; % Let's plot a single trajectory
linkaxes([ax1 ax2], 'x')
legend('actual','decoded');

subplot(3,1,3)
plot(PARAMS.data.ca_time,velocity,'LineWidth',1)
title 'Speed'
xlabel 'Time (s)'
ylabel 'Velocity'
ax2 = gca;
ax2.XLim = [515 521]; % Let's plot a single trajectory
linkaxes([ax1 ax2], 'x')



% Plotting decoded position during post-REM
    
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
    

% Plotting decoded position during pre-REM
    
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


% plotting for replay
figure;
scatter(sampling_percentage,Replay_score_actual,'*');
ylim([0 1]);
title('scatter plot of sampling percentage and replay score');
xlabel('Sampling Percentage');
ylabel('Replay score')
xline(0.7,'r','LineWidth',2);

%plotting slope histogram
figure;
histogram(slope, 'EdgeAlpha',0);
xlabel('Slope')
ylabel('Frequency')
box off




% script for Fig 1


% script for fig 2




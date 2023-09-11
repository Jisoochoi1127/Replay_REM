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
% DONE_separate the selection part from replay - from line 227
% DONE_save all the variables from replay script for plotting figures
% DONE_check params is loaded???- no ; should save in decoding.PARAMS
% DONE_ treaming script - save only necessary variables.
% check jumpiness - what happens with nan (mas(diff(abs..)
%check decoding.PARAMS woudl be ideal? 
% check actual jumpiness part (newly added)

% for slection script
% DONE_Include selection script in the master Monday
% DONE_include slope threshold
% DONE_remove plotting from replay script.
% Run selction script through master script
% sampling_percentage_time ; less than threshold... - that cause error

% etc
% Add place cell script before decoding
% Add theta amplitude/frequency during replay script
% Add Bias index script 
% Check the selction of replay _by slope/ jumpiness

% add plotting repaly script
%2023.09.08
% DONE_modify decoding binning - using absolute value
% DONE_place cell script - change training_ts to whole frames
% DONE cell_used ; include place cell only, non place cell only
% DONE include place cell output in decoding part - line 194
%Decoding error in hat- match sampling between closed and open arms

%2023.09.11
% check fnames, fnames_PC works in the master script
% Bayesian_JS2 - check if place cell, non place cells _ cell used works


%% set path for data
data_dir = '/Users/jisoo/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/Data';
inter_dir = '/Users/jisoo/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/inter/inter_data';
mkdir(inter_dir);
decoding_dir = '/Users/jisoo/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/inter/decoding';
mkdir(decoding_dir);
replay_dir = '/Users/jisoo/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/inter/replay';
mkdir(replay_dir);
selection_dir='/Users/jisoo/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/inter/selection_replay';
mkdir(selection_dir);


% Path for Comp Canada

% data_dir = '/lustre06/project/6064766/datasets/Jisoo/data';
% inter_dir = '/lustre06/project/6064766/datasets/Jisoo/inter';
% mkdir(inter_dir);
% decoding_dir = '/lustre06/project/6064766/datasets/Jisoo/inter/decoding';
% mkdir(decoding_dir);
% replay_dir = '/lustre06/project/6064766/datasets/Jisoo/inter/replay';
% mkdir(replay_dir);
% selection_dir = '/lustre06/project/6064766/datasets/Jisoo/inter/selection_replay';
% mkdir(selection_dir);
% 
% addpath('/home/ecar/Github/Replay_REM')
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
PARAMS.decoding.numshuffles = 1000; % calculate average decoding error with 1000 shuffled data

% set the parameters for replay
PARAMS.replay.step_size=5; %step_size ; overlapped portion. 
PARAMS.replay.windowsize=14; %windowsize=duration of time window.
PARAMS.replay.numshuffles = 1000; 
PARAMS.replay.sampling_threshold=0.7;%how many data points for each time window
PARAMS.replay.jumpiness=50;
PARAMS.replay.min_slope=1;
PARAMS.replay.max_slope=10;


%% Place cell script here
cd(inter_dir);

fnames = dir('*_data.mat');


for iF = 1:length(fnames)
    warning off
    load(fnames(iF).name)
    warning on
    fprintf('Extracting place cells for: %s   %s....', info.subject, info.session)
    tic
    
    % set the parameters for data

    PARAMS.data.ca_time= ms.time/1000;
    PARAMS.data.ca_data=ms.RawTraces ;
    PARAMS.data.behav_time=behav.time/1000;
    PARAMS.data.behav_vec=behav.position(:,1);
    PARAMS.data.num_surrogates=1000;
    
    [PCs_properties] = extract_place_cells(inter_dir, info, PARAMS, ms, behav);
    toc
    
    save([inter_dir filesep info.subject '_' info.session '_PCs.mat'], 'PCs_properties')
    
    fprintf('done\n')
end



%% Decoding
cd(inter_dir);

fnames = dir('*_data.mat');
fnames_PC=dir('*_PCs.mat');


for iF = 1:length(fnames)
    warning off
    load(fnames(iF).name)
    load(fnames_PC(iF).name)
    
    %load place cell output here
    warning on
    fprintf('Generating TCs for: %s   %s....', info.subject, info.session)
    tic
    
    % set the parameters for data
    PARAMS.data.ca_time= ms.time/1000;
    PARAMS.data.ca_data=ms.RawTraces ;
    PARAMS.data.behav_time=behav.time/1000;
    PARAMS.data.behav_vec=behav.position(:,1);
    
    [decoding] = Bayesian_JS2(decoding_dir, info, PARAMS, ms, behav, all_binary_pre_REM, all_binary_post_REM,PCs_properties);
    toc
    
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
   
    [Replay] = Replay_Fidelity_linear_regression(PARAMS,replay_dir,decoding);
    toc
            
 
    save([replay_dir filesep Replay.info.subject '_' Replay.info.session '_Replay.mat'], 'Replay')
  
    clear decoding 
    
    fprintf('done\n')
end



%% Replay selection

cd(replay_dir);
fnames = dir('*_Replay.mat');
for iF = 1:length(fnames)

    warning off
    load(fnames(iF).name)
    warning on
    fprintf('Generating TCs for: %s   %s....', Replay.info.subject, Replay.info.session)
    tic
   
    [Selected_replay]=Replay_selection(PARAMS, [], Replay)
    toc

    
    save([selection_dir filesep Replay.info.subject '_' Replay.info.session '_selected_replay.mat'], 'Selected_replay')
        
     clear Replay
     
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


%% plotting for replay

figure;
scatter(Replay.sampling_per_actual,Replay.Replay_score_actual,'*');
ylim([0 1]);
title('scatter plot of sampling percentage and replay score');
xlabel('Sampling Percentage');
ylabel('Replay score')
xline(0.7,'r','LineWidth',2);

%plotting slope histogram
figure;
histogram(Replay.Replay_slope_actual, 'EdgeAlpha',0);
xlabel('Slope')
ylabel('Frequency')
box off



%% plotting CI95, 99
% Getting confidence Intervals using position shuffle 
x1=shuffle_1;
figure;
histogram(x1);
CIFcn = @(x1,p)prctile(x1,abs([0,100]-(100-p)/2));
CI_95_Shuffle_position = CIFcn(x1,95); 
arrayfun(@(x1)xline(x1,'-m','95 prctile'),CI_95_Shuffle_position(2));
xlabel('Replay score')
ylabel('Frequency')
legend('Position shuffle');
box off

% Getting confidence Intervals using time shuffle
x2=shuffle_2;
figure;
histogram(x2);
CIFcn = @(x2,p)prctile(x2,abs([0,100]-(100-p)/2)); % this doesn't need to be normal disgtribution
CI_95_Shuffle_time = CIFcn(x2,95); 
arrayfun(@(x2)xline(x2,'-m','95 prctile'),CI_95_Shuffle_time(2));
xlabel('Replay score')
ylabel('Frequency')
legend('Time shuffle');
box off

x1=shuffle_p_score;
figure;
histogram(x1);
CIFcn = @(x1,p)prctile(x1,abs([0,100]-(100-p)/2));
CI_99_Shuffle_position = CIFcn(x1,99); 
arrayfun(@(x1)xline(x1,'-m','99 prctile'),CI_99_Shuffle_position(2));
xlabel('Replay score')
ylabel('Frequency')
legend('Position shuffle');

box off
% Getting confidence Intervals using time shuffle
x2=shuffle_t_score;
figure;
histogram(x2);
CIFcn = @(x2,p)prctile(x2,abs([0,100]-(100-p)/2)); % this doesn't need to be normal disgtribution
CI_99_Shuffle_time = CIFcn(x2,99); 
arrayfun(@(x2)xline(x2,'-m','99 prctile'),CI_99_Shuffle_time(2));
xlabel('Replay score')
ylabel('Frequency')
legend('Time shuffle');
box off
% script for Fig 1


% script for fig 2




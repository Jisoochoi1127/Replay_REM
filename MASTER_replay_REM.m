%%% MASTER_REM_replay

%% Jisoo to do list

% for master script
% DONE_ collect data first as ms is needed for params
% DONE_ Set up parameters -- May 11th 
% need to add line(line 45) for passing some sessions doens't have sleep data
%(ex)pv1254-Dark transparet

% for Bayesian_JS2 script
% DONE_ change Bayesian JS2 variable using params -- May 11th 
% Separate plotting section from Bayesian JS2
% Remove all_binary_pre_SW,all_binary_post_SW in the analysis

% for Replay script
% separate the selection part from replay
% input jumpiness in replay script 



%% set path for data
data_dir = '/Users/jisoo/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/Data';
inter_dir = '/Users/jisoo/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/inter/inter_data';
mkdir(inter_dir);
decoding_dir = '/Users/jisoo/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/inter/decoding';
mkdir(decoding_dir);


%% collect data and generate intermediate files.

cd(data_dir)
sub_list = dir('pv*');

for iSub = length(sub_list):-1:1
    cd([data_dir filesep sub_list(iSub).name]);
    sess_list = dir('*D*');
    
    for iS  = length(sess_list):-1:1
        
        cd([data_dir filesep sub_list(iSub).name filesep sess_list(iS).name]);
        
        % need to add line for pass this if the session doesn't have this
        % structure
        
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
        save([inter_dir filesep info.subject '_' info.session '_data.mat'], 'ms', 'behav', 'all_binary_pre_REM', 'all_binary_pre_SW', 'all_binary_post_REM', 'all_binary_post_SW', 'info')
        warning on
    end
end

%% PARAMS

% set RNG for consistency
rng(1234, 'twister'); 

% set the parameters for data
PARAMS.data.ca_time= ms.time/1000;
PARAMS.data.ca_data=ms.RawTraces ;
PARAMS.data.behav_time=behav.time/1000;
PARAMS.data.behav_vec=behav.position(:,1);

% set the parameters for decoding

PARAMS.decoding.bin_size = 3;
PARAMS.decoding.cell_used = logical(ones(size(ca_data,2),1)); % Use every cell
PARAMS.decoding.sampling_frequency = 30; % This data set has been sampled at 30 images per second
PARAMS.decoding.z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
PARAMS.decoding.min_speed_threshold = 5; % 2 cm.s-1


% set the parameters for replay

PARAMS.replay.step_size=10;
PARAMS.replay.windowsize=29;
PARAMS.replay.numshuffles = 1000; 
PARAMS.replay.sampling_threshold=0.7;%how many data points for each time window
%PARAMS.replay.slope
PARAMS.replay.jumpiness=50;

% style parameters for figures

%% Generate tunning curves
cd(inter_dir);

fnames = dir('*_data.mat');

for iF = 1:length(fnames)
    warning off
    load(fnames(iF).name)
    warning on
    fprintf('Generating TCs for: %s   %s....', info.subject, info.session)
    tic
    [decoding] = Bayesian_JS2(PARAMS, ms, behav, all_binary_pre_REM, all_binary_pre_SW, all_binary_post_SW, all_binary_post_REM);
    toc
    
    % where and how do we save this decoding?
    save([decoding_dir filesep info.subject '_' info.session '_data.mat'], 'decoding')
    
    
    clear decoding ms behav info all_binary_pre_REM all_binary_pre_SW all_binary_post_REM all_binary_post_SW
    
    
    fprintf('done\n')
end


% function for generating TCs goes here. Give it the Inter_dir as an input
% along with the paramters and let it go.


%% Decoding



%% Replay




%% Figures


% script for Fig 1


% script for fig 2




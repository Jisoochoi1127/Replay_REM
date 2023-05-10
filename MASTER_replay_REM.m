%%% MASTER_REM_replay


%% PARAMS

% put all the parameters for the entire project here. Edit as needed for
% different analyses
data_dir = '/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/Data';
inter_dir = '/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/inter/inter_data';
mkdir(inter_dir);
decoding_dir = '/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/inter/decoding';
mkdir(decoding_dir);

% set RNG for consistency
rng(1234, 'twister'); 

% set the parameters 


% decoding paramters
PARAMS.decoding.bin_size = 3;

% style parameters for figures


%% collect data and generate intermediate files.
cd(data_dir)
sub_list = dir('pv*');

for iSub = length(sub_list):-1:1
    cd([data_dir filesep sub_list(iSub).name]);
    sess_list = dir('*D*');
    
    for iS  = length(sess_list):-1:1
        
        cd([data_dir filesep sub_list(iSub).name filesep sess_list(iS).name]);
        
        
        warning off
        load('ms.mat', 'ms')
        load('behav.mat')
        load('all_binary_pre_REM.mat')
        load('all_binary_pre_SW.mat')
        load('all_binary_post_REM.mat')
        load('all_binary_post_SW.mat')
        
        info = [];
        info.subject = sub_list(iSub).name;
        info.session = sess_list(iS).name;
        info.nCells = ms.numNeurons;
        
        
        fprintf('Saving %s   %s ...\n', info.subject, info.session)
        save([inter_dir filesep info.subject '_' info.session '_data.mat'], 'ms', 'behav', 'all_binary_pre_REM', 'all_binary_pre_SW', 'all_binary_post_REM', 'all_binary_post_SW', 'info')
        warning on
    end
end



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




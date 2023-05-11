%For Bayesian_JS2 script


%path=/Users/jisoo/Williams Lab Dropbox/Williams Lab Team Folder/Jisoo/Manuscript/Data;


% Add id, data and loop.



%loading input
load 'ms'
load 'behav'
load 'all_binary_post_REM'
load 'all_binary_post_SW'
load 'all_binary_pre_REM'
load 'all_binary_pre_SW'

% Parameter
ca_time= ms.time/1000;
ca_data=ms.RawTraces ;
behav_time=behav.time/1000;
behav_vec=behav.position(:,1);


%% decoding _ wake, REM, SWS
bin_size = 3; %bin size for decoding
cell_used = logical(ones(size(ca_data,2),1)); % Use every cell

[decoding]=Bayesian_JS2(ms,behav,all_binary_post_REM,all_binary_post_SW,all_binary_pre_REM,all_binary_pre_SW);

%% Replay Line fitting

% Parameter
step_size=10;
windowsize=29;
numshuffles = 1000; %for line fitting shuffle
sampling_threshold=0.7; %how many data points for each time window

% this part, we need to change
Final_Replay_idx=find(Sig_max_jumpiness_99 <50 & abs(Sig_slope_99)>1 & abs(Sig_slope_99)<10)

[out]=Replay_Fidelity_linear_regression(decoding)

  
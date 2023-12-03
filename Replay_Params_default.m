function PARAMS = Replay_Params_default

   
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
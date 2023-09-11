function [PC_properties]= Bayesian_JS2(inter_dir, info, PARAMS, ms, behav)
rng(PARAMS.rng, 'twister'); 

%% Interpolate behav and calcium
[unique_behavtime, IAbehav, ICbehav]=unique(PARAMS.data.behav_time);
[unique_mstime, IAms, ICms]=unique(PARAMS.data.ca_time);
interp_behav_vec = interp1(IAbehav,PARAMS.data.behav_vec(IAbehav,1),IAms);
interp_behav_vec(end) = interp_behav_vec(end-1);

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[velocity] = extract_velocity(interp_behav_vec, PARAMS.data.ca_time);

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis

running_ts = velocity > PARAMS.decoding.min_speed_threshold;

%% Compute occupancy and joint probabilities
% Make sure that your binning vector includes every data point of
% interp_behav_vec using the min/max function:
bin_vector = 0:PARAMS.decoding.bin_size:100+PARAMS.decoding.bin_size; % start : bin_size : end
bin_centers_vector = bin_vector + PARAMS.decoding.bin_size/2;
bin_centers_vector(end) = [];

% First let's binarize traces from all cells #TODO serialize w/ decoding
binarized_data = zeros(size(PARAMS.data.ca_data));
for cell_i = 1:size(PARAMS.data.ca_data,2)
    binarized_data(:,cell_i) = extract_binary(PARAMS.data.ca_data(:,cell_i), PARAMS.decoding.sampling_frequency, PARAMS.decoding.z_threshold);
end

%% Shuffle and establish place cells
numShuffles=PARAMS.data.num_surrogates
p_value_threshold = 0.05 % Threshold to be considered a place cell. TODO move to params

for cell_i = 1:size(binarized_data,2)
    [PC_properties.MI(cell_i), ~, ~, PC_properties.marginal_likelihood(cell_i), tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, bin_vector, logical(ones(size(binarized_data,1),1)));
    PC_properties.peak_rate(cell_i), PC_properties.peak_loc(cell_i) = max(tuning_curve_data(:,cell_i))

    PC_properties.isPC(cell_i) = 1 % By default, neuron will be a place cell unless does not pass significance test

    % Shuffle
    shuffled_MIs = zeros(numShuffles);
    num_shuffled_above_actual = 0

    for k = 1:numShuffles
        random_ts = ceil(rand*length(PARAMS.data.ca_time));
        shuffled_binarized = zeros(length(binarized_trace),1);

        % Permute the trace
        shuffled_binarized(1:random_ts) = binarized_trace(end-random_ts+1:end);
        shuffled_binarized(random_ts+1:end) = binarized_trace(1:end-random_ts);
        
        % Compute tuning curve
        [shuffled_MI, ~, ~, ~, ~] = extract_1D_information(shuffled_binarized, interp_behav_vec, bin_vector, running_ts);
        if shuffled_MI>MI(cell_i)
            num_shuffled_above_actual = num_shuffled_above_actual+1
        end

        if num_shuffled_above_actual/numShuffles > p_value_threshold
            PC_properties.isPC(cell_i) = 0 % Too many surrogates outperformed actual data. This is not a place cell
            break % Since this is likely not a place cell, no need to shuffle further
        end
    end 
end

end

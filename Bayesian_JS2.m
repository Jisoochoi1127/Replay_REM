%add Params. input 

function [decoding]= Bayesian_JS2(decoding_dir, info,PARAMS, ms, behav, all_binary_pre_REM, all_binary_post_REM)
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
bin_vector = min(interp_behav_vec):PARAMS.decoding.bin_size:max(interp_behav_vec)+PARAMS.decoding.bin_size; % start : bin_size : end
bin_centers_vector = bin_vector + PARAMS.decoding.bin_size/2;
bin_centers_vector(end) = [];

%% Decoding
% First let's binarize traces from all cells
binarized_data = zeros(size(PARAMS.data.ca_data));
for cell_i = 1:size(PARAMS.data.ca_data,2)
    binarized_data(:,cell_i) = extract_binary(PARAMS.data.ca_data(:,cell_i), PARAMS.decoding.sampling_frequency, PARAMS.decoding.z_threshold);
end

%% Select random frames to train decoder 

if contains(PARAMS.decoding.training_set_creation_method,'Block')
    
    [training_ts] =Analysis_Create_training_decoding_alternate(PARAMS.data.ca_time,PARAMS.decoding.len);
    
   training_ts(running_ts == 0) = 0; % Exclude periods of immobility from the traing set

   
else
    training_ts = create_training_set(PARAMS.data.ca_time,'random',0.5);
    training_ts(running_ts == 0) = 0;
    % Exclude periods of immobility from the traing set
end


%% Create tuning curves for every cell
    for cell_i = 1:size(binarized_data,2)
        [MI(cell_i), PDF(:,cell_i), occupancy_vector, prob_being_active(cell_i), tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, bin_vector, training_ts);
    end
    
    Occupancy_bin=occupancy_vector;
    

%% Decode position

% First, let us establish the timestamps used for decoding.
decoding_ts = ~training_ts; % Training timestamps are excluded
decoding_ts(running_ts == 0) = 0; % Periods of immobility are excluded

% Minimal a priori (use to remove experimental a priori)
occupancy_vector = occupancy_vector./occupancy_vector*(1/length(occupancy_vector));


    
%% Extract  decoded position during Track

if PARAMS.decoding.cell_used == 'Whole_cell'; %Whole_cell, Place_cell, Non_Place_cell
    
    PARAMS.decoding.cell_used=logical(ones(size(PARAMS.data.ca_data,2),1)); % Use every cell
    
end
% else if PARAMS.decoding.cell_used = 'Place_cell'
%         PARAMS.decoding.cell_used= 
%         
%     end
% 
% else if PARAMS.decoding.cell_used = 'Non_Place_cell'
%         PARAMS.decoding.cell_used=
%         
%     end
% end
[WAKE_decoded_probabilities] = bayesian_decode1D(binarized_data, occupancy_vector, prob_being_active, tuning_curve_data, PARAMS.decoding.cell_used);

[max_decoded_prob, decoded_bin] = max(WAKE_decoded_probabilities,[],1);
decoded_position = bin_centers_vector(decoded_bin);

% added by jisoo
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



%% Remove training timestamps to assess decoding error rate
decoded_bin(~decoding_ts) = nan;
decoded_position(~decoding_ts) = nan;
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
mean_decoding_error = mean(abs_decoding_error, 'omitnan');

%% To calculate decoding error in Half left/right in HAT
Half_left_actual=find(actual_position<=50);
Half_right_actual=find(actual_position>50);
Half_left_decoded=decoded_position(Half_left_actual);
Half_right_decoded=decoded_position(Half_right_actual);
decoding_error_half_left = actual_position(Half_left_actual) - Half_left_decoded;
mean_decoding_error_half_left = mean(abs(decoding_error_half_left), 'omitnan');

decoding_error_half_right = actual_position(Half_right_actual) - Half_right_decoded;
mean_decoding_error_half_right = mean(abs(decoding_error_half_right), 'omitnan');



%% Compute confusion matrix
confusion_matrix = zeros(length(bin_centers_vector),length(bin_centers_vector));

for actual_i = 1:length(bin_centers_vector)
   for decoded_i = 1:length(bin_centers_vector)
       confusion_matrix(actual_i,decoded_i) = sum(decoded_bin == decoded_i & actual_bin == actual_i)./sum(actual_bin == actual_i);
   end
end


%% Create tuning curves with shuffled binarized
mean_decoding_error_shuffle=[];
for k=PARAMS.decoding.numshuffles:-1:1
    
    Shuffled_binarized_data =zeros(size(binarized_data));
    random_ts=[];
    for ii=size(binarized_data,2):-1:1
     
        random_ts(ii) = randi(size(binarized_data,1));
        Shuffled_binarized_data(:,ii)=circshift(binarized_data(:,ii),random_ts(ii),1);
      
    end
    
    
for cell_i = size(Shuffled_binarized_data,2):-1:1
    [MI_S(cell_i), PDF_S(:,cell_i), occupancy_vector_S, prob_being_active_S(cell_i), tuning_curve_data_shuffled(:,cell_i) ] = extract_1D_information(Shuffled_binarized_data(:,cell_i), interp_behav_vec, bin_vector, training_ts);
end


[Shuffled_decoded_probabilities] = bayesian_decode1D(Shuffled_binarized_data, occupancy_vector_S, prob_being_active_S, tuning_curve_data_shuffled, PARAMS.decoding.cell_used);

[max_shuffled_decoded_prob, shuffled_decoded_bin] = max(Shuffled_decoded_probabilities,[],1);
decoded_position_shuffle = bin_centers_vector(shuffled_decoded_bin);

%% Remove training timestamps to assess decoding error rate
shuffled_decoded_bin(~decoding_ts) = nan;
decoded_position_shuffle(~decoding_ts) = nan;
Shuffled_decoded_probabilities(:,~decoding_ts) = nan;


%% Compute decoding error
decoding_error_shuffle = abs(actual_position - decoded_position_shuffle);
mean_decoding_error_shuffle(k) = mean(decoding_error_shuffle, 'omitnan');

end
Avg_total_shuffle_error=mean(mean_decoding_error_shuffle);

%% Extract decoded position during post REM sleep

    [REM_decoded_probabilities] = bayesian_decode1D(all_binary_post_REM, occupancy_vector, prob_being_active, tuning_curve_data, PARAMS.decoding.cell_used);
    
    % Maximum a posteriori
    [REM_max_prob, REM_decoded_bin] = max(REM_decoded_probabilities,[],1);
    REM_decoded_position = bin_centers_vector(REM_decoded_bin);
    
    REM_decoded_bin(isnan(REM_max_prob)) = nan;

    REM_decoded_position(isnan(REM_decoded_bin)) = nan;



  
%% Extract decoded position during pre REM sleep

[pre_REM_decoded_probabilities] = bayesian_decode1D(all_binary_pre_REM, occupancy_vector, prob_being_active, tuning_curve_data, PARAMS.decoding.cell_used);

% Maximum a posteriori
[pre_REM_max_prob, pre_REM_decoded_bin] = max(pre_REM_decoded_probabilities,[],1);
pre_REM_decoded_position = bin_centers_vector(pre_REM_decoded_bin);

pre_REM_decoded_bin(isnan(pre_REM_max_prob)) = nan;

pre_REM_decoded_position(isnan(pre_REM_decoded_bin)) = nan;

 

%% Save the output of decoded results

decoding.tuning_curve_data=tuning_curve_data;
decoding.Velocity_bin=Velocity_bin;
decoding.confusion_matrix=confusion_matrix;

decoding.bin_centers_vector=bin_centers_vector;
decoding.WAKE_decoded_probabilities =whole_Wake_decoded_probabilities;
decoding.WAKE_decoded_position=whole_Wake_decoded_position;
decoding.MI=MI;
decoding.PDF=PDF;

decoding.Testing_wake_decoded_probabilities=WAKE_decoded_probabilities;
decoding.Testing_wake_decoded_position=decoded_position;


decoding.REM_decoded_probabilities=REM_decoded_probabilities;
decoding.REM_max_prob=REM_max_prob;
decoding.REM_decoded_position=REM_decoded_position;

decoding.pre_REM_decoded_probabilities=pre_REM_decoded_probabilities;
decoding.pre_REM_max_prob=pre_REM_max_prob;
decoding.pre_REM_decoded_position=pre_REM_decoded_position;


decoding.decoding_error=decoding_error;
decoding.abs_decoding_error=abs_decoding_error;
decoding.mean_decoding_error=mean_decoding_error;
decoding.mean_decoding_error_shuffle=mean_decoding_error_shuffle;
decoding.Avg_total_shuffle_error=Avg_total_shuffle_error;


 end



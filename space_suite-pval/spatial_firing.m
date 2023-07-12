function [ place_field_data ] = spatial_firing(trace, behav, cell_ID, varargin)
%spatial_firing Returns structure with spatial information analysis for
%calcium.
%INPUTS
%   trace    vector with behavioural calcium data
%   behav   struct with track info
%   cell_ID    cell#
% <OPTIONALS> 'params' structure 
% dVector logical directionality vector for laps to be considered
% analysis_method   1 binarized trace, 2 raw trace. default 1   
% min_speed   minimum running speed, default 2cm/s.
% bin_size    size of bin in cm to bin the area. default 4
% smoothing   logical. If we want to smooth the place field or nor. Default True    
% place_field_threshold Threshold in SD over suffled distribution from which place fields are selected. Default 2  
% min_place_field_stability Stability between the split trace from which place fields are selected [0 to 1] Default 0.5
%'numShuffles' % of shuffles to create the random distribution default: 200.  
% OUTPUTS
% place_field_data, binarized or raw for a cell
% FIELDS:
%     IsPlaceCell
%     numPlaceFields
%     Params
%     Cell_ID
%     Z_map
%     PlaceField
%     PlaceFieldPeakJointProbability
%     PlaceFieldCentroid
%     PlaceFieldArea in cm2
%     PlaceFieldStability value ranges 0 to 1
%     InFieldActivityRatio
%     PlaceFieldMask
%     SI spatian information bits per spike
%     sparsity
%     background(if the behavioural track)
%     CellFiringProb  total firing probability of the recording
%     spatial_matrix calcium corrected by time spent in bin
%

%   Based on Guillaume Etter's script
%Author: E Vico

%% Parameter parsing
    p = inputParser;
    default_dVector.diVector= 1+zeros(1,length(trace));
    default_dVector.laps= 0;
    addParameter(p,'dVector',default_dVector)
    addParameter(p,'analysis_method',1,@isnumeric)
    addParameter(p,'min_speed',5,@isnumeric)
    addParameter(p,'bin_size',4,@isnumeric)
    addParameter(p,'smoothing',true,@islogical) 
    addParameter(p,'place_field_threshold',.05,@isnumeric)
    addParameter(p,'min_place_field_stability',0.5,@isnumeric)
    addParameter(p,'numShuffles',1000,@isnumeric)
    
    parse(p,varargin{:})
    dVector = p.Results.dVector;    
    analysis_method = p.Results.analysis_method;
    min_speed = p.Results.min_speed;
    bin_size = p.Results.bin_size;
    smoothing = p.Results.smoothing ;
    place_field_threshold = p.Results.place_field_threshold;
    min_place_field_stability = p.Results.min_place_field_stability;
    numShuffles = p.Results.numShuffles;
    
    X_bin_vector=behav.X_bin_vector;    
    interpolated_X=behav.interpolated_X;
    interpolated_speed=behav.interpolated_speed;
%%
%start analysis      
    spatial_matrix_early = zeros(1,length(X_bin_vector));
    spatial_matrix_late = zeros(1,length(X_bin_vector));
    spatial_matrix_full = zeros(1,length(X_bin_vector));
    S_map = zeros(1,length(X_bin_vector));
    time_in_bin = zeros(1,length(X_bin_vector));
    visited_bins=0;
    
    
    for i = 1:length(X_bin_vector)        
        position_idx = find(interpolated_X>X_bin_vector(i)-bin_size & interpolated_X < X_bin_vector(i));
        position_idx=position_idx(dVector.diVector(position_idx)==1);
        early_position_idx = [];
        late_position_idx = [];
        if ~isempty(position_idx)
            position_idx(interpolated_speed(position_idx)<min_speed)=[]; % Remove indices for particular speeds
            visited_bins=visited_bins+1;

            if length(position_idx) > 2
                early_position_idx = position_idx(1:ceil(length(position_idx)./2));
                late_position_idx = position_idx(ceil(length(position_idx)./2):end);
            elseif length(position_idx) == 2
                early_position_idx = position_idx(1);
                late_position_idx = position_idx(2);
            elseif length(position_idx) == 1
                early_position_idx = [];
                late_position_idx = [];
            end            

            time_in_bin(1,i)= length(position_idx); 
            if analysis_method == 1 % Using binary information                  
                firing_in_bin_idx = find(trace(position_idx) == 1);

                early_firing_in_bin_idx = find(trace(early_position_idx) == 1);
                late_firing_in_bin_idx = find(trace(late_position_idx) == 1);             

                S_map(1,i)= length(firing_in_bin_idx);               
                spatial_matrix_full(1,i) = length(firing_in_bin_idx)/length(position_idx);               
                spatial_matrix_early(1,i) = length(early_firing_in_bin_idx)/length(early_position_idx);
                spatial_matrix_late(1,i) = length(late_firing_in_bin_idx)/length(late_position_idx);

            elseif analysis_method == 2 % Using raw calcium traces

                S_map(1,i)= sum(trace(position_idx));        
                spatial_matrix_full(1,i) = sum(trace(position_idx))./length(position_idx);
                spatial_matrix_early(1,i) = sum(trace(early_position_idx))./length(early_position_idx);
                spatial_matrix_late(1,i) = sum(trace(late_position_idx))./length(late_position_idx); 
            end
        else 
            % if the behavior doesn't reach to the first bin.% added by
            % jisoo for pv1254_HATD1
            early_position_idx = [];
            late_position_idx = [];
            time_in_bin(1,i)= [];
            S_map(1,i)=[];
            spatial_matrix_full(1,i)=[];
            spatial_matrix_early(1,i) =[];
            spatial_matrix_late(1,i)=[];
            
            
        end
        
    end

%% Smoothing of place fields

    if smoothing
        spatial_matrix_full=smooth_mat(spatial_matrix_full, bin_size,'gauss', 1); 
        spatial_matrix_early=smooth_mat(spatial_matrix_early,bin_size,'gauss'); 
        spatial_matrix_late=smooth_mat(spatial_matrix_late, bin_size,'gauss'); 
        S_map=smooth_mat(S_map, bin_size,'gauss'); 
    else
        spatial_matrix_early(isnan(spatial_matrix_early))=0;
        spatial_matrix_late(isnan(spatial_matrix_late))=0;
    end

    %% Compute stability
    if sum(abs(spatial_matrix_early))  ==0 || sum(abs(spatial_matrix_late))  ==0
        place_field_stability=0;
    else
        place_field_stability = abs(corr2(spatial_matrix_early,spatial_matrix_late));
    end

    %% Perform shuffling
    shuffled_matrix = zeros(numShuffles, length(X_bin_vector));

    if numShuffles > 0
        for shuffle_i = 1:numShuffles
            random_frame = randi([300 length(trace)-300]);    %random frame value between 1 and input
            for i = 1:length(X_bin_vector)                
                position_idx = find(interpolated_X>X_bin_vector(i)-bin_size & interpolated_X < X_bin_vector(i));
                position_idx=position_idx(dVector.diVector(position_idx)==1);

                if ~isempty(position_idx)
                    position_idx(interpolated_speed(position_idx)<min_speed)=[]; % Remove indices for particular speeds

                    shuffled_trace = zeros(1,length(trace));
                    shuffled_trace = circshift(trace,random_frame);   
                    if analysis_method == 1 % Using binary information                        
                        firing_in_bin_idx = find(shuffled_trace(position_idx) == 1);
                        shuffled_matrix(shuffle_i,i) = length(firing_in_bin_idx)/length(position_idx);

                    elseif analysis_method == 2 % Using raw calcium traces
                        shuffled_matrix(shuffle_i,i) = sum(shuffled_trace(position_idx))./length(position_idx);
                    end
                
                else
                    % if the behavior doesn't reach to the first bin.% added by
                    % jisoo for pv1254_HATD1
                   shuffled_matrix(shuffle_i,i)=[];
                end
            end
        

    %% Smoothing of place fields
            if smoothing
                shuffled_matrix(shuffle_i,:)=smooth_mat(shuffled_matrix(shuffle_i,:), bin_size,'gauss',1); 
            end

        end

    %% Find the p values
        p_map=1 - sum(spatial_matrix_full > shuffled_matrix)/numShuffles;
%p<.05 95 percentile
    else
        p_map =1- (spatial_matrix_full- nanmean(spatial_matrix_full))/nanstd(spatial_matrix_full,0);
    end

    %% Place field analysis

    place_field_mask = p_map <= place_field_threshold & place_field_stability > min_place_field_stability;
    
    place_field = spatial_matrix_full;
    place_field(~place_field_mask) = 0;

    avg_activity_in_place_field = mean(spatial_matrix_full(place_field_mask));
    peak_joint_probability = max(place_field(:));
    avg_activity_out_place_field = nanmean(spatial_matrix_full(~place_field_mask));
    
    if analysis_method==1
        firing_prob = sum(trace)./(size(trace,1)/30);
        if isnan(avg_activity_out_place_field) || avg_activity_out_place_field ==0
            infield_activity_ratio = avg_activity_in_place_field; 
        else
            infield_activity_ratio = avg_activity_in_place_field./avg_activity_out_place_field; %doesnt handle method 2
        end
    else
        firing_prob = NaN;
        if isnan(avg_activity_out_place_field) || avg_activity_out_place_field ==0
            infield_activity_ratio = avg_activity_in_place_field; 
        else
            infield_activity_ratio=abs(avg_activity_in_place_field-avg_activity_out_place_field/avg_activity_out_place_field);
        end%increase calcium from outoffield
    end
    if sum(place_field_mask(:)) ~= 0
        is_place_cell = true;

        reg_props = regionprops(place_field_mask,'Centroid','Area','Eccentricity');
        numPlaceFields = length(reg_props);

        for field_i = 1:numPlaceFields % for each field
            place_field_centroid{field_i} = (reg_props(field_i).Centroid);
            place_field_area{field_i} = (reg_props(field_i).Area)*bin_size;
        end

    else
        is_place_cell = false;
        numPlaceFields = 0;
        place_field_centroid = NaN;
        place_field_area = NaN;
    end

    %% SI & sparsity calculation
    %occupancy= (time_in_bin/time_lap/visited_bins);
    occupancy= time_in_bin/nansum(time_in_bin);
    if analysis_method ==2
        S_map= S_map -min(S_map)+.001;
    end
    S_map=S_map/nanmean(S_map);
    SI = nansum((occupancy.*S_map).*log2(S_map));
    sparsity=(nansum(((occupancy.*S_map)))^2)./nansum((occupancy.*(S_map.^2)));
    
    %divide laps for plotting.
    lap_mat=zeros(dVector.laps,length(X_bin_vector));
    if dVector.laps ~=0
        for d=1:dVector.laps
            ix1=dVector.times(d,1);
            ix2=dVector.times(d,2);
            lap_v=trace(ix1:ix2);
            position_idx=dsearchn(X_bin_vector',(interpolated_X(ix1:ix2)));
            dd=position_idx(logical([1;diff(position_idx)~=0]));            
            for d2=1:length(dd)
                lap_mat(d,d2)=mean(lap_v(position_idx==dd(d2)));
            end
            
        end
        lap_mat(d,:)=smooth_mat(lap_mat(d,:), bin_size,'gauss',1); 
        if dd(1)> dd(5)
            lap_mat=fliplr(lap_mat);
        end            
    end

    %% Output the data
    place_field_data.IsPlaceCell = is_place_cell;
    place_field_data.numPlaceFields = numPlaceFields;
    place_field_data.Params = p.Results;
    place_field_data.Cell_ID = cell_ID;
    place_field_data.p_map = p_map;
    place_field_data.spatial_matrix=spatial_matrix_full;
    place_field_data.PlaceField = place_field;
    place_field_data.PlaceFieldPeakJointProbability = peak_joint_probability;
    place_field_data.PlaceFieldCentroid = place_field_centroid;
    place_field_data.PlaceFieldArea = place_field_area;
    place_field_data.PlaceFieldStability = place_field_stability;
    place_field_data.InFieldActivityRatio = infield_activity_ratio;
    place_field_data.PlaceFieldMask = place_field_mask;
    place_field_data.SI=SI;
    place_field_data.sparsity=sparsity;
    place_field_data.background=behav.background;
    place_field_data.CellFiringProb=firing_prob;
    place_field_data.occupancy=occupancy;
    place_field_data.unc_spatial_matrix=S_map;
    place_field_data.lap_mat=lap_mat;

end




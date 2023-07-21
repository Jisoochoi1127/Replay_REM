

% check Replay.PARAMS.replay.sampling_threshold;
%DONE change variables matched with linear regression script
% Do we need Final_Replay_actual ? or just actual? - maybe remove actual
% and just use final_replay_actual
%maybe remove 95% because we only use 99%

%% Selecting replay with parameters

function [Selected_replay]=Replay_selection(PARAMS,selection_dir, Replay)

%Selection of replay events with sampling threshold, jumpiness, and slope.
%% Selection of replay using shuffled data
% Remove nan, inf, values less than sampling threshold.
Selected_position=Replay.sampling_percentage_position>Replay.PARAMS.replay.sampling_threshold;
Selected_position_shuffle=Replay.score_shuffle_position(Selected_position);
Selected_position_jumpiness=Replay.Shuffle_position_jumpiness(Selected_position);
Selected_position_slope=Replay.slope_position(Selected_position);

nan_inf_position=~isfinite(Selected_position_shuffle);
Selected_position_shuffle(nan_inf_position)=[];
Selected_position_jumpiness(nan_inf_position)=[];
Selected_position_slope(nan_inf_position)=[];

Selected_time=Replay.sampling_percentage_time>Replay.PARAMS.replay.sampling_threshold;
Selected_time_shuffle=Replay.score_shuffle_time(Selected_time);
Selected_time_jumpiness=Replay.Shuffle_time_jumpiness(Selected_time);
Selected_time_slope=Replay.slope_time(Selected_time);

nan_inf_position=~isfinite(Selected_time_shuffle);
Selected_time_shuffle(nan_inf_position)=[];
Selected_time_jumpiness(nan_inf_position)=[];
Selected_time_slope(nan_inf_position)=[];

%score that exceed sampling rate from position and time shuffles

shuffle_p_score=Selected_position_shuffle; 
shuffle_p_jumpiness=Selected_position_jumpiness;
shuffle_p_slope=Selected_position_slope;

shuffle_t_score=Selected_time_shuffle;
shuffle_t_jumpiness=Selected_time_jumpiness;
shuffle_t_slope=Selected_time_slope;


%% Selection of replay using actual data

Selected_actual=Replay.sampling_per_actual>Replay.PARAMS.replay.sampling_threshold;
actual=Replay.Replay_score_actual(Selected_actual);
nan_inf_actual=~isfinite(actual);
actual(nan_inf_actual)=[];

% don't know why i added this but actual = Final_Replay_actual same.
% for actual data
Replay_actual=isfinite(Replay.Replay_score_actual); % excluding nan, infinite values
Selected_Replay_actual=find(Replay_actual==1 & Replay.sampling_per_actual>Replay.PARAMS.replay.sampling_threshold );
Final_Replay_actual=Replay.Replay_score_actual(Selected_Replay_actual);



%% calculate pvalue 
x1=shuffle_p_score;
CIFcn = @(x1,p)prctile(x1,abs([0,100]-(100-p)/2));
CI_95_Shuffle_position = CIFcn(x1,95); 
CI_99_Shuffle_position = CIFcn(x1,99); 

x2=shuffle_t_score;
CIFcn = @(x2,p)prctile(x2,abs([0,100]-(100-p)/2)); % this doesn't need to be normal disgtribution
CI_95_Shuffle_time = CIFcn(x2,95); 
CI_99_Shuffle_time = CIFcn(x2,99); 


%% Find significant replay using 95 percentile 
%significant_idx_95; which window has significant replay events 
Significant_95=find(Final_Replay_actual>CI_95_Shuffle_position(2) &Final_Replay_actual>CI_95_Shuffle_time(2)); 
Significant_idx_95=Selected_Replay_actual(Significant_95);% now we can find which window has significant replay

Significant_99=find(Final_Replay_actual>CI_99_Shuffle_position(2) &Final_Replay_actual>CI_99_Shuffle_time(2)); 
Significant_idx_99=Selected_Replay_actual(Significant_99);% now we can find which window has significant replay


%% for shuffle

Final_shuffle_position_idx=find(shuffle_p_jumpiness<PARAMS.replay.jumpiness & abs(shuffle_p_slope)>PARAMS.replay.min_slope & abs(shuffle_p_slope)<PARAMS.replay.max_slope );
Final_shuffle_time_idx=find(shuffle_t_jumpiness<PARAMS.replay.jumpiness & abs(shuffle_t_slope)>PARAMS.replay.min_slope  & abs(shuffle_t_slope)<PARAMS.replay.max_slope );

% Score that satisfy jumpiness<50, 1<slope<10 & exceed 99% of shuffle
Final_shuffle_p=shuffle_p_score(Final_shuffle_position_idx);
Numb_shuffle_p_sig=find(Final_shuffle_p>CI_99_Shuffle_position(2) &Final_shuffle_p>CI_99_Shuffle_time(2));

Final_shuffle_p=shuffle_t_score(Final_shuffle_time_idx);
Numb_shuffle_t_sig=find(Final_shuffle_p>CI_99_Shuffle_position(2) &Final_shuffle_p>CI_99_Shuffle_time(2));


    
%% Profile of significant replay 

%for 95% CI
Sig_num_replay_95=length(Significant_idx_95);
Sig_replay_score_95=Replay_score_actual(Significant_idx_95);
Sig_avg_replay_score_95=mean(Sig_replay_score_95)
Sig_slope_95=slope(Significant_idx_95);
Sig_avg_slope_95=mean(Sig_slope_95);

%for 99% CI
Sig_num_replay_99=length(Significant_idx_99);
Sig_replay_score_99=Replay_score_actual(Significant_idx_99);
Sig_avg_replay_score_99=mean(Sig_replay_score_99)
Sig_slope_99=slope(Significant_idx_99);
Sig_avg_slope_99=mean(Sig_slope_99);


%% getting decoded occupancy, extend, jumpiness for each 95%, 99%CI


Sig_location_95=[];

for i=1:length(Significant_idx_95);
   
    time= windowStartPoints(Significant_idx_95(i)):windowStartPoints(Significant_idx_95(i))+PARAMS.replay.windowsize;
   
    Sig_start_frame_95(i)=windowStartPoints(Significant_idx_95(i));
    y=decoding.REM_decoded_position(time);
    y(find(isnan(y)))=[];
    
    Sig_extent_95(i)= abs(y(end)-y(1)); 
    Sig_avg_jumpiness_95(i)=mean(abs(diff(y)));
    Sig_max_jumpiness_95(i)=max(abs(diff(y)));; 
    location=y;
    Sig_location_95=[Sig_location_95,location];
    
   
end

Sig_location_99=[];

for i=1:length(Significant_idx_99);
   
    time= windowStartPoints(Significant_idx_99(i)):windowStartPoints(Significant_idx_99(i))+PARAMS.replay.windowsize;
   
    Sig_start_frame_99(i)=windowStartPoints(Significant_idx_99(i));
    y=decoding.REM_decoded_position(time);
    y(find(isnan(y)))=[];
    
    Sig_extent_99(i)= abs(y(end)-y(1)); % can't work if nan is at the end.
    Sig_avg_jumpiness_99(i)=mean(abs(diff(y)));
    Sig_max_jumpiness_99(i)=max(abs(diff(y))); % average or max?
    location=y;
    Sig_location_99=[Sig_location_99,location]; 
    % To check the extent and location of the replay 
    Sig_location_replay{i}=location;
        
   
end

%% Selection - Threshold for jumpiness, slope for final replay events.
Final_Replay_idx=find(Sig_max_jumpiness_99 <PARAMS.replay.jumpiness & abs(Sig_slope_99)>PARAMS.replay.min_slope & abs(Sig_slope_99)<PARAMS.replay.max_slope)
Final_start_frame=Sig_start_frame_99(Final_Replay_idx);
Final_Replay_score=Sig_replay_score_99(Final_Replay_idx);
Final_Replay_slope=Sig_slope_99(Final_Replay_idx);
Final_Replay_extent=Sig_extent_99(Final_Replay_idx);
Final_Replay_avg_jumpiness=Sig_avg_jumpiness_99(Final_Replay_idx);
Final_Replay_max_jumpiness=Sig_max_jumpiness_99(Final_Replay_idx);
%Final_Replay_location=Sig_location_99(Final_Replay_idx);
for i=1:length(Final_Replay_idx);
Final_Sig_location_replay(i)=Sig_location_replay(Final_Replay_idx(i));
end


%% Calculate ratio
sig_shuff_p=length(Numb_shuffle_p_sig)/length(shuffle_p_score)*100;
sig_shuff_t=length(Numb_shuffle_t_sig)/length(shuffle_t_score)*100;
sig_actual=length(Final_Replay_score)/length(Final_Replay_actual)*100;

%% save all the variables

%check which variable i need for plotting and further analysis.
Selected_replay.Final_Replay_actual=Final_Replay_actual;

Selected_replay.Final_start_frame=Final_start_frame;
Selected_replay.Final_Replay_score=Final_Replay_score;
Selected_replay.Final_Replay_slope=Final_Replay_slope;
Selected_replay.Final_Replay_extent=Final_Replay_extent;
Selected_replay.Final_Replay_avg_jumpiness=Final_Replay_avg_jumpiness;
Selected_replay.Final_Replay_max_jumpiness=Final_Replay_max_jumpiness;
Selected_replay.Final_Sig_location_replay=Final_Sig_location_replay;

Selected_replay.sig_shuff_p=sig_shuff_p;
Selected_replay.sig_shuff_t=sig_shuff_t;
Selected_replay.sig_actual=sig_actual;



end
function [Selected_replay]=Replay_selection(PARAMS,out_replay)

tic
%Selection of replay events with sampling threshold, jumpiness, and slope.
%% Selection of replay using shuffled data
% Remove nan, inf, values less than sampling threshold.
Selected_position=sampling_percentage_position>PARAMS.replay.sampling_threshold;
Selected_position_shuffle=Replay_score_shuffle_position(Selected_position);
Selected_position_jumpiness=Shuffle_position_jumpiness(Selected_position);
Selected_position_slope=slope_position(Selected_position);

nan_inf_position=~isfinite(Selected_position_shuffle);
Selected_position_shuffle(nan_inf_position)=[];
Selected_position_jumpiness(nan_inf_position)=[];
Selected_position_slope(nan_inf_position)=[];

Selected_time=sampling_percentage_time>PARAMS.replay.sampling_threshold;
Selected_time_shuffle=Replay_score_shuffle_time(Selected_time);
Selected_time_jumpiness=Shuffle_time_jumpiness(Selected_time);
Selected_time_slope=slope_time(Selected_time);

nan_inf_position=~isfinite(Selected_time_shuffle);
Selected_time_shuffle(nan_inf_position)=[];
Selected_time_jumpiness(nan_inf_position)=[];
Selected_time_slope(nan_inf_position)=[];
%shuffle_1; score from position shuffled after selection
%shuffle_2; score from time shuffled after selection

shuffle_1=Selected_position_shuffle; %score that exceed sampling rate
shuffle_1_jumpiness=Selected_position_jumpiness;
shuffle_1_slope=Selected_position_slope;

shuffle_2=Selected_time_shuffle;
shuffle_2_jumpiness=Selected_time_jumpiness;
shuffle_2_slope=Selected_time_slope;




%% Selection of replay using actual data

Selected_actual=sampling_percentage>PARAMS.replay.sampling_threshold;
actual=Replay_score_actual(Selected_actual);
nan_inf_actual=~isfinite(actual);
actual(nan_inf_actual)=[];

% don't know why i added this but actual = Final_Replay_actual same.
% for actual data
Replay_actual=isfinite(Replay_score_actual); %1 - excluding nan, infinite values
Selected_Replay_actual=find(Replay_actual==1 & sampling_percentage>PARAMS.replay.sampling_threshold );
% ( if slope has also threshold) abs(slope) >min_slope & abs(slope) < max_slope
Final_Replay_actual=Replay_score_actual(Selected_Replay_actual);



%% calculate pvalue 

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

x1=shuffle_1;
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
x2=shuffle_2;
figure;
histogram(x2);
CIFcn = @(x2,p)prctile(x2,abs([0,100]-(100-p)/2)); % this doesn't need to be normal disgtribution
CI_99_Shuffle_time = CIFcn(x2,99); 
arrayfun(@(x2)xline(x2,'-m','99 prctile'),CI_99_Shuffle_time(2));
xlabel('Replay score')
ylabel('Frequency')
legend('Time shuffle');
box off

%% Find significant replay using 95 percentile 
%significant_idx_95; which window has significant replay events 
Significant_95=find(Final_Replay_actual>CI_95_Shuffle_position(2) &Final_Replay_actual>CI_95_Shuffle_time(2)); 
Significant_idx_95=Selected_Replay_actual(Significant_95);% now we can find which window has significant replay

Significant_99=find(Final_Replay_actual>CI_99_Shuffle_position(2) &Final_Replay_actual>CI_99_Shuffle_time(2)); 
Significant_idx_99=Selected_Replay_actual(Significant_99);% now we can find which window has significant replay


%% for shuffle
Final_shuffle_position_idx=find(shuffle_1_jumpiness<50 & abs(shuffle_1_slope)>1 & abs(shuffle_1_slope)<10 );
Final_shuffle_time_idx=find(shuffle_2_jumpiness<50 & abs(shuffle_2_slope)>1 & abs(shuffle_2_slope)<10 );

% Score that satisfy jumpiness<50, 1<slope<10 & exceed 99% of shuffle
Final_shuffle_1=shuffle_1(Final_shuffle_position_idx);
Numb_shuffle_1_sig=find(Final_shuffle_1>CI_99_Shuffle_position(2) &Final_shuffle_1>CI_99_Shuffle_time(2));

Final_shuffle_2=shuffle_2(Final_shuffle_time_idx);
Numb_shuffle_2_sig=find(Final_shuffle_2>CI_99_Shuffle_position(2) &Final_shuffle_2>CI_99_Shuffle_time(2));


    
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
Final_Replay_idx=find(Sig_max_jumpiness_99 <50 & abs(Sig_slope_99)>1 & abs(Sig_slope_99)<10)
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
out.sig_shuff1=length(Numb_shuffle_1_sig)/length(shuffle_1)*100;
out.sig_shuff2=length(Numb_shuffle_2_sig)/length(shuffle_2)*100;
out.sig_actual=length(Final_Replay_score)/length(Final_Replay_actual)*100;

toc
end
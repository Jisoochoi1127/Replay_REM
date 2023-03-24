%% 1. Binning the time
tic
load decoding

%step_size ; overlapped portion. 
%windowsize=duration of time window.

step_size=5;  %overlap 
 %step_size=2
 %step_size=10;
windowsize=14; % duration for calculating R2
%windowsize=7
% windowsize=29;
%step_size=2;  %overlap 
%windowsize=9; % duration for calculating R2

% To calulcate fidelity with wake

Input=decoding.pre_REM_decoded_position    ;
Input_Prob=decoding.pre_REM_decoded_probabilities  ;

%Input=decoding.WAKE_decoded_position;
%Input_Prob=decoding.WAKE_decoded_probabilities;

% To calulcate fidelity with preREM
%Input=decoding.pre_REM_decoded_position  ;
%Input_Prob=decoding.pre_REM_decoded_probabilities  ;

inputDataTimeVec=1:length(Input);
refTimeVec = 1:step_size:length(inputDataTimeVec);
windowStartPoints = dsearchn(inputDataTimeVec',refTimeVec');
numWindows=length(windowStartPoints);


numshuffles = 1000;
sampling_threshold=0.7; %how many data points for each time window


%% 2. Apply linear regression and calculate R square

for window = 1:numWindows
    
    if  windowStartPoints(window)<=(length(inputDataTimeVec)-windowsize);
        
        time= windowStartPoints(window):windowStartPoints(window)+windowsize;
        
        % if window size exceed the length of the inputdata
    else windowStartPoints(window)>(length(inputDataTimeVec)-windowsize);
        %time= windowStartPoints(window):inputDataTimeVec(end);
        time=[]; 
        
    end
    
    x=time;
    y=Input(time);
    %y=decoded_position(time);
    
    %Ignoring nan values
    nanIndex=isnan(y);
    y(nanIndex) = [];
    x = x(~nanIndex);
    
    % check the smapling point
    sampling_percentage(window)=length(y)/length(time);
    
    
    %calculate Rsquare
    p = polyfit(x,y,1);
    
    slope(window)=p(1);
    intercept(window)=p(2);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq(window) = 1 - SSresid/SStotal;
    Replay_score_actual=rsq;
    
    %plotting (Comment off only to see few examples)
%     bincenter=1:size(decoding.REM_decoded_probabilities,1);
%     figure;
%     scatter(time,decoding.REM_decoded_position(time));
%     % imagesc(time,bincenter,decoding.REM_decoded_probabilities(:,time));
%     caxis([0,0.3]);
%     hold on;
%     %plot(x,yfit,'r','LineWidth',2);
%     plot(time,polyval(p,time),'r','LineWidth',2);
%     xlabel('frames');
%     ylabel('decoded position(cm)');
%     title('linear fitting to decoded position during REM');
%     
end

%plotting sampling percentage for each window time
figure;
scatter(sampling_percentage,Replay_score_actual,'*');
ylim([0 1]);
title('scatter plot of sampling percentage and replay score');
xlabel('Sampling Percentage');
ylabel('Replay score')
yline(0.7,'r','LineWidth',2);

%plotting slope histogram
figure;
histogram(slope, 'EdgeAlpha',0);
xlabel('Slope')
ylabel('Frequency')
box off
%% Shuffling

%input=decoding.REM_decoded_probabilities  ;
Replay_score_shuffle_position=[];
Replay_score_shuffle_time=[];


for k=1:numshuffles;
    % Use decoded probabilities - shuffle- get maximum value. 
    % 
    %  shifting each colum randomly to get shuffled decoded position
    tic
    
    Shuffled_decoded_prob =zeros(size(Input_Prob));
    random_ts=[];
    for i=1:size(Input_Prob,2);
        random_ts(i) = ceil(rand*size(Input_Prob,1));
        
        Shuffled_decoded_prob(:,i)=circshift(Input_Prob(:,i),random_ts(i),1);
      
    end
    
    [Input_max_prob, Input_decoded_bin] = max(Shuffled_decoded_prob,[],1);
    Shuffled_decoded_position = decoding.bin_centers_vector(Input_decoded_bin);
    
    Input_decoded_bin(isnan(Input_max_prob)) = nan;
    % added by jisoo _ nan value for position as well
    
    Shuffled_decoded_position(isnan(Input_decoded_bin)) = nan;


   
    % for time shuffling
    %randperm(n) returns a row vector containing a random permutation of the integers from 1 to n without repeating elements.
    time_shuffle=randperm(size(Input,2));
    Shuffled_time=Input(:,time_shuffle);
    
    
 %% calculate replay score, R2 using shuffled data
    
  
for window = 1:numWindows
    
    if  windowStartPoints(window)<=(length(inputDataTimeVec)-windowsize);
        
        time= windowStartPoints(window):windowStartPoints(window)+windowsize;
        
        % if window size exceed the length of the inputdata
    else windowStartPoints(window)>(length(inputDataTimeVec)-windowsize);
       % time= windowStartPoints(window):inputDataTimeVec(end); 
       time=[]; % if it's not emptied, it cause time shuffling problem that cause high score.
        
    end
    
    x_shuffle_time=time;
    y_shuffle_position=Shuffled_decoded_position(time);
    
    

   % Ignoring nan values
    nanIndex=isnan(y_shuffle_position);
    y_shuffle_position(nanIndex) = [];
    x_shuffle_time = x_shuffle_time(~nanIndex);
    
    sampling_percentage_position(k,window)=length(y_shuffle_position)/length(time);

    
    if ~isempty(x_shuffle_time) & ~isempty(y_shuffle_position)
        p = polyfit(x_shuffle_time,y_shuffle_position,1);
        slope_position(k,window)=p(1);
        yfit_shuffle_position = polyval(p,x_shuffle_time);
        yresid_shuffle_position = y_shuffle_position - yfit_shuffle_position;
        SSresid_shuffle_position = sum(yresid_shuffle_position.^2);
        SStotal_shuffle_position = (length(y_shuffle_position)-1) * var(y_shuffle_position);
        rsq_shuffle_position(window) = 1 - SSresid_shuffle_position/SStotal_shuffle_position;
        
          % Recently added by jisoo _2022.03.14
        if ~isempty(abs(diff(y_shuffle_position)))
        Shuffle_position_jumpiness(k,window)=max(abs(diff(y_shuffle_position)));
        end
        
    else
        
        % remove empty window
        p=[];
        yfit_shuffle_position = [];
        yresid_shuffle_position = [];
        SSresid_shuffle_position = [];
        SStotal_shuffle_position = [];
        rsq_shuffle_position(window)=nan;
         slope_position(k,window)=nan;
        
        
    end
    
    % Shuffling for time
    
     x_shuffle_time=time;
    y_shuffle_time= Shuffled_time(time);
    
    %Ignoring nan values
    nanIndex=isnan(y_shuffle_time);
    y_shuffle_time(nanIndex) = [];
    x_shuffle_time = x_shuffle_time(~nanIndex);
    
    sampling_percentage_time(k,window)=length(y_shuffle_time)/length(time);

    
    if ~isempty(x_shuffle_time) & ~isempty(y_shuffle_time)
        
        p = polyfit(x_shuffle_time,y_shuffle_time,1);
        slope_time(k,window)=p(1);
        yfit_shuffle_time = polyval(p,x_shuffle_time);
        yresid_shuffle_time = y_shuffle_time - yfit_shuffle_time;
        SSresid_shuffle_time = sum(yresid_shuffle_time.^2);
        SStotal_shuffle_time = (length(y_shuffle_time)-1) * var(y_shuffle_time);
        rsq_shuffle_time(window) = 1 - SSresid_shuffle_time/SStotal_shuffle_time;
        
         % Recently added by jisoo _2022.03.14
        if ~isempty(abs(diff(y_shuffle_time)))
            
            Shuffle_time_jumpiness(k,window)=max(abs(diff(y_shuffle_time)));
        end
%         if rsq_shuffle_time(window)==1 & sampling_percentage_time(k,window)>0.7
%           figure;scatter(time, Shuffled_time(time));
% 
%         end
        
    else
        
        p=[];
        yfit_shuffle_time = [];
        yresid_shuffle_time = [];
        SSresid_shuffle_time = [];
        SStotal_shuffle_time = [];
        rsq_shuffle_time(window)=nan;
        slope_time(k,window)=nan;
        
    end
    
end

%Caculate replay score (R2) using position shuffle
shuffle_position=rsq_shuffle_position;
Replay_score_shuffle_position=[Replay_score_shuffle_position;shuffle_position]

%Caculate replay score (R2) using time shuffle
shuffle_time=rsq_shuffle_time;
Replay_score_shuffle_time=[Replay_score_shuffle_time;shuffle_time]
fprintf('Shuff #%0.0f  took ', k)
toc
fprintf('\n')

end
% figure;
% scatter(rsq_shuffle_position, sampling_percentage_position,'*');
% xlim([0 1]);
% title('scatter plot of sampling percentage and replay score(shuffled position)');
% xlabel('Replay score');
% ylabel('sampling_percentage')

% figure;
% scatter(rsq_shuffle_time, sampling_percentage_time,'*');
% xlim([0 1]);
% title('scatter plot of sampling percentage and replay score (shuffled time)');
% xlabel('Replay score');
% ylabel('sampling_percentage')

%% Selection of replay
% Remove nan, inf, values less than sampling threshold.
Selected_position=sampling_percentage_position>sampling_threshold;
Selected_position_shuffle=Replay_score_shuffle_position(Selected_position);
Selected_position_jumpiness=Shuffle_position_jumpiness(Selected_position);
Selected_position_slope=slope_position(Selected_position);


nan_inf_position=~isfinite(Selected_position_shuffle);
Selected_position_shuffle(nan_inf_position)=[];
Selected_position_jumpiness(nan_inf_position)=[];
Selected_position_slope(nan_inf_position)=[];


Selected_time=sampling_percentage_time>sampling_threshold;
Selected_time_shuffle=Replay_score_shuffle_time(Selected_time);
Selected_time_jumpiness=Shuffle_time_jumpiness(Selected_time);
Selected_time_slope=slope_time(Selected_time);

nan_inf_position=~isfinite(Selected_time_shuffle);
Selected_time_shuffle(nan_inf_position)=[];
Selected_time_jumpiness(nan_inf_position)=[];
Selected_time_slope(nan_inf_position)=[];

%shuffle_1; score from position shuffled after selection
%shuffle_2; score from time shuffled after selection



shuffle_1=Selected_position_shuffle;
shuffle_1_jumpiness=Selected_position_jumpiness;
shuffle_1_slope=Selected_position_slope;

shuffle_2=Selected_time_shuffle;
shuffle_2_jumpiness=Selected_time_jumpiness;
shuffle_2_slope=Selected_time_slope;


Selected_actual=sampling_percentage>sampling_threshold;
actual=Replay_score_actual(Selected_actual);
nan_inf_actual=~isfinite(actual);
actual(nan_inf_actual)=[];

% don't know why i added this but actual = Final_Replay_actual same.
% for actual data
Replay_actual=isfinite(Replay_score_actual); %1 - excluding nan, infinite values
Selected_Replay_actual=find(Replay_actual==1 & sampling_percentage>sampling_threshold );
% ( if slope has also threshold) abs(slope) >min_slope & abs(slope) < max_slope
Final_Replay_actual=Replay_score_actual(Selected_Replay_actual);


%Plot histogram of shuffled, actual distribution after selection
figure;
histogram(shuffle_1,'BinWidth',0.025,'EdgeAlpha',0);
hold on;
histogram(shuffle_2,'BinWidth',0.025,'EdgeAlpha',0);
hold on;
histogram(actual,'BinWidth',0.025,'EdgeAlpha',0);
legend('Space shuffle','Time shuffle','Actual data');
xlim([0 1]);
xlabel('Replay score');
ylabel('Probability');
title('Replay score excludidng nan and inf values');


%plotting final replay with actual data
figure; 
histogram(Final_Replay_actual,'BinWidth',0.025,'EdgeAlpha',0)
xlabel('Replay score')
ylabel('Frequency')
title('Replay score histogram with selected actual data')

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
   
    time= windowStartPoints(Significant_idx_95(i)):windowStartPoints(Significant_idx_95(i))+windowsize;
   
    Sig_start_frame_95(i)=windowStartPoints(Significant_idx_95(i));
    y=Input(time);
    y(find(isnan(y)))=[];
    
    Sig_extent_95(i)= abs(y(end)-y(1)); 
    Sig_avg_jumpiness_95(i)=mean(abs(diff(y)));
    Sig_max_jumpiness_95(i)=max(abs(diff(y)));; 
    location=y;
    Sig_location_95=[Sig_location_95,location];
    
   
end

Sig_location_99=[];

for i=1:length(Significant_idx_99);
   
    time= windowStartPoints(Significant_idx_99(i)):windowStartPoints(Significant_idx_99(i))+windowsize;
   
    Sig_start_frame_99(i)=windowStartPoints(Significant_idx_99(i));
    y=Input(time);
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

%% Plotting normalized shuffle, actual replay score

histogram(shuffle_1,'Normalization','probability','BinWidth',0.01,'DisplayStyle','stairs','LineWidth',2,'BinLimits',[0,1]);
hold on;
histogram(shuffle_2,'Normalization','probability','BinWidth',0.01,'DisplayStyle','stairs','LineWidth',2,'BinLimits',[0,1]);
hold on;
histogram(Final_Replay_actual,'Normalization','probability','BinWidth',0.01,'DisplayStyle','stairs','LineWidth',2,'BinLimits',[0,1]);
xlabel('Replay score')
ylabel('Normalized count')
legend('position shuffle','time shuffle','actual data', 'Location','northwest')
box off
%% Calculate ratio
out.sig_shuff1=length(Numb_shuffle_1_sig)/length(shuffle_1)*100;
out.sig_shuff2=length(Numb_shuffle_2_sig)/length(shuffle_2)*100;
out.sig_actual=length(Final_Replay_score)/length(Final_Replay_actual)*100;


%%  With different jumpiness threshold
% Final_Replay_idx=find(Sig_max_jumpiness_99 <40 & abs(Sig_slope_99)>1 & abs(Sig_slope_99)<10)
% Final_start_frame=Sig_start_frame_99(Final_Replay_idx);
% Final_Replay_score=Sig_replay_score_99(Final_Replay_idx);
% Final_Replay_slope=Sig_slope_99(Final_Replay_idx);
% Final_Replay_extent=Sig_extent_99(Final_Replay_idx);
% Final_Replay_avg_jumpiness=Sig_avg_jumpiness_99(Final_Replay_idx);
% Final_Replay_max_jumpiness=Sig_max_jumpiness_99(Final_Replay_idx);
% Final_Replay_location=Sig_location_99(Final_Replay_idx);


toc
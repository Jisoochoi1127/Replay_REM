                                                                               %% calculating replay fidelity using linear regression
% made by jisoo
%April,22 2021

% To do list (May 10, 2021)
% Done-     make time window script _Done
% Done-     Understand what p indicates?
% Done-     plotting decoded position with linear plot.
% Done-     check how R2 match with linear fit. ; 
% Summary ; There are cases that nan value more than ~.check the jumpiness
 % should add criteria for replay - consequtive. how nan value? shorten the
 % time winodw?? (For loop to get maximum R2 value)
 
 
% To do list (May 11, 2021)
%Ongoing-    Check the criteria for chossing candidate events
%Done-     check the paper for shuffle(Davidson2009,Drieu 2019
%           nest,Zielinski 2021
%Done-     script for shuffle data (Column shuffle)
%Done-     calculate R2 for shuffled data

% To do list (May 12, 2021)
%Done-      check polyfit - optimizing?(nope, yresid?) 
%Done-      Check the criteria for chossing candidate events
%Done-      fix -shuffle; location of decoded position.     
%Done-      Check the shuffle part - with data. small number of shuffle.         
%Done-      Watch pvalue youtube video
%Done-      Think about how to set criteria for replay selection
%            (use diff for check jumpiness, check % of nan in window)
% summary - issue ; nan value, Inf, 1 value. - narrow timewindow, criteria?

% To do list (May 13, 2021)
%Delayed      check when replay score reach 1.
%Done       - Ask Eric about pvalue
%Done       - plot distribution of Shuffle, actual replay score
%           Calculate Pvalue , 
%          Check significant replay compared to shuffled .
%          Plot the significant replay (95%, 99%...)

% To do list (May 14, 2021)
%Done        check when replay score reach 1,- inf, -0.667..)
%           Calculate Pvalue , 
%          Check significant replay compared to shuffled .
%          Plot the significant replay (95%, 99%...)
%Replay ? 1. Slope, 2.. Extent (distance) 3.#of cells involved 4. place cell/anxiety cell proportion.
% Add function for unit identity shuffle
% input, output
% Replay score with different time window 
%check replay is not detected between different episodes.

%% To do list (May 17, 2021)
% Done         check the criteria of sampling point & histogram
% Done         scatter plot - number of sample, replay score
% Done         Selection - sampling >50% point (nan excluding), slope (min,
% Done         max value), for shuffle as well?
% Done          add more shuffling - time random 
% Done           check time shuffle with data

% To do list (May 18, 2021)
% Done         selection and then pvalue? check order in the paper
% shuffle - criteria ? - slope ,sampling percentage, jumpiness
% plot the final replay after selection


% To do list (May 19, 2021)
% Done        calculate percentile . 99, 95
% Done        plot 99, 95 in distribution
%              shuffle - slope ,sampling percentage, jumpiness
% Done        Check significant replay compared to shuffled using 95 CI.
% Done         plot significant replay(pv1043 LTD1- 12 , 627 best.
% Done          pv1060_HATD1 - best ; 41, 53
% Done        profile number of replay, slope, mean score, distribtuion
%            think the way to plot the singificant replay 



% To do list (May 20, 2021)
%Done       extent, location , jumpiness of replay
%           Check paper-extent, location , jumpiness of replay
%Done       shuffling K part check. 
% sampling percentage - 50% check,



%% To do list (May 24,2021)
%Done       shuffle polyfiterror- when both inputs are empty
%Done       shuffle polyfit =[], when x, y are both empty..
%Done       ask Eric - shuffling - average or use the whole scores
% if i use whole distribution - there is still high score - 

%summary - set threshold for both shuffle, actual for sampling percentage.

%% To do list (May 21, 2021)
%
%         colum - change shuffling - each colum randomly shift.


% Done       Slope - check how it change with LT, HAT
% same procedure for timewindow. 
%plotting significant replay
% Decoding - lower the threshold and then try decoding during running again
%check window overlapped portion - other paper. 

% check replay is not occur between episodes
% replay at the end of the time window - not same length?


%Replay - when it's occur during each episode
% Replay - with theta power

%% rerun decoding -to get bincenter vector.
% why sampling_percentage_position all same? - make sense; because it's
% same.
%Done     Select the shuffle only from (sampling percentage>50)? -position
%Done     Select the shuffle only from (sampling percentage>50)? -time
%Done     Check histogram of shuffle -position,time after selection
%     Actual - have the original form to find significant replay

%% June, 1 2021

%Done     Check the whole script part by part
%Done     Check shuffle is working well (time, position)
%Done     plotting after shuffling to verify shuffle
%Done     Include both CI 95, 99 results 
%Done     Plotting decoded position that has significant replay

%% June 2 , 2021
% Done     plot 99 Percentile replay - need window time and decoded prob.
% Done     Run 1000 shuffles. - took 2hrs....(pv1060_LTD1),
% HATD5,HATDSwitch- switch has error?
% Done  Create plotting script for significant replay events

%% June 3, 2021

%       use different time window (larger 1s- 30frames)with overlapped- 
%  Some session - TIme shuffle has still high score after removing sampling
%  threshold.... - why?
% Check other replay paper.

%% June 7, 2021
%Done   Time shuffle - why it leads high score? for LTD5
%Done   Whenever, it shuffle 10 times, it gives me different 99CI. 
%Done    need more  shuffling- Try pv1060_LTD5 1000shuffle

%Done    check other shuffle methods - other literatures.
%  Time window - why we decide this window

%% June 8, 2021
%Done    stop at the replay score 1 in shuffled time. (it's the end frames)
%Done    check how replay looks like
%Done    solve time shuffling issue by empty the last window.



%% June 9, 2021
%Done   check other session . HATD5, HATDSwitch- pv1060
% - seems like replay decreases with Habituation
%Done   check other mice pv 1069 - time shuffling too low?

%% June 10, 2021
%Done  Check Replay during day 1 versus day 5 pv1043
%Done  Check Replay during day 1 versus day 5 pv1060
%Done  Check Replay during day 1 versus day 5 pv1069
%summary ; replay events day1 >day5, but jumpiness higher?

%% Jun 11, 2021
%Done       Plot replay events at day 5 pv1060
%Done       running 1000 shuffling _pv1060 LTD1
%Done       running 1000 shuffling_pv1043_LTD1
%Done       add colormap(hot) in plotting script 
% Check replay events with Theta amplitude/frequency
% Run shuffle 1000
%before run, maybe check slope again.

% Remove jumpiness>50
% think about the slope again

%% June 18, 2021
% Modified   Sig_avg_jumpiness_99(i)=mean(abs(diff(y)));
    %Sig_max_jumpiness_99(i)=max(abs(diff(y)));
% without abs, it cause error to find max jumpiness in negative value.
%% 1. Binning the time
tic
load decoding
%load decoding_training90.mat


%step_size ; overlapped portion. 
%windowsize=duration of time window.

%step_size=5;  %overlap 
 %windowsize=14; % duration for calculating R2
 step_size=10;
  windowsize=29;
  
%   step_size=20;
%   windowsize=59;
%step_size=2;  %overlap 
%windowsize=9; % duration for calculating R2
inputDataTimeVec=1:length(decoding.REM_decoded_position);
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
    y=decoding.REM_decoded_position(time);
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
xline(0.7,'r','LineWidth',2);

%plotting slope histogram
figure;
histogram(slope, 'EdgeAlpha',0);
xlabel('Slope')
ylabel('Frequency')
box off
%% Shuffling

%input=decoding.REM_decoded_probabilities  ;
input=decoding.REM_decoded_position;
Replay_score_shuffle_position=[];
Replay_score_shuffle_time=[];
Shuffle_position_jumpiness=NaN(numshuffles,length(Replay_score_actual));
Shuffle_time_jumpiness=NaN(numshuffles,length(Replay_score_actual));

warning off
for k=1:numshuffles;
    tic
    % Use decoded probabilities - shuffle- get maximum value. 
    % 
    %  shifting each colum randomly to get shuffled decoded position
    
    Shuffled_decoded_prob =zeros(size(decoding.REM_decoded_probabilities));
    random_ts=[];
    for i=1:size(decoding.REM_decoded_probabilities,2);
        random_ts(i) = ceil(rand*size(decoding.REM_decoded_probabilities,1));
        
        Shuffled_decoded_prob(:,i)=circshift(decoding.REM_decoded_probabilities(:,i),random_ts(i),1);
      
    end
    
    [REM_max_prob, REM_decoded_bin] = max(Shuffled_decoded_prob,[],1);
    Shuffled_decoded_position = decoding.bin_centers_vector(REM_decoded_bin);
    
    REM_decoded_bin(isnan(REM_max_prob)) = nan;
    % added by jisoo _ nan value for position as well
    
    Shuffled_decoded_position(isnan(REM_decoded_bin)) = nan;


   
    % for time shuffling
    %randperm(n) returns a row vector containing a random permutation of the integers from 1 to n without repeating elements.
    time_shuffle=randperm(size(input,2));
    Shuffled_time=input(:,time_shuffle);
    
    
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
%         Shuffle_position_jumpiness(k,window)=nan;
        
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
      %  Shuffle_time_jumpiness(k,window)=nan;
        
    end
    
end

%Caculate replay score (R2) using position shuffle
shuffle_position=rsq_shuffle_position;
Replay_score_shuffle_position=[Replay_score_shuffle_position;shuffle_position];

%Caculate replay score (R2) using time shuffle
shuffle_time=rsq_shuffle_time;
Replay_score_shuffle_time=[Replay_score_shuffle_time;shuffle_time];
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

shuffle_1=Selected_position_shuffle; %score that exceed sampling rate
shuffle_1_jumpiness=Selected_position_jumpiness;
shuffle_1_slope=Selected_position_slope;

shuffle_2=Selected_time_shuffle;
shuffle_2_jumpiness=Selected_time_jumpiness;
shuffle_2_slope=Selected_time_slope;




%%
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
   
    time= windowStartPoints(Significant_idx_99(i)):windowStartPoints(Significant_idx_99(i))+windowsize;
   
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

%% to do statistics for above plot
% a=isnan()
% x=[shuffle_1,Final_Replay_actual]
% p = kruskalwallis(x)



%%  With different jumpiness threshold
% Final_Replay_idx=find(Sig_max_jumpiness_99 <40 & abs(Sig_slope_99)>1 & abs(Sig_slope_99)<10)
% Final_start_frame=Sig_start_frame_99(Final_Replay_idx);
% Final_Replay_score=Sig_replay_score_99(Final_Replay_idx);
% Final_Replay_slope=Sig_slope_99(Final_Replay_idx);
% Final_Replay_extent=Sig_extent_99(Final_Replay_idx);
% Final_Replay_avg_jumpiness=Sig_avg_jumpiness_99(Final_Replay_idx);
% Final_Replay_max_jumpiness=Sig_max_jumpiness_99(Final_Replay_idx);
% Final_Replay_location=Sig_location_99(Final_Replay_idx);

%% Replya_extent=min, max of   
%% plotting significant replay events - % collect p for each window.
        
%      time= windowStartPoints(window):windowStartPoints(window)+windowsize;
%     figure;
%     scatter(time,decoding.REM_decoded_position(time));
%     hold on; 
%     plot(time,polyval(p,time),'r','LineWidth',2);

%% After getting significant value, Selection process

% Continuity of position
% Continuity of decoding 
% d<...? yresid

% unit identity shuffle
toc
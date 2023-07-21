                                                                               %% calculating replay fidelity using linear regression
% check decoding.PARAMS. 
%% 1. Binning the time
function [Replay]=Replay_Fidelity_linear_regression(PARAMS,replay_dir,decoding);

rng(PARAMS.rng, 'twister');

inputDataTimeVec=1:length(decoding.REM_decoded_position);
refTimeVec = 1:PARAMS.replay.step_size:length(inputDataTimeVec);
windowStartPoints = dsearchn(inputDataTimeVec',refTimeVec');
numWindows=length(windowStartPoints);

%% 2. Apply linear regression and calculate R square

for window = 1:numWindows
    
    if  windowStartPoints(window)<=(length(inputDataTimeVec)-PARAMS.replay.windowsize);
        
        time= windowStartPoints(window):windowStartPoints(window)+PARAMS.replay.windowsize;
        
        % if window size exceed the length of the inputdata
    else windowStartPoints(window)>(length(inputDataTimeVec)-PARAMS.replay.windowsize);
        %time= windowStartPoints(window):inputDataTimeVec(end);
        time=[]; 
        
    end
    
    x=time;
    y=decoding.REM_decoded_position(time);
    
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
    
    if ~isempty(abs(diff(y)))
        actual_jumpiness(window)=max(abs(diff(y)));
        end
end


%% Shuffling

input=decoding.REM_decoded_position;
Replay_score_shuffle_position=[];
Replay_score_shuffle_time=[];
Shuffle_position_jumpiness=NaN(PARAMS.replay.numshuffles,length(Replay_score_actual));
Shuffle_time_jumpiness=NaN(PARAMS.replay.numshuffles,length(Replay_score_actual));

warning off
for k=1:PARAMS.replay.numshuffles;
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
    
    %REM max prob nan but it finds value in REM decoded bin...?
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
    
    if  windowStartPoints(window)<=(length(inputDataTimeVec)-PARAMS.replay.windowsize);
        
        time= windowStartPoints(window):windowStartPoints(window)+PARAMS.replay.windowsize;
        
        % if window size exceed the length of the inputdata
    else windowStartPoints(window)>(length(inputDataTimeVec)-PARAMS.replay.windowsize);
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
        rsq_shuffle_position(k,window) = 1 - SSresid_shuffle_position/SStotal_shuffle_position;
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
        rsq_shuffle_position(k,window)=nan;
        slope_position(k,window)=nan;
        Shuffle_position_jumpiness(k,window)=nan;
        
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
        rsq_shuffle_time(k,window) = 1 - SSresid_shuffle_time/SStotal_shuffle_time;
        % Recently added by jisoo _2022.03.14
        if ~isempty(abs(diff(y_shuffle_time)))
            
            Shuffle_time_jumpiness(k,window)=max(abs(diff(y_shuffle_time)));
        end
        
        
    else
        
        p=[];
        yfit_shuffle_time = [];
        yresid_shuffle_time = [];
        SSresid_shuffle_time = [];
        SStotal_shuffle_time = [];
        rsq_shuffle_time(k,window)=nan;
        slope_time(k,window)=nan;
       Shuffle_time_jumpiness(k,window)=nan;
        
    end
    
end

%Caculate replay score(R2), slope, jumpiness, sampling perc using position, time shuffle
Replay_score_shuffle_position=rsq_shuffle_position;
Replay_score_shuffle_time=rsq_shuffle_time;
Replay_slope_shuffle_position=slope_position;
Replay_jumpiness_shuffle_position=Shuffle_position_jumpiness;
Replay_slope_shuffle_time=slope_time;
Replay_jumpiness_shuffle_position=Shuffle_time_jumpiness;
sampling_per_shuffle_positon=sampling_percentage_position;
sampling_per_shuffle_time=sampling_percentage_time;

fprintf('Shuff #%0.0f  took ', k)
toc
fprintf('\n')
end

%Save variables
Replay.PARAMS=PARAMS;
Replay.info=decoding.info;
Replay.windowStartPoints=windowStartPoints;
Replay.numWindows=numWindows;

Replay.sampling_per_actual=sampling_percentage;
Replay.Replay_score_actual=Replay_score_actual;
Replay.Replay_slope_actual=slope;
Replay.actual_jumpiness=actual_jumpiness;

Replay.sampling_percentage_position=sampling_percentage_position;
Replay.score_shuffle_position=rsq_shuffle_position;
Replay.slope_position=slope_position;
Replay.Shuffle_position_jumpiness=Shuffle_position_jumpiness;

Replay.sampling_percentage_time=sampling_percentage_time;
Replay.score_shuffle_time=rsq_shuffle_time;
Replay.slope_time=slope_time;
Replay.Shuffle_time_jumpiness=Shuffle_time_jumpiness;
Replay.decoding.REM_decoded_position=decoding.REM_decoded_position;

end
function [ms] = msExtractBinary_detrendTraces(ms,z_threshold);
%MSEXTRACTBINARY Summary of this function goes here
%   Detailed explanation goes here

%% Parameters
    z_threshold = 2;
 

    [bFilt,aFilt] = butter(2,  2/(30/2), 'low');

    for trace_i = 1:ms.numNeurons;
        detrend_raw_trace=detrend(ms.RawTraces(:,trace_i),2);  %2;  d removes the quadratic trend. but didn't work?
        filt_trace = zscore(filtfilt(bFilt,aFilt,detrend_raw_trace));

        d1_trace = diff(filt_trace);
        d1_trace(end+1) = 0;

        binary_trace = filt_trace*0;
        binary_trace(filt_trace>z_threshold & d1_trace>0) = 1;

        ms.detrendRaw(:,trace_i)=detrend_raw_trace;
        ms.Binary(:,trace_i) = binary_trace;
       
        

    end

 
%         ms.Binary=ms.Binary;
%         ms.Sum=sum(ms.Binary);
%         ms.mean=mean(ms.Sum);

end

% sum1_5=sum(Threshold_1_5);
% sum2=sum(Threshold_2);
% sum3=sum(Threshold_3);
% 
% mean1_5=mean(sum1_5);
% mean2=mean(sum2);
% mean3=mean(sum3);

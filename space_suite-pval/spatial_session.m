function spatial_analysis=spatial_session(ms, behav, varargin)
%session_space preps the behavioural session for further analysis
%INPUTS
%   ms calcium data structures
%   behav   struct with track info
%   cell_ID    cell#
% <OPTIONALS> 'params' structure 
% method   1 binarized trace, 2 raw trace, or 'both' default both  
% min_speed   Rt01425: 200. 
%endportion proportion of track at each end which will be excluded. default
%is .1 (10%)
%OUTPUT
%   spatial_analysis: Full spatial analysis for linear track, 
%   binarized spatial_analysis.bin and/or RAW spatial_analysis.raw

%Author: E Vico

%% parse the inputs
    p = inputParser;
    addParameter(p,'method','both')
    addParameter(p,'bin_size',4,@isnumeric)
    addParameter(p,'min_speed',5,@isnumeric)
    addParameter(p,'smoothing',true,@islogical) 
    addParameter(p,'place_field_threshold',.05,@isnumeric)
    addParameter(p,'min_place_field_stability',0.5,@isnumeric)
    addParameter(p,'numShuffles',1000,@isnumeric)
    addParameter(p, 'endportion',0.1,@isnumeric)
    
    parse(p,varargin{:})
    method = p.Results.method;
    bin_size = p.Results.bin_size;
    min_speed = p.Results.min_speed;
    bin_size = p.Results.bin_size;
    smoothing = p.Results.smoothing ;
    place_field_threshold = p.Results.place_field_threshold;
    min_place_field_stability = p.Results.min_place_field_stability;
    numShuffles = p.Results.numShuffles;
    endportion= p.Results.endportion;
    params=p.Results;
    params=rmfield(params,'method');
    params=rmfield(params,'endportion');
    
    if method==1
        am=1:1;
    elseif method==2
        am=2:2;
    else
        am=1:2; 
    end      
        
    % Binning euclidian space 1D
    pr_behav.X_bin_vector =  (endportion*behav.trackLength)+bin_size:bin_size: (behav.trackLength- endportion*behav.trackLength);%taking the edge bins out

    %% Make sure time vectors contain unique time values
    [ubehav, IAbehav, ~]=unique(behav.time);
    [ums, IAms, ~]=unique(ms.time);

   
    % Interpolate the position of the mouse on the linear track
    pr_behav.interpolated_X=interp1(ubehav,behav.position(IAbehav,1),ums); 
    pr_behav.interpolated_Y=interp1(ubehav,behav.position(IAbehav,2),ums);  
    pr_behav.interpolated_speed=interp1(ubehav,behav.speed(IAbehav),ums);
    pr_behav.background=behav.background;
    
    %calculate directions    
    rVector= direction_laps(pr_behav, 'r');
    lVector= direction_laps(pr_behav, 'l');
    
    %% Get binarized trace
    [ms] = msExtractBinary_detrendTraces(ms,2);
    
    %% Calculating spatial firing
    binned_spatial={};
    raw_spatial={};
    for analysis_method=am
        counter=1;
        for cell_ID=1:ms.numNeurons;             
            if analysis_method ==1
                trace = ms.Binary(IAms,cell_ID);                   
            else
                trace = ms.RawTraces(IAms,cell_ID);
                trace = zscore(trace);                     
            end                 

            right_field_data  = spatial_firing(trace, pr_behav, cell_ID, 'dVector', rVector, 'analysis_method', analysis_method, params);
            left_field_data  = spatial_firing(trace, pr_behav, cell_ID, 'dVector',  lVector,  'analysis_method',analysis_method, params);
            all_field_data  = spatial_firing(trace, pr_behav, cell_ID, 'analysis_method',analysis_method, params);     
                        
            %%
%             %plotting
            if right_field_data.IsPlaceCell==1 && left_field_data.IsPlaceCell==0 && counter<100
                plot_place(right_field_data, ms, pr_behav, behav);
                counter=counter+1;
                if analysis_method==1
                    name=sprintf('SF_right_cell%d.fig', cell_ID);
                else
                    name=sprintf('SF_raw_right_cell%d.fig', cell_ID);
                end
               % saveas(gcf, name); close
            elseif (right_field_data.IsPlaceCell==0 && left_field_data.IsPlaceCell==1) && counter<100
               plot_place(left_field_data, ms, pr_behav, behav);
                counter=counter+1;
                if analysis_method==1
                    name=sprintf('SF_left_cell%d.fig', cell_ID);
                  %  saveas(gcf, name); close
                end
            elseif all_field_data.IsPlaceCell==1 && counter<100
                plot_place(all_field_data, ms, pr_behav, behav);
                counter=counter+1;
                if analysis_method==1
                    name=sprintf('SF_bi_cell%d.fig', cell_ID);
                  %  saveas(gcf, name); close
                end
            end    
            
        
            if analysis_method ==1
                binned_spatial(cell_ID,:)={right_field_data, left_field_data, all_field_data};
            else
                raw_spatial(cell_ID,:)={right_field_data, left_field_data, all_field_data};
            end
        end
    end
    spatial_analysis.bin=binned_spatial;
    spatial_analysis.raw=raw_spatial;

end
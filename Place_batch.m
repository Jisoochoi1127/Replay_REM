function [PCs_properties] = Place_batch(PARAMS, fname, inter_dir)
%% wrapper function to allow for parallel generation place cells properties. 


        warning off
        load(fname)
        
        %load place cell output here
        warning on
        fprintf('Generating PCs for: %s   %s....', info.subject, info.session)
        tic
        
        % set the parameters for data
        PARAMS.data.ca_time= ms.time/1000;
        PARAMS.data.ca_data=ms.RawTraces ;
        PARAMS.data.behav_time=behav.time/1000;
        PARAMS.data.behav_vec=behav.position(:,1);
        PARAMS.data.num_surrogates=1000;
 
        [PCs_properties] = extract_place_cells(inter_dir, info, PARAMS, ms, behav);
        toc
        
        save([inter_dir filesep info.subject '_' info.session '_PCs.mat'], 'PCs_properties')
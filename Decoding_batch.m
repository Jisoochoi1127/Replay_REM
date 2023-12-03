function [decoding] = Decoding_batch(PARAMS, fname,fname_PC, decoding_dir)
%% wrapper function to allow for parallel generation of tuning curves. 


         warning off
        load(fname)
        load(fname_PC)
        
        %load place cell output here
        warning on
        fprintf('Generating TCs for: %s   %s....', info.subject, info.session)
        tic
        
        % set the parameters for data
        PARAMS.data.ca_time= ms.time/1000;
        PARAMS.data.ca_data=ms.RawTraces ;
        PARAMS.data.behav_time=behav.time/1000;
        PARAMS.data.behav_vec=behav.position(:,1);
        
        [decoding] = Bayesian_JS2(decoding_dir, info, PARAMS, ms, behav, all_binary_pre_REM, all_binary_post_REM,PCs_properties);
        toc
        
        save([decoding_dir filesep info.subject '_' info.session '_decoding.mat'], 'decoding')
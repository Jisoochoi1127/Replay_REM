%starting in animal direcctory
fs=dir();

for i=3:length(fs)
   
    cd (fs(i).name)
    tic
    %1st classification of folders. Will have to do hypno for more careful analysis
    %2 is sw, 3 rem, 1 beh, 4nouse/further process
    %first take only folders with H
    ca_names = dir('H*');
    folders = {ca_names([ca_names.isdir]).name};
    c_folders=folders;
    %we have to order the folders
    for f=1: length(folders)
        times{1,f} = sscanf(folders{f},'H%d_M%d_S%d*');
    end
    wrong_digit=cellfun(@(x) x<10, times, 'UniformOutput', false);
    wrong_digit_folder=times(logical(cellfun(@sum, wrong_digit)));
    replacing_ind=find(cellfun(@sum, wrong_digit));

    for ri=1: length (replacing_ind)
        if sum(wrong_digit{replacing_ind(ri)})==1;
            if wrong_digit{replacing_ind(ri)}(1)
                folders{replacing_ind(ri)} = sprintf('H0%d_M%d_S%d',wrong_digit_folder{ri}(1),wrong_digit_folder{ri}(2),wrong_digit_folder{ri}(3));
            elseif wrong_digit{replacing_ind(ri)}(2)
                folders{replacing_ind(ri)} = sprintf('H%d_M0%d_S%d',wrong_digit_folder{ri}(1),wrong_digit_folder{ri}(2),wrong_digit_folder{ri}(3));
            else
                folders{replacing_ind(ri)} = sprintf('H%d_M%d_S0%d',wrong_digit_folder{ri}(1),wrong_digit_folder{ri}(2),wrong_digit_folder{ri}(3));
            end
        elseif sum(wrong_digit{replacing_ind(ri)})==2;
             if wrong_digit{replacing_ind(ri)}(1) && wrong_digit{replacing_ind(ri)}(2)
                folders{replacing_ind(ri)} = sprintf('H0%d_M0%d_S%d',wrong_digit_folder{ri}(1),wrong_digit_folder{ri}(2),wrong_digit_folder{ri}(3));
             elseif wrong_digit{replacing_ind(ri)}(1) && wrong_digit{replacing_ind(ri)}(3)
                folders{replacing_ind(ri)} = sprintf('H0%d_M%d_S0%d',wrong_digit_folder{ri}(1),wrong_digit_folder{ri}(2),wrong_digit_folder{ri}(3));
             else
                folders{replacing_ind(ri)} = sprintf('H%d_M0%d_S0%d',wrong_digit_folder{ri}(1),wrong_digit_folder{ri}(2),wrong_digit_folder{ri}(3));
             end
        else
            folders{replacing_ind(ri)} = sprintf('H0%d_M0%d_S0%d',wrong_digit_folder{ri}(1),wrong_digit_folder{ri}(2),wrong_digit_folder{ri}(3));
        end
    end

    [~, correct_index]=sort(folders);
    correct_folders=c_folders(correct_index);
    folder_class=4+zeros(length(folders),1);
    for f=1: length(folders)
        if contains(correct_folders{f},'mix','IgnoreCase',true);
            folder_class(f)= 4;
        elseif contains(correct_folders{f},'swrem','IgnoreCase',true);
            folder_class(f)= 4;
        elseif contains(correct_folders{f},'remsw','IgnoreCase',true);
        folder_class(f)= 4;
        elseif contains(correct_folders{f},'rem','IgnoreCase',true);
           folder_class(f)= 3;
       elseif contains(correct_folders{f}, 'sw', 'IgnoreCase',true);
           folder_class(f)= 2;
       elseif contains(correct_folders{f}, 'l', 'IgnoreCase',true);
           folder_class(f)= 1;
       end
    end
    %%
    try
        load('ms.mat');
        %let's look at the behav segment
        inx=find(folder_class==1);
        load('behav.mat');
    catch
        cd ..
        clearvars -except i fs
        continue
    end
    %ms to pass
    msB=[];
    %correct TS
    tf=sum(ms.timestamps1(1:inx-1));  
    msB.time=ms.time(tf+1:tf+ms.timestamps1(inx))-ms.time(tf+1); %timestamps0 for non-beh

    %slice frames
    msB.RawTraces=ms.RawTraces(tf+1:tf+ms.timestamps1(inx), :);%timestamps0 for non-beh
    msB.numFrames=tf+1+ms.timestamps1(inx)-(tf+1);%timestamps0 for non-beh
    msB.numNeurons=ms.numNeurons;
    msB= msExtractBinary_detrendTraces(msB);
    %%
    %per cell loop
    
    try 
        cd ("spatial_analysis")
    catch
        mkdir spatial_analysis    
        cd ("spatial_analysis")
    end

    spatial_analysis=spatial_session(msB, behav);
    spatial_analysis.speed=mean(behav.speed);
    save("spatial_analysis.mat", 'spatial_analysis');
    cd ..
    toc
    cd ..
    clearvars -except i fs
end


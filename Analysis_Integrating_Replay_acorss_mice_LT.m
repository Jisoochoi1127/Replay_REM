%% Integrate Replay events across mice
%June 21,2021
%Done    create script to integrate whole replay events
%Done    add data day to select the day to integrate data.
%Done    add plotting section for whole replay events across mice
%        compare LTD1 and LTD5
%Done    Make separate script for HAT
%Done    Running HATD5, HATDSwitch as well

%% June, 30, 2021
%Done   Added num of replay, sum. 
%Done   modified figure for subplot
%   finished figure for Science paper fig 2. 

%%

Data_day='LTD1'
% For day 1
%for 1043
%replay_dir = ['/Users/jisoo/Dropbox (Williams Lab)/Jisoo/JisooProject2020/2020_Results_aftercutting/Decoding/pv1043' filesep Data_day];
replay_dir = ['/Users/jisoo/Dropbox (Williams Lab)/Jisoo/JisooProject2020/2020_Results_aftercutting/Decoding/pv1043' filesep Data_day];

cd(replay_dir)

load 1000shuffling_1s

Whole_score=Final_Replay_score;
Whole_slope=Final_Replay_slope;
Whole_location_occupancy=Final_Sig_location_replay;
Whole_avg_jump=Final_Replay_avg_jumpiness;
Whole_max_jump=Final_Replay_max_jumpiness;
Whole_extent=Final_Replay_extent;
Whole_num_replay=length(Final_Replay_score);

replay_dir = ['/Users/jisoo/Dropbox (Williams Lab)/Jisoo/JisooProject2020/2020_Results_aftercutting/Decoding/pv1060' filesep Data_day];
cd(replay_dir)

load 1000shuffling_1s

Whole_score=[Whole_score, Final_Replay_score];
Whole_slope=[Whole_slope, Final_Replay_slope];
Whole_location_occupancy=[Whole_location_occupancy,Final_Sig_location_replay];
Whole_avg_jump=[Whole_avg_jump,Final_Replay_avg_jumpiness];
Whole_max_jump=[Whole_max_jump,Final_Replay_max_jumpiness];
Whole_extent=[Whole_extent,Final_Replay_extent];
Whole_num_replay=[Whole_num_replay;length(Final_Replay_score)];

replay_dir = ['/Users/jisoo/Dropbox (Williams Lab)/Jisoo/JisooProject2020/2020_Results_aftercutting/Decoding/pv1069' filesep Data_day];
cd(replay_dir)

load 1000shuffling_1s

Whole_score=[Whole_score,Final_Replay_score];
Whole_slope=[Whole_slope,Final_Replay_slope];
Whole_location_occupancy=[Whole_location_occupancy,Final_Sig_location_replay];
Whole_avg_jump=[Whole_avg_jump,Final_Replay_avg_jumpiness];
Whole_max_jump=[Whole_max_jump,Final_Replay_max_jumpiness];
Whole_extent=[Whole_extent,Final_Replay_extent];
Whole_num_replay=[Whole_num_replay;length(Final_Replay_score)];
Whole_sum_num_replay=sum(Whole_num_replay);

if contains(Data_day,'LTD5')
replay_dir = ['/Users/jisoo/Dropbox (Williams Lab)/Jisoo/JisooProject2020/2020_Results_aftercutting/Decoding/pv1252' filesep Data_day];
cd(replay_dir)

load 1000shuffling_1s

Whole_score=[Whole_score,Final_Replay_score];
Whole_slope=[Whole_slope,Final_Replay_slope];
Whole_location_occupancy=[Whole_location_occupancy,Final_Sig_location_replay];
Whole_avg_jump=[Whole_avg_jump,Final_Replay_avg_jumpiness];
Whole_max_jump=[Whole_max_jump,Final_Replay_max_jumpiness];
Whole_extent=[Whole_extent,Final_Replay_extent];
Whole_num_replay=[Whole_num_replay;length(Final_Replay_score)];
Whole_sum_num_replay=sum(Whole_num_replay);

end
    
replay_dir = ['/Users/jisoo/Dropbox (Williams Lab)/Jisoo/JisooProject2020/2020_Results_aftercutting/Decoding/pv1254' filesep Data_day];
cd(replay_dir)

load 1000shuffling_1s

Whole_score=[Whole_score,Final_Replay_score];
Whole_slope=[Whole_slope,Final_Replay_slope];
Whole_location_occupancy=[Whole_location_occupancy,Final_Sig_location_replay];
Whole_avg_jump=[Whole_avg_jump,Final_Replay_avg_jumpiness];
Whole_max_jump=[Whole_max_jump,Final_Replay_max_jumpiness];
Whole_extent=[Whole_extent,Final_Replay_extent];
Whole_num_replay=[Whole_num_replay;length(Final_Replay_score)];
Whole_sum_num_replay=sum(Whole_num_replay);


%% Mean 

Location_occupancy_whole=cell2mat(Whole_location_occupancy);
mean_location_ocuupancy=mean(Location_occupancy_whole);
mean_score=mean(Whole_score);
mean_slope=mean(Whole_slope);
mean_avg_jump=mean(Whole_avg_jump);
mean_max_jump=mean(Whole_max_jump);
mean_extent=mean(Whole_extent);


%% plotting output
save_dir=['/Users/jisoo/Dropbox (Williams Lab)/Jisoo/JisooProject2020/2020_Results_aftercutting/Decoding/Replay_across_Mice' filesep Data_day];
cd(save_dir)

subplot(2,3,1);
histogram(Whole_score,'BinWidth',0.025,'EdgeAlpha',0.1,'FaceColor','k');
xlabel('Replay Score')
ylabel('Frequency')

subplot(2,3,2);
histogram(Whole_slope,'BinWidth',3, 'EdgeAlpha',0.5,'FaceColor','k');
xlabel('Slope')
ylabel('Frequency')


subplot(2,3,3);
histogram(Whole_extent,'BinWidth',5,'EdgeAlpha',0.1,'FaceColor','k');
xlabel('Whole extent')
ylabel('Frequency')

subplot(2,3,4);
histogram(Whole_avg_jump,'BinWidth',5,'EdgeAlpha',0.5,'FaceColor','k');
xlabel('Average Jumpiness')
ylabel('Frequency')

subplot(2,3,5);
histogram(Whole_max_jump,'BinWidth',5,'EdgeAlpha',0.5,'FaceColor','k');
xlabel('Max Jumpiness')
ylabel('Frequency')


subplot(2,3,6);
 Location_occupancy_whole=cell2mat(Whole_location_occupancy);
histogram( Location_occupancy_whole,'BinWidth',5,'EdgeAlpha',0.5,'FaceColor','k');
xlabel('Location occupancy')
ylabel('Frequency')

saveas(gcf,'Replay_across_mice_1s.jpg');

%% Output
Whole_Replay.Whole_score=Whole_score';
Whole_Replay.Whole_slope=Whole_slope';
Whole_Replay.Whole_location_occupancy=Whole_location_occupancy';
Whole_Replay.Whole_avg_jump=Whole_avg_jump';
Whole_Replay.Whole_max_jump=Whole_max_jump';
Whole_Replay.Whole_extent=Whole_extent';
Whole_Replay.Whole_num_replay=Whole_num_replay; %num of replay from individual
Whole_Replay.Whole_sum_num_replay=Whole_sum_num_replay;

Whole_Replay.Mean_Whole_score=mean_score
Whole_Replay.Mean_Whole_slope=mean_slope;
Whole_Replay.Mean_Whole_location_occupancy=mean_location_ocuupancy;
Whole_Replay.Mean_Whole_avg_jump=mean_avg_jump;
Whole_Replay.Mean_Whole_max_jump=mean_max_jump;
Whole_Replay.Mean_Whole_extent=mean_extent;



save ('Whole_Replay_1s.mat','Whole_Replay')

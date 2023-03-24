%% horizontal line plot

%dir='/Users/jisoo/Williams Lab Dropbox/Guillaume Etter/Jisoo/JisooProject2020/2020_Results_aftercutting/Decoding'
dir='/Users/jisoo/Williams Lab Dropbox/Guillaume Etter/Jisoo/JisooProject2020/2020_Results_aftercutting/Makr_comment'
id='pv1254'
data='HATDSwitch'


path=[dir filesep id filesep data];
cd(path)
load 1000shuffling_1s

x=Final_Sig_location_replay;
y = 1:length(x);
z = Final_Replay_extent;

%c_ord = linspecer(length(x));  
c_ord = parula(length(x));  


[~, s_idx] = sort(abs(z));
x = x(:,s_idx);



for k=1:length(Final_Replay_score);
    
    location=Final_Sig_location_replay{k};
    left=length(find(location<50))/length(location);
    right=length(find(location>=50))/length(location);
    
    % other bias indexing by checking the distance
    
   a=location(end)-50;
   b=location(1)-50;
   
    %for pv1043,1060,1192- HATD1,5 - left ; open 
    % for pv1191- HATD1 - right ; open 

    
    %for pv1043,1060,1192- HATDS - right; open 
    
    %Bias_idx(ii)=left-right; %open-closed arm 
    
    if contains(id,'pv1191') | contains(data,'HATDSwitch')

    Bias_idx(k)= right-left;   %open-closed arm for pv1191_HATd1
    Bias_distance(k)=(a-b)/100; % if it's positive; right side bias

    else
        
         Bias_idx(k)= left-right;
         Bias_distance(k)=(a-b)/100; 
    end
   
    
end
figure

for ii = 1:length(x)
   % hold on
    
    location=x{ii}
    y=[ii*ones(1,length(location))];
    
%     left=length(find(location<50))/length(location);
%     right=length(find(location>=50))/length(location);
%     
%     % other bias indexing by checking the distance
%     
%    a=location(end)-50;
%    b=location(1)-50;
%    
   
    
    % for pv1043,1060,1192- HATD1,5 - left ; open 
    % for pv1191- HATD1 - right ; open 

    
    %for pv1043,1060,1192- HATDS - right; open 
    
    %Bias_idx(ii)=left-right; %open-closed arm 
    
%     if contains(id,'pv1191') | contains(data,'HATDSwitch')
% 
%     Bias_idx(ii)= right-left;   %open-closed arm for pv1191_HATd1
%     Bias_distance(ii)=(a-b)/100; % if it's positive; right side bias
% 
%     else
%         
%          Bias_idx(ii)= left-right;
%          Bias_distance(ii)=(a-b)/100; 
%     end
    
%plot(x{ii}, y(:,ii), '-', 'color', c_ord(:,ii))
plot(location, y, '-','LineWidth',1)
xlabel('Location of replay');
ylabel('Number of replay');
xlim([0 100]);

hold on


end
box off
Average_bias_ratio=mean(Bias_idx);
Average_bias_distance=mean(Bias_distance);
%% Out put
Bias.Bias_idx=Bias_idx';
Bias.Average_bias_ratio=Average_bias_ratio;

save Bias_idx
saveas(gcf,'Replay_bias.jpg');
save('Bias.mat','Bias')



%Replay plotting

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
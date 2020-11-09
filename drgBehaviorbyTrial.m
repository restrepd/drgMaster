function drgBehaviorbyTrial(handles)


[perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);


try
    close 1
catch
end

hFig1 = figure(1);

%Plot the percent correct
% subplot(3,1,1)
trials=1:length(perCorr);

%Plot in different colors
plot(trials,perCorr,'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
hold on

if sum(encoding_trials)>0
    plot(trials(encoding_trials),perCorr(encoding_trials),'o','MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0 114/255 178/255] )
end

if sum(retrieval_trials)>0
    plot(trials(retrieval_trials),perCorr(retrieval_trials),'o','MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [158/255 31/255 99/255] )
end

set(gca, 'box', 'off')

%Plot black
% plot(trials,perCorr,'ok','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7)
% set(gca,'FontSize',20,'LineWidth',1.25)


ylim([0 110]);
xlabel('Trial No for Odor On events')
ylabel('Percent correct')
title('% correct vs trial number, blue=encoding, red=retrieval')





function drgBehaviorbyTrial(handles)


[perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);


try
    close 1
catch
end

hFig1 = figure(1);
% set(hFig1, 'units','normalized','position',[.25 .25 .23 .65])

%Plot the percent correct
% subplot(3,1,1)
trials=1:length(perCorr);

%Plot in different colors
plot(trials,perCorr,'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
hold on
plot(trials(encoding_trials),perCorr(encoding_trials),'ob')
plot(trials(retrieval_trials),perCorr(retrieval_trials),'or')

%Plot black
% plot(trials,perCorr,'ok','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7)
% set(gca,'FontSize',20,'LineWidth',1.25)


ylim([0 110]);
xlabel('Trial No for Odor On events')
ylabel('Percent correct')
title('% correct vs trial number, blue=encoding, red=retrieval')





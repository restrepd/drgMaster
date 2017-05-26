function drgBehaviorbyBlock(handles)
disp('justin debug');

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

n = 20;
blperf = perCorr(1:n:length(perCorr));
nobl = 1:length(handles.drg.draq_d.blocks);
plot(1:length(blperf),blperf,'-or','MarkerEdgeColor','r');

% plot(trials,perCorr,'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
hold on
% plot(trials(encoding_trials),perCorr(encoding_trials),'ob')
% plot(trials(retrieval_trials),perCorr(retrieval_trials),'or')
ylim([40 100]);
exnobl = nobl(end)+1;
xlim([0 exnobl]);
xlabel('Block number')
ylabel('Percent correct')
title('% correct vs block number')






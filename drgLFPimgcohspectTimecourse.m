function drgLFPimgcohspectTimecourse(handles)
%Generates a timecourse of the imaginary coherence between two LFP channels

[t,f, all_Cxy_timecourse, this_trialNo]=drgGetLFPimgCoherenceForThisEvTypeNo(handles);


freq=f';

Cxy_timecourse(:,:)=mean(all_Cxy_timecourse,1);

if handles.autoscale==1
    maxCxy=prctile(Cxy_timecourse(:),99);
    minCxy=prctile(Cxy_timecourse(:),1);
else
    maxCxy=handles.maxLogP;
    minCxy=handles.minLogP;
end

if minCxy==maxCxy
    minCxy=maxCxy-0.01;
end

%Plot the timecourse
try
    close(1)
catch
end
hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.07 .1 .75 .3])


drg_pcolor(repmat(t',length(freq),1)',repmat(f,length(t),1),Cxy_timecourse(1:length(freq),1:length(t))')


colormap fire
shading interp
caxis([minCxy maxCxy]);
xlabel('Time (sec)')
ylabel('Frequency (Hz)');
title(['Coherence timecourse ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])

try
    close 2
catch
end

hFig2 = figure(2);
set(hFig2, 'units','normalized','position',[.83 .1 .05 .3])

prain=[minCxy:(maxCxy-minCxy)/99:maxCxy];
drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
colormap fire
shading interp
ax=gca;
set(ax,'XTickLabel','')
% end

%This code is here for Daniels' Figure 1
%Plot the timecourse
try
    close(3)
catch
end
hFig3 = figure(3);
set(hFig3, 'units','normalized','position',[.07 .1 .75 .3])
plot(t',mean(Cxy_timecourse((f>=6)&(f<=14),1:length(t)))','-k','LineWidth',3)
ylim([0 1])
xlabel('Time (sec)')
ylabel('Theta coherence')
pffft=1

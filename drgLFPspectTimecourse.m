function drgLFPspectTimecourse(handles)
%Generates a timecourse of the LFP power in decibels 10*log10(Power)



[t,f,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo]=drgGetLFPPowerForThisEvTypeNo(handles);



freq=handles.burstLowF:1:handles.burstHighF;


%Timecourse doing average after log
%Get max and min
if handles.subtractRef==0
    log_P_timecourse=zeros(length(f),length(t));
    log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),1);
    if handles.autoscale==1
        maxLogP=prctile(log_P_timecourse(:),99);
        minLogP=prctile(log_P_timecourse(:),1);
    else
        maxLogP=handles.maxLogP;
        minLogP=handles.minLogP;
    end
else
    log_P_timecourse=zeros(length(f),length(t));
    log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),1);
    log_P_timecourse_ref=zeros(length(f),length(t));
    log_P_timecourse_ref(:,:)=repmat(mean(10*log10(all_Power_ref),1)',1,length(t));
    if handles.autoscale==1
        maxLogP=prctile(log_P_timecourse(:)-log_P_timecourse_ref(:),99);
        minLogP=prctile(log_P_timecourse(:)-log_P_timecourse_ref(:),1);
    else
        maxLogP=handles.maxLogP;
        minLogP=handles.minLogP;
    end
end

%Note: Diego added this on purpose to limit the range to 10 dB
%This results in emphasizing changes in the top 10 dB
if maxLogP-minLogP>10
    minLogP=maxLogP-10;
end

try
    close 1
catch
end

%Plot the timecourse
hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.07 .1 .75 .3])

%     if handles.subtractRef==0
%         pcolor(repmat(t,length(freq),1)',repmat(freq,length(t),1),P_timecourse')
%     else
%         pcolor(repmat(t,length(freq),1)',repmat(freq,length(t),1),P_timecourse'-P_timecourse_ref')
%     end

if handles.subtractRef==0
    drg_pcolor(repmat(t,length(freq),1)',repmat(freq,length(t),1),log_P_timecourse')
else
    %pcolor(repmat(t,length(f),1)',repmat(f,length(t),1),10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
    drg_pcolor(repmat(t,length(freq),1)',repmat(freq,length(t),1),log_P_timecourse'-log_P_timecourse_ref')
    %imagesc(t,f,10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
end

colormap jet
shading interp
caxis([minLogP maxLogP]);
xlabel('Time (sec)')
ylabel('Frequency (Hz)');
title(['Power (dB) timecourse ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])

try
    close 2
catch
end

hFig2 = figure(2);
set(hFig2, 'units','normalized','position',[.83 .1 .05 .3])

prain=[minLogP:(maxLogP-minLogP)/99:maxLogP];
drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
colormap jet
shading interp
ax=gca;
set(ax,'XTickLabel','')
% end
pffft=1

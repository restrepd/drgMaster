function drgLFPxspectTimecourse(handles)
%Generates a timecourse of the LFP power in decibels 10*log10(Power)


[t,f,all_xspec_ref, all_xspec_timecourse]=drgGetLFPxspectrogramForThisEvTypeNo(handles);



freq=f';


%Timecourse doing average after log
%Get max and min
if handles.subtractRef==0
    log_xspec_timecourse=zeros(length(f),length(t));
    log_xspec_timecourse(:,:)=mean(20*log10(all_xspec_timecourse),1);
    if handles.autoscale==1
        maxLog_xspec=prctile(log_xspec_timecourse(:),99);
        minLog_xspec=prctile(log_xspec_timecourse(:),1);
    else
        maxLog_xspec=handles.maxLogP;
        minLog_xspec=handles.minLogP;
    end
else
    log_xspec_timecourse=zeros(length(f),length(t));
    log_xspec_timecourse(:,:)=mean(20*log10(all_xspec_timecourse),1);
    log_xspec_timecourse_ref=zeros(length(f),length(t));
    log_xspec_timecourse_ref(:,:)=repmat(mean(20*log10(all_xspec_ref),1)',1,length(t));
    if handles.autoscale==1
        maxLog_xspec=prctile(log_xspec_timecourse(:)-log_xspec_timecourse_ref(:),99);
        minLog_xspec=prctile(log_xspec_timecourse(:)-log_xspec_timecourse_ref(:),1);
    else
        maxLog_xspec=handles.maxLogP;
        minLog_xspec=handles.minLogP;
    end
end
 
%Note: Diego added this on purpose to limit the range to 10 dB
%This results in emphasizing changes in the top 10 dB
if handles.autoscale==1
    if maxLog_xspec-minLog_xspec>10
        minLog_xspec=maxLog_xspec-10;
    end
end

try
    close 1
catch
end

%Plot the timecourse
hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.07 .1 .75 .3])

if handles.subtractRef==0
    drg_pcolor(repmat(t,length(freq),1)',repmat(freq',length(t),1),log_xspec_timecourse')
else
    drg_pcolor(repmat(t,length(freq),1)',repmat(freq',length(t),1),log_xspec_timecourse'-log_xspec_timecourse_ref')
end

colormap jet
shading interp
caxis([minLog_xspec maxLog_xspec]);
xlabel('Time (sec)')
ylabel('Frequency (Hz)');
title(['Cross-spectrogram for power (dB) timecourse for ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo} ' LFP' num2str(handles.peakLFPNo) 'xLFP' num2str(handles.burstLFPNo) ])

try
    close 2
catch
end

hFig2 = figure(2);
set(hFig2, 'units','normalized','position',[.83 .1 .05 .3])

prain=[minLog_xspec:(maxLog_xspec-minLog_xspec)/99:maxLog_xspec];
drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
colormap jet
shading interp
ax=gca;
set(ax,'XTickLabel','')


function drgLFPwaveTimecourse(handles)
%Generates a timecourse of the LFP power in decibels 10*log10(Power)


[t,freq,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo]=drgGetLFPwavePowerForThisEvTypeNo(handles);

%Timecourse doing average after log
%Get max and min
if handles.subtractRef==0
    log_P_timecourse=zeros(length(freq),length(t));
    log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),1);
    
    %Per trial power plot
    log_P_per_trial_timecourse=zeros(length(freq)*length(this_trialNo),length(t));
    y_shift=0;
    for trialNo=1:length(this_trialNo)
        this_log_P_timecourse=zeros(length(freq),length(t));
        this_log_P_timecourse(:,:)=10*log10(all_Power_timecourse(trialNo,:,:));
        log_P_per_trial_timecourse(y_shift+1:y_shift+length(freq),:)=this_log_P_timecourse;
        y_shift=y_shift+length(freq);
    end
    
    if handles.autoscale==1
        maxLogPper=prctile(log_P_timecourse(:),99);
        minLogPper=prctile(log_P_timecourse(:),1);
        %Note: Diego added this on purpose to limit the range to 10 dB
        %This results in emphasizing changes in the top 10 dB
        if maxLogPper-minLogPper>12
            minLogPper=maxLogPper-12;
        end
    else
        maxLogPper=handles.maxLogP;
        minLogPper=handles.minLogP;
    end
else
    log_P_timecourse=zeros(length(freq),length(t));
    log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),1);
    log_P_timecourse_ref=zeros(length(freq),length(t));
    log_P_timecourse_ref(:,:)=repmat(mean(10*log10(all_Power_ref),1)',1,length(t));
    
    %Per trial power plot
    log_P_per_trial_timecourse_sub=zeros(length(freq)*length(this_trialNo),length(t));
    y_shift=0;
    sy_shift=0;
    shifted_freq=[];
    for trialNo=1:length(this_trialNo)
        this_log_P_timecourse=zeros(length(freq),length(t));
        this_log_P_timecourse(:,:)=10*log10(all_Power_timecourse(trialNo,:,:));
        this_log_P_timecourse_ref=zeros(length(freq),length(t));
        this_log_P_timecourse_ref(:,:)=repmat(mean(10*log10(all_Power_ref(trialNo,:)),1)',1,length(t));
        log_P_per_trial_timecourse_sub(y_shift+1:y_shift+length(freq),:)=this_log_P_timecourse-this_log_P_timecourse_ref;
        shifted_freq(1,y_shift+1:y_shift+length(freq))=freq+(trialNo-1)*freq(end);
        y_shift=y_shift+length(freq);
    end
    
    max_delta=16;
    if handles.autoscale==1

        deltaLogP=log_P_timecourse'-log_P_timecourse_ref';
        maxLogPper=prctile(deltaLogP(:),99);
        minLogPper=prctile(deltaLogP(:),1);
        %Note: Diego added this on purpose to limit the range to 10 dB
        %This results in emphasizing changes in the top 10 dB
        if maxLogPper-minLogPper>max_delta
            minLogPper=maxLogPper-max_delta;
        end
        
    else
        maxLogPper=handles.maxLogP;
        minLogPper=handles.minLogP;
    end
end
 
%Calculate the licks
[lick_freq,times_lick_freq,lick_traces,CIlickf,lick_trace_times,stamped_lick_ii,these_stamped_lick_times,no_trials,trials_included]=drgGetLicks(handles);


if ~isempty(this_trialNo)
    
    for ii=1:12
        try
            close(ii)
        catch
        end
    end
    
    try
        close 1
    catch
    end
    
    %Plot the timecourse
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.07 .05 .75 .3])
    if handles.subtractRef==0
        drg_pcolor(repmat(t,length(freq),1)',repmat(freq,length(t),1),log_P_timecourse')
    else
        %pcolor(repmat(t,length(f),1)',repmat(f,length(t),1),10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
        drg_pcolor(repmat(t,length(freq),1)',repmat(freq,length(t),1),log_P_timecourse'-log_P_timecourse_ref')
        %imagesc(t,f,10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
    end
    
    colormap jet
    shading interp
    caxis([minLogPper maxLogPper]);
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)');
    title(['Power (dB, wavelet) timecourse ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
    
    %If this is a single trial show the licks
    if length(this_trialNo)==1
        
       
        %PLot the licks
        hold on
        
        %Clear the top strip
        ffrom=freq(1)+0.95*(freq(end)-freq(1));
        
        
        for t1=t(1):t(2)-t(1):t(end)
            plot([t1 t1],[ffrom freq(end)],'-w','LineWidth',3)
        end
        
        
        for ilick=1:stamped_lick_ii
            plot([these_stamped_lick_times(ilick) these_stamped_lick_times(ilick)],[ffrom freq(end)],'-k','LineWidth',3)
        end
    end
    
    try
        close 2
    catch
    end
    
    hFig2 = figure(2);
    set(hFig2, 'units','normalized','position',[.83 .1 .05 .3])
    
    prain=[minLogPper:(maxLogPper-minLogPper)/99:maxLogPper];
    drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    colormap jet
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')
    
    
    try
        close 5
    catch
    end
    
    %Plot the per-trial timecourse
    hFig5 = figure(5);
    set(hFig5, 'units','normalized','position',[.07 .1 .75 .3])
    
    if handles.subtractRef==0
        drg_pcolor(repmat(t,length(freq)*length(this_trialNo),1)',repmat(shifted_freq,length(t),1),log_P_per_trial_timecourse')
    else
        %pcolor(repmat(t,length(f),1)',repmat(f,length(t),1),10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
        drg_pcolor(repmat(t,length(freq)*length(this_trialNo),1)',repmat(shifted_freq,length(t),1),log_P_per_trial_timecourse_sub')
        %imagesc(t,f,10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
    end
    
    colormap jet
    shading interp
    caxis([minLogPper maxLogPper]);
    xlabel('Time (sec)')
    ylabel('Frequency*trialNo');
    title(['Power (dB, wavelet) timecourse per trial ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
    
    try
        close 6
    catch
    end
    
    hFig6 = figure(6);
    set(hFig6, 'units','normalized','position',[.83 .1 .05 .3])
    
    prain=[minLogPper:(maxLogPper-minLogPper)/99:maxLogPper];
    drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    colormap jet
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')
    
    %Plot the lick traces
    try
        close 7
    catch
    end
    
    hFig7 = figure(7);
    set(hFig7, 'units','normalized','position',[.07 .4 .75 .3])
    
    hold on
    
    per99=prctile(lick_traces(:),99.9);
    per1=prctile(lick_traces(:),1);
    
    mean_licks=zeros(1,length(lick_trace_times));
    
    
    y_shift=0;

    %Plot lick traces
    for ii=1:no_trials
        plot(lick_trace_times,lick_traces(ii,:)+y_shift,'-r')
        y_shift=y_shift+1.5*(per99-per1);
    end
    
    y_shift=y_shift+1.5*(per99-per1);
    ylim([0 y_shift])
    xlim([t(1) t(end)])
    xlabel('time(sec)')
    title('Lick traces')
    
    
     %Plot the lick frequency
    try
        close 8
    catch
    end
     
    hFig8 = figure(8);
    set(hFig8, 'units','normalized','position',[.07 .7 .75 .3])
    
    
    if length(this_trialNo)==1
        plot(times_lick_freq, lick_freq)
    else
        [hl1, hp1] = boundedline(times_lick_freq',lick_freq', CIlickf', 'r');
    end
    
    ylim([0 1.2*max(lick_freq)])
    xlabel('Time (sec)')
    ylabel('frequency (Hz)')
    title('Lick frequency')
    
    figure(7)
    figure(8)
    figure(1)
    figure(2)
end

pfft=1;



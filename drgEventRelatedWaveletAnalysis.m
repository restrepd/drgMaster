function [log_P_t,no_trials_w_event,which_event,f,out_times,times,phase_per_trial,no_trials,no_events_per_trial,t_per_event_per_trial,trial_map,perCorrERP,no_ref_evs_per_trial]=drgEventRelatedWaveletAnalysis(handles)
%Performs an event-related analysis. The event is signaled by a sharp chane
%in the reference voltage. This is used to analyze lick-related changes in
%LFP

[perCorr, encoding_trials, retrieval_trials,encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);

odorOn=2;

if handles.displayData==1
    fprintf(1, '\n');
end

which_event=[];
angleLFP = [];
no_events_per_trial=[];
no_ref_events_per_trial=[];
t_per_event_per_trial=[];
time_per_event=[];
perCorrERP=[];


%Generates a trial per trial phase histogram
sessionNo=handles.sessionNo;
Fs=floor(handles.drg.session(sessionNo).draq_p.ActualRate);
lowF1=handles.peakLowF;
lowF2=handles.peakHighF;
highF1=handles.burstLowF;
highF2=handles.burstHighF;
pad_time=handles.time_pad;
n_phase_bins=handles.n_phase_bins;

dec_n=fix(handles.drg.session(sessionNo).draq_p.ActualRate/1000);

%Setup the wavelet scales
%   scales = helperCWTTimeFreqVector(minfreq,maxfreq,f0,dt,NumVoices)
%   f0 - center frequency of the wavelet in cycles/unit time
%   dt - sampling interval
%   NumVoices - number of voices per octave

NumVoices=5;
minfreq=handles.burstLowF;
maxfreq=handles.burstHighF;
dt=1/Fs;
f0=5/(2*pi);

a0 = 2^(1/NumVoices);
minscale = f0/(maxfreq*dt);
maxscale = f0/(minfreq*dt);
minscale = floor(NumVoices*log2(minscale));
maxscale = ceil(NumVoices*log2(maxscale));
scales = a0.^(minscale:maxscale).*dt;

window=round(handles.window*handles.drg.draq_p.ActualRate);
noverlap=round(0.975*handles.window*handles.drg.draq_p.ActualRate);

no_time_pts=floor(handles.window*handles.drg.session(sessionNo).draq_p.ActualRate)+1;
times=[1:no_time_pts]/handles.drg.session(sessionNo).draq_p.ActualRate;
times=times-(handles.window/2);



%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;

%Calculate the threshold value to detect a lick
all_refs=[];
for trNo=firstTr:lastTr
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            %Note: handles.peakLFPNo is the reference LFP
            [referenceLFP, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.startRef, handles.time_end);
             
            if (can_read1==1)
                all_refs=[all_refs referenceLFP];
            end
        end
    end
end

thershold_ref=prctile(all_refs,1)+((prctile(all_refs,99)-prctile(all_refs,1))/2);

%First find the time range for the spectrogram
if handles.subtractRef==1
    if handles.time_start+handles.time_pad<handles.startRef+handles.time_pad
        min_t=handles.time_start+handles.time_pad-(handles.window/2);
    else
        min_t=handles.startRef+handles.time_pad-(handles.window/2);
    end
    
    if handles.time_end-handles.time_pad>handles.endRef-handles.time_pad
        max_t=handles.time_end-handles.time_pad+handles.window;
    else
        max_t=handles.endRef-handles.time_pad+handles.window;
    end
else
    min_t=handles.time_start+handles.time_pad-(handles.window/2);
    max_t=handles.time_end-handles.time_pad+handles.window;
end


%Now get the LFP phase of the events

events=[];
phase=[];
time=[];
no_events=0;
no_trials=0;
no_trials_w_event=0;

ERLFP=[];
phase_per_trial=[];
trial_map=[];

log_P_t=[];

inter_lick_intervals_ref=[];
ii_ili_ref=0;

inter_lick_intervals=[];
ii_ili=0;

for trNo=firstTr:lastTr
    
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            %Note: handles.peakLFPNo is the reference LFP
            [referenceLFP, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, min_t, max_t);
            [LFP, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, min_t, max_t);
            
            if (can_read1==1)&(can_read2==1)
                no_trials=no_trials+1;
                trial_map(no_trials)=trNo;
                no_events_per_trial(no_trials)=0;
                no_ref_events_per_trial(no_trials)=0;
                t_per_event_per_trial(no_trials,1)=0;
                time(no_trials)=handles.drg.session(sessionNo).trial_start(trialNo);
                which_trial(no_trials)=1;
                perCorr_per_histo(no_trials)=50;
                perCorrERP(no_trials)=perCorr(drgFindEvNo(handles,trialNo,sessionNo,odorOn));

                
                
                %Get LFP phase
                bpFiltLFP = designfilt('bandpassiir','FilterOrder',20, ...
                    'HalfPowerFrequency1',highF1,'HalfPowerFrequency2',highF2, ...
                    'SampleRate',Fs);
                thfiltLFP=filtfilt(bpFiltLFP,LFP);
                thisangleLFP = angle(hilbert(thfiltLFP)); % LFP phase
                angleLFP = [angleLFP thisangleLFP];
                
                
                %Now do the wavelet transform
                decLFP=decimate(LFP,dec_n);
                decFs=Fs/dec_n;
               
                cwtLFP = cwtft({detrend(double(decLFP)),1/decFs},'wavelet','morl','scales',scales);
                Prev=abs(cwtLFP.cfs).^2;
                P=Prev(end:-1:1,:);
                DT=1/decFs;
                t = 0:DT:(numel(decLFP)*DT)-DT;
                frev=cwtLFP.frequencies;
                f=frev(end:-1:1)';
                
                %Calculate the wings for the calculation of the event-triggered spectrogram
                times_spec=t+min_t;
                wing_ii=int64((handles.window/2)/(times_spec(2)-times_spec(1)));
                dt= times_spec(2)- times_spec(1);
                out_times=double([-wing_ii:wing_ii])*dt;
                
                %Get events in the reference time window if nescessary
                if handles.subtractRef==1
                    
                    %Trim off the time pads
                    ii_start_ref=floor((handles.window/2)*handles.drg.session(sessionNo).draq_p.ActualRate);
                    delta_ii_end_ref=floor((handles.window)*handles.drg.session(sessionNo).draq_p.ActualRate);
                    end_ref=floor((handles.endRef-handles.time_pad-(handles.startRef-handles.time_pad)+(handles.window/2))*handles.drg.session(sessionNo).draq_p.ActualRate);
                    refLFP_ref=referenceLFP(ii_start_ref:end_ref-delta_ii_end_ref-1);
                    thaLFP_ref=thisangleLFP(ii_start_ref:end_ref-delta_ii_end_ref-1);
                    
                    if refLFP_ref(1)>thershold_ref
                        ii=find(refLFP_ref<thershold_ref,1,'first');
                    else
                        ii=1;
                    end
                    
                    the_end=0;
                    
                    ref_power_these_events=[];
                    no_ref_evs_this_trial=0;
                    
                    this_lick_ii=0;
                    these_lick_times=[];
                    
                    %Find the events (licks)
                    while the_end==0
                        next_event=find(refLFP_ref(ii:end)>thershold_ref,1,'first');
                        if isempty(next_event)
                            the_end=1;
                        else
                            
                            ii=ii+next_event-1;
                            
                            %Exclude if the inter event interval is too
                            %small due to noise in the lick signal
                            
                            %Find the inter lick interval
                            this_lick_ii=this_lick_ii+1;
                            these_lick_times(this_lick_ii)=(ii/handles.drg.session(sessionNo).draq_p.ActualRate);
                            if this_lick_ii>1
                                %record the inter lick interval
                                ii_ili_ref=ii_ili_ref+1;
                                inter_lick_intervals_ref(ii_ili_ref)=these_lick_times(this_lick_ii)-these_lick_times(this_lick_ii-1);
                            end
                            
                            %Enter the event (lick) in the timecourse only if it is
                            %not within a burst of high frequency noise
                            if ii_ili_ref>0
                                if inter_lick_intervals_ref(ii_ili_ref)>handles.smallest_inter_lick_interval
                                    
                                    %Make sure that the array is large enough
                                    this_time=handles.startRef+pad_time+(ii/handles.drg.session(sessionNo).draq_p.ActualRate);
                                    [mint,mint_ii]=min(abs(times_spec-this_time));
                                    
                                    if (mint_ii+wing_ii<=length(times_spec))&(mint_ii-wing_ii>=1)
                                        no_ref_evs_this_trial=no_ref_evs_this_trial+1;
                                        no_ref_events_per_trial(no_trials)=no_ref_evs_this_trial;
                                        lot=length(out_times);
                                        ref_Power_these_events(no_ref_evs_this_trial,1:length(f),1:length(out_times))=P(:,mint_ii-wing_ii:mint_ii+wing_ii);
                                    end
                                    
                                end
                            end
                            
                            end_event=find(refLFP_ref(ii:end)<thershold_ref,1,'first');
                            if isempty(end_event)
                                the_end=1;
                            else
                                ii=ii+end_event-1;
                            end
                        end
                    end 
                else
                    no_ref_evs_this_trial=1;
                end
                
                no_ref_evs_per_trial(no_trials)=no_ref_evs_this_trial;
                
                %Get events
                
                %Trim off the time pads
                ii_start=floor(((handles.time_start-handles.startRef)+(handles.window/2))*handles.drg.session(sessionNo).draq_p.ActualRate);
                delta_ii_end=floor((handles.window)*handles.drg.session(sessionNo).draq_p.ActualRate);
                ref=referenceLFP(ii_start:end-delta_ii_end-1);
                thaLFP=thisangleLFP(ii_start:end-delta_ii_end-1);
                
                if ref(1)>thershold_ref
                    ii=find(ref<thershold_ref,1,'first');
                else
                    ii=1;
                end
                
                the_end=0;
                
                phase_this_trial=[];
                all_Power_these_events=[];
                ERLFP_this_trial=[];
                no_evs_this_trial=0;
                
                
                this_lick_ii=0;
                these_lick_times=[];
                
                while the_end==0
                    next_event=find(ref(ii:end)>thershold_ref,1,'first');
                    if isempty(next_event)
                        the_end=1;
                    else
                        
                        ii=ii+next_event-1;
                        
                        
                        %Exclude if the inter event interval is too
                        %small due to noise in the lick signal
                        
                        %Find the inter lick interval
                        this_lick_ii=this_lick_ii+1;
                        these_lick_times(this_lick_ii)=(ii/handles.drg.session(sessionNo).draq_p.ActualRate);
                        if this_lick_ii>1
                            %record the inter lick interval
                            ii_ili=ii_ili+1;
                            inter_lick_intervals(ii_ili)=these_lick_times(this_lick_ii)-these_lick_times(this_lick_ii-1);
                        end
                        
                        %Enter the event (lick) in the timecourse only if it is
                        %not within a burst of high frequency noise
                        if ii_ili>0
                            if inter_lick_intervals(ii_ili)>handles.smallest_inter_lick_interval
                                
                                %Make sure that the array is large enough
                                this_time=handles.time_start+pad_time+(ii/handles.drg.session(sessionNo).draq_p.ActualRate);
                                [mint,mint_ii]=min(abs(times_spec-this_time));
                                
                                if (mint_ii+wing_ii<=length(times_spec))&(mint_ii-wing_ii>=1)
                                    no_events=no_events+1;
                                    events(no_events)=ii;
                                    time_per_event(no_events)=handles.time_start+pad_time+(ii/handles.drg.session(sessionNo).draq_p.ActualRate);
                                    phase(no_events)=thaLFP(ii);
                                    no_evs_this_trial=no_evs_this_trial+1;
                                    phase_this_trial(no_evs_this_trial)=thaLFP(ii);
                                    
                                    t_per_event_per_trial(no_trials,no_evs_this_trial)=time_per_event(no_events);
                                    
                                    ERLFP(no_events,:)=LFP(1,floor(ii_start+ii-(handles.window/2)*handles.drg.session(sessionNo).draq_p.ActualRate):...
                                        floor(ii_start+ii+(handles.window/2)*handles.drg.session(sessionNo).draq_p.ActualRate));
                                    
                                    ERLFP_this_trial(no_evs_this_trial,:)=LFP(1,floor(ii_start+ii-(handles.window/2)*handles.drg.session(sessionNo).draq_p.ActualRate):...
                                        floor(ii_start+ii+(handles.window/2)*handles.drg.session(sessionNo).draq_p.ActualRate));
                                    
                                    lot=length(out_times);
                                    all_Power_these_events(no_evs_this_trial,1:length(f),1:length(out_times))=P(:,mint_ii-wing_ii:mint_ii+wing_ii);
                                    
                                end
                                
                            end
                        end
                        
                        end_event=find(ref(ii:end)<thershold_ref,1,'first');
                        if isempty(end_event)
                            the_end=1;
                        else
                            ii=ii+end_event-1;
                        end
                    end
                end

                no_events_per_trial(no_trials)=no_evs_this_trial;
                
                if (no_evs_this_trial>0)&(no_evs_this_trial<handles.max_events_per_sec*(handles.time_end-handles.time_start-2*pad_time))...
                        &(no_ref_evs_this_trial>0)&(no_ref_evs_this_trial<handles.max_events_per_sec*(handles.time_end-handles.time_start-2*pad_time))
                    no_trials_w_event=no_trials_w_event+1;
                    ERLFP_per_trial(no_trials_w_event,:)=mean(ERLFP_this_trial,1);
                    %Per trial event related spectrogram
                    %Event-related spectrogram
                    %Timecourse doing average after log
                    %Get max and min
                    if handles.subtractRef==0
                        log_P_timecourse=zeros(length(f),length(out_times));
                        log_P_timecourse(:,:)=mean(10*log10(all_Power_these_events),1);
                        log_P_t(no_trials,1:length(f),1:length(out_times))=log_P_timecourse(:,:);
                    else
                        log_P_timecourse=zeros(length(f),length(out_times));
                        log_P_timecourse(:,:)=mean(10*log10(all_Power_these_events),1);
                        log_P_timecourse_ref=zeros(length(f),length(out_times));
                        log_P_timecourse_ref(:,:)=mean(10*log10(ref_Power_these_events),1);
                        log_P_t(no_trials,1:length(f),1:length(out_times))=log_P_timecourse(:,:)-log_P_timecourse_ref(:,:);
                    end
                    phase_per_trial(no_trials)=circ_mean(phase_this_trial');
                    
                    %for debugging get theta
                    this_lpt=zeros(1,1);
                    this_lpt=mean(log_P_t(no_trials,(f>=6)&(f<=12),floor(length(out_times)/2)+1),2);
                else
                    log_P_t(no_trials,1:length(f),1:length(out_times))=zeros(length(f),length(out_times));
                    phase_per_trial(no_trials)=0;
                    this_lpt=0;
                end
               
                if handles.displayData==1
                    fprintf(1, 'Trial No: %d, no of events: %d, logP theta=%d\n',no_trials,no_evs_this_trial,this_lpt);
                end

                if (no_evs_this_trial>handles.max_events_per_sec*(handles.time_end-handles.time_start-2*pad_time))
                    no_events=no_events-no_evs_this_trial;
                    events=events(1:end-no_evs_this_trial);
                    time_per_event=time_per_event(1:end-no_evs_this_trial);
                    phase=phase(1:end-no_evs_this_trial);
                    ERLFP=ERLFP(1:end-no_evs_this_trial,:);
                end
                
                if isfield(handles,'drgbchoices')
                    for evTypeNo=1:length(handles.drgbchoices.evTypeNos)
                        switch handles.evTypeNo
                            case 1
                                %tstart is the reference event
                                if handles.drgbchoices.evTypeNos(evTypeNo)==1
                                    %This is tstart
                                    if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(handles.drgbchoices.referenceEvent).times(evNo))>0
                                        which_event(evTypeNo,no_trials)=1;
                                    else
                                        which_event(evTypeNo,no_trials)=0;
                                    end
                                else
                                    %These are not tstart, and the time
                                    %should be compared at OdorOn
                                    %This is tstart
                                    if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(2).times(evNo))>0
                                        which_event(evTypeNo,no_trials)=1;
                                    else
                                        which_event(evTypeNo,no_trials)=0;
                                    end
                                end
                            otherwise
                                %OdorOn is the reference event
                                if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(handles.drgbchoices.referenceEvent).times(evNo))>0
                                    which_event(evTypeNo,no_trials)=1;
                                else
                                    which_event(evTypeNo,no_trials)=0;
                                end
                        end
                        
                        
                        
                    end
                end
                
                
                
                
                
            end
        end
    end
    
    %end %if eventstamps...
end %for evNo

if handles.displayData==1
    fprintf(1, '\n');
    
    %Avearge event-related filtered LFP
    try
        close 4
    catch
    end
    
    hFig4 = figure(4);
    set(hFig4, 'units','normalized','position',[.0 .5 .25 .25])
    
    no_time_pts=floor(handles.window*handles.drg.session(sessionNo).draq_p.ActualRate)+1;
    times=[1:no_time_pts]/handles.drg.session(sessionNo).draq_p.ActualRate;
    times=times-(handles.window/2);
    
    ERLFP_per_trial=ERLFP_per_trial-mean(mean(ERLFP_per_trial,1));
    meanERLFP_per_trial=mean(ERLFP_per_trial,1);
    CI = bootci(1000, @mean, ERLFP_per_trial);
    
    mean_ERLFP=[meanERLFP_per_trial;meanERLFP_per_trial];
    ERLFP_CI=CI-mean_ERLFP;
    ERLFP_CI(1,:)=-ERLFP_CI(1,:);
    
    [hl1, hp1] = boundedline(times',meanERLFP_per_trial', ERLFP_CI', 'b');
    
    
    title('Event-related LFP')
    xlim([-0.5 0.5])
    
    pct1=prctile(meanERLFP_per_trial,1);
    pct99=prctile(meanERLFP_per_trial,99);
    max_mean=pct99+0.1*(pct99-pct1);
    min_mean=pct1-0.1*(pct99-pct1);
    ylim([min_mean max_mean])
    
    ylabel('uV')
    xlabel('Time (s)')
    set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        
    
    
    %Event-related spectrogram
    %Timecourse doing average after log
 
    
    try
        close 1
    catch
    end
    
    %Plot the ERP timecourse
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.07 .1 .55 .3])
    
    %Calculate the mean ERP timecourse
    mean_log_P_t=zeros(length(f),length(out_times));
    mean_log_P_t(:,:)=mean(log_P_t,1);
    
     %Get max and min
    if handles.autoscale==1
        maxLogP=prctile(mean_log_P_t(:),99);
        minLogP=prctile(mean_log_P_t(:),1);
    else
        maxLogP=handles.maxLogP;
        minLogP=handles.minLogP;
    end
    
    
    %Note: Diego added this on purpose to limit the range to 10 dB
    %This results in emphasizing changes in the top 10 dB
    if maxLogP-minLogP>10
        minLogP=maxLogP-10;
    end
    
    
    drg_pcolor(repmat(out_times-mean(out_times),length(f),1)',repmat(f,1,length(out_times))',mean_log_P_t')
    
    
    colormap jet
    shading interp
    caxis([minLogP maxLogP]);
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)');
    title(['Event-related power wavelet spectrogram (dB) for ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
    
    try
        close 2
    catch
    end
    
    hFig2 = figure(2);
    set(hFig2, 'units','normalized','position',[.63 .1 .05 .3])
    
    prain=[minLogP:(maxLogP-minLogP)/99:maxLogP];
    drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    colormap jet
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')
    ylabel('dB')
    
    %Phase rose plot
    try
        close 5
    catch
    end
    
    hFig5 = figure(5);
    set(hFig5, 'units','normalized','position',[.69 .1 .3 .3])
    
    if no_events>0
        polarhistogram(pi*phase/180,6)
        title('LFP phase of events')
    end
    
    try
        close 3
    catch
    end
    
    hFig3 = figure(3);
    set(hFig3, 'units','normalized','position',[.25 .5 .25 .25])
    edges=[handles.time_start+handles.time_pad:0.1:handles.time_end-handles.time_pad];
    histogram(time_per_event,edges)
    title('Histogram for time of events')
    xlabel('Time (sec)')
    pffft=1;
end




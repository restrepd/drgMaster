function [log_P_t,no_trials_w_event,which_event,f,out_times,times,phase_per_trial,no_trials,no_events_per_trial,t_per_event_per_trial,trial_map,perCorrERP,no_ref_evs_per_trial]=drgERWAtimecourse(handles)
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
time_per_event_vetted=[];
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


phase_per_trial=[];
trial_map=[];

log_P_t=[];

inter_lick_intervals_ref=[];
ii_ili_ref=0;

inter_lick_intervals=[];
ii_ili=0;

ERLFP_per_event=[];

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
                            else
                                ii_ili_ref=ii_ili_ref+1;
                                inter_lick_intervals_ref(ii_ili_ref)=these_lick_times(this_lick_ii);
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
                time_per_these_events=[];
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
                        else
                            ii_ili=ii_ili+1;
                            inter_lick_intervals(ii_ili)=these_lick_times(this_lick_ii);
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
                                    time_per_these_events(no_evs_this_trial)=handles.time_start+pad_time+(ii/handles.drg.session(sessionNo).draq_p.ActualRate);
                                    
                                    t_per_event_per_trial(no_trials,no_evs_this_trial)=time_per_event(no_events);
                                    
                                    ERLFP_per_event(no_events,:)=LFP(1,floor(ii_start+ii-(handles.window/2)*handles.drg.session(sessionNo).draq_p.ActualRate):...
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
                        log_P_t_per_event(no_events-no_evs_this_trial+1:no_events,1:length(f),1:length(out_times))=10*log10(all_Power_these_events);
                    else
                        log_P_timecourse=zeros(length(f),length(out_times));
                        log_P_timecourse(:,:)=mean(10*log10(all_Power_these_events),1);
                        log_P_timecourse_ref=zeros(length(f),length(out_times));
                        log_P_timecourse_ref(:,:)=mean(10*log10(ref_Power_these_events),1);
                        log_P_t(no_trials,1:length(f),1:length(out_times))=log_P_timecourse(:,:)-log_P_timecourse_ref(:,:);
                        log_P_t_ref=zeros(1,length(f),length(out_times));
                        log_P_t_ref(1,:,:)=log_P_timecourse_ref;
                        log_P_t_per_event(no_events-no_evs_this_trial+1:no_events,1:length(f),1:length(out_times))=...
                            10*log10(all_Power_these_events)-repmat(log_P_t_ref,no_evs_this_trial,1,1);
                        
                    end
                    
                    phase_per_trial(no_trials)=circ_mean(phase_this_trial');
                    time_per_event_vetted=[time_per_event_vetted time_per_these_events];
                    
                    %for debugging get theta
                    this_lpt=zeros(1,1);
                    this_lpt=mean(log_P_t(no_trials,(f>=6)&(f<=12),floor(length(out_times)/2)+1),2);
                    
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
                    
                else
                    %Do not use this trial
                    no_trials=no_trials-1;
                    this_lpt=0;
                    trial_map=trial_map(1:no_trials);
                    no_events_per_trial=no_events_per_trial(1:no_trials);
                    no_ref_events_per_trial=no_ref_events_per_trial(1:no_trials);
                    t_per_event_per_trial=t_per_event_per_trial(1:no_trials);
                    time=time(1:no_trials);
                    which_trial=which_trial(1:no_trials);
                    perCorr_per_histo=perCorr_per_histo(1:no_trials);
                    perCorrERP=perCorrERP(1:no_trials);
                    phase_per_trial=phase_per_trial(1:no_trials);
                    log_P_t=log_P_t(1:no_trials,1:length(f),1:length(out_times));
                    no_events_per_trial=no_events_per_trial(1:no_trials);
                end
                
                if handles.displayData==1
                    fprintf(1, 'Trial No: %d, no of events: %d, logP theta=%d\n',no_trials,no_evs_this_trial,this_lpt);
                end
                
                if (no_evs_this_trial>=handles.max_events_per_sec*(handles.time_end-handles.time_start-2*pad_time))||...
                        (no_ref_evs_this_trial>=handles.max_events_per_sec*(handles.time_end-handles.time_start-2*pad_time))
                    no_events=no_events-no_evs_this_trial;
                    events=events(1:end-no_evs_this_trial);
                    time_per_event=time_per_event(1:end-no_evs_this_trial);
                    phase=phase(1:end-no_evs_this_trial);
                    ERLFP_per_event=ERLFP_per_event(1:end-no_evs_this_trial,:);
                end
                
                
                
                
                
                
                
            end
        end
    end
    
    %end %if eventstamps...
end %for evNo

if handles.displayData==1
    
    fprintf(1, '\n');
    
    %Note: close all closes drgMaster!
    for ii=1:12
        try
            close(ii)
        catch
        end
    end
    
    %Initialize variables needed to generate the figures
    t_bin=0.05;
    t_start=handles.time_start+handles.time_pad;
    t_end=handles.time_end-handles.time_pad;
    times=(t_start:t_bin:t_end-t_bin)+t_bin/2;
    
    %Timecourse for event-related filtered LFP
    try
        close 1
    catch
    end
    
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.05 .1 .65 .3])
    
    
    no_ERLFP_time_pts=floor(handles.window*handles.drg.session(sessionNo).draq_p.ActualRate)+1;
    
    %Calculate the ERLFP per time bin
    ERLFP_per_t_bin=zeros((t_end-t_start)/t_bin,no_ERLFP_time_pts);
    
    ii=0;
    for t=t_start:t_bin:t_end-t_bin
        ii=ii+1;
        t_mask=(time_per_event>=t)&(time_per_event<t+t_bin);
        if sum(t_mask)>=1
            ERLFP_per_t_bin(ii,1:no_ERLFP_time_pts)=mean(ERLFP_per_event(t_mask,:),1)';
        else
            if ii>1
                ERLFP_per_t_bin(ii,1:no_ERLFP_time_pts)=ERLFP_per_t_bin(ii-1,1:no_ERLFP_time_pts);
            end
        end
    end
    
    ERLFP_per_t_bin=ERLFP_per_t_bin-repmat(mean(ERLFP_per_t_bin,2),1,no_ERLFP_time_pts);
    
    delta_times=[1:no_ERLFP_time_pts]/handles.drg.session(sessionNo).draq_p.ActualRate;
    
    drg_pcolor(repmat(times',1,no_ERLFP_time_pts),repmat(delta_times,length(times),1),ERLFP_per_t_bin)
    
    minERLFP=prctile(ERLFP_per_t_bin(:),1);
    maxERLFP=prctile(ERLFP_per_t_bin(:),99);
    
    colormap jet
    shading interp
    caxis([minERLFP maxERLFP]);
    xlabel('Time (sec)')
    ylabel('Delta time (sec)');
    title(['Mean event-related LFP timecourse for ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
    set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
    
    %Figure 2 rainbow
    
    
    
    try
        close 4
    catch
    end
    
    hFig4 = figure(4);
    set(hFig4, 'units','normalized','position',[.072 .12 .55 .3])
    edges=[t_start:t_bin:t_end];
    h1=histogram(time_per_event,edges);
    licks_per_sec=h1.Values*(1/t_bin)*(1/no_trials);
    bar(edges(1:length(licks_per_sec)),licks_per_sec)
    title('Histogram for time of events')
    xlabel('Time (sec)')
    ylabel('Licks per sec')
    xlim([t_start t_end])
    %Event-related timecourse spectrograms
    
    %Calculate the log Power per time bin
    
    
    log_P_t__per_t_bin=zeros((t_end-t_start)/t_bin,length(f),length(out_times));
    
    
    ii=0;
    for t=t_start:t_bin:t_end-t_bin
        ii=ii+1;
        t_mask=(time_per_event_vetted>=t)&(time_per_event_vetted<t+t_bin);
        if sum(t_mask)>=1
            log_P_t__per_t_bin(ii,1:length(f),1:length(out_times))=mean(log_P_t_per_event(t_mask,1:length(f),1:length(out_times)),1);
        else
            if ii>1
                log_P_t__per_t_bin(ii,1:length(f),1:length(out_times))=log_P_t__per_t_bin(ii-1,1:length(f),1:length(out_times));
            end
        end
    end
    
        %Get max and min
    if handles.autoscale==1
        mindB=prctile(log_P_t__per_t_bin(:),1);
        maxdB=prctile(log_P_t__per_t_bin(:),99);
    else
        maxdB=handles.maxLogP;
        mindB=handles.minLogP;
    end
    
    
    %Plot timecourse for different time shifts
    figNo=4;
    
    
    
    for t_shift=-0.2:0.1:0.2
        figNo=figNo+1;
        try
            close figNo
        catch
        end
        
        hFig = figure(figNo);
        set(hFig, 'units','normalized','position',[.07+0.02*(figNo-2) .1+0.02*(figNo-2) .55 .3])
        
        [minabs,ii_t_shift]=min(abs(out_times-t_shift));
        
        this_log_P_t__per_t_bin=zeros((t_end-t_start)/t_bin,length(f));
        this_log_P_t__per_t_bin(:,:)=log_P_t__per_t_bin(:,:,ii_t_shift);
        
        drg_pcolor(repmat(times',1,length(f)),repmat(f',length(times),1),this_log_P_t__per_t_bin)
        
        
    
        colormap jet
        shading interp
        caxis([mindB maxdB]);
        xlabel('Time (sec)')
        ylabel('Frequency (Hz)');
        title(['Mean power event-related LFP timecourse (dB) for ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo} ' time shift= ' num2str(t_shift)])
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        
    end
    
    %Rainbow for dB
    try
        close 3
    catch
    end
    
    hFig3 = figure(3);
    set(hFig3, 'units','normalized','position',[.8 .1 .05 .3])
    
    
    
    prain=[mindB:(maxdB-mindB)/99:maxdB];
    drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    colormap jet
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')
    ylabel('dB')
    
    
    pffft=1;
    
    
end




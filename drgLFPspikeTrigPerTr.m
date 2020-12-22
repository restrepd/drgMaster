function drgLFPspikeTrigPerTr(handles)
%Generates a spike-triggered LFP in a trial range

sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
freq=handles.burstLowF:ceil((handles.burstHighF-handles.burstLowF)/60):handles.burstHighF;
highF1=handles.burstLowF;
highF2=handles.burstHighF;
pad_time=handles.time_pad;
n_phase_bins=handles.n_phase_bins;
window=round(handles.window*handles.drg.draq_p.ActualRate); 
noverlap=round(handles.noverlap*handles.drg.draq_p.ActualRate); 

%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;


%Butterworth filter for spike-triggered LFP
fpass=[handles.burstLowF handles.burstHighF];


bpFilt = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',fpass(1),'HalfPowerFrequency2',fpass(2), ...
    'SampleRate',handles.drg.session(sessionNo).draq_p.ActualRate);

%Now get the phase of gamma bursts witin theta
no_trials=0;
no_spikes=0;
no_spikes_ref=0;
no_spikes_perm=0;

spike_tr_LFP=[];
filt_spike_tr_LFP=[];

spike_tr_LFPperm=[];
filt_spike_tr_LFPperm=[];

window_dt=window/Fs;

valid_events=[];
no_valid=0;

for trNo=firstTr:lastTr
    
    if handles.save_drgb==0
        trial_no=trNo
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    
    if evNo~=-1
        
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            %Read the LFP spanning from the beggining to the end of the times
            %looked at
            if (handles.time_start-(window_dt/2))<(handles.startRef-(window_dt/2))
                this_start=handles.time_start-(window_dt/2);
            else
                this_start=handles.startRef-(window_dt/2);
            end
            
            if (handles.time_end+(window_dt/2))>(handles.endRef+(window_dt/2))
                this_end=handles.time_end+(window_dt/2);
            else
                this_end=handles.endRef+(window_dt/2);
            end
            
            [LFP, trialNo, can_read] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, this_start, this_end);
            ffLFP=filtfilt(bpFilt,LFP);
            delta_ii=ceil((window_dt/2)*handles.drg.session(sessionNo).draq_p.ActualRate);
            
            if (can_read==1)
                
                
                this_event_time=handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo);
                 
                %Get the spikes
                ii_spikes=find((handles.drg.unit(handles.unitNo).spike_times>=this_event_time+handles.time_start+handles.time_pad)&...
                    (handles.drg.unit(handles.unitNo).spike_times<=this_event_time+handles.time_end-handles.time_pad));
                
                ii_spikes_ref=find((handles.drg.unit(handles.unitNo).spike_times>=this_event_time+handles.startRef+handles.time_pad)&...
                    (handles.drg.unit(handles.unitNo).spike_times<=this_event_time+handles.endRef-handles.time_pad));
                
                
                no_trials=no_trials+1;
                
                if length(ii_spikes>0)&length(ii_spikes_ref>0)
                    no_valid=no_valid+1;
                    valid_events(no_valid)=evNo;
                    no_ii_spikes(no_valid)=length(ii_spikes);
                    no_ii_spikes_ref(no_valid)=length(ii_spikes_ref);
                    
                    for kk=1:length(ii_spikes)
                        no_spikes=no_spikes+1;
                        time_this_spike=handles.drg.unit(handles.unitNo).spike_times(ii_spikes(kk))-handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo);
                        ii_this_spike=ceil(1+(time_this_spike-this_start)*handles.drg.session(sessionNo).draq_p.ActualRate);
                        thisffLFP=[];
                        thisffLFP=ffLFP(ii_this_spike-delta_ii:ii_this_spike+delta_ii);
                        filt_spike_tr_LFP(no_spikes,1:length(thisffLFP))=thisffLFP;
                    end
                    
                    for kk=1:length(ii_spikes_ref)
                        no_spikes_ref=no_spikes_ref+1;
                        time_this_spike_ref=handles.drg.unit(handles.unitNo).spike_times(ii_spikes_ref(kk))-handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo);
                        ii_this_spike_ref=ceil(1+(time_this_spike_ref-this_start)*handles.drg.session(sessionNo).draq_p.ActualRate);
                        thisffLFP_ref=[];
                        thisffLFP_ref=ffLFP(ii_this_spike_ref-delta_ii:ii_this_spike_ref+delta_ii);
                        filt_spike_tr_LFP_ref(no_spikes_ref,1:length(thisffLFP_ref))=thisffLFP_ref;
                    end
                    
                end
                
            end
        end
    end
end



try
    close 1
catch
end

%Plot the timecourse
hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.02 .4 .5 .5])


min_points=min([size(filt_spike_tr_LFP_ref,2) size(filt_spike_tr_LFP,2)]);
time=[1:min_points]/Fs;
time=time-window_dt/2;


%Plot the filtered spike-triggered LFP
%subplot(2,1,1)
shadedErrorBar(time,mean(filt_spike_tr_LFP(:,1:min_points),1),std(filt_spike_tr_LFP(:,1:min_points),0,1)/sqrt(size(filt_spike_tr_LFP,1)),'-b')
hold on

mint=-3/fpass(1);
maxt=3/fpass(1);
first_ii=find(time>mint,1);
if isempty(first_ii)
    first_ii=1;
end
if time(first_ii)>mint
    mint=time(first_ii);
end
last_ii=find(time>maxt,1);
if isempty(last_ii)
    last_ii=length(time);
end
if time(last_ii)<maxt
    maxt=time(last_ii);
end

meanf=mean(filt_spike_tr_LFP,1)';
pct99=prctile(meanf(first_ii:last_ii),99);
pct1=prctile(meanf(first_ii:last_ii),1);
plot([0 0],[pct1-0.3*abs(pct99-pct1) pct99+0.3*abs(pct99-pct1)],'-r')
ylim([pct1-0.3*abs(pct99-pct1) pct99+0.3*abs(pct99-pct1)])
xlim([mint maxt])
title('Filtered spike-triggered LFP vs time')
xlabel('Time (sec)')


%Find the phase for the spike

angle_meanf = angle(hilbert(meanf)); % phase modulation of amplitude
[minf,ii_zero]=min(abs(time));
spike_phase=(360/(2*pi()))*(angle_meanf(ii_zero)+pi())
pffft=1
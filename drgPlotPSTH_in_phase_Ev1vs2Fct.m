function drgPlotPSTH_in_phase_Ev1vs2Fct(handles)
%PSTH plot for spikes in phase with the spike histogram calculated by drgSpikePhase

unitNo=handles.unitNo;
drg=handles.drg;
sessionNo=drg.unit(unitNo).sessionNo;

%Butterworth filter for spike-triggered LFP
fpass=[handles.burstLowF handles.burstHighF];


bpFilt = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',fpass(1),'HalfPowerFrequency2',fpass(2), ...
    'SampleRate',handles.drg.session(sessionNo).draq_p.ActualRate);
no_bins=length(handles.norm_phase_histo)-1;

try
    close 1
catch
end



evTypeNo=handles.evTypeNo;
evTypeNo2=handles.evTypeNo2;
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;

%Enter unit and event
%unitNo=2;



%Enter the event type
%   Events 1 through 6
%     'TStart'    'OdorOn'    'Hit'    'HitE'    'S+'    'S+E'
%   Events 7 through 13
%     'Miss'    'MissE'    'CR'    'CRE'    'S-'    'S-E'    'FA'
%   Events 14 through 17
%     'FAE'    'Reinf'    'L+'    'L-'
%evTypeNo=1;


%bin_size=0.10;
bin_size=0.05;
%bin_size=0.005


textout='drgPlotPSTH'


nobins=fix((handles.time_end-handles.time_start-2*handles.time_pad)/bin_size);

%Calculate the PSTH for event 2
PSTH2=zeros(1,nobins);
time=handles.time_start+handles.time_pad:(handles.time_end-handles.time_start-2*handles.time_pad)/nobins:handles.time_end-handles.time_pad;
time=time(1:end-1);

trial_num2=0;
for trNo=firstTr:lastTr
    
    evNo = drgFindEvNo(handles,trNo,sessionNo,evTypeNo2);
    if evNo~=-1
        excludeTrial=drgExcludeTrial(drg,drg.unit(unitNo).channel,drg.session(sessionNo).events(evTypeNo2).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            trial_num2=trial_num2+1;
        end
    end
    %end %if eventstamps...
end %for evNo

PSTH2_per_trial=zeros(trial_num2,nobins);

noTrials2=0;
spike_times=[];
spike_times=drg.unit(unitNo).spike_times;

for trNo=firstTr:lastTr
    
    if handles.save_drgb==0
        trialNo=trNo
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo,evTypeNo2);
    if evNo~=-1
        
        excludeTrialSp=drgExcludeTrial(drg,drg.unit(unitNo).channel,drg.session(sessionNo).events(evTypeNo2).times(evNo),sessionNo);
        excludeTrialLFP=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo2).times(evNo),sessionNo);

        
        if (excludeTrialSp==0)&(excludeTrialLFP==0)
            

            [LFP, trialNo, can_read] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, evTypeNo2, handles.time_start+handles.time_pad, handles.time_end-handles.time_pad);
        
            if (can_read==1)
                
                ffLFP=filtfilt(bpFilt,LFP);
                angle_ffLFP = angle(hilbert(ffLFP)); % phase
                
                noTrials2=noTrials2+1;
                these_spikes2=(spike_times>drg.session(sessionNo).events(evTypeNo2).times(evNo)+handles.time_start+handles.time_pad)&...
                    (spike_times<=drg.session(sessionNo).events(evTypeNo2).times(evNo)+handles.time_end-handles.time_pad);
                these_spike_times2=spike_times(these_spikes2)-(drg.session(sessionNo).events(evTypeNo2).times(evNo)+handles.time_start+handles.time_pad);
                
                for spk=1:length(these_spike_times2)
                    this_LFP_ii=ceil(these_spike_times2(spk)*handles.drg.session(sessionNo).draq_p.ActualRate);
                    this_phase=angle_ffLFP(this_LFP_ii);
                    this_phase_bin=ceil((this_phase+pi)/(2*pi/no_bins));
                    this_bin=ceil(these_spike_times2(spk)/bin_size);
                    PSTH2(1,this_bin)=PSTH2(1,this_bin)+handles.norm_phase_histo(this_phase_bin);
                    PSTH2_per_trial(noTrials2,this_bin)=PSTH2_per_trial(noTrials2,this_bin)+handles.norm_phase_histo(this_phase_bin);
                end %for spk
                
            end
        end
    end
    %end %if eventstamps...
end %for evNo

number_of_trials_included2=noTrials2

PSTH2=PSTH2/(noTrials2*bin_size);
PSTH2_per_trial=PSTH2_per_trial/bin_size;

%Now plot the PSTH for event 2
figure(1)
shadedErrorBar(time,mean(PSTH2_per_trial,1),std(PSTH2_per_trial,0,1)/sqrt(noTrials2),'-b')

%Calculate the PSTH for event 1
PSTH=zeros(1,nobins);


trial_num=0;
for trNo=firstTr:lastTr
    
    evNo = drgFindEvNo(handles,trNo,sessionNo,evTypeNo);
    
    if evNo~=-1
        excludeTrial=drgExcludeTrial(drg,drg.unit(unitNo).channel,drg.session(sessionNo).events(evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            trial_num=trial_num+1;
        end
    end
    %end %if eventstamps...
end %for evNo

PSTH_per_trial=zeros(trial_num,nobins);

noTrials=0;
spike_times=[];
spike_times=drg.unit(unitNo).spike_times;

for trNo=firstTr:lastTr
    
    if handles.save_drgb==0
        trialNo=trNo
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo,evTypeNo);
    if evNo~=-1
        excludeTrialSp=drgExcludeTrial(drg,drg.unit(unitNo).channel,drg.session(sessionNo).events(evTypeNo).times(evNo),sessionNo);
        excludeTrialLFP=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);

        if (excludeTrialSp==0)&(excludeTrialLFP==0)
            
            [LFP, trialNo, can_read] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start+handles.time_pad, handles.time_end-handles.time_pad);
            
            
            if (can_read==1)
                
                ffLFP=filtfilt(bpFilt,LFP);
                angle_ffLFP = angle(hilbert(ffLFP)); % phase
                
                noTrials=noTrials+1;
                these_spikes=(spike_times>drg.session(sessionNo).events(evTypeNo).times(evNo)+handles.time_start+handles.time_pad)&...
                    (spike_times<=drg.session(sessionNo).events(evTypeNo).times(evNo)+handles.time_end-handles.time_pad);
                these_spike_times=spike_times(these_spikes)-(drg.session(sessionNo).events(evTypeNo).times(evNo)+handles.time_start+handles.time_pad);
                
                for spk=1:length(these_spike_times)
                    this_LFP_ii=ceil(these_spike_times(spk)*handles.drg.session(sessionNo).draq_p.ActualRate);
                    this_phase=angle_ffLFP(this_LFP_ii);
                    this_phase_bin=ceil((this_phase+pi)/(2*pi/no_bins));
                    this_bin=ceil(these_spike_times(spk)/bin_size);
                    PSTH(1,this_bin)=PSTH(1,this_bin)+handles.norm_phase_histo(this_phase_bin);
                    PSTH_per_trial(noTrials,this_bin)=PSTH_per_trial(noTrials,this_bin)+handles.norm_phase_histo(this_phase_bin);
                end %for spk
            end
        end
    end
    %end %if eventstamps...
end %for evNo

number_of_trials_included=noTrials

PSTH=PSTH/(noTrials*bin_size);
PSTH_per_trial=PSTH_per_trial/bin_size;

%Now plot the PSTH

hold on
shadedErrorBar(time,mean(PSTH_per_trial,1),std(PSTH_per_trial,0,1)/sqrt(noTrials),'-r')
title(['PSTH for ' drg.session.eventlabels{evTypeNo} ' (red) vs. ' drg.session.eventlabels{evTypeNo2} ' (blue)'])
ylabel('Frequency (Hz)')
xlabel('Time (sec)')



pfft=1;
%save('/Users/restrepd/Documents/Grants/Complex odor 2017/Figures/Ming/1017ch1u1splusPSTH.mat','PSTH_per_trial','noTrials','time')





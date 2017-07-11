function handles=drgPlotPSTHFct(handles)
%PSTH plot for a range of trials



drg=handles.drg;
unitNo=handles.unitNo;
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;

%Enter unit and event
%unitNo=2;

sessionNo=drg.unit(unitNo).sessionNo;

%Enter the event type
%   Events 1 through 6
%     'TStart'    'OdorOn'    'Hit'    'HitE'    'S+'    'S+E'
%   Events 7 through 13
%     'Miss'    'MissE'    'CR'    'CRE'    'S-'    'S-E'    'FA' 
%   Events 14 through 17
%     'FAE'    'Reinf'    'L+'    'L-'
%evTypeNo=1;


bin_size=0.05;

nobins=ceil((handles.time_end-handles.time_start-2*handles.time_pad)/bin_size); 

PSTH=zeros(1,nobins);
time=handles.time_start+handles.time_pad:(handles.time_end-handles.time_start-2*handles.time_pad)/nobins:handles.time_end-handles.time_pad;
time=time(1:end-1);

[perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);

trial_num=0;
for trNo=firstTr:lastTr
    
    evNo = drgFindEvNo(handles,trNo,sessionNo,handles.evTypeNo);

    if evNo~=-1
        excludeTrial=drgExcludeTrial(drg,drg.unit(unitNo).channel,drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
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
    
     evNo = drgFindEvNo(handles,trNo,sessionNo,handles.evTypeNo);
    
     if evNo~=-1
        excludeTrial=drgExcludeTrial(drg,drg.unit(unitNo).channel,drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            noTrials=noTrials+1;
            these_spikes=(spike_times>drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+handles.time_start+handles.time_pad)&...
                (spike_times<=drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+handles.time_end-handles.time_pad);
            these_spike_times=spike_times(these_spikes)-(drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+handles.time_start+handles.time_pad);
            
            for spk=1:length(these_spike_times)
                this_bin=ceil(these_spike_times(spk)/bin_size);
                PSTH(1,this_bin)=PSTH(1,this_bin)+1;
                PSTH_per_trial(noTrials,this_bin)=PSTH_per_trial(noTrials,this_bin)+1;
            end %for spk
            
            if handles.save_drgb==1
                for eventTypeNo=1:length(handles.drgbchoices.evTypeNos)
                    if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(eventTypeNo)).times==handles.drg.session(1).events(handles.drgbchoices.referenceEvent).times(evNo))>0
                        handles.drgb.unit(handles.drgb.unit_no).which_event(eventTypeNo,noTrials)=1;
                    else
                        handles.drgb.unit(handles.drgb.unit_no).which_event(eventTypeNo,noTrials)=0;
                    end
                end
                
                
                handles.drgb.unit(handles.drgb.unit_no).perCorr(noTrials)=perCorr(find(abs(handles.drg.session(1).events(handles.evTypeNo).times(evNo)-handles.drg.session(1).events(2).times)...
                    <= handles.max_dt_between_events,1,'first'));
                
            end
        end
     end
    %end %if eventstamps...
end %for evNo

number_of_trials_included=noTrials;

PSTH=PSTH/(noTrials*bin_size);
PSTH_per_trial=PSTH_per_trial/bin_size;


%Now plot the PSTH
if handles.displayData==1
    try
        close 1
    catch
    end
    figure(1)
    %bar(time,PSTH,'b');
    shadedErrorBar(time,mean(PSTH_per_trial,1),std(PSTH_per_trial,0,1)/sqrt(noTrials),'-b')
    title(['PSTH for ' drg.session.eventlabels{handles.evTypeNo}])
    ylabel('Frequency (Hz)')
    xlabel('Time (sec)')
end


%Save the data if requested by drgb
if handles.save_drgb==1
    handles.drgb.unit(handles.drgb.unit_no).PSTH_per_trial=PSTH_per_trial;
    handles.drgb.time=time;    
end



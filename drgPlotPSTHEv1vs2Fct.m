function drgPlotPSTHEv1vs2Fct(handles)
%PSTH plot for a range of trials

try
    close 1
catch
end

drg=handles.drg;
unitNo=handles.unitNo;
evTypeNo=handles.evTypeNo;
evTypeNo2=handles.evTypeNo2;
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
    
    evNo = drgFindEvNo(handles,trNo,sessionNo,evTypeNo2);
    if evNo~=-1
        excludeTrial=drgExcludeTrial(drg,drg.unit(unitNo).channel,drg.session(sessionNo).events(evTypeNo2).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            noTrials2=noTrials2+1;
            these_spikes2=(spike_times>drg.session(sessionNo).events(evTypeNo2).times(evNo)+handles.time_start+handles.time_pad)&...
                (spike_times<=drg.session(sessionNo).events(evTypeNo2).times(evNo)+handles.time_end-handles.time_pad);
            these_spike_times2=spike_times(these_spikes2)-(drg.session(sessionNo).events(evTypeNo2).times(evNo)+handles.time_start+handles.time_pad);
            
            for spk=1:length(these_spike_times2)
                this_bin=ceil(these_spike_times2(spk)/bin_size);
                PSTH2(1,this_bin)=PSTH2(1,this_bin)+1;
                PSTH2_per_trial(noTrials2,this_bin)=PSTH2_per_trial(noTrials2,this_bin)+1;
            end %for spk
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
    
    evNo = drgFindEvNo(handles,trNo,sessionNo,evTypeNo);
    if evNo~=-1
        excludeTrial=drgExcludeTrial(drg,drg.unit(unitNo).channel,drg.session(sessionNo).events(evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            noTrials=noTrials+1;
            these_spikes=(spike_times>drg.session(sessionNo).events(evTypeNo).times(evNo)+handles.time_start+handles.time_pad)&...
                (spike_times<=drg.session(sessionNo).events(evTypeNo).times(evNo)+handles.time_end-handles.time_pad);
            these_spike_times=spike_times(these_spikes)-(drg.session(sessionNo).events(evTypeNo).times(evNo)+handles.time_start+handles.time_pad);
            
            for spk=1:length(these_spike_times)
                this_bin=ceil(these_spike_times(spk)/bin_size);
                PSTH(1,this_bin)=PSTH(1,this_bin)+1;
                PSTH_per_trial(noTrials,this_bin)=PSTH_per_trial(noTrials,this_bin)+1;
            end %for spk
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





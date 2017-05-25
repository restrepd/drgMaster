function drgPlotPSTHFct(handles)
%PSTH plot for a range of trials

try
    close 1
catch
end

drg=handles.drg;
unitNo=handles.unitNo;
evTypeNo=handles.evTypeNo;
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

PSTH=zeros(1,nobins);
time=handles.time_start+handles.time_pad:(handles.time_end-handles.time_start-2*handles.time_pad)/nobins:handles.time_end-handles.time_pad;
time=time(1:end-1);

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

number_of_trials_included=noTrials;

PSTH=PSTH/(noTrials*bin_size);
PSTH_per_trial=PSTH_per_trial/bin_size;

%Now plot the PSTH
figure(1)
%bar(time,PSTH,'b');
shadedErrorBar(time,mean(PSTH_per_trial,1),std(PSTH_per_trial,0,1)/sqrt(noTrials),'-b')
title(['PSTH for ' drg.session.eventlabels{evTypeNo}])
ylabel('Frequency (Hz)')
xlabel('Time (sec)')


%save('/Users/restrepd/Documents/Grants/Complex odor 2017/Figures/Ming/1017ch1u1splusPSTH.mat','PSTH_per_trial','noTrials','time')





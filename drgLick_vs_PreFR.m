function [rholickFR,pvallickFR,percent_var_explained]=drgLick_vs_PreFR(handles)

%
%   Finds the relationship between FR and percent lick interval and event chosen by the user
%


%Fist get the licks
sessionNo=handles.sessionNo;
lfpElectrode=19; %19 is the lick

%Enter the event type
%   Events 1 through 6
%     'TStart'    'OdorOn'    'Hit'    'HitE'    'S+'    'S+E'
%   Events 7 through 13
%     'Miss'    'MissE'    'CR'    'CRE'    'S-'    'S-E'    'FA'
%   Events 14 through 19
%     'FAE'    'Reinf'    'L+'    'L-' 'S+TStart' 'S-TStart'
%   'S+TStart' = 18

%Enter trials
firstTr=1;
lastTr=handles.drg.session(sessionNo).events(handles.evTypeNo).noTimes;


allnoEvs1=0;
licks=[];

for evNo=firstTr:lastTr
    
    %Note: the same trials are excluded as those excluded for this unit in
    %the PSTH
    excludeTrial=drgExcludeTrial(handles.drg,handles.drg.unit(handles.unitNo).channel,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
    
    if excludeTrial==0
        
        thisLFP=[];
        [thisLFP, trialNo, can_read] = drgGetTrialLFPData(handles, lfpElectrode, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
        allnoEvs1=allnoEvs1+1;
        if (can_read==1)
            licks(1:length(thisLFP),allnoEvs1)=thisLFP;
        else
%             szLFP=size(licks);
%             licks(szLFP(1),allnoEvs1)=zeros(szLFP(1),1);
        end
    end
end

%Now calculate percent lick
szlicks=size(licks);



skip_artifact_n=ceil(handles.time_pad*handles.drg.session(sessionNo).draq_p.ActualRate); %Need to skip a large jump at the start due to the filtering
times=[1:(szlicks(1)-2*skip_artifact_n)+1]/handles.drg.session(sessionNo).draq_p.ActualRate;
times=times+handles.time_start+handles.time_pad;
thresholded_licks=zeros(szlicks(2),length(licks(skip_artifact_n:end-skip_artifact_n,1))');

for noTr=1:szlicks(2)
    threshold=((prctile(licks(:,noTr),99.5)-prctile(licks(:,noTr),0.5))/2)+prctile(licks(:,noTr),0.5);
    thresholded_licks(noTr,:)=(licks(skip_artifact_n:end-skip_artifact_n,noTr)>=threshold)';
end

thresholded_licks=double(thresholded_licks);
trials=1:szlicks(2);
try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.05 .15 .85 .3])
drg_pcolor(repmat(times,szlicks(2),1),repmat(trials',1,length(times)),double(thresholded_licks))
colormap jet
shading flat
xlabel('Time (sec)')
ylabel('Trial No');
title(['Licks per trial ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
caxis([0 1]);

%Now get the preFR


%Enter the event type
%   Events 1 through 6
%     'TStart'    'OdorOn'    'Hit'    'HitE'    'S+'    'S+E'
%   Events 7 through 13
%     'Miss'    'MissE'    'CR'    'CRE'    'S-'    'S-E'    'FA' 
%   Events 14 through 17
%     'FAE'    'Reinf'    'L+'    'L-'
%evTypeNo=1;


%bin_size=0.10;
bin_size=0.02;
%bin_size=0.005


textout='drgPlotPSTH'

time_pre=handles.time_start+handles.time_pad;
time_post=handles.time_end-handles.time_pad;

nobins=fix((time_post-time_pre)/bin_size); 

PSTH=zeros(1,nobins);
itime=1:nobins;
itime=itime+fix(time_pre/bin_size);
time=double(itime)*bin_size;

noTrials=0;
spike_times=[];
spike_times=handles.drg.unit(handles.unitNo).spike_times;
no_spikes_per_trial=[];

for evNo=firstTr:lastTr
    
    %evNo
    excludeTrial=drgExcludeTrial(handles.drg,handles.drg.unit(handles.unitNo).channel,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
    
    
    if excludeTrial==0
        
        noTrials=noTrials+1;
        these_spikes=(spike_times>handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+time_pre)&...
            (spike_times<=handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+time_post);
        these_spike_times=spike_times(these_spikes)-(handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+time_pre);
        no_spikes_per_trial(noTrials)=length(these_spike_times);
        for spk=1:length(these_spike_times)
            this_bin=ceil(these_spike_times(spk)/bin_size);
            PSTH(1,this_bin)=PSTH(1,this_bin)+1;
        end %for spk
    end
        %end
    %end %if eventstamps...
end %for evNo

number_of_trials_included=noTrials;



PSTH=PSTH/(noTrials*bin_size);

%Now plot the PSTH
try
    close 2
catch
end

hFig2 = figure(2);
set(hFig2, 'units','normalized','position',[.05 .55 .4 .3])
bar(time,PSTH,'b');
title(['PSTH for ' handles.drg.session.eventlabels{handles.evTypeNo}])
ylabel('Frequency (Hz)')
xlabel('Time (sec)')

%Now plot the percent lick vs. FR

try
    close 3
catch
end

hFig3 = figure(3);
set(hFig3, 'units','normalized','position',[.55 .55 .4 .3])

FR=(no_spikes_per_trial/(time_post-time_pre))';
percent_lick=100*sum(thresholded_licks,2)/length(times);
minPer=min(percent_lick);
maxPer=max(percent_lick);

beta_FR_per=nlinfit(percent_lick,FR,@dr_fitline,[0 1])

plot(percent_lick,FR,'ob')
hold on
per=[minPer-0.05*(maxPer-minPer) maxPer+0.05*(maxPer-minPer)];
plot(per, per*beta_FR_per(1)+beta_FR_per(2),'-b')
ylabel('Firing rate (Hz)')
xlabel('Percent lick')
xlim(per)
title('Firing rate vs. fpercent lick')

[rholickFR,pvallickFR] = corr(percent_lick,FR)

FR_var=var(FR);
remaining_var=sum((FR-(percent_lick*beta_FR_per(1)+beta_FR_per(2))).^2)/(length(FR)-1);
percent_var_explained=100*(FR_var-remaining_var)/FR_var

pfffft=1

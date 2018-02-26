function handles=drgLickPerTrialTimecourse(handles)

%
%   Finds the percent lick in the 0.5 to 2.5 sec interval for all OdorOn events
%


%Enter LFP tetrode and event
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

%     excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
%     
%     if excludeTrial==0
        
        thisLFP=[];
        [thisLFP, trialNo, can_read] = drgGetTrialLFPData(handles, lfpElectrode, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
        allnoEvs1=allnoEvs1+1;
        if (can_read==1)
            licks(1:length(thisLFP),allnoEvs1)=thisLFP;
        else
%             szLFP=size(licks);
%             licks(szLFP(1),allnoEvs1)=zeros(szLFP(1),1);
        end
%     end
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

pfffft=1

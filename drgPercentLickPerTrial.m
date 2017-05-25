function handles=drgPercentLickPerTrial(handles)

%
%   Finds the percent lick in the 0.5 to 2.5 sec interval for all OdorOn events
%


%Enter LFP tetrode and event
sessionNo=handles.sessionNo;
lfpElectrode=19; %19 is the lick
time_start=0.5;
time_end=2.5;

%Enter the event type
%   Events 1 through 6
%     'TStart'    'OdorOn'    'Hit'    'HitE'    'S+'    'S+E'
%   Events 7 through 13
%     'Miss'    'MissE'    'CR'    'CRE'    'S-'    'S-E'    'FA'
%   Events 14 through 19
%     'FAE'    'Reinf'    'L+'    'L-' 'S+TStart' 'S-TStart'
%   'S+TStart' = 18
odorOn=2; %Defaults to OdorOn, this will not work in some experiments where OdorOn is not defined as 2 
this_ev_t_no=handles.evTypeNo;
handles.evTypeNo=odorOn;

%Enter trials
firstTr=1;
lastTr=handles.drg.session(sessionNo).noTrials;


allnoEvs1=0;
licks=[];

for evNo=firstTr:lastTr
    
    %     excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
    %
    %     if excludeTrial==0
    
    thisLFP=[];
    [thisLFP, trialNo, can_read] = drgGetTrialLFPData(handles, lfpElectrode, evNo, odorOn, time_start, time_end);
    allnoEvs1=allnoEvs1+1;
    if (can_read==1)
        licks(1:length(thisLFP),allnoEvs1)=thisLFP;
    end
    %     end
end

%Now calculate percent lick
szlicks=size(licks);
threshold=((prctile(licks(:),99.5)-prctile(licks(:),0.5))/2)+prctile(licks(:),0.5);
for evNo=firstTr:lastTr
    handles.drg.session(sessionNo).percent_lick(evNo)=100*sum(licks(:,evNo)>threshold)/szlicks(1);
end
handles.evTypeNo=this_ev_t_no;
%pfffft=1

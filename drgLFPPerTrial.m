function drgLFPPerTrial(handles)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   plots a number of LFP trials for a certain event type
%


%Enter LFP tetrode and event
sessionNo=handles.sessionNo;
lfpElectrode=handles.peakLFPNo;

%Enter the event type
%   Events 1 through 6
%     'TStart'    'OdorOn'    'Hit'    'HitE'    'S+'    'S+E'
%   Events 7 through 13
%     'Miss'    'MissE'    'CR'    'CRE'    'S-'    'S-E'    'FA'
%   Events 14 through 19
%     'FAE'    'Reinf'    'L+'    'L-' 'S+TStart' 'S-TStart'
%   'S+TStart' = 18
evTypeNo=handles.evTypeNo;

%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;

fpass=[handles.burstLowF handles.burstHighF];


bpFilt = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',fpass(1),'HalfPowerFrequency2',fpass(2), ...
    'SampleRate',handles.drg.session(sessionNo).draq_p.ActualRate);


allnoEvs1=0;
LFP=[];

for trNo=firstTr:lastTr
    if handles.save_drgb==0
        trial_no=trNo
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);

    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            thisLFP=[];
            [thisLFP, trialNo, can_read] = drgGetTrialLFPData(handles, lfpElectrode, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            
            if (can_read==1)
                allnoEvs1=allnoEvs1+1;
                LFP(1:length(thisLFP),allnoEvs1)=thisLFP;
            end
        end
    end
end



szdfc=size(LFP);

skip_artifact_n=ceil(handles.time_pad*handles.drg.session(sessionNo).draq_p.ActualRate); %Need to skip a large jump at the start due to the filtering
times=[1:(szdfc(1)-2*skip_artifact_n)+1]/handles.drg.session(sessionNo).draq_p.ActualRate;
times=times+handles.time_start+handles.time_pad;

ffLFP=[];
this_ff=[];

if lfpElectrode<18
    %This is an LFP
    for ii=1:allnoEvs1
        this_ff=filtfilt(bpFilt,LFP(:,ii));
        
        ffLFP(:,ii)=this_ff;
    end
else
    %This is >=trigger
    ffLFP=LFP;
end

sdff=std(ffLFP(1:end));
deltaOne=prctile(ffLFP(1:end),99)-prctile(ffLFP(1:end),1);


yfac=6;


    
try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.05 .25 .43 .65])

for ii=1:allnoEvs1
%for ii=1:10
    plot(times, ffLFP(skip_artifact_n:end-skip_artifact_n,ii)+4*deltaOne*ii,'-b');
    hold on
end



pffft=1




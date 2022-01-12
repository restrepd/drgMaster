function handles=drgLickTraces(handles)

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
firstEvNo=handles.firstEvNo;
lastEvNo=handles.lastEvNo;

Fr=handles.drg.session(sessionNo).draq_p.ActualRate;

%Get traces for Ev1
allnoEvs1=0;
ev1_licks=[];

for evNo=firstEvNo:lastEvNo

%     excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
%     
%     if excludeTrial==0
        
        thisLFP=[];
        [thisLFP, trialNo, can_read] = drgGetTrialLFPData(handles, lfpElectrode, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
        allnoEvs1=allnoEvs1+1;
        if (can_read==1)
            ev1_licks(1:length(thisLFP),allnoEvs1)=thisLFP;
        else
%             szLFP=size(licks);
%             licks(szLFP(1),allnoEvs1)=zeros(szLFP(1),1);
        end
%     end
end

%Get traces for Ev2
allnoEvs2=0;
ev2_licks=[];

for evNo=firstEvNo:lastEvNo

%     excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
%     
%     if excludeTrial==0
        
        thisLFP=[];
        [thisLFP, trialNo, can_read] = drgGetTrialLFPData(handles, lfpElectrode, evNo, handles.evTypeNo2, handles.time_start, handles.time_end);
        allnoEvs2=allnoEvs2+1;
        if (can_read==1)
            ev2_licks(1:length(thisLFP),allnoEvs2)=thisLFP;
        else
%             szLFP=size(licks);
%             licks(szLFP(1),allnoEvs1)=zeros(szLFP(1),1);
        end
%     end
end

%Plot the lick traces
try
    close 1
catch
end

hFig1 = figure(1);
ax=gca;ax.LineWidth=3;

hold on

skip_artifact_n=ceil(handles.time_pad*handles.drg.session(sessionNo).draq_p.ActualRate); %Need to skip a large jump at the start due to the filtering

times=[1:(size(ev2_licks,1)-2*skip_artifact_n)+1]/handles.drg.session(sessionNo).draq_p.ActualRate;

hipctile=prctile(ev1_licks(:),99);
lowpctile=prctile(ev1_licks(:),1);
delta_trial=2*(hipctile-lowpctile);
y_offset=delta_trial;


%S- first
for trNo=1:size(ev2_licks,2)
   plot(times,ev2_licks(handles.time_pad*Fr:handles.time_pad*Fr+length(times)-1,trNo)+y_offset,'-','Color',[0 114/255 178/255], 'LineWidth',1.5) 
   y_offset=y_offset+delta_trial;
end

y_offset=y_offset+delta_trial*3;

%S+ next
for trNo=1:size(ev1_licks,2)
   plot(times,ev1_licks(handles.time_pad*Fr:handles.time_pad*Fr+length(times)-1, trNo)+y_offset,'-','Color',[158/255 31/255 99/255],'LineWidth',1.5) 
   y_offset=y_offset+delta_trial;
end


xlabel('Time (sec)')
title(['Lick traces'])

%Calculate the binary licks
lick_threshold=prctile(ev1_licks(:),1)+((prctile(ev1_licks(:),99)-prctile(ev1_licks(:),1))/2);
  
dt=handles.dt_tPRP;
t_pac=[handles.time_start+handles.time_pad:dt:handles.time_end-handles.time_pad];



lick_trace_times=times-2;

%Genetrate the licks timecourse for ev1
binary_lick_per_t_ev1=zeros(size(ev1_licks,2),length(t_pac));
for trNo=1:size(ev1_licks,2)
    for ii_t=1:length(t_pac)
        if sum(ev1_licks((lick_trace_times>(t_pac(ii_t)-((t_pac(2)-t_pac(1))/2)))&(lick_trace_times<=(t_pac(ii_t)+((t_pac(2)-t_pac(1))/2))),trNo)>lick_threshold)>0
            binary_lick_per_t_ev1(trNo,ii_t)=1;
        end
    end
end

binary_lick_per_t_ev2=zeros(size(ev2_licks,2),length(t_pac));
for trNo=1:size(ev2_licks,2)
    for ii_t=1:length(t_pac)
        if sum(ev2_licks((lick_trace_times>(t_pac(ii_t)-((t_pac(2)-t_pac(1))/2)))&(lick_trace_times<=(t_pac(ii_t)+((t_pac(2)-t_pac(1))/2))),trNo)>lick_threshold)>0
            binary_lick_per_t_ev2(trNo,ii_t)=1;
        end
    end
end

%Calculate p values
p_vals_licks=zeros(1,length(t_pac));
for ii_t=1:length(t_pac)
    p_vals_licks(ii_t)=ranksum(binary_lick_per_t_ev1(:,ii_t),binary_lick_per_t_ev2(:,ii_t));
end

try
    close(2)
catch
end
hFig=figure(2);

% set(hFig, 'units','normalized','position',[.1 .1 .7 .7])
%             subplot(2,1,1)
ax=gca;ax.LineWidth=3;

hold on
plot(t_pac,log10(p_vals_licks),'-k','LineWidth',3)
plot([t_pac(1) t_pac(end)],[log10(0.05) log10(0.05)],'-r','LineWidth', 2)

title('log10(p value) for licks')
ylabel('log10(p value)')
xlabel('Time (sec)')
pfft=1;


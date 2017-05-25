function handles=drgITI(handles)

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
odorOn=2; %Defaults to OdorOn, this will not work in some experiments where OdorOn is not defined as 2 
times=handles.drg.session(sessionNo).events(odorOn).times;
t0=times(1:end-1);
ITI=times(2:end)-times(1:end-1);

try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.05 .15 .5 .3])
plot(t0,ITI,'ob')
xlabel('Time (sec)')
ylabel('ITI (sec)');
title('Inter trial interval (sec)')


pfffft=1

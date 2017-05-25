function drgPlotPSTHBlockFct(drg, unitNo, evTypeNo)
%PSTH plot for a range of trials



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
bin_size=0.02;
%bin_size=0.005

%Enter trials



textout='drgPlotPSTH'


nobins=fix((drg.time_post-drg.time_pre)/bin_size);

PSTH=zeros(1,nobins);
itime=1:nobins;
itime=itime+fix(drg.time_pre/bin_size);
time=double(itime)*bin_size;


spike_times=[];
spike_times=drg.unit(unitNo).spike_times;

no_blocks=ceil(drg.session(sessionNo).events(2).noTimes/20);
PSTH=zeros(no_blocks,length(time));
 
for block=1:no_blocks
    
    noTrials(block)=0;
    firstTr=find(drg.session(sessionNo).events(evTypeNo).times>=drg.session(sessionNo).blocks(block,1),1,'first');
    lastTr=find(drg.session(sessionNo).events(evTypeNo).times<=drg.session(sessionNo).blocks(block,2),1,'last');
    
    for evNo=firstTr:lastTr
        
        %evNo
        excludeTrial=drgExcludeTrial(drg,drg.unit(unitNo).channel,drg.session(sessionNo).events(evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            noTrials(block)=noTrials(block)+1;
            these_spikes=(spike_times>drg.session(sessionNo).events(evTypeNo).times(evNo)+drg.time_pre)&...
                (spike_times<=drg.session(sessionNo).events(evTypeNo).times(evNo)+drg.time_post);
            these_spike_times=spike_times(these_spikes)-(drg.session(sessionNo).events(evTypeNo).times(evNo)+drg.time_pre);
            
            for spk=1:length(these_spike_times)
                this_bin=ceil(these_spike_times(spk)/bin_size);
                PSTH(block,this_bin)=PSTH(block,this_bin)+1;
            end %for spk
        end
        %end
        %end %if eventstamps...
    end %for evNo
    
    PSTH(block,:)=PSTH(block,:)/(noTrials(block)*bin_size);
 
end

maxPSTH=max(max(PSTH));

try
    close 1
catch
end

figure(1)

for block=1:no_blocks
    
    %Now plot the PSTH
    subplot(no_blocks,1,block)
    bar(time,PSTH(block,:),'b');
    ylim([0 1.2*maxPSTH]);
    ylabel(num2str(block))
    
    if block==1
        title(['PSTH for ' drg.session.eventlabels{evTypeNo}])
    end
    
    if block==no_blocks
        xlabel('Time (sec)')
    end
    
end

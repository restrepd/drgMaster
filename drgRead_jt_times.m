function drgRead_jt_times(drg_directory,jt_times_file)

drg.drg_directory=drg_directory;

drg.jt_times_file=jt_times_file;

drg.drg_output_file=[jt_times_file(10:end-4) '_drg.mat'];

spikes_processed=0;

%These are variables that are initialized or hard coded, do not change unless you really want to do something
%different

no_singles=0;
no_singles_Fee=0;
drg.noUnits=0;

%Pre-sets
drg.time_pre=-2.5;
drg.time_post=3;        %Note that this used to be 5 seconds and I changed it because of acquisition problems in Anan's data
drg.block_choice=1;     %0:nsampler, 1:spm; 2:cm; 3:osampler or spmult
drg.bin=0.1;            %Time bin for PSTH per trial

drg.excludeBadTraces=1;

drg.maxLFP=3900;
drg.minLFP=100;
drg.delta_max_min_out=0.6;

drg.bin_size=0.10;
drg.nobins=fix((drg.time_post-drg.time_pre)/drg.bin_size)+1;
drg.start_pre=1;
drg.end_pre=24;

% fid=fopen([drg.drg_directory drg.drg_file '.txt']);
% %fid=fopen([drg.drg_directory 'drg_files.txt']);
% %fid=fopen([drg.drg_directory 'drs_files_gogo_sniff.txt']);
%
% num_files=textscan(fid,'%d',1);
% files=textscan(fid,'%s',num_files{1});
% drg.fls=files{1};
%
% for (filNum=1:num_files{1})
%     file_number=filNum
%Load Wilder's data
filNum=1;
load([drg_directory jt_times_file])


drg.noSessions=filNum;

drg.draq_d=draq_d;
drg.draq_p=draq_p;
drg.drta_p=drta_p;

%Change in drta_p the name and directory to the current 
%jt_times name and directory
%Note: this is here to fix a bug that happens
%when the name of the file and jt_times are changed after
%running drta

drg.drta_p.PathName=[drg.drg_directory];
%Find if the file is .rhd or .dg
if strcmp(drg.drta_p.fullName(end-3:end-3),'.')
    %This is .rhd
    drg.drta_p.FileName=[drg.jt_times_file(10:end-4) '.rhd'];
    drg.drta_p.fullName=[drg.drg_directory drg.jt_times_file(10:end-4) '.rhd'];
else
    %This is .dg
    drg.drta_p.FileName=[drg.jt_times_file(10:end-4) '.dg'];
    drg.drta_p.fullName=[drg.drg_directory drg.jt_times_file(10:end-4) '.dg'];
end


%Enter events in drg
drg.nEventTypes=draq_d.nEventTypes;
nEventTypes=draq_d.nEventTypes;
drg.session(filNum).eventlabels=draq_d.eventlabels;
noEvs=zeros(1,nEventTypes);

for evNo=1:draq_d.noEvents
    noEvs(draq_d.eventType(evNo))=noEvs(draq_d.eventType(evNo))+1;
    drg.session(filNum).events(draq_d.eventType(evNo)).times(noEvs(draq_d.eventType(evNo)))=draq_d.events(evNo);
    drg.session(filNum).events(draq_d.eventType(evNo)).noTimes=noEvs(draq_d.eventType(evNo));
end

%Save trial starts
drg.session(filNum).trial_start=draq_d.t_trial;

%Now add events
%Adding events should be moved to drta and a patch program!


%If there are Hit, CR, etc then add lick+ and lick-: L+ and L-
if (~isempty(find(strcmp('Hit',drg.session(filNum).eventlabels), 1)))||...
        (~isempty(find(strcmp('FA',drg.session(filNum).eventlabels), 1)))||...
        (~isempty(find(strcmp('CR',drg.session(filNum).eventlabels), 1)))||...
        (~isempty(find(strcmp('Miss',drg.session(filNum).eventlabels), 1)))
    %Add L+
    drg.session(filNum).eventlabels{length(drg.session(filNum).eventlabels)+1}='L+';
    drg.session(filNum).events(length(drg.session(filNum).eventlabels)).times=sort([drg.session(filNum).events(strcmp(drg.session(filNum).eventlabels,'Hit')).times ...
        drg.session(filNum).events(strcmp(drg.session(filNum).eventlabels,'FA')).times]);
    drg.session(filNum).events(length(drg.session(filNum).eventlabels)).noTimes=length(sort([drg.session(filNum).events(strcmp(drg.session(filNum).eventlabels,'Hit')).times ...
        drg.session(filNum).events(strcmp(drg.session(filNum).eventlabels,'FA')).times]));
    
    %Add L-
    drg.session(filNum).eventlabels{length(drg.session(filNum).eventlabels)+1}='L-';
    drg.session(filNum).events(length(drg.session(filNum).eventlabels)).times=sort([drg.session(filNum).events(strcmp(drg.session(filNum).eventlabels,'CR')).times ...
        drg.session(filNum).events(strcmp(drg.session(filNum).eventlabels,'Miss')).times]);
    drg.session(filNum).events(length(drg.session(filNum).eventlabels)).noTimes=length(sort([drg.session(filNum).events(strcmp(drg.session(filNum).eventlabels,'CR')).times ...
        drg.session(filNum).events(strcmp(drg.session(filNum).eventlabels,'Miss')).times]));
end

%Now if there are S+ and S- add 'S+TStart' and 'S-TStart'
%S+TStart
if ~isempty(find(strcmp('S+',drg.session(filNum).eventlabels), 1))
    drg.session(filNum).eventlabels{length(drg.session(filNum).eventlabels)+1}='S+TStart';
    drg.session(filNum).events(length(drg.session(filNum).eventlabels)).noTimes=drg.session(filNum).events(find(strcmp('S+',drg.session(filNum).eventlabels))).noTimes;
    for ii=1:drg.session(filNum).events(find(strcmp('S+',drg.session(filNum).eventlabels))).noTimes
        [mindeltat,minjj]=min(abs(drg.session(filNum).events(1).times-drg.session(filNum).events(find(strcmp('S+',drg.session(filNum).eventlabels))).times(ii)));
        drg.session(filNum).events(length(drg.session(filNum).eventlabels)).times(ii)=drg.session(filNum).events(1).times(minjj);
    end
end

if ~isempty(find(strcmp('S-',drg.session(filNum).eventlabels), 1))
    %S-TStart
    drg.session(filNum).eventlabels{length(drg.session(filNum).eventlabels)+1}='S-TStart';
    drg.session(filNum).events(length(drg.session(filNum).eventlabels)).noTimes=drg.session(filNum).events(find(strcmp('S-',drg.session(filNum).eventlabels))).noTimes;
    for ii=1:drg.session(filNum).events(find(strcmp('S-',drg.session(filNum).eventlabels))).noTimes
        [mindeltat,minjj]=min(abs(drg.session(filNum).events(1).times-drg.session(filNum).events(find(strcmp('S-',drg.session(filNum).eventlabels))).times(ii)));
        drg.session(filNum).events(length(drg.session(filNum).eventlabels)).times(ii)=drg.session(filNum).events(1).times(minjj);
    end
end

%Add Odor1or2 of cm
switch drg.block_choice
    case 1
        
    case 2
        %These are cm trials. Add Odor1or2-S+
        drg.session(filNum).eventlabels{length(drg.session(filNum).eventlabels)+1}='Odor1or2-S+';
        drg.session(filNum).events(length(drg.session(filNum).eventlabels)).times=sort([drg.session(filNum).events(strcmp(drg.session(filNum).eventlabels,'Odor1-S+')).times...
            drg.session(filNum).events(strcmp(drg.session(filNum).eventlabels,'Odor2-S+')).times]);
        
    case 3
end

%Now redo blocks.
drg.session(filNum).blocks=draq_d.blocks;
switch drg.block_choice
    case 1
        %Leave the blocks alone
        
    case 2
        %Redo block using cm block logic (20 S- trials, and 10 S+ trials
        %for two odors)
        drg.session(filNum).blocks(1,1)=drg.session(filNum).trial_start(1)-100;
        indx=20;
        inbl=2;
        while (indx<length(drg.session(filNum).events(strcmp(drg.session(filNum).eventlabels,'S-')==1).times))
            drg.session(filNum).blocks(inbl,1)=0.001+(drg.session(filNum).events(strcmp(drg.session(filNum).eventlabels,'S-')==1).times(indx)...
                +drg.session(filNum).events(strcmp(drg.session(filNum).eventlabels,'S-')==1).times(indx+1))/2;
            indx=indx+20;
            inbl=inbl+1;
        end
        indx=20;
        inbl=1;
        while (indx<length(drg.session(filNum).events(strcmp(drg.session(filNum).eventlabels,'S-')==1).times))
            drg.session(filNum).blocks(inbl,2)=-0.001+(drg.session(filNum).events(strcmp(drg.eventlabels,'S-')==1).times(indx)...
                +drg.session(filNum).events(strcmp(drg.session(filNum).eventlabels,'S-')==1).times(indx+1))/2;
            indx=indx+20;
            inbl=inbl+1;
        end
        drg.session(filNum).blocks(inbl,2)=drg.session(filNum).trial_start(end)+100;
        
    case 3
        %Sort into a single block
        drg.session(filNum).blocks(1,1)=drg.session(filNum).trial_start(1)-100;
        drg.session(filNum).blocks(1,2)=max(drg.session(filNum).events(1).times)+100;
end


%If the draq_p.dgordra does not exist this is a Wilder file (dra)
try
    drg.session(filNum).dgordra=draq_p.dgordra;
    drg.dgordra=draq_p.dgordra;
catch
    drg.session(filNum).dgordra=1;
    drg.dgordra=1;
end

drg.session(filNum).draq_p=draq_p;
drg.session(filNum).draq_d=draq_d;

numUnits=0;

switch drg.session(filNum).dgordra
    case {1,4}
        %this is dra or Plexon
        drg.draq_p.no_chans=16;
        drg.session(filNum).noUnits=0;
        for chNo=1:length(noSpikes)
            if drta_p.ch_processed(chNo)==1
                for clusNo=1:max(cluster_class_per_file(offset_for_chan(chNo)+1:offset_for_chan(chNo)+noSpikes(chNo)));
                    this_clus=find(cluster_class_per_file(offset_for_chan(chNo)+1:offset_for_chan(chNo)+noSpikes(chNo))==clusNo);
                    if ~isempty(this_clus)
                        szts=size(this_clus);
                        numUnits=numUnits+1;
                        drg.session(filNum).noUnits=drg.session(filNum).noUnits+1;
                        %                             drg.spikes(numUnits+drg.noUnits).time=all_timestamp_per_file(this_clus+offset_for_chan(chNo));
                        %                             drg.ch_un{numUnits+drg.noUnits}=['ch' num2str(chNo) 'u' num2str(clusNo)];
                        %                             drg.channel(numUnits+drg.noUnits)=chNo;
                        %                             drg.unitinch(numUnits+drg.noUnits)=clusNo;
                        %                             drg.array(numUnits+drg.noUnits)=floor((chNo-1)/8);       %Assumes the arrays are eight channel arrays
                        %                             drg.sessionNo(numUnits+drg.noUnits)=filNum;
                        
                        %Really, this should be specified per unit
                        drg.unit(numUnits+drg.noUnits).spike_times=all_timestamp_per_file(this_clus+offset_for_chan(chNo));
                        drg.unit(numUnits+drg.noUnits).ch_un=['ch' num2str(chNo) 'u' num2str(clusNo)];
                        drg.unit((numUnits+drg.noUnits)).channel=chNo;
                        drg.unit(numUnits+drg.noUnits).unitinch=clusNo;
                        drg.unit(numUnits+drg.noUnits).array=floor((chNo-1)/8);       %Assumes the arrays are eight channel arrays
                        drg.unit(numUnits+drg.noUnits).sessionNo=filNum;
                        drg.unit(numUnits+drg.noUnits).dgordra=1;
                        drg.unit(numUnits+drg.noUnits).blocks=draq_d.blocks;
                        
                    end
                    %end
                end
            end
        end
        
    case {2,3}
        %this is dg or rhd
        drg.session(filNum).noUnits=0;
        for chNo=1:4
            if isfield(drta_p,'tetr_processed')
                spikes_processed=1;
                if drta_p.tetr_processed(chNo)==1
                    max_clusNo=max(cluster_class_per_file(offset_for_chan(chNo)+1:offset_for_chan(chNo)+noSpikes(chNo)));
                    for clusNo=1:max_clusNo
                        this_clus=find(cluster_class_per_file(offset_for_chan(chNo)+1:offset_for_chan(chNo)+noSpikes(chNo))==clusNo);
                        if ~isempty(this_clus)
                            szts=size(this_clus);
                            numUnits=numUnits+1;
                            drg.session(filNum).noUnits=drg.session(filNum).noUnits+1;
                            
                            
                            %Really, this should be specified per unit
                            drg.unit(numUnits+drg.noUnits).spike_times=all_timestamp_per_file(this_clus+offset_for_chan(chNo));
                            drg.unit(numUnits+drg.noUnits).ch_un=['ch' num2str(chNo) 'u' num2str(clusNo)];
                            drg.unit((numUnits+drg.noUnits)).channel=chNo;
                            drg.unit(numUnits+drg.noUnits).unitinch=clusNo;
                            drg.unit(numUnits+drg.noUnits).array=floor((chNo-1)/8);       %Assumes the arrays are eight channel arrays
                            drg.unit(numUnits+drg.noUnits).sessionNo=filNum;
                            drg.unit(numUnits+drg.noUnits).dgordra=0;
                            drg.unit(numUnits+drg.noUnits).blocks=draq_d.blocks;
                            try
                                if (clusNo<=max_clusNo)&exist('units_per_tet','var')
                                    drg.unit(numUnits+drg.noUnits).Lratio=units_per_tet(chNo).Lratio(clusNo+1);
                                    drg.unit(numUnits+drg.noUnits).IsolDist=units_per_tet(chNo).IsolDist(clusNo+1);
                                else
                                    drg.unit(numUnits+drg.noUnits).Lratio=NaN;
                                    drg.unit(numUnits+drg.noUnits).IsolDist=NaN;
                                end
                            catch
                            end
                            
                        end
                        %end
                    end
                end
            end
        end
    
        
end

%Save the user choices for which trials, channels, trialsxchannel to
%process
drg.session(filNum).channels_processed=drta_p.ch_processed;
drg.session(filNum).trials_x_channel_processed=drta_p.trial_ch_processed;
drg.session(filNum).trials_processed=drta_p.trial_allch_processed;
if isfield(drta_p,'doSubtract')
    drg.session(filNum).doSubtract=drta_p.doSubtract;
    drg.session(filNum).subtractCh=drta_p.subtractCh;
else
    drg.session(filNum).doSubtract=0;
end
%Save the total number of seconds per trial and the start time for each
%trial
drg.session(filNum).start_times=draq_d.t_trial;
drg.session(filNum).sec_per_trial=draq_p.sec_per_trigger;

%Number of trials per session
drg.session(filNum).noTrials=drg.draq_d.noTrials;
drg.session(filNum).no_chans=drg.draq_p.no_chans;

%Calculate inter spike intervals (ISI) and classify each unit as a single or multi
%unit



%     for unitNo=1:numUnits
%         ISIs=[];
%         ISIs=drg.spikes(unitNo+drg.noUnits).time(2:end)-drg.spikes(unitNo+drg.noUnits).time(1:end-1);
%         drg.perViol(unitNo+drg.noUnits)=100*sum(ISIs<refractPer)/length(ISIs);
%         if drg.perViol(unitNo+drg.noUnits)<maxPercentViol
%             drg.SingleUnit(unitNo+drg.noUnits)=1;
%         else
%             drg.SingleUnit(unitNo+drg.noUnits)=0;
%         end
%     end




%Extract the PSTH per trial, making sure to exclude the trials excluded by
%the user
%Do all OdorOn events (event 2), centered on odor on
%Also get inter trial intervals and classify unit as single or multi unit
%     evTypeNo=2;
%     noEvents=length(drg.session(filNum).events(evTypeNo).times);
%     szblocks=size(drg.session(filNum).blocks);
for unitNo=1:numUnits
    drg.unit(drg.noUnits+unitNo).session=filNum;
    
    
    drg.unit(drg.noUnits+unitNo).ISIs=[];
    drg.unit(drg.noUnits+unitNo).noISIs=0;
    %         spikes1=[];
    %         spikes1=drg.unit(drg.noUnits+unitNo).spike_times;
    spike_times=drg.unit(drg.noUnits+unitNo).spike_times;
    drg.unit(drg.noUnits+unitNo).noTrials=0;
    %Now do all trials
    no_evs_included=0;
    basalFR=[];
    
    
    %         for evNo=1:noEvents
    %             %Find out if this trial should be used
    %
    %             if drgExcludeTrial(drg,drg.unit((drg.noUnits+unitNo)).channel,drg.session(filNum).events(evTypeNo).times(evNo),filNum)==0
    %                 drg.unit(drg.noUnits+unitNo).noTrials=drg.unit(drg.noUnits+unitNo).noTrials+1;
    %                 perieventindx=[];
    %                 perieventSpikes=[];
    %                 %perieventindx=(spikes1>(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre))&(spikes1<(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_post));
    %                 perieventindx=(spikes1>(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre))&(spikes1<(drg.session(filNum).events(evTypeNo).times(evNo)));
    %                 perieventSpikes=spikes1(perieventindx)-drg.session(filNum).events(evTypeNo).times(evNo);
    %                 drg.unit(drg.noUnits+unitNo).ISIs=[drg.unit(drg.noUnits+unitNo).ISIs perieventSpikes(2:end)-perieventSpikes(1:end-1)];
    %                 drg.unit(drg.noUnits+unitNo).noISIs=drg.unit(drg.noUnits+unitNo).noISIs+sum(perieventindx)-1;
    %
    %                 basal_indx=(spikes1>(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre))&(spikes1<drg.session(filNum).events(evTypeNo).times(evNo));
    %                 no_evs_included=no_evs_included+1;
    %                 basalFR(no_evs_included)=sum(basal_indx)/(-drg.time_pre);
    %             end
    %
    %         end
    
    %Do OdorOn
    evTypeNo=2; %This is odorOn
    
    for evNo=1:drg.session(filNum).events(evTypeNo).noTimes
        
        %Get an S+ trial
        excludeTrial=drgExcludeTrial(drg,drg.unit((drg.noUnits+unitNo)).channel,drg.session(filNum).events(evTypeNo).times(evNo),filNum);
        
        if excludeTrial==0
            
            drg.unit(drg.noUnits+unitNo).noTrials=drg.unit(drg.noUnits+unitNo).noTrials+1;
            %                 perieventindx=[];
            %                 perieventSpikes=[];
            %                 %perieventindx=(spikes1>(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre))&(spikes1<(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_post));
            %                 perieventindx=(spikes1>(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre))&(spikes1<(drg.session(filNum).events(evTypeNo).times(evNo)));
            %                 perieventSpikes=spikes1(perieventindx)-drg.session(filNum).events(evTypeNo).times(evNo);
            %                 drg.unit(drg.noUnits+unitNo).ISIs=[drg.unit(drg.noUnits+unitNo).ISIs perieventSpikes(2:end)-perieventSpikes(1:end-1)];
            %                 drg.unit(drg.noUnits+unitNo).noISIs=drg.unit(drg.noUnits+unitNo).noISIs+sum(perieventindx)-1;
            %                 basal_indx=[];
            %                 basal_indx=(spikes1>(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre))&(spikes1<drg.session(filNum).events(evTypeNo).times(evNo));
            %                 no_evs_included=no_evs_included+1;
            %                 basalFR(no_evs_included)=sum(basal_indx)/(-drg.time_pre);
            noTrials=drg.unit(drg.noUnits+unitNo).noTrials;
            
            these_spikes=[];
            these_spikes=(spike_times>drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre)&...
                (spike_times<=drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_post);
            these_spike_times=[];
            these_spike_times=spike_times(these_spikes)-(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre);
            last_pre=find(these_spike_times<-drg.time_pre,1,'last');
            
            
            PSTHSminus=zeros(1,drg.nobins);
            for spk=1:length(these_spike_times)
                this_bin=ceil(these_spike_times(spk)/drg.bin_size);
                PSTHSminus(1,this_bin)=PSTHSminus(1,this_bin)+1;
            end %for spk
            
            
            PSTHSminus=PSTHSminus/drg.bin_size;
            
            %Now enter the BFR per unit
            basalFR=[basalFR mean(PSTHSminus(1,drg.start_pre:drg.end_pre))];
            
            if ~isempty(last_pre)
                this_ISI=[];
                this_ISI=these_spike_times(2:last_pre)-these_spike_times(1:last_pre-1);
                drg.unit(drg.noUnits+unitNo).ISIs=[drg.unit(drg.noUnits+unitNo).ISIs this_ISI(this_ISI>0)];
                drg.unit(drg.noUnits+unitNo).noISIs=drg.unit(drg.noUnits+unitNo).noISIs+sum(this_ISI>0);
            end
            
        end
        
    end
   
    %
    %         %Do S+
    %         evTypeNo=5; %This is S+
    %         for evNo=1:drg.session(filNum).events(evTypeNo).noTimes
    %
    %             %Get an S+ trial
    %             excludeTrial=drgExcludeTrial(drg,drg.unit((drg.noUnits+unitNo)).channel,drg.session(filNum).events(evTypeNo).times(evNo),filNum);
    %
    %             if excludeTrial==0
    %
    %                 drg.unit(drg.noUnits+unitNo).noTrials=drg.unit(drg.noUnits+unitNo).noTrials+1;
    %                 %                 perieventindx=[];
    %                 %                 perieventSpikes=[];
    %                 %                 %perieventindx=(spikes1>(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre))&(spikes1<(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_post));
    %                 %                 perieventindx=(spikes1>(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre))&(spikes1<(drg.session(filNum).events(evTypeNo).times(evNo)));
    %                 %                 perieventSpikes=spikes1(perieventindx)-drg.session(filNum).events(evTypeNo).times(evNo);
    %                 %                 drg.unit(drg.noUnits+unitNo).ISIs=[drg.unit(drg.noUnits+unitNo).ISIs perieventSpikes(2:end)-perieventSpikes(1:end-1)];
    %                 %                 drg.unit(drg.noUnits+unitNo).noISIs=drg.unit(drg.noUnits+unitNo).noISIs+sum(perieventindx)-1;
    %                 %                 basal_indx=[];
    %                 %                 basal_indx=(spikes1>(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre))&(spikes1<drg.session(filNum).events(evTypeNo).times(evNo));
    %                 %                 no_evs_included=no_evs_included+1;
    %                 %                 basalFR(no_evs_included)=sum(basal_indx)/(-drg.time_pre);
    %                 noTrials=drg.unit(drg.noUnits+unitNo).noTrials;
    %
    %                 these_spikes=[];
    %                 these_spikes=(spike_times>drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre)&...
    %                     (spike_times<=drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_post);
    %                 these_spike_times=[];
    %                 these_spike_times=spike_times(these_spikes)-(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre);
    %                 last_pre=find(these_spike_times<-drg.time_pre,1,'last');
    %
    %
    %                 PSTHSminus=zeros(1,drg.nobins);
    %                 for spk=1:length(these_spike_times)
    %                     this_bin=ceil(these_spike_times(spk)/drg.bin_size);
    %                     PSTHSminus(1,this_bin)=PSTHSminus(1,this_bin)+1;
    %                 end %for spk
    %
    %
    %                 PSTHSminus=PSTHSminus/drg.bin_size;
    %
    %                 %Now enter the BFR per unit
    %                 basalFR=[basalFR mean(PSTHSminus(1,drg.start_pre:drg.end_pre))];
    %
    %                 if ~isempty(last_pre)
    %                     this_ISI=[];
    %                     this_ISI=these_spike_times(2:last_pre)-these_spike_times(1:last_pre-1);
    %                     drg.unit(drg.noUnits+unitNo).ISIs=[drg.unit(drg.noUnits+unitNo).ISIs this_ISI(this_ISI>0)];
    %                     drg.unit(drg.noUnits+unitNo).noISIs=drg.unit(drg.noUnits+unitNo).noISIs+sum(this_ISI>0);
    %                 end
    %
    %             end
    %
    %         end
    %
    %         %Do S-
    %         evTypeNo=11; %This is S-
    %         for evNo=1:drg.session(filNum).events(evTypeNo).noTimes
    %
    %             %Get an S- trial
    %             excludeTrial=drgExcludeTrial(drg,drg.unit((drg.noUnits+unitNo)).channel,drg.session(filNum).events(evTypeNo).times(evNo),filNum);
    %
    %             if excludeTrial==0
    %                 drg.unit(drg.noUnits+unitNo).noTrials=drg.unit(drg.noUnits+unitNo).noTrials+1;
    %                 %                 perieventindx=[];
    %                 %                 perieventSpikes=[];
    %                 %                 %perieventindx=(spikes1>(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre))&(spikes1<(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_post));
    %                 %                 perieventindx=(spikes1>(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre))&(spikes1<(drg.session(filNum).events(evTypeNo).times(evNo)));
    %                 %                 perieventSpikes=spikes1(perieventindx)-drg.session(filNum).events(evTypeNo).times(evNo);
    %                 %                 drg.unit(drg.noUnits+unitNo).ISIs=[drg.unit(drg.noUnits+unitNo).ISIs perieventSpikes(2:end)-perieventSpikes(1:end-1)];
    %                 %                 drg.unit(drg.noUnits+unitNo).noISIs=drg.unit(drg.noUnits+unitNo).noISIs+sum(perieventindx)-1;
    %                 %                 basal_indx=[];
    %                 %                 basal_indx=(spikes1>(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre))&(spikes1<drg.session(filNum).events(evTypeNo).times(evNo));
    %                 %                 no_evs_included=no_evs_included+1;
    %                 %                 basalFR(no_evs_included)=sum(basal_indx)/(-drg.time_pre);
    %                 noTrials=drg.unit(drg.noUnits+unitNo).noTrials;
    %
    %                 these_spikes=[];
    %                 these_spikes=(spike_times>drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre)&...
    %                     (spike_times<=drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_post);
    %                 these_spike_times=[];
    %                 these_spike_times=spike_times(these_spikes)-(drg.session(filNum).events(evTypeNo).times(evNo)+drg.time_pre);
    %                 last_pre=find(these_spike_times<-drg.time_pre,1,'last');
    %
    %
    %                 PSTHSminus=zeros(1,drg.nobins);
    %                 for spk=1:length(these_spike_times)
    %                     this_bin=ceil(these_spike_times(spk)/drg.bin_size);
    %                     PSTHSminus(1,this_bin)=PSTHSminus(1,this_bin)+1;
    %                 end %for spk
    %
    %
    %                 PSTHSminus=PSTHSminus/drg.bin_size;
    %
    %                 %Now enter the BFR per unit
    %                 basalFR=[basalFR mean(PSTHSminus(1,drg.start_pre:drg.end_pre))];
    %
    %                 if ~isempty(last_pre)
    %                     this_ISI=[];
    %                     this_ISI=these_spike_times(2:last_pre)-these_spike_times(1:last_pre-1);
    %                     drg.unit(drg.noUnits+unitNo).ISIs=[drg.unit(drg.noUnits+unitNo).ISIs this_ISI(this_ISI>0)];
    %                     drg.unit(drg.noUnits+unitNo).noISIs=drg.unit(drg.noUnits+unitNo).noISIs+sum(this_ISI>0);
    %                 end
    %
    %             end
    %
    %         end
    
    drg.unit(drg.noUnits+unitNo).basalFR=mean(basalFR);
    drg.basal_firing_rate(drg.noUnits+unitNo)=mean(basalFR);
    
    %Determine whether these are single or multiunits
    %Method 0
    %These are the old values, taken from Quian Quiroga
    %         refractPer=0.002;       %Refractory period in seconds. We use the value from Quian Quiroga
    %         maxPercentViol=3;       %Maximum percent ISI violations
    %
    %         drg.unit(drg.noUnits+unitNo).perViol=100*sum((drg.unit(drg.noUnits+unitNo).ISIs<refractPer)/drg.unit(drg.noUnits+unitNo).noISIs);
    %
    %         %These are the new values (the calculation is different, see below
    %         %Method 1
    %         refractPer=0.001; %Refractory period in seconds.
    %         minRefPeriod=0.0015;
    %         maxRefPeriod=0.0035;
    %         maxPercentViol=10;     %Maximum percent ISI violations
    %
    %         drg.unit(drg.noUnits+unitNo).perViol=100*sum((drg.unit(drg.noUnits+unitNo).ISIs<refractPer)/refractPer)/...
    %             (sum((drg.unit(drg.noUnits+unitNo).ISIs>minRefPeriod)&(drg.unit(drg.noUnits+unitNo).ISIs<=maxRefPeriod))/(maxRefPeriod-minRefPeriod));
    %
    
    % Refractory period and percent violation
    drg.refractPer=0.001;       %Refractory period in seconds.
    drg.maxPercentViol=0.75;       %Maximum percent ISI violations
    
    drg.unit(drg.noUnits+unitNo).perViol=100*sum((drg.unit(drg.noUnits+unitNo).ISIs<drg.refractPer)/drg.unit(drg.noUnits+unitNo).noISIs);
    
    
    if drg.unit(drg.noUnits+unitNo).perViol<=drg.maxPercentViol
        drg.unit(drg.noUnits+unitNo).SingleUnit=1;
        drg.single_or_multi(drg.noUnits+unitNo)=1;
        no_singles=no_singles+1;
    else
        drg.unit(drg.noUnits+unitNo).SingleUnit=0;
        drg.single_or_multi(drg.noUnits+unitNo)=0;
    end
    
    %Now calculate R_2_20 as in Fee et al J. Neurosci Meth. 69:175,
    %1996. We use 20 msec because of Stark et al Neuron 80:1263, 2013
    drg.Fee_refractPer=0.002;
    drg.Fee_no_detection=0.0005;
    drg.Fee_twenty=0.02;
    drg.Fee_criterion=0.4;
    
    F_2=sum((drg.unit(drg.noUnits+unitNo).ISIs<drg.Fee_refractPer)&(drg.unit(drg.noUnits+unitNo).ISIs>drg.Fee_no_detection));
    F_20=sum((drg.unit(drg.noUnits+unitNo).ISIs>drg.Fee_refractPer)&(drg.unit(drg.noUnits+unitNo).ISIs<drg.Fee_twenty));
    drg.unit(drg.noUnits+unitNo).R_2_20_Fee=((drg.Fee_twenty-drg.Fee_refractPer)/(drg.Fee_refractPer-drg.Fee_no_detection))*F_2/F_20;
    
    if drg.unit(drg.noUnits+unitNo).R_2_20_Fee<=drg.Fee_criterion
        drg.unit(drg.noUnits+unitNo).SingleUnit_Fee=1;
        drg.single_or_multi_Fee(drg.noUnits+unitNo)=1;
        no_singles_Fee=no_singles_Fee+1;
    else
        drg.unit(drg.noUnits+unitNo).SingleUnit_Fee=0;
        drg.single_or_multi_Fee(drg.noUnits+unitNo)=0;
    end
end

drg.noUnits=drg.noUnits+numUnits;
%end
percent_single_units=100*no_singles/drg.noUnits;
no_units=drg.noUnits;

% %Plot BFR histogram
% figure(4)
% subplot(2,1,1);
% edges=[0:5:300];
% hist(drg.basal_firing_rate(drg.single_or_multi==1),edges)
% mean_single=mean(drg.basal_firing_rate(drg.single_or_multi==1));
% SD_single=std(drg.basal_firing_rate(drg.single_or_multi==1));
% num_single=sum(drg.single_or_multi==1);
% xlabel('Frequency (Hz)')
% title('Basal firing rate for singles (Hz)')
% xlim([0 300])
% subplot(2,1,2);
% hist(drg.basal_firing_rate(drg.single_or_multi==0),edges)
% mean_multi=mean(drg.basal_firing_rate(drg.single_or_multi==0));
% SD_multi=std(drg.basal_firing_rate(drg.single_or_multi==0));
% num_multi=sum(drg.single_or_multi==0);
% title ('Basal firing rate for multis (Hz)')
% xlabel('Frequency (Hz)')
% xlim([0 300])
% singles_per_multi=mean_multi/mean_single;
% 
% %Plot BFR histogram Fee
% figure(5)
% subplot(2,1,1);
% edges=[0:5:300];
% hist(drg.basal_firing_rate(drg.single_or_multi_Fee==1),edges)
% mean_single_Fee=mean(drg.basal_firing_rate(drg.single_or_multi_Fee==1));
% SD_single_Fee=std(drg.basal_firing_rate(drg.single_or_multi_Fee==1));
% num_single_Fee=sum(drg.single_or_multi_Fee==1);
% xlabel('Frequency (Hz)')
% title('Basal firing rate for singles (Hz), R Fee')
% xlim([0 300])
% subplot(2,1,2);
% hist(drg.basal_firing_rate(drg.single_or_multi_Fee==0),edges)
% mean_multi_Fee=mean(drg.basal_firing_rate(drg.single_or_multi_Fee==0));
% SD_multi_Fee=std(drg.basal_firing_rate(drg.single_or_multi_Fee==0));
% num_multi_Fee=sum(drg.single_or_multi_Fee==0);
% title ('Basal firing rate for multis (Hz), R Fee')
% xlabel('Frequency (Hz)')
% xlim([0 300])
% singles_per_multi_Fee=mean_multi_Fee/mean_single_Fee;
if spikes_processed==0
   drg.unit(1).sessionNo=1; 
end

%Save drg
save([drg.drg_directory drg.drg_output_file],'drg','-v7.3')


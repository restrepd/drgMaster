function drgPlotPSTHEncRetrFct(handles)
%PSTH plot the encoding and retrieval segments


[perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);
%Enter unit and event
%unitNo=2;

try
    close 1
catch
end
 
figure(1)
 
%Plot the percent correct
subplot(3,1,1)
trials=1:length(perCorr);
plot(trials,perCorr,'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
hold on
plot(trials(encoding_trials),perCorr(encoding_trials),'ob')
plot(trials(retrieval_trials),perCorr(retrieval_trials),'or')
ylim([40 110]);
xlabel('Trial #')
ylabel('Percent correct')
title('% correct vs trial number, blue=encoding, red=retrieval')
 

%Do encoding
sessionNo=handles.drg.unit(handles.unitNo).sessionNo;
bin_size=0.02;
nobins=fix((handles.drg.time_post-handles.drg.time_pre)/bin_size); 
itime=1:nobins;
itime=itime+fix(handles.drg.time_pre/bin_size);
time=double(itime)*bin_size;
noTrials_encoding=0;
spike_times=[];
spike_times=handles.drg.unit(handles.unitNo).spike_times;
PSTHencoding=zeros(1,nobins);


for evNo=encoding_this_evTypeNo
    
    %evNo
    excludeTrial=drgExcludeTrial(handles.drg,handles.drg.unit(handles.unitNo).channel,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
    
    if excludeTrial==0
        
        noTrials_encoding=noTrials_encoding+1;
        these_spikes=(spike_times>handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+handles.drg.time_pre)&...
            (spike_times<=handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+handles.drg.time_post);
        these_spike_times=spike_times(these_spikes)-(handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+handles.drg.time_pre);
        
        for spk=1:length(these_spike_times)
            this_bin=ceil(these_spike_times(spk)/bin_size);
            PSTHencoding(1,this_bin)=PSTHencoding(1,this_bin)+1;
        end %for spk
    end
        %end
    %end %if eventstamps...
end %for evNo

number_of_trials_encoding=noTrials_encoding

PSTHencoding=PSTHencoding/(noTrials_encoding*bin_size);

%Now plot the PSTHencoding
subplot(3,1,2)
bar(time,PSTHencoding,'b');
title(['PSTH, encoding segment, for ' handles.drg.session.eventlabels{handles.evTypeNo}])
ylabel('Frequency (Hz)')
xlabel('Time (sec)')

%Do retrieval
sessionNo=handles.drg.unit(handles.unitNo).sessionNo;
bin_size=0.02;
nobins=fix((handles.drg.time_post-handles.drg.time_pre)/bin_size); 
itime=1:nobins;
itime=itime+fix(handles.drg.time_pre/bin_size);
time=double(itime)*bin_size;
noTrials_retrieval=0;
spike_times=[];
spike_times=handles.drg.unit(handles.unitNo).spike_times;
PSTHretrieval=zeros(1,nobins);


for evNo=retrieval_this_evTypeNo
    
    %evNo
    excludeTrial=drgExcludeTrial(handles.drg,handles.drg.unit(handles.unitNo).channel,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
    
    if excludeTrial==0
        
        noTrials_retrieval=noTrials_retrieval+1;
        these_spikes=(spike_times>handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+handles.drg.time_pre)&...
            (spike_times<=handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+handles.drg.time_post);
        these_spike_times=spike_times(these_spikes)-(handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+handles.drg.time_pre);
        
        for spk=1:length(these_spike_times)
            this_bin=ceil(these_spike_times(spk)/bin_size);
            PSTHretrieval(1,this_bin)=PSTHretrieval(1,this_bin)+1;
        end %for spk
    end
        %end
    %end %if eventstamps...
end %for evNo

number_retreival_trials=noTrials_retrieval

PSTHretrieval=PSTHretrieval/(noTrials_retrieval*bin_size);

%Now plot the PSTHretrieval
subplot(3,1,3)
bar(time,PSTHretrieval,'b');
title(['PSTH, retrieval segment, for ' handles.drg.session.eventlabels{handles.evTypeNo}])
ylabel('Frequency (Hz)')
xlabel('Time (sec)')
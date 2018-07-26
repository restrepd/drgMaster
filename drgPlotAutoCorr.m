function drgPlotAutoCorr(handles)
%Plot autocorrelaogram

noShuffles=5;
sessionNo=handles.drg.unit(handles.unitNo).sessionNo;

%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;

bin_size=0.001;
auto_width=0.02;
nobins=fix(2*(auto_width/bin_size));
delta_times=[-auto_width+bin_size/2:bin_size:auto_width-bin_size/2];
noTrials=0;
spike_times=[];
spike_times=handles.drg.unit(handles.unitNo).spike_times;
Auto=zeros(1,nobins);
no_comp_spikes=0;
shAuto=zeros(1,nobins);
no_comp_spikes_sh=0;
time_start=handles.time_start+handles.time_pad;
time_end=handles.time_end-handles.time_pad;

for trNo=firstTr:lastTr
    
    if handles.save_drgb==0
        trial_no=trNo
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    
    if evNo~=-1
        
        %evNo
        excludeTrial=drgExcludeTrial(handles.drg,handles.drg.unit(handles.unitNo).channel,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            noTrials=noTrials+1;
            these_spikes=(spike_times>handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+time_start+auto_width)&...
                (spike_times<=handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+time_end-auto_width);
            these_spike_times=spike_times(these_spikes)-(handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)+time_start);
            
            for spkref=1:length(these_spike_times)
                for spk=1:length(these_spike_times)
                    if spkref~=spk
                        deltat=these_spike_times(spk)-these_spike_times(spkref);
                        if abs(deltat)<auto_width
                            this_bin=fix(deltat/bin_size)+floor(auto_width/bin_size);
                            Auto(1,this_bin)=Auto(1,this_bin)+1;
                        end
                    end
                end
            end %for spkref
            no_comp_spikes=no_comp_spikes+length(these_spike_times)-1;
            
            %Now do random shuffled spikes
            for noS=1:noShuffles
                shuf_spike_times=(time_end-time_start)*rand(1,length(these_spike_times));
                for spkref=1:length(these_spike_times)
                    for spk=1:length(these_spike_times)
                        if spkref~=spk
                            deltat=shuf_spike_times(spk)-shuf_spike_times(spkref);
                            if abs(deltat)<auto_width
                                this_bin=fix(deltat/bin_size)+floor(auto_width/bin_size);
                                shAuto(1,this_bin)=shAuto(1,this_bin)+1;
                            end
                        end
                    end
                end %for spkref
                no_comp_spikes_sh=no_comp_spikes_sh+length(these_spike_times)-1;
            end
        end
        %end
        %end %if eventstamps...
    end %if evNo
end %for trNo=

Auto=Auto/no_comp_spikes;
shAuto=shAuto/no_comp_spikes_sh;

%Now plot the autocorrelogram
try
    close 1
catch
end

%Plot the timecourse
hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.02 .4 .5 .5])

subplot(2,1,1)
bar(delta_times,Auto,'b');
title(['Autocorrelogram for ' handles.drg.session.eventlabels{handles.evTypeNo}])
ylabel('Correlation coefficient')
xlabel('delta time (sec)')

%Now plot the autocorrelogram - random
subplot(2,1,2)
bar(delta_times,Auto-shAuto,'b');
title(['Autocorrelogram -random for ' handles.drg.session.eventlabels{handles.evTypeNo}])
ylabel('Correlation coefficient')
xlabel('delta time (sec)')
y_max=max(Auto-shAuto);
ylim([0 1.2*y_max])

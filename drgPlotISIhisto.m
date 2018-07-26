function drgPlotISIhisto(handles)
%Plot autocorrelaogram


sessionNo=handles.drg.unit(handles.unitNo).sessionNo;

%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;

skip_ISI=0;


auto_width=0.05;
noTrials=0;
spike_times=handles.drg.unit(handles.unitNo).spike_times;
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
            
            at_end=0;
            
            %This routine produces a histogram discarding close spikes (for
            %use with complex cerebellar spikes)
            
            ii=0;
            ISIs=[];
            while at_end==0
                
                if ii+1<length(these_spike_times)-1
                    ii=ii+1;
                    jj=ii+1;
                    %find the next ISI (if there is one)
                    jj_found=0;
                    while (at_end==0)&(jj_found==0)

                        if jj<=length(these_spike_times)
                            if these_spike_times(jj)-these_spike_times(ii)>skip_ISI
                                ISIs=[ISIs these_spike_times(jj)-these_spike_times(ii)];
                                jj_found=1;
                                ii=jj-1;
                            else
                                jj=jj+1;
                            end
                        else
                            at_end=1;
                        end
                        
                    end
                else
                    at_end=1;
                end
                
            end
            
        end
        %end
        %end %if eventstamps...
    end %if evNo
end %for trNo=



%Now plot the autocorrelogram
try
    close 1
catch
end


hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.02 .4 .5 .5])
edges=[0:0.2:10];
histogram(ISIs*1000,edges)

title('ISI histogram')
xlabel('ISI (ms)')
ylabel('Counts/bin)')

ISImedian=median(ISIs*1000)

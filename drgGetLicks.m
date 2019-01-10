function [lick_freq,times_lick_freq,lick_traces,CIlickf,lick_trace_times,stamped_lick_ii,these_stamped_lick_times,no_trials,trials_included]=drgGetLicks(handles)

%Get licks for those trials that have an LFP that can be read
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;

%First find the time range for the spectrogram, calculate using
%handles.time_pad to exclude artifacts at the end of
min_t=handles.time_start;
max_t=handles.time_end;

%Calculate the threshold value to detect a lick
sessionNo=1;
all_licks=[];
evtNo=handles.evTypeNo;
handles.evTypeNo=2;
hpNo=handles.peakLFPNo;
handles.peakLFPNo=19;
for trNo=1:handles.drg.session(1).draq_d.noTrials
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            %Note: handles.peakLFPNo is the reference LFP
            [lickLFP, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.startRef, handles.time_end);
            
            if (can_read1==1)
                all_licks=[all_licks lickLFP];
            end
        end
    end
end
handles.evTypeNo=evtNo;

thershold_licks=prctile(all_licks,1)+((prctile(all_licks,99)-prctile(all_licks,1))/2);


no_trials=0;
trials_included=[];
temp_stamped_lick_ii=zeros(1,lastTr-firstTr+1);
these_temp_stamped_lick_times=ones(lastTr-firstTr+1,250);
lick_traces=zeros(lastTr-firstTr+1,handles.drg.session(sessionNo).draq_p.ActualRate*(min_t-max_t));  %On purpose this is larger than the total number of points
dt_licks=handles.dt_lick;
lick_freq=zeros(1,ceil((handles.time_end-handles.time_start-2*handles.time_pad)/dt_licks)-1);
lick_freq_per_trial=zeros(lastTr-firstTr+1,ceil((handles.time_end-handles.time_start-2*handles.time_pad)/dt_licks)-1);
times_lick_freq=([1:(ceil((handles.time_end-handles.time_start-2*handles.time_pad)/dt_licks)-1)]*dt_licks)-(dt_licks/2)+min_t+handles.time_pad;

for trNo=firstTr:lastTr
    
    if handles.save_drgb==0
        trialNo=trNo
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    
    if evNo~=-1
 
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            %Can the LFP be read?
            [LFP, trialNo, can_read] = drgGetTrialLFPData(handles, hpNo, evNo, handles.evTypeNo, min_t, max_t);
            
            if (can_read==1)
                
                no_trials=no_trials+1;
                trials_included(no_trials)=trNo;
                
                %Get the licks
                [lickLFP, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, min_t, max_t);
                
                ii_start_licks=handles.time_pad*handles.drg.session(sessionNo).draq_p.ActualRate;
                ii_end_licks=length(lickLFP)-ii_start_licks;
                lick_trace=[];
                lick_trace=lickLFP(int32(ii_start_licks):int32(ii_end_licks));
                lick_traces(no_trials,1:length(lick_trace))=lick_trace;
            
                
                inter_lick_intervals_licks=[];
                ii_ili_licks=0;
                
                if lick_trace(1)>thershold_licks
                    ii=find(lick_trace<thershold_licks,1,'first');
                else
                    ii=1;
                end
                
                the_end=0;
                ref_power_these_events=[];
                no_licks_evs_this_trial=0;
                
                this_lick_ii=0;
                these_lick_times=[];
                
                
                %Find the events (licks)
                first_lick=1;
                while the_end==0
                    next_event=find(lick_trace(ii:end)>thershold_licks,1,'first');
                    if isempty(next_event)
                        the_end=1;
                    else
                        
                        ii=ii+next_event-1;
                        
                        %Exclude if the inter event interval is too
                        %small due to noise in the lick signal
                        
                        this_lick_ii=this_lick_ii+1;
                        these_lick_times(this_lick_ii)=(ii/handles.drg.session(sessionNo).draq_p.ActualRate);
                        if this_lick_ii>1
                            %record the inter lick interval
                            ii_ili_licks=ii_ili_licks+1;
                            inter_lick_intervals_licks(ii_ili_licks)=these_lick_times(this_lick_ii)-these_lick_times(this_lick_ii-1);
                        else
                            ii_ili_licks=ii_ili_licks+1;
                            inter_lick_intervals_licks(ii_ili_licks)=these_lick_times(this_lick_ii);
                        end
                        
                        %Enter the event (lick) in the timecourse only if it is
                        %not within a burst of high frequency noise
                        if inter_lick_intervals_licks(ii_ili_licks)>handles.smallest_inter_lick_interval
                            %stamp the lick
                            temp_stamped_lick_ii(no_trials)=temp_stamped_lick_ii(no_trials)+1;
                            these_temp_stamped_lick_times(no_trials,temp_stamped_lick_ii(no_trials))=(ii/handles.drg.session(sessionNo).draq_p.ActualRate)+min_t+handles.time_pad;
                        end

                        end_event=find(lick_trace(ii:end)<thershold_licks,1,'first');
                        if isempty(end_event)
                            the_end=1;
                        else
                            ii=ii+end_event-1;
                        end
                    end
                end
                
                %Find the lick frequency
                for jj=1:temp_stamped_lick_ii(no_trials)
                    if (these_temp_stamped_lick_times(no_trials,jj)>times_lick_freq(1)-(dt_licks/2))&...
                           (these_temp_stamped_lick_times(no_trials,jj)<times_lick_freq(end)+(dt_licks/2)) 
                        jj_lick_time=ceil((these_temp_stamped_lick_times(no_trials,jj)-(min_t+handles.time_pad))/dt_licks);
                        lick_freq(1,jj_lick_time)=lick_freq(1,jj_lick_time)+1;
                        lick_freq_per_trial(no_trials,jj_lick_time)=lick_freq_per_trial(no_trials,jj_lick_time)+1;
                    end
                end
                
                
            end
        end
    end
end
 
%Convolve lick_freq using a window
conv_win=ones(1,handles.window/dt_licks);
lick_freq=conv(lick_freq,conv_win,'same');
lick_freq=lick_freq/(no_trials*handles.window);

for ii=1:lastTr-firstTr+1
    lick_freq_per_trial(ii,:)=conv(lick_freq_per_trial(ii,:),conv_win,'same');
end
lick_freq_per_trial=lick_freq_per_trial/(dt_licks*handles.window);
lick_freq_per_trial=lick_freq_per_trial(1:no_trials,:);

CIlickf = bootci(1000, @mean, lick_freq_per_trial);
CIlickf(1,:)=mean(lick_freq_per_trial)-CIlickf(1,:);
CIlickf(2,:)=CIlickf(2,:)-mean(lick_freq_per_trial);

lick_traces=lick_traces(1:no_trials,1:ceil((handles.time_end-handles.time_start-2*handles.time_pad)*handles.drg.session(sessionNo).draq_p.ActualRate));
lick_trace_times=([1:ceil((handles.time_end-handles.time_start-2*handles.time_pad)*handles.drg.session(sessionNo).draq_p.ActualRate)]/handles.drg.session(sessionNo).draq_p.ActualRate)+min_t+handles.time_pad;
handles.peakLFPNo=hpNo;

stamped_lick_ii=zeros(1,no_trials);
these_stamped_lick_times=ones(no_trials,250);

stamped_lick_ii(1,:)=temp_stamped_lick_ii(1,1:no_trials);
these_stamped_lick_times(:,:)=these_temp_stamped_lick_times(1:no_trials,:);

pffft=1;



    
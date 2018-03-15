function [per_corr_per_trial, lick_timecourse,which_event,no_trials,inter_lick_intervals]=drgLickTimecourseThisEv(handles)
%Performs lick analysis. The lick event is signaled by a sharp change
%in the reference voltage. Events with small inter lick intervals are due to noise and are rejected

[perCorr, encoding_trials, retrieval_trials,encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);


if handles.displayData==1
    fprintf(1, '\n');
end

odorOn=2;
lickLFPNo=19;
which_event=[];
inter_lick_intervals=[];
ii_ili=0;
no_lick_dt=floor(((handles.time_end-handles.time_pad)-(handles.time_start+handles.time_pad))/handles.dt_lick);

%Generates a trial per trial phase histogram
sessionNo=handles.sessionNo;
no_time_pts=floor(handles.window*handles.drg.session(sessionNo).draq_p.ActualRate)+1;
lick_timecourse=[];
per_corr_per_trial=[];
max_events_per_trial=((handles.time_end-handles.time_pad)-(handles.time_start+handles.time_pad))*handles.max_events_per_sec;


%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;

%Calculate the threshold value to detect a lick
all_refs=[];
for trNo=firstTr:lastTr
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            %Note: handles.peakLFPNo is the reference LFP
            [referenceLFP, trialNo, can_read1] = drgGetTrialLFPData(handles, lickLFPNo, evNo, handles.evTypeNo, handles.startRef, handles.time_end);
            
            if (can_read1==1)
                all_refs=[all_refs referenceLFP];
            end
        end
    end
end

thershold_ref=prctile(all_refs,1)+((prctile(all_refs,99)-prctile(all_refs,1))/2);

%First find the time range for the spectrogram
if handles.subtractRef==1
    if handles.time_start+handles.time_pad<handles.startRef+handles.time_pad
        min_t=handles.time_start+handles.time_pad-(handles.window/2);
    else
        min_t=handles.startRef+handles.time_pad-(handles.window/2);
    end
    
    if handles.time_end-handles.time_pad>handles.endRef-handles.time_pad
        max_t=handles.time_end-handles.time_pad+handles.window;
    else
        max_t=handles.endRef-handles.time_pad+handles.window;
    end
else
    min_t=handles.time_start+handles.time_pad-(handles.window/2);
    max_t=handles.time_end-handles.time_pad+handles.window;
end


%Now get the LFP phase of the events

no_trials=0;

for trNo=firstTr:lastTr
    
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            %Note: handles.peakLFPNo is the reference LFP
            [ref, trialNo, can_read1] = drgGetTrialLFPData(handles, lickLFPNo, evNo, handles.evTypeNo, handles.time_start+handles.time_pad, handles.time_end-handles.time_pad);
            
            if (can_read1==1)
                
                %Get the events
                lick_timecourse_this_trial(1,1:no_lick_dt)=zeros(1,no_lick_dt);
                
                
                if ref(1)>thershold_ref
                    ii=find(ref<thershold_ref,1,'first');
                else
                    ii=1;
                end
                
                the_end=0;
                
                this_lick_ii=0;
                these_lick_times=[];
                while the_end==0
                    next_event=find(ref(ii:end)>thershold_ref,1,'first');
                    if isempty(next_event)
                        the_end=1;
                    else
                        
                        ii=ii+next_event-1;
                        
                        this_jj=floor((ii/handles.drg.session(sessionNo).draq_p.ActualRate)/handles.dt_lick)+1;
                        if this_jj<=no_lick_dt
                            
                            %Find the inter lick interval
                            this_lick_ii=this_lick_ii+1;
                            these_lick_times(this_lick_ii)=(ii/handles.drg.session(sessionNo).draq_p.ActualRate);
                            if this_lick_ii>1
                                %record the inter lick interval
                                ii_ili=ii_ili+1;
                                inter_lick_intervals(ii_ili)=these_lick_times(this_lick_ii)-these_lick_times(this_lick_ii-1);
                            end
                            
                            %Enter the lick in the timecourse only if it is
                            %not within a burst of high frequency noise
                            if ii_ili>0
                                if inter_lick_intervals(ii_ili)>handles.smallest_inter_lick_interval
                                    lick_timecourse_this_trial(1,this_jj)=lick_timecourse_this_trial(1,this_jj)+1;
                                end
                            end
                        end
                        
                        
                        end_event=find(ref(ii:end)<thershold_ref,1,'first');
                        if isempty(end_event)
                            the_end=1;
                        else
                            ii=ii+end_event-1;
                        end
                    end
                end
                
%                 for kk=1:no_lick_dt
%                    if lick_timecourse_this_trial(kk)>3
%                        lick_timecourse_this_trial(kk)=3;
%                    end
%                 end
                
                if sum(lick_timecourse_this_trial)<max_events_per_trial
                    
                    no_trials=no_trials+1;
                    lick_timecourse(no_trials,1:no_lick_dt)=lick_timecourse_this_trial(1,:);
                    per_corr_per_trial(no_trials)=perCorr(drgFindEvNo(handles,trialNo,sessionNo,odorOn));
                    
                    
                    if isfield(handles,'drgbchoices')
                        for evTypeNo=1:length(handles.drgbchoices.evTypeNos)
                            switch handles.evTypeNo
                                case 1
                                    %tstart is the reference event
                                    if handles.drgbchoices.evTypeNos(evTypeNo)==1
                                        %This is tstart
                                        if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(handles.drgbchoices.referenceEvent).times(evNo))>0
                                            which_event(evTypeNo,no_trials)=1;
                                        else
                                            which_event(evTypeNo,no_trials)=0;
                                        end
                                    else
                                        %These are not tstart, and the time
                                        %should be compared at OdorOn
                                        %This is tstart
                                        if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(2).times(evNo))>0
                                            which_event(evTypeNo,no_trials)=1;
                                        else
                                            which_event(evTypeNo,no_trials)=0;
                                        end
                                    end
                                otherwise
                                    %OdorOn is the reference event
                                    if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(handles.drgbchoices.referenceEvent).times(evNo))>0
                                        which_event(evTypeNo,no_trials)=1;
                                    else
                                        which_event(evTypeNo,no_trials)=0;
                                    end
                            end
                            
                            
                            
                        end
                    end
                end
                
                
                
                
                
            end
        end
    end
    
    %end %if eventstamps...
end %for evNo



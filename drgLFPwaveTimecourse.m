function drgLFPwaveTimecourse(handles)
%Generates a timecourse of the LFP power in decibels 10*log10(Power)


[t,freq,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo]=drgGetLFPwavePowerForThisEvTypeNo(handles);



%Timecourse doing average after log
%Get max and min
if handles.subtractRef==0
    log_P_timecourse=zeros(length(freq),length(t));
    log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),1);
    if handles.autoscale==1
        maxLogP=prctile(log_P_timecourse(:),99);
        minLogP=prctile(log_P_timecourse(:),1);
    else
        maxLogP=handles.maxLogP;
        minLogP=handles.minLogP;
    end
else
    log_P_timecourse=zeros(length(freq),length(t));
    log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),1);
    log_P_timecourse_ref=zeros(length(freq),length(t));
    log_P_timecourse_ref(:,:)=repmat(mean(10*log10(all_Power_ref),1)',1,length(t));
    if handles.autoscale==1
        maxLogP=prctile(log_P_timecourse(:)-log_P_timecourse_ref(:),99);
        minLogP=prctile(log_P_timecourse(:)-log_P_timecourse_ref(:),1);
    else
        maxLogP=handles.maxLogP;
        minLogP=handles.minLogP;
    end
end

%Note: Diego added this on purpose to limit the range to 10 dB
%This results in emphasizing changes in the top 10 dB
if maxLogP-minLogP>12
    minLogP=maxLogP-12;
end

try
    close 1
catch
end

%Plot the timecourse
hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.07 .1 .75 .3])

%     if handles.subtractRef==0
%         pcolor(repmat(t,length(freq),1)',repmat(freq,length(t),1),P_timecourse')
%     else
%         pcolor(repmat(t,length(freq),1)',repmat(freq,length(t),1),P_timecourse'-P_timecourse_ref')
%     end

if ~isempty(this_trialNo)
    
    if handles.subtractRef==0
        drg_pcolor(repmat(t,length(freq),1)',repmat(freq,length(t),1),log_P_timecourse')
    else
        %pcolor(repmat(t,length(f),1)',repmat(f,length(t),1),10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
        drg_pcolor(repmat(t,length(freq),1)',repmat(freq,length(t),1),log_P_timecourse'-log_P_timecourse_ref')
        %imagesc(t,f,10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
    end
    
    colormap jet
    shading interp
    caxis([minLogP maxLogP]);
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)');
    title(['Power (dB, wavelet) timecourse ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
    
    %If this is a single trial show the licks
    if length(this_trialNo)==1
        
        if handles.subtractRef==1
            if handles.time_start<handles.startRef
                min_t=handles.time_start;
            else
                min_t=handles.startRef;
            end
            
            if handles.time_end>handles.endRef
                max_t=handles.time_end;
            else
                max_t=handles.endRef;
            end
        else
            min_t=handles.time_start;
            max_t=handles.time_end;
        end
        
        
        %Calculate the threshold value to detect a lick
        sessionNo=1;
        all_refs=[];
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
                    [referenceLFP, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.startRef, handles.time_end);
                    
                    if (can_read1==1)
                        all_refs=[all_refs referenceLFP];
                    end
                end
            end
        end
        handles.evTypeNo=evtNo;
        
        thershold_ref=prctile(all_refs,1)+((prctile(all_refs,99)-prctile(all_refs,1))/2);
        
        %Get the licks
        
        evNo = drgFindEvNo(handles,handles.trialNo,sessionNo);
        
        [referenceLFP, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, min_t, max_t);
        handles.peakLFPNo=hpNo;
        
        ii_start_ref=handles.time_pad*handles.drg.session(sessionNo).draq_p.ActualRate;
        ii_end_ref=length(referenceLFP)-ii_start_ref;
        refLFP_ref=referenceLFP(ii_start_ref:ii_end_ref);
        
        inter_lick_intervals_ref=[];
        ii_ili_ref=0;
        
        if refLFP_ref(1)>thershold_ref
            ii=find(refLFP_ref<thershold_ref,1,'first');
        else
            ii=1;
        end
        
        the_end=0;
        ref_power_these_events=[];
        no_ref_evs_this_trial=0;
        
        this_lick_ii=0;
        stamped_lick_ii=0;
        these_lick_times=[];
        these_stamped_lick_times=[];
        
        %Find the events (licks)
        while the_end==0
            next_event=find(refLFP_ref(ii:end)>thershold_ref,1,'first');
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
                    ii_ili_ref=ii_ili_ref+1;
                    inter_lick_intervals_ref(ii_ili_ref)=these_lick_times(this_lick_ii)-these_lick_times(this_lick_ii-1);
                end
                
                %Enter the event (lick) in the timecourse only if it is
                %not within a burst of high frequency noise
                if ii_ili_ref>0
                    if inter_lick_intervals_ref(ii_ili_ref)>handles.smallest_inter_lick_interval
                        %stamp the lick
                        stamped_lick_ii=stamped_lick_ii+1;
                        these_stamped_lick_times(stamped_lick_ii)=(ii/handles.drg.session(sessionNo).draq_p.ActualRate);
                    end
                end
                
                end_event=find(refLFP_ref(ii:end)<thershold_ref,1,'first');
                if isempty(end_event)
                    the_end=1;
                else
                    ii=ii+end_event-1;
                end
            end
        end
        
        %PLot the licks
        hold on
        
        %Clear the top strip
        ffrom=freq(1)+0.95*(freq(end)-freq(1));
        
        
        for t1=t(1):t(2)-t(1):t(end)
            plot([t1 t1],[ffrom freq(end)],'-w','LineWidth',3)
        end
        
        
        for ilick=1:stamped_lick_ii
            plot([these_stamped_lick_times(ilick)+t(1) these_stamped_lick_times(ilick)+t(1)],[ffrom freq(end)],'-k','LineWidth',3)
        end
    end
    
    try
        close 2
    catch
    end
    
    hFig2 = figure(2);
    set(hFig2, 'units','normalized','position',[.83 .1 .05 .3])
    
    prain=[minLogP:(maxLogP-minLogP)/99:maxLogP];
    drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    colormap jet
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')
end


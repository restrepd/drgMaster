function [log_P_t,which_event,freq,out_times,perCorrERP,no_events_per_trial,trials_with_wavelick,f_lick,mean_lick_phase]=drgLFP_ERWA(handles)
%Generates a event related wavelet power for licks

odorOn=2;
sessionNo=handles.sessionNo;

%ERWA timing
wing_dt=0.2;
Fs=handles.drg.session(1).draq_p.ActualRate;
dec_n=fix(Fs/1000);
decFs=Fs/dec_n;
dt=1/decFs;
wing_ii=fix(wing_dt/dt);
out_times=double([-wing_ii:wing_ii])*dt;

%These are the inputs from the pad
time_start=handles.time_start;
time_end=handles.time_end;
startRef=handles.startRef;
endRef=handles.endRef;

%Add a wing to the start and end times
if (time_start-wing_dt)>(startRef-wing_dt)
    handles.time_start=startRef-wing_dt;
else
    handles.time_start=time_start-wing_dt;
end
if (time_end+wing_dt)>(endRef+wing_dt)
    handles.time_end=time_end+wing_dt;
else
    handles.time_end=endRef+wing_dt;
end

[t_apt,freq,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo]=drgGetLFPwavePowerForThisEvTypeNo(handles);

%Calculate the licks (no wing nescessary)
[lick_freq,times_lick_freq,lick_traces,CIlickf,lick_trace_times,stamped_lick_ii,these_stamped_lick_times,no_trials,trials_included,lick_threshold]=drgGetLicks(handles);
 
%Get PAC
handles=drgThetaAmpPhaseTrialRange(handles);
perCorr=handles.drgb.PAC.perCorr;

%Find ERWA for licks
wing_dt=0.2;
dt=t_apt(2)-t_apt(1);
wing_ii=fix(wing_dt/dt);
out_times=double([-wing_ii:wing_ii])*dt;
envelope_times=handles.drgb.PAC.out_times_env+handles.time_start+handles.time_pad;

total_no_licks=0;
total_trNo=0;
for trNo=1:length(this_trialNo)
    this_lick_trNo=find(this_trialNo(trNo)==trials_included);
    if exist('this_lick_trNo')~=0
        %Calculate lick-referenced power within the measurement window
        licks_for_this_trial=[];
        licks_for_this_trial=these_stamped_lick_times(this_lick_trNo,1:stamped_lick_ii(this_lick_trNo));
        no_licks=sum((licks_for_this_trial>time_start+handles.time_pad)&(licks_for_this_trial<=time_end-handles.time_pad));
        if no_licks>0
            total_no_licks=total_no_licks+no_licks;
            total_trNo=total_trNo+1;
        end
    end
end


no_trs=0;
log_P_t=zeros(total_trNo,length(freq),length(out_times));
perCorrERP=zeros(1,total_trNo);
no_events_per_trial=zeros(1,total_trNo);
trials_with_wavelick=zeros(1,total_trNo);
mean_lick_phase=zeros(1,total_trNo);
f_lick=zeros(1,total_trNo);
first_which=0;
%log_P_timecourse=zeros(total_trNo,length(freq),length(t_apt));

for trNo=1:length(this_trialNo)
    
    this_lick_trNo=find(this_trialNo(trNo)==trials_included);
    if exist('this_lick_trNo')~=0
        
        %Calculate lick-referenced power within the measurement window
        licks_for_this_trial=[];
        licks_for_this_trial=these_stamped_lick_times(this_lick_trNo,1:stamped_lick_ii(this_lick_trNo));
        no_licks=sum((licks_for_this_trial>time_start+handles.time_pad)&(licks_for_this_trial<=time_end-handles.time_pad));
        these_angles=zeros(1,no_licks);
        if no_licks>0
            %Calculate the mean phase for licks
            jj_l=0;
            for kk_licks=1:length(licks_for_this_trial)
                if (licks_for_this_trial(kk_licks)>time_start+handles.time_pad)&(licks_for_this_trial(kk_licks)<=time_end-handles.time_pad)
                   [mint this_ii_lick]=min(abs(envelope_times-licks_for_this_trial(kk_licks)));
                   jj_l=jj_l+1;
                   these_angles(jj_l)=handles.drgb.PAC.PACtimecourse(trNo).decanglethetaLFP(this_ii_lick);
                end
            end
            mean_lick_phase(trNo)=(180/pi)*circ_axial(circ_mean(these_angles')');
            
            %Calculate log power
            if handles.subtractRef==0
                this_log_P_timecourse=zeros(length(freq),length(t_apt));
                this_log_P_timecourse(:,:)=10*log10(all_Power_timecourse(trNo,:,:));
            else
                this_log_P_timecourse=zeros(length(freq),length(t_apt));
                this_log_P_timecourse(:,:)=10*log10(all_Power_timecourse(trNo,:,:));
                this_log_P_timecourse_ref=zeros(length(freq),length(t_apt));
                this_log_P_timecourse_ref(:,:)=repmat(mean(10*log10(all_Power_ref(trNo,:)),1)',1,length(t_apt));
                this_log_P_timecourse=this_log_P_timecourse-this_log_P_timecourse_ref;
            end
            no_trs=no_trs+1;
            no_events_per_trial(no_trs)=no_licks;
            f_lick(no_trs)=no_licks/(time_end-time_start-2*handles.time_pad);
            %log_P_timecourse(no_trs,:,:)=this_log_P_timecourse;
            
            these_log_P_ts=zeros(no_licks,length(freq),length(out_times));
            jj_lick=0;
            for ii_lick=1:length(licks_for_this_trial)
                if (licks_for_this_trial(ii_lick)>time_start+handles.time_pad)&(licks_for_this_trial(ii_lick)<=time_end-handles.time_pad)
                    [min_dt,min_ii]=min(abs(t_apt-licks_for_this_trial(ii_lick)));
                    jj_lick=jj_lick+1;
                    these_log_P_ts(jj_lick,:,:)=this_log_P_timecourse(:,min_ii(1)-wing_ii:min_ii(1)+wing_ii);
                end
            end
            log_P_t(no_trs,:,:)=mean(these_log_P_ts,1);
            
            perCorrERP(no_trs)=perCorr(trNo);
            trials_with_wavelick(no_trs)=trials_included(trNo);
            
            if isfield(handles,'drgbchoices')
                if first_which==0
                    which_event=zeros(length(handles.drgbchoices.evTypeNos),total_trNo);
                    first_which=1;
                end
                for evTypeNo=1:length(handles.drgbchoices.evTypeNos)
                    switch handles.evTypeNo
                        case 1
                            %tstart is the reference event
                            if handles.drgbchoices.evTypeNos(evTypeNo)==1
                                %This is tstart
                                if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(handles.drgbchoices.referenceEvent).times(evNo))>0
                                    which_event(evTypeNo,no_trs)=1;
                                else
                                    which_event(evTypeNo,no_trs)=0;
                                end
                            else
                                %These are not tstart, and the time
                                %should be compared at OdorOn
                                %This is tstart
                                if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(2).times(evNo))>0
                                    which_event(evTypeNo,no_trs)=1;
                                else
                                    which_event(evTypeNo,no_trs)=0;
                                end
                            end
                        otherwise
                            %OdorOn is the reference event
                            if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(handles.drgbchoices.referenceEvent).times(evNo))>0
                                which_event(evTypeNo,no_trs)=1;
                            else
                                which_event(evTypeNo,no_trs)=0;
                            end
                    end
                    
                end
            end
        end
    end
end


no_trials=length(this_trialNo);

handles.time_start=time_start;
handles.time_end=time_end;
handles.startRef=startRef;
handles.endRef=endRef;

if handles.displayData==1
    if ~isempty(this_trialNo)
        
        
        
        %Lick-related spectrogram
        %Timecourse doing average after log
        try
            close 15
        catch
        end
        
        %Plot the ERP timecourse
        hFig1 = figure(15);
        set(hFig1, 'units','normalized','position',[.04 .1 .55 .3])
        
        %Calculate the mean ERP timecourse
        mean_log_P_t=zeros(length(freq),length(out_times));
        mean_log_P_t(:,:)=mean(log_P_t,1);
        
        %Get max and min
        if handles.autoscale==1
            maxLogP=prctile(mean_log_P_t(:),99);
            minLogP=prctile(mean_log_P_t(:),1);
            
            %Note: Diego added this on purpose to limit the range to 10 dB
            %This results in emphasizing changes in the top 10 dB
            if maxLogP-minLogP>10
                minLogP=maxLogP-10;
            end
        else
            maxLogP=handles.maxLogP;
            minLogP=handles.minLogP;
        end
        
        
        drg_pcolor(repmat(out_times-mean(out_times),length(freq),1)',repmat(freq,length(out_times),1),mean_log_P_t')
        
        
        colormap jet
        shading interp
        caxis([minLogP maxLogP]);
        xlabel('Time (sec)')
        ylabel('Frequency (Hz)');
        title(['Lick-related power wavelet spectrogram (dB) for ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
        
        try
            close 16
        catch
        end
        
        hFig2 = figure(16);
        set(hFig2, 'units','normalized','position',[.6 .1 .05 .3])
        
        prain=[minLogP:(maxLogP-minLogP)/99:maxLogP];
        drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
        colormap jet
        shading interp
        ax=gca;
        set(ax,'XTickLabel','')
        ylabel('dB')
        
        try
            close 17
        catch
        end
        
        hFig2 = figure(17);
        set(hFig2, 'units','normalized','position',[.67 .1 .3 .3])
        polarhistogram(pi*mean_lick_phase/180,6)
        title(['Mean angle for per trial phase of the licks  for ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
        
    end
end
pfft=1;



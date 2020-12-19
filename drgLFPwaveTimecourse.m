function handles=drgLFPwaveTimecourse(handles)
%Generates a timecourse of the LFP power in decibels 10*log10(Power)

handles.drgb.PACwave=[];

[t_apt,freq,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo]=drgGetLFPwavePowerForThisEvTypeNo(handles);

handles.drgb.PACwave.trialNos_PRP=this_trialNo;

%Calculate the licks
[lick_freq,times_lick_freq,lick_traces,CIlickf,lick_trace_times,stamped_lick_ii,these_stamped_lick_times,no_trials,...
    trials_included,lick_threshold, lick_freq_per_trial,trials_included_per_trial]=drgGetLicks(handles);

handles.drgb.PACwave.times_lick_freq=times_lick_freq;
handles.drgb.PACwave.lick_trials_included=trials_included;
handles.drgb.PACwave.lick_freq_per_trial=lick_freq_per_trial;
handles.drgb.PACwave.lick_freq_per_trial=lick_freq_per_trial;
handles.drgb.PACwave.stamped_lick_ii=stamped_lick_ii;
handles.drgb.PACwave.stamped_lick_ii=stamped_lick_ii;
handles.drgb.PACwave.these_stamped_lick_times=these_stamped_lick_times;

%Get PAC
handles=drgThetaAmpPhaseTrialRange(handles);

handles.drgb.PACwave.trialNos_PAC=handles.drgb.PAC.this_trialNo;

%If there was a problem with signal saturation for this LFP the number of
%trials is zero, and we should skip further analysis
if handles.drgb.PAC.no_trials>0
    %Please note that more trials are excluded from the PAC analysis than from
    %the lick analysis
    
    %Find the wavelet power at the peak and trough
    dt=handles.window-handles.noverlap;
    t_pac=[handles.time_start+handles.time_pad:dt:handles.time_end-handles.time_pad];
    handles.drgb.PACwave.t_pac=t_pac;
    peakPower=zeros(handles.drgb.PAC.no_trials,length(t_pac));
    %Find the value of gamma power at each point
    out_times_env=handles.drgb.PAC.out_times_env;
    out_times_env=out_times_env+handles.time_start+handles.time_pad;
    meanPeakAngle=(handles.drgb.PAC.peakAngleForPower*(pi/180))-pi;
    meanTroughAngle=(handles.drgb.PAC.troughAngleForPower*(pi/180))-pi;
    
    
    %Find peak wavelet power
    peakPowerSpectrum=[];
    troughPowerSpectrum=[];
    
    for trNum=1:handles.drgb.PAC.no_trials
        this_peakPower=[];
        this_peakPower_spectrum=[];
        this_peakPower_times=[];
        ii=1;
        jj=0;
        at_end=0;
        this_angleTetaLFP=handles.drgb.PAC.PACtimecourse(trNum).decanglethetaLFP;
        this_LFPenv=handles.drgb.PAC.PACtimecourse(trNum).decLFPgenv;
        while at_end==0
            ii_next=find(this_angleTetaLFP(ii:end)>=meanPeakAngle,1,'first');
            if (~isempty(ii_next))&(ii+ii_next-1<=length(t_apt))&(ii+ii_next-1<=length(out_times_env))
                jj=jj+1;
                this_peakPower(jj)=mean(10*log10(all_Power_timecourse(trNum,:,ii+ii_next-1)),2);
                this_peakPower_spectrum(jj,:)=10*log10(all_Power_timecourse(trNum,:,ii+ii_next-1));
                this_peakPower_times(jj)=out_times_env(ii+ii_next-1);
                ii=ii+ii_next;
                ii_next=find(this_angleTetaLFP(ii:end)<meanPeakAngle,1,'first');
                if ~isempty(ii_next)
                    ii=ii+ii_next;
                else
                    at_end=1;
                end
                
            else
                at_end=1;
            end
        end
        
        d_t_pac=t_pac(2)-t_pac(1);
        
        %         for ii_t=1:length(t_pac)
        %             %If this point is before any peakPower times enter the first
        %             %peakPower
        %             if t_pac(ii_t)+(d_t_pac/2)<=this_peakPower_times(1)
        %                 peakPower(trNum,ii_t)=this_peakPower(1);
        %                 %this_peakPower_t_pac(ii_t)=this_peakPower(1);
        %             else
        %                 %If this point is after any peakPower times enter the last
        %                 %peakPower
        %                 if t_pac(ii_t)-(d_t_pac/2)>=this_peakPower_times(end)
        %                     peakPower(trNum,ii_t)=this_peakPower(end);
        %                     %this_peakPower_t_pac(ii_t)=this_peakPower(end);
        %                 else
        %                     %Now make the power value equal to the mean within the
        %                     %bin
        %                     ii_pt=find((this_peakPower_times>=t_pac(ii_t)-d_t_pac/2)&(this_peakPower_times<=t_pac(ii_t)+d_t_pac/2));
        %                     if ~isempty(ii_pt)
        %                         peakPower(trNum,ii_t)=mean(this_peakPower(ii_pt));
        %                     else
        %                         peakPower(trNum,ii_t-1);
        %                     end
        %                 end
        %             end
        %         end
        %
        
        
        for ii_t=1:length(t_pac)
            if t_pac(ii_t)<=this_peakPower_times(1)
                peakPower(trNum,ii_t)=this_peakPower(1);
                %this_peakPower_t_pac(ii_t)=this_peakPower(1);
            else
                if t_pac(ii_t)>=this_peakPower_times(end)
                    peakPower(trNum,ii_t)=this_peakPower(end);
                    %this_peakPower_t_pac(ii_t)=this_peakPower(end);
                else
                    ii_pt=find(this_peakPower_times>=t_pac(ii_t),1,'first');
                    peakPower(trNum,ii_t)=this_peakPower(ii_pt-1)+(t_pac(ii_t)-this_peakPower_times(ii_pt-1))*((this_peakPower(ii_pt)-this_peakPower(ii_pt-1))/(this_peakPower_times(ii_pt)-this_peakPower_times(ii_pt-1)));
                end
            end
        end
        
        if handles.subtractRef==1
            peakPower(trNum,:)=peakPower(trNum,:)-mean(peakPower(trNum,(t_pac>=handles.startRef+handles.time_pad)&(t_pac<=handles.endRef-handles.time_pad)));
        end
        
        handles.drgb.PACwave.peakPowerSpectrum(trNum,:)=mean(this_peakPower_spectrum,1);
        handles.drgb.PACwave.PACtimecourse(trNum).peakPower=peakPower(trNum,:);
        handles.drgb.PACwave.meanPeakPower(trNum)=mean(peakPower(trNum,:));
    end
    
    %Find trough power
    troughPower=zeros(handles.drgb.PAC.no_trials,length(t_pac));
    for trNum=1:handles.drgb.PAC.no_trials
        this_troughPower=[];
        this_troughPower_spectrum=[];
        this_troughPower_times=[];
        ii=1;
        jj=0;
        at_end=0;
        this_angleTetaLFP=handles.drgb.PAC.PACtimecourse(trNum).decanglethetaLFP;
        this_LFPenv=handles.drgb.PAC.PACtimecourse(trNum).decLFPgenv;
        while at_end==0
            ii_next=find(this_angleTetaLFP(ii:end)>=meanTroughAngle,1,'first');
            if (~isempty(ii_next))&(ii+ii_next-1<=length(t_apt))&(ii+ii_next-1<=length(out_times_env))
                jj=jj+1;
                this_troughPower(jj)=mean(10*log10(all_Power_timecourse(trNum,:,ii+ii_next-1)),2);
                this_troughPower_spectrum(jj,:)=10*log10(all_Power_timecourse(trNum,:,ii+ii_next-1));
                this_troughPower_times(jj)=out_times_env(ii+ii_next-1);
                ii=ii+ii_next;
                ii_next=find(this_angleTetaLFP(ii:end)<meanTroughAngle,1,'first');
                if ~isempty(ii_next)
                    ii=ii+ii_next;
                else
                    at_end=1;
                end
            else
                at_end=1;
            end
        end
        
        d_t_pac=t_pac(2)-t_pac(1);
        
%         for ii_t=1:length(t_pac)
%             %If this point is before any troughPower times enter the first
%             %troughPower
%             if t_pac(ii_t)+(d_t_pac/2)<=this_troughPower_times(1)
%                 troughPower(trNum,ii_t)=this_troughPower(1);
%                 %this_troughPower_t_pac(ii_t)=this_troughPower(1);
%             else
%                 %If this point is after any troughPower times enter the last
%                 %troughPower
%                 if t_pac(ii_t)-(d_t_pac/2)>=this_troughPower_times(end)
%                     troughPower(trNum,ii_t)=this_troughPower(end);
%                     %this_troughPower_t_pac(ii_t)=this_troughPower(end);
%                 else
%                     %Now make the power value equal to the mean within the
%                     %bin
%                     ii_pt=find((this_troughPower_times>=t_pac(ii_t)-d_t_pac/2)&(this_troughPower_times<=t_pac(ii_t)+d_t_pac/2));
%                     if ~isempty(ii_pt)
%                         troughPower(trNum,ii_t)=mean(this_troughPower(ii_pt));
%                     else
%                         troughPower(trNum,ii_t-1);
%                     end
%                 end
%             end
%         end
        
                for ii_t=1:length(t_pac)
                    if t_pac(ii_t)<=this_troughPower_times(1)
                        troughPower(trNum,ii_t)=this_troughPower(1);
                        %this_troughPower_t_pac(ii_t)=this_troughPower(1);
                    else
                        if t_pac(ii_t)>=this_troughPower_times(end)
                            troughPower(trNum,ii_t)=this_troughPower(end);
                            %this_troughPower_t_pac(ii_t)=this_troughPower(end);
                        else
                            ii_pt=find(this_troughPower_times>=t_pac(ii_t),1,'first');
                            troughPower(trNum,ii_t)=this_troughPower(ii_pt-1)+(t_pac(ii_t)-this_troughPower_times(ii_pt-1))*((this_troughPower(ii_pt)-this_troughPower(ii_pt-1))/(this_troughPower_times(ii_pt)-this_troughPower_times(ii_pt-1)));
                        end
                    end
                end
        
        if handles.subtractRef==1
            troughPower(trNum,:)=troughPower(trNum,:)-mean(troughPower(trNum,(t_pac>=handles.startRef+handles.time_pad)&(t_pac<=handles.endRef-handles.time_pad)));
        end
        
        handles.drgb.PACwave.troughPowerSpectrum(trNum,:)=mean(this_troughPower_spectrum,1);
        handles.drgb.PACwave.PACtimecourse(trNum).troughPower=troughPower(trNum,:);
        handles.drgb.PACwave.meanTroughPower(trNum)=mean(troughPower(trNum,:));
    end
    
    %     %Calculate mean power
    %     for trNum=1:handles.drgb.PAC.no_trials
    %
    %         for ii_t=1:length(t_pac)
    %             if ii_t==1
    %                 meanPower(trNum,ii_t)=mean(10*log10(all_Power_timecourse(trNum,:,1)),2);
    %             else
    %                 if ii_t==length(t_pac)
    %                     meanPower(trNum,ii_t)=mean(10*log10(all_Power_timecourse(trNum,:,end)),2);
    %                 else
    %                     meanPower(trNum,ii_t)=mean(mean(10*log10(all_Power_timecourse(trNum,:,(t_apt>=t_pac(ii_t)-(dt/2))&(t_apt<t_pac(ii_t)+(dt/2)))),2),3);
    %                 end
    %             end
    %         end
    %         if handles.subtractRef==1
    %             meanPower(trNum,:)=meanPower(trNum,:)-mean(meanPower(trNum,(t_pac>=handles.startRef+handles.time_pad)&(t_pac<=handles.endRef-handles.time_pad)));
    %         end
    %
    %         handles.drgb.PACwave.PACtimecourse(trNum).meanPower=meanPower(trNum,:);
    %         handles.drgb.PACwave.meanPower(trNum)=mean(meanPower(trNum,:));
    %     end
    
    %Place lick per trial
    ntrs_out=0;
    n_lpw_out=0;
    %Let's space the lick-triggered power by 70 msec
    delta_ii=floor(0.07/(t_apt(2)-t_apt(1)));
    ii_span=4;
    
    d_t_pac=t_pac(2)-t_pac(1);
    
    for trNum=1:length(trials_included_per_trial)
        ii_tr=find(handles.drgb.PAC.this_trialNo==trials_included_per_trial(trNum),1);
        if ~isempty(ii_tr)
            ntrs_out=ntrs_out+1;
            if handles.subtractRef==1
                lick_freq_per_trial(trNum,:)=lick_freq_per_trial(trNum,:)-mean(lick_freq_per_trial(trNum,(t_pac>=handles.startRef+handles.time_pad)&(t_pac<=handles.endRef-handles.time_pad)));
            end
            
            handles.drgb.PACwave.trNo_lick_and_PRP(ntrs_out)=trials_included_per_trial(trNum);
            handles.drgb.PACwave.no_lick_trials=ntrs_out;
            handles.drgb.PACwave.lick_timecourse(ntrs_out).lick_f=lick_freq_per_trial(trNum,:);
            handles.drgb.PACwave.lick_timecourse(ntrs_out).times_lick_freq=times_lick_freq;
            handles.drgb.PACwave.mean_lick_freq(ntrs_out)=mean(lick_freq_per_trial(trNum,:));
            
            %Now enter peakPower and troughPower for trials that include licks
            handles.drgb.PACwave.meanPeakPower_per_lick_trial(ntrs_out)=handles.drgb.PACwave.meanPeakPower(ii_tr);
            handles.drgb.PACwave.meanTroughPower_per_lick_trial(ntrs_out)=handles.drgb.PACwave.meanTroughPower(ii_tr);
            
            handles.drgb.PACwave.lick_triggered_wpower_no(ntrs_out)=0;
            this_power=zeros(1,size(all_Power_timecourse,2));
            this_lickPower_times=[];
            this_lickPower=[];
            for jj=1:stamped_lick_ii(trNum)
                
                %Find the t_apt
                [min_t_apt,ii_min_t_apt]=min(abs(t_apt-these_stamped_lick_times(trNum,jj)));
                
                %Log in power stamp if the data are in the array
                if ((ii_min_t_apt-delta_ii*ii_span)>0)&((ii_min_t_apt+delta_ii*ii_span)<length(t_apt))
                    handles.drgb.PACwave.lick_triggered_wpower_no(ntrs_out)= handles.drgb.PACwave.lick_triggered_wpower_no(ntrs_out)+1;
                    handles.drgb.PACwave.lick_triggered_wpower_t(ntrs_out,handles.drgb.PACwave.lick_triggered_wpower_no(ntrs_out))=these_stamped_lick_times(trNum,jj);
                    for iis=1:2*ii_span+1
                        this_power=zeros(1,size(all_Power_timecourse,2));
                        this_power(1,:)=all_Power_timecourse(ii_tr,:,ii_min_t_apt-delta_ii*ii_span+delta_ii*(iis-1));
                        handles.drgb.PACwave.lick_triggered_wpower(ntrs_out,handles.drgb.PACwave.lick_triggered_wpower_no(ntrs_out),iis,1:size(all_Power_timecourse,2))=...
                            this_power(1,:);
                        if iis==ii_span+1
                            this_lickPower(handles.drgb.PACwave.lick_triggered_wpower_no(ntrs_out))=mean(10*log10(this_power(1,:)));
                        end
                    end
                    this_lickPower_times(handles.drgb.PACwave.lick_triggered_wpower_no(ntrs_out))=handles.drgb.PACwave.lick_triggered_wpower_t(ntrs_out,handles.drgb.PACwave.lick_triggered_wpower_no(ntrs_out));
                    
                end
                
            end
            
            %Now enter the lick powers
            
            %Note: Lick power is only computed if there are lick times
            if ~isempty(this_lickPower_times)
                n_lpw_out=n_lpw_out+1;
                handles.drgb.PACwave.no_lickPower_trials=n_lpw_out;
                handles.drgb.PACwave.lickPower_trials(n_lpw_out)=trials_included_per_trial(trNum);
%                 for ii_t=1:length(t_pac)
%                     %If this point is before any lickPower times enter the first
%                     %lickPower
%                     if t_pac(ii_t)+(d_t_pac/2)<=this_lickPower_times(1)
%                         lickPower(n_lpw_out,ii_t)=this_lickPower(1);
%                         %this_lickPower_t_pac(ii_t)=this_lickPower(1);
%                     else
%                         %If this point is after any lickPower times enter the last
%                         %lickPower
%                         if t_pac(ii_t)-(d_t_pac/2)>=this_lickPower_times(end)
%                             lickPower(n_lpw_out,ii_t)=this_lickPower(end);
%                             %this_lickPower_t_pac(ii_t)=this_lickPower(end);
%                         else
%                             %Now make the power value equal to the mean within the
%                             %bin
%                             ii_pt=find((this_lickPower_times>=t_pac(ii_t)-d_t_pac/2)&(this_lickPower_times<=t_pac(ii_t)+d_t_pac/2));
%                             if ~isempty(ii_pt)
%                                 lickPower(n_lpw_out,ii_t)=mean(this_lickPower(ii_pt));
%                             else
%                                 lickPower(n_lpw_out,ii_t)=lickPower(n_lpw_out,ii_t-1);
%                             end
%                         end
%                     end
%                 end
                
                
                for ii_t=1:length(t_pac)
                    if t_pac(ii_t)<=this_lickPower_times(1)
                        lickPower(n_lpw_out,ii_t)=this_lickPower(1);
                        %this_troughPower_t_pac(ii_t)=this_troughPower(1);
                    else
                        if t_pac(ii_t)>=this_troughPower_times(end)
                            lickPower(n_lpw_out,ii_t)=this_lickPower(end);
                            %this_troughPower_t_pac(ii_t)=this_troughPower(end);
                        else
                            ii_pt=find(this_lickPower_times>=t_pac(ii_t),1,'first');
                            if ~isempty(ii_pt)
                                lickPower(n_lpw_out,ii_t)=this_lickPower(ii_pt-1)+(t_pac(ii_t)-this_lickPower_times(ii_pt-1))*((this_lickPower(ii_pt)-this_lickPower(ii_pt-1))/(this_lickPower_times(ii_pt)-this_lickPower_times(ii_pt-1)));
                            else
                                lickPower(n_lpw_out,ii_t)=this_lickPower(end);
                            end
                        end
                    end
                end
                
%                          for ii_t=1:length(t_pac)
%                     if t_pac(ii_t)<=this_troughPower_times(1)
%                         troughPower(trNum,ii_t)=this_troughPower(1);
%                         %this_troughPower_t_pac(ii_t)=this_troughPower(1);
%                     else
%                         if t_pac(ii_t)>=this_troughPower_times(end)
%                             troughPower(trNum,ii_t)=this_troughPower(end);
%                             %this_troughPower_t_pac(ii_t)=this_troughPower(end);
%                         else
%                             ii_pt=find(this_troughPower_times>=t_pac(ii_t),1,'first');
%                             troughPower(trNum,ii_t)=this_troughPower(ii_pt-1)+(t_pac(ii_t)-this_troughPower_times(ii_pt-1))*((this_troughPower(ii_pt)-this_troughPower(ii_pt-1))/(this_troughPower_times(ii_pt)-this_troughPower_times(ii_pt-1)));
%                         end
%                     end
%                 end
                
                if handles.subtractRef==1
                    lickPower(n_lpw_out,:)=lickPower(n_lpw_out,:)-mean(lickPower(n_lpw_out,(t_pac>=handles.startRef+handles.time_pad)&(t_pac<=handles.endRef-handles.time_pad)));
                end
                
                handles.drgb.PACwave.PACtimecourse(n_lpw_out).lickPower=lickPower(n_lpw_out,:);
                handles.drgb.PACwave.meanLickPower(n_lpw_out)=mean(lickPower(n_lpw_out,:));
            end
        end
    end
    
    
    
    
    
    
    
    if handles.displayData==1
        if ~isempty(this_trialNo)
            
            %Timecourse doing average after log
            %Get max and min
            if handles.subtractRef==0
                log_P_timecourse=zeros(length(freq),length(t_apt));
                log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),1);
                
                %Per trial power plot
                log_P_per_trial_timecourse=zeros(length(freq)*length(this_trialNo),length(t_apt));
                y_shift=0;
                for trialNo=1:length(this_trialNo)
                    this_log_P_timecourse=zeros(length(freq),length(t_apt));
                    this_log_P_timecourse(:,:)=10*log10(all_Power_timecourse(trialNo,:,:));
                    log_P_per_trial_timecourse(y_shift+1:y_shift+length(freq),:)=this_log_P_timecourse;
                    shifted_freq(1,y_shift+1:y_shift+length(freq))=freq+(trialNo-1)*freq(end);
                    y_shift=y_shift+length(freq);
                end
                
                if handles.autoscale==1
                    maxLogPper=prctile(log_P_timecourse(:),99);
                    minLogPper=prctile(log_P_timecourse(:),1);
                    %Note: Diego added this on purpose to limit the range to 10 dB
                    %This results in emphasizing changes in the top 10 dB
                    if maxLogPper-minLogPper>12
                        minLogPper=maxLogPper-12;
                    end
                else
                    maxLogPper=handles.maxLogP;
                    minLogPper=handles.minLogP;
                end
            else
                log_P_timecourse=zeros(length(freq),length(t_apt));
                log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),1);
                log_P_timecourse_ref=zeros(length(freq),length(t_apt));
                log_P_timecourse_ref(:,:)=repmat(mean(10*log10(all_Power_ref),1)',1,length(t_apt));
                
                %Per trial power plot
                log_P_per_trial_timecourse_sub=zeros(length(freq)*length(this_trialNo),length(t_apt));
                y_shift=0;
                sy_shift=0;
                shifted_freq=[];
                for trialNo=1:length(this_trialNo)
                    this_log_P_timecourse=zeros(length(freq),length(t_apt));
                    this_log_P_timecourse(:,:)=10*log10(all_Power_timecourse(trialNo,:,:));
                    this_log_P_timecourse_ref=zeros(length(freq),length(t_apt));
                    this_log_P_timecourse_ref(:,:)=repmat(mean(10*log10(all_Power_ref(trialNo,:)),1)',1,length(t_apt));
                    log_P_per_trial_timecourse_sub(y_shift+1:y_shift+length(freq),:)=this_log_P_timecourse-this_log_P_timecourse_ref;
                    shifted_freq(1,y_shift+1:y_shift+length(freq))=freq+(trialNo-1)*freq(end);
                    y_shift=y_shift+length(freq);
                end
                
                max_delta=16;
                if handles.autoscale==1
                    
                    deltaLogP=log_P_timecourse'-log_P_timecourse_ref';
                    maxLogPper=prctile(deltaLogP(:),99);
                    minLogPper=prctile(deltaLogP(:),1);
                    %Note: Diego added this on purpose to limit the range to 10 dB
                    %This results in emphasizing changes in the top 10 dB
                    if maxLogPper-minLogPper>max_delta
                        minLogPper=maxLogPper-max_delta;
                    end
                    
                else
                    maxLogPper=handles.maxLogP;
                    minLogPper=handles.minLogP;
                end
            end
            
            figNo=8;
            
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            %Plot the timecourse
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.07 .05 .75 .3])
            if handles.subtractRef==0
                drg_pcolor(repmat(t_apt,length(freq),1)',repmat(freq,length(t_apt),1),log_P_timecourse')
            else
                %pcolor(repmat(t,length(f),1)',repmat(f,length(t),1),10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
                drg_pcolor(repmat(t_apt,length(freq),1)',repmat(freq,length(t_apt),1),log_P_timecourse'-log_P_timecourse_ref')
                %imagesc(t,f,10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
            end
            
            colormap fire
            shading interp
            caxis([minLogPper maxLogPper]);
            xlabel('Time (sec)')
            ylabel('Frequency (Hz)');
            title(['Power (dB, wavelet) timecourse ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
            
            %If this is a single trial show the licks
            if length(this_trialNo)==1
                
                
                %PLot the licks
                hold on
                
                %Clear the top strip
                ffrom=freq(1)+0.95*(freq(end)-freq(1));
                
                
                for t1=t_apt(1):t_apt(2)-t_apt(1):t_apt(end)
                    plot([t1 t1],[ffrom freq(end)],'-w','LineWidth',3)
                end
                
                
                for ilick=1:stamped_lick_ii
                    plot([these_stamped_lick_times(ilick) these_stamped_lick_times(ilick)],[ffrom freq(end)],'-k','LineWidth',3)
                end
            end
            
            %Plot the spectrum
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            %Plot the timecourse
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.1 .1 .4 .4])
            hold on
            
            
            this_mean_dbWB=zeros(1,length(freq));
            this_mean_dbWB(1,:)=mean(handles.drgb.PACwave.troughPowerSpectrum,1);
            
            try
                CI=[];
                CI = bootci(1000, {@mean, handles.drgb.PACwave.troughPowerSpectrum})';
                CI(:,1)= this_mean_dbWB'-CI(:,1);
                CI(:,2)=CI(:,2)- this_mean_dbWB';
                
                
                [hlCR, hpCR] = boundedline(freq',this_mean_dbWB', CI, 'b');
            catch
                plot(freq',this_mean_dbWB', 'b');
            end
            
            this_mean_dbWB=zeros(1,length(freq));
            this_mean_dbWB(1,:)=mean(handles.drgb.PACwave.peakPowerSpectrum,1);
            
            CI=[];
            CI = bootci(1000, {@mean, handles.drgb.PACwave.peakPowerSpectrum})';
            CI(:,1)= this_mean_dbWB'-CI(:,1);
            CI(:,2)=CI(:,2)- this_mean_dbWB';
            
            
            [hlCR, hpCR] = boundedline(freq',this_mean_dbWB', CI, 'r');
            
            title('Wavelet power spectrum')
            xlabel('Frequency (Hz)')
            ylabel('dB')
            
            %Uncomment if you want to save power spectra
            peakPowerSpectrum=handles.drgb.PACwave.peakPowerSpectrum;
            troughPowerSpectrum=handles.drgb.PACwave.troughPowerSpectrum;
            save('C:\Users\Diego Restrepo\OneDrive - The University of Colorado Denver\CaMKII Paper\Figure 3 PRP\elec5_splus_proficient.mat','peakPowerSpectrum',...
                'troughPowerSpectrum','freq')
            
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.83 .1 .05 .3])
            
            prain=[minLogPper:(maxLogPper-minLogPper)/99:maxLogPper];
            drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
            colormap fire
            shading interp
            ax=gca;
            set(ax,'XTickLabel','')
            
            
             figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
          
            %Plot the per-trial timecourse
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.07 .1 .75 .3])
            
            if handles.subtractRef==0
                drg_pcolor(repmat(t_apt,length(freq)*length(this_trialNo),1)',repmat(shifted_freq,length(t_apt),1),log_P_per_trial_timecourse')
            else
                %pcolor(repmat(t,length(f),1)',repmat(f,length(t),1),10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
                drg_pcolor(repmat(t_apt,length(freq)*length(this_trialNo),1)',repmat(shifted_freq,length(t_apt),1),log_P_per_trial_timecourse_sub')
                %imagesc(t,f,10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
            end
            
            colormap fire
            shading interp
            caxis([minLogPper maxLogPper]);
            xlabel('Time (sec)')
            ylabel('Frequency*trialNo');
            title(['Power (dB, wavelet) timecourse per trial ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
            
            figNo=figNo+1;
  
            try
                close(figNo)
            catch
            end
            
          
            %Plot the per-trial timecourse
            hFig = figure(figNo);            
            set(hFig, 'units','normalized','position',[.83 .1 .05 .3])
            
            prain=[minLogPper:(maxLogPper-minLogPper)/99:maxLogPper];
            drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
            colormap fire
            shading interp
            ax=gca;
            set(ax,'XTickLabel','')
            
            %Plot the lick traces
            figNo=figNo+1;
  
            try
                close(figNo)
            catch
            end
            
          
            %Plot the per-trial timecourse
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.07 .4 .6 .23])
            
            hold on
            
            per99=prctile(lick_traces(:),99.9);
            per1=prctile(lick_traces(:),1);
            
            mean_licks=zeros(1,length(lick_trace_times));
            
            
            y_shift=0;
            
            %Plot lick traces
            for ii=1:no_trials
                if sum(trials_included(ii)==handles.drgb.PAC.this_trialNo)>0
                    plot(lick_trace_times,lick_traces(ii,:)+y_shift,'-r')
                    y_shift=y_shift+1.5*(per99-per1);
                end
            end
            
            y_shift=y_shift+1.5*(per99-per1);
            ylim([0 y_shift])
            xlim([t_apt(1) t_apt(end)])
            xlabel('time(sec)')
            title('Lick traces')
            
            
            %Plot the lick frequency
            figNo=figNo+1;
  
            try
                close(figNo)
            catch
            end
            
          
            %Plot the per-trial timecourse
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.07 .7 .7 .23])
            
            
            if length(this_trialNo)==1
                plot(times_lick_freq, lick_freq)
            else
                [hl1, hp1] = boundedline(times_lick_freq',lick_freq', CIlickf', 'r');
            end
            
            ylim([0 1.2*max(lick_freq)])
            xlabel('Time (sec)')
            ylabel('frequency (Hz)')
            title('Lick frequency')
            
            %Plot the timecourse for power
            figNo=figNo+1;
  
            try
                close(figNo)
            catch
            end
            
          
            %Plot the per-trial timecourse
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.15 .6 .7 .23])
            
            hold on
            
          
            mean_lickPower=mean(lickPower,1)';
            if no_trials>2
                CI_lickPower = bootci(1000, {@mean, lickPower})';
                CI_lickPower(:,1)=mean_lickPower-CI_lickPower(:,1);
                CI_lickPower(:,2)=CI_lickPower(:,2)-mean_lickPower;
            end
            
            mean_peakPower=mean(peakPower,1)';
            if no_trials>2
                CI_peakPower = bootci(1000, {@mean, peakPower})';
                CI_peakPower(:,1)=mean_peakPower-CI_peakPower(:,1);
                CI_peakPower(:,2)=CI_peakPower(:,2)-mean_peakPower;
            end
            
            mean_troughPower=mean(troughPower,1)';
            if no_trials>2
                CI_troughPower = bootci(1000, {@mean, troughPower})';
                CI_troughPower(:,1)=mean_troughPower-CI_troughPower(:,1);
                CI_troughPower(:,2)=CI_troughPower(:,2)-mean_troughPower;
            end
            
            if no_trials>2
%                 [hllick, hplick]=boundedline(t_pac,mean_lickPower, CI_lickPower, 'g');
                
                [hltrough, hptrough]=boundedline(t_pac,mean_troughPower, CI_troughPower, 'b');
                [hlpeak, hppeak]=boundedline(t_pac,mean_peakPower, CI_peakPower, 'r');
            else
                
%                 plot(t_pac,mean_lickPower, 'g');
                
                plot(t_pac,mean_troughPower,  'b');
                plot(t_pac,mean_peakPower,  'r');
            end
            xlabel('Time(sec)')
            ylabel('Wavelet power (dB)')
            
            if no_trials>2
                legend([hltrough hlpeak],{'Trough','Peak'})
            else
                legend('Mean','Trough')
            end
            
            title(['Power (dB, wavelet)  ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
            
            if handles.autoscale==0
                ylim([handles.minLogP handles.maxLogP])
            end
            
             %Plot the timecourse for power
           figNo=figNo+1;
  
            try
                close(figNo)
            catch
            end
            
          
            %Plot the per-trial timecourse
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.15 .6 .7 .23])
            
            hold on
            
          
            mean_lickPower=mean(lickPower,1)';
            if no_trials>2
                CI_lickPower = bootci(1000, {@mean, lickPower})';
                CI_lickPower(:,1)=mean_lickPower-CI_lickPower(:,1);
                CI_lickPower(:,2)=CI_lickPower(:,2)-mean_lickPower;
            end
             
            if no_trials>2
                
                [hllick, hplick]=boundedline(t_pac,mean_lickPower, CI_lickPower, 'g');

            else
                
                plot(t_pac,mean_lickPower, 'g');
                
              
            end
            xlabel('Time(sec)')
            ylabel('Wavelet power (dB)')
            
          
            
            title(['Lick-referenced power (dB, wavelet)  ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
            
            if handles.autoscale==0
                ylim([handles.minLogP handles.maxLogP])
            end
            
        end
    end
end
pfft=1;



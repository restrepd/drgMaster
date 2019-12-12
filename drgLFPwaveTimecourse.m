function handles=drgLFPwaveTimecourse(handles)
%Generates a timecourse of the LFP power in decibels 10*log10(Power)


[t_apt,freq,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo]=drgGetLFPwavePowerForThisEvTypeNo(handles);



%Calculate the licks
[lick_freq,times_lick_freq,lick_traces,CIlickf,lick_trace_times,stamped_lick_ii,these_stamped_lick_times,no_trials,...
    trials_included,lick_threshold, lick_freq_per_trial,trials_included_per_trial]=drgGetLicks(handles);

%Get PAC
handles=drgThetaAmpPhaseTrialRange(handles);

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
    for trNum=1:handles.drgb.PAC.no_trials
        this_peakPower=[];
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
        
        handles.drgb.PACwave.PACtimecourse(trNum).peakPower=peakPower(trNum,:);
        handles.drgb.PACwave.meanPeakPower(trNum)=mean(peakPower(trNum,:));
    end
    
    %Find trough power
    troughPower=zeros(handles.drgb.PAC.no_trials,length(t_pac));
    for trNum=1:handles.drgb.PAC.no_trials
        this_troughPower=[];
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
    for trNum=1:length(trials_included_per_trial)
        ii_tr=find(handles.drgb.PAC.this_trialNo==trials_included_per_trial(trNum),1);
        if ~isempty(ii_tr)
            ntrs_out=ntrs_out+1;
            if handles.subtractRef==1
                lick_freq_per_trial(trNum,:)=lick_freq_per_trial(trNum,:)-mean(lick_freq_per_trial(trNum,(t_pac>=handles.startRef+handles.time_pad)&(t_pac<=handles.endRef-handles.time_pad)));
            end
            
            handles.drgb.PACwave.no_lick_trials=ntrs_out;
            handles.drgb.PACwave.lick_timecourse(ntrs_out).lick_f=lick_freq_per_trial(trNum,:);
            handles.drgb.PACwave.lick_timecourse(ntrs_out).times_lick_freq=times_lick_freq;
            handles.drgb.PACwave.mean_lick_freq(ntrs_out)=mean(lick_freq_per_trial(trNum,:));
            
            %Now enter peakPower and troughPower for trials that include licks
            handles.drgb.PACwave.meanPeakPower_per_lick_trial(ntrs_out)=handles.drgb.PACwave.meanPeakPower(ii_tr);
            handles.drgb.PACwave.meanTroughPower_per_lick_trial(ntrs_out)=handles.drgb.PACwave.meanTroughPower(ii_tr);
            
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
            
            try
                close 9
            catch
            end
            
            %Plot the timecourse
            hFig9 = figure(9);
            set(hFig9, 'units','normalized','position',[.07 .05 .75 .3])
            if handles.subtractRef==0
                drg_pcolor(repmat(t_apt,length(freq),1)',repmat(freq,length(t_apt),1),log_P_timecourse')
            else
                %pcolor(repmat(t,length(f),1)',repmat(f,length(t),1),10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
                drg_pcolor(repmat(t_apt,length(freq),1)',repmat(freq,length(t_apt),1),log_P_timecourse'-log_P_timecourse_ref')
                %imagesc(t,f,10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
            end
            
            colormap jet
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
            
            try
                close 10
            catch
            end
            
            hFig10 = figure(10);
            set(hFig10, 'units','normalized','position',[.83 .1 .05 .3])
            
            prain=[minLogPper:(maxLogPper-minLogPper)/99:maxLogPper];
            drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
            colormap jet
            shading interp
            ax=gca;
            set(ax,'XTickLabel','')
            
            
            try
                close 11
            catch
            end
            
            %Plot the per-trial timecourse
            hFig11 = figure(11);
            set(hFig11, 'units','normalized','position',[.07 .1 .75 .3])
            
            if handles.subtractRef==0
                drg_pcolor(repmat(t_apt,length(freq)*length(this_trialNo),1)',repmat(shifted_freq,length(t_apt),1),log_P_per_trial_timecourse')
            else
                %pcolor(repmat(t,length(f),1)',repmat(f,length(t),1),10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
                drg_pcolor(repmat(t_apt,length(freq)*length(this_trialNo),1)',repmat(shifted_freq,length(t_apt),1),log_P_per_trial_timecourse_sub')
                %imagesc(t,f,10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
            end
            
            colormap jet
            shading interp
            caxis([minLogPper maxLogPper]);
            xlabel('Time (sec)')
            ylabel('Frequency*trialNo');
            title(['Power (dB, wavelet) timecourse per trial ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
            
            try
                close 12
            catch
            end
            
            hFig12 = figure(12);
            set(hFig12, 'units','normalized','position',[.83 .1 .05 .3])
            
            prain=[minLogPper:(maxLogPper-minLogPper)/99:maxLogPper];
            drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
            colormap jet
            shading interp
            ax=gca;
            set(ax,'XTickLabel','')
            
            %Plot the lick traces
            try
                close 13
            catch
            end
            
            hFig13 = figure(13);
            set(hFig13, 'units','normalized','position',[.07 .4 .6 .23])
            
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
            try
                close 14
            catch
            end
            
            hFig14 = figure(14);
            set(hFig14, 'units','normalized','position',[.07 .7 .7 .23])
            
            
            if length(this_trialNo)==1
                plot(times_lick_freq, lick_freq)
            else
                [hl1, hp1] = boundedline(times_lick_freq',lick_freq', CIlickf', 'r');
            end
            
            ylim([0 1.2*max(lick_freq)])
            xlabel('Time (sec)')
            ylabel('frequency (Hz)')
            title('Lick frequency')
            
            %Plot the timecourse for the phase
            try
                close 15
            catch
            end
            
            hFig15 = figure(15);
            set(hFig15, 'units','normalized','position',[.15 .6 .7 .23])
            
            hold on
            
            try
                mean_calc=1;
                mean_meanPower=mean(meanPower,1)';
                if no_trials>2
                    CI_meanPower = bootci(1000, {@mean, meanPower})';
                    CI_meanPower(:,1)=mean_meanPower-CI_meanPower(:,1);
                    CI_meanPower(:,2)=CI_meanPower(:,2)-mean_meanPower;
                end
            catch
                mean_calc=0;
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
                if mean_calc==1
                    [hlmean, hpmean]=boundedline(t_pac,mean_meanPower, CI_meanPower, 'g');
                end
                [hltrough, hptrough]=boundedline(t_pac,mean_troughPower, CI_troughPower, 'b');
                [hlpeak, hppeak]=boundedline(t_pac,mean_peakPower, CI_peakPower, 'r');
            else
                if mean_calc==1
                    plot(t_pac,mean_meanPower, 'g');
                end
                plot(t_pac,mean_troughPower,  'b');
                plot(t_pac,mean_peakPower,  'r');
            end
            xlabel('Time(sec)')
            ylabel('Wavelet power (dB)')
            if mean_calc==1
                if no_trials>2
                    legend([hltrough hlpeak hlmean],{'Trough','Peak','Mean'})
                else
                    legend('Mean','Trough','Peak')
                end
            else
                if no_trials>2
                    legend([hltrough hlpeak],{'Trough','Peak'})
                else
                    legend('Trough','Peak')
                end
            end
            title(['Power (dB, wavelet)  ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
            
            if handles.autoscale==0
                ylim([handles.minLogP handles.maxLogP])
            end
            
            
            
        end
    end
end
pfft=1;



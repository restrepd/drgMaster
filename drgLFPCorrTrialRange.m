function handles=drgLFPCorrTrialRange(handles)

%Generates a trial per trial phase histogram
odorOn=2;

%Note that if event 5 exists the phase is referenced to event 5 (splus
%or non-match). Ohterwise phase is referenced to itself
if length(handles.drg.session(1).events)>=5
    splus=5;
else
    splus=handles.evTypeNo;
end
 
sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
lowF1=handles.peakLowF;
lowF2=handles.peakHighF;
highF1=handles.burstLowF;
highF2=handles.burstHighF;
pad_time=handles.time_pad;
n_phase_bins=handles.n_phase_bins;

%Empty vectors
handles.drgb.LFPcorr.no_trials=0;
handles.drgb.LFPcorr.rho=[];
handles.drgb.LFPcorr.t_lag=[];
handles.drgb.LFPcorr.trialNo=[];

%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;


[perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);
%handles=drgPercentLickPerTrial(handles);

% if handles.displayData==1
%     try
%         close 1
%     catch
%     end
%
%     hFig1 = figure(1);
%     set(hFig1, 'units','normalized','position',[.55 .27 .35 .15])
% end

 
%Now get the phase of gamma bursts witin theta
% no_encoding_trials=0;
% no_retrieval_trials=0;
no_trials=0;
% enc_phase_histo=[];
% retr_phase_histo=[];
% all_phase_histo=[];
% all_out_times=[];
% all_theta_wave=[];
% MI_enc=[];
% MI_retr=[];
% which_event=[];
% all_out_time_PAChisto=[];
% spm=[];
trials_attempted=0;
 
for trNo=firstTr:lastTr
     
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    
    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
         
        if excludeTrial==0
             
            
            
            trials_attempted=trials_attempted+1;
            [LFP1, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            %This excludes flat lick recordings
            if handles.peakLFPNo==19
               if sum(LFP1)==0
                   can_read1=0;
               end
            end
            [LFP2, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
             
            if (can_read1==1)&(can_read2==1)
                
                no_trials=no_trials+1;
                handles.drgb.LFPcorr.trialNo(no_trials)=trNo;
                
                if handles.displayData==1
                    fprintf(1, ['Trial No %d, event No %d, processed trial no %d\n'], trNo,evNo,no_trials);
                end
                
                %                 %If this is being used for batch processing find out whether
                %                 %this is an S+
                %                 eventNo = drgFindEvNo(handles,trNo,sessionNo,splus);
                %                 if eventNo~=-1
                %                     spm(no_trials)=1;
                %                 else
                %                     spm(no_trials)=0;
                %                 end
                
                time(no_trials)=handles.drg.session(sessionNo).trial_start(trNo);
                handles.drgb.LFPcorr.this_trialNo(no_trials)=trNo;
                if trNo==1
                    handles.drgb.LFPcorr.delta_t_trial(no_trials)=100; %There are no trials before the first trial
                else
                    handles.drgb.LFPcorr.delta_t_trial(no_trials)=handles.drg.session(sessionNo).trial_start(trNo)-handles.drg.session(sessionNo).trial_start(trNo-1);
                end
                
                [rho, t_lag]=drgGetLFPcorr(LFP1,LFP2,Fs,lowF1,lowF2,pad_time);
                
                weighted_mean_t_lag=sum(t_lag.*rho)/sum(rho);
                [max_rho, max_ii]=max(rho);
                max_rho_t_lag=t_lag(max_ii);
                
                %Save the output
                handles.drgb.LFPcorr.no_trials=no_trials;
                handles.drgb.LFPcorr.rho(no_trials,1:length(t_lag))=rho;
                handles.drgb.LFPcorr.weighted_mean_t_lag(no_trials)=weighted_mean_t_lag;
                handles.drgb.LFPcorr.max_rho_t_lag(no_trials)=max_rho_t_lag;
                handles.drgb.LFPcorr.t_lag=t_lag;
                
                
            else
                if handles.displayData==1
                    fprintf(1, ['Trial No %d, event No %d, LFP could not be read\n'], trNo,evNo);
                end
            end
        else
            if handles.displayData==1
                fprintf(1, ['Trial No %d, event No %d, trial excluded\n'], trNo,evNo);
            end
        end
    else
        if handles.displayData==1
            fprintf(1, ['Trial No %d, event No %d\n'], trNo,evNo);
        end
    end
    %end
    %end %if eventstamps...
end %for evNo





%Calculate correlation for shuffled trials
perm_trNo = randperm(length(handles.drgb.LFPcorr.trialNo));
no_sh_trials=0;
for trNo=firstTr:lastTr
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    
    
    
    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            
            if sum(perm_trNo==trNo)>0
                
                [LFP1, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
                %This excludes flat lick recordings
                if handles.peakLFPNo==19
                    if sum(LFP1)==0
                        can_read1=0;
                    end
                end
                
                evNo_sh = drgFindEvNo(handles,perm_trNo(trNo),sessionNo);
                [LFP2, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo_sh, handles.evTypeNo, handles.time_start, handles.time_end);
                
                if (can_read1==1)&(can_read2==1)
                    
                    no_sh_trials=no_sh_trials+1;
                    handles.drgb.LFPcorr.sh_trialNo(no_sh_trials)=trNo;
                    
                    if handles.displayData==1
                        fprintf(1, ['Trial No %d, event No %d, processed trial no %d\n'], trNo,evNo,no_sh_trials);
                    end
                    
                    %                 %If this is being used for batch processing find out whether
                    %                 %this is an S+
                    %                 eventNo = drgFindEvNo(handles,trNo,sessionNo,splus);
                    %                 if eventNo~=-1
                    %                     spm(no_sh_trials)=1;
                    %                 else
                    %                     spm(no_sh_trials)=0;
                    %                 end
                    
                    %                     time(no_sh_trials)=handles.drg.session(sessionNo).trial_start(trNo);
                    %                     handles.drgb.LFPcorr.this_trialNo(no_sh_trials)=trNo;
                    %                     if trNo==1
                    %                         handles.drgb.LFPcorr.delta_t_trial(no_sh_trials)=100; %There are no trials before the first trial
                    %                     else
                    %                         handles.drgb.LFPcorr.delta_t_trial(no_sh_trials)=handles.drg.session(sessionNo).trial_start(trNo)-handles.drg.session(sessionNo).trial_start(trNo-1);
                    %                     end
                    
                    
                    [rho_sh, t_lag_sh]=drgGetLFPcorr(LFP1,LFP2,Fs,lowF1,lowF2,pad_time);
                    
                    
                    weighted_mean_t_lag_sh=sum(t_lag_sh.*rho_sh)/sum(rho_sh);
                    [max_rho_sh, max_ii_sh]=max(rho_sh);
                    max_rho_t_lag_sh=t_lag_sh(max_ii_sh);
                    
                    %Save the output
                    handles.drgb.LFPcorr.no_sh_trials=no_sh_trials;
                    handles.drgb.LFPcorr.rho_sh(no_sh_trials,1:length(t_lag_sh))=rho_sh;
                    handles.drgb.LFPcorr.weighted_mean_t_lag_sh(no_sh_trials)=weighted_mean_t_lag_sh;
                    handles.drgb.LFPcorr.max_rho_t_lag_sh(no_sh_trials)=max_rho_t_lag_sh;
                    handles.drgb.LFPcorr.t_lag_sh=t_lag_sh;
                    
                    
                else
                    if handles.displayData==1
                        fprintf(1, ['Trial No %d, event No %d, LFP could not be read\n'], trNo,evNo);
                    end
                end
            end
        else
            if handles.displayData==1
                fprintf(1, ['Trial No %d, event No %d, trial excluded\n'], trNo,evNo);
            end
        end
    else
        if handles.displayData==1
            fprintf(1, ['Trial No %d, event No %d\n'], trNo,evNo);
        end
    end
    %end
    %end %if eventstamps...
end %for evNo



%Calculate the weigted mean t_lag
rho=handles.drgb.LFPcorr.rho;
mean_rho=mean(rho);

if handles.displayData==1
    fprintf(1, ['\nLFP crosscorrelation processed for %d out of %d trials \n\n'], no_trials,trials_attempted);
    
    try
        close 1
    catch
    end
    
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.05 .2 .2 .7])
    
   
    
    min_prob=prctile(rho(:),5);
    max_prob=prctile(rho(:),95);
    
    trials=[1:no_trials];
    %pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
    drg_pcolor(repmat(t_lag,length(trials),1),repmat(trials',1,length(t_lag)),rho)
    %         colormap jet
    colormap fire
    shading flat
    % min_prob=0.0113;
    % max_prob=0.0314;
    caxis([min_prob    max_prob])
    xlabel('Lag (sec)')
    ylabel('Trial No');
    title(['LFP cross-correlation'])
    
    
    try
        close 2
    catch
    end
    
    hFig2 = figure(2);
    set(hFig2, 'units','normalized','position',[.25 .65 .05 .3])
    
    
    prain=[min_prob:(max_prob-min_prob)/99:max_prob];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    %             colormap jet
    colormap fire
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')
    

    try
        close 3
    catch
    end
    
    hFig3 = figure(3);
    set(hFig3, 'units','normalized','position',[.27 .12 .27 .27])
    

    max_rho_t_lag=handles.drgb.LFPcorr.max_rho_t_lag;
     
    edges=[-0.0525:0.005:0.0525];
    histogram(max_rho_t_lag,edges)
    
    %Show mean rho
    this_ylim=ylim;
    hold on
    plot([mean(max_rho_t_lag) mean(max_rho_t_lag)],[0 this_ylim(2)],'-r')
    
    text(-0.05,this_ylim(1)+0.8*(this_ylim(2)-this_ylim(1)),['Mean rho (ms) = ' num2str(mean(max_rho_t_lag)*1000)])
    title('max rho lag time')
    xlabel('lag time (sec)')
    ylabel('Counts')
    
       try
        close 4
    catch
    end
    
    hFig4 = figure(4);
    set(hFig4, 'units','normalized','position',[.27 .12 .27 .27])
    

    max_rho_t_lag_sh=handles.drgb.LFPcorr.max_rho_t_lag_sh;
     
    edges=[-0.0525:0.005:0.0525];
    histogram(max_rho_t_lag_sh,edges)
    
    %Show mean rho
    this_ylim=ylim;
    hold on
    plot([mean(max_rho_t_lag_sh) mean(max_rho_t_lag_sh)],[0 this_ylim(2)],'-r')
    
    text(-0.05,this_ylim(1)+0.8*(this_ylim(2)-this_ylim(1)),['Mean rho (ms) = ' num2str(mean(max_rho_t_lag)*1000)])
    title('max rho lag time, shifted trials')
    xlabel('lag time (sec)')
    ylabel('Counts')
    
    p=ranksum(max_rho_t_lag,max_rho_t_lag_sh)
    
end

pfffft=1;






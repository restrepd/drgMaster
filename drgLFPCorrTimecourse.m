function handles=drgLFPCorrTimecourse(handles)

if ~isfield(handles,'randpermLFP')
    handles.randpermLFP=0;
end

%Generates a trial per trial lag histogram
%following the procedure of Adhikari et al. 2010
%https://doi.org/10.1016/j.jneumeth.2010.06.019
%Figures 4F,G

odorOn=2;
 
sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
lowF1=handles.peakLowF;
lowF2=handles.peakHighF;
highF1=handles.burstLowF;
highF2=handles.burstHighF;
pad_time=handles.time_pad;
n_phase_bins=handles.n_phase_bins;
decimation_factor=40;

%Empty vectors
handles.drgb.LFPcorr=[];

%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;

[perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);

no_trials=0;
trials_attempted=0;

vetted_trNos=[];
vetted_evNos=[];
ii_vet=0;
for trNo=firstTr:lastTr
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            [LFP1, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            
            [LFP2, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            
            LFP1=[];
            LFP2=[];
            
            if (can_read1==1)&(can_read2==1)
                ii_vet=ii_vet+1;
                vetted_trNos(ii_vet)=trNo;
                vetted_evNos(ii_vet)=evNo;
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

if handles.randpermLFP==1
    perm_ii_vet = randperm(length(vetted_trNos));
end

for ii_vet=1:length(vetted_trNos)
    
    trNo=vetted_trNos(ii_vet);
    evNo=vetted_evNos(ii_vet);
    
    no_trials=no_trials+1;
    
    handles.drgb.LFPcorr.no_trials=no_trials;
    handles.drgb.LFPcorr.perCorr(no_trials)=perCorr(drgFindEvNo(handles,trNo,sessionNo,odorOn));
    handles.drgb.LFPcorr.trial(no_trials).trialNo=trNo;
    if handles.save_drgb==1
        for evTypeNo=1:length(handles.drgbchoices.evTypeNos)
            if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(handles.drgbchoices.referenceEvent).times(evNo))>0
                handles.drgb.LFPcorr.which_event(evTypeNo,no_trials)=1;
            else
                handles.drgb.LFPcorr.which_event(evTypeNo,no_trials)=0;
            end
        end
    end
    
    if handles.displayData==1
        fprintf(1, ['Trial No %d, event No %d, processed trial no %d\n'], trNo,evNo,no_trials);
    end
    
    %Now parse the LFP timecourse
    dt=0.1;
    dt_window=2; %2 sec works much better than 1.5 sec
    time_start=handles.time_start;
    time_end=handles.time_end;
    no_steps=0;
    while handles.time_start+dt_window+2*handles.time_pad<time_end
        
        handles.time_end=handles.time_start+dt_window+2*handles.time_pad;
        [LFP1, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
        
        if handles.randpermLFP==1
            evNo_LFP2=vetted_evNos(perm_ii_vet(ii_vet));
        else
            evNo_LFP2=evNo;
        end
        
        [LFP2, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo_LFP2, handles.evTypeNo, handles.time_start, handles.time_end);
        
        LFP1=decimate(LFP1,decimation_factor);
        LFP2=decimate(LFP2,decimation_factor);
        no_steps=no_steps+1;
        time(no_steps)=handles.time_start+handles.time_pad+(dt_window/2);
        
        [rho, t_lag]=drgGetLFPcorr(LFP1,LFP2,Fs/decimation_factor,lowF1,lowF2,pad_time);
        
        weighted_mean_t_lag=sum(t_lag.*rho)/sum(rho);
        [max_rho, max_ii]=max(rho);
        max_rho_t_lag=t_lag(max_ii);
        
        
        %Save the output
        handles.drgb.LFPcorr.trial(no_trials).no_steps=no_steps;
        handles.drgb.LFPcorr.trial(no_trials).time=time;
        
        %handles.drgb.LFPcorr.trial(no_trials).rho(no_steps,1:length(t_lag))=rho;
        %I am not saving rho to save memory
        handles.drgb.LFPcorr.trial(no_trials).weighted_mean_t_lag(no_steps)=weighted_mean_t_lag;
        handles.drgb.LFPcorr.trial(no_trials).max_rho_t_lag(no_steps)=max_rho_t_lag;
        handles.drgb.LFPcorr.trial(no_trials).max_rho(no_steps)=max_rho;
        handles.drgb.LFPcorr.trial(no_trials).t_lag=t_lag;
        
        handles.time_start=handles.time_start+dt;
        
    end
    handles.time_start=time_start;
    handles.time_end=time_end; 
    
end %for ii_vet

handles.drgb.LFPcorr.time=time;

%Calculate the weigted mean t_lag

if handles.displayData==1

    max_rho_t_lag=zeros(no_trials,no_steps);
    max_rho=zeros(no_trials,no_steps);
    for ii=1:no_trials
        this_max_rho_t_lag=handles.drgb.LFPcorr.trial(ii).max_rho_t_lag;
        this_max_rho=handles.drgb.LFPcorr.trial(ii).max_rho;
        for jj=1:no_steps-1
            if abs(this_max_rho_t_lag(jj))>0.08
                if jj==1
                    this_max_rho_t_lag(jj)=0;
                else
                    this_max_rho_t_lag(jj)=this_max_rho_t_lag(jj-1);
                end
            end
        end

        max_rho_t_lag(ii,:)=this_max_rho_t_lag;
        max_rho(ii,:)=this_max_rho;
    end
    

    
%     %Plot lag time histogram
%     try
%         close 3
%     catch
%     end
%     
%     
%     hFig3 = figure(3);
%     set(hFig3, 'units','normalized','position',[.05 .33 .25 .45])
%     
%     plot(time,max_rho_t_lag)
    
    try
        close 4
    catch
    end
    
    
    hFig4 = figure(4);
    set(hFig4, 'units','normalized','position',[.05 .33 .25 .45])
    
    subplot(2,1,1)
    mean_max_rho_t_lag=mean(max_rho_t_lag);
    tempCI = bootci(1000, {@mean, max_rho_t_lag})';
    CI(:,1)=mean_max_rho_t_lag'-tempCI(:,1);
    CI(:,2)=tempCI(:,2)-mean_max_rho_t_lag';
    [hlCR, hpCR] = boundedline(time,mean_max_rho_t_lag, CI, 'b');
    
        xlabel('Time (sec)')
    ylabel('Lag time max rho')
    title('Lag timecourse')
    
  
    subplot(2,1,2)
        mean_max_rho=mean(max_rho);
    tempCI = bootci(1000, {@mean, max_rho})';
    CI(:,1)=mean_max_rho'-tempCI(:,1);
    CI(:,2)=tempCI(:,2)-mean_max_rho';
    [hlCR, hpCR] = boundedline(time,mean_max_rho, CI, 'b');
    
    xlabel('Time (sec)')
    ylabel('Max rho')
    title('Max rho timecourse')
    
    try
        close 1
    catch
    end
    
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.05 .05 .4 .6])

    
    
    min_prob=prctile(max_rho_t_lag(:),5);
    max_prob=prctile(max_rho_t_lag(:),95);
    
   
    %pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
    trials=1:no_trials;
    drg_pcolor(repmat(time,no_trials,1),repmat(trials',1,no_steps),max_rho_t_lag)
    %         colormap jet
    colormap fire
    shading flat
    % min_prob=0.0113;
    % max_prob=0.0314;
    caxis([min_prob    max_prob])
    xlabel('Time (sec)')
    ylabel('Trial');
    title(['Lag (sec)'])
    
    
    try
        close 2
    catch
    end
    
    hFig2 = figure(2);
    set(hFig2, 'units','normalized','position',[.45 .05 .05 .3])
    
    
    prain=[min_prob:(max_prob-min_prob)/99:max_prob];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    %             colormap jet
    colormap fire
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')
    
   
    
end
   
pfffft=1;






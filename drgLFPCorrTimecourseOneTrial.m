function handles=drgLFPCorrTimecourseOneTrial(handles)

if ~isfield(handles,'flipLFP')
    handles.flipLFP=0;
end

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
 
%Find the first instance of the event the user chose
trNo=firstTr-1;
foundTrial=0;   
while foundTrial==0
    trNo=trNo+1;
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
                foundTrial=1;
%                 
%                 no_trials=no_trials+1;
%                 handles.drgb.LFPcorr.trialNo(no_trials)=trNo;
%                 
%                 if handles.displayData==1
%                     fprintf(1, ['Trial No %d, event No %d, processed trial no %d\n'], trNo,evNo,no_trials);
%                 end
%                 
%                 %                 %If this is being used for batch processing find out whether
%                 %                 %this is an S+
%                 %                 eventNo = drgFindEvNo(handles,trNo,sessionNo,splus);
%                 %                 if eventNo~=-1
%                 %                     spm(no_trials)=1;
%                 %                 else
%                 %                     spm(no_trials)=0;
%                 %                 end
%                 
%                 time(no_trials)=handles.drg.session(sessionNo).trial_start(trNo);
%                 handles.drgb.LFPcorr.this_trialNo(no_trials)=trNo;
%                 if trNo==1
%                     handles.drgb.LFPcorr.delta_t_trial(no_trials)=100; %There are no trials before the first trial
%                 else
%                     handles.drgb.LFPcorr.delta_t_trial(no_trials)=handles.drg.session(sessionNo).trial_start(trNo)-handles.drg.session(sessionNo).trial_start(trNo-1);
%                 end
%                 
%                 [rho, t_lag]=drgGetLFPcorr(LFP1,LFP2,Fs,lowF1,lowF2,pad_time);
%                 
%                 weighted_mean_t_lag=sum(t_lag.*rho)/sum(rho);
%                 [max_rho, max_ii]=max(rho);
%                 max_rho_t_lag=t_lag(max_ii);
%                 
%                 %Save the output
%                 handles.drgb.LFPcorr.no_trials=no_trials;
%                 handles.drgb.LFPcorr.rho(no_trials,1:length(t_lag))=rho;
%                 handles.drgb.LFPcorr.weighted_mean_t_lag(no_trials)=weighted_mean_t_lag;
%                 handles.drgb.LFPcorr.max_rho_t_lag(no_trials)=max_rho_t_lag;
%                 handles.drgb.LFPcorr.max_rho(no_trials)=max_rho;
%                 handles.drgb.LFPcorr.t_lag=t_lag;
%                 
%                 
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


end

no_trials=no_trials+1;
handles.drgb.LFPcorr_tc.trialNo=trNo;

if handles.displayData==1
    fprintf(1, ['Trial No %d, event No %d, processed trial no %d\n'], trNo,evNo,no_trials);
end

%Now parse the LFP timecourse
dt=0.05;
dt_window=2; %2 sec works much better than 1.5 sec
time_start=handles.time_start;
time_end=handles.time_end;
no_steps=0;
while handles.time_start+dt_window+2*handles.time_pad<time_end
    
    handles.time_end=handles.time_start+dt_window+2*handles.time_pad;
    [LFP1, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
    
    %This excludes flat lick recordings
    if handles.peakLFPNo==19
        if sum(LFP1)==0
            can_read1=0;
        end
    end
    
    [LFP2, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
    
    if handles.flipLFP==1
        LFP2=LFP2(end:-1:1);
    end
    
    if (can_read1==1)&(can_read2==1)
        no_steps=no_steps+1;
        time(no_steps)=handles.time_start+handles.time_pad+(dt_window/2);
        
        [rho, t_lag]=drgGetLFPcorr(LFP1,LFP2,Fs,lowF1,lowF2,pad_time);
         
        weighted_mean_t_lag=sum(t_lag.*rho)/sum(rho);
        [max_rho, max_ii]=max(rho);
        max_rho_t_lag=t_lag(max_ii);

        
        %Save the output
        handles.drgb.LFPcorr_tc.no_steps=no_steps;
        handles.drgb.LFPcorr_tc.time=time;
        
        handles.drgb.LFPcorr_tc.rho(no_steps,1:length(t_lag))=rho;
        handles.drgb.LFPcorr.weighted_mean_t_lag(no_steps)=weighted_mean_t_lag;
        handles.drgb.LFPcorr.max_rho_t_lag(no_steps)=max_rho_t_lag;
        handles.drgb.LFPcorr.max_rho(no_steps)=max_rho;
        handles.drgb.LFPcorr.t_lag=t_lag;
        

    end
    
    handles.time_start=handles.time_start+dt;
    
end
      
handles.drgb.LFPcorr.time=time;

%Calculate the weigted mean t_lag

if handles.displayData==1

    %Plot lag time histogram
    try
        close 3
    catch
    end
    
    
    hFig3 = figure(3);
    set(hFig3, 'units','normalized','position',[.05 .33 .25 .45])
    
    subplot(2,1,1)
    max_rho_t_lag=handles.drgb.LFPcorr.max_rho_t_lag;
    plot(time(abs(max_rho_t_lag)~=0.1),max_rho_t_lag(abs(max_rho_t_lag)~=0.1),'-ob')
    xlabel('Time (sec)')
    ylabel('Lag (sec)')
    title('Lag timecourse')
    
    subplot(2,1,2)
    max_rho=handles.drgb.LFPcorr.max_rho;
    plot(time(abs(max_rho_t_lag)~=0.1),max_rho(abs(max_rho_t_lag)~=0.1),'-ob')
    xlabel('Time (sec)')
    ylabel('Max rho')
    title('Max rho timecourse')
    
    try
        close 1
    catch
    end
    
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.05 .05 .6 .25])

    rho=handles.drgb.LFPcorr_tc.rho;
    
    min_prob=prctile(rho(:),5);
    max_prob=prctile(rho(:),95);
    
   
    %pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
    drg_pcolor(repmat(time,length(t_lag),1),repmat(t_lag',1,no_steps),rho')
    %         colormap jet
    colormap fire
    shading flat
    % min_prob=0.0113;
    % max_prob=0.0314;
    caxis([min_prob    max_prob])
    xlabel('Time (sec)')
    ylabel('Lag (sec)');
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
    
   
    
end
   
pfffft=1;






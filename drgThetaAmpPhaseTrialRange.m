function handles=drgThetaAmpPhaseTrialRange(handles)
 
%Generates a trial per trial phase histogram
odorOn=2;
sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
lowF1=handles.peakLowF;
lowF2=handles.peakHighF;
highF1=handles.burstLowF;
highF2=handles.burstHighF;
pad_time=handles.time_pad;
n_phase_bins=handles.n_phase_bins;

%Empty vectors
handles.drgb.PAC.no_trials=0;
handles.drgb.PAC.meanVectorLength=[];
handles.drgb.PAC.meanVectorAngle=[];
handles.drgb.PAC.peakAngle=[];
handles.drgb.PAC.mod_indx=[];
handles.drgb.PAC.this_trialNo=[];
handles.drgb.PAC.all_phase_histo=[];
handles.drgb.PAC.all_theta_wave=[];
handles.drgb.PAC.perCorr=[];
handles.drgb.PAC.which_event=[];

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
no_encoding_trials=0;
no_retrieval_trials=0;
no_trials=0;
enc_phase_histo=[];
retr_phase_histo=[];
all_phase_histo=[];
all_out_times=[];
all_theta_wave=[];
MI_enc=[];
MI_retr=[];
which_event=[];


for trNo=firstTr:lastTr
    
    if handles.displayData==1
        trial_no=trNo
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);

    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            [LFPlow, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            [LFPhigh, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            
            if (can_read1==1)&(can_read2==1)
                no_trials=no_trials+1;
                time(no_trials)=handles.drg.session(sessionNo).trial_start(trialNo);
                handles.drgb.PAC.this_trialNo(no_trials)=trNo;
                perCorr_per_histo(no_trials)=50;
                if handles.peakLFPNo==18
                    %This is the sniff
                    [meanVectorLength(no_trials), meanVectorAngle(no_trials), peakAngle(no_trials), mod_indx(no_trials), phase, phase_histo, theta_wave]=drgGetThetaAmpPhaseSniff(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,pad_time,n_phase_bins,handles.which_method);
                else
                    [meanVectorLength(no_trials), meanVectorAngle(no_trials), peakAngle(no_trials), mod_indx(no_trials), phase, phase_histo, theta_wave, meanPeakAngle, out_times, out_phase, out_time_PAChisto]=drgGetThetaAmpPhase(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,pad_time,n_phase_bins,handles.which_method);
                end
                
                out_times=out_times+handles.time_start+handles.time_pad;
                
                %Save the output
                handles.drgb.PAC.no_trials=no_trials;
                handles.drgb.PAC.meanVectorLength(no_trials)=meanVectorLength(no_trials);
                handles.drgb.PAC.meanVectorAngle(no_trials)=meanVectorAngle(no_trials);
                handles.drgb.PAC.peakAngle(no_trials)=peakAngle(no_trials);
                handles.drgb.PAC.mod_indx(no_trials)=mod_indx(no_trials);
                handles.drgb.PAC.all_phase_histo(no_trials,1:n_phase_bins+1)=phase_histo;
                handles.drgb.PAC.all_theta_wave(no_trials,1:n_phase_bins+1)=theta_wave;
                handles.drgb.PAC.perCorr(no_trials)=perCorr(drgFindEvNo(handles,trialNo,sessionNo,odorOn));
                handles.drgb.PAC.meanPeakAngle(no_trials)=meanPeakAngle;
                handles.drgb.PAC.PACtimecourse(no_trials).out_times=out_times;
                handles.drgb.PAC.PACtimecourse(no_trials).out_phase=out_phase;
                handles.drgb.PAC.PACtimecourse(no_trials).out_time_PAChisto=out_time_PAChisto;
   
                all_out_phase(no_trials,1:length(out_times))=out_phase;
                if no_trials==1
                    all_out_time_PAChisto=[];
                    all_out_time_PAChisto=out_time_PAChisto;
                else
                all_out_time_PAChisto=all_out_time_PAChisto+out_time_PAChisto;
                end
                all_phase_histo(no_trials,1:n_phase_bins+1)=phase_histo;
                all_theta_wave(no_trials,1:n_phase_bins+1)=theta_wave;
                no_encoding_trials=no_encoding_trials+1;
                enc_phase_histo(no_encoding_trials,1:n_phase_bins+1)=phase_histo;
                MI_enc=[MI_enc mod_indx(no_trials)];
                
                if handles.save_drgb==1
                    for evTypeNo=1:length(handles.drgbchoices.evTypeNos)
                        if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(handles.drgbchoices.referenceEvent).times(evNo))>0
                            handles.drgb.PAC.which_event(evTypeNo,no_trials)=1;
                        else
                            handles.drgb.PAC.which_event(evTypeNo,no_trials)=0;
                        end
                    end
                end
            end
        end
    end
    %end
    %end %if eventstamps...
end %for evNo

all_out_time_PAChisto=all_out_time_PAChisto/no_trials;

mean_MI_enc=mean(MI_enc);

if handles.displayData==1
    
     try
        close 1
    catch
    end
    
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.15 .6 .8 .23])
    
    plot(out_times,circ_rad2ang(circ_mean(circ_ang2rad(all_out_phase)-180))+180,'ob')
    ylim([0 360])
    
    try
        close 6
    catch
    end
    
    hFig6 = figure(6);
    set(hFig6, 'units','normalized','position',[.15 .2 .8 .23])
    
    min_prob=prctile(all_out_time_PAChisto(:),5);
    max_prob=prctile(all_out_time_PAChisto(:),95);
    
  
    %pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
    drg_pcolor(repmat(out_times,length(phase(1:end-1)),1),repmat(phase(1:end-1)',1,length(out_times)),all_out_time_PAChisto')
    colormap jet
    shading flat
    % min_prob=0.0113;
    % max_prob=0.0314;
    caxis([min_prob    max_prob])
    xlabel('Time (sec)')
    ylabel('Phase for low freq oscillation (deg)');
    title(['Phase-amplitude coupling timecourse ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
   
%     (180/pi)*circ_mean(all_out_phase'*pi/180)'
    try
        close 2
    catch
    end
    
    hFig2 = figure(2);
    set(hFig2, 'units','normalized','position',[.25 .1 .23 .8])
    
    %Plot the modulation index as a function of trial number
    subplot(4,1,1)
    trials=1:no_trials;
    plot(trials,mod_indx,'o-')
    xlabel('Trial No')
    ylabel('Modulation index')
    title('Modulation index vs trial number')
    
    % figure(11)
    % plot(trials,mod_indx,'o-k','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7,'LineWidth',2)
    % xlabel('Trial No')
    % ylabel('Modulation index')
    % title('Modulation index vs trial number')
    % ylim([0 0.07])
    % set(gca,'FontSize',20,'LineWidth',1.25)
    
    
    
    %Now plot the mean theta waveform
    subplot(4,1,2)
    shadedErrorBar(phase,mean(all_theta_wave,1),std(all_theta_wave,0,1),'-b')
    xlim([0 360])
    title('Mean low frequency waveform')
    xlabel('Degrees')
    
    %Now plot the encoding theta phase histogram for gamma
    subplot(4,1,3)
    shadedErrorBar(phase,mean(enc_phase_histo,1),std(enc_phase_histo,0,1),'-b')
    xlim([0 360])
    title('Phase-amplitude coupling')
    ylabel('Probability')
    xlabel('Phase for low freq oscillation (deg)')
    yl_enc=ylim;
    text(20,yl_enc(1)+0.9*(yl_enc(2)-yl_enc(1)),['MI= ' num2str(mean_MI_enc)])
    
    
    
    try
        close 3
    catch
    end
    
    hFig3 = figure(3);
    set(hFig3, 'units','normalized','position',[.01 .1 .23 .8])
    
    min_prob=prctile(all_phase_histo(:),5);
    max_prob=prctile(all_phase_histo(:),95);
    
    trials=1:no_trials;
    %pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
    drg_pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
    colormap jet
    shading flat
    % min_prob=0.0113;
    % max_prob=0.0314;
    caxis([min_prob    max_prob])
    xlabel('Phase for low freq oscillation (deg)')
    ylabel('Trial No');
    title(['Phase-amplitude coupling per trial' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
    
    try
        close 4
    catch
    end
    
    hFig4 = figure(4);
    set(hFig4, 'units','normalized','position',[.49 .1 .05 .3])
    
    
    prain=[min_prob:(max_prob-min_prob)/99:max_prob];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    colormap jet
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')
    
    try
        close 5
    catch
    end
    
    hFig5 = figure(5);
    set(hFig5, 'units','normalized','position',[.59 .05 .35 .35])
    
    if no_encoding_trials>0
        rose(pi*meanVectorAngle/180,12)
        title('Mean phase')
    end
end

if handles.displayData==1
    fprintf(1, ['\nPCA processed for %d out of %d trials \n\n'], no_trials,lastTr);
end

pfffft=1;






function handles=drgGetThetaAmpPhaseAllTrials(handles)
%Gets PAC for all trials using either single electrode LFP or or all electrode LFPs per tetrode\

sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
lowF1=handles.peakLowF;
lowF2=handles.peakHighF;
highF1=handles.burstLowF;
highF2=handles.burstHighF;
pad_time=handles.time_pad;
n_phase_bins=handles.n_phase_bins;
odorOn=2;

%First LFP in this tetrode
firstLFP=4*ceil(handles.peakLFPNo/4)-3;

%Get the information for percent correct
[perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);


handles=drgPercentLickPerTrial(handles);


%Now get the phase of gamma bursts witin theta
no_trials=0;
all_phase_histo=[];
all_theta_wave=[];
which_trial=[];
which_event=[];



for evNo=1:handles.drg.session(handles.sessionNo).events(handles.evTypeNo).noTimes
    
    if handles.save_drgb==0
        event_no=evNo
    end
    
    if handles.perTetrode==1
        noLFPs=0;
        meanVectorLength=[];
        meanVectorAngle=[];
        peakAngle=[];
        mod_indx=[];
        all_LFPs_phase_histo=[];
        all_LFPs_theta_wave=[];
        for thisLFP=firstLFP:firstLFP+3
            excludeTrial=drgExcludeTrialLFP(handles.drg,thisLFP,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
            if excludeTrial==0
                
                [LFPlow, trialNo, can_read1] = drgGetTrialLFPData(handles, thisLFP, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
                [LFPhigh, trialNo, can_read2] = drgGetTrialLFPData(handles, thisLFP, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
                
                if (can_read1==1)&(can_read2==1)
                    noLFPs=noLFPs+1;
                    [meanVectorLength, meanVectorAngle, peakAngle,  mod_indx, handles.drgb.PAC.phase, phase_histo, theta_wave]=drgGetThetaAmpPhase(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,pad_time,n_phase_bins,handles.which_method);
                    all_LFPs_phase_histo(noLFPs,1:n_phase_bins+1)=phase_histo;
                    all_LFPs_theta_wave(noLFPs,1:n_phase_bins+1)=theta_wave;
                end
            end
        end
        
        if noLFPs>=1
            no_trials=no_trials+1;
            handles.drgb.PAC.no_trials=no_trials;
            
            phase=handles.drgb.PAC.phase;
            phase_histo=mean(all_LFPs_phase_histo,1);
            
            phaseAngle=pi*phase/180;
            meanVectorLength=sqrt((mean(phase_histo.*sin(phaseAngle)))^2+(mean(phase_histo.*cos(phaseAngle)))^2);
            meanX=mean(phase_histo.*cos(phaseAngle));
            meanY=mean(phase_histo.*sin(phaseAngle));
            if meanY>0
                meanVectorAngle=(180/pi)*acos(meanX/sqrt(meanY^2+meanX^2));
            else
                meanVectorAngle=360-(180/pi)*acos(meanX/sqrt(meanY^2+meanX^2));
            end
            handles.drgb.PAC.meanVectorLength(no_trials)=meanVectorLength;
            handles.drgb.PAC.meanVectorAngle(no_trials)=meanVectorAngle;
            
            
            [max_hist peak_bin]=max(phase_histo);
            handles.drgb.PAC.peakAngle(no_trials)=phase(peak_bin);
            
            mean_prob=mean(phase_histo)*ones(1,length(phase_histo));
            DKL=sum(phase_histo(1:end-1).*log(phase_histo(1:end-1)./mean_prob(1:end-1)));
            MI_Tort=DKL/log(n_phase_bins);
            handles.drgb.PAC.mod_indx(no_trials)=MI_Tort;
            
            handles.drgb.PAC.all_phase_histo(no_trials,1:n_phase_bins+1)=mean(all_LFPs_phase_histo,1);
            handles.drgb.PAC.all_theta_wave(no_trials,1:n_phase_bins+1)=mean(all_LFPs_theta_wave,1);
            handles.drgb.PAC.perCorr(no_trials)=perCorr(drgFindEvNo(handles,trialNo,sessionNo,odorOn));
            handles.drgb.PAC.percent_lick(no_trials)=handles.drg.session(sessionNo).percent_lick(trialNo);
            
            if handles.displayData==0
                for evTypeNo=1:length(handles.drgbchoices.evTypeNos)
                    if sum(handles.drg.session(1).events(handles.drgbchoices.evTypeNos(evTypeNo)).times==handles.drg.session(1).events(handles.drgbchoices.referenceEvent).times(evNo))>0
                        handles.drgb.PAC.which_event(evTypeNo,no_trials)=1;
                    else
                        handles.drgb.PAC.which_event(evTypeNo,no_trials)=0;
                    end
                end
            end
        end
        
    else
        
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            [LFPlow, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            [LFPhigh, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            
            if (can_read1==1)&(can_read2==1)
                
                [meanVectorLength, meanVectorAngle, peakAngle,  mod_indx, handles.drgb.PAC.phase, phase_histo, theta_wave]=drgGetThetaAmpPhase(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,pad_time,n_phase_bins,handles.which_method);
                no_trials=no_trials+1;
                handles.drgb.PAC.no_trials=no_trials;
                handles.drgb.PAC.meanVectorLength(no_trials)=meanVectorLength;
                handles.drgb.PAC.meanVectorAngle(no_trials)=meanVectorAngle;
                handles.drgb.PAC.peakAngle(no_trials)=peakAngle;
                handles.drgb.PAC.mod_indx(no_trials)=mod_indx;
                handles.drgb.PAC.all_phase_histo(no_trials,1:n_phase_bins+1)=phase_histo;
                handles.drgb.PAC.all_theta_wave(no_trials,1:n_phase_bins+1)=theta_wave;
                
                handles.drgb.PAC.perCorr(no_trials)=perCorr(abs(handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo)-handles.drg.session(sessionNo).events(2).times)...
                            <= handles.max_dt_between_events);
                handles.drgb.PAC.percent_lick(no_trials)=handles.drg.session(sessionNo).percent_lick(trialNo);
                
                if handles.displayData==0
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






function handles=drgClosedLoopPhase(handles)

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
handles.drgb.CL.no_trials=0;
handles.drgb.CL.laser_phases=[];
handles.drgb.CL.all_phase_histo=[];
handles.drgb.CL.all_theta_wave=[];


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
all_laser_phases=[];
all_theta_wave=[];
MI_enc=[];
MI_retr=[];
which_event=[];
all_out_time_PAChisto=[];
spm=[];
trials_attempted=0;

for trNo=firstTr:lastTr
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    
    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            trials_attempted=trials_attempted+1;
            [LFPlow, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            %This excludes flat lick recordings
            if handles.peakLFPNo==19
                if sum(LFPlow)==0
                    can_read1=0;
                end
            end
            
            %Closed loop laser triggers are in channel 18
            
            this_burstLFPNo=handles.burstLFPNo;
            handles.burstLFPNo=18;
            [LFPhigh, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            handles.burstLFPNo=this_burstLFPNo;
            
            if (can_read1==1)&(can_read2==1)
                
                no_trials=no_trials+1;
                
                if handles.displayData==1
                    fprintf(1, ['Trial No %d, event No %d, processed trial no %d\n'], trNo,evNo,no_trials);
                end
                
                %If this is being used for batch processing find out whether
                %this is an S+
                eventNo = drgFindEvNo(handles,trNo,sessionNo,splus);
                if eventNo~=-1
                    spm(no_trials)=1;
                else
                    spm(no_trials)=0;
                end
                
                time(no_trials)=handles.drg.session(sessionNo).trial_start(trNo);
                handles.drgb.CL.this_trialNo(no_trials)=trNo;
                if trNo==1
                    handles.drgb.CL.delta_t_trial(no_trials)=100; %There are no trials before the first trial
                else
                    handles.drgb.CL.delta_t_trial(no_trials)=handles.drg.session(sessionNo).trial_start(trNo)-handles.drg.session(sessionNo).trial_start(trNo-1);
                end
                perCorr_per_histo(no_trials)=50;
                
                
                [meanVectorLength(no_trials), meanVectorAngle(no_trials), peakAngle(no_trials), mod_indx(no_trials), phase,...
                    phase_histo, theta_wave, meanPeakAngle, out_times, out_phase, out_time_PAChisto, decLFPgenv, decanglethetaLFP, out_times_env,laser_phases]...
                    =drgGetClosedLoopPhase(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,pad_time,n_phase_bins,handles.which_method);
                
                
                out_times=out_times+handles.time_start+handles.time_pad;
                
                %Save the output
                handles.drgb.CL.no_trials=no_trials;
                handles.drgb.CL.laser_phases=[handles.drgb.CL.laser_phases laser_phases];
                handles.drgb.CL.all_phase_histo(no_trials,1:n_phase_bins+1)=phase_histo;
                handles.drgb.CL.all_theta_wave(no_trials,1:n_phase_bins+1)=theta_wave;
                handles.drgb.CL.phase=phase;
                
                
                all_phase_histo(no_trials,1:n_phase_bins+1)=phase_histo;
                all_theta_wave(no_trials,1:n_phase_bins+1)=theta_wave;
                all_laser_phases=[all_laser_phases laser_phases];
                
%                 try
%                     close(23)
%                 catch
%                 end
%                 hFig=figure(23);
%                 histogram(all_laser_phases)
                
                
                pffft=1;
                
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


if handles.displayData==1
    
    
    try
        close 2
    catch
    end
    
    hFig2 = figure(2);
    set(hFig2, 'units','normalized','position',[.25 .1 .23 .8])
    
    
    %Now plot the mean theta waveform
    subplot(2,1,1)
    shadedErrorBar(phase,mean(all_theta_wave,1),std(all_theta_wave,0,1),'-b')
    xlim([0 360])
    title('Mean low frequency waveform')
    xlabel('Degrees')
    
    %Now plot the encoding theta phase histogram for gamma
    subplot(2,1,2)
    edges=[0:10:360];
    histogram(all_laser_phases,edges)
    xlim([0 360])
    title('Phase for laser pulse')
    ylabel('Number')
    xlabel('Phase (deg)')
    
    pffft=1;
end







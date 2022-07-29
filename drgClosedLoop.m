function handles=drgClosedLoop(handles)

%Generates a trial per trial phase histogram
odorOn=2;

% laser_method=1 for data acquired with DT9812
% laser_method=2 for data acquired with Arduino UNO

laser_method=2;

all_Power_timecourse=[];

%Parameters for laser power
sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;

%LFP power parameters
handles.LFPPowerSpectrumLowF=4;
handles.LFPPowerSpectrumHighF=95;

handles.no_lfp_bands=4;

%Theta
handles.lowF(1)=6;
handles.highF(1)=14;

%Beta
handles.lowF(2)=15;
handles.highF(2)=30;

%Low gamma
handles.lowF(3)=35;
handles.highF(3)=55;

%High gamma
handles.lowF(4)=65;
handles.highF(4)=95;

handles.f_label{1}='Theta';
handles.f_label{2}='Beta';
handles.f_label{3}='Low gamma';
handles.f_label{4}='High gamma';


dec_n=fix(handles.drg.session(sessionNo).draq_p.ActualRate/1000);

%Setup the wavelet scales
%   scales = helperCWTTimeFreqVector(minfreq,maxfreq,f0,dt,NumVoices)
%   f0 - center frequency of the wavelet in cycles/unit time
%   dt - sampling interval
%   NumVoices - number of voices per octave

NumVoices=5;
minfreq=handles.LFPPowerSpectrumLowF;
maxfreq=handles.LFPPowerSpectrumHighF;
dt=1/Fs;
f0=5/(2*pi);

a0 = 2^(1/NumVoices);
minscale = f0/(maxfreq*dt);
maxscale = f0/(minfreq*dt);
minscale = floor(NumVoices*log2(minscale));
maxscale = ceil(NumVoices*log2(maxscale));
scales = a0.^(minscale:maxscale).*dt;


%Note that if event 5 exists the phase is referenced to event 5 (splus
%or non-match). Ohterwise phase is referenced to itself
if length(handles.drg.session(1).events)>=5
    splus=5;
else
    splus=handles.evTypeNo;
end



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
            [LFPhigh, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
             
            %Read the laser
            switch laser_method
                case 1
                    laserLFPNo=18;
                case 2
                    laserLFPNo=21;
            end
            
            [LFPlaser, trialNo, can_read3] = drgGetTrialLFPData(handles, laserLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            
            if (can_read1==1)&(can_read2==1)&(can_read3==1)
                
                
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
                handles.drgb.PAC.this_trialNo(no_trials)=trNo;
                if trNo==1
                    handles.drgb.PAC.delta_t_trial(no_trials)=100; %There are no trials before the first trial
                else
                    handles.drgb.PAC.delta_t_trial(no_trials)=handles.drg.session(sessionNo).trial_start(trNo)-handles.drg.session(sessionNo).trial_start(trNo-1);
                end
                perCorr_per_histo(no_trials)=50;
               
                switch laser_method
                    case 1
                        if handles.peakLFPNo==18
                            %This is the sniff
                            [meanVectorLength(no_trials), meanVectorAngle(no_trials), peakAngle(no_trials), mod_indx(no_trials), phase, phase_histo, theta_wave]=drgGetThetaAmpPhaseSniff(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,pad_time,n_phase_bins,handles.which_method);
                        else
                            [meanVectorLength(no_trials), meanVectorAngle(no_trials), peakAngle(no_trials), mod_indx(no_trials), phase,...
                                phase_histo, theta_wave, meanPeakAngle, out_times, out_phase, out_time_PAChisto, decLFPgenv, decanglethetaLFP, out_times_env]...
                                =drgGetThetaAmpPhase(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,pad_time,n_phase_bins,handles.which_method);
                        end
                    case 2
                        [meanVectorLength(no_trials), meanVectorAngle(no_trials), peakAngle(no_trials), mod_indx(no_trials), phase,...
                            phase_histo, theta_wave, meanPeakAngle, out_times, out_phase, out_time_PAChisto, decLFPgenv, decanglethetaLFP, out_times_env]...
                            =drgGetThetaAmpPhase(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,pad_time,n_phase_bins,handles.which_method);
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
                handles.drgb.PAC.perCorr(no_trials)=perCorr(drgFindEvNo(handles,trNo,sessionNo,odorOn));
                handles.drgb.PAC.meanPeakAngle(no_trials)=meanPeakAngle;
                handles.drgb.PAC.PACtimecourse(no_trials).out_times=out_times;
                handles.drgb.PAC.PACtimecourse(no_trials).out_phase=out_phase;
                handles.drgb.PAC.PACtimecourse(no_trials).out_time_PAChisto=out_time_PAChisto;
                handles.drgb.PAC.PACtimecourse(no_trials).decLFPgenv=decLFPgenv;
                handles.drgb.PAC.PACtimecourse(no_trials).decanglethetaLFP=decanglethetaLFP;
                handles.drgb.PAC.out_times_env=out_times_env;
                handles.drgb.PAC.phase=phase;
                
                all_out_phase(no_trials,1:length(out_times))=out_phase;
                if no_trials==1
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
                 
                %Now figure out the number and phase of laser pulses
                
              
                
                switch laser_method
                    
                    case 1
                        dt_all_pulses=40;
                        dt_down=10;
                        dt_next_pulse=27;
                        at_end=0;
                        ii=1;
                        handles.drgb.PAC.laser(no_trials).laser_times=[];
                        handles.drgb.PAC.laser(no_trials).laser_phase=[];
                        ii_laser=0;
                        thisLFPlaser=LFPlaser(pad_time*Fs:pad_time*Fs+out_times_env(end)*Fs-1);
                        Fs_out_times_env=1/(out_times_env(2)-out_times_env(1));
                        dec_factor=floor(Fs/Fs_out_times_env);
                        thisLFPlaser_dec=decimate(thisLFPlaser,dec_factor);
                        
                        if thisLFPlaser_dec(1)>1
                            this_laser_down_jj=find(thisLFPlaser_dec(ii:end)<1,1,'first');
                            ii=ii+this_laser_down_jj;
                        end
%                         figure(11)
%                         plot(thisLFPlaser_dec)
%                         figure(12)
%                         plot(decanglethetaLFP)
                        while at_end==0
                            this_laser_jj=find(thisLFPlaser_dec(ii:end)>1,1,'first');
                            if ~isempty(this_laser_jj)
                                
                                if (ii+this_laser_jj-1 +dt_all_pulses) < length(out_times_env)
                                    if (thisLFPlaser_dec(ii+this_laser_jj-1+dt_down)<1)&(thisLFPlaser_dec(ii+this_laser_jj-1+dt_next_pulse)>1)
                                        %This is a double pulse
                                        ii_laser=ii_laser+1;
                                        ii_out=(ii+this_laser_jj-1);
                                        handles.drgb.PAC.laser(no_trials).laser_times(ii_laser)=out_times_env(ii_out);
                                        handles.drgb.PAC.laser(no_trials).laser_phase(ii_laser)=decanglethetaLFP(ii_out);
                                        
                                    end
                                    
                                    
                                    ii=ii+this_laser_jj-1+dt_all_pulses;
                                    
                                else
                                    at_end=1;
                                end
                            else
                                at_end=1;
                            end
                        end
                    case 2
                        
                        if trNo>33
                            pffft=1;
                        end
                        
                        dt_down=10;
                        dt_next_pulse=27;
                        at_end=0;
                        ii=1;
                        handles.drgb.PAC.laser(no_trials).laser_times=[];
                        handles.drgb.PAC.laser(no_trials).laser_phase=[];
                        ii_laser=0;
                        handles.drgb.PAC.laser(no_trials).ii_laser=0;
                        thisLFPlaser=LFPlaser(pad_time*Fs:pad_time*Fs+out_times_env(end)*Fs-1);
                        Fs_out_times_env=1/(out_times_env(2)-out_times_env(1));
                        dec_factor=floor(Fs/Fs_out_times_env);
                        thisLFPlaser_dec=decimate(thisLFPlaser,dec_factor);
                         
                        if thisLFPlaser_dec(1)>1
                            this_laser_down_jj=find(thisLFPlaser_dec(ii:end)<1,1,'first');
                            ii=ii+this_laser_down_jj;
                        end
%                         figure(11)
%                         plot(thisLFPlaser_dec)
%                         figure(12)
%                         plot(decanglethetaLFP)
                        laser_times=[];
                        while at_end==0
                            this_laser_jj=find(thisLFPlaser_dec(ii:end)>1,1,'first');
                            if ~isempty(this_laser_jj)
                                
                                if (ii+this_laser_jj-1 +dt_next_pulse) < length(out_times_env)
                                    
                                   
                                    ii_laser=ii_laser+1;
                                    ii_out=(ii+this_laser_jj-1);
                                    handles.drgb.PAC.laser(no_trials).ii_laser=ii_laser;
                                    handles.drgb.PAC.laser(no_trials).laser_times(ii_laser)=out_times_env(ii_out)+handles.time_start+handles.time_pad;
                                    laser_times(ii_laser)=out_times_env(ii_out);
                                    handles.drgb.PAC.laser(no_trials).laser_phase(ii_laser)=decanglethetaLFP(ii_out);
     
                                    ii=ii+this_laser_jj-1+dt_next_pulse;
                                    
                                else
                                    at_end=1;
                                end
                            else
                                at_end=1;
                            end
                        end
                end
                
                pffft=1;
                

              
                
                
                  
                %Now do the wavelet transform
                handles.peakLowF=handles.LFPPowerSpectrumLowF;
                handles.peakHighF=handles.LFPPowerSpectrumHighF;
                handles.burstLowF=handles.LFPPowerSpectrumLowF;
                handles.burstHighF=handles.LFPPowerSpectrumHighF;
                
                [LFP, trialNo, can_read] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            
                handles.peakLowF=lowF1;
                handles.peakHighF=lowF2;
                handles.burstLowF=highF1;
                handles.burstHighF=highF2;
                
                decLFP=decimate(LFP,dec_n);
                decFs=Fs/dec_n;
                
                cwtLFP = cwtft({detrend(double(decLFP)),1/decFs},'wavelet','morl','scales',scales);
                Prev=abs(cwtLFP.cfs).^2;
                P=Prev(end:-1:1,:);
                DT=1/decFs;
                t = 0:DT:(numel(decLFP)*DT)-DT;
                frev=cwtLFP.frequencies;
                f=frev(end:-1:1);

                all_times=t+handles.time_start;
%                 out_times=all_times(:,(all_times>=handles.time_start+handles.time_pad)&(all_times<=handles.time_end-handles.time_pad));
                ii_from=find(all_times>=handles.time_start+handles.time_pad,1,'first');
                this_all_Power_timecourse(1:length(f),1:length(out_times_env))=P(:,ii_from:ii_from+length(out_times_env)-1);
                this_LFP(1,1:length(out_times_env))=decLFP(:,ii_from:ii_from+length(out_times_env)-1);
                all_Power_timecourse(no_trials,1:length(f),1:length(out_times_env))=this_all_Power_timecourse;
                handles.drgb.PAC.laser(no_trials).all_Power_timecourse=this_all_Power_timecourse;
                handles.drgb.PAC.laser(no_trials).out_times=out_times_env;
                if handles.subtractRef==1
                    P_ref=[];
                    P_ref=P(:,(all_times>=handles.startRef+handles.time_pad)&(all_times<=handles.endRef-handles.time_pad));
                    all_Power_ref(no_trials,1:length(f))=mean(P_ref,2);
                end
                
                %Now pull out the laser-referenced power windows
%                 dt_window=0.2; %Seconds for laser window
%                 dt_window_ii=ceil(dt_window/(out_times_env(2)-out_times_env(1)));
%                 half_dt_window_ii=ceil(dt_window_ii/2);
                dt_before=0.01;
                ii_dt_before=ceil(dt_before/(out_times_env(2)-out_times_env(1)));
                dt_after=0.02;
                ii_dt_after=ceil(dt_after/(out_times_env(2)-out_times_env(1)));
                dt_window_ii=ii_dt_before+ii_dt_after+1;
                ii_ref=2;
                
                for ii_l=1:length(laser_times)
                    [min_dt,min_ii]=min(abs(laser_times(ii_l)-out_times_env));
                    if (min_ii>ii_dt_before)&(min_ii+ii_dt_after<length(out_times_env))
                        length_snips=size(this_all_Power_timecourse(:,min_ii-ii_dt_before:min_ii+ii_dt_after),2);
                        handles.drgb.PAC.laser(no_trials).power_snips(ii_l).all_Power_timecourse_snip=this_all_Power_timecourse(:,min_ii-ii_dt_before:min_ii+ii_dt_after);
                        handles.drgb.PAC.laser(no_trials).power_snips(ii_l).LFP_snip=this_LFP(1,min_ii-ii_dt_before:min_ii+ii_dt_after);
                        handles.drgb.PAC.laser(no_trials).power_snips(ii_l).snip_calculated=1;
                    else
                        handles.drgb.PAC.laser(no_trials).power_snips(ii_l).all_Power_timecourse_snip=[];
                        handles.drgb.PAC.laser(no_trials).power_snips(ii_l).LFP_snip=[];
                        handles.drgb.PAC.laser(no_trials).power_snips(ii_l).snip_calculated=0;
                    end
                end
                
                
                
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

if ~isempty(spm)
    %Now figure out the gamma power at the peak
    if handles.displayData==0
        if handles.use_peakAngle==0
            %With no display (i.e., for batch) do this calculation only for S+
            %Note that I use phase, mean(enc_phase_histo), as opposed to
            %meanPeakAngle because meanPeakAngle may switch spuriously from
            %trial to trial when there are two peaks
            
            [max_hist max_ii]=max(mean(enc_phase_histo(logical(spm),:),1));
            meanPeakAngle=(phase(max_ii)*(pi/180))-pi;
            handles.drgb.PAC.peakAngleForPower=phase(max_ii);
            
            [min_hist min_ii]=min(mean(enc_phase_histo(logical(spm),:),1));
            meanTroughAngle=(phase(min_ii)*(pi/180))-pi;
            handles.drgb.PAC.troughAngleForPower=phase(min_ii);
        else
            meanPeakAngle=(handles.peakAngle_for_power*(pi/180))-pi;
            handles.drgb.PAC.peakAngleForPower=handles.peakAngle_for_power;
            meanTroughAngle=(handles.troughAngle_for_power*(pi/180))-pi;
            handles.drgb.PAC.troughAngleForPower=handles.troughAngle_for_power;
        end
    else
        if handles.use_peakAngle==0
            %Note that I use phase, mean(enc_phase_histo), as opposed to
            %meanPeakAngle because meanPeakAngle may switch spuriously from
            %trial to trial when there are two peaks
            [max_hist max_ii]=max(mean(enc_phase_histo));
            meanPeakAngle=(phase(max_ii)*(pi/180))-pi;
            handles.drgb.PAC.peakAngleForPower=phase(max_ii);
            if handles.displayData==1
                handles.peakAngle_for_power=phase(max_ii);
                set(handles.set_peakAngle,'String',num2str(phase(max_ii)));
            end
            [min_hist min_ii]=min(mean(enc_phase_histo));
            meanTroughAngle=(phase(min_ii)*(pi/180))-pi;
            handles.drgb.PAC.troughAngleForPower=phase(min_ii);
            if handles.displayData==1
                handles.troughAngle_for_power=phase(min_ii);
                set(handles.set_troughAngle,'String',num2str(phase(min_ii)));
            end
        else
            handles.peakAngle_for_power=str2num(get(handles.set_peakAngle,'String'));
            handles.troughAngle_for_power=str2num(get(handles.set_troughAngle,'String'));
            meanPeakAngle=(handles.peakAngle_for_power*(pi/180))-pi;
            handles.drgb.PAC.peakAngleForPower=handles.peakAngle_for_power;
            meanTroughAngle=(handles.troughAngle_for_power*(pi/180))-pi;
            handles.drgb.PAC.troughAngleForPower=handles.troughAngle_for_power;
        end
    end
    
    dt=handles.window-handles.noverlap;
    t=[handles.time_start+handles.time_pad:dt:handles.time_end-handles.time_pad];
    handles.drgb.PAC.t_power=t;
    peakPower=zeros(no_trials,length(t));
    %Find the value of gamma power at each point
    out_times_env=out_times_env+handles.time_start+handles.time_pad;
    
    %Find peak power
    for trNum=1:no_trials
        peakAmp=[];
        peakAmp_times=[];
        ii=1;
        jj=0;
        at_end=0;
        this_angleTetaLFP=handles.drgb.PAC.PACtimecourse(trNum).decanglethetaLFP;
        this_LFPenv=handles.drgb.PAC.PACtimecourse(trNum).decLFPgenv;
        while at_end==0
            ii_next=find(this_angleTetaLFP(ii:end)>=meanPeakAngle,1,'first');
            if (~isempty(ii_next))&(ii+ii_next-1<=length(out_times_env))
                jj=jj+1;
                peakAmp(jj)=this_LFPenv(ii+ii_next-1);
                peakAmp_times(jj)=out_times_env(ii+ii_next-1);
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
        
        for ii_t=1:length(t)
            if t(ii_t)<=peakAmp_times(1)
                peakPower(trNum,ii_t)=10*log10((peakAmp(1)^2)/2);
                %peakAmp_t(ii_t)=peakAmp(1);
            else
                if t(ii_t)>=peakAmp_times(end)
                    peakPower(trNum,ii_t)=10*log10((peakAmp(end)^2)/2);
                    %peakAmp_t(ii_t)=peakAmp(end);
                else
                    ii_pt=find(peakAmp_times>=t(ii_t),1,'first');
                    this_amp=peakAmp(ii_pt-1)+(t(ii_t)-peakAmp_times(ii_pt-1))*((peakAmp(ii_pt)-peakAmp(ii_pt-1))/(peakAmp_times(ii_pt)-peakAmp_times(ii_pt-1)));
                    %peakAmp_t(ii_t)=this_amp;
                    peakPower(trNum,ii_t)=10*log10((this_amp^2)/2);
                end
            end
        end
        
        if handles.subtractRef==1
            peakPower(trNum,:)=peakPower(trNum,:)-mean(peakPower(trNum,(t>=handles.startRef+handles.time_pad)&(t<=handles.endRef-handles.time_pad)));
        end
        
        handles.drgb.PAC.PACtimecourse(trNum).peakPower=peakPower(trNum,:);
        handles.drgb.PAC.meanPeakPower(trNum)=mean(peakPower(trNum,:));
    end
    
    %Find trough power
    troughPower=zeros(no_trials,length(t));
    for trNum=1:no_trials
        troughAmp=[];
        troughAmp_times=[];
        ii=1;
        jj=0;
        at_end=0;
        this_angleTetaLFP=handles.drgb.PAC.PACtimecourse(trNum).decanglethetaLFP;
        this_LFPenv=handles.drgb.PAC.PACtimecourse(trNum).decLFPgenv;
        while at_end==0
            ii_next=find(this_angleTetaLFP(ii:end)>=meanTroughAngle,1,'first');
            if (~isempty(ii_next))&(ii+ii_next-1<=length(out_times_env))
                jj=jj+1;
                troughAmp(jj)=this_LFPenv(ii+ii_next-1);
                troughAmp_times(jj)=out_times_env(ii+ii_next-1);
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
        
        for ii_t=1:length(t)
            if t(ii_t)<=troughAmp_times(1)
                troughPower(trNum,ii_t)=10*log10((troughAmp(1)^2)/2);
                %troughAmp_t(ii_t)=troughAmp(1);
            else
                if t(ii_t)>=troughAmp_times(end)
                    troughPower(trNum,ii_t)=10*log10((troughAmp(end)^2)/2);
                    %troughAmp_t(ii_t)=troughAmp(end);
                else
                    ii_pt=find(troughAmp_times>=t(ii_t),1,'first');
                    this_amp=troughAmp(ii_pt-1)+(t(ii_t)-troughAmp_times(ii_pt-1))*((troughAmp(ii_pt)-troughAmp(ii_pt-1))/(troughAmp_times(ii_pt)-troughAmp_times(ii_pt-1)));
                    %troughAmp_t(ii_t)=this_amp;
                    troughPower(trNum,ii_t)=10*log10((this_amp^2)/2);
                end
            end
        end
        
        if handles.subtractRef==1
            troughPower(trNum,:)=troughPower(trNum,:)-mean(troughPower(trNum,(t>=handles.startRef+handles.time_pad)&(t<=handles.endRef-handles.time_pad)));
        end
        
        handles.drgb.PAC.PACtimecourse(trNum).troughPower=troughPower(trNum,:);
        handles.drgb.PAC.meanTroughPower(trNum)=mean(troughPower(trNum,:));
    end
    
    
    if handles.displayData==1
        
        all_out_time_PAChisto=all_out_time_PAChisto/no_trials;
        mean_MI_enc=mean(MI_enc);
        
        %Plot the timecourse for the phase
        try
            close 1
        catch
        end
        
        hFig1 = figure(1);
        set(hFig1, 'units','normalized','position',[.15 .6 .8 .23])
        
        %     plot(out_times,circ_rad2ang(circ_mean(circ_ang2rad(all_out_phase)-180))+180,'ob')
        plot(out_times,(180/pi)*circ_axial(circ_mean(all_out_phase*pi/180)'),'ob')
        ylim([0 360])
        
        
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
        
        
        try
            close 7
        catch
        end
        
        hFig7 = figure(7);
        set(hFig7, 'units','normalized','position',[.25 .2 .7 .23])
        hold on
        
        if no_trials>2
            [hltrough, hptrough]=boundedline(t,mean_troughPower, CI_troughPower, 'b');
            [hlpeak, hppeak]=boundedline(t,mean_peakPower, CI_peakPower, 'r');
        else
            plot(t,mean_troughPower, 'b')
            plot(t,mean_peakPower, 'r')
        end
        xlabel('Time(sec)')
        ylabel('Peak power (dB)')
        
        if no_trials>2
            legend([hltrough hlpeak],{'Trough','Peak'})
        else
            legend('Through','Peak')
        end
        title(['Power (dB, PAC high frequency envelope)  ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
        
        
        if handles.autoscale==0
            ylim([handles.minLogP handles.maxLogP]);
        end
        
        try
            close 6
        catch
        end
        
        hFig6 = figure(6);
        set(hFig6, 'units','normalized','position',[.25 .2 .7 .23])
        
        min_prob=prctile(all_out_time_PAChisto(:),5);
        max_prob=prctile(all_out_time_PAChisto(:),95);
        
        
        %pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),all_phase_histo)
        drg_pcolor(repmat(out_times,length(phase(1:end-1)),1),repmat(phase(1:end-1)',1,length(out_times)),all_out_time_PAChisto')
%         colormap jet
        colormap fire
        shading flat
        % min_prob=0.0113;
        % max_prob=0.0314;
        caxis([min_prob    max_prob])
        xlabel('Time (sec)')
        ylabel('Phase for low freq oscillation (deg)');
        title(['Phase-amplitude coupling timecourse ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
        
        %     (180/pi)*circ_mean(all_out_phase'*pi/180)'
        try
            close 8
        catch
        end
        
        hFig8=figure(8);
        hold on
        
        %Naive
        if sum((handles.drgb.PAC.perCorr>=45)&(handles.drgb.PAC.perCorr<=65))>2
            [f_MI,x_MI] = drg_ecdf(mod_indx((handles.drgb.PAC.perCorr>=45)&(handles.drgb.PAC.perCorr<=65)));
            plot(x_MI,f_MI,'-b','LineWidth',3);
        end
        
        %Proficient
        if sum(handles.drgb.PAC.perCorr>=80)>2
            [f_MI,x_MI] = drg_ecdf(mod_indx(handles.drgb.PAC.perCorr>=80));
            plot(x_MI,f_MI,'-r','LineWidth',3);
        end
        
        xlabel('MI')
        ylabel('Probability')
        title('Cumulative histogram for MI')
        legend('Naive','Proficient')
        
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
        hold on
        if length(mod_indx)>20
            no_conv_points=6;
            conv_win=ones(1,no_conv_points);
            mod_indx_extend=[mean(mod_indx(1:no_conv_points/2))*ones(1,no_conv_points) mod_indx mean(mod_indx(end-(no_conv_points/2)+1:end))*ones(1,no_conv_points)];
            conv_mod_indx_extend=conv(mod_indx_extend,conv_win)/no_conv_points;
            conv_mod_indx=conv_mod_indx_extend(no_conv_points+1:no_conv_points+length(mod_indx));
            plot(trials,conv_mod_indx,'-b','LineWidth',3)
        end
        xlabel('Trial No')
        ylabel('Modulation index')
        title('Modulation index vs trial number')
        xlim([0 no_trials])
        
        
        %Plot the mean vector length as a function of trial number
        subplot(4,1,2)
        trials=1:no_trials;
        plot(trials,meanVectorLength,'o-')
        hold on
        if length(meanVectorLength)>20
            no_conv_points=6;
            conv_win=ones(1,no_conv_points);
            meanVectorLength_extend=[mean(meanVectorLength(1:no_conv_points/2))*ones(1,no_conv_points) meanVectorLength mean(meanVectorLength(end-(no_conv_points/2)+1:end))*ones(1,no_conv_points)];
            conv_meanVectorLength_extend=conv(meanVectorLength_extend,conv_win)/no_conv_points;
            conv_meanVectorLength=conv_meanVectorLength_extend(no_conv_points+1:no_conv_points+length(meanVectorLength));
            plot(trials,conv_meanVectorLength,'-b','LineWidth',3)
        end
        xlabel('Trial No')
        ylabel('Mean vector length')
        title('Mean vector length vs trial number')
        xlim([0 no_trials]);
        
        % figure(11)
        % plot(trials,mod_indx,'o-k','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7,'LineWidth',2)
        % xlabel('Trial No')
        % ylabel('Modulation index')
        % title('Modulation index vs trial number')
        % ylim([0 0.07])
        % set(gca,'FontSize',20,'LineWidth',1.25)
        
        
        
        %Now plot the mean theta waveform
        subplot(4,1,3)
        shadedErrorBar(phase,mean(all_theta_wave,1),std(all_theta_wave,0,1),'-b')
        xlim([0 360])
        title('Mean low frequency waveform')
        xlabel('Degrees')
        
        %Now plot the encoding theta phase histogram for gamma
        subplot(4,1,4)
        shadedErrorBar(phase,mean(enc_phase_histo,1),std(enc_phase_histo,0,1),'-b')
        xlim([0 360])
        title('Phase-amplitude coupling')
        ylabel('Probability')
        xlabel('Phase for low freq oscillation (deg)')
        yl_enc=ylim;
        text(20,yl_enc(1)+0.9*(yl_enc(2)-yl_enc(1)),['MI= ' num2str(mean_MI_enc)])
        
        
        if no_trials>2
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
%             colormap jet
            colormap fire
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
%             colormap jet
            colormap fire
            shading interp
            ax=gca;
            set(ax,'XTickLabel','')
        end
        
        try
            close 5
        catch
        end
        
        hFig5 = figure(5);
        set(hFig5, 'units','normalized','position',[.59 .05 .35 .35])
        
        if no_encoding_trials>0
            polarhistogram(pi*peakAngle/180,12)
            title('Peak angle')
        end
        
        all_laser_phase=[];
        laser_pulse_no=[];
        for trNo=1:no_trials
            all_laser_phase=[all_laser_phase handles.drgb.PAC.laser(trNo).laser_phase];
            laser_pulse_no(trNo)=handles.drgb.PAC.laser(trNo).ii_laser;
        end
        
        
         try
            close 9
        catch
        end
        
        hFig9 = figure(9);
        ax=gca;ax.LineWidth=3;
        set(hFig9, 'units','normalized','position',[.59 .05 .35 .35])
        
        if length(all_laser_phase)>0
            polarhistogram(all_laser_phase+pi,12)
            title('Laser phase')
        end
        
         try
            close 10
        catch
        end
        
        hFig10 = figure(10);
        ax=gca;ax.LineWidth=3;
        set(hFig10, 'units','normalized','position',[.22 .05 .35 .35])
        
        if length(all_laser_phase)>0
            plot((180*all_laser_phase/pi)+180,'.k')
            title('Laser phase')
            ylabel('phase')
        end
        
        laser_phase=(180*all_laser_phase/pi)+180;
        save([handles.fullName(1:end-7) 'laser.mat'],'laser_phase')
         
         try
            close 11
        catch
        end
        
        hFig11 = figure(11);
        set(hFig11, 'units','normalized','position',[.0 .05 .35 .35])
        
        if length(laser_pulse_no)>0
            plot(laser_pulse_no,'ok')
            title('Number of laser pulses')
            ylabel('No')
        end
        
        %Show a pseudocolor of the average power timecourse
        
        %Timecourse doing average after log
        freq=f;
        t_apt=out_times_env;
        if handles.subtractRef==0
            log_P_timecourse=zeros(length(freq),length(t_apt));
            log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),1);
            
            %Per trial power plot
            log_P_per_trial_timecourse=zeros(length(freq)*trialNo,length(t_apt));
            y_shift=0;
            log_P_snip_per_trial_timecourse_sub=[];
            ii_snip_sub=0;
            
            for trialNo=1:length(no_trials)
                this_log_P_timecourse=zeros(length(freq),length(t_apt));
                this_log_P_timecourse(:,:)=10*log10(all_Power_timecourse(trialNo,:,:));
                log_P_per_trial_timecourse(y_shift+1:y_shift+length(freq),:)=this_log_P_timecourse;
                shifted_freq(1,y_shift+1:y_shift+length(freq))=freq+(trialNo-1)*freq(end);
                y_shift=y_shift+length(freq);
                
                 for ii_snips=1:length(handles.drgb.PAC.laser(trialNo).power_snips)
                    if (handles.drgb.PAC.laser(trialNo).power_snips(ii_snips).snip_calculated==1)
                        this_Power_timecourse_snip=handles.drgb.PAC.laser(trialNo).power_snips(ii_snips).all_Power_timecourse_snip;
                        this_log_P_snip_timecourse(:,:)=10*log10(this_Power_timecourse_snip);
                        ii_snip_sub=ii_snip_sub+1;
                        log_P_snip_per_trial_timecourse_sub(ii_snip_sub,:,:)=this_log_P_snip_timecourse;
                    end
                end
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
            log_P_per_trial_timecourse_sub=zeros(length(freq)*length(no_trials),length(t_apt));
            log_P_snip_per_trial_timecourse_sub=[];
            LFP_snip_per_trial_timecourse=[];
            logP_trials=[];
            ii_snip_sub=0;
            
            y_shift=0;
            sy_shift=0;
            shifted_freq=[];
            for trialNo=1:no_trials
                this_log_P_timecourse=zeros(length(freq),length(t_apt));
                this_log_P_timecourse(:,:)=10*log10(all_Power_timecourse(trialNo,:,:));
                this_log_P_timecourse_ref=zeros(length(freq),length(t_apt));
                this_log_P_timecourse_ref(:,:)=repmat(mean(10*log10(all_Power_ref(trialNo,:)),1)',1,length(t_apt));
                this_log_P_timecourse_ref_snip=[];
                this_log_P_timecourse_ref_snip(:,:)=repmat(mean(10*log10(all_Power_ref(trialNo,:)),1)',1,dt_window_ii);
                log_P_per_trial_timecourse_sub(y_shift+1:y_shift+length(freq),:)=this_log_P_timecourse-this_log_P_timecourse_ref;
                shifted_freq(1,y_shift+1:y_shift+length(freq))=freq+(trialNo-1)*freq(end);
                y_shift=y_shift+length(freq);
                
                for ii_snips=1:length(handles.drgb.PAC.laser(trialNo).power_snips)
                    if handles.drgb.PAC.laser(trialNo).power_snips(ii_snips).snip_calculated==1
                        this_Power_timecourse_snip=handles.drgb.PAC.laser(trialNo).power_snips(ii_snips).all_Power_timecourse_snip;
                        this_log_P_snip_timecourse(:,:)=10*log10(this_Power_timecourse_snip);
                        ii_snip_sub=ii_snip_sub+1;
                        log_P_snip_per_trial_timecourse_sub(ii_snip_sub,:,:)=this_log_P_snip_timecourse-this_log_P_timecourse_ref_snip;
                        LFP_snip_per_trial_timecourse(ii_snip_sub,:)=handles.drgb.PAC.laser(trialNo).power_snips(ii_snips).LFP_snip;
                        logP_trials(ii_snip_sub)=trialNo;
                    end
                end
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
        
        %Plot the average wavelet power timecourse
        figNo=12;
        
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
        
        
        
        %Show a pseudocolor of average the laser snips
        figNo=13;
        
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        avg_log_P_snip_per_trial_timecourse_sub=zeros(size(log_P_snip_per_trial_timecourse_sub,2),size(log_P_snip_per_trial_timecourse_sub,3));
        avg_log_P_snip_per_trial_timecourse_sub(:,:)=mean(log_P_snip_per_trial_timecourse_sub,1);
        t_snips=[0:out_times_env(2)-out_times_env(1):dt_before+dt_after]-dt_before;
        t_snips=[t_snips t_snips(end)+out_times_env(2)-out_times_env(1)];
        
        %Plot the timecourse
        hFig = figure(figNo);
        set(hFig, 'units','normalized','position',[.07 .05 .3 .3])
         
        %pcolor(repmat(t,length(f),1)',repmat(f,length(t),1),10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
        time_length=size(avg_log_P_snip_per_trial_timecourse_sub,2);
        if time_length<length(t_snips)
            t_snips=t_snips(1:time_length);
        end
        drg_pcolor(repmat(t_snips,length(freq),1)',repmat(freq,length(t_snips),1),avg_log_P_snip_per_trial_timecourse_sub')
        %imagesc(t,f,10*log10(P_timecourse')-10*log10(P_timecourse_ref'))
        
        
        colormap fire
        shading interp
        caxis([minLogPper maxLogPper]);
        xlabel('Time (sec)')
        ylabel('Frequency (Hz)');
        title(['Laser-referenced power (dB, wavelet) timecourse ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
        
        %Plot the laser-referenced power for each bandwidth
         figNo=14;
         
         figNo=figNo+1;
         try
             close(figNo)
         catch
         end
         
         hFig = figure(figNo);
         set(hFig, 'units','normalized','position',[.07 .05 .4 .4])
         
      
          
         for ii_bw=1:length(handles.lowF)
             subplot(2,2,ii_bw)
             hold on
             this_bw_log_P_snip_per_trial_timecourse_sub=zeros(size(log_P_snip_per_trial_timecourse_sub,1),size(log_P_snip_per_trial_timecourse_sub,3));
             this_bw_log_P_snip_per_trial_timecourse_sub(:,:)=mean(log_P_snip_per_trial_timecourse_sub(:,(f>=handles.lowF(ii_bw))&(f>=handles.highF(ii_bw)),:),2);
             mean_this_bw_log_P_snip_per_trial_timecourse_sub=mean(this_bw_log_P_snip_per_trial_timecourse_sub(:,ii_dt_before-ii_ref:ii_dt_before+ii_ref),2);
             this_bw_log_P_snip_per_trial_timecourse_sub(:,:)=this_bw_log_P_snip_per_trial_timecourse_sub(:,:)-repmat(mean_this_bw_log_P_snip_per_trial_timecourse_sub,1,size(this_bw_log_P_snip_per_trial_timecourse_sub,2));
             CIsp = bootci(1000, @mean, this_bw_log_P_snip_per_trial_timecourse_sub);
             meansp=mean(this_bw_log_P_snip_per_trial_timecourse_sub,1);
             CIsp(1,:)=meansp-CIsp(1,:);
             CIsp(2,:)=CIsp(2,:)-meansp;

             [hlsp, hpsp] = boundedline(t_snips',mean(this_bw_log_P_snip_per_trial_timecourse_sub,1)', CIsp', 'cmap',[80/255 194/255 255/255]);

             this_yl=ylim;

             plot([0 0],this_yl)

             %              plot(t_snips',this_bw_log_P_snip_per_trial_timecourse_sub','-k')

             title(handles.f_label{ii_bw})
             xlabel('Time (sec)')
             ylabel('dB')
         end

         %Plot the laser-referenced LFP
         figNo=15;

         figNo=figNo+1;
         try
             close(figNo)
         catch
         end

         hFig = figure(figNo);
         set(hFig, 'units','normalized','position',[.07 .05 .4 .4])



         mean_LFP_snip_per_trial_timecourse=mean(LFP_snip_per_trial_timecourse(:,ii_dt_before-ii_ref:ii_dt_before+ii_ref),2);
         LFP_snip_per_trial_timecourse(:,:)=LFP_snip_per_trial_timecourse(:,:)-repmat(mean_LFP_snip_per_trial_timecourse,1,size(this_bw_log_P_snip_per_trial_timecourse_sub,2));
            

         hold on
         
         CIsp = bootci(1000, @mean, LFP_snip_per_trial_timecourse);
         meansp=mean(LFP_snip_per_trial_timecourse,1);
         CIsp(1,:)=meansp-CIsp(1,:);
         CIsp(2,:)=CIsp(2,:)-meansp;

         [hlsp, hpsp] = boundedline(t_snips',mean(LFP_snip_per_trial_timecourse,1)', CIsp', 'cmap',[80/255 194/255 255/255]);

         this_yl=ylim;

         plot([0 0],this_yl)

         %              plot(t_snips',this_bw_log_P_snip_per_trial_timecourse_sub','-k')

         title('Laser-referenced LFP')
         xlabel('Time (sec)')
         ylabel('uV')


         try
             close 17
         catch
         end

         hFig2 = figure(17);
         ax=gca;ax.LineWidth=3;
         set(hFig10, 'units','normalized','position',[.22 .05 .35 .35])

         %Now plot the mean theta waveform
         shadedErrorBar(phase,mean(all_theta_wave,1),std(all_theta_wave,0,1),'-b')
         xlim([0 360])
         title('Mean low frequency waveform')
         xlabel('Degrees')

%          %Plot the laser-referenced power for one bandwidth for the grant
%          figNo=17;
% 
%          figNo=figNo+1;
%          try
%              close(figNo)
%          catch
%          end
% 
%          hFig = figure(figNo);
%          ax=gca;ax.LineWidth=3;
%          set(hFig10, 'units','normalized','position',[.22 .05 .35 .35])
% 
%          hold on
% 
%          ii_bw=4;
% 
%          hold on
%          this_bw_log_P_snip_per_trial_timecourse_sub=zeros(size(log_P_snip_per_trial_timecourse_sub,1),size(log_P_snip_per_trial_timecourse_sub,3));
%          this_bw_log_P_snip_per_trial_timecourse_sub(:,:)=mean(log_P_snip_per_trial_timecourse_sub(:,(f>=handles.lowF(ii_bw))&(f>=handles.highF(ii_bw)),:),2);
%          mean_this_bw_log_P_snip_per_trial_timecourse_sub=mean(this_bw_log_P_snip_per_trial_timecourse_sub(:,ii_dt_before-ii_ref:ii_dt_before+ii_ref),2);
%          this_bw_log_P_snip_per_trial_timecourse_sub(:,:)=this_bw_log_P_snip_per_trial_timecourse_sub(:,:)-repmat(mean_this_bw_log_P_snip_per_trial_timecourse_sub,1,size(this_bw_log_P_snip_per_trial_timecourse_sub,2));
% 
%          %Plot before laser
%          pre_bw_log_P_snip_per_trial_timecourse_sub=this_bw_log_P_snip_per_trial_timecourse_sub((logP_trials>=20)&(logP_trials<=40),:);
%          CIsp = bootci(1000, @mean, pre_bw_log_P_snip_per_trial_timecourse_sub);
%          meansp=mean(pre_bw_log_P_snip_per_trial_timecourse_sub,1);
%          CIsp(1,:)=meansp-CIsp(1,:);
%          CIsp(2,:)=CIsp(2,:)-meansp;
% 
%          [hlsp, hpsp] = boundedline(t_snips',mean(pre_bw_log_P_snip_per_trial_timecourse_sub,1)', CIsp', 'cmap',[0/255 0/255 0/255]);
% 
%          %Plot after laser
%          post_bw_log_P_snip_per_trial_timecourse_sub=this_bw_log_P_snip_per_trial_timecourse_sub((logP_trials>=42)&(logP_trials<=80),:);
%          CIsp = bootci(1000, @mean, post_bw_log_P_snip_per_trial_timecourse_sub);
%          meansp=mean(post_bw_log_P_snip_per_trial_timecourse_sub,1);
%          CIsp(1,:)=meansp-CIsp(1,:);
%          CIsp(2,:)=CIsp(2,:)-meansp;
% 
%          [hlsp, hpsp] = boundedline(t_snips',mean(post_bw_log_P_snip_per_trial_timecourse_sub,1)', CIsp', 'cmap',[80/255 194/255 255/255]);
% 
%          this_yl=ylim;
% 
%          plot([0 0],this_yl)
% 
%          %              plot(t_snips',this_bw_log_P_snip_per_trial_timecourse_sub','-k')
% 
%          title(handles.f_label{ii_bw})
%          xlabel('Time (sec)')
%          ylabel('dB')
%  

         pffft=1;
    end
end

if handles.displayData==1
    fprintf(1, ['\nPCA processed for %d out of %d trials \n\n'], no_trials,trials_attempted);
end

pfffft=1;






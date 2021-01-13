function handles=drgSpikePhase(handles)

%Generates a phase histogram
sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
%window=round(handles.window*handles.drg.draq_p.ActualRate); 


%Time shift the LFP if nescessary
if isfield(handles,'time_shift')
    time_shift=handles.time_shift;
else
    time_shift=0;
end

%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;

%Butterworth filter for spike-triggered LFP
fpass=[handles.burstLowF handles.burstHighF];

bpFilt = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',fpass(1),'HalfPowerFrequency2',fpass(2), ...
    'SampleRate',handles.drg.session(sessionNo).draq_p.ActualRate);
 
%Now get the phase of gamma bursts witin theta
no_trials=0;
valid_events=[];
no_valid=0;
no_spikes=0;

no_bins=handles.n_phase_bins/2;
theta_wave=zeros(1,no_bins+1);
n_wave=zeros(1,no_bins+1);
phase_histo=zeros(1,no_bins+1);
phase=[0:360/no_bins:360];

ii_pad=ceil(handles.time_pad*Fs);

for trNo=firstTr:lastTr
    
    if handles.save_drgb==0
        trial_no=trNo;
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
     
    if evNo~=-1
        
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
 
            
            [LFP, trialNo, can_read] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start+time_shift, handles.time_end+time_shift);
             
            
            if (can_read==1)
                
                this_event_time=handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo);
                
                ffLFP=filtfilt(bpFilt,LFP);
                angle_ffLFP = angle(hilbert(ffLFP)); % phase
                
                
                for ii=ii_pad:length(LFP)-ii_pad
                    n_wave(ceil((angle_ffLFP(ii)+pi)/(2*pi/no_bins)))=n_wave(ceil((angle_ffLFP(ii)+pi)/(2*pi/no_bins)))+1;
                    theta_wave(ceil((angle_ffLFP(ii)+pi)/(2*pi/no_bins)),n_wave(ceil((angle_ffLFP(ii)+pi)/(2*pi/no_bins))))=ffLFP(ii);
                end
                

                %Get the spikes
                
                
                ii_spikes=find((handles.drg.unit(handles.unitNo).spike_times>=this_event_time+handles.time_start+handles.time_pad)&...
                    (handles.drg.unit(handles.unitNo).spike_times<=this_event_time+handles.time_end-handles.time_pad));
                
              
                
                no_trials=no_trials+1;
                
                if length(ii_spikes>0)
                    no_valid=no_valid+1;
                    valid_events(no_valid)=evNo;
                    no_ii_spikes(no_valid)=length(ii_spikes);

                    for kk=1:length(ii_spikes)
                        no_spikes=no_spikes+1;
                        time_this_spike=handles.drg.unit(handles.unitNo).spike_times(ii_spikes(kk))-handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo);
                        ii_this_spike=ceil(1+(time_this_spike-handles.time_start)*handles.drg.session(sessionNo).draq_p.ActualRate);
                        phase_histo(ceil((angle_ffLFP(ii_this_spike)+pi)/(2*pi/no_bins)))=phase_histo(ceil((angle_ffLFP(ii_this_spike)+pi)/(2*pi/no_bins)))+1;
                    end
  
                end
                
            end
        end
    end
end

%Calculate the modulation index
norm_phase_histo=phase_histo/sum(phase_histo);
handles.sp_phase=phase;

mean_prob=mean(norm_phase_histo)*ones(1,length(norm_phase_histo));
DKLii=norm_phase_histo(1:end-1).*log(norm_phase_histo(1:end-1)./mean_prob(1:end-1));
DKL=sum(DKLii(~isnan(DKLii)));
% DKL=sum(norm_phase_histo(1:end-1).*log(norm_phase_histo(1:end-1)./mean_prob(1:end-1)));
MI_Tort=DKL/log(no_bins);

%Peak angle
[max_hist peak_bin]=max(phase_histo);
peakAngle=phase(peak_bin);

%Mean vector length and vector angle
phaseAngle=pi*phase/180;
meanVectorLength=sqrt((mean(norm_phase_histo.*sin(phaseAngle)))^2+(mean(norm_phase_histo.*cos(phaseAngle)))^2);
meanX=mean(norm_phase_histo.*cos(phaseAngle));
meanY=mean(norm_phase_histo.*sin(phaseAngle));
if meanY>0
    meanVectorAngle=(180/pi)*acos(meanX/sqrt(meanY^2+meanX^2));
else
    meanVectorAngle=360-(180/pi)*acos(meanX/sqrt(meanY^2+meanX^2));
end

%MLR taken from Adhikari et al. Neuron 65:257, 2010
phase_rad=circ_ang2rad(phase(1:end-1));
d=phase_rad(2)-phase_rad(1);
[pval z] = circ_rtest(phase_rad, phase_histo(1:end-1), d);
MRL=sqrt(z/sum(phase_histo));

if handles.displayData==1
    fprintf(1, 'Bandpass filter from %d to %d Hz\n', handles.burstLowF,handles.burstHighF);
    fprintf(1, 'Time shift (sec)= %d\n', time_shift);
    fprintf(1, 'MI Tort= %d\n', MI_Tort);
    fprintf(1, 'Peak angle= %d\n', peakAngle);
    fprintf(1, 'Mean vector angle= %d\n', meanVectorAngle);
    fprintf(1, 'meanVectorLength= %d\n', meanVectorLength);
    fprintf(1, 'MRL= %d\n', MRL);
    fprintf(1, 'Number of spikes= %d\n', sum(phase_histo));
end


phase_histo(no_bins+1)=phase_histo(1);
norm_phase_histo(no_bins+1)=norm_phase_histo(1);
handles.norm_phase_histo=norm_phase_histo;
for ii=1:no_bins
    mean_theta_wave(ii)=mean(theta_wave(ii,1:n_wave(ii)),2);
    std_theta_wave(ii)=std(theta_wave(ii,1:n_wave(ii)),0,2);
end
mean_theta_wave(no_bins+1)=mean_theta_wave(1);
std_theta_wave(no_bins+1)=std_theta_wave(1);

%Output structure
handles.drgb.SpikePhase.no_trials=no_trials;
handles.drgb.SpikePhase.MI_Tort=MI_Tort;
handles.drgb.SpikePhase.peakAngle=peakAngle;
handles.drgb.SpikePhase.meanVectorAngle=meanVectorAngle;
handles.drgb.SpikePhase.meanVectorLength=meanVectorLength;
handles.drgb.SpikePhase.MRL=MRL;
handles.drgb.SpikePhase.noSpikes=sum(phase_histo);



if handles.displayData==1
    try
        close 1
    catch
    end
    
    %Plot the histogram
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.02 .4 .5 .5])
    
    subplot(2,1,2)
    bar(phase,norm_phase_histo)
    xlim([0 360])
    xlabel('degrees')
    ylabel('Probability')
    title('Spike phase')
    
    
    subplot(2,1,1)
    shadedErrorBar(phase,mean_theta_wave,std_theta_wave,'-b')
    xlim([0 360])
    xlabel('degrees')
    ylabel('uv')
    title('LFP')
end



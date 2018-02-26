function [meanVectorLength, meanVectorAngle, peakAngle, MI_Tort, phase, phase_histo, theta_wave]=drgGetThetaAmpPhase(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,time_pad,no_bins,method)
%Generates the phase histogram for the emvelope and pac
%function [pac_value, mod_indx, phase, phase_histo, theta_wave]=drgGetThetaAmpPhase(LFP,Fs,lowF1,lowF2,highF1,highF2,time_pad,no_bins)


%Time pad is used to exclude filter artifacts at the end
Fs=floor(Fs);

switch method
    case 1
        %Butterworth
        ii_pad=ceil(time_pad*Fs);
        
        %Filter gamma
        bpFiltgamma = designfilt('bandpassiir','FilterOrder',20, ...
            'HalfPowerFrequency1',highF1,'HalfPowerFrequency2',highF2, ...
            'SampleRate',Fs);
        filtLFPgamma=filtfilt(bpFiltgamma,LFPhigh);
    case 2
        %Tort eegfilt
        
        %First resample data for fast eegfilt
        rsFs=1024;
        LFPlowrs=resample(LFPlow,rsFs,Fs);    %Note: It is assumed that the high rate will not be larger than 250 Hz
        
        %Time pad is used to exclude filter artifacts at the end
        ii_pad=ceil(time_pad*rsFs);
        
        %Filter
        filtLFPgamma = eegfilt(LFPlowrs, rsFs, highF1, highF2);
end

LFPgenv = abs(hilbert(filtLFPgamma)); % gamma envelope
angleGammaEnv = angle(hilbert(filtLFPgamma)); % phase modulation of amplitude

% Now filter the theta
switch method
    case 1
        bpFilttheta = designfilt('bandpassiir','FilterOrder',20, ...
            'HalfPowerFrequency1',lowF1,'HalfPowerFrequency2',lowF2, ...
            'SampleRate',Fs);
        thfiltLFPgenv=filtfilt(bpFilttheta,LFPgenv);
    case 2
        thfiltLFPgenv = eegfilt(LFPgenv, rsFs, lowF1, lowF2);
end
   
angleThFlGammaEnv = angle(hilbert(thfiltLFPgenv)); % phase modulation of amplitude
% vect_length_Env = abs(hilbert(thfiltLFPgenv)); % vector length of the envelope filtered at theta
% meanVectorLength=sqrt((mean(vect_length_Env.*sin(angleThFlGammaEnv)))^2+(mean(vect_length_Env.*cos(angleThFlGammaEnv)))^2)/mean(LFPgenv);




%Filter LFP theta
switch method
    case 1
        thfiltLFP=filtfilt(bpFilttheta,LFPlow);
    case 2
        thfiltLFP = eegfilt(LFPlowrs, rsFs, lowF1, lowF2);
end

%thfiltLFP=LFPlow;

anglethetaLFP = angle(hilbert(thfiltLFP)); % phase modulation of amplitude

%This is PAC as defined by Bradley Voytek
pac_value = abs(sum(exp(1i * (anglethetaLFP - angleThFlGammaEnv)), 'double')) / length(LFPlow);






%Now generate the histogram
phase_histo=zeros(1,no_bins);
theta_wave=zeros(1,no_bins);
envelope_wave=zeros(1,no_bins);
amplitude_wave=zeros(1,no_bins);
phase=[0:360/no_bins:360];
sum_env=sum(LFPgenv(ii_pad:end-ii_pad));
for ii=ii_pad:length(LFPgenv)-ii_pad
    phase_histo(ceil((anglethetaLFP(ii)+pi)/(2*pi/no_bins)))=phase_histo(ceil((anglethetaLFP(ii)+pi)/(2*pi/no_bins)))+LFPgenv(ii); 
    theta_wave(ceil((anglethetaLFP(ii)+pi)/(2*pi/no_bins)))=theta_wave(ceil((anglethetaLFP(ii)+pi)/(2*pi/no_bins)))+thfiltLFP(ii);
    %envelope_wave(ceil((anglethetaLFP(ii)+pi)/(2*pi/no_bins)))=envelope_wave(ceil((anglethetaLFP(ii)+pi)/(2*pi/no_bins)))+thfiltLFPgenv(ii);
end
phase_histo(no_bins+1)=phase_histo(1);
theta_wave(no_bins+1)=theta_wave(1);
phase_histo=phase_histo/sum_env;
[max_hist peak_bin]=max(phase_histo);
peakAngle=phase(peak_bin);

%Calculate the modulation index defined by Tort et al J Neurophysiol 104: 1195?1210, 2010
%Note that the pvalue for Tort et al is the same as phase_histo
mean_prob=mean(phase_histo)*ones(1,length(phase_histo));
DKL=sum(phase_histo(1:end-1).*log(phase_histo(1:end-1)./mean_prob(1:end-1)));
MI_Tort=DKL/log(no_bins);

%This creates a polar plot of the gamma LFP with the theta angle


% figure(11)
% polar(pi*phase/180,phase_histo)

phaseAngle=pi*phase/180;
meanVectorLength=sqrt((mean(phase_histo.*sin(phaseAngle)))^2+(mean(phase_histo.*cos(phaseAngle)))^2);
meanX=mean(phase_histo.*cos(phaseAngle));
meanY=mean(phase_histo.*sin(phaseAngle));
if meanY>0
    meanVectorAngle=(180/pi)*acos(meanX/sqrt(meanY^2+meanX^2));
else
    meanVectorAngle=360-(180/pi)*acos(meanX/sqrt(meanY^2+meanX^2));
end


% 
% % %Plot for troubleshooting
% figure(10)
% subplot(3,1,1)
% plot(LFPgenv)
% title('Envelope for high frequency LFP')
% 
% subplot(3,1,2)
% plot(anglethetaLFP)
% title('Angle for low frequency LFP/sniff')
% 
% subplot(3,1,3)
% plot(LFPlow)
% title('Low frequency LFP/sniff')
% 
% pffft=1;





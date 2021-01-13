function [rho, t_lag]=drgGetLFPcorr(LFP1,LFPshifted,Fs,F1,F2,time_pad)
%Generates the phase histogram for the envelope and pac
%function [pac_value, mod_indx, phase, phase_histo, theta_wave]=drgGetThetaAmpPhase(LFP,Fs,lowF1,lowF2,highF1,highF2,time_pad,no_bins)
 
 
%Time pad is used to exclude filter artifacts at the end
Fs=floor(Fs);
 

%Butterworth
ii_pad=ceil(time_pad*Fs);

%Filter gamma
bpFilt = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',F1,'HalfPowerFrequency2',F2, ...
    'SampleRate',Fs);

filtLFP1=filtfilt(bpFilt,LFP1);
LFP1env = abs(hilbert(filtLFP1)); % envelope
powerLFP1filt=LFP1env.^2;
angleLFP1 = angle(hilbert(filtLFP1)); % phase modulation of amplitude

filtLFPshifted=filtfilt(bpFilt,LFPshifted);
LFPshiftedenv = abs(hilbert(filtLFPshifted)); % envelope
powerLFPshiftedfilt=LFPshiftedenv.^2;
angleLFPshifted = angle(hilbert(filtLFPshifted)); % phase modulation of amplitude

max_t_shift=(1/mean([F1,F2]))/2;
max_ii_shift=max_t_shift*Fs;

n_bins=41;
t_lag=[];
rho=[];
for ii_shift=1:n_bins
    ii=ceil(-max_ii_shift+(ii_shift-1)*(max_ii_shift/((n_bins-1)/2)));
    t_lag(ii_shift)=ii/Fs;
    if ii<0
        powerLFPshiftedshifted=powerLFPshiftedfilt(-ii:end);
        powerLFP1shifted=powerLFP1filt(1:end+ii+1);
    else
        powerLFPshiftedshifted=powerLFPshiftedfilt(1:end-ii);
        powerLFP1shifted=powerLFP1filt(ii+1:end);
    end
    if ii_shift==13
       pffft=1; 
    end
    rho(ii_shift)=corr(powerLFP1shifted',powerLFPshiftedshifted');
end

pffft=1;

%Uncomment to plot figure
% try
%     close 3
% catch
% end
% 
% hFig3 = figure(3);
% set(hFig3, 'units','normalized','position',[.3 .3 .3 .3])
% 
% plot(t_lag,rho,'-b')
% xlabel('Lag (sec)')
% ylabel('Rho')



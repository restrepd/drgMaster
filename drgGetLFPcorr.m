function [rho, t_lag]=drgGetLFPcorr(LFP1,LFPshifted,Fs,F1,F2,time_pad)
%Generates the phase histogram for the envelope and pac
%function [pac_value, mod_indx, phase, phase_histo, theta_wave]=drgGetThetaAmpPhase(LFP,Fs,lowF1,lowF2,highF1,highF2,time_pad,no_bins)
 
%See Adhikari et al J. Neurosci Meth. 191:191 and Adhikari et al Neuron
%65:257, 2010
 
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
%Note: I tried to use power (LFP squared), it did not make a difference
% powerLFP1filt=LFP1env.^2;
% powerLFP1filt=LFP1env;
LFP1env=LFP1env-mean(LFP1env);
angleLFP1 = angle(hilbert(filtLFP1)); % phase modulation of amplitude

filtLFPshifted=filtfilt(bpFilt,LFPshifted);
LFPshiftedenv = abs(hilbert(filtLFPshifted)); % envelope
% powerLFPshiftedfilt=LFPshiftedenv.^2;
% powerLFPshiftedfilt=LFPshiftedenv;
LFPshiftedenv=LFPshiftedenv-mean(LFPshiftedenv);
angleLFPshifted = angle(hilbert(filtLFPshifted)); % phase modulation of amplitude

max_t_shift=(1/mean([F1,F2]));
max_ii_shift=ceil(max_t_shift*Fs);


if length(LFP1env)~=length(LFPshiftedenv)
    min_length=min([length(LFP1env) length(LFPshiftedenv)]);
    LFP1env=LFP1env(1:min_length);
    LFPshiftedenv=LFPshiftedenv(1:min_length);
end

[rho, lags]=xcorr(LFP1env,LFPshiftedenv,max_ii_shift,'coeff');
t_lag=lags/Fs;


% %How about using the angle?
% %This did not work
% if length(angleLFP1)~=length(angleLFPshifted)
%     min_length=min([length(powerLFP1filt) length(angleLFPshifted)]);
%     angleLFP1=angleLFP1(1:min_length);
%     angleLFPshifted=angleLFPshifted(1:min_length);
% end
% 
% [rho_angle,lags_angle]=xcorr(angleLFP1,angleLFPshifted,max_ii_shift,'coeff');
% t_lag_angle=lags/Fs;

pffft=1;

% figure(1)
% sub_phase=angleLFPshifted-angleLFP1;
% sub_phase_shifted=
% plot(angleLFPshifted-angleLFP1,'-b')

% 
% n_bins=81;
% t_lag=[];
% rho=[];
% for ii_shift=1:n_bins
%     ii=ceil(-max_ii_shift+(ii_shift-1)*(max_ii_shift/((n_bins-1)/2)));
%     t_lag(ii_shift)=ii/Fs;
%     if ii<0
%         powerLFPshiftedshifted=powerLFPshiftedfilt(-ii:end);
%         powerLFP1shifted=powerLFP1filt(1:end+ii+1);
%     else
%         powerLFPshiftedshifted=powerLFPshiftedfilt(1:end-ii);
%         powerLFP1shifted=powerLFP1filt(ii+1:end);
%     end
% 
% 
%     if length(powerLFP1shifted)==length(powerLFPshiftedshifted)
%         rho(ii_shift)=corr(powerLFP1shifted',powerLFPshiftedshifted');
%     else
%         if length(powerLFP1shifted)>length(powerLFPshiftedshifted)
%             rho(ii_shift)=corr(powerLFP1shifted(1:length(powerLFPshiftedshifted))',powerLFPshiftedshifted');
%         else
%             rho(ii_shift)=corr(powerLFP1shifted',powerLFPshiftedshifted(1:length(powerLFP1shifted))');
%         end
%     end
%     
% 
% end
% 
% pffft=1;
% 
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
% 
% try
%     close 4
% catch
% end
% 
% hFig4 = figure(4);
% set(hFig4, 'units','normalized','position',[.3 .3 .3 .3])
% 
% plot(LFP1env,'-b')
% hold on
% plot(LFPshiftedenv,'-r')
% 
% 
% 
% pffft=1
% 
% try
%     close 5
% catch
% end
% 
% hFig5 = figure(5);
% set(hFig5, 'units','normalized','position',[.3 .3 .3 .3])
% 
% plot(t_lag,rho,'-b')
% xlabel('Lag (sec)')
% ylabel('Rho')
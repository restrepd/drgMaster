function [meanVectorLength, meanVectorAngle, peakAngle, peakAngles, MI_Tort, phase, phase_histo, theta_wave, meanPeakAngle, out_times, out_phase, out_time_PAChisto, decLFPgenv, decanglethetaLFP, out_times_env]=drgGetThetaAmpPhase(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,time_pad,no_bins,method)
%Generates the phase histogram for the envelope and pac
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

anglethetaLFP = angle(hilbert(thfiltLFP)); % theta phase

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

%Find the timecourse of the phase
min_dt=(1/lowF2); %Skip peaks if they are closer than this dt

%Find the peaks of the phase that define each cycle
[peaks,loc]=findpeaks(anglethetaLFP,Fs,'MinPeakHeight',2,'MinPeakDistance',min_dt);
time_histos=zeros(length(loc)-1,no_bins);
tstart=zeros(length(loc)-1,1);
tend=zeros(length(loc)-1,1);
for ii=1:length(loc)-1
    if (loc(ii)*Fs>ii_pad)&(loc(ii)*Fs<length(anglethetaLFP)-ii_pad)
        
        %Find the peak for the envelope
        this_LFPgenv=LFPgenv(ceil(loc(ii)*Fs):ceil(loc(ii+1)*Fs));
        this_anglethetaLFP=anglethetaLFP(ceil(loc(ii)*Fs):ceil(loc(ii+1)*Fs));
        [maxamp,max_ii]=max(this_LFPgenv);
        
        peak_amp_phase_timecourse_t(ii)=((max_ii+ceil(loc(ii)*Fs))/Fs);
        peak_amp_phase_timecourse_phase(ii)=this_anglethetaLFP(max_ii);
        
        for jj=1:no_bins
            this_angle=-pi+jj*(2*pi/no_bins)-0.5*(2*pi/no_bins);
            [minkk,kk]=min(abs(anglethetaLFP(ceil(loc(ii)*Fs):ceil(loc(ii+1)*Fs))-this_angle));
            time_histos(ii,jj)=LFPgenv(ceil(loc(ii)*Fs)+kk);
            tstart(jj)=loc(ii)-time_pad;
            tend(jj)=loc(ii+1)-time_pad-1/Fs;
        end
        
        
    end
end


%Normalize the time histo
for ii=1:length(loc)-1
    time_histos(ii,:)=time_histos(ii,:)/sum(time_histos(ii,:));
end

meanPeakAngle=(360/(2*pi))*(circ_mean(peak_amp_phase_timecourse_phase')+pi);

%Create an output evenly spaced time vector for peak of phase
out_times=[0.05:0.1:(length(anglethetaLFP)/Fs)-2*time_pad];
out_phase=-1000*ones(1,length(out_times));

for jj=1:round(((length(anglethetaLFP)/Fs)-2*time_pad)/0.1)
    ii_this_dt=find(((peak_amp_phase_timecourse_t-time_pad)>out_times(jj)-0.05)&((peak_amp_phase_timecourse_t-time_pad)<=out_times(jj)+0.05));
    if ~isempty(ii_this_dt)
        out_phase(jj)=circ_mean(peak_amp_phase_timecourse_phase(((peak_amp_phase_timecourse_t-time_pad)>out_times(jj)-0.05)&((peak_amp_phase_timecourse_t-time_pad)<=out_times(jj)+0.05))');
    end
end
   
jj=0;
at_end=0;
to_average=[];
while at_end==0
    ii_no_phase=[];
    jj=jj+1;

    if (out_phase(jj)==-1000)
        this_out_phase=-1000;
        while (this_out_phase==-1000)
            ii_no_phase=[ii_no_phase jj];
            jj=jj+1;
            if jj>length(out_phase)
                at_end=1;
                this_out_phase=0;
                for kk=1:length(ii_no_phase)
                    out_phase(ii_no_phase(kk))=circ_mean(to_average');
                end
            else
               this_out_phase=out_phase(jj);
            end
        end
        if at_end==0
            to_average=[to_average out_phase(jj)];
        end
        for kk=1:length(ii_no_phase)
            out_phase(ii_no_phase(kk))=circ_mean(to_average');
        end
        if at_end==0
            to_average=out_phase(jj);
        end
        if (jj==length(out_phase))
            at_end=1;
        end
    else
        to_average=out_phase(jj);
        if jj==length(out_phase)
            at_end=1;
        end
    end
end

out_phase=(360/(2*pi))*(out_phase+pi);

out_time_PAChisto=zeros(length(out_times),no_bins);

for iit=1:length(out_times)
    iiloc_bef=find(loc<out_times(iit)+time_pad,1,'last');
    if iiloc_bef==length(loc)
        iiloc_bef=length(loc)-1;
    end
    
    if ~isempty(iiloc_bef)
        out_time_PAChisto(iit,:)=time_histos(iiloc_bef,:);
    else
        out_time_PAChisto(iit,:)=time_histos(1,:);
    end
%     out_time_PAChisto(iit,:)=time_histos(iiloc_bef,:);
end

dec_n=fix(Fs/1000);
decLFPgenv=decimate(LFPgenv(int64(time_pad*Fs):end-int64(time_pad*Fs)-1),dec_n);
decanglethetaLFP=decimate(anglethetaLFP(int64(time_pad*Fs):end-int64(time_pad*Fs)-1),dec_n);
% decLFPgenv=decimate(LFPgenv(time_pad*Fs:end-time_pad*Fs-1),20);
% decanglethetaLFP=decimate(anglethetaLFP(time_pad*Fs:end-time_pad*Fs-1),20);
out_times_env=[1:length(decLFPgenv)]*(1/(Fs/dec_n));

peakAngles=peak_amp_phase_timecourse_phase;

%Note peak_amp_phase_timecourse_t goes from zero to
%(time_end-time_start)+2*time_pad


% % %Plot for troubleshooting
% figure(10)
% subplot(3,1,1)
% plot(LFPgenv)
% hold on
% for ii=1:length(loc)-1
%     plot([peak_amp_phase_timecourse_t(ii)*Fs peak_amp_phase_timecourse_t(ii)*Fs],[min(LFPgenv) max(LFPgenv)],'-r')
% end
% title('Envelope for high frequency LFP')
% 
% %Plot the theta phase
% subplot(3,1,2)
% plot(anglethetaLFP)
% title('Theta phase')
% 
% subplot(3,1,3)
% plot(anglethetaLFP)
% hold on
% for ii=1:length(loc)-1
%     plot(peak_amp_phase_timecourse_t(ii)*Fs,peak_amp_phase_timecourse_phase(ii),'or')
% end
% plot(loc*Fs,(max(anglethetaLFP)+0.05*(max(anglethetaLFP)-min(anglethetaLFP))*ones(1,length(loc))),'ob')
% plot((out_times+time_pad)*Fs,out_phase,'ok')
% title('Angle for low frequency LFP/sniff')

% %Plot for Justin's Fig. 4 (example of PAC wavelet power)
% %
% figure(16)
% ii_pad=ceil(time_pad*Fs);
% trimmed_anglethetaLFP=anglethetaLFP(ii_pad:end-ii_pad);
% timeAngle=([1:length(trimmed_anglethetaLFP)]/Fs)-2;
% plot(timeAngle,(180/pi)*(trimmed_anglethetaLFP+pi),'-b')
% xlim([-2 5])

% 
% % Plots for Daniel's Fig. 1
% % % Raw data
% % set(gcf,'renderer','zbuffer')
% % set(gcf,'renderer','painters');
% 
% try
%     close(17)
% catch
% end
% hFig=figure(17);
% set(hFig, 'units','normalized','position',[.02 .05 .5 .2])
% time=[1:length(LFPhigh(time_pad*Fs+1:end-time_pad*Fs))]/Fs;
% plot(time,LFPhigh(time_pad*Fs+1:end-time_pad*Fs),'-b','LineWidth',0.2)
% xlabel('Time (sec)')
% ylim([-5000 5000])
% ylabel('uv')
% 
%  
% try
%     close(18)
% catch
% end
% hFig=figure(18);
% set(hFig, 'units','normalized','position',[.02 .05 .5 .2])
% time=[1:length(LFPhigh(time_pad*Fs+1:end-time_pad*Fs))]/Fs;
% plot(time,LFPhigh(time_pad*Fs+1:end-time_pad*Fs),'-b','LineWidth',0.2)
% xlabel('Time (sec)')
% ylim([-5000 5000])
% xlim([3 4])
% ylabel('uv')
% 
% %Theta
% try
%     close(19)
% catch
% end
% hFig=figure(19);
% set(hFig, 'units','normalized','position',[.02 .25 .5 .2])
% time=[1:length(thfiltLFP(time_pad*Fs+1:end-time_pad*Fs))]/Fs;
% plot(time,thfiltLFP(time_pad*Fs+1:end-time_pad*Fs),'-b')
% xlabel('Time (sec)')
% ylabel('uv')
% ylim([-7000 7000])
% xlim([3 4])
% 
% 
% 
% %Theta phase
% try
%     close(20)
% catch
% end
% hFig=figure(20);
% set(hFig, 'units','normalized','position',[.02 .45 .5 .2])
% time=[1:length(anglethetaLFP(time_pad*Fs+1:end-time_pad*Fs))]/Fs;
% plot(time,anglethetaLFP(time_pad*Fs+1:end-time_pad*Fs)*(180/pi),'-b')
% xlabel('Time (sec)')
% ylabel('degrees')
% ylim([-200 200])
% xlim([3 4])
% 
% 
% %Gamma + envelope
% try
%     close(21)
% catch
% end
% hFig=figure(21);
% set(hFig, 'units','normalized','position',[.02 .65 .5 .2])
% hold on
% time=[1:length(LFPgenv(time_pad*Fs+1:end-time_pad*Fs))]/Fs;
% plot(time,filtLFPgamma(time_pad*Fs+1:end-time_pad*Fs),'-b')
% plot(time,LFPgenv(time_pad*Fs+1:end-time_pad*Fs),'-r')
% plot(time,-LFPgenv(time_pad*Fs+1:end-time_pad*Fs),'-r')
% xlabel('Time (sec)')
% ylabel('uv')
% xlim([3 4])
% ylim([-500 500])
% 
%Plots for Daniel's new Figure 5 (licks vs theta)
%These are used with LFP 1
% 
% 
% %Theta
% try
%     close(21)
% catch
% end
% hFig=figure(21);
% set(hFig, 'units','normalized','position',[.02 .37 .5 .2])
% hold on
% ax=gca;ax.LineWidth=3;
% time=[1:length(thfiltLFP(time_pad*Fs+1:end-time_pad*Fs))]/Fs;
% plot(time,thfiltLFP(time_pad*Fs+1:end-time_pad*Fs),'-b','LineWidth',2)
% plot(time,30*LFPhigh(time_pad*Fs+1:end-time_pad*Fs)*(180/pi)-3000,'-','LineWidth',1,'Color',[0.7 0.7 0.7])
% xlabel('Time (sec)')
% ylabel('uv')
% ylim([-4500 4500])
% 
% 
% %Licks
% try
%     close(22)
% catch
% end
% hFig=figure(22);
% set(hFig, 'units','normalized','position',[.02 .65 .5 .2])
% ax=gca;ax.LineWidth=3;
% time=[1:length(anglethetaLFP(time_pad*Fs+1:end-time_pad*Fs))]/Fs;
% plot(time,LFPhigh(time_pad*Fs+1:end-time_pad*Fs)*(180/pi),'-k','LineWidth',2)
% xlabel('Time (sec)')
% % ylabel('degrees')
% ylim([-50 250])
% 
% %Theta
% try
%     close(19)
% catch
% end
% hFig=figure(19);
% set(hFig, 'units','normalized','position',[.02 .37 .5 .2])
% hold on
% ax=gca;ax.LineWidth=3;
% time=[1:length(thfiltLFP(time_pad*Fs+1:end-time_pad*Fs))]/Fs;
% plot(time,thfiltLFP(time_pad*Fs+1:end-time_pad*Fs),'-b','LineWidth',2)
% plot(time,30*LFPhigh(time_pad*Fs+1:end-time_pad*Fs)*(180/pi)-3000,'-','LineWidth',1,'Color',[0.7 0.7 0.7])
% xlabel('Time (sec)')
% ylabel('uv')
% ylim([-4500 4500])
% xlim([2.5 4.5])
% 
% %Licks
% try
%     close(20)
% catch
% end
% hFig=figure(20);
% set(hFig, 'units','normalized','position',[.02 .65 .5 .2])
% ax=gca;ax.LineWidth=3;
% time=[1:length(anglethetaLFP(time_pad*Fs+1:end-time_pad*Fs))]/Fs;
% plot(time,LFPhigh(time_pad*Fs+1:end-time_pad*Fs)*(180/pi),'-k','LineWidth',2)
% xlabel('Time (sec)')
% % ylabel('degrees')
% ylim([-50 250])
% xlim([2.5 4.5])


% %Theta phase
% try
%     close(20)
% catch
% end
% hFig=figure(20);
% set(hFig, 'units','normalized','position',[.02 .09 .5 .2])
% ax=gca;ax.LineWidth=3;
% hold on
% time=[1:length(anglethetaLFP(time_pad*Fs+1:end-time_pad*Fs))]/Fs;
% plot(time,anglethetaLFP(time_pad*Fs+1:end-time_pad*Fs)*(180/pi),'-b','LineWidth',3)
% plot(time,2.2*LFPhigh(time_pad*Fs+1:end-time_pad*Fs)*(180/pi)-240,'-','LineWidth',1,'Color',[0.7 0.7 0.7])
% xlabel('Time (sec)')
% ylabel('degrees')
% ylim([-250 250])
% xlim([3 4])


% %Gamma + envelope
% try
%     close(21)
% catch
% end
% hFig=figure(21);
% set(hFig, 'units','normalized','position',[.02 .65 .5 .2])
% ax=gca;ax.LineWidth=3;
% hold on
% time=[1:length(LFPgenv(time_pad*Fs+1:end-time_pad*Fs))]/Fs;
% plot(time,filtLFPgamma(time_pad*Fs+1:end-time_pad*Fs),'-b','LineWidth',3)
% plot(time,LFPgenv(time_pad*Fs+1:end-time_pad*Fs),'-r','LineWidth',3)
% plot(time,-LFPgenv(time_pad*Fs+1:end-time_pad*Fs),'-r','LineWidth',3)
% xlabel('Time (sec)')
% ylabel('uv')
% xlim([3 4])
% ylim([-1200 1200])

% %Licks
% try
%     close(22)
% catch
% end
% hFig=figure(22);
% set(hFig, 'units','normalized','position',[.02 .65 .5 .2])
% ax=gca;ax.LineWidth=3;
% time=[1:length(anglethetaLFP(time_pad*Fs+1:end-time_pad*Fs))]/Fs;
% plot(time,LFPhigh(time_pad*Fs+1:end-time_pad*Fs)*(180/pi),'-k','LineWidth',3)
% xlabel('Time (sec)')
% % ylabel('degrees')
% ylim([-50 250])
% xlim([3 4])

       
% %Uncomment if you want to save the theta LFP to test with closed loop
% try
%     load('C:\Users\Diego Restrepo\Desktop\theta.mat')
%     theta_out=[theta_out thfiltLFP];
% catch
%     theta_out=thfiltLFP;
% end
% save('C:\Users\Diego Restrepo\Desktop\theta.mat','theta_out')
 
pffft=1;





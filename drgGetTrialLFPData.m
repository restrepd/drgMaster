function [dataforch,trialNo,can_read] = drgGetTrialLFPData(handles, lfpElectrode, evNo, evTypeNo, time_start, time_end)
%This function either extracts the LFP data from drg.data_dg
%or generates simulated data with specific theta/gamma oscillations
%evNo is the event number for the evTypeNo

%If testLFP==0 the function reads the trial
%If testLFP~=0 the function simulates data

can_read=1;


testLFP=handles.data_vs_simulate;
% if testLFP>0
%     if time_start==handles.startRef+handles.time_pad
%         if time_end==handles.endRef-handles.time_pad
%             testLFP=4;  %This is the reference period in a simulation trial
%         end
%     end
% end

sessionNo=handles.sessionNo;
LFPtheta=[];
LFPtheta2=[];
LFPbursts=[];

try
    trialNo=find(handles.drg.session(sessionNo).trial_start<handles.drg.session(sessionNo).events(evTypeNo).times(evNo),1,'last');
catch
    %pffft=1
end


data_allch = drgGetThisTrial(handles,evNo);

% data_allch=zeros(floor(handles.drg.session(sessionNo).draq_p.ActualRate*handles.drg.session(sessionNo).draq_p.sec_per_trigger),handles.drg.session(sessionNo).draq_p.no_chans);
% for ii=1:handles.drg.session(sessionNo).draq_p.no_chans
%     data_allch(1:end-2000,ii)=data_this_trial(floor((ii-1)*handles.drg.session(sessionNo).draq_p.ActualRate*handles.drg.session(sessionNo).draq_p.sec_per_trigger)+1:...
%         floor((ii-1)*handles.drg.session(sessionNo).draq_p.ActualRate*handles.drg.session(sessionNo).draq_p.sec_per_trigger)...
%         +floor(handles.drg.session(sessionNo).draq_p.ActualRate*handles.drg.session(sessionNo).draq_p.sec_per_trigger)-2000);
% end
data=[];
data=data_allch(:,lfpElectrode);


if handles.notch60==1
    notch60HzFilt = designfilt('bandstopiir','FilterOrder',2, ...
        'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
        'DesignMethod','butter','SampleRate',floor(handles.drg.session(sessionNo).draq_p.ActualRate));
    data=filtfilt(notch60HzFilt,data);
end

%Get the data
ii_from=floor(((handles.drg.session(sessionNo).events(evTypeNo).times(evNo)-handles.drg.session(sessionNo).trial_start(trialNo))+time_start)...
    *handles.drg.draq_p.ActualRate+1);
ii_to=floor(((handles.drg.session(sessionNo).events(evTypeNo).times(evNo)-handles.drg.session(sessionNo).trial_start(trialNo))...
    +time_end)*handles.drg.draq_p.ActualRate);
ii_zero1=floor(((handles.drg.session(sessionNo).events(evTypeNo).times(evNo)-handles.drg.session(sessionNo).trial_start(trialNo)))...
    *handles.drg.draq_p.ActualRate+1);

dataforch=[];
if ii_to<=length(data)
    if (ii_from<1)
%         can_read=0;
        dataforch=[zeros(-ii_from+1,1); data(1:ii_to)];
    else
        dataforch=data(ii_from:ii_to);
    end
else
    %data is too short
%     can_read=0;
    dataforch=[data(ii_from:end);zeros(ii_to-length(data)+1,1)];
end

if testLFP==0
    dataforch=dataforch';
else
    dataforch=[];
    %Setup a signal with high and low frequency oscillations
    Fs=handles.drg.draq_p.ActualRate;
    t = 0:1/Fs:time_end-time_start;
    t_tort=1:1:length(t);
    
    %Theta
    Ftheta=handles.phaseHz;
    LFPtheta = cos(2*pi*Ftheta*t);
    
    %This is a Gaussian high gamma burst
    FWHMdt=(1/Ftheta)/3; %0.0440 msec FWHM burst, 25 Hz
    FWHM=ceil(Fs*FWHMdt);
    X=[1:3*FWHM]; %Burst is calculated for 3*FWHM samples with peak in the middle
    k=normpdf(X,3*FWHM/2,FWHM/2.35482); %normpdf(x,mu,sigma), FWHM=2*sqrt(2*ln(2))=2.3548*sigma
    k=k/max(k); %Normalized to 1 at peak
    Fgamma=handles.amplitudeHz; %75
    g_burst=k.*cos(2*pi*Fgamma*t(1:length(k)));
    
    
    LFPbursts=zeros(1,length(LFPtheta));
    
    %Important: It is key to have noise, otherwise the fourier is broad
    noiseFact=0.2;   %0.2
    
    dataforch=zeros(1,length(LFPtheta));
    time=time_start:(time_end-time_start)/length(dataforch):time_end;
    time=time(1:end-1);
    
    odorLFP=~((time<0)|(time>2.5));
    no_odorLFP=((time<0)|(time>2.5));
    
    switch testLFP
        
        case 1
            %LFP gamma burst in trough, high MI
            gammaFact=0.5;
            shift_burst=-ceil(length(k)/2);
            shift_trough=0.5*(1/Ftheta)*Fs; %trough
            for ii=1:floor((time_end-time_start)*Ftheta)
                LFPbursts(floor(ii*(Fs/Ftheta)+shift_burst+shift_trough)+1:floor(ii*(Fs/Ftheta)+shift_trough+shift_burst)+length(k))=gammaFact*g_burst;
            end
            dataforch=(LFPtheta+LFPbursts(1:length(LFPtheta))).*(0.05*no_odorLFP+odorLFP)+noiseFact*randn(1,length(LFPtheta));
        case 2
            %LFP gamma burst in peak, high MI
            gammaFact=0.5;
            shift_burst=-ceil(length(k)/2);
            shift_trough=0; %peak
            for ii=1:floor((time_end-time_start)*Ftheta)
                LFPbursts(floor(ii*(Fs/Ftheta)+shift_burst+shift_trough)+1:floor(ii*(Fs/Ftheta)+shift_trough+shift_burst)+length(k))=gammaFact*g_burst;
            end
            dataforch=(LFPtheta+LFPbursts(1:length(LFPtheta))).*(0.05*no_odorLFP+odorLFP)+noiseFact*randn(1,length(LFPtheta));
        case 3
            %LFP gamma burst in trough, low MI
            gammaFact=0.1;
            shift_burst=-ceil(length(k)/2);
            shift_trough=0.5*(1/Ftheta)*Fs; %trough
            for ii=1:floor((time_end-time_start)*Ftheta)
                LFPbursts(floor(ii*(Fs/Ftheta)+shift_burst+shift_trough)+1:floor(ii*(Fs/Ftheta)+shift_trough+shift_burst)+length(k))=gammaFact*g_burst;
            end
            dataforch=(LFPtheta+LFPbursts(1:length(LFPtheta))).*(0.05*no_odorLFP+odorLFP)+noiseFact*randn(1,length(LFPtheta));
        case 4
            %LFP gamma burst in peak, low MI
            gammaFact=0.1;
            shift_burst=-ceil(length(k)/2);
            shift_trough=0; %peak
            for ii=1:floor((time_end-time_start)*Ftheta)
                LFPbursts(floor(ii*(Fs/Ftheta)+shift_burst+shift_trough)+1:floor(ii*(Fs/Ftheta)+shift_trough+shift_burst)+length(k))=gammaFact*g_burst;
            end
            dataforch=(LFPtheta+LFPbursts(1:length(LFPtheta))).*(0.05*no_odorLFP+odorLFP)+noiseFact*randn(1,length(LFPtheta));
        case 5
            %Steady gamma
%             gammaFact=2;
%             dataforch = gammaFact*cos(2*pi*Fgamma*t)+noiseFact*randn(1,length(LFPtheta));
            
            dataforch=sin(2*pi*t_tort*Fgamma/Fs);
            dataforch= dataforch+noiseFact*randn(1,length(dataforch));
            dataforch(1:-time_start*handles.drg.draq_p.ActualRate)=0;
            dataforch(((-time_start+2.5)*handles.drg.draq_p.ActualRate):end)=0;
        case 6
            %Noise
            dataforch=noiseFact*randn(1,length(LFPtheta));
        case 7
            %Tort's PAC simulation
%             nonmodulatedamplitude=2;
%             dataforch=(0.2*(sin(2*pi*t_tort*Ftheta/Fs)+1)+nonmodulatedamplitude*0.1).*sin(2*pi*t_tort*Fgamma/Fs)+sin(2*pi*t_tort*Ftheta/Fs);
%             dataforch=dataforch+1*randn(1,length(dataforch));
            
            %Tort spectrogram simulation
            dataforch=sin(2*pi*t_tort*Ftheta/Fs);
            dataforch(1:-time_start*handles.drg.draq_p.ActualRate)=0;
            dataforch(((-time_start+2.5)*handles.drg.draq_p.ActualRate):end)=0;
    end
    
    %     try
    %         close 4
    %     catch
    %     end
    %
    %     hFig4 = figure(4);
    %     set(hFig4, 'units','normalized','position',[.55 .15 .35 .25])
    %     plot(time,dataforch);
    
    
    
    %pffft=1
end




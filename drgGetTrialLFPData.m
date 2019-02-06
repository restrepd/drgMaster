function [dataforch,trialNo,can_read] = drgGetTrialLFPData(handles, lfpElectrode, evNo, evTypeNo, time_start, time_end)
%This function either extracts the LFP data from drg.data_dg
%or generates simulated data with specific theta/gamma oscillations
%evNo is the event number for the evTypeNo

%If testLFP==0 the function reads the trial
%If testLFP~=0 the function simulates data assuming that the trace is noise
%before t=0, an oscillatory signal from 0 to 2.5 and noise from 2.5 onwards

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
    

    data_allch = drgGetThisTrial(handles,evNo,evTypeNo);
    
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
    
    try
        if handles.notch50==1
            notch50HzFilt = designfilt('bandstopiir','FilterOrder',2, ...
                'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
                'DesignMethod','butter','SampleRate',floor(handles.drg.session(sessionNo).draq_p.ActualRate));
            data=filtfilt(notch50HzFilt,data);
        end
    catch
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
        
        %Generate an oscillatory signal for the entire length (this is
        %trimmed below)
        t = 0:1/Fs:time_end-time_start;
        t_tort=1:1:length(t);
        
        %Theta
        Ftheta=handles.phaseHz;
        LFPtheta = cos(2*pi*Ftheta*t);
        
        %This is a Gaussian high gamma burst
        FWHMdt=(1/Ftheta)/2; %0.0440 msec FWHM burst, 25 Hz
        FWHM=ceil(Fs*FWHMdt);
        X=[1:3*FWHM]; %Burst is calculated for 3*FWHM samples with peak in the middle
        k=normpdf(X,3*FWHM/2,FWHM/2.35482); %normpdf(x,mu,sigma), FWHM=2*sqrt(2*ln(2))=2.3548*sigma
        k=k/max(k); %Normalized to 1 at peak
        Fgamma=handles.amplitudeHz; %75
        g_burst=k.*cos(2*pi*Fgamma*t(1:length(k)));
        
        
        LFPbursts=zeros(1,length(LFPtheta));
        
        %Important: It is key to have noise, otherwise the fourier is broad
        noiseFact=0.1;   %0.2
        
        dataforch=zeros(1,length(LFPtheta));
        time=time_start:(time_end-time_start)/length(dataforch):time_end;
        time=time(1:end-1);
        
        odorLFP=~((time<0)|(time>2.5));
        no_odorLFP=((time<0)|(time>2.5));
        
        %Genetrate pink noise
        % https://www.dsprelated.com/freebooks/sasp/Example_Synthesis_1_F_Noise.html
        % Pink noise``1/f noise''
        Nx = length(LFPtheta);  % number of samples to synthesize
        B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
        A = [1 -2.494956002   2.017265875  -0.522189400];
        nT60 = round(log(1000)/(1-max(abs(roots(A))))); % T60 est.
        v = randn(1,Nx+nT60); % Gaussian white noise: N(0,1)
        pink_noise = filter(B,A,v);    % Apply 1/F roll-off to PSD
        pink_noise = pink_noise(nT60+1:end);    % Skip transient response
        pink_n_fact=2000;
        
        %Generate licks at the troughs of the theta LFP
        [peaks,loc]=findpeaks(-LFPtheta,1,'MinPeakHeight',0.5,'MinPeakDistance',0);
        licks=zeros(1,length(LFPtheta));
        for peakNo=1:length(loc)
            licks(loc(peakNo):loc(peakNo)+200)=500;
        end
         
        switch testLFP
            
            case 1
                %LFP gamma burst at trough of theta, zero degrees, high MI
                gammaFact=15;
                shift_burst=-ceil(length(k)/2);
                shift_trough=0.5*(1/Ftheta)*Fs; %trough
                for ii=1:floor((time_end-time_start)*Ftheta)
                    LFPbursts(floor(ii*(Fs/Ftheta)+shift_burst+shift_trough)+1:floor(ii*(Fs/Ftheta)+shift_trough+shift_burst)+length(k))=gammaFact*g_burst;
                end
                if lfpElectrode==19
                    dataforch=licks;
                else
                    dataforch=500*(LFPtheta+0.03*LFPbursts(1:length(LFPtheta))).*(0.02*no_odorLFP+odorLFP)+pink_n_fact*pink_noise;
%                     dataforch=500*(LFPtheta+0.015*LFPbursts(1:length(LFPtheta))).*(0.02*no_odorLFP+odorLFP)+pink_n_fact*pink_noise;
                end
                pffft=1;
                %                 figure(10)
                %                 plot(dataforch)
            case 2
                %LFP gamma burst at 180 degrees, peak of theta high MI
                gammaFact=15;
                shift_burst=-ceil(length(k)/2);
                shift_trough=0; %peak
                for ii=1:floor((time_end-time_start)*Ftheta)
                    LFPbursts(floor(ii*(Fs/Ftheta)+shift_burst+shift_trough)+1:floor(ii*(Fs/Ftheta)+shift_trough+shift_burst)+length(k))=gammaFact*g_burst;
                end
                if lfpElectrode==19
                    dataforch=licks;
                else
                    dataforch=500*(LFPtheta+0.03*LFPbursts(1:length(LFPtheta))).*(0.02*no_odorLFP+odorLFP)+pink_n_fact*pink_noise;
                end
            case 3
                %LFP gamma burst at zero degrees, low MI
                gammaFact=10;
                shift_burst=-ceil(length(k)/2);
                shift_trough=0.5*(1/Ftheta)*Fs; %trough
                for ii=1:floor((time_end-time_start)*Ftheta)
                    LFPbursts(floor(ii*(Fs/Ftheta)+shift_burst+shift_trough)+1:floor(ii*(Fs/Ftheta)+shift_trough+shift_burst)+length(k))=gammaFact*g_burst;
                end
                if lfpElectrode==19
                    dataforch=licks;
                else
%                     dataforch=500*(LFPtheta+0.015*LFPbursts(1:length(LFPtheta))).*(0.02*no_odorLFP+odorLFP)+pink_n_fact*pink_noise;
                     dataforch=500*(LFPtheta+0.015*LFPbursts(1:length(LFPtheta))).*(0.02*no_odorLFP+odorLFP)+pink_n_fact*pink_noise;
                end
                
            case 4
                %LFP gamma burst at 180 degrees, low MI
                gammaFact=10;
                shift_burst=-ceil(length(k)/2);
                shift_trough=0; %peak
                for ii=1:floor((time_end-time_start)*Ftheta)
                    LFPbursts(floor(ii*(Fs/Ftheta)+shift_burst+shift_trough)+1:floor(ii*(Fs/Ftheta)+shift_trough+shift_burst)+length(k))=gammaFact*g_burst;
                end
                if lfpElectrode==19
                    dataforch=licks;
                else
                    dataforch=500*(LFPtheta+0.0125*LFPbursts(1:length(LFPtheta))).*(0.02*no_odorLFP+odorLFP)+pink_n_fact*pink_noise;
                end
            case 5
                
                %10 Hz and 80Hz 10 dB
                dataforch=0.002*(6*sin(2*pi*t_tort*Ftheta/Fs)+6*sin(2*pi*t_tort*Fgamma/Fs));
                dataforch= dataforch+noiseFact*randn(1,length(dataforch));
                dataforch=dataforch.*(0.02*no_odorLFP+odorLFP)+noiseFact*randn(1,length(dataforch));
            case 6
                %10 Hz and 80Hz 20 dB
                dataforch=0.005*(6*sin(2*pi*t_tort*Ftheta/Fs)+6*sin(2*pi*t_tort*Fgamma/Fs));
                dataforch= dataforch+noiseFact*randn(1,length(dataforch));
                dataforch=dataforch.*(0.02*no_odorLFP+odorLFP)+noiseFact*randn(1,length(dataforch));
            case 7
                %10 and 80Hz centered on the lick
                if lfpElectrode==19
                    %Licks starting at -1 and 1.5
                    licks=(~((time<1.5)|(time>2)))|(~((time<-1)|(time>-0.5)));
                    dataforch=1000*licks+2000;
                else
                    %0.4sec 80 Hz sine wave burst centered at the lick
                    LFPburst=(~((time<1.3)|(time>1.7)));
                    noLFPburst=~LFPburst;
                    dataforch=0.005*(6*sin(2*pi*t_tort*Ftheta/Fs)+6*sin(2*pi*t_tort*Fgamma/Fs));
                    dataforch=dataforch.*(0.05*noLFPburst+LFPburst)+noiseFact*randn(1,length(dataforch));
                end
            case 8
                %LFP gamma burst at zero degrees, low MI, low power
                gammaFact=1;
                shift_burst=-ceil(length(k)/2);
                shift_trough=0.5*(1/Ftheta)*Fs; %trough
                for ii=1:floor((time_end-time_start)*Ftheta)
                    LFPbursts(floor(ii*(Fs/Ftheta)+shift_burst+shift_trough)+1:floor(ii*(Fs/Ftheta)+shift_trough+shift_burst)+length(k))=gammaFact*g_burst;
                end
                if lfpElectrode==19
                    dataforch=licks;
                else
%                     dataforch=500*(LFPtheta+0.015*LFPbursts(1:length(LFPtheta))).*(0.02*no_odorLFP+odorLFP)+pink_n_fact*pink_noise;
                     dataforch=100*(LFPtheta+0.015*LFPbursts(1:length(LFPtheta))).*(0.02*no_odorLFP+odorLFP)+pink_n_fact*pink_noise;
                end
                
                % The cases below are used by batch processing
            case 9
                
                %if it is S-handles.drgbchoices.evTypeNos=[3 5 7 9 11 13 16:21];
   
                %First find out if this is S+ or S-
                
                if sum(handles.drg.session(1).events(5).times==handles.drg.session(1).events(evTypeNo).times(evNo))>0
                    %This is S+
                    %LFP gamma burst at 180 degrees, high MI
                    gammaFact=15;
                    shift_burst=-ceil(length(k)/2);
                    shift_trough=0; %peak
                    for ii=1:floor((time_end-time_start)*Ftheta)
                        LFPbursts(floor(ii*(Fs/Ftheta)+shift_burst+shift_trough)+1:floor(ii*(Fs/Ftheta)+shift_trough+shift_burst)+length(k))=gammaFact*g_burst;
                    end
                    if lfpElectrode==19
                        dataforch=licks;
                    else
                        dataforch=500*(LFPtheta+0.0125*LFPbursts(1:length(LFPtheta))).*(0.02*no_odorLFP+odorLFP)+pink_n_fact*pink_noise;
                    end
                    
                else
                    %This is S-
                    %LFP gamma burst at 180 degrees, low MI
                    gammaFact=10;
                    shift_burst=-ceil(length(k)/2);
                    shift_trough=0; %peak
                    for ii=1:floor((time_end-time_start)*Ftheta)
                        LFPbursts(floor(ii*(Fs/Ftheta)+shift_burst+shift_trough)+1:floor(ii*(Fs/Ftheta)+shift_trough+shift_burst)+length(k))=gammaFact*g_burst;
                    end
                    
                    if lfpElectrode==19
                        dataforch=licks;
                    else
                        dataforch=500*(LFPtheta+0.0125*LFPbursts(1:length(LFPtheta))).*(0.02*no_odorLFP+odorLFP)+pink_n_fact*pink_noise;
                    end
                end
                
            case 10
                
                %if it is S-handles.drgbchoices.evTypeNos=[3 5 7 9 11 13 16:21];
  
                %First find out if this is S+ or S-
                
                if sum(handles.drg.session(1).events(5).times==handles.drg.session(1).events(evTypeNo).times(evNo))>0
                    %This is S+
                   %LFP gamma burst at 180 degrees, high MI
                    gammaFact=15;
                    shift_burst=-ceil(length(k)/2);
                    shift_trough=0; %peak
                    for ii=1:floor((time_end-time_start)*Ftheta)
                        LFPbursts(floor(ii*(Fs/Ftheta)+shift_burst+shift_trough)+1:floor(ii*(Fs/Ftheta)+shift_trough+shift_burst)+length(k))=gammaFact*g_burst;
                    end
                    
                    if lfpElectrode==19
                        dataforch=licks;
                    else
                        dataforch=500*(LFPtheta+0.0125*LFPbursts(1:length(LFPtheta))).*(0.02*no_odorLFP+odorLFP)+pink_n_fact*pink_noise;
                    end
                else
                    %This is S-
                    %LFP gamma burst at zero degrees, high MI
                    gammaFact=15;
                    shift_burst=-ceil(length(k)/2);
                    shift_trough=0.5*(1/Ftheta)*Fs; %trough
                    for ii=1:floor((time_end-time_start)*Ftheta)
                        LFPbursts(floor(ii*(Fs/Ftheta)+shift_burst+shift_trough)+1:floor(ii*(Fs/Ftheta)+shift_trough+shift_burst)+length(k))=gammaFact*g_burst;
                    end
                    
                    
                    if lfpElectrode==19
                        dataforch=licks;
                    else
                        dataforch=500*(LFPtheta+0.015*LFPbursts(1:length(LFPtheta))).*(0.02*no_odorLFP+odorLFP)+pink_n_fact*pink_noise;
                    end
                end
        end
        

    end
    
catch
    can_read=0;
    dataforch=[];
    trialNo=-1;
end




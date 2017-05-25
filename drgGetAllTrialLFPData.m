function [dataforch,trialNo, ii_zero] = drgGetAllTrialLFPData(handles, lfpElectrode, evNo)
%This function either extracts the LFP data from drg.data_dg
%or generates simulated data with specific theta/gamma oscillations

%If testLFP==0 the function reads the trial
%If testLFP~=0 the function simulates data

testLFP=handles.data_vs_simulate;
evTypeNo=handles.evTypeNo;
sessionNo=handles.sessionNo;

trialNo=find(handles.drg.session(sessionNo).trial_start<handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),1,'last');

ii_zero=floor(((handles.drg.session(sessionNo).events(evTypeNo).times(evNo)-handles.drg.session(sessionNo).trial_start(trialNo)))...
    *handles.drg.draq_p.ActualRate+1);
    
if testLFP==0
    
    %Get data for this trial
    
    data_allch = drgGetThisTrial(handles,evNo);
    
%     data_allch=zeros(floor(handles.drg.session(sessionNo).draq_p.ActualRate*handles.drg.session(sessionNo).draq_p.sec_per_trigger),handles.drg.session(sessionNo).draq_p.no_chans);
%     for ii=1:handles.drg.session(sessionNo).draq_p.no_chans
%         data_allch(1:end-2000,ii)=data_this_trial(floor((ii-1)*handles.drg.session(sessionNo).draq_p.ActualRate*handles.drg.session(sessionNo).draq_p.sec_per_trigger)+1:...
%             floor((ii-1)*handles.drg.session(sessionNo).draq_p.ActualRate*handles.drg.session(sessionNo).draq_p.sec_per_trigger)...
%             +floor(handles.drg.session(sessionNo).draq_p.ActualRate*handles.drg.session(sessionNo).draq_p.sec_per_trigger)-2000);
%     end
    dataforch=[];
    dataforch=data_allch(:,lfpElectrode);
    %             for ii=1:handles.draq_p.no_chans
    %                 fseek(handles.p.fid, (ii-1)*size_per_ch_bytes+trial_offset, 'bof');
    %                 data(:,ii)=fread(handles.p.fid,no_unit16_per_ch,'uint16');
    %             end
    
    
    %Get the data
    
else
    %Setup a signal with high and low frequency oscillations
    Fs=handles.drg.draq_p.ActualRate;
    t = 0:1/Fs:handles.drg.session(sessionNo).draq_p.sec_per_trigger;
    
    %Theta
    Ftheta=8;
    LFPtheta = cos(2*pi*Ftheta*t);
    
    %This is a Gaussian gamma burst
    dt=0.04;
    FWHM=ceil(Fs*dt);
    X=[1:3*FWHM];
    k=normpdf(X,3*FWHM/2,FWHM/2.35482);
    k=k/max(k);
    Fgamma=75;
    g_burst=k.*cos(2*pi*Fgamma*t(1:length(k)));
    LFPgamma=zeros(1,length(LFPtheta));
    gammaFact=0.2;
    
    if (testLFP==1)||(testLFP==2)
        
        if (testLFP==1)
            shift_burst=-ceil(length(k)/2); %Trough
        else
            shift_burst=0; %Peak
        end
        
        %LFP gamma bursts
        for ii=0:ceil(2*Ftheta)-1
            LFPgamma(ceil((Fs/(2*Ftheta))+ii*(Fs/Ftheta)+shift_burst+1):ceil((Fs/(2*Ftheta))+ii*(Fs/Ftheta)+shift_burst+length(k)))=gammaFact*g_burst;
        end
        
    else
        %Continuous LFP
        LFPgamma = gammaFact*cos(2*pi*Fgamma*t);
    end
    
    LFPgamma=LFPgamma(1:length(LFPtheta));
    
    noiseFact=0.2;
    
    %LFP
    dataforch=LFPtheta+LFPgamma+noiseFact*randn(1,length(LFPtheta));
    
%     try
%         close 1
%     catch
%     end
%     
%     figure(1)
%     subplot(5,1,1)
%     plot(t,dataforch)
%     xlabel('Seconds')
%     title('LFP')
%     
end




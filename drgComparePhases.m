<<<<<<< HEAD
function handles=drgComparePhases(handles)

anglerefLFP = [];
anglesniff = [];
delta_phase_timecourse=[];

%Generates a trial per trial phase histogram
sessionNo=handles.sessionNo;
Fs=floor(handles.drg.session(sessionNo).draq_p.ActualRate);
lowF1=handles.peakLowF;
lowF2=handles.peakHighF;
highF1=handles.burstLowF;
highF2=handles.burstHighF;
pad_time=handles.time_pad;
n_phase_bins=handles.n_phase_bins;

%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;


%Now get the phase of gamma bursts witin theta
no_encoding_trials=0;
no_retrieval_trials=0;
no_trials=0;
enc_phase_histo=[];
retr_phase_histo=[];
all_phase_histo=[];
all_theta_wave=[];
MI_enc=[];
MI_retr=[];
which_event=[];

for trNo=firstTr:lastTr
    
    if handles.displayData==1
        trial_no=trNo
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            %Note: handles.peakLFPNo is the reference LFP
            [LFPref, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            [sniff, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);

            if (can_read1==1)&(can_read2==1)
                no_trials=no_trials+1;
                time(no_trials)=handles.drg.session(sessionNo).trial_start(trialNo);
                which_trial(no_trials)=1;
                perCorr_per_histo(no_trials)=50;
                
                %Get theta phase
                bpFilttheta = designfilt('bandpassiir','FilterOrder',20, ...
                    'HalfPowerFrequency1',lowF1,'HalfPowerFrequency2',lowF2, ...
                    'SampleRate',Fs);
                thfiltLFP=filtfilt(bpFilttheta,LFPref);
                thisanglerefLFP = angle(hilbert(thfiltLFP)); % phase modulation of theta amplitude
                anglerefLFP = [anglerefLFP thisanglerefLFP];
                
                %Get sniff phase
                bpFiltsniff = designfilt('bandpassiir','FilterOrder',20, ...
                    'HalfPowerFrequency1',highF1,'HalfPowerFrequency2',highF2, ...
                    'SampleRate',Fs);
                thfiltsniff=filtfilt(bpFiltsniff,sniff);
                thisanglesniff = angle(hilbert(thfiltsniff)); % sniff phase
                anglesniff = [anglesniff thisanglesniff];
                
                %Delta phase
                delta_phase_timecourse(no_trials,:)=thisanglesniff-thisanglerefLFP;
               
                
                %Plot each trial for troubleshooting
%                 try
%                     close 1
%                 catch
%                 end
%                 
%                 hFig1 = figure(1);
%                 set(hFig1, 'units','normalized','position',[.25 .35 .5 .2])
%                 
%                 plot(thfiltLFP,'-r')
%                 
%                 try
%                     close 2
%                 catch
%                 end
%                 
%                 hFig2 = figure(2);
%                 set(hFig2, 'units','normalized','position',[.25 .1 .5 .2])
%                 
%                 plot(sniff,'-b')
%                  hold on
%                  [yupper,yunder]=envelope(sniff,10000,'peak');
%                 plot(yupper,'-r')
%                 plot(yunder,'-r')
%                 
%                  try
%                     close 3
%                 catch
%                 end
%                 
%                 hFig3 = figure(3);
%                 set(hFig3, 'units','normalized','position',[.25 .6 .5 .2])
%                 
%                 plot(thisanglesniff,'-b')
%                 
%                 hold on
%                 
%                 maxsniff=max(sniff);
%                 minsniff=min(sniff);
%                 %plot(thisanglerefLFP,'-r')
%                 plot(2*pi()*(sniff-minsniff)/(maxsniff-minsniff) -pi(),'-r')
                
%                 pffft=1
            end
        end
    end
    
    %end %if eventstamps...
end %for evNo

if handles.displayData==1
    try
        close 2
    catch
    end
    
    hFig2 = figure(2);
    set(hFig2, 'units','normalized','position',[.05 .5 .3 .35])
    edges=[-2*pi():4*pi()/50:2*pi()];
    histogram(anglesniff-anglerefLFP,edges,'Normalization','probability')
    title('Histogram of the difference in phase between sniff and theta')
    ylabel('Probability')
    xlabel('Delta phase sniff - LFP')
    
    try
        close 3
    catch
    end
    
    hFig3 = figure(3);
    set(hFig3, 'units','normalized','position',[.05 .1 .7 .3])
    
    delta_t=0.2;
    time=handles.time_start+pad_time;
    delta_phase_t_hist=[];
    ii=0;
    
    while time<handles.time_end-pad_time
        ii=ii+1;
        t(ii)=time;
        deltat=time-handles.time_start+pad_time;
        these_phases=delta_phase_timecourse(:,floor(deltat*Fs):floor((deltat+delta_t)*Fs));
        delta_phase_t_hist(ii,1:50)=histcounts(these_phases(:),edges)/length(length(these_phases(:)));
        time=time+delta_t;
    end
    
    phase=[-2*pi()+(4*pi()/100):4*pi()/50:2*pi()-(4*pi()/100)]';
    drg_pcolor(repmat(t,50,1),repmat(phase,1,length(t)),delta_phase_t_hist')
    colormap jet
    shading flat
    
    min_prob=prctile(delta_phase_t_hist(:),1);
    max_prob=prctile(delta_phase_t_hist(:),99);
    caxis([min_prob    max_prob])
    
    xlabel('Time (sec)')
    ylabel('Radians');
    title('Timecourse for histogram of delta phase sniff - LFP (deg)')
    
    try
        close 4
    catch
    end
    
    hFig4 = figure(4);
    set(hFig4, 'units','normalized','position',[.76 .1 .05 .3])
    
    prain=[min_prob:(max_prob-min_prob)/99:max_prob];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    colormap jet
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')
    
    try
        close 1
    catch
    end
    
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.06 .45 .3 .35])
    
    plot(anglerefLFP,anglesniff,'ob','MarkerSize',0.6)
    xlabel(['Voltage electrode #' num2str(handles.peakLFPNo)])
    ylabel(['Voltage electrode #' num2str(handles.burstLFPNo)])
    xlim([-pi() pi()])
    ylim([-pi() pi()])
%     [rho,p_val]=corr(anglerefLFP,anglesniff);
%    title(['rho= ' num2str(rho) ' p_value= ' num2str(p_val)])
    
    pffft=1;
end





=======
function handles=drgComparePhases(handles)

anglerefLFP = [];
anglesniff = [];
delta_phase_timecourse=[];

%Generates a trial per trial phase histogram
sessionNo=handles.sessionNo;
Fs=floor(handles.drg.session(sessionNo).draq_p.ActualRate);
lowF1=handles.peakLowF;
lowF2=handles.peakHighF;
highF1=handles.burstLowF;
highF2=handles.burstHighF;
pad_time=handles.time_pad;
n_phase_bins=handles.n_phase_bins;

%Enter trials
firstTr=handles.trialNo;
lastTr=handles.lastTrialNo;


%Now get the phase of gamma bursts witin theta
no_encoding_trials=0;
no_retrieval_trials=0;
no_trials=0;
enc_phase_histo=[];
retr_phase_histo=[];
all_phase_histo=[];
all_theta_wave=[];
MI_enc=[];
MI_retr=[];
which_event=[];

for trNo=firstTr:lastTr
    
    if handles.displayData==1
        trial_no=trNo
    end
    
    evNo = drgFindEvNo(handles,trNo,sessionNo);
    if evNo~=-1
        excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
        
        if excludeTrial==0
            
            %Note: handles.peakLFPNo is the reference LFP
            [LFPref, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
            [sniff, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);

            if (can_read1==1)&(can_read2==1)
                no_trials=no_trials+1;
                time(no_trials)=handles.drg.session(sessionNo).trial_start(trialNo);
                which_trial(no_trials)=1;
                perCorr_per_histo(no_trials)=50;
                
                %Get theta phase
                bpFilttheta = designfilt('bandpassiir','FilterOrder',20, ...
                    'HalfPowerFrequency1',lowF1,'HalfPowerFrequency2',lowF2, ...
                    'SampleRate',Fs);
                thfiltLFP=filtfilt(bpFilttheta,LFPref);
                thisanglerefLFP = angle(hilbert(thfiltLFP)); % phase modulation of theta amplitude
                anglerefLFP = [anglerefLFP thisanglerefLFP];
                
                %Get sniff phase
                bpFiltsniff = designfilt('bandpassiir','FilterOrder',20, ...
                    'HalfPowerFrequency1',highF1,'HalfPowerFrequency2',highF2, ...
                    'SampleRate',Fs);
                thfiltsniff=filtfilt(bpFiltsniff,sniff);
                thisanglesniff = angle(hilbert(thfiltsniff)); % sniff phase
                anglesniff = [anglesniff thisanglesniff];
                
                %Delta phase
                delta_phase_timecourse(no_trials,:)=thisanglesniff-thisanglerefLFP;
               
                
                %Plot each trial for troubleshooting
%                 try
%                     close 1
%                 catch
%                 end
%                 
%                 hFig1 = figure(1);
%                 set(hFig1, 'units','normalized','position',[.25 .35 .5 .2])
%                 
%                 plot(thfiltLFP,'-r')
%                 
%                 try
%                     close 2
%                 catch
%                 end
%                 
%                 hFig2 = figure(2);
%                 set(hFig2, 'units','normalized','position',[.25 .1 .5 .2])
%                 
%                 plot(sniff,'-b')
%                  hold on
%                  [yupper,yunder]=envelope(sniff,10000,'peak');
%                 plot(yupper,'-r')
%                 plot(yunder,'-r')
%                 
%                  try
%                     close 3
%                 catch
%                 end
%                 
%                 hFig3 = figure(3);
%                 set(hFig3, 'units','normalized','position',[.25 .6 .5 .2])
%                 
%                 plot(thisanglesniff,'-b')
%                 
%                 hold on
%                 
%                 maxsniff=max(sniff);
%                 minsniff=min(sniff);
%                 %plot(thisanglerefLFP,'-r')
%                 plot(2*pi()*(sniff-minsniff)/(maxsniff-minsniff) -pi(),'-r')
                
%                 pffft=1
            end
        end
    end
    
    %end %if eventstamps...
end %for evNo

if handles.displayData==1
    try
        close 2
    catch
    end
    
    hFig2 = figure(2);
    set(hFig2, 'units','normalized','position',[.05 .5 .3 .35])
    edges=[-2*pi():4*pi()/50:2*pi()];
    histogram(anglesniff-anglerefLFP,edges,'Normalization','probability')
    title('Histogram of the difference in phase between sniff and theta')
    ylabel('Probability')
    xlabel('Delta phase sniff - LFP')
    
    try
        close 3
    catch
    end
    
    hFig3 = figure(3);
    set(hFig3, 'units','normalized','position',[.05 .1 .7 .3])
    
    delta_t=0.2;
    time=handles.time_start+pad_time;
    delta_phase_t_hist=[];
    ii=0;
    
    while time<handles.time_end-pad_time
        ii=ii+1;
        t(ii)=time;
        deltat=time-handles.time_start+pad_time;
        these_phases=delta_phase_timecourse(:,floor(deltat*Fs):floor((deltat+delta_t)*Fs));
        delta_phase_t_hist(ii,1:50)=histcounts(these_phases(:),edges)/length(length(these_phases(:)));
        time=time+delta_t;
    end
    
    phase=[-2*pi()+(4*pi()/100):4*pi()/50:2*pi()-(4*pi()/100)]';
    drg_pcolor(repmat(t,50,1),repmat(phase,1,length(t)),delta_phase_t_hist')
    colormap jet
    shading flat
    
    min_prob=prctile(delta_phase_t_hist(:),1);
    max_prob=prctile(delta_phase_t_hist(:),99);
    caxis([min_prob    max_prob])
    
    xlabel('Time (sec)')
    ylabel('Radians');
    title('Timecourse for histogram of delta phase sniff - LFP (deg)')
    
    try
        close 4
    catch
    end
    
    hFig4 = figure(4);
    set(hFig4, 'units','normalized','position',[.76 .1 .05 .3])
    
    prain=[min_prob:(max_prob-min_prob)/99:max_prob];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    colormap jet
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')
    
    try
        close 1
    catch
    end
    
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.06 .45 .3 .35])
    
    plot(anglerefLFP,anglesniff,'ob','MarkerSize',0.6)
    xlabel(['Voltage electrode #' num2str(handles.peakLFPNo)])
    ylabel(['Voltage electrode #' num2str(handles.burstLFPNo)])
    xlim([-pi() pi()])
    ylim([-pi() pi()])
%     [rho,p_val]=corr(anglerefLFP,anglesniff);
%    title(['rho= ' num2str(rho) ' p_value= ' num2str(p_val)])
    
    pffft=1;
end





>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855

function drgPhaseTimecourse(handles)

%Generates a trial per trial phase histogram
sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
lowF1=handles.peakLowF;
lowF2=handles.peakHighF;
highF1=handles.burstLowF;
highF2=handles.burstHighF;
pad_time=handles.time_pad;
n_phase_bins=handles.n_phase_bins;
phase_window=1;
phase_step=0.1; %0.25 sec

[perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);

try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.05 .5 .35 .4])


%Plot the percent correct
subplot(3,1,1)
trials=1:length(perCorr);
plot(trials,perCorr,'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
hold on
plot(trials(encoding_trials),perCorr(encoding_trials),'ob')
plot(trials(retrieval_trials),perCorr(retrieval_trials),'or')
ylim([40 110]);
xlabel('Trial No')
ylabel('Percent correct')
title('% correct vs trial number, blue=encoding, red=retrieval')

%Now get the phase of gamma bursts witin theta
no_encoding_trials=0;
no_retrieval_trials=0;
enc_phase_histo=[];
retr_phase_histo=[];
all_phase_histo=[];
all_theta_wave=[];
MI_enc=[]
MI_retr=[];
time_stamps_enc=[];
time_stamps_retr=[];


for evNo=1:handles.drg.session(handles.sessionNo).events(handles.evTypeNo).noTimes
    
    event_no=evNo
    
    excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(handles.evTypeNo).times(evNo),sessionNo);
    
    if excludeTrial==0
        
        %Can we read all the data?
        times=handles.time_start-phase_window/2:phase_step:handles.time_end-phase_window/2;
        [LFPlow, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, times(1), times(end)+phase_window);
        [LFPhigh, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, times(1), times(end)+phase_window);
        
        if (can_read1==1)&(can_read2==1)
            ii_t=0;
            for time=times
                
                ii_t=ii_t+1;
                
                [LFPlow, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, time, time+phase_window);
                [LFPhigh, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, time, time+phase_window);
                
                
                [meanVectorLength, meanVectorAngle, peakAngle, mod_indx, phase, phase_histo, theta_wave]=drgGetThetaAmpPhase(LFPlow,LFPhigh,Fs,lowF1,lowF2,highF1,highF2,pad_time,n_phase_bins,handles.which_method);
                
                if isnan(mod_indx)
                    pfffft=1
                end
                
                if encoding_trials(trialNo)==1
                    if ii_t==1
                        no_encoding_trials=no_encoding_trials+1;
                    end
                    enc_phase_histo(no_encoding_trials,ii_t,1:n_phase_bins+1)=phase_histo;
                    time_stamps_enc(ii_t,1:n_phase_bins+1)=time-phase_window/2;
                    MI_enc(no_encoding_trials,ii_t)=mod_indx;
                end
                
                if retrieval_trials(trialNo)==1
                    if ii_t==1
                        no_retrieval_trials=no_retrieval_trials+1;
                    end
                    retr_phase_histo(no_retrieval_trials,ii_t,1:n_phase_bins+1)=phase_histo;
                    MI_retr(no_retrieval_trials,ii_t)=mod_indx;
                    time_stamps_retr(ii_t,1:n_phase_bins+1)=phase_window/2;
                end
                
            end
        end
    end
    
end %for evNo

mean_enc_phase_histo=zeros(size(enc_phase_histo,2),size(enc_phase_histo,3));
mean_enc_phase_histo(:,:)=mean(enc_phase_histo,1);

mean_retr_phase_histo=zeros(size(retr_phase_histo,2),size(retr_phase_histo,3));
mean_retr_phase_histo(:,:)=mean(retr_phase_histo,1);

maxphase=prctile([mean_enc_phase_histo(:)' mean_retr_phase_histo(:)'],98);
minphase=prctile([mean_enc_phase_histo(:)' mean_retr_phase_histo(:)'],2);

try
    close 2
catch
end

hFig2 = figure(2);
set(hFig2, 'units','normalized','position',[.05 .05 .75 .35])


drg_pcolor(time_stamps_enc,repmat(phase,size(time_stamps_enc,1),1),mean_enc_phase_histo)
colormap jet
shading flat
caxis([minphase maxphase]);
xlabel('Time (sec)')
ylabel('Phase (deg)');
title(['Phase timecourse for encoding in ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])

try
    close 3
catch
end

hFig3 = figure(3);
set(hFig3, 'units','normalized','position',[.1 .06 .75 .35])


drg_pcolor(time_stamps_retr,repmat(phase,size(time_stamps_retr,1),1),mean_retr_phase_histo)
colormap jet
shading flat
caxis([minphase maxphase]);
xlabel('Time (sec)')
ylabel('Phase (deg)');
title(['Phase timecourse for retrieval in ' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])


try
    close 4
catch
end

hFig4 = figure(4);
set(hFig4, 'units','normalized','position',[.90 .1 .05 .3])

prain=[minphase:(maxphase-minphase)/99:maxphase];
pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
colormap jet
shading interp
ax=gca;
set(ax,'XTickLabel','')

figure(1)
subplot(3,1,2)
shadedErrorBar(time_stamps_enc(:,1),mean(MI_enc,1),std(MI_enc,0,1)/sqrt(size(MI_enc,1)))
title('MI for encoding')
xlabel('Time(sec)')
ylabel('MI')

subplot(3,1,3)
shadedErrorBar(time_stamps_retr(:,1),mean(MI_retr,1),std(MI_retr,0,1)/sqrt(size(MI_retr,1)))
title('MI for retreival')
xlabel('Time(sec)')
ylabel('MI')

pffft=1


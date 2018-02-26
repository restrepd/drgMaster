function drgMIvsHiLoPerTrial(handles)

deltaLowF_PAC=handles.deltaLowF_PAC;
deltaHighF_PAC=handles.deltaHighF_PAC;
bandwidth_lowF=handles.bandwidth_lowF;
bandwidth_highF=handles.bandwidth_highF;

%Generates a comodulogram for the encoding and retrieval segments
sessionNo=handles.sessionNo;
Fs=handles.drg.session(sessionNo).draq_p.ActualRate;
lowF1=handles.peakLowF-(handles.bandwidth_lowF/2);
lowF2=handles.peakHighF-(handles.bandwidth_lowF/2);
highF1=handles.burstLowF-(handles.bandwidth_highF/2);
highF2=handles.burstHighF-(handles.bandwidth_highF/2);
pad_time=handles.time_pad;
n_phase_bins=handles.n_phase_bins;

[perCorr, encoding_trials, retrieval_trials, encoding_this_evTypeNo,retrieval_this_evTypeNo]=drgFindEncRetr(handles);


try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.55 .27 .35 .15])

%Plot the percent correct
trials=1:length(perCorr);
plot(trials,perCorr,'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
hold on
plot(trials(encoding_trials),perCorr(encoding_trials),'ob')
plot(trials(retrieval_trials),perCorr(retrieval_trials),'or')
plot(trials(handles.trialNo),perCorr(handles.trialNo),'om','MarkerSize',10,'MarkerEdgeColor','m','MarkerFaceColor','m')
ylim([40 110]);
xlabel('Trial #')
ylabel('Percent correct')
title('% correct vs trial number, blue=encoding, red=retrieval')

%Now get the phase of gamma bursts witin theta
comodulogram=[];

evNo=handles.trialNo;
evTypeNo=handles.evTypeNo;

excludeTrial=drgExcludeTrialLFP(handles.drg,handles.peakLFPNo,handles.drg.session(sessionNo).events(evTypeNo).times(evNo),sessionNo);

if excludeTrial==0
    
    
    [LFPlow, trialNo, can_read1] = drgGetTrialLFPData(handles, handles.peakLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
    [LFPhigh, trialNo, can_read2] = drgGetTrialLFPData(handles, handles.burstLFPNo, evNo, handles.evTypeNo, handles.time_start, handles.time_end);
    
    if (can_read1==1)&(can_read2==1)
        lowFii=0;
        for lowF=lowF1:deltaLowF_PAC:lowF2
            lowFii=lowFii+1
            highFii=0;
            for highF=highF1:deltaHighF_PAC:highF2
                [meanVectorLength, meanVectorAngle, peakAngle, mod_indx, phase, phase_histo, theta_wave]=drgGetThetaAmpPhase(LFPlow, LFPhigh,Fs,lowF,lowF+bandwidth_lowF,highF,highF+bandwidth_highF,pad_time,n_phase_bins,handles.which_method);
                highFii=highFii+1;
                comodulogram(lowFii,highFii)=mod_indx;
            end
        end
    end
end
%end
%end %if eventstamps...

%IMPORTANT: pcolor does not plot the last column/row
%Because of this we have to add one row/column
lowFii=0;
for lowF=lowF1:deltaLowF_PAC:lowF2+deltaLowF_PAC
    lowFii=lowFii+1;
    highFii=0;
    for highF=highF1:deltaHighF_PAC:highF2+deltaHighF_PAC
        highFii=highFii+1;
        comLowF(lowFii,highFii)=lowF+(bandwidth_lowF/2)-(deltaLowF_PAC/2);
        comHighF(lowFii,highFii)=highF+(bandwidth_highF/2);
    end
end


try
    close 2
catch
end

hFig2 = figure(2);
set(hFig2, 'units','normalized','position',[.25 .1 .23 .8])

szcom=size(comodulogram);


max_com=max(max(comodulogram))
min_com=min(min(comodulogram))

show_com=zeros(szcom(1)+1,szcom(2)+1);
show_com(1:szcom(1),1:szcom(2))=comodulogram;

pcolor(comLowF,comHighF,show_com)
colormap jet
caxis([min_com max_com]);
shading interp
xlabel('Frequency for phase (Hz)')
ylabel('Frequency for amplitude (Hz)');
title(['Phase-amplitude MI comodulogram trial No ' num2str(handles.trialNo) ' ' drgGetEventID(handles,handles.trialNo)])


try
    close 3
catch
end

hFig3 = figure(3);
set(hFig3, 'units','normalized','position',[.49 .1 .05 .3])

prain=[min_com:(max_com-min_com)/99:max_com];
pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
colormap jet
shading interp
ax=gca;
set(ax,'XTickLabel','')
pffft=1





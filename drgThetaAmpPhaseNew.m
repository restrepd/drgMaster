function handles=drgThetaAmpPhaseNew(handles)

%Get PAC
handles=drgGetThetaAmpPhaseAllTrials(handles);

%Now display the data
if handles.displayData==1
    
    try
        close 1
    catch
    end
    
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.55 .27 .35 .15])
    
    
    %Plot the percent correct
    perCorr=handles.drgb.PAC.perCorr;
    encoding_trials=perCorr<65;
    retrieval_trials=perCorr>=80;
    trials=1:length(perCorr);
    plot(trials,perCorr,'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
    hold on
    plot(trials(encoding_trials),perCorr(encoding_trials),'ob')
    plot(trials(retrieval_trials),perCorr(retrieval_trials),'or')
    ylim([40 110]);
    xlabel('Trial No')
    ylabel('Percent correct')
    title('% correct vs trial number, blue=learning, red=retrieval')
    
    try
        close 2
    catch
    end
    
    hFig2 = figure(2);
    set(hFig2, 'units','normalized','position',[.25 .67 .23 .22])
    
    %Plot the modulation index as a function of trial number
    subplot(2,1,1)
    no_trials=handles.drgb.PAC.no_trials;
    trials=1:no_trials;
    mod_indx=handles.drgb.PAC.mod_indx;
    plot(trials,mod_indx,'o-')
    xlabel('Trial No')
    ylabel('Modulation index')
    title('Modulation index vs trial number')
    
    %Plot the vector length as a function of trial number
    subplot(2,1,2)
    meanVectorLength=handles.drgb.PAC.meanVectorLength;
    plot(trials,meanVectorLength,'o-')
    xlabel('Trial No')
    ylabel('Vector length')
    title('Vector length vs trial number')
    
    try
        close 6
    catch
    end
    
    hFig6 = figure(6);
    set(hFig6, 'units','normalized','position',[.25 .08 .23 .5])
    
    %Now plot the mean theta waveform
    subplot(4,1,1)
    phase=handles.drgb.PAC.phase;
    all_theta_wave=handles.drgb.PAC.all_theta_wave;
    shadedErrorBar(phase,mean(all_theta_wave,1),std(all_theta_wave,0,1),'-b')
    xlim([0 360])
    title('Mean low frequency waveform')
    xlabel('Degrees')
    
    %Now plot the encoding theta phase histogram for gamma
    no_encoding_trials=sum(encoding_trials);
    mean_MI_enc=mean(mod_indx(encoding_trials));
    if no_encoding_trials>0
        subplot(4,1,2)
        shadedErrorBar(phase,mean(handles.drgb.PAC.all_phase_histo(encoding_trials,:),1),std(handles.drgb.PAC.all_phase_histo(encoding_trials,:),0,1),'-b')
        xlim([0 360])
        title('Phase-amplitude coupling for encoding segment')
        ylabel('Probability')
        xlabel('Phase for low freq oscillation (deg)')
        yl_enc=ylim;
        text(20,yl_enc(1)+0.9*(yl_enc(2)-yl_enc(1)),['MI= ' num2str(mean_MI_enc)])
    end
    
    
    %Now plot the middle theta phase histogram for gamma
    mid_trials=(perCorr>=65)&(perCorr<80);
    no_mid_trials=sum(mid_trials);
    mean_MI_mid=mean(mod_indx(mid_trials));
    if no_mid_trials>0
        subplot(4,1,3)
        shadedErrorBar(phase,mean(handles.drgb.PAC.all_phase_histo(mid_trials,:),1),std(handles.drgb.PAC.all_phase_histo(mid_trials,:),0,1),'-b')
        xlim([0 360])
        title('Phase-amplitude coupling for middle segment')
        ylabel('Probability')
        xlabel('Phase for low freq oscillation (deg)')
        yl_mid=ylim;
        text(20,yl_mid(1)+0.9*(yl_mid(2)-yl_mid(1)),['MI= ' num2str(mean_MI_mid)])
    end
    
    %Now plot the retrieval theta phase histogram for gamma
    no_retrieval_trials=sum(retrieval_trials);
    mean_MI_retr=mean(mod_indx(retrieval_trials));
    if no_retrieval_trials>0
        subplot(4,1,4)
        shadedErrorBar(phase,mean(handles.drgb.PAC.all_phase_histo(retrieval_trials,:),1),std(handles.drgb.PAC.all_phase_histo(retrieval_trials,:),0,1),'-b')
        xlim([0 360])
        title('Phase-amplitude coupling for retreival segment')
        ylabel('Probability')
        xlabel('Phase for low freq oscillation (deg)')
        yl_retr=ylim;
        text(20,yl_retr(1)+0.9*(yl_retr(2)-yl_retr(1)),['MI= ' num2str(mean_MI_retr)])
    end
    
    try
        close 3
    catch
    end
    
    hFig3 = figure(3);
    set(hFig3, 'units','normalized','position',[.01 .1 .23 .8])
    
    
    drg_pcolor(repmat(phase,length(trials),1),repmat(trials',1,length(phase)),handles.drgb.PAC.all_phase_histo)
    colormap jet
    shading flat
    xlabel('Phase for low freq oscillation (deg)')
    ylabel('Trial No');
    title(['Phase-amplitude coupling per trial' handles.drg.session(1).draq_d.eventlabels{handles.evTypeNo}])
    z_per=prctile(handles.drgb.PAC.all_phase_histo(:),95);
    z_min=prctile(handles.drgb.PAC.all_phase_histo(:),5);
    caxis([z_min z_per]);
    
    try
        close 4
    catch
    end
    
    hFig4 = figure(4);
    set(hFig4, 'units','normalized','position',[.5 .6 .05 .3])
    
    
    prain=[z_min:(1.2*z_per-z_min)/99:1.2*z_per];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    colormap jet
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')
    
    try
        close 5
    catch
    end
    
    hFig5 = figure(5);
    set(hFig5, 'units','normalized','position',[.49 .05 .5 .35])
    
    meanVectorAngle=handles.drgb.PAC.meanVectorAngle;
    meanVectorAngle_enc=meanVectorAngle(encoding_trials);
    if no_encoding_trials>0
        subplot(1,3,1)
        rose(pi*meanVectorAngle_enc/180,12)
        title('Mean angle, encoding')
    end
    
    meanVectorAngle_mid=meanVectorAngle(mid_trials);
    if no_mid_trials>0
        subplot(1,3,2)
        rose(pi*meanVectorAngle_mid/180,12)
        title('Mean angle, middle')
    end
    
    meanVectorAngle_retr=meanVectorAngle(retrieval_trials);
    if no_retrieval_trials>0
        subplot(1,3,3)
        rose(pi*meanVectorAngle_retr/180,12)
        title('Mean angle, retrieval')
    end
    
end






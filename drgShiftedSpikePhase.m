function handles=drgShiftedSpikePhase(handles)

max_shift=1/mean([handles.burstLowF handles.burstHighF]);
no_bins=20;
time_shift=[];
MI_Tort=[];
meanVectorLength=[];
MRL=[];

for ii_shift=1:no_bins
    handles.time_shift=-max_shift+(ii_shift-1)*(2*max_shift/no_bins)+ (2*max_shift/no_bins)/2;
    time_shift(ii_shift)=handles.time_shift;
    handles=drgSpikePhase(handles);
    MI_Tort(ii_shift)=handles.drgb.SpikePhase.MI_Tort;
    meanVectorLength(ii_shift)=handles.drgb.SpikePhase.meanVectorLength;
    MRL(ii_shift)=handles.drgb.SpikePhase.MRL;
end

handles.time_shift=0;

if handles.displayData==1
    try
        close 1
    catch
    end
    
    %Plot the histogram
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.02 .4 .5 .5])
    
    bar(time_shift,MI_Tort)

    xlabel('Time shift (sec)')
    ylabel('MI')
    title('Modulation index vs. LFP shift')
    
     try
        close 2
    catch
    end
    
    %Plot the histogram
    hFig2 = figure(2);
    set(hFig2, 'units','normalized','position',[.02 .4 .5 .5])
    
    bar(time_shift,meanVectorLength)

    xlabel('Time shift (sec)')
    ylabel('MVL')
    title('Mean vector length vs. LFP shift')
    
        try
        close 3
    catch
    end
    
    %Plot the histogram
    hFig3 = figure(3);
    set(hFig3, 'units','normalized','position',[.02 .4 .5 .5])
    
    bar(time_shift,MRL)

    xlabel('Time shift (sec)')
    ylabel('MRL')
    title('MRL vs. LFP shift')
end
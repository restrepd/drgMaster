function handles=drgLFPCorrAllPairs(handles)

if ~isfield(handles,'drgbchoices')
    handles.drgbchoices.which_electrodes_pref=[1:4,13:16]; %Prefrontal
    handles.drgbchoices.which_electrodes_hipp=[5:12]; %Hippocampus
end

handles.displayData=0;
handles.corr_out=[];

%Between brain regions
%handles.burstLFPNo is hippocampus
%handles.peakLFPNo is prefrontal
handles.corr_out.no_between_pairs=0;
for burstLFPNo=handles.drgbchoices.which_electrodes_hipp
    for peakLFPNo=handles.drgbchoices.which_electrodes_pref
        handles.burstLFPNo=burstLFPNo;
        handles.peakLFPNo=peakLFPNo;
        handles=drgLFPCorrTrialRange(handles);
        handles.corr_out.no_between_pairs=handles.corr_out.no_between_pairs+1;
        handldes.corr_out.between.LFPcorr(handles.corr_out.no_between_pairs).burstLFPNo=handles.burstLFPNo;
        handldes.corr_out.between.LFPcorr(handles.corr_out.no_between_pairs).peakLFPNo=handles.peakLFPNo;
        handldes.corr_out.between.LFPcorr(handles.corr_out.no_between_pairs).LFPcorr=handles.drgb.LFPcorr;
        fprintf(1, ['Hippocampus electrode %d, prefrontal electrode %d\n'], handles.burstLFPNo,handles.peakLFPNo);
    end
end

%Within hippocampus
%handles.burstLFPNo is hippocampus
%handles.peakLFPNo is prefrontal
handles.corr_out.no_hipp_pairs=0;
for burstLFPii=1:length(handles.drgbchoices.which_electrodes_hipp)
    for peakLFPii=burstLFPii+1:length(handles.drgbchoices.which_electrodes_hipp)
        handles.burstLFPNo=handles.drgbchoices.which_electrodes_hipp(burstLFPii);
        handles.peakLFPNo=handles.drgbchoices.which_electrodes_hipp(peakLFPii);
        handles=drgLFPCorrTrialRange(handles);
        handles.corr_out.no_hipp_pairs=handles.corr_out.no_hipp_pairs+1;
        handldes.corr_out.hipp.LFPcorr(handles.corr_out.no_hipp_pairs).burstLFPNo=handles.burstLFPNo;
        handldes.corr_out.hipp.LFPcorr(handles.corr_out.no_hipp_pairs).peakLFPNo=handles.peakLFPNo;
        handldes.corr_out.hipp.LFPcorr(handles.corr_out.no_hipp_pairs).LFPcorr=handles.drgb.LFPcorr;
        fprintf(1, ['Hippocampus electrode %d, hippocampus electrode %d\n'], handles.burstLFPNo,handles.peakLFPNo);
    end
end

%Within prefrontal
%handles.burstLFPNo is hippocampus
%handles.peakLFPNo is prefrontal
handles.corr_out.no_pref_pairs=0;
for burstLFPii=1:length(handles.drgbchoices.which_electrodes_pref)
    for peakLFPii=burstLFPii+1:length(handles.drgbchoices.which_electrodes_pref)
        handles.burstLFPNo=handles.drgbchoices.which_electrodes_pref(burstLFPii);
        handles.peakLFPNo=handles.drgbchoices.which_electrodes_pref(peakLFPii);
        handles=drgLFPCorrTrialRange(handles);
        handles.corr_out.no_pref_pairs=handles.corr_out.no_pref_pairs+1;
        handldes.corr_out.pref.LFPcorr(handles.corr_out.no_pref_pairs).burstLFPNo=handles.burstLFPNo;
        handldes.corr_out.pref.LFPcorr(handles.corr_out.no_pref_pairs).peakLFPNo=handles.peakLFPNo;
        handldes.corr_out.pref.LFPcorr(handles.corr_out.no_pref_pairs).LFPcorr=handles.drgb.LFPcorr;
        fprintf(1, ['Prefrontal electrode %d, prefrontal electrode %d\n'], handles.burstLFPNo,handles.peakLFPNo);
    end
end

%Plot the 

%Between
between_mean_max_rho_t_lag=[];
for ii_pairs=1:handles.corr_out.no_between_pairs
    between_mean_max_rho_t_lag(ii_pairs)=mean(hanldes.corr_out.between.LFPcorr(ii_pairs).LFPcorr.max_rho_t_lag);
end

%Hippocampus
hipp_mean_max_rho_t_lag=[];
for ii_pairs=1:handles.corr_out.no_hipp_pairs
    hipp_mean_max_rho_t_lag(ii_pairs)=mean(handldes.corr_out.hipp.LFPcorr(ii_pairs).LFPcorr.max_rho_t_lag);
end

%Prefrontal
pref_mean_max_rho_t_lag=[];
for ii_pairs=1:handles.corr_out.no_pref_pairs
    pref_mean_max_rho_t_lag(ii_pairs)=mean(handldes.corr_out.pref.LFPcorr(ii_pairs).LFPcorr.max_rho_t_lag);
end

try
    close 4
catch
end

hFig4 = figure(4);
set(hFig4, 'units','normalized','position',[.27 .12 .27 .27])




[f_mrtlag,x_mrtlag] = drg_ecdf(between_mean_max_rho_t_lag);
plot(x_mrtlag,f_mrtlag,'-r','LineWidth',3);

hold on

[f_mrtlag,x_mrtlag] = drg_ecdf(hipp_mean_max_rho_t_lag);
plot(x_mrtlag,f_mrtlag,'-b','LineWidth',3);

[f_mrtlag,x_mrtlag] = drg_ecdf(pref_mean_max_rho_t_lag);
plot(x_mrtlag,f_mrtlag,'-k','LineWidth',3);

plot([0 0],[0 1],'-k')
title('max rho lag time')
xlabel('lag time (sec)')
ylabel('Probability')
    
pffft=1;
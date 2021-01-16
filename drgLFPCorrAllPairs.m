function handles=drgLFPCorrAllPairs(handles)

if ~isfield(handles,'drgbchoices')
    handles.drgbchoices.which_electrodes_pref=[1:4,13:16]; %Prefrontal
    handles.drgbchoices.which_electrodes_hipp=[5:12]; %Hippocampus
%     
%     handles.drgbchoices.which_electrodes_pref=[1:8]; %Prefrontal
%     handles.drgbchoices.which_electrodes_hipp=[9:16]; %Hippocampus
    
%     handles.drgbchoices.which_electrodes_hipp=[1 3 5 7 9 11 13 15]; %Prefrontal
%     handles.drgbchoices.which_electrodes_pref=[2 4 6 8 10 12 14 16]; %Hippocampus
end

handles.displayData=0;
handles.corr_out=[];

gcp

%Between brain regions
%handles.burstLFPNo is hippocampus
%handles.peakLFPNo is prefrontal

handles.corr_out.no_between_pairs=0;
par_out_LFPcorr=[];
for burstLFPNo=handles.drgbchoices.which_electrodes_hipp
    for peakLFPNo=handles.drgbchoices.which_electrodes_pref
        handles.corr_out.no_between_pairs=handles.corr_out.no_between_pairs+1;
        handles.corr_out.between.burstLFPNo(handles.corr_out.no_between_pairs)=burstLFPNo;
        handles.corr_out.between.peakLFPNo(handles.corr_out.no_between_pairs)=peakLFPNo;
        par_out_LFPcorr(handles.corr_out.no_between_pairs).LFPcorr=[];
    end
end

no_between_pairs=handles.corr_out.no_between_pairs;
parfor ii_between_pairs=1:no_between_pairs
    handlespf=struct();
    handlespf=handles;
    handlespf.burstLFPNo=handlespf.corr_out.between.burstLFPNo(ii_between_pairs);
    handlespf.peakLFPNo=handlespf.corr_out.between.peakLFPNo(ii_between_pairs);
    handlespf=drgLFPCorrTrialRange(handlespf);
    par_out_LFPcorr(ii_between_pairs).LFPcorr=handlespf.drgb.LFPcorr;
    fprintf(1, ['Hippocampus electrode %d, prefrontal electrode %d\n'], handlespf.burstLFPNo,handlespf.peakLFPNo);
end

%Save the output to handles
for ii_between_pairs=1:no_between_pairs
    handles.corr_out.between.LFPcorr(ii_between_pairs).LFPcorr=par_out_LFPcorr(ii_between_pairs).LFPcorr;
end


%Within hippocampus
handles.corr_out.no_hipp_pairs=0;
par_out_LFPcorr=[];
for burstLFPii=1:length(handles.drgbchoices.which_electrodes_hipp)
    for peakLFPii=burstLFPii+1:length(handles.drgbchoices.which_electrodes_hipp)
        handles.corr_out.no_hipp_pairs=handles.corr_out.no_hipp_pairs+1;
        handles.corr_out.hipp.burstLFPNo(handles.corr_out.no_hipp_pairs)=handles.drgbchoices.which_electrodes_hipp(burstLFPii);
        handles.corr_out.hipp.peakLFPNo(handles.corr_out.no_hipp_pairs)=handles.drgbchoices.which_electrodes_hipp(peakLFPii);
        par_out_LFPcorr(handles.corr_out.no_hipp_pairs).LFPcorr=[];
    end
end

no_hipp_pairs=handles.corr_out.no_hipp_pairs;
parfor ii_hipp_pairs=1:no_hipp_pairs
    handlespf=struct();
    handlespf=handles;
    handlespf.burstLFPNo=handlespf.corr_out.hipp.burstLFPNo(ii_hipp_pairs);
    handlespf.peakLFPNo=handlespf.corr_out.hipp.peakLFPNo(ii_hipp_pairs);
    handlespf=drgLFPCorrTrialRange(handlespf);
    par_out_LFPcorr(ii_hipp_pairs).LFPcorr=handlespf.drgb.LFPcorr;
    fprintf(1, ['Hippocampus electrode %d, hippocampus electrode %d\n'], handlespf.burstLFPNo,handlespf.peakLFPNo);
end

%Save the output to handles
for ii_hipp_pairs=1:no_hipp_pairs
    handles.corr_out.hipp.LFPcorr(ii_hipp_pairs).LFPcorr=par_out_LFPcorr(ii_hipp_pairs).LFPcorr;
end



%Within prefrontal
handles.corr_out.no_pref_pairs=0;
par_out_LFPcorr=[];
for burstLFPii=1:length(handles.drgbchoices.which_electrodes_pref)
    for peakLFPii=burstLFPii+1:length(handles.drgbchoices.which_electrodes_pref)
        handles.corr_out.no_pref_pairs=handles.corr_out.no_pref_pairs+1;
        handles.corr_out.pref.burstLFPNo(handles.corr_out.no_pref_pairs)=handles.drgbchoices.which_electrodes_pref(burstLFPii);
        handles.corr_out.pref.peakLFPNo(handles.corr_out.no_pref_pairs)=handles.drgbchoices.which_electrodes_pref(peakLFPii);
        par_out_LFPcorr(handles.corr_out.no_pref_pairs).LFPcorr=[];
    end
end

no_pref_pairs=handles.corr_out.no_pref_pairs;
parfor ii_pref_pairs=1:no_pref_pairs
    handlespf=struct();
    handlespf=handles;
    handlespf.burstLFPNo=handlespf.corr_out.pref.burstLFPNo(ii_pref_pairs);
    handlespf.peakLFPNo=handlespf.corr_out.pref.peakLFPNo(ii_pref_pairs);
    handlespf=drgLFPCorrTrialRange(handlespf);
    par_out_LFPcorr(ii_pref_pairs).LFPcorr=handlespf.drgb.LFPcorr;
    fprintf(1, ['Prefrontal electrode %d, Prefrontal electrode %d\n'], handlespf.burstLFPNo,handlespf.peakLFPNo);
end

%Save the output to handles
for ii_pref_pairs=1:no_pref_pairs
    handles.corr_out.pref.LFPcorr(ii_pref_pairs).LFPcorr=par_out_LFPcorr(ii_pref_pairs).LFPcorr;
end

%Plot the histograms

%Between
between_mean_max_rho_t_lag=[];
between_max_rho_t_lag=[];
for ii_pairs=1:handles.corr_out.no_between_pairs
    between_mean_max_rho_t_lag(ii_pairs)=mean(handles.corr_out.between.LFPcorr(ii_pairs).LFPcorr.max_rho_t_lag);
    between_max_rho_t_lag=[between_max_rho_t_lag handles.corr_out.between.LFPcorr(ii_pairs).LFPcorr.max_rho_t_lag];
end

%Hippocampus
hipp_mean_max_rho_t_lag=[];
for ii_pairs=1:handles.corr_out.no_hipp_pairs
    hipp_mean_max_rho_t_lag(ii_pairs)=mean(handles.corr_out.hipp.LFPcorr(ii_pairs).LFPcorr.max_rho_t_lag);
end

%Prefrontal
pref_mean_max_rho_t_lag=[];
for ii_pairs=1:handles.corr_out.no_pref_pairs
    pref_mean_max_rho_t_lag(ii_pairs)=mean(handles.corr_out.pref.LFPcorr(ii_pairs).LFPcorr.max_rho_t_lag);
end

try
    close 3
catch
end

hFig3 = figure(3);
set(hFig3, 'units','normalized','position',[.27 .12 .27 .27])


edges=[-0.02:0.001:0.02];
histogram(between_mean_max_rho_t_lag,edges)

hold on
thisylm=ylim;
plot([mean(between_mean_max_rho_t_lag), mean(between_mean_max_rho_t_lag)],[thisylm(1) thisylm(2)],'-k')
this_ylim=ylim;
text(0.005,this_ylim(1)+0.8*(this_ylim(2)-this_ylim(1)),['Mean lag= ', num2str(1000*mean(between_mean_max_rho_t_lag)) ' ms'])
p_between=ranksum(between_mean_max_rho_t_lag,zeros(1,length(between_mean_max_rho_t_lag)))

title('Between')

try
    close 4
catch
end

hFig4 = figure(4);
set(hFig4, 'units','normalized','position',[.27 .12 .27 .27])


edges=[-0.02:0.001:0.02];
histogram(hipp_mean_max_rho_t_lag,edges)

hold on
thisylm=ylim;
plot([mean(hipp_mean_max_rho_t_lag), mean(hipp_mean_max_rho_t_lag)],[thisylm(1) thisylm(2)],'-k')
this_ylim=ylim;
text(0.005,this_ylim(1)+0.8*(this_ylim(2)-this_ylim(1)),['Mean lag= ', num2str(1000*mean(hipp_mean_max_rho_t_lag)) ' ms'])
p_hippocampus=ranksum(hipp_mean_max_rho_t_lag,zeros(1,length(hipp_mean_max_rho_t_lag)))

title('Hippocampus')

try
    close 5
catch
end

hFig5 = figure(5);
set(hFig5, 'units','normalized','position',[.27 .12 .27 .27])

edges=[-0.02:0.001:0.02];
histogram(pref_mean_max_rho_t_lag,edges)

hold on
thisylm=ylim;
plot([mean(pref_mean_max_rho_t_lag), mean(pref_mean_max_rho_t_lag)],[thisylm(1) thisylm(2)],'-k')
this_ylim=ylim;
text(0.005,this_ylim(1)+0.8*(this_ylim(2)-this_ylim(1)),['Mean lag= ', num2str(1000*mean(pref_mean_max_rho_t_lag)) ' ms'])
title('Prefrontal')

p_prefrontal=ranksum(pref_mean_max_rho_t_lag,zeros(1,length(pref_mean_max_rho_t_lag)))
    



    
%Between
between_mean_max_rho_t_lag_sh=[];
between_max_rho_t_lag_sh=[];
for ii_pairs=1:handles.corr_out.no_between_pairs
    between_mean_max_rho_t_lag_sh(ii_pairs)=mean(handles.corr_out.between.LFPcorr(ii_pairs).LFPcorr.max_rho_t_lag_sh);
    between_max_rho_t_lag_sh=[between_max_rho_t_lag_sh handles.corr_out.between.LFPcorr(ii_pairs).LFPcorr.max_rho_t_lag_sh];
end

%Hippocampus
hipp_mean_max_rho_t_lag_sh=[];
for ii_pairs=1:handles.corr_out.no_hipp_pairs
    hipp_mean_max_rho_t_lag_sh(ii_pairs)=mean(handles.corr_out.hipp.LFPcorr(ii_pairs).LFPcorr.max_rho_t_lag_sh);
end

%Prefrontal
pref_mean_max_rho_t_lag_sh=[];
for ii_pairs=1:handles.corr_out.no_pref_pairs
    pref_mean_max_rho_t_lag_sh(ii_pairs)=mean(handles.corr_out.pref.LFPcorr(ii_pairs).LFPcorr.max_rho_t_lag_sh);
end

try
    close 7
catch
end

hFig7 = figure(7);
set(hFig7, 'units','normalized','position',[.27 .12 .27 .27])

edges=[-0.02:0.001:0.02];
histogram(between_mean_max_rho_t_lag_sh,edges)

hold on
thisylm=ylim;
plot([mean(between_mean_max_rho_t_lag_sh), mean(between_mean_max_rho_t_lag_sh)],[thisylm(1) thisylm(2)],'-k')
this_ylim=ylim;
text(0.005,this_ylim(1)+0.8*(this_ylim(2)-this_ylim(1)),['Mean lag= ', num2str(1000*mean(between_mean_max_rho_t_lag_sh)) ' ms'])
title('Between shuffled')
p_between_shuffled=ranksum(between_mean_max_rho_t_lag_sh,zeros(1,length(between_mean_max_rho_t_lag_sh)))
  

try
    close 8
catch
end

hFig8 = figure(8);
set(hFig8, 'units','normalized','position',[.27 .12 .27 .27])

edges=[-0.02:0.001:0.02];
histogram(hipp_mean_max_rho_t_lag_sh,edges)

hold on
thisylm=ylim;
plot([mean(hipp_mean_max_rho_t_lag_sh), mean(hipp_mean_max_rho_t_lag_sh)],[thisylm(1) thisylm(2)],'-k')
this_ylim=ylim;
text(0.005,this_ylim(1)+0.8*(this_ylim(2)-this_ylim(1)),['Mean lag= ', num2str(1000*mean(hipp_mean_max_rho_t_lag_sh)) ' ms'])
title('Hippocampus shuffled')
p_hippocampus_shuffled=ranksum(hipp_mean_max_rho_t_lag_sh,zeros(1,length(hipp_mean_max_rho_t_lag_sh)))

try
    close 9
catch
end

hFig9 = figure(9);
set(hFig9, 'units','normalized','position',[.27 .12 .27 .27])

edges=[-0.02:0.001:0.02];
histogram(pref_mean_max_rho_t_lag_sh,edges)

hold on
thisylm=ylim;
plot([mean(pref_mean_max_rho_t_lag_sh), mean(pref_mean_max_rho_t_lag_sh)],[thisylm(1) thisylm(2)],'-k')
this_ylim=ylim;
text(0.005,this_ylim(1)+0.8*(this_ylim(2)-this_ylim(1)),['Mean lag= ', num2str(1000*mean(pref_mean_max_rho_t_lag_sh)) ' ms'])
title('Prefrontal shuffled')
p_prefrontal_shuffled=ranksum(pref_mean_max_rho_t_lag_sh,zeros(1,length(pref_mean_max_rho_t_lag_sh)))
pffft=1;
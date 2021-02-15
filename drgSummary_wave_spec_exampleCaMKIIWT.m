function drgSummary_wave_spec_exampleCaMKIIWT
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'drgbDiscPar' file listing 'Discriminant_*.mat' output files
%from drgLFPDiscriminantBatch
%
%Performs summary analises for LDA and  PAC analysis performed with case 19
%of drgAnalysisBatchLFP

warning('off')

close all
clear all

freq=[4:(95-4)/100:95];

PACnames{1}='Beta';
PACnames{2}='Low gamma';
PACnames{3}='High gamma';

prof_naive_leg{1}='Proficient';
prof_naive_leg{2}='Naive';

group_legend{1}='WT';
group_legend{2}='Het';
group_legend{3}='KO';

peak_legend{1}='Peak';
peak_legend{2}='Trough';

evTypeLabels{1}='S+';
evTypeLabels{2}='S-';

%Location of files
% PathName='C:\Users\Diego Restrepo\OneDrive - The University of Colorado Denver\CaMKII Paper\Figure 3 PRP\';
PathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Figure_3_PRP/';

% handles_pars.which_electrodes=[5:12]; %Hippocampus
% handles_pars.which_electrodes=[1:4,13:16]; %Prefrontal

%The numbers above are incorrect, here is take 2
% handles.drgbchoices.which_electrodes=[1:8]; %Hippocampus 2
% handles.drgbchoices.which_electrodes=[9:16]; %Prefrontal 2



figNo=0;





%Now plot for the hippocampus the dependence of delta LFP to freq

glm_dB_hipp=[];
glm_ii_hipp=0;

glm_dB_pre=[];
glm_ii_pre=0;

glm_dB_both=[];
glm_ii_both=0;

maxdB=-1000000;
mindB=1000000;

%Let's do hippocampus
for peak=1:2
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .35 .3])
    hold on
    
    these_dBWB=[];
    
    for evNo=2:-1:1
        
        for per_ii=2:-1:1
            
            %Load the data
            if peak==1
                if evNo==1
                    if per_ii==1
                        %S+ proficient
                        load([PathName '6251853931_elec5_splus_proficient.mat'])
                        these_dBWB=peakPowerSpectrum;
                    else
                        %S+ naive
                        load([PathName '6151853931_elec5_splus_naive.mat'])
                        these_dBWB=peakPowerSpectrum;
                    end
                else
                     if per_ii==1
                        %S- proficient
                        load([PathName '6251853931_elec5_sminus_proficient.mat'])
                        these_dBWB=peakPowerSpectrum;
                    else
                        %S- naive
                        load([PathName '6151853931_elec5_sminus_naive.mat'])
                        these_dBWB=peakPowerSpectrum;
                    end
                end
            else
                if evNo==1
                    if per_ii==1
                        %S+ proficient
                        load([PathName '6251853931_elec5_splus_proficient.mat'])
                        these_dBWB=troughPowerSpectrum;
                    else
                        %S+ naive
                        load([PathName '6151853931_elec5_splus_naive.mat'])
                        these_dBWB=troughPowerSpectrum;
                    end
                else
                     if per_ii==1
                        %S- proficient
                        load([PathName '6251853931_elec5_sminus_proficient.mat'])
                        these_dBWB=troughPowerSpectrum;
                    else
                        %S- naive
                        load([PathName '6151853931_elec5_sminus_naive.mat'])
                        these_dBWB=troughPowerSpectrum;
                    end
                end
            end
            
            
            
            this_mean_dbWB=zeros(1,length(freq));
            this_mean_dbWB(1,:)=mean(these_dBWB,1);
            
            CI=[];
            CI = bootci(1000, {@mean, these_dBWB})';
            CI(:,1)= this_mean_dbWB'-CI(:,1);
            CI(:,2)=CI(:,2)- this_mean_dbWB';
            
            if evNo==1
                if per_ii==1
                    %S+ Proficient
                    [hlCR, hpCR] = boundedline(freq',this_mean_dbWB', CI, 'cmap',[158/255 31/255 99/255]);
                else
                    %S+ Naive
                    [hlCR, hpCR] = boundedline(freq',this_mean_dbWB', CI, 'cmap',[238/255 111/255 179/255]);
                end
            else
                if per_ii==1
                    %S- Proficient
                    [hlCR, hpCR] = boundedline(freq',this_mean_dbWB', CI, 'cmap',[0 114/255 178/255]);
                else
                    %S- naive
                    [hlCR, hpCR] = boundedline(freq',this_mean_dbWB', CI, 'cmap',[80/255 194/255 255/255]);
                end
            end
            
            no_dbWB=size(these_dBWB,1);
            for iif=1:length(freq)
                glm_dB_hipp.data(glm_ii_hipp+1:glm_ii_hipp+no_dbWB)=these_dBWB(:,iif);
                glm_dB_hipp.freq(glm_ii_hipp+1:glm_ii_hipp+no_dbWB)=freq(iif);
                glm_dB_hipp.per_ii(glm_ii_hipp+1:glm_ii_hipp+no_dbWB)=per_ii;
                glm_dB_hipp.evNo(glm_ii_hipp+1:glm_ii_hipp+no_dbWB)=evNo;
                glm_dB_hipp.peak(glm_ii_hipp+1:glm_ii_hipp+no_dbWB)=peak;
                glm_ii_hipp=glm_ii_hipp+no_dbWB;
                
                glm_dB_both.data(glm_ii_both+1:glm_ii_both+no_dbWB)=these_dBWB(:,iif);
                glm_dB_both.freq(glm_ii_both+1:glm_ii_both+no_dbWB)=freq(iif);
                glm_dB_both.per_ii(glm_ii_both+1:glm_ii_both+no_dbWB)=per_ii;
                glm_dB_both.evNo(glm_ii_both+1:glm_ii_both+no_dbWB)=evNo;
                glm_dB_both.peak(glm_ii_both+1:glm_ii_both+no_dbWB)=peak;
                glm_dB_both.tissue(glm_ii_both+1:glm_ii_both+no_dbWB)=1;
                glm_ii_both=glm_ii_both+no_dbWB;
            end
            
            maxdB=max([maxdB max(this_mean_dbWB)]);
            mindB=min([mindB min(this_mean_dbWB)]);
            
            
        end
    end
    
    title([peak_legend{peak} ' wavelet LFP power vs freq for hippocampus'])
    
    
    %Proficient/Naive annotations
    annotation('textbox',[0.75 0.85 0.3 0.1],'String','S+ Proficient','FitBoxToText','on','Color',[158/255 31/255 99/255],'LineStyle','none');
    annotation('textbox',[0.75 0.80 0.3 0.1],'String','S+ Naive','FitBoxToText','on','Color',[238/255 111/255 179/255],'LineStyle','none');
    annotation('textbox',[0.75 0.75 0.3 0.1],'String','S- Proficient','FitBoxToText','on','Color',[[0 114/255 178/255]],'LineStyle','none');
    annotation('textbox',[0.75 0.70 0.3 0.1],'String','S- Naive','FitBoxToText','on','Color',[80/255 194/255 255/255],'LineStyle','none');
    
    
    xlim([4 95])
    
    xlabel('freq (Hz)')
    
    ylabel('dB')
    
end

%Let's do prefrontal
for peak=1:2
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .35 .3])
    hold on
    
    these_dBWB=[];
    
    for evNo=2:-1:1
        
        for per_ii=2:-1:1
            
            %Load the data
            if peak==1
                if evNo==1
                    if per_ii==1
                        %S+ proficient
                        load([PathName '6251853931_elec12_splus_proficient.mat'])
                        these_dBWB=peakPowerSpectrum;
                    else
                        %S+ naive
                        load([PathName '6151853931_elec12_splus_naive.mat'])
                        these_dBWB=peakPowerSpectrum;
                    end
                else
                     if per_ii==1
                        %S- proficient
                        load([PathName '6251853931_elec12_sminus_proficient.mat'])
                        these_dBWB=peakPowerSpectrum;
                    else
                        %S- naive
                        load([PathName '6151853931_elec12_sminus_naive.mat'])
                        these_dBWB=peakPowerSpectrum;
                    end
                end
            else
                if evNo==1
                    if per_ii==1
                        %S+ proficient
                        load([PathName '6251853931_elec12_splus_proficient.mat'])
                        these_dBWB=troughPowerSpectrum;
                    else
                        %S+ naive
                        load([PathName '6151853931_elec12_splus_naive.mat'])
                        these_dBWB=troughPowerSpectrum;
                    end
                else
                     if per_ii==1
                        %S- proficient
                        load([PathName '6251853931_elec12_sminus_proficient.mat'])
                        these_dBWB=troughPowerSpectrum;
                    else
                        %S- naive
                        load([PathName '6151853931_elec12_sminus_naive.mat'])
                        these_dBWB=troughPowerSpectrum;
                    end
                end
            end
            
            
            
            this_mean_dbWB=zeros(1,length(freq));
            this_mean_dbWB(1,:)=mean(these_dBWB,1);
            
            CI=[];
            CI = bootci(1000, {@mean, these_dBWB})';
            CI(:,1)= this_mean_dbWB'-CI(:,1);
            CI(:,2)=CI(:,2)- this_mean_dbWB';
            
            if evNo==1
                if per_ii==1
                    %S+ Proficient
                    [hlCR, hpCR] = boundedline(freq',this_mean_dbWB', CI, 'cmap',[158/255 31/255 99/255]);
                else
                    %S+ Naive
                    [hlCR, hpCR] = boundedline(freq',this_mean_dbWB', CI, 'cmap',[238/255 111/255 179/255]);
                end
            else
                if per_ii==1
                    %S- Proficient
                    [hlCR, hpCR] = boundedline(freq',this_mean_dbWB', CI, 'cmap',[0 114/255 178/255]);
                else
                    %S- naive
                    [hlCR, hpCR] = boundedline(freq',this_mean_dbWB', CI, 'cmap',[80/255 194/255 255/255]);
                end
            end
            
            no_dbWB=size(these_dBWB,1);
            for iif=1:length(freq)
                glm_dB_pre.data(glm_ii_pre+1:glm_ii_pre+no_dbWB)=these_dBWB(:,iif);
                glm_dB_pre.freq(glm_ii_pre+1:glm_ii_pre+no_dbWB)=freq(iif);
                glm_dB_pre.per_ii(glm_ii_pre+1:glm_ii_pre+no_dbWB)=per_ii;
                glm_dB_pre.evNo(glm_ii_pre+1:glm_ii_pre+no_dbWB)=evNo;
                glm_dB_pre.peak(glm_ii_pre+1:glm_ii_pre+no_dbWB)=peak;
                glm_ii_pre=glm_ii_pre+no_dbWB;
                
                glm_dB_both.data(glm_ii_both+1:glm_ii_both+no_dbWB)=these_dBWB(:,iif);
                glm_dB_both.freq(glm_ii_both+1:glm_ii_both+no_dbWB)=freq(iif);
                glm_dB_both.per_ii(glm_ii_both+1:glm_ii_both+no_dbWB)=per_ii;
                glm_dB_both.evNo(glm_ii_both+1:glm_ii_both+no_dbWB)=evNo;
                glm_dB_both.peak(glm_ii_both+1:glm_ii_both+no_dbWB)=peak;
                glm_dB_both.tissue(glm_ii_both+1:glm_ii_both+no_dbWB)=2;
                glm_ii_both=glm_ii_both+no_dbWB;
            end
            
            maxdB=max([maxdB max(this_mean_dbWB)]);
            mindB=min([mindB min(this_mean_dbWB)]);
            
            
        end
    end
    
    title([peak_legend{peak} ' wavelet LFP power vs freq for prefrontal'])
    
    
    %Proficient/Naive annotations
    annotation('textbox',[0.75 0.85 0.3 0.1],'String','S+ Proficient','FitBoxToText','on','Color',[158/255 31/255 99/255],'LineStyle','none');
    annotation('textbox',[0.75 0.80 0.3 0.1],'String','S+ Naive','FitBoxToText','on','Color',[238/255 111/255 179/255],'LineStyle','none');
    annotation('textbox',[0.75 0.75 0.3 0.1],'String','S- Proficient','FitBoxToText','on','Color',[[0 114/255 178/255]],'LineStyle','none');
    annotation('textbox',[0.75 0.70 0.3 0.1],'String','S- Naive','FitBoxToText','on','Color',[80/255 194/255 255/255],'LineStyle','none');
    
    
    xlim([4 95])
    
    xlabel('freq (Hz)')
    
    ylabel('dB')
    
end

for fn=1:figNo
    figure(fn)
    ylim([mindB-0.05*(maxdB-mindB) maxdB+0.05*(maxdB-mindB)])
end

%Perform the glm
fprintf(1, ['\n\nglm for delta dB calculated for hippocampus\n'])
tbl = table(glm_dB_hipp.data',glm_dB_hipp.freq',glm_dB_hipp.evNo',glm_dB_hipp.per_ii',glm_dB_hipp.peak',...
    'VariableNames',{'dB','frequency','spm','naive_proficient','peak'});
mdl = fitglm(tbl,'dB~frequency+spm+naive_proficient+peak+frequency*spm*naive_proficient*peak'...
    ,'CategoricalVars',[3,4,5])

%Perform the glm
fprintf(1, ['\n\nglm for delta dB calculated for prefrontal\n'])
tbl = table(glm_dB_pre.data',glm_dB_pre.freq',glm_dB_pre.evNo',glm_dB_pre.per_ii',glm_dB_pre.peak',...
    'VariableNames',{'dB','frequency','spm','naive_proficient','peak'});
mdl = fitglm(tbl,'dB~frequency+spm+naive_proficient+peak+frequency*spm*naive_proficient*peak'...
    ,'CategoricalVars',[3,4,5])


%Perform the glm mi for both brain regions
fprintf(1, ['\n\nglm for delta dB for prefrontal and hippocampus\n'])
tbl = table(glm_dB_both.data',glm_dB_both.freq',glm_dB_both.evNo',glm_dB_both.per_ii',glm_dB_both.peak',glm_dB_both.tissue',...
    'VariableNames',{'dB','frequency','spm','naive_proficient','peak','tissue'});
mdl = fitglm(tbl,'dB~frequency+spm+naive_proficient+peak+tissue+frequency*spm*naive_proficient*peak*tissue'...
    ,'CategoricalVars',[3,4,5,6])



pffft=1

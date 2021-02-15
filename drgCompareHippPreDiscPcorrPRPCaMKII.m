function drgCompareHippPreDiscPcorrPRPCaMKII
%Analyzes the linear discriminant analysis performed by drgAnalyzeLFPDiscriminantBatchCaMKIIGrp
%Takes as a in input the 'drgbDiscPar' file listing 'Discriminant_*.mat' output files
%from drgLFPDiscriminantBatch
%
%Performs summary analises for LDA and  PAC analysis performed with case 19
%of drgAnalysisBatchLFP

warning('off')

close all
clear all


figNo=0;

PACnames{1}='Beta';
PACnames{2}='gamma';


prof_naive_leg{1}='Proficient';
prof_naive_leg{2}='Naive';

fwd_rev_leg{1}='Forward';
fwd_rev_leg{2}='Reverse';

group_legend{1}='WT';
group_legend{2}='Het';
group_legend{3}='KO';

evTypeLabels{1}='S+';
evTypeLabels{2}='S-';

peak_label{1}='Trough';
peak_label{2}='Peak';


%Location of files
% hippPathName='F:\Datos summary CaMKII111720\PRP drgAnalysisBatchLFPCaMKII case 24 output for summary\';
hippPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Discriminant/';

%Hippocampus
hippFileName='pcorr_Discriminant_CaMKIIpzz1paea_disc_01302021_hippo2.mat';

% prePathName='F:\Datos summary CaMKII111720\PRP drgAnalysisBatchLFPCaMKII case 24 output for summary\';
prePathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Discriminant/';

%Prefrontal
preFileName='pcorr_Discriminant_CaMKIIpzz1paea_disc_02012021_pref2.mat';


%Now process the hippocampus

%Load data


load([hippPathName hippFileName])
hippo_pcorr_out=pcorr_out;

load([prePathName preFileName])
pre_pcorr_out=pcorr_out;

for PACii=1:2
    
    glm_pcorr_dt_peak=[];
    glm_ii_pcdt_p=0;
    glm_pcorr_dt_trough=[];
    glm_ii_pcdt_t=0;
    
    p_corr_dt_peak_stats=[];
    ii_p_stats=0;
    p_corr_dt_trough_stats=[];
    ii_t_stats=0;
    
    %Bar graph plot for peak percent correct
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
    hold on
    
    
    bar_offset = 0;
    
    per_ii=1;  %We only do proficient
    
    for grNo=1:3
        
        %Hippocampus
        bar_offset = bar_offset +1;
        
        
        switch grNo
            case 1
                bar(bar_offset,mean(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data),'g','LineWidth', 3,'EdgeColor','none')
            case 2
                bar(bar_offset,mean(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data),'b','LineWidth', 3,'EdgeColor','none')
            case 3
                bar(bar_offset,mean(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data),'y','LineWidth', 3,'EdgeColor','none')
        end
        
        
        CI = bootci(1000, {@mean, hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data},'type','cper');
        plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        plot(bar_offset*ones(1,length(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data)),hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        
        data=hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data;
        glm_pcorr_dt_peak.data(glm_ii_pcdt_p+1:glm_ii_pcdt_p+length(data))=data;
        glm_pcorr_dt_peak.group(glm_ii_pcdt_p+1:glm_ii_pcdt_p+length(data))=grNo;
        glm_pcorr_dt_peak.brainArea(glm_ii_pcdt_p+1:glm_ii_pcdt_p+length(data))=1;
        glm_ii_pcdt_p=glm_ii_pcdt_p+length(data);
        
        ii_p_stats=ii_p_stats+1;
        p_corr_dt_peak_stats(ii_p_stats).data=data;
        p_corr_dt_peak_stats(ii_p_stats).description=[group_legend{grNo} ' ' ...
            'Hippocampus'];
        
        
        %Prefrontal
        bar_offset = bar_offset +1;
        
        
        switch grNo
            case 1
                bar(bar_offset,mean(pre_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data),'g','LineWidth', 3,'EdgeColor','none')
            case 2
                bar(bar_offset,mean(pre_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data),'b','LineWidth', 3,'EdgeColor','none')
            case 3
                bar(bar_offset,mean(pre_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data),'y','LineWidth', 3,'EdgeColor','none')
        end
        
        
        CI = bootci(1000, {@mean, pre_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data},'type','cper');
        plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        plot(bar_offset*ones(1,length(pre_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data)),pre_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        
        bar_offset = bar_offset + 1;
        
        data=pre_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).data;
        glm_pcorr_dt_peak.data(glm_ii_pcdt_p+1:glm_ii_pcdt_p+length(data))=data;
        glm_pcorr_dt_peak.group(glm_ii_pcdt_p+1:glm_ii_pcdt_p+length(data))=grNo;
        glm_pcorr_dt_peak.brainArea(glm_ii_pcdt_p+1:glm_ii_pcdt_p+length(data))=2;
        glm_ii_pcdt_p=glm_ii_pcdt_p+length(data);
        
        ii_p_stats=ii_p_stats+1;
        p_corr_dt_peak_stats(ii_p_stats).data=data;
        p_corr_dt_peak_stats(ii_p_stats).description=[group_legend{grNo} ' ' ...
            'Prefrontal'];
    end
    
    
    
    
    ylim([45 100])
    
    title(['LDA percent correct for peak theta/' PACnames{PACii} ])
    
    
    xticks([1 2 3 5 6 7])
    xticklabels({'hwt', 'hH', 'hKO', 'pwt', 'pH', 'pKO'})
    
    
    ylabel('Percent correct')
    
    %Perform the glm
    fprintf(1, ['\n\nglm for mean percent correct peak Theta/' PACnames{PACii} '\n'])
    tbl = table(glm_pcorr_dt_peak.data',glm_pcorr_dt_peak.group',glm_pcorr_dt_peak.brainArea',...
        'VariableNames',{'mean_pc','group','hippo_pre'});
    mdl = fitglm(tbl,'mean_pc~group+hippo_pre+group*hippo_pre'...
        ,'CategoricalVars',[2,3])
    
    %Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for mean percent correct peak Theta/' PACnames{PACii} '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(p_corr_dt_peak_stats);
        fprintf(1, '\n\n')
    catch
    end
    
    %Bar graph plot for trough percent correct
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
    hold on
    
    
    bar_offset = 0;
    
    per_ii=1;  %We only do proficient
    
    for grNo=1:3
        
        %Hippocampus
        bar_offset = bar_offset +1;
        
        
        switch grNo
            case 1
                bar(bar_offset,mean(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).data),'g','LineWidth', 3,'EdgeColor','none')
            case 2
                bar(bar_offset,mean(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).data),'b','LineWidth', 3,'EdgeColor','none')
            case 3
                bar(bar_offset,mean(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).data),'y','LineWidth', 3,'EdgeColor','none')
        end
        
        
        CI = bootci(1000, {@mean, hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).data},'type','cper');
        plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        plot(bar_offset*ones(1,length(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).data)),hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        
        data=hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).data;
        glm_pcorr_dt_trough.data(glm_ii_pcdt_t+1:glm_ii_pcdt_t+length(data))=data;
        glm_pcorr_dt_trough.group(glm_ii_pcdt_t+1:glm_ii_pcdt_t+length(data))=grNo;
        glm_pcorr_dt_trough.brainArea(glm_ii_pcdt_t+1:glm_ii_pcdt_t+length(data))=1;
        glm_ii_pcdt_t=glm_ii_pcdt_t+length(data);
        
        ii_t_stats=ii_t_stats+1;
        p_corr_dt_trough_stats(ii_t_stats).data=data;
        p_corr_dt_trough_stats(ii_t_stats).description=[group_legend{grNo} ' ' ...
            'Hippocampus'];
        
        bar_offset = bar_offset +1;
        
        
        switch grNo
            case 1
                bar(bar_offset,mean(pre_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).data),'g','LineWidth', 3,'EdgeColor','none')
            case 2
                bar(bar_offset,mean(pre_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).data),'b','LineWidth', 3,'EdgeColor','none')
            case 3
                bar(bar_offset,mean(pre_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).data),'y','LineWidth', 3,'EdgeColor','none')
        end
        
        
        CI = bootci(1000, {@mean, pre_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).data},'type','cper');
        plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        plot(bar_offset*ones(1,length(pre_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).data)),pre_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        
        bar_offset = bar_offset + 1;
        
        data=pre_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).data;
        glm_pcorr_dt_trough.data(glm_ii_pcdt_t+1:glm_ii_pcdt_t+length(data))=data;
        glm_pcorr_dt_trough.group(glm_ii_pcdt_t+1:glm_ii_pcdt_t+length(data))=grNo;
        glm_pcorr_dt_trough.brainArea(glm_ii_pcdt_t+1:glm_ii_pcdt_t+length(data))=2;
        glm_ii_pcdt_t=glm_ii_pcdt_t+length(data);
        
        ii_t_stats=ii_t_stats+1;
        p_corr_dt_trough_stats(ii_t_stats).data=data;
        p_corr_dt_trough_stats(ii_t_stats).description=[group_legend{grNo} ' ' ...
            'Prefrontal'];
    end
    
    
    
    
    ylim([45 100])
    
    title(['LDA percent correct for trough theta/' PACnames{PACii} ])
    
    
    xticks([1 2 3 5 6 7])
    xticklabels({'hwt', 'hH', 'hKO', 'pwt', 'pH', 'pKO'})
    
    
    ylabel('Percent correct')
    
    %Perform the glm
    fprintf(1, ['\n\nglm for mean percent correct trough Theta/' PACnames{PACii} '\n'])
    tbl = table(glm_pcorr_dt_trough.data',glm_pcorr_dt_trough.group',glm_pcorr_dt_trough.brainArea',...
        'VariableNames',{'mean_pc','group','hippo_pre'});
    mdl = fitglm(tbl,'mean_pc~group+hippo_pre+group*hippo_pre'...
        ,'CategoricalVars',[2,3])
    
    %Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for mean percent correct trough Theta/' PACnames{PACii} '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(p_corr_dt_trough_stats);
        fprintf(1, '\n\n')
    catch
    end
    
    
end

pfft=1;
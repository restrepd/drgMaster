function drgSummaryDiscPcorrPRPCaMKIIin
%Analyzes the linear discriminant analysis performed by drgAnalyzeLFPDiscriminantBatchCaMKIIGrp

%
%Performs summary analises for LDA and  PAC analysis performed with case 19
%of drgAnalysisBatchLFP

warning('off')

close all
clear all

edges=[30:0.5:100];
rand_offset=0.5;

figNo=0;

PACnames{1}='Beta';
PACnames{2}='gamma';


prof_naive_leg{1}='Proficient';
prof_naive_leg{2}='Naive';


group_legend{1}='KN92';
group_legend{2}='KN93';

evTypeLabels{1}='S+';
evTypeLabels{2}='S-';

peak_label{1}='Trough';
peak_label{2}='Peak';



%Location of files for proficient=[80 100]

%Path names
hippPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Drug infusion/';
prePathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Drug infusion/';

%Text file for statistical output
fileID = fopen([hippPathName 'drgSummaryDiscPcorrPRPCaMKIIstatsKN93.txt'],'w');

%Files


%aceto
hippFileName{2}='pcorr_Discriminant_CaMKII_druginf_AcetoKN93_04182022.mat';
file_legend{2}='aceto';

%ethylben
hippFileName{1}='pcorr_Discriminant_CaMKII_druginf_ethylbenzoateKN93_04192022.mat';
file_legend{1}='ethylben';



%Load the table of mouse numbers
mouse_no_table='/Users/restrepd/Documents/Projects/CaMKII_analysis/Reply_to_reviewers/camkii_mice_per_odor_pair_for_discriminant_inf.xlsx';
T_mouse_no = readtable(mouse_no_table);


%Now do the calculations per mouse per odor pair
%Separate calculations for each odor pair


for PACii=1:2
    
    
    %Now process percent correct for the hippocampus, proficient, peak
    
    glm_disc=[];
    glm_ii_disc=0;
    
    p_corr_stats=[];
    ii_stats=0;
    
    %Bar graph plot for peak percent correct
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .2 .3 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_offset = 0;
    
    per_ii=1;  %We only do proficient
    
    all_ps=[];
    
    for grNo=1:2
        these_percorr=[];
        
        for fileNo=1:2
            if ~isempty(hippFileName{fileNo})
                load([hippPathName hippFileName{fileNo}])
                hippo_pcorr_out=pcorr_out;
                
                these_percorr=[these_percorr hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data];
            end
        end
         
        %Hippocampus
        bar_offset = bar_offset +1;
        
        switch grNo
            case 1
                bar(bar_offset,mean(these_percorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 223/255])
            case 2
                bar(bar_offset,mean(these_percorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[204/255 121/255 167/255])
            case 3
                bar(bar_offset,mean(these_percorr),'y','LineWidth', 3,'EdgeColor','none')
        end
        
        %Violin plot
        
        [mean_out, CIout]=drgViolinPoint(these_percorr,edges,bar_offset,rand_offset,'k','k',7);
        
        %         CI = bootci(1000, {@mean, hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data},'type','cper');
        %         plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        %         plot(bar_offset*ones(1,length(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data)),hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        %
        glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_percorr))=these_percorr;
        glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_percorr))=grNo*ones(1,length(these_percorr));
        glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_percorr))=ones(1,length(these_percorr));
        glm_ii_disc=glm_ii_disc+length(these_percorr);
        
        
        ii_stats=ii_stats+1;
        p_corr_stats(ii_stats).data=these_percorr;
        p_corr_stats(ii_stats).description=[group_legend{grNo} ' peak'];
        
    end
    
    yticks([50 60 70 80 90 100])
    ylim([50 100])
    
    xlim([0.3 2.8])
    
    xticks([1 2])
    xticklabels({'KN92', 'KN93'})
    
    title(['LDA percent correct for hippocampus for peak per mouse per odor pair theta/' PACnames{PACii} ])
    
    
    ylabel('Percent correct')

    
    %Bar graph plot for trough percent correct
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .2 .3 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_offset = 0;
    
    per_ii=1;  %We only do proficient
    
    all_ps=[];
    
    for grNo=1:2
        these_percorr=[];
        
        for fileNo=1:2
            if ~isempty(hippFileName{fileNo})
                load([hippPathName hippFileName{fileNo}])
                hippo_pcorr_out=pcorr_out;
                
                these_percorr=[these_percorr hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data];
            end
        end
        
        %Hippocampus
        bar_offset = bar_offset +1;
        
        switch grNo
            case 1
                bar(bar_offset,mean(these_percorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 223/255])
            case 2
                bar(bar_offset,mean(these_percorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[204/255 121/255 167/255])
            case 3
                bar(bar_offset,mean(these_percorr),'y','LineWidth', 3,'EdgeColor','none')
        end
        
        %Violin plot
        
        [mean_out, CIout]=drgViolinPoint(these_percorr,edges,bar_offset,rand_offset,'k','k',7);
        
        %         CI = bootci(1000, {@mean, hippo_pcorr_out.PACii(PACii).troughs.pcorr(per_ii).group(grNo).odor_data},'type','cper');
        %         plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        %         plot(bar_offset*ones(1,length(hippo_pcorr_out.PACii(PACii).troughs.pcorr(per_ii).group(grNo).odor_data)),hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        %
        glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_percorr))=these_percorr;
        glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_percorr))=grNo*ones(1,length(these_percorr));
        glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_percorr))=zeros(1,length(these_percorr));
        glm_ii_disc=glm_ii_disc+length(these_percorr);
        
        
        ii_stats=ii_stats+1;
        p_corr_stats(ii_stats).data=these_percorr;
        p_corr_stats(ii_stats).description=[group_legend{grNo} ' trough'];
        
    end
    
    yticks([50 60 70 80 90 100])
    ylim([50 100])
    
    xlim([0.3 2.8])
    
    xticks([1 2])
    xticklabels({'KN92', 'KN93'})
    
    title(['LDA percent correct for hippocampus for trough per mouse per odor pair theta/' PACnames{PACii} ])
    
    
    ylabel('Percent correct')
    
    %Perform the glm
    fprintf(1, ['glm for LDA accuracy per mouse per odor pair for theta/' PACnames{PACii} '\n'])
    fprintf(fileID, ['glm for LDA accuracy per mouse per odor pair for theta/' PACnames{PACii} '\n']);
    tbl = table(glm_disc.data',glm_disc.group',glm_disc.peak',...
        'VariableNames',{'performance','KN92_vs_KN93','peak_vs_trough'});
    mdl = fitglm(tbl,'performance~KN92_vs_KN93+peak_vs_trough'...
        ,'CategoricalVars',[2 3])
    
      txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
    fprintf(fileID,'%s\n', txt);
    
    %Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for LDA accuracy per mouse per electrode for theta' PACnames{PACii} '\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for LDA accuracy per mouse per electrode for theta' PACnames{PACii} '\n']);
   
    [output_data] = drgMutiRanksumorTtest(p_corr_stats, fileID);
    
    
%     %Nested ANOVAN
%     %https://www.mathworks.com/matlabcentral/answers/491365-within-between-subjects-in-anovan
%     nesting=[0 0 0; ... % This line indicates that group factor is not nested in any other factor.
%         0 0 0; ... % This line indicates that perCorr is not nested in any other factor.
% %         0 0 0 0; ... % This line indicates that event is not nested in any other factor.
%         1 1 0];    % This line indicates that mouse_no (the last factor) is nested under group, perCorr and event
%     % (the 1 in position 1 on the line indicates nesting under the first factor).
%     figNo=figNo+1;
%     
%     [p anovanTbl stats]=anovan(glm_disc.data,{glm_disc.group glm_disc.peak glm_disc.mouse_no},...
%         'model','interaction',...
%         'nested',nesting,...
%         'varnames',{'KN92_vs_KN93', 'peak_vs_tough','mouse_no'});
%     
%     fprintf(fileID, ['\n\nNested ANOVAN for LDA accuracy per mouse per odor pair for '  bandwidth_names{pacii} '\n'])
%     drgWriteANOVANtbl(anovanTbl,fileID);
%     fprintf(fileID, '\n\n');
%     
end


fclose(fileID);
pfft=1;
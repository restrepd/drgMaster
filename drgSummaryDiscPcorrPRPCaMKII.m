function drgSummaryDiscPcorrPRPCaMKII
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

fwd_rev_leg{1}='Forward';
fwd_rev_leg{2}='Reverse';

group_legend{1}='WT';
group_legend{2}='Het';
group_legend{3}='KO';

evTypeLabels{1}='S+';
evTypeLabels{2}='S-';

peak_label{1}='Trough';
peak_label{2}='Peak';


% %Location of files for proficient = [85 100]
% 
% %Path names
% hippPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Discriminant/';
% prePathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Discriminant/';
% 
% %Files
% 
% %Hippocampus
% 
% %pzz1PAEA
% hippFileName{8}='pcorr_Discriminant_CaMKIIpzz1paea_disc_01302021_hippo2.mat';
% file_legend{8}='pzz1PAEA';
% 
% %pzz1EAPA
% hippFileName{7}='pcorr_Discriminant_CaMKIIpzz1ethylace_disc_PRPhipp2_01222020.mat';
% file_legend{7}='pzz1EAPA';
% 
% %pz1PAEA
% hippFileName{6}='pcorr_Discriminant_CaMKIIpz1paea_disc_02042021_hipp2.mat';
% file_legend{6}='pz1PAEA';
% 
% %pz1EAPA
% hippFileName{5}='pcorr_Discriminant_CaMKIIPZ1EAPA_disc_PRPhippo2_0212021.mat';
% file_legend{5}='pz1EAPA';
% 
% %PAEA
% hippFileName{4}='pcorr_Discriminant_CaMKIIPAEA_disc_PRPhipp2_01222020.mat';
% file_legend{4}='PAEA';
% 
% %EBAP
% hippFileName{2}='pcorr_Discriminant_CaMKIIEBAP_disc_PRPhipp2_01222020.mat';
% file_legend{2}='EBAP';
% 
% %EAPA
% hippFileName{3}='pcorr_Discriminant_CaMKIIEAPA_disc_PRPhippo2_02082021.mat';
% file_legend{3}='EAPA';
% 
% %APEB aceto
% hippFileName{1}='pcorr_Discriminant_CaMKIIaceto_disc_PRPhippo2_0272021.mat';
% file_legend{1}='APEB';
% 
% %Prefrontal files
% 
% %pzz1PAEA
% preFileName{8}='pcorr_Discriminant_CaMKIIpzz1paea_disc_02012021_pref2.mat';
% 
% %pzz1EAPA
% preFileName{7}='pcorr_Discriminant_CaMKIIpzz1ethylace_disc_PRPpref2_01262020.mat';
% 
% %pz1PAEA
% preFileName{6}='pcorr_Discriminant_CaMKIIpz1paea_disc_02072021_pref2.mat';
% 
% %pz1EAPA
% preFileName{5}='pcorr_Discriminant_CaMKIIPZ1EAPA_disc_PRPprefront2_01302021.mat';
% 
% %PAEA
% preFileName{4}='pcorr_Discriminant_CaMKIIPAEA_disc_PRPprefront2_01252021.mat';
% 
% %EBAP
% preFileName{2}='pcorr_Discriminant_CaMKIIEBAP_disc_PRPprefront2_01262021.mat';
% 
% %EAPA
% preFileName{3}='pcorr_Discriminant_CaMKIIEAPA_disc_PRPrefront2_01232021.mat';
% 
% %APEB aceto
% preFileName{1}='pcorr_Discriminant_CaMKIIaceto_disc_PRPprefront2_01302021.mat';


%Location of files for proficient=[80 100]

%Path names
hippPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Discriminant new 80/';
prePathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Discriminant new 80/';

%Files

%Hippocampus

%pzz1PAEA
hippFileName{8}='pcorr_Discriminant_CaMKIIpzz1paea_disc80_04052021_hipp.mat';
file_legend{8}='pzz1PAEA';

%pzz1EAPA
hippFileName{7}='pcorr_Discriminant_CaMKIIpzz1ethylace_disc_80hipp2_04082021.mat';
file_legend{7}='pzz1EAPA';

%pz1PAEA
hippFileName{6}='pcorr_Discriminant_CaMKIIpz1paea_disc_04112021_80hipp.mat';
file_legend{6}='pz1PAEA';

%pz1EAPA
hippFileName{5}='pcorr_Discriminant_CaMKIIpz1eapa_disc_80_04142021_hipp.mat';
file_legend{5}='pz1EAPA';

%PAEA
hippFileName{4}='pcorr_Discriminant_CaMKIIPAEA_disc80_hipp2_04162021.mat';
file_legend{4}='PAEA';

%EBAP
hippFileName{2}='pcorr_Discriminant_CaMKIIEBAP_disc_80_hipp04182021.mat';
file_legend{2}='EBAP';

%EAPA
hippFileName{3}='pcorr_Discriminant_CaMKIIEAPA_dis80c_hipp2_04212021.mat';
file_legend{3}='EAPA';

%APEB aceto
hippFileName{1}='pcorr_Discriminant_CaMKIIAPEB_disc_PRP80_hippo_04142021.mat';
file_legend{1}='APEB';

%Prefrontal files

%pzz1PAEA
preFileName{8}='pcorr_Discriminant_CaMKIIpzz1paea_disc80_04052021_pref.mat';

%pzz1EAPA
preFileName{7}='pcorr_Discriminant_CaMKIIpzz1ethylace_disc_80pref2_04072021.mat';

%pz1PAEA
preFileName{6}='pcorr_Discriminant_CaMKIIpz1paea_disc_04112021_80pref.mat';

%pz1EAPA
preFileName{5}='pcorr_Discriminant_CaMKIIEAPA_dis80c_pre2_04212021.mat';

%PAEA
preFileName{4}='pcorr_Discriminant_CaMKIIPAEA_disc80_pre2_04162021.mat';

%EBAP
preFileName{2}='pcorr_Discriminant_CaMKIIEBAP_disc_80_pre04182021.mat';

%EAPA
preFileName{3}='pcorr_Discriminant_CaMKIIEAPA_dis80c_pre2_04212021.mat';

%APEB aceto
preFileName{1}='pcorr_Discriminant_CaMKIIAPEB_disc_PRP80_pre_04232021.mat';


%Separate calculations for each odor pair
%Now process percent correct for the hippocampus, proficient, peak

for PACii=1:2
    
    %     glm_pcorr_dt_peak=[];
    %     glm_ii_pcdt_p=0;
    %     glm_pcorr_dt_trough=[];
    %     glm_ii_pcdt_t=0;
    
    
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
    
    all_ps=[];
    
    for fileNo=1:8
        
        if ~isempty(hippFileName{fileNo})
            load([hippPathName hippFileName{fileNo}])
            hippo_pcorr_out=pcorr_out;
            
            p_corr_dt_peak_stats=[];
            ii_p_stats=0;
            
            for grNo=1:3
                
                %Hippocampus
                bar_offset = bar_offset +1;
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data),'y','LineWidth', 3,'EdgeColor','none')
                end
                
                
                CI = bootci(1000, {@mean, hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data)),hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                ii_p_stats=ii_p_stats+1;
                p_corr_dt_peak_stats(ii_p_stats).data=hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data;
                p_corr_dt_peak_stats(ii_p_stats).description=[file_legend{fileNo} ' ' group_legend{grNo} ' ' ...
                    'Hippocampus'];
            end
            
            bar_offset = bar_offset +1;
            
            %Do ranksum/t test
            [output_data] = drgMutiRanksumorTtestNoFDR(p_corr_dt_peak_stats);
            all_ps=[all_ps output_data.p];
            
        else
            bar_offset = bar_offset +7;
        end
        
    end
    pFDR=drsFDRpval(all_ps);
    fprintf(1, ['\n\npFDR = %d \n\n'],pFDR)
    
    ylim([50 100])
    
    title(['LDA percent correct for hippocampus for peak theta/' PACnames{PACii} ])
    
    
    
    
    ylabel('Percent correct')
end




%Now process percent correct for the hippocampus, proficient, trough

for PACii=1:2
    
    %     glm_pcorr_dt_trough=[];
    %     glm_ii_pcdt_p=0;
    %     glm_pcorr_dt_trough=[];
    %     glm_ii_pcdt_t=0;
    
    
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
    
    all_ps=[];
    
    for fileNo=1:8
        
        if ~isempty(hippFileName{fileNo})
            load([hippPathName hippFileName{fileNo}])
            hippo_pcorr_out=pcorr_out;
            
            p_corr_dt_trough_stats=[];
            ii_p_stats=0;
            
            for grNo=1:3
                
                %Hippocampus
                bar_offset = bar_offset +1;
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data),'y','LineWidth', 3,'EdgeColor','none')
                end
                
                
                CI = bootci(1000, {@mean, hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data)),hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                ii_p_stats=ii_p_stats+1;
                p_corr_dt_trough_stats(ii_p_stats).data=hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data;
                p_corr_dt_trough_stats(ii_p_stats).description=[file_legend{fileNo} ' ' group_legend{grNo} ' ' ...
                    'Hippocampus'];
            end
            
            bar_offset = bar_offset +1;
            
            %Do ranksum/t test
            [output_data] = drgMutiRanksumorTtestNoFDR(p_corr_dt_trough_stats);
            all_ps=[all_ps output_data.p];
            
        else
            bar_offset = bar_offset +7;
        end
        
    end
    pFDR=drsFDRpval(all_ps);
    fprintf(1, ['\n\npFDR = %d \n\n'],pFDR)
    
    ylim([50 100])
    
    title(['LDA percent correct for hippocampus for trough theta/' PACnames{PACii} ])
    
    
    
    
    ylabel('Percent correct')
end


%Now process percent correct for the prefrontal, proficient, peak

for PACii=1:2
    
    %     glm_pcorr_dt_peak=[];
    %     glm_ii_pcdt_p=0;
    %     glm_pcorr_dt_trough=[];
    %     glm_ii_pcdt_t=0;
    
    
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
    
    all_ps=[];
    
    for fileNo=1:8
        
        if ~isempty(preFileName{fileNo})
            load([prePathName preFileName{fileNo}])
            pref_pcorr_out=pcorr_out;
            
            p_corr_dt_peak_stats=[];
            ii_p_stats=0;
            
            for grNo=1:3
                
                %Hippocampus
                bar_offset = bar_offset +1;
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data),'y','LineWidth', 3,'EdgeColor','none')
                end
                
                
                CI = bootci(1000, {@mean, pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data)),pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                ii_p_stats=ii_p_stats+1;
                p_corr_dt_peak_stats(ii_p_stats).data=pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data;
                p_corr_dt_peak_stats(ii_p_stats).description=[file_legend{fileNo} ' ' group_legend{grNo} ' ' ...
                    'Hippocampus'];
            end
            
            bar_offset = bar_offset +1;
            
            %Do ranksum/t test
            [output_data] = drgMutiRanksumorTtestNoFDR(p_corr_dt_peak_stats);
            all_ps=[all_ps output_data.p];
            
        else
            bar_offset = bar_offset +7;
        end
        
    end
    pFDR=drsFDRpval(all_ps);
    fprintf(1, ['\n\npFDR = %d \n\n'],pFDR)
    
    ylim([50 100])
    
    xticks([2 6 10 14 18 22 26 30])
    xticklabels({file_legend{1}, file_legend{2},file_legend{3}, file_legend{4},file_legend{5}, file_legend{6},file_legend{7}, file_legend{8}})
    
    title(['LDA percent correct for prefrontal for peak theta/' PACnames{PACii} ])
    
    
    
    
    ylabel('Percent correct')
end




%Now process percent correct for the prefrontal, proficient, trough

for PACii=1:2
    
    %     glm_pcorr_dt_trough=[];
    %     glm_ii_pcdt_p=0;
    %     glm_pcorr_dt_trough=[];
    %     glm_ii_pcdt_t=0;
    
    
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
    
    all_ps=[];
    
    for fileNo=1:8
        
        if ~isempty(preFileName{fileNo})
            load([prePathName preFileName{fileNo}])
            pref_pcorr_out=pcorr_out;
            
            p_corr_dt_trough_stats=[];
            ii_p_stats=0;
            
            for grNo=1:3
                
                %Hippocampus
                bar_offset = bar_offset +1;
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data),'y','LineWidth', 3,'EdgeColor','none')
                end
                
                
                CI = bootci(1000, {@mean, pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data)),pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                ii_p_stats=ii_p_stats+1;
                p_corr_dt_trough_stats(ii_p_stats).data=pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data;
                p_corr_dt_trough_stats(ii_p_stats).description=[file_legend{fileNo} ' ' group_legend{grNo} ' ' ...
                    'Hippocampus'];
            end
            
            bar_offset = bar_offset +1;
            
            %Do ranksum/t test
            [output_data] = drgMutiRanksumorTtestNoFDR(p_corr_dt_trough_stats);
            all_ps=[all_ps output_data.p];
            
        else
            bar_offset = bar_offset +7;
        end
        
    end
    pFDR=drsFDRpval(all_ps);
    fprintf(1, ['\n\npFDR = %d \n\n'],pFDR)
    
    ylim([50 100])
    
    xticks([2 6 10 14 18 22 26 30])
    xticklabels({file_legend{1}, file_legend{2},file_legend{3}, file_legend{4},file_legend{5}, file_legend{6},file_legend{7}, file_legend{8}})
    
    title(['LDA percent correct for prefrontal for trough theta/' PACnames{PACii} ])
    
    
    
    
    ylabel('Percent correct')
end

%Do a GLM to compare groups/brain regions for each odor pair trough
for PACii=1:2
    
    %     glm_pcorr_dt_trough=[];
    %     glm_ii_pcdt_p=0;
    %     glm_pcorr_dt_trough=[];
    %     glm_ii_pcdt_t=0;
    
    
    
    per_ii=1;  %We only do proficient
    
    
    
    for fileNo=1:8
        
        glm_disc=[];
        glm_ii_disc=0;
        
        if ~isempty(hippFileName{fileNo})
            
            load([prePathName preFileName{fileNo}])
            pref_pcorr_out=pcorr_out;
            
            load([hippPathName hippFileName{fileNo}])
            hippo_pcorr_out=pcorr_out;
            
            
            for grNo=1:3
                
                %Hippocampus
                glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data))=hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data;
                glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data))=grNo*ones(1,length(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data));
                glm_disc.brain_area(glm_ii_disc+1:glm_ii_disc+length(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data))=zeros(1,length(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data));
                glm_ii_disc=glm_ii_disc+length(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data);
                
                %Prefrontal
                glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data))=pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data;
                glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data))=grNo*ones(1,length(pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data));
                glm_disc.brain_area(glm_ii_disc+1:glm_ii_disc+length(pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data))=ones(1,length(pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data));
                glm_ii_disc=glm_ii_disc+length(pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data);
                
                
            end
            
            
            %Perform the glm
            fprintf(1, ['glm for performance for trough theta/' PACnames{PACii} '\n'])
            tbl = table(glm_disc.data',glm_disc.group',glm_disc.brain_area',...
                'VariableNames',{'performance','group','brain_area'});
            mdl = fitglm(tbl,'performance~group+brain_area+group*brain_area'...
                ,'CategoricalVars',[2,3])
            
        else
            bar_offset = bar_offset +7;
        end
        
    end
    
    
end



%Do a GLM to compare groups/brain regions for each odor pair peak
for PACii=1:2
    
    %     glm_pcorr_dt_peak=[];
    %     glm_ii_pcdt_p=0;
    %     glm_pcorr_dt_peak=[];
    %     glm_ii_pcdt_t=0;
    
    
    
    per_ii=1;  %We only do proficient
    
    
    
    for fileNo=1:8
        
        glm_disc=[];
        glm_ii_disc=0;
        
        if ~isempty(hippFileName{fileNo})
            
            load([prePathName preFileName{fileNo}])
            pref_pcorr_out=pcorr_out;
            
            load([hippPathName hippFileName{fileNo}])
            hippo_pcorr_out=pcorr_out;
            
            
            for grNo=1:3
                
                %Hippocampus
                glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data))=hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data;
                glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data))=grNo*ones(1,length(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data));
                glm_disc.brain_area(glm_ii_disc+1:glm_ii_disc+length(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data))=zeros(1,length(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data));
                glm_ii_disc=glm_ii_disc+length(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data);
                
                %Prefrontal
                glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data))=pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data;
                glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data))=grNo*ones(1,length(pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data));
                glm_disc.brain_area(glm_ii_disc+1:glm_ii_disc+length(pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data))=ones(1,length(pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data));
                glm_ii_disc=glm_ii_disc+length(pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data);
                
                
            end
            
            
            %Perform the glm
            fprintf(1, ['glm for performance for peak theta/' PACnames{PACii} '\n'])
            tbl = table(glm_disc.data',glm_disc.group',glm_disc.brain_area',...
                'VariableNames',{'performance','group','brain_area'});
            mdl = fitglm(tbl,'performance~group+brain_area+group*brain_area'...
                ,'CategoricalVars',[2,3])
            
        else
            bar_offset = bar_offset +7;
        end
        
    end
    
    
end


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
    
    for grNo=1:3
        these_percorr=[];
        
        for fileNo=1:8
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
                bar(bar_offset,mean(these_percorr),'g','LineWidth', 3,'EdgeColor','none')
            case 2
                bar(bar_offset,mean(these_percorr),'b','LineWidth', 3,'EdgeColor','none')
            case 3
                bar(bar_offset,mean(these_percorr),'y','LineWidth', 3,'EdgeColor','none')
        end
        
        %Violin plot
        
        [mean_out, CIout]=drgViolinPoint(these_percorr,edges,bar_offset,rand_offset,'k','k',3);
        
        %         CI = bootci(1000, {@mean, hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data},'type','cper');
        %         plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        %         plot(bar_offset*ones(1,length(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data)),hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        %
        glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_percorr))=these_percorr;
        glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_percorr))=grNo*ones(1,length(these_percorr));
        glm_ii_disc=glm_ii_disc+length(these_percorr);
        
        
        ii_stats=ii_stats+1;
        p_corr_stats(ii_stats).data=these_percorr;
        
        switch grNo
            case 1
                p_corr_stats(ii_stats).description=['WT'];
            case 2
                p_corr_stats(ii_stats).description=['Het'];
            case 3
                p_corr_stats(ii_stats).description=['KO'];
        end
        
    end
    
    yticks([50 60 70 80 90 100])
    ylim([50 100])
    
    xlim([0.3 3.8])
    
    xticks([1 2 3])
    xticklabels({'WT', 'Het', 'KO'})
    
    title(['LDA percent correct for hippocampus for peak per mouse per odor pair theta/' PACnames{PACii} ])
    

    ylabel('Percent correct')
    
       %Perform the glm
    fprintf(1, ['glm for percent correct vs genotype for hippocampus for peak per mouse per odor pair for theta/' PACnames{PACii} '\n'])
    tbl = table(glm_disc.data',glm_disc.group',...
        'VariableNames',{'performance','group'});
    mdl = fitglm(tbl,'performance~group'...
        ,'CategoricalVars',[2])
    
    %Do ranksum/t test
    [output_data] = drgMutiRanksumorTtest(p_corr_stats);
    
    
    
    %Now process percent correct for the hippocampus, proficient, trough
    
    glm_disc=[];
    glm_ii_disc=0;
    
    p_corr_stats=[];
    ii_stats=0;
    
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
    
    for grNo=1:3
        these_percorr=[];
        
        for fileNo=1:8
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
                bar(bar_offset,mean(these_percorr),'g','LineWidth', 3,'EdgeColor','none')
            case 2
                bar(bar_offset,mean(these_percorr),'b','LineWidth', 3,'EdgeColor','none')
            case 3
                bar(bar_offset,mean(these_percorr),'y','LineWidth', 3,'EdgeColor','none')
        end
        
        %Violin plot
        
        [mean_out, CIout]=drgViolinPoint(these_percorr,edges,bar_offset,rand_offset,'k','k',3);
        
        %         CI = bootci(1000, {@mean, hippo_pcorr_out.PACii(PACii).troughs.pcorr(per_ii).group(grNo).odor_data},'type','cper');
        %         plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        %         plot(bar_offset*ones(1,length(hippo_pcorr_out.PACii(PACii).troughs.pcorr(per_ii).group(grNo).odor_data)),hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        %
        glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_percorr))=these_percorr;
        glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_percorr))=grNo*ones(1,length(these_percorr));
        glm_ii_disc=glm_ii_disc+length(these_percorr);
        
        
        ii_stats=ii_stats+1;
        p_corr_stats(ii_stats).data=these_percorr;
        
        switch grNo
            case 1
                p_corr_stats(ii_stats).description=['WT'];
            case 2
                p_corr_stats(ii_stats).description=['Het'];
            case 3
                p_corr_stats(ii_stats).description=['KO'];
        end
        
    end
    
    yticks([50 60 70 80 90 100])
    ylim([50 100])
    
    xlim([0.3 3.8])
    
    xticks([1 2 3])
    xticklabels({'WT', 'Het', 'KO'})
    
    title(['LDA percent correct for hippocampus for trough per mouse per odor pair theta/' PACnames{PACii} ])
    

    ylabel('Percent correct')
    
       %Perform the glm
    fprintf(1, ['glm for percent correct vs genotype for hippocampus for trough per mouse per odor pair for theta/' PACnames{PACii} '\n'])
    tbl = table(glm_disc.data',glm_disc.group',...
        'VariableNames',{'performance','group'});
    mdl = fitglm(tbl,'performance~group'...
        ,'CategoricalVars',[2])
    
    %Do ranksum/t test
    [output_data] = drgMutiRanksumorTtest(p_corr_stats);
    
    
     %Now process percent correct for the prefrontal, proficient, peak
    
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
    
    for grNo=1:3
        these_percorr=[];
        
        for fileNo=1:8
            if ~isempty(preFileName{fileNo})
                load([prePathName preFileName{fileNo}])
                pre_pcorr_out=pcorr_out;
                
                these_percorr=[these_percorr pre_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data];
            end
        end
        
        %Hippocampus
        bar_offset = bar_offset +1;
        
        switch grNo
            case 1
                bar(bar_offset,mean(these_percorr),'g','LineWidth', 3,'EdgeColor','none')
            case 2
                bar(bar_offset,mean(these_percorr),'b','LineWidth', 3,'EdgeColor','none')
            case 3
                bar(bar_offset,mean(these_percorr),'y','LineWidth', 3,'EdgeColor','none')
        end
        
        %Violin plot
        
        [mean_out, CIout]=drgViolinPoint(these_percorr,edges,bar_offset,rand_offset,'k','k',3);
        
        %         CI = bootci(1000, {@mean, hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data},'type','cper');
        %         plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        %         plot(bar_offset*ones(1,length(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data)),hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        %
        glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_percorr))=these_percorr;
        glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_percorr))=grNo*ones(1,length(these_percorr));
        glm_ii_disc=glm_ii_disc+length(these_percorr);
        
        
        ii_stats=ii_stats+1;
        p_corr_stats(ii_stats).data=these_percorr;
        
        switch grNo
            case 1
                p_corr_stats(ii_stats).description=['WT'];
            case 2
                p_corr_stats(ii_stats).description=['Het'];
            case 3
                p_corr_stats(ii_stats).description=['KO'];
        end
        
    end
    
    yticks([50 60 70 80 90 100])
    ylim([50 100])
    
    xlim([0.3 3.8])
    
    xticks([1 2 3])
    xticklabels({'WT', 'Het', 'KO'})
    
    title(['LDA percent correct for prefrontal for peak per mouse per odor pair theta/' PACnames{PACii} ])
    

    ylabel('Percent correct')
    
       %Perform the glm
    fprintf(1, ['glm for percent correct vs genotype for prefrontal for peak per mouse per odor pair for theta/' PACnames{PACii} '\n'])
    tbl = table(glm_disc.data',glm_disc.group',...
        'VariableNames',{'performance','group'});
    mdl = fitglm(tbl,'performance~group'...
        ,'CategoricalVars',[2])
    
    %Do ranksum/t test
    [output_data] = drgMutiRanksumorTtest(p_corr_stats);
    
    
    
    %Now process percent correct for the prefrontal, proficient, trough
    
    glm_disc=[];
    glm_ii_disc=0;
    
    p_corr_stats=[];
    ii_stats=0;
    
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
    
    for grNo=1:3
        these_percorr=[];
        
        for fileNo=1:8
            if ~isempty(preFileName{fileNo})
                load([prePathName preFileName{fileNo}])
                pre_pcorr_out=pcorr_out;
                
                these_percorr=[these_percorr pre_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data];
            end
        end
        
        %Hippocampus
        bar_offset = bar_offset +1;
        
        switch grNo
            case 1
                bar(bar_offset,mean(these_percorr),'g','LineWidth', 3,'EdgeColor','none')
            case 2
                bar(bar_offset,mean(these_percorr),'b','LineWidth', 3,'EdgeColor','none')
            case 3
                bar(bar_offset,mean(these_percorr),'y','LineWidth', 3,'EdgeColor','none')
        end
        
        %Violin plot
        
        [mean_out, CIout]=drgViolinPoint(these_percorr,edges,bar_offset,rand_offset,'k','k',3);
        
        %         CI = bootci(1000, {@mean, hippo_pcorr_out.PACii(PACii).troughs.pcorr(per_ii).group(grNo).odor_data},'type','cper');
        %         plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        %         plot(bar_offset*ones(1,length(hippo_pcorr_out.PACii(PACii).troughs.pcorr(per_ii).group(grNo).odor_data)),hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        %
        glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_percorr))=these_percorr;
        glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_percorr))=grNo*ones(1,length(these_percorr));
        glm_ii_disc=glm_ii_disc+length(these_percorr);
        
        
        ii_stats=ii_stats+1;
        p_corr_stats(ii_stats).data=these_percorr;
        
        switch grNo
            case 1
                p_corr_stats(ii_stats).description=['WT'];
            case 2
                p_corr_stats(ii_stats).description=['Het'];
            case 3
                p_corr_stats(ii_stats).description=['KO'];
        end
        
    end
    
    yticks([50 60 70 80 90 100])
    ylim([50 100])
    
    xlim([0.3 3.8])
    
    xticks([1 2 3])
    xticklabels({'WT', 'Het', 'KO'})
    
    title(['LDA percent correct for prefrontal for trough per mouse per odor pair theta/' PACnames{PACii} ])
    

    ylabel('Percent correct')
    
       %Perform the glm
    fprintf(1, ['glm for percent correct vs genotype for prefrontal for trough per mouse per odor pair for theta/' PACnames{PACii} '\n'])
    tbl = table(glm_disc.data',glm_disc.group',...
        'VariableNames',{'performance','group'});
    mdl = fitglm(tbl,'performance~group'...
        ,'CategoricalVars',[2])
    
    %Do ranksum/t test
    [output_data] = drgMutiRanksumorTtest(p_corr_stats);
    pffft=1;
    
end


%Now do the plots for the WT only


%Now process percent correct for the prefrontal and hippocampus, proficient, peak

for PACii=1:2
    
    %     glm_pcorr_dt_peak=[];
    %     glm_ii_pcdt_p=0;
    %     glm_pcorr_dt_trough=[];
    %     glm_ii_pcdt_t=0;
    
    
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
    
    all_ps=[];
    
    p_corr_dt_peak_stats=[];
    ii_p_stats=0;
    
    grNo=1;
    
    for fileNo=1:8
        
        
        %Prefrontal
        if ~isempty(preFileName{fileNo})
            load([prePathName preFileName{fileNo}])
            pref_pcorr_out=pcorr_out;
            
            bar_offset = bar_offset +1;
            bar(bar_offset,mean(pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
            
            CI = bootci(1000, {@mean, pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data)),pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            ii_p_stats=ii_p_stats+1;
            p_corr_dt_peak_stats(ii_p_stats).data=pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data;
            p_corr_dt_peak_stats(ii_p_stats).description=[file_legend{fileNo} ' ' ...
                'Prefrontal'];
            
            bar_offset = bar_offset +1;
        else
            bar_offset = bar_offset +2;
        end
        
        %Hippocampus
        if ~isempty(hippFileName{fileNo})
            load([hippPathName hippFileName{fileNo}])
            hippo_pcorr_out=pcorr_out;
            
            bar_offset = bar_offset +1;
            bar(bar_offset,mean(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
            
            CI = bootci(1000, {@mean, hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data)),hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            ii_p_stats=ii_p_stats+1;
            p_corr_dt_peak_stats(ii_p_stats).data=hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data;
            p_corr_dt_peak_stats(ii_p_stats).description=[file_legend{fileNo} ' ' ...
                'Prefrontal'];
            
            bar_offset = bar_offset +2;
        else
            bar_offset = bar_offset +3;
        end
        
        
    end
    
    %Do ranksum/t test
    [output_data] = drgMutiRanksumorTtest(p_corr_dt_peak_stats);
    all_ps=[all_ps output_data.p];
    
    
    ylim([50 100])
    
    xticks([2 6 10 14 18 22 26 30])
    xticklabels({file_legend{1}, file_legend{2},file_legend{3}, file_legend{4},file_legend{5}, file_legend{6},file_legend{7}, file_legend{8}})
    
    text(0,98,'Hippocampus','Color',[0 114/255 178/255])
    text(0,96,'Prefrontal','Color',[158/255 31/255 99/255])
    
    title(['LDA percent correct for peak theta/' PACnames{PACii} ])
    
    ylabel('Percent correct')
end




%Now process percent correct for the prefrontal, proficient, trough

for PACii=1:2
    
    %     glm_pcorr_dt_peak=[];
    %     glm_ii_pcdt_p=0;
    %     glm_pcorr_dt_trough=[];
    %     glm_ii_pcdt_t=0;
    
    
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
    
    all_ps=[];
    
    p_corr_dt_trough_stats=[];
    ii_p_stats=0;
    
    grNo=1;
    
    for fileNo=1:8
        
        
        %Prefrontal
        if ~isempty(preFileName{fileNo})
            load([prePathName preFileName{fileNo}])
            pref_pcorr_out=pcorr_out;
            
            bar_offset = bar_offset +1;
            bar(bar_offset,mean(pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
            
            CI = bootci(1000, {@mean, pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data)),pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            ii_p_stats=ii_p_stats+1;
            p_corr_dt_trough_stats(ii_p_stats).data=pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data;
            p_corr_dt_trough_stats(ii_p_stats).description=[file_legend{fileNo} ' ' ...
                'Prefrontal'];
            
            bar_offset = bar_offset +1;
        else
            bar_offset = bar_offset +2;
        end
        
        %Hippocampus
        if ~isempty(hippFileName{fileNo})
            load([hippPathName hippFileName{fileNo}])
            hippo_pcorr_out=pcorr_out;
            
            bar_offset = bar_offset +1;
            bar(bar_offset,mean(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
            
            CI = bootci(1000, {@mean, hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data)),hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            ii_p_stats=ii_p_stats+1;
            p_corr_dt_trough_stats(ii_p_stats).data=hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data;
            p_corr_dt_trough_stats(ii_p_stats).description=[file_legend{fileNo} ' ' ...
                'Prefrontal'];
            
            bar_offset = bar_offset +2;
        else
            bar_offset = bar_offset +3;
        end
        
        
    end
    
    %Do ranksum/t test
    [output_data] = drgMutiRanksumorTtest(p_corr_dt_trough_stats);
    all_ps=[all_ps output_data.p];
    
    
    ylim([50 100])
    
    xticks([2 6 10 14 18 22 26 30])
    xticklabels({file_legend{1}, file_legend{2},file_legend{3}, file_legend{4},file_legend{5}, file_legend{6},file_legend{7}, file_legend{8}})
    
    text(0,98,'Hippocampus','Color',[0 114/255 178/255])
    text(0,96,'Prefrontal','Color',[158/255 31/255 99/255])
    
    title(['LDA percent correct for trough theta/' PACnames{PACii} ])
    
    ylabel('Percent correct')
end

%Now do the wild type analysis per odor pair

%Now process percent correct for the hippocampus, peak

for PACii=1:2
    
    %     glm_pcorr_dt_peak=[];
    %     glm_ii_pcdt_p=0;
    %     glm_pcorr_dt_trough=[];
    %     glm_ii_pcdt_t=0;
    
    
    %Bar graph plot for peak percent correct
    
    %Hippocampus
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
    hold on
    
    
    mean_hippo_peak_naive=[];
    mean_hippo_peak_proficient=[];
    all_hippo_peak_shuffled=[];
    all_hippo_peak_naive=[];
    all_hippo_peak_proficient=[];
    bar_offset = 0;
    
    per_ii=1;  %We only do proficient
    
    all_ps=[];
    
    glm_disc=[];
    glm_ii_disc=0;
    
    p_corr_stats=[];
    ii_stats=0;
    
    
    for fileNo=1:8
        
        if ~isempty(hippFileName{fileNo})
            load([hippPathName hippFileName{fileNo}])
            hippo_pcorr_out=pcorr_out;
            
            
            grNo=1;
            
             
            
            %Shuffled
            bar_offset = bar_offset +1;
            per_ii=1;
            these_pcorr=hippo_pcorr_out.PACii(PACii).shuffled.pcorr(per_ii).group(grNo).odor_data;
            per_ii=2;
            these_pcorr=hippo_pcorr_out.PACii(PACii).shuffled.pcorr(per_ii).group(grNo).odor_data;
            all_hippo_peak_shuffled=[all_hippo_peak_shuffled these_pcorr];
            
            bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[120/255 120/255 120/255])
            
            
            CI = bootci(1000, {@mean, these_pcorr},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
            glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=1*ones(1,length(these_pcorr));
            glm_disc.odor_pair(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=fileNo*ones(1,length(these_pcorr));
            glm_ii_disc=glm_ii_disc+length(these_pcorr);
            
            
            ii_stats=ii_stats+1;
            p_corr_stats(ii_stats).data=these_pcorr;
            p_corr_stats(ii_stats).description=[file_legend{fileNo} ' ' ...
                'Shuffled'];
            
            %Naive
            bar_offset = bar_offset +1;
            per_ii=2;
            these_pcorr=hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data;
            mean_hippo_peak_naive(fileNo)=mean(these_pcorr);
            all_hippo_peak_naive=[all_hippo_peak_naive these_pcorr];
            
            bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
            
            
            CI = bootci(1000, {@mean, these_pcorr},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
            glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
            glm_disc.odor_pair(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=fileNo*ones(1,length(these_pcorr));
            glm_ii_disc=glm_ii_disc+length(these_pcorr);
            
            
            ii_stats=ii_stats+1;
            p_corr_stats(ii_stats).data=these_pcorr;
            p_corr_stats(ii_stats).description=[file_legend{fileNo} ' ' ...
                'Naive'];
            
            
            %Proficient
            bar_offset = bar_offset +1;
            per_ii=1;
            these_pcorr=hippo_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data;
            mean_hippo_peak_proficient(fileNo)=mean(these_pcorr);
            all_hippo_peak_proficient=[all_hippo_peak_proficient these_pcorr];
            
            bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
            
            
            CI = bootci(1000, {@mean, these_pcorr},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
            glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=3*ones(1,length(these_pcorr));
            glm_disc.odor_pair(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=fileNo*ones(1,length(these_pcorr));
            glm_ii_disc=glm_ii_disc+length(these_pcorr);
            
            
            ii_stats=ii_stats+1;
            p_corr_stats(ii_stats).data=these_pcorr;
            p_corr_stats(ii_stats).description=[file_legend{fileNo} ' ' ...
                'Proficient'];
            
            bar_offset = bar_offset +1;
            
            
            
        else
            bar_offset = bar_offset +7;
        end
        
    end
    
    plot([0 32],[50 50],'-k')
    
    ylim([40 100])
    
    xticks([2 6 10 14 18 22 26 30])
    xticklabels({file_legend{1}, file_legend{2},file_legend{3}, file_legend{4},file_legend{5}, file_legend{6},file_legend{7}, file_legend{8}})
    
    text(0,98,'Shuffled','Color',[120/255 120/255 120/255])
    text(0,96,'Naive','Color',[238/255 111/255 179/255])
    text(0,96,'Proficient','Color',[158/255 31/255 99/255])
    
    title(['LDA percent correct for hippocampus for peak theta/' PACnames{PACii} ])
    ylabel('Percent correct')
    
    %Perform the glm
    fprintf(1, ['glm for performance for hippocampus peak theta/' PACnames{PACii} '\n'])
    tbl = table(glm_disc.data',glm_disc.group',glm_disc.odor_pair',...
        'VariableNames',{'performance','group','odor_pair'});
    mdl = fitglm(tbl,'performance~group+odor_pair+group*odor_pair'...
        ,'CategoricalVars',[2,3])
    
    %Do ranksum/t test
    [output_data] = drgMutiRanksumorTtest(p_corr_stats);
    
    
    %     glm_pcorr_dt_peak=[];
    %     glm_ii_pcdt_p=0;
    %     glm_pcorr_dt_trough=[];
    %     glm_ii_pcdt_t=0;
    
    
    %Bar graph plot for peak percent correct
    
    %Hippocampus
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
    
    all_ps=[];
    
    glm_disc=[];
    glm_ii_disc=0;
    
    p_corr_stats=[];
    ii_stats=0;
    
    mean_hippo_trough_naive=[];
    mean_hippo_trough_proficient=[];
    
    all_hippo_trough_shuffled=[];
    all_hippo_trough_naive=[];
    all_hippo_trough_proficient=[];
    
    
    for fileNo=1:8
        
        if ~isempty(hippFileName{fileNo})
            load([hippPathName hippFileName{fileNo}])
            hippo_pcorr_out=pcorr_out;
            
            
            grNo=1;
            
            
            
            %Shuffled
            bar_offset = bar_offset +1;
            per_ii=1;
            these_pcorr=hippo_pcorr_out.PACii(PACii).shuffled.pcorr(per_ii).group(grNo).odor_data;
            per_ii=2;
            these_pcorr=hippo_pcorr_out.PACii(PACii).shuffled.pcorr(per_ii).group(grNo).odor_data;
            all_hippo_trough_shuffled=[all_hippo_trough_shuffled these_pcorr];
            
            bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[120/255 120/255 120/255])
            
            
            CI = bootci(1000, {@mean, these_pcorr},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
            glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=1*ones(1,length(these_pcorr));
            glm_disc.odor_pair(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=fileNo*ones(1,length(these_pcorr));
            glm_ii_disc=glm_ii_disc+length(these_pcorr);
            
            
            ii_stats=ii_stats+1;
            p_corr_stats(ii_stats).data=these_pcorr;
            p_corr_stats(ii_stats).description=[file_legend{fileNo} ' ' ...
                'Shuffled'];
            
            %Naive
            bar_offset = bar_offset +1;
            per_ii=2;
            these_pcorr=hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data;
            mean_hippo_trough_naive(fileNo)=mean(these_pcorr);
            all_hippo_trough_naive=[all_hippo_trough_naive these_pcorr];
            
            bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
            
            
            CI = bootci(1000, {@mean, these_pcorr},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
            glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
            glm_disc.odor_pair(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=fileNo*ones(1,length(these_pcorr));
            glm_ii_disc=glm_ii_disc+length(these_pcorr);
            
            
            ii_stats=ii_stats+1;
            p_corr_stats(ii_stats).data=these_pcorr;
            p_corr_stats(ii_stats).description=[file_legend{fileNo} ' ' ...
                'Naive'];
            
            
            %Proficient
            bar_offset = bar_offset +1;
            per_ii=1;
            these_pcorr=hippo_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data;
            mean_hippo_trough_proficient(fileNo)=mean(these_pcorr);
            all_hippo_trough_proficient=[all_hippo_trough_proficient these_pcorr];
            
            bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
            
            
            CI = bootci(1000, {@mean, these_pcorr},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
            glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=3*ones(1,length(these_pcorr));
            glm_disc.odor_pair(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=fileNo*ones(1,length(these_pcorr));
            glm_ii_disc=glm_ii_disc+length(these_pcorr);
            
            
            ii_stats=ii_stats+1;
            p_corr_stats(ii_stats).data=these_pcorr;
            p_corr_stats(ii_stats).description=[file_legend{fileNo} ' ' ...
                'Proficient'];
            
            bar_offset = bar_offset +1;
            
            
            
        else
            bar_offset = bar_offset +7;
        end
        
    end
    
    plot([0 32],[50 50],'-k')
    
    ylim([40 100])
    
    xticks([2 6 10 14 18 22 26 30])
    xticklabels({file_legend{1}, file_legend{2},file_legend{3}, file_legend{4},file_legend{5}, file_legend{6},file_legend{7}, file_legend{8}})
    
    text(0,98,'Shuffled','Color',[120/255 120/255 120/255])
    text(0,96,'Naive','Color',[238/255 111/255 179/255])
    text(0,96,'Proficient','Color',[158/255 31/255 99/255])
    
    title(['LDA percent correct for hippocampus for trough theta/' PACnames{PACii} ])
    ylabel('Percent correct')
    
    %Perform the glm
    fprintf(1, ['glm for performance for hippocampus trough theta/' PACnames{PACii} '\n'])
    tbl = table(glm_disc.data',glm_disc.group',glm_disc.odor_pair',...
        'VariableNames',{'performance','group','odor_pair'});
    mdl = fitglm(tbl,'performance~group+odor_pair+group*odor_pair'...
        ,'CategoricalVars',[2,3])
    
    %Do ranksum/t test
    [output_data] = drgMutiRanksumorTtest(p_corr_stats);
    
    
    %     glm_pcorr_dt_peak=[];
    %     glm_ii_pcdt_p=0;
    %     glm_pcorr_dt_trough=[];
    %     glm_ii_pcdt_t=0;
    
    
    %Bar graph plot for peak percent correct
    
    %Prefrontal
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
    
    all_ps=[];
    
    glm_disc=[];
    glm_ii_disc=0;
    
    p_corr_stats=[];
    ii_stats=0;
    
    mean_pref_peak_naive=[];
    mean_pref_peak_proficient=[];
    all_pref_peak_naive=[];
    all_pref_peak_proficient=[];
    all_pref_peak_shuffled=[];
    
    for fileNo=1:8
        
        if ~isempty(preFileName{fileNo})
            load([prePathName preFileName{fileNo}])
            pref_pcorr_out=pcorr_out;
            
            
            grNo=1;
            
            
            
            %Shuffled
            bar_offset = bar_offset +1;
            per_ii=1;
            these_pcorr=pref_pcorr_out.PACii(PACii).shuffled.pcorr(per_ii).group(grNo).odor_data;
            per_ii=2;
            these_pcorr=pref_pcorr_out.PACii(PACii).shuffled.pcorr(per_ii).group(grNo).odor_data;
            all_pref_peak_shuffled=[all_pref_peak_shuffled these_pcorr];
                
            bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[120/255 120/255 120/255])
            
            
            CI = bootci(1000, {@mean, these_pcorr},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
            glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=1*ones(1,length(these_pcorr));
            glm_disc.odor_pair(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=fileNo*ones(1,length(these_pcorr));
            glm_ii_disc=glm_ii_disc+length(these_pcorr);
            
            
            ii_stats=ii_stats+1;
            p_corr_stats(ii_stats).data=these_pcorr;
            p_corr_stats(ii_stats).description=[file_legend{fileNo} ' ' ...
                'Shuffled'];
            
            %Naive
            bar_offset = bar_offset +1;
            per_ii=2;
            these_pcorr=pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data;
            mean_pref_peak_naive(fileNo)=mean(these_pcorr);
            all_pref_peak_naive=[all_pref_peak_naive these_pcorr];
            
            bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
            
            
            CI = bootci(1000, {@mean, these_pcorr},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
            glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
            glm_disc.odor_pair(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=fileNo*ones(1,length(these_pcorr));
            glm_ii_disc=glm_ii_disc+length(these_pcorr);
            
            
            ii_stats=ii_stats+1;
            p_corr_stats(ii_stats).data=these_pcorr;
            p_corr_stats(ii_stats).description=[file_legend{fileNo} ' ' ...
                'Naive'];
            
            
            %Proficient
            bar_offset = bar_offset +1;
            per_ii=1;
            these_pcorr=pref_pcorr_out.PACii(PACii).peaks.pcorr(per_ii).group(grNo).odor_data;
            mean_pref_peak_proficient(fileNo)=mean(these_pcorr);
            all_pref_peak_proficient=[all_pref_peak_proficient these_pcorr];
            
            bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
            
            
            CI = bootci(1000, {@mean, these_pcorr},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
            glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=3*ones(1,length(these_pcorr));
            glm_disc.odor_pair(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=fileNo*ones(1,length(these_pcorr));
            glm_ii_disc=glm_ii_disc+length(these_pcorr);
            
            
            ii_stats=ii_stats+1;
            p_corr_stats(ii_stats).data=these_pcorr;
            p_corr_stats(ii_stats).description=[file_legend{fileNo} ' ' ...
                'Proficient'];
            
            bar_offset = bar_offset +1;
            
            
            
        else
            bar_offset = bar_offset +7;
        end
        
    end
    
    plot([0 32],[50 50],'-k')
    
    ylim([40 100])
    
    xticks([2 6 10 14 18 22 26 30])
    xticklabels({file_legend{1}, file_legend{2},file_legend{3}, file_legend{4},file_legend{5}, file_legend{6},file_legend{7}, file_legend{8}})
    
    text(0,98,'Shuffled','Color',[120/255 120/255 120/255])
    text(0,96,'Naive','Color',[238/255 111/255 179/255])
    text(0,96,'Proficient','Color',[158/255 31/255 99/255])
    
    title(['LDA percent correct for prefrontal for peak theta/' PACnames{PACii} ])
    ylabel('Percent correct')
    
    %Perform the glm
    fprintf(1, ['glm for performance for prefrontal peak theta/' PACnames{PACii} '\n'])
    tbl = table(glm_disc.data',glm_disc.group',glm_disc.odor_pair',...
        'VariableNames',{'performance','group','odor_pair'});
    mdl = fitglm(tbl,'performance~group+odor_pair+group*odor_pair'...
        ,'CategoricalVars',[2,3])
    
    %Do ranksum/t test
    [output_data] = drgMutiRanksumorTtest(p_corr_stats);
    
    
    
    %     glm_pcorr_dt_peak=[];
    %     glm_ii_pcdt_p=0;
    %     glm_pcorr_dt_trough=[];
    %     glm_ii_pcdt_t=0;
    
    
    %Bar graph plot for peak percent correct
    
    %Prefrontal trough
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
    
    all_ps=[];
    
    glm_disc=[];
    glm_ii_disc=0;
    
    p_corr_stats=[];
    ii_stats=0;
    
    mean_pref_trough_naive=[];
    mean_pref_trough_proficient=[];
    all_pref_trough_naive=[];
    all_pref_trough_proficient=[];
    all_pref_trough_shuffled=[];
    
    for fileNo=1:8
        
        if ~isempty(preFileName{fileNo})
            load([prePathName preFileName{fileNo}])
            pref_pcorr_out=pcorr_out;
            
            
            grNo=1;
            
            
            
            %Shuffled
            bar_offset = bar_offset +1;
            per_ii=1;
            these_pcorr=pref_pcorr_out.PACii(PACii).shuffled.pcorr(per_ii).group(grNo).odor_data;
            per_ii=2;
            these_pcorr=pref_pcorr_out.PACii(PACii).shuffled.pcorr(per_ii).group(grNo).odor_data;
            all_pref_trough_shuffled=[all_pref_trough_shuffled these_pcorr];
            
            bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[120/255 120/255 120/255])
            
            
            CI = bootci(1000, {@mean, these_pcorr},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
            glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=1*ones(1,length(these_pcorr));
            glm_disc.odor_pair(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=fileNo*ones(1,length(these_pcorr));
            glm_ii_disc=glm_ii_disc+length(these_pcorr);
            
            
            ii_stats=ii_stats+1;
            p_corr_stats(ii_stats).data=these_pcorr;
            p_corr_stats(ii_stats).description=[file_legend{fileNo} ' ' ...
                'Shuffled'];
            
            %Naive
            bar_offset = bar_offset +1;
            per_ii=2;
            these_pcorr=pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data;
            mean_pref_trough_naive(fileNo)=mean(these_pcorr);
            all_pref_trough_naive=[all_pref_trough_naive these_pcorr];
            
            
            bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
            
            
            CI = bootci(1000, {@mean, these_pcorr},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
            glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
            glm_disc.odor_pair(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=fileNo*ones(1,length(these_pcorr));
            glm_ii_disc=glm_ii_disc+length(these_pcorr);
            
            
            ii_stats=ii_stats+1;
            p_corr_stats(ii_stats).data=these_pcorr;
            p_corr_stats(ii_stats).description=[file_legend{fileNo} ' ' ...
                'Naive'];
            
            
            %Proficient
            bar_offset = bar_offset +1;
            per_ii=1;
            these_pcorr=pref_pcorr_out.PACii(PACii).trough.pcorr(per_ii).group(grNo).odor_data;
            mean_pref_trough_proficient(fileNo)=mean(these_pcorr);
            all_pref_trough_proficient=[all_pref_trough_proficient these_pcorr];
            
            bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
            
            
            CI = bootci(1000, {@mean, these_pcorr},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
            glm_disc.group(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=3*ones(1,length(these_pcorr));
            glm_disc.odor_pair(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=fileNo*ones(1,length(these_pcorr));
            glm_ii_disc=glm_ii_disc+length(these_pcorr);
            
            
            ii_stats=ii_stats+1;
            p_corr_stats(ii_stats).data=these_pcorr;
            p_corr_stats(ii_stats).description=[file_legend{fileNo} ' ' ...
                'Proficient'];
            
            bar_offset = bar_offset +1;
            
            
            
        else
            bar_offset = bar_offset +7;
        end
        
    end
    
    plot([0 32],[50 50],'-k')
    
    ylim([40 100])
    
    xticks([2 6 10 14 18 22 26 30])
    xticklabels({file_legend{1}, file_legend{2},file_legend{3}, file_legend{4},file_legend{5}, file_legend{6},file_legend{7}, file_legend{8}})
    
    text(0,98,'Shuffled','Color',[120/255 120/255 120/255])
    text(0,96,'Naive','Color',[238/255 111/255 179/255])
    text(0,96,'Proficient','Color',[158/255 31/255 99/255])
    
    title(['LDA percent correct for prefrontal for trough theta/' PACnames{PACii} ])
    ylabel('Percent correct')
    
    %Perform the glm
    fprintf(1, ['glm for performance for prefrontal trough theta/' PACnames{PACii} '\n'])
    tbl = table(glm_disc.data',glm_disc.group',glm_disc.odor_pair',...
        'VariableNames',{'performance','group','odor_pair'});
    mdl = fitglm(tbl,'performance~group+odor_pair+group*odor_pair'...
        ,'CategoricalVars',[2,3])
    
    %Do ranksum/t test
    [output_data] = drgMutiRanksumorTtest(p_corr_stats);
    
    
    %Now plot the bar graphs for the mean percent correct per mouse
    
    glm_disc=[];
    glm_ii_disc=0;
    
    p_corr_stats=[];
    ii_stats=0;
    
    %Hippocampus peak
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .2 .3 .4])
    hold on
    
    
    bar_offset = 0;
    
    
    %Naive
    bar_offset = bar_offset +1;
    
    
    these_pcorr=mean_hippo_peak_naive;
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
    
    
    CI = bootci(1000, {@mean, these_pcorr},'type','cper');
    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
    plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
    
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=zeros(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Hippo peak naive'];
    
    
    %Proficient
    bar_offset = bar_offset +1;
    these_pcorr=mean_hippo_peak_proficient;
    
    
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
    
    
    CI = bootci(1000, {@mean, these_pcorr},'type','cper');
    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
    plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
    
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Hippo peak proficient'];
    
    
    plot([0 3],[50 50],'-k')
    
    ylim([40 100])
    xlim([0.3 2.8])
    
    xticks([1 2])
    xticklabels({'Naive', 'Proficient'})
    
    title(['LDA mean percent correct for hippocampus for peak theta/' PACnames{PACii} ])
    ylabel('Percent correct')
    
    %Hippocampus trough
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .2 .3 .4])
    hold on
    
    
    bar_offset = 0;
    
    
    %Naive
    bar_offset = bar_offset +1;
    
    
    these_pcorr=mean_hippo_trough_naive;
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
    
    
    CI = bootci(1000, {@mean, these_pcorr},'type','cper');
    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
    plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
    
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=zeros(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=zeros(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Hippo trough naive'];
    
    
    %Proficient
    bar_offset = bar_offset +1;
    these_pcorr=mean_hippo_trough_proficient;
    
    
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
    
    
    CI = bootci(1000, {@mean, these_pcorr},'type','cper');
    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
    plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
    
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=zeros(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Hippo trough proficient'];
    
    
    plot([0 3],[50 50],'-k')
    
    ylim([40 100])
    xlim([0.3 2.8])
    
    xticks([1 2])
    xticklabels({'Naive', 'Proficient'})
    
    title(['LDA mean percent correct for hippocampus for trough theta/' PACnames{PACii} ])
    ylabel('Percent correct')
    
    %Prefrontal peak
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .2 .3 .4])
    hold on
    
    
    bar_offset = 0;
    
    
    %Naive
    bar_offset = bar_offset +1;
    
    
    these_pcorr=mean_pref_peak_naive;
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
    
    
    CI = bootci(1000, {@mean, these_pcorr},'type','cper');
    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
    plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
    
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=zeros(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Hippo peak naive'];
    
    
    %Proficient
    bar_offset = bar_offset +1;
    these_pcorr=mean_pref_peak_proficient;
    
    
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
    
    
    CI = bootci(1000, {@mean, these_pcorr},'type','cper');
    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
    plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
    
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Hippo peak proficient'];
    
    
    plot([0 3],[50 50],'-k')
    
    ylim([40 100])
    xlim([0.3 2.8])
    
    xticks([1 2])
    xticklabels({'Naive', 'Proficient'})
    
    title(['LDA mean percent correct for prefrontal for peak theta/' PACnames{PACii} ])
    ylabel('Percent correct')
    
    %Prefrontal trough
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .2 .3 .4])
    hold on
    
    
    bar_offset = 0;
    
    
    %Naive
    bar_offset = bar_offset +1;
    
    
    these_pcorr=mean_pref_trough_naive;
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
    
    
    CI = bootci(1000, {@mean, these_pcorr},'type','cper');
    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
    plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
    
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=zeros(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=zeros(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Hippo trough naive'];
    
    
    %Proficient
    bar_offset = bar_offset +1;
    these_pcorr=mean_pref_trough_proficient;
    
    
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
    
    
    CI = bootci(1000, {@mean, these_pcorr},'type','cper');
    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
    plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
    
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=zeros(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Hippo trough proficient'];
    
    
    plot([0 3],[50 50],'-k')
    
    ylim([40 100])
    xlim([0.3 2.8])
    
    xticks([1 2])
    xticklabels({'Naive', 'Proficient'})
    
    title(['LDA mean percent correct for prefrontal for trough theta/' PACnames{PACii} ])
    ylabel('Percent correct')
    
    
    %Perform the glm
    fprintf(1, ['glm for mean performance for theta/' PACnames{PACii} '\n'])
    tbl = table(glm_disc.data',glm_disc.peak',glm_disc.proficient',glm_disc.brain_region',...
        'VariableNames',{'performance','peak','proficient','brain_region'});
    mdl = fitglm(tbl,'performance~peak+proficient+brain_region+peak*proficient*brain_region'...
        ,'CategoricalVars',[2,3,4])
    
    %Do ranksum/t test
    [output_data] = drgMutiRanksumorTtest(p_corr_stats);
    
    %Finally, do the plots for percent correct per odor pair per mouse
    %Now plot the bar graphs for the mean percent correct per mouse
    
    glm_disc=[];
    glm_ii_disc=0;
    
    p_corr_stats=[];
    ii_stats=0;
    
    
    
    %Hippocampus peak
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
    
    %Shuffled
    bar_offset = bar_offset +1;
    
    
    these_pcorr=all_hippo_peak_shuffled;
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[120/255 120/255 120/255])
     
    [mean_out, CIout]=drgViolinPoint(these_pcorr,edges,bar_offset,rand_offset,'k','k',3);
    
%     CI = bootci(1000, {@mean, these_pcorr},'type','cper');
%     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%     plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%     
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=3*ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Hippo peak shuffled'];
    
    %Naive
    bar_offset = bar_offset +1;
    
    
    these_pcorr=all_hippo_peak_naive;
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
    
    [mean_out, CIout]=drgViolinPoint(these_pcorr,edges,bar_offset,rand_offset,'k','k',3);
    
%     CI = bootci(1000, {@mean, these_pcorr},'type','cper');
%     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%     plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%     
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Hippo peak naive'];
    
    
    %Proficient
    bar_offset = bar_offset +1;
    these_pcorr=all_hippo_peak_proficient;
    
    
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
    
    
    [mean_out, CIout]=drgViolinPoint(these_pcorr,edges,bar_offset,rand_offset,'k','k',3);
%     CI = bootci(1000, {@mean, these_pcorr},'type','cper');
%     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%     plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%     
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Hippo peak proficient'];
    
    
    plot([0 4],[50 50],'-k')
    
    ylim([40 100])
    xlim([0.3 3.8])
    
    xticks([1 2 3])
    xticklabels({'Shuffled', 'Naive', 'Proficient'})
    
    title(['LDA percent correct for hippocampus per odor pair per mouse for peak theta/' PACnames{PACii} ])
    ylabel('Percent correct')
    
    %Hippocampus trough
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .2 .3 .4])
    hold on
    
    
    bar_offset = 0;
    
    ax=gca;ax.LineWidth=3;
    
    %Shuffled
    bar_offset = bar_offset +1;
    
    
    these_pcorr=all_hippo_trough_shuffled;
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[120/255 120/255 120/255])
    
    [mean_out, CIout]=drgViolinPoint(these_pcorr,edges,bar_offset,rand_offset,'k','k',3);
    
%     CI = bootci(1000, {@mean, these_pcorr},'type','cper');
%     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%     plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%     
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=zeros(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=3*ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Hippo trough shuffled'];
    
    %Naive
    bar_offset = bar_offset +1;
    
    
    these_pcorr=all_hippo_trough_naive;
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
    
    [mean_out, CIout]=drgViolinPoint(these_pcorr,edges,bar_offset,rand_offset,'k','k',3);
    
%     CI = bootci(1000, {@mean, these_pcorr},'type','cper');
%     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%     plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%     
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=zeros(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Hippo trough naive'];
    
    
    %Proficient
    bar_offset = bar_offset +1;
    these_pcorr=all_hippo_trough_proficient;
    
    
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
    
    [mean_out, CIout]=drgViolinPoint(these_pcorr,edges,bar_offset,rand_offset,'k','k',3);
     
%     CI = bootci(1000, {@mean, these_pcorr},'type','cper');
%     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%     plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
    
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=zeros(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Hippo trough proficient'];
    
    
    plot([0 4],[50 50],'-k')
    
    ylim([40 100])
    xlim([0.3 3.8])
    
    xticks([1 2 3])
    xticklabels({'Shuffled', 'Naive', 'Proficient'})
    
    title(['LDA percent correct for hippocampus per mouse per odor pair for trough theta/' PACnames{PACii} ])
    ylabel('Percent correct')
    
    %Prefrontal peak
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.1 .2 .3 .4])
    hold on
    
    
    bar_offset = 0;
    
     ax=gca;ax.LineWidth=3;
     
     
    %Shuffled
    bar_offset = bar_offset +1;
    
    
    these_pcorr=all_pref_peak_shuffled;
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[120/255 120/255 120/255])
    
    
    [mean_out, CIout]=drgViolinPoint(these_pcorr,edges,bar_offset,rand_offset,'k','k',3);
    
%     CI = bootci(1000, {@mean, these_pcorr},'type','cper');
%     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%     plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
    
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=3*ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Pref peak shuffled'];
    
     
    %Naive
    bar_offset = bar_offset +1;
    
    
    these_pcorr=all_pref_peak_naive;
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
    
    [mean_out, CIout]=drgViolinPoint(these_pcorr,edges,bar_offset,rand_offset,'k','k',3);
    
%     CI = bootci(1000, {@mean, these_pcorr},'type','cper');
%     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%     plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%     
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Pref peak naive'];
    
    
    %Proficient
    bar_offset = bar_offset +1;
    these_pcorr=all_pref_peak_proficient;
    
    
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
    
    [mean_out, CIout]=drgViolinPoint(these_pcorr,edges,bar_offset,rand_offset,'k','k',3);
    
    
%     CI = bootci(1000, {@mean, these_pcorr},'type','cper');
%     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%     plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
    
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Pref peak proficient'];
    
    
    plot([0 4],[50 50],'-k')
    
    ylim([40 100])
    xlim([0.3 3.8])
    
    xticks([1 2 3])
    xticklabels({'Shuffled','Naive', 'Proficient'})
    
    title(['LDA percent correct per mouse per odor pair for prefrontal for peak theta/' PACnames{PACii} ])
    ylabel('Percent correct')
    
    %Prefrontal trough
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
    
    %Naive
    bar_offset = bar_offset +1;
    
    
    these_pcorr=all_pref_trough_shuffled;
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[120/255 120/255 120/255])
    
    
    [mean_out, CIout]=drgViolinPoint(these_pcorr,edges,bar_offset,rand_offset,'k','k',3);
    
%     CI = bootci(1000, {@mean, these_pcorr},'type','cper');
%     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%     plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%     
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=zeros(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=3*ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Pref trough shuffled'];
    
    %Naive
    bar_offset = bar_offset +1;
    
    
    these_pcorr=all_pref_trough_naive;
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
    
     [mean_out, CIout]=drgViolinPoint(these_pcorr,edges,bar_offset,rand_offset,'k','k',3);
     
%     CI = bootci(1000, {@mean, these_pcorr},'type','cper');
%     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%     plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
    
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=zeros(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Pref trough naive'];
    
    
    %Proficient
    bar_offset = bar_offset +1;
    these_pcorr=all_pref_trough_proficient;
    
    
    bar(bar_offset,mean(these_pcorr),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
    
    [mean_out, CIout]=drgViolinPoint(these_pcorr,edges,bar_offset,rand_offset,'k','k',3);
     
%     
%     CI = bootci(1000, {@mean, these_pcorr},'type','cper');
%     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%     plot(bar_offset*ones(1,length(these_pcorr)),these_pcorr,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%     
    glm_disc.data(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=these_pcorr;
    glm_disc.peak(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=zeros(1,length(these_pcorr));
    glm_disc.proficient(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=ones(1,length(these_pcorr));
    glm_disc.brain_region(glm_ii_disc+1:glm_ii_disc+length(these_pcorr))=2*ones(1,length(these_pcorr));
    glm_ii_disc=glm_ii_disc+length(these_pcorr);
    
    
    ii_stats=ii_stats+1;
    p_corr_stats(ii_stats).data=these_pcorr;
    p_corr_stats(ii_stats).description=['Pref trough proficient'];
    
    
    plot([0 4],[50 50],'-k')
    
    ylim([40 100])
    xlim([0.3 3.8])
    
    xticks([1 2 3])
    xticklabels({'Shuffled', 'Naive', 'Proficient'})
    
    title(['LDApercent correct for prefrontal per mouse per odor pair for trough theta/' PACnames{PACii} ])
    ylabel('Percent correct')
    
    
    %Perform the glm
    fprintf(1, ['glm for percent correct per mouse per odor pair for theta/' PACnames{PACii} '\n'])
    tbl = table(glm_disc.data',glm_disc.peak',glm_disc.proficient',glm_disc.brain_region',...
        'VariableNames',{'performance','peak','proficient','brain_region'});
    mdl = fitglm(tbl,'performance~peak+proficient+brain_region+peak*proficient*brain_region'...
        ,'CategoricalVars',[2,3,4])
    
    %Do ranksum/t test
    [output_data] = drgMutiRanksumorTtest(p_corr_stats);
    
    pffft=1;
end
pfft=1;
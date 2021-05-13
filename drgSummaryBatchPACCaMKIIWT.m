function drgSummaryBatchPACCaMKIIWT
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'drgbDiscPar' file listing 'Discriminant_*.mat' output files
%from drgLFPDiscriminantBatch
%
%Performs summary analises for LDA and  PAC analysis performed with case 19
%of drgAnalysisBatchLFP

warning('off')

close all
clear all

%If you want statistics to be done with the value for each odor pair for
%each mouse make this variable 1
mouse_op=1;

PACnames{1}='Beta';
PACnames{2}='Low gamma';
PACnames{3}='High gamma';

prof_naive_leg{1}='Proficient';
prof_naive_leg{2}='Naive';

group_legend{1}='WT';
group_legend{2}='Het';
group_legend{3}='KO';

evTypeLabels{1}='S+';
evTypeLabels{2}='S-';

% %Location of files for proficient 85-100
% PathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/LFPPAC drgAnalysisBatchLFPCaMKII case 19 output take 2/';
%
% %Hippocampus
% hippFileName{1}='CaMKIIacetoLFPPACall121319_hippocampusLFP2.mat';
% hippFileName{2}='CaMKIIethylbenLFPPACal222019_hippocampusLFP2.mat';
% hippFileName{3}='CaMKIIethylaceLFPPACall12219_hippocampusLFP2.mat';
% hippFileName{4}='CaMKIIpropylaceLFPPACall1220_hippocampusLFP2.mat';
% hippFileName{5}='CaMKIIpointzero1ethylacepropylaceLFPPACall122019_hippocampusLFP2.mat';
% hippFileName{6}='CaMKIIpointzero1propylaceLFPPACall11620_hippocampusLFP2.mat';
% hippFileName{7}='CaMKIIpzz1ethylaceLFPPACall10719_hippocampusLFP2.mat';
% hippFileName{8}='CaMKIIpzz1propylaceLFPPACall121919_hippocampusLFP2.mat';
%
% %Prefrontal
% preFileName{1}='CaMKIIacetoLFPPACall121319_prefrontalLFP2.mat';
% preFileName{2}='CaMKIIethylbenLFPPACal222019_prefrontalLFP2.mat';
% preFileName{3}='CaMKIIethylaceLFPPACall12219_prefrontalLFP2.mat';
% preFileName{4}='CaMKIIpropylaceLFPPACall1220_prefrontalLFP2.mat';
% preFileName{5}='CaMKIIpointzero1ethylacepropylaceLFPPACall122019_prefrontalLFP2.mat';
% preFileName{6}='CaMKIIpointzero1propylaceLFPPACall11620_prefrontalLFP2.mat';
% preFileName{7}='CaMKIIpzz1ethylaceLFPPACall10719_prefrontalLFP2.mat';
% preFileName{8}='CaMKIIpzz1propylaceLFPPACall121919_prefrontalLFP2.mat';

%Location of files for proficient 80-100
PathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/LFPPAC drgAnalysisBatchLFPCaMKII case 19 output take 2 80/';

%Hippocampus
hippFileName{1}='CaMKIIacetoLFPPACall121319_hippocampusLFP80.mat';
hippFileName{2}='CaMKIIethylbenLFPPACal222019_hippocampusLFP80.mat';
hippFileName{3}='CaMKIIethylaceLFPPACall12219_hippocampusLFP80.mat';
hippFileName{4}='CaMKIIpropylaceLFPPACall1220_hippocampusLFP80.mat';
hippFileName{5}='CaMKIIpointzero1ethylacepropylaceLFPPACall122019_hippocampusLFP80.mat';
hippFileName{6}='CaMKIIpointzero1propylaceLFPPACall11620_hippocampusLFP80.mat';
hippFileName{7}='CaMKIIpzz1ethylaceLFPPACall10719_hippocampusLFP80.mat';
hippFileName{8}='CaMKIIpzz1propylaceLFPPACall121919_hippocampusLFP80.mat';

%Prefrontal
preFileName{1}='CaMKIIacetoLFPPACall121319_prefrontalLFP80.mat';
preFileName{2}='CaMKIIethylbenLFPPACal222019_prefrontalLFP80.mat';
preFileName{3}='CaMKIIethylaceLFPPACall12219_prefrontalLFP80.mat';
preFileName{4}='CaMKIIpropylaceLFPPACall1220_prefrontalLFP80.mat';
preFileName{5}='CaMKIIpointzero1ethylacepropylaceLFPPACall122019_prefrontalLFP80.mat';
preFileName{6}='CaMKIIpointzero1propylaceLFPPACall11620_prefrontalLFP80.mat';
preFileName{7}='CaMKIIpzz1ethylaceLFPPACall10719_prefrontalLFP80.mat';
preFileName{8}='CaMKIIpzz1propylaceLFPPACall121919_prefrontalLFP80.mat';


%Load data hippocampus
all_hippo=[];

for ii=1:length(hippFileName)
    load([PathName hippFileName{ii}])
    all_hippo(ii).handles_out=handles_out;
end

%Load data prefrontal
all_pre=[];

for ii=1:length(preFileName)
    load([PathName preFileName{ii}])
    all_pre(ii).handles_out=handles_out;
end

figNo=0;




for pacii=[1 3]    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    %Now plot for the hippocampus the average MI for each electrode calculated per mouse
    %(including all sessions for each mouse)
    
    glm_mi_hipp=[];
    glm_ii_hipp=0;
    
    glm_mi_both=[];
    glm_mi_ii_both=0;
    
    glm_PA_both=[];
    glm_PA_ii_both=0;
    
    id_ii=0;
    input_data=[];
    
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    %             try
    %                 close(figNo+pacii)
    %             catch
    %             end
    %             hFig=figure(figNo+pacii);
    
    set(hFig, 'units','normalized','position',[.1 .5 .3 .4])
    hold on
    
    ax=gca;ax.LineWidth=3;
    
    bar_lab_loc=[];
    no_ev_labels=0;
    ii_gr_included=0;
    bar_offset = 0;
    edges=[0:0.001:0.05];
    rand_offset=0.5;
    
    for evNo=1:2
        
        for per_ii=2:-1:1
            
            grNo=1;
            
            bar_offset = bar_offset + 1;
            
            %Get these MI values
            these_mi=[];
            ii_mi=0;
            for ii=1:length(hippFileName)
                this_jj=[];
                for jj=1:all_hippo(ii).handles_out.mi_ii
                    if all_hippo(ii).handles_out.mi_values(jj).pacii==pacii
                        if all_hippo(ii).handles_out.mi_values(jj).evNo==evNo
                            if all_hippo(ii).handles_out.mi_values(jj).per_ii==per_ii
                                if all_hippo(ii).handles_out.mi_values(jj).groupNo==grNo
                                    this_jj=jj;
                                end
                                
                            end
                            
                        end
                    end
                    
                end
                if ~isempty(this_jj)
                    if mouse_op==1
                        these_mi(ii_mi+1:ii_mi+length(all_hippo(ii).handles_out.mi_values(this_jj).MI_per_mouse))=all_hippo(ii).handles_out.mi_values(this_jj).MI_per_mouse;
                        ii_mi=ii_mi+length(all_hippo(ii).handles_out.mi_values(this_jj).MI_per_mouse);
                    else
                        ii_mi=ii_mi+1;
                        these_mi(ii_mi)=all_hippo(ii).handles_out.mi_values(this_jj).MI;
                    end
                end
            end
            
            if evNo==2
                if per_ii==1
                    %S+ Proficient
                    bar(bar_offset,mean(these_mi),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                else
                    %S- Naive
                    bar(bar_offset,mean(these_mi),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                end
            else
                if per_ii==1
                    %S- Proficient
                    bar(bar_offset,mean(these_mi),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                else
                    %S- naive
                    bar(bar_offset,mean(these_mi),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                end
            end
            
            
            %Violin plot
            
            [mean_out, CIout]=drgViolinPoint(these_mi,edges,bar_offset,rand_offset,'k','k',3);
            %             CI = bootci(1000, {@mean, these_mi},'type','cper');
            %             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            %             plot(bar_offset*ones(1,length(these_mi)),these_mi,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            %
            
            %                                 %Save data for glm and ranksum
            
            glm_mi_hipp.data(glm_ii_hipp+1:glm_ii_hipp+length(these_mi))=these_mi;
            glm_mi_hipp.perCorr(glm_ii_hipp+1:glm_ii_hipp+length(these_mi))=per_ii*ones(1,length(these_mi));
            glm_mi_hipp.event(glm_ii_hipp+1:glm_ii_hipp+length(these_mi))=evNo*ones(1,length(these_mi));
            glm_ii_hipp=glm_ii_hipp+length(these_mi);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=these_mi;
            input_data(id_ii).description=[evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
            
            
            glm_mi_both.data(glm_mi_ii_both+1:glm_mi_ii_both+length(these_mi))=these_mi;
            glm_mi_both.brain_region(glm_mi_ii_both+1:glm_mi_ii_both+length(these_mi))=ones(1,length(these_mi));
            glm_mi_both.perCorr(glm_mi_ii_both+1:glm_mi_ii_both+length(these_mi))=per_ii*ones(1,length(these_mi));
            glm_mi_both.event(glm_mi_ii_both+1:glm_mi_ii_both+length(these_mi))=evNo*ones(1,length(these_mi));
            glm_mi_ii_both=glm_mi_ii_both+length(these_mi);
            
            %             end
            
            
        end
        bar_offset = bar_offset + 1;
        
    end
    
    title(['Average MI per mouse per odor pair for PAC theta/' PACnames{pacii} ' hippocampus'])
    
    
    %Proficient/Naive annotations
    annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color',[158/255 31/255 99/255],'LineStyle','none');
    annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color',[0 114/255 178/255],'LineStyle','none');
    
    
    xticks([1 2 4 5])
    xticklabels({'nS+', 'pS+','nS-', 'pS-'})
    
    yticks([0 0.01 0.02])
    
    ylim([0 0.02])
    
    ylabel('Modulation Index')
    
    %Perform the glm
    fprintf(1, ['glm for average MI per mouse per odor pair for PAC theta' PACnames{pacii} ' hippocampus\n'])
    
    fprintf(1, ['\n\nglm for MI for Theta/' PACnames{pacii} '\n'])
    tbl = table(glm_mi_hipp.data',glm_mi_hipp.perCorr',glm_mi_hipp.event',...
        'VariableNames',{'MI','perCorr','event'});
    mdl = fitglm(tbl,'MI~perCorr+event+perCorr*event'...
        ,'CategoricalVars',[2,3])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for average MI per mouse per odor pair for PAC theta' PACnames{pacii} ' hippocampus\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
    
    
    %Now plot the average PA variance for the hippocampus for each electrode calculated per mouse
    %(including all sessions for each mouse)
    
    
    
    
    glm_PA_hipp=[];
    glm_ii=0;
    
    id_ii=0;
    input_data=[];
    
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .3 .4])
    hold on
    
    bar_lab_loc=[];
    no_ev_labels=0;
    ii_gr_included=0;
    bar_offset = 0;
    edges=[0:20:1000];
    rand_offset=0.5;
    
    %Subtract the average for S+ proficient
    %Get these MI values
    ref_PA=[];
    ii_PA=0;
    evNo=1;
    per_ii=1;
    for ii=1:length(hippFileName)
        this_jj=[];
        for jj=1:all_hippo(ii).handles_out.PA_ii
            if all_hippo(ii).handles_out.PA_values(jj).pacii==pacii
                if all_hippo(ii).handles_out.PA_values(jj).evNo==evNo
                    if all_hippo(ii).handles_out.PA_values(jj).per_ii==per_ii
                        if all_hippo(ii).handles_out.PA_values(jj).groupNo==grNo
                            this_jj=jj;
                        end
                        
                    end
                    
                end
            end
            
        end
        
        if ~isempty(this_jj)
            if mouse_op==1
                ref_PA(ii_PA+1:ii_PA+length(all_hippo(ii).handles_out.PA_values(this_jj).PA_var_per_mouse))=all_hippo(ii).handles_out.PA_values(this_jj).PA_var_per_mouse;
                ii_PA=ii_PA+length(all_hippo(ii).handles_out.PA_values(this_jj).PA_var_per_mouse);
            else
                ii_PA=ii_PA+1;
                ref_PA(ii_PA)=all_hippo(ii).handles_out.PA_values(this_jj).PA_var;
            end
        end
    end
    
    
    for evNo=1:2
        
        for per_ii=2:-1:1
            
            grNo=1;
            bar_offset = bar_offset +1;
            
            %Get these MI values
            these_PA=[];
            ii_PA=0;
            for ii=1:length(hippFileName)
                this_jj=[];
                for jj=1:all_hippo(ii).handles_out.PA_ii
                    if all_hippo(ii).handles_out.PA_values(jj).pacii==pacii
                        if all_hippo(ii).handles_out.PA_values(jj).evNo==evNo
                            if all_hippo(ii).handles_out.PA_values(jj).per_ii==per_ii
                                if all_hippo(ii).handles_out.PA_values(jj).groupNo==grNo
                                    this_jj=jj;
                                end
                                
                            end
                            
                        end
                    end
                    
                end
                if ~isempty(this_jj)
                    if mouse_op==1
                        these_PA(ii_PA+1:ii_PA+length(all_hippo(ii).handles_out.PA_values(this_jj).PA_var_per_mouse))=all_hippo(ii).handles_out.PA_values(this_jj).PA_var_per_mouse;
                        ii_PA=ii_PA+length(all_hippo(ii).handles_out.PA_values(this_jj).PA_var_per_mouse);
                    else
                        ii_PA=ii_PA+1;
                        these_PA(ii_PA)=all_hippo(ii).handles_out.PA_values(this_jj).PA_var;
                    end
                end
            end
            
            these_PA=these_PA-mean(ref_PA);
            
            if evNo==2
                if per_ii==1
                    %S+ Proficient
                    bar(bar_offset,mean(these_PA),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                else
                    %S- Naive
                    bar(bar_offset,mean(these_PA),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                end
            else
                if per_ii==1
                    %S- Proficient
                    bar(bar_offset,mean(these_PA),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                else
                    %S- naive
                    bar(bar_offset,mean(these_PA),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                end
            end
            
            
            %Violin plot
            
            [mean_out, CIout]=drgViolinPoint(these_PA,edges,bar_offset,rand_offset,'k','k',3);
            %             CI = bootci(1000, {@mean, these_PA},'type','cper');
            %             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            %             plot(bar_offset*ones(1,length(these_PA)),these_PA,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            %
            %                                 %Save data for glm and ranksum
            
            glm_PA_hipp.data(glm_ii+1:glm_ii+length(these_PA))=these_PA;
            glm_PA_hipp.perCorr(glm_ii+1:glm_ii+length(these_PA))=per_ii*ones(1,length(these_PA));
            glm_PA_hipp.event(glm_ii+1:glm_ii+length(these_PA))=evNo*ones(1,length(these_PA));
            glm_ii=glm_ii+length(these_PA);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=these_PA;
            input_data(id_ii).description=[evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
            
            glm_PA_both.data(glm_PA_ii_both+1:glm_PA_ii_both+length(these_PA))=these_PA;
            glm_PA_both.brain_region(glm_PA_ii_both+1:glm_PA_ii_both+length(these_mi))=ones(1,length(these_mi));
            glm_PA_both.perCorr(glm_PA_ii_both+1:glm_PA_ii_both+length(these_PA))=per_ii*ones(1,length(these_PA));
            glm_PA_both.event(glm_PA_ii_both+1:glm_PA_ii_both+length(these_PA))=evNo*ones(1,length(these_PA));
            glm_PA_ii_both=glm_PA_ii_both+length(these_PA);
            
        end
        bar_offset = bar_offset + 1;
        
    end
    
    title(['Delta PA variance per mouse per odor pair for PAC theta/' PACnames{pacii} ' hippocampus'])
    
    
    
    %Proficient/Naive annotations
    annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color',[158/255 31/255 99/255],'LineStyle','none');
    annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color',[0 114/255 178/255],'LineStyle','none');
    
    plot([0 6],[0 0],'-k','LineWidth',2)
    xlim([0 6])
    
    xticks([1 2 4 5])
    xticklabels({'nS+', 'pS+','nS-', 'pS-'})
    
    yticks([-1000 0 1000])
    
    ylim([-1500  1000])
    
    ylabel('Delta PA variance')
    
    %Perform the glm
    fprintf(1, ['glm for delta PA variance per mouse per odor pair for PAC theta' PACnames{pacii} ' hippocampus\n'])
    
    fprintf(1, ['\n\nglm for delta PA variance for Theta/' PACnames{pacii} '\n'])
    tbl = table(glm_PA_hipp.data',glm_PA_hipp.perCorr',glm_PA_hipp.event',...
        'VariableNames',{'MI','perCorr','event'});
    mdl = fitglm(tbl,'MI~perCorr+event+perCorr*event'...
        ,'CategoricalVars',[2,3])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for delta PA variance per mouse per odor pair for PAC theta' PACnames{pacii} ' hippocampus\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
    
    %Now plot the average MI for prefrontal for each electrode calculated per mouse
    %(including all sessions for each mouse)
    
    
    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_mi_pre=[];
    glm_ii=0;
    
    id_ii=0;
    input_data=[];
    
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    %             try
    %                 close(figNo+pacii)
    %             catch
    %             end
    %             hFig=figure(figNo+pacii);
    
    set(hFig, 'units','normalized','position',[.1 .5 .3 .4])
    hold on
    
    bar_lab_loc=[];
    no_ev_labels=0;
    ii_gr_included=0;
    bar_offset = 0;
    
    edges=[0:0.001:0.05];
    rand_offset=0.5;
    
    for evNo=1:2
        
        for per_ii=2:-1:1
            
            grNo=1;
            
            bar_offset = bar_offset + 1;
            
            %Get these MI values
            these_mi=[];
            ii_mi=0;
            for ii=1:length(hippFileName)
                this_jj=[];
                for jj=1:all_pre(ii).handles_out.mi_ii
                    if all_pre(ii).handles_out.mi_values(jj).pacii==pacii
                        if all_pre(ii).handles_out.mi_values(jj).evNo==evNo
                            if all_pre(ii).handles_out.mi_values(jj).per_ii==per_ii
                                if all_pre(ii).handles_out.mi_values(jj).groupNo==grNo
                                    this_jj=jj;
                                end
                                
                            end
                            
                        end
                    end
                    
                end
                if ~isempty(this_jj)
                    if mouse_op==1
                        these_mi(ii_mi+1:ii_mi+length(all_pre(ii).handles_out.mi_values(this_jj).MI_per_mouse))=all_pre(ii).handles_out.mi_values(this_jj).MI_per_mouse;
                        ii_mi=ii_mi+length(all_pre(ii).handles_out.mi_values(this_jj).MI_per_mouse);
                    else
                        ii_mi=ii_mi+1;
                        these_mi(ii_mi)=all_pre(ii).handles_out.mi_values(this_jj).MI;
                    end
                end
            end
            
            
            
            if evNo==2
                if per_ii==1
                    %S+ Proficient
                    bar(bar_offset,mean(these_mi),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                else
                    %S- Naive
                    bar(bar_offset,mean(these_mi),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                end
            else
                if per_ii==1
                    %S- Proficient
                    bar(bar_offset,mean(these_mi),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                else
                    %S- naive
                    bar(bar_offset,mean(these_mi),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                end
            end
            
            
            
            %Violin plot
            
            [mean_out, CIout]=drgViolinPoint(these_mi,edges,bar_offset,rand_offset,'k','k',3);
            %             CI = bootci(1000, {@mean, these_mi},'type','cper');
            %             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            %             plot(bar_offset*ones(1,length(these_mi)),these_mi,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            glm_mi_pre.data(glm_ii+1:glm_ii+length(these_mi))=these_mi;
            glm_mi_pre.brain_region(glm_ii+1:glm_ii+length(these_mi))=2*ones(1,length(these_mi));
            glm_mi_pre.perCorr(glm_ii+1:glm_ii+length(these_mi))=per_ii*ones(1,length(these_mi));
            glm_mi_pre.event(glm_ii+1:glm_ii+length(these_mi))=evNo*ones(1,length(these_mi));
            glm_ii=glm_ii+length(these_mi);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=these_mi;
            input_data(id_ii).description=[evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
            
            glm_mi_both.data(glm_mi_ii_both+1:glm_mi_ii_both+length(these_mi))=these_mi;
            glm_mi_both.brain_region(glm_mi_ii_both+1:glm_mi_ii_both+length(these_mi))=2*ones(1,length(these_mi));
            glm_mi_both.perCorr(glm_mi_ii_both+1:glm_mi_ii_both+length(these_mi))=per_ii*ones(1,length(these_mi));
            glm_mi_both.event(glm_mi_ii_both+1:glm_mi_ii_both+length(these_mi))=evNo*ones(1,length(these_mi));
            glm_mi_ii_both=glm_mi_ii_both+length(these_mi);
            
            
            
        end
        bar_offset = bar_offset + 1;
        
    end
    
    title(['Average MI per mouse per odor pair for PAC theta/' PACnames{pacii} ' prefrontal'])
    
    
    %Annotations identifying groups
    x_interval=0.8/ii_gr_included;
    for ii=1:ii_gr_included
        annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
    end
    
    %Proficient/Naive annotations
    annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color',[158/255 31/255 99/255],'LineStyle','none');
    annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color',[0 114/255 178/255],'LineStyle','none');
    
    
    xticks([1 2 4 5])
    xticklabels({'nS+', 'pS+','nS-', 'pS-'})
    
    yticks([0 0.01 0.02])
    
    ylim([0 0.02])
    
    ylabel('Modulation Index')
    
    %Perform the glm
    fprintf(1, ['glm for average MI per mouse per odor pair for PAC theta' PACnames{pacii} ' prefrontal\n'])
    
    fprintf(1, ['\n\nglm for MI for Theta/' PACnames{pacii} '\n'])
    tbl = table(glm_mi_pre.data',glm_mi_pre.perCorr',glm_mi_pre.event',...
        'VariableNames',{'MI','perCorr','event'});
    mdl = fitglm(tbl,'MI~perCorr+event+perCorr*event'...
        ,'CategoricalVars',[2,3])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for average MI per mouse per odor pair for PAC theta' PACnames{pacii} ' prefrontal\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
    
    %Now plot the average PA variance for prefrontal for each electrode calculated per mouse
    %(including all sessions for each mouse)
    
    
    
    
    glm_PA_pre=[];
    glm_ii=0;
    
    id_ii=0;
    input_data=[];
    
    %Plot the average
    figNo = figNo +1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    %             try
    %                 close(figNo+pacii)
    %             catch
    %             end
    %             hFig=figure(figNo+pacii);
    
    set(hFig, 'units','normalized','position',[.1 .5 .3 .4])
    hold on
    
    bar_lab_loc=[];
    no_ev_labels=0;
    ii_gr_included=0;
    bar_offset = 0;
    
    edges=[0:20:1000];
    rand_offset=0.5;
    
    %Subtract the average for S+ proficient
    %Get these MI values
    ref_PA=[];
    ii_PA=0;
    evNo=1;
    per_ii=1;
    for ii=1:length(hippFileName)
        this_jj=[];
        for jj=1:all_pre(ii).handles_out.PA_ii
            if all_pre(ii).handles_out.PA_values(jj).pacii==pacii
                if all_pre(ii).handles_out.PA_values(jj).evNo==evNo
                    if all_pre(ii).handles_out.PA_values(jj).per_ii==per_ii
                        if all_pre(ii).handles_out.PA_values(jj).groupNo==grNo
                            this_jj=jj;
                        end
                        
                    end
                    
                end
            end
            
        end
        if ~isempty(this_jj)
            if mouse_op==1
                ref_PA(ii_PA+1:ii_PA+length(all_pre(ii).handles_out.PA_values(this_jj).PA_var_per_mouse))=all_pre(ii).handles_out.PA_values(this_jj).PA_var_per_mouse;
                ii_PA=ii_PA+length(all_pre(ii).handles_out.PA_values(this_jj).PA_var_per_mouse);
            else
                ii_PA=ii_PA+1;
                ref_PA(ii_PA)=all_pre(ii).handles_out.PA_values(this_jj).PA_var;
            end
        end
    end
    
    for evNo=1:2
        
        for per_ii=2:-1:1
            
            grNo=1;
            bar_offset = bar_offset +1;
            
            
            %Get these MI values
            these_PA=[];
            ii_PA=0;
            for ii=1:length(hippFileName)
                this_jj=[];
                for jj=1:all_pre(ii).handles_out.PA_ii
                    if all_pre(ii).handles_out.PA_values(jj).pacii==pacii
                        if all_pre(ii).handles_out.PA_values(jj).evNo==evNo
                            if all_pre(ii).handles_out.PA_values(jj).per_ii==per_ii
                                if all_pre(ii).handles_out.PA_values(jj).groupNo==grNo
                                    this_jj=jj;
                                end
                                
                            end
                            
                        end
                    end
                    
                end
                if ~isempty(this_jj)
                    if mouse_op==1
                        these_PA(ii_PA+1:ii_PA+length(all_pre(ii).handles_out.PA_values(this_jj).PA_var_per_mouse))=all_pre(ii).handles_out.PA_values(this_jj).PA_var_per_mouse;
                        ii_PA=ii_PA+length(all_pre(ii).handles_out.PA_values(this_jj).PA_var_per_mouse);
                    else
                        ii_PA=ii_PA+1;
                        these_PA(ii_PA)=all_pre(ii).handles_out.PA_values(this_jj).PA_var;
                    end
                end
            end
            
            these_PA=these_PA-mean(ref_PA);
            
            if evNo==2
                if per_ii==1
                    %S+ Proficient
                    bar(bar_offset,mean(these_PA),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                else
                    %S- Naive
                    bar(bar_offset,mean(these_PA),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                end
            else
                if per_ii==1
                    %S- Proficient
                    bar(bar_offset,mean(these_PA),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                else
                    %S- naive
                    bar(bar_offset,mean(these_PA),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                end
            end
            
            
            %Violin plot
            
            [mean_out, CIout]=drgViolinPoint(these_PA,edges,bar_offset,rand_offset,'k','k',3);
            %             CI = bootci(1000, {@mean, these_PA},'type','cper');
            %             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            %             plot(bar_offset*ones(1,length(these_PA)),these_PA,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            %                                 %Save data for glm and ranksum
            
            
            glm_PA_pre.data(glm_ii+1:glm_ii+length(these_PA))=these_PA;
            glm_PA_pre.perCorr(glm_ii+1:glm_ii+length(these_PA))=per_ii*ones(1,length(these_PA));
            glm_PA_pre.brain_region(glm_ii+1:glm_ii+length(these_mi))=2*ones(1,length(these_mi));
            glm_PA_pre.event(glm_ii+1:glm_ii+length(these_PA))=evNo*ones(1,length(these_PA));
            glm_ii=glm_ii+length(these_PA);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=these_PA;
            input_data(id_ii).description=[evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
            
            
            glm_PA_both.data(glm_PA_ii_both+1:glm_PA_ii_both+length(these_PA))=these_PA;
            glm_PA_both.brain_region(glm_PA_ii_both+1:glm_PA_ii_both+length(these_mi))=2*ones(1,length(these_mi));
            glm_PA_both.perCorr(glm_PA_ii_both+1:glm_PA_ii_both+length(these_PA))=per_ii*ones(1,length(these_PA));
            glm_PA_both.event(glm_PA_ii_both+1:glm_PA_ii_both+length(these_PA))=evNo*ones(1,length(these_PA));
            glm_PA_ii_both=glm_PA_ii_both+length(these_PA);
            
            
        end
        bar_offset = bar_offset + 1;
        
    end
    
    plot([0 6],[0 0],'-k','LineWidth',2)
    xlim([0 6])
    
    title(['Delta PA variance per mouse per odor pair for PAC theta/' PACnames{pacii} ' prefrontal'])
    
    
    
    
    %Proficient/Naive annotations
    annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color',[158/255 31/255 99/255],'LineStyle','none');
    annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color',[0 114/255 178/255],'LineStyle','none');
    
    
    xticks([1 2 4 5])
    xticklabels({'nS+', 'pS+','nS-', 'pS-'})
    
    yticks([-1000 0 1000])
    
    ylim([-1500 1000])
    
    ylabel('Delta PA variance')
    
    
    
    %Perform the glm
    fprintf(1, ['glm for delta PA variance per mouse per odor pair for PAC theta' PACnames{pacii} ' prefrontal\n'])
    
    fprintf(1, ['\n\nglm for delta PA variance for Theta/' PACnames{pacii} '\n'])
    tbl = table(glm_PA_pre.data',glm_PA_pre.perCorr',glm_PA_pre.event',...
        'VariableNames',{'MI','perCorr','event'});
    mdl = fitglm(tbl,'MI~perCorr+event+perCorr*event'...
        ,'CategoricalVars',[2,3])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for delta PA variance per mouse per odor pair for PAC theta' PACnames{pacii} ' prefrontal\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
    %Perform the glm mi for prefrontal
    fprintf(1, ['glm for average MI per mouse per odor pair for PAC theta' PACnames{pacii} ' prefrontal and hippocampus\n'])
    
    fprintf(1, ['\n\nglm for MI for Theta/' PACnames{pacii} '\n'])
    tbl = table(glm_mi_pre.data',glm_mi_pre.perCorr',glm_mi_pre.event',...
        'VariableNames',{'MI','perCorr','event'});
    mdl = fitglm(tbl,'MI~perCorr+event+perCorr*event'...
        ,'CategoricalVars',[2,3])
    
    %Perform the glm PA for both brain regions
    fprintf(1, ['glm for delta PA variance per mouse per odor pair for PAC theta' PACnames{pacii} ' prefrontal and hippocampus\n'])
    
    fprintf(1, ['\n\nglm for PA variance for Theta/' PACnames{pacii} '\n'])
    tbl = table(glm_PA_both.data',glm_PA_both.brain_region',glm_PA_both.perCorr',glm_PA_both.event',...
        'VariableNames',{'PA','brain_region','perCorr','event'});
    mdl = fitglm(tbl,'PA~perCorr+brain_region+event+perCorr*brain_region*event'...
        ,'CategoricalVars',[2,3,4])
    
end
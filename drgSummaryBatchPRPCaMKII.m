function drgSummaryBatchPRPCaMKII
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

peak_label{1}='Trough';
peak_label{2}='Peak';

% %Location of files
% % hippPathName='F:\Datos summary CaMKII111720\PRP drgAnalysisBatchLFPCaMKII case 24 output for summary\';
% hippPathName='/Users/restrepd/Documents/Projects/CaMKII analysis/PRP old/';
% %Hippocampus
% hippFileName{1}='spm_LFP_acetowavephasepower32620_hippoPRPnew.mat';
% hippFileName{2}='spm_LFP_acetowavephasepower32620_hippoPRPnewpropyl.mat';
% hippFileName{3}='spm_LFP_ethylbenwavephasepower41420_hippoPRPnew.mat';
% hippFileName{4}='spm_LFP_ethylwavephasepower3520_hippoPRPnew.mat';
% hippFileName{5}='spm_LFP_pz1ethyllwavephasepower0213020_hippoPRPnew.mat';
% hippFileName{6}='spm_LFP_pz1propylwavephasepower013020_hippoPRPnew.mat';
% hippFileName{7}='spm_LFP_pzz1ethyllwavephasepower043020_hippoPRPnew.mat';
% hippFileName{8}='spm_LFP_pzz1propylwavephasepower071220_hippoPRPnew.mat';
%
% % prePathName='F:\Datos summary CaMKII111720\PRP drgAnalysisBatchLFPCaMKII case 24 output for summary\';
% prePathName='/Users/restrepd/Documents/Projects/CaMKII analysis/PRP old/';
%
% %Prefrontal
% preFileName{1}='spm_LFP_acetowavephasepower32620_prefrontPRPnew.mat';
% preFileName{2}='spm_LFP_acetowavephasepower32620_prefrontPRPnewpropyl.mat';
% preFileName{3}='spm_LFP_ethylbenwavephasepower41420_prefrontPRPnew.mat';
% preFileName{4}='spm_LFP_ethylwavephasepower3520_prefrontPRPnew.mat';
% preFileName{5}='spm_LFP_pz1ethyllwavephasepower0213020_prefrontPRPnew.mat';
% preFileName{6}='spm_LFP_pz1propylwavephasepower013020_prefrontPRPnew.mat';
% preFileName{7}='spm_LFP_pzz1ethyllwavephasepower043020_prefrontPRPnew.mat';
% preFileName{8}='spm_LFP_pzz1propylwavephasepower071220_prefrontPRPnew.mat';
%
%

%Location of files
% hippPathName='F:\Datos summary CaMKII111720\PRP drgAnalysisBatchLFPCaMKII case 24 output for summary\';
hippPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/PRP drgAnalysisBatchLFPCaMKII case 24 output for 80/';
% fwd_rev_hippo=[

% %Hippocampus
% hippFileName{1}='spm_LFP_acetowavephasepower32620_hippoPRPnew.mat';
% hippFileName{2}='spm_LFP_acetowavephasepower32620_hippoPRPnewpropyl.mat';
% hippFileName{3}='spm_LFP_ethylbenwavephasepower41420_hippoPRPnew.mat';
% hippFileName{4}='spm_LFP_ethylwavephasepower3520_hippoPRPnew.mat';
% hippFileName{5}='spm_LFP_pz1ethyllwavephasepower0213020_hippoPRPnew.mat';
% hippFileName{6}='spm_LFP_pz1propylwavephasepower013020_hippoPRPnew.mat';
% hippFileName{7}='spm_LFP_pzz1ethyllwavephasepower043020_hippoPRPnew.mat';
% hippFileName{8}='spm_LFP_pzz1propylwavephasepower071220_hippoPRPnew.mat';

% %Hippocampus
% hippFileName{1}='spm_LFP_acetowavephasepower32620_hippocampusLFP2.mat';
% hippFileName{2}='spm_LFP_ethylbenwavephasepower41420_hippocampusLFP2.mat';
% hippFileName{3}='spm_LFP_ethylwavephasepower3520_hippocampusLFP2.mat';
% hippFileName{4}='spm_LFP_acetowavephasepower32620 2_hippocampusLFP2.mat';
% hippFileName{5}='spm_LFP_pz1ethyllwavephasepower0213020_hippocampusLFP2.mat';
% hippFileName{6}='spm_LFP_pz1propylwavephasepower013020_hippocampusLFP2.mat';
% hippFileName{7}='spm_LFP_pzz1ethyllwavephasepower043020_hippocampusLFP2.mat';
% hippFileName{8}='spm_LFP_pzz1propylwavephasepower071220_hippocampusLFP2.mat';

%Hippocampus proficient =[80 100]
hippFileName{1}='spm_LFP_acetowavephasepower32620_hippocampusLFP80.mat';
hippFileName{2}='spm_LFP_ethylbenwavephasepower41420_hippocampusLFP80.mat';
hippFileName{3}='spm_LFP_ethylwavephasepower3520_hippocampusLFP80.mat';
hippFileName{4}='spm_LFP_acetowavephasepower32620 2_hippocampusLFP80.mat';
hippFileName{5}='spm_LFP_pz1ethyllwavephasepower0213020_hippocampusLFP80.mat';
hippFileName{6}='spm_LFP_pz1propylwavephasepower013020_hippocampusLFP80.mat';
hippFileName{7}='spm_LFP_pzz1ethyllwavephasepower043020_hippocampusLFP80.mat';
hippFileName{8}='spm_LFP_pzz1propylwavephasepower071220_hippocampusLFP80.mat';

% prePathName='F:\Datos summary CaMKII111720\PRP drgAnalysisBatchLFPCaMKII case 24 output for summary\';
prePathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/PRP drgAnalysisBatchLFPCaMKII case 24 output for 80/';

% %Prefrontal
% preFileName{1}='spm_LFP_acetowavephasepower32620_prefrontPRPnew.mat';
% preFileName{2}='spm_LFP_acetowavephasepower32620_prefrontPRPnewpropyl.mat';
% preFileName{3}='spm_LFP_ethylbenwavephasepower41420_prefrontPRPnew.mat';
% preFileName{4}='spm_LFP_ethylwavephasepower3520_prefrontPRPnew.mat';
% preFileName{5}='spm_LFP_pz1ethyllwavephasepower0213020_prefrontPRPnew.mat';
% preFileName{6}='spm_LFP_pz1propylwavephasepower013020_prefrontPRPnew.mat';
% preFileName{7}='spm_LFP_pzz1ethyllwavephasepower043020_prefrontPRPnew.mat';
% preFileName{8}='spm_LFP_pzz1propylwavephasepower071220_prefrontPRPnew.mat';

% %Prefrontal
% preFileName{1}='spm_LFP_acetowavephasepower32620_prefrontalLFP2.mat';
% preFileName{2}='spm_LFP_ethylbenwavephasepower41420_prefrontalLFP2.mat';
% preFileName{3}='spm_LFP_ethylwavephasepower3520_prefrontalLFP2.mat';
% preFileName{4}='spm_LFP_acetowavephasepower32620 2_prefrontalLFP2.mat';
% preFileName{5}='spm_LFP_pz1ethyllwavephasepower0213020_prefrontalLFP2.mat';
% preFileName{6}='spm_LFP_pz1propylwavephasepower013020_prefrontalLFP2.mat';
% preFileName{7}='spm_LFP_pzz1ethyllwavephasepower043020_prefrontalLFP2.mat';
% preFileName{8}='spm_LFP_pzz1propylwavephasepower071220_prefrontalLFP2.mat';

%Prefrontal proficient =[80 100]
preFileName{1}='spm_LFP_acetowavephasepower32620_prefrontalLFP80.mat';
preFileName{2}='spm_LFP_ethylbenwavephasepower41420_prefrontalLFP80.mat';
preFileName{3}='spm_LFP_ethylwavephasepower3520_prefrontalLFP80.mat';
preFileName{4}='spm_LFP_acetowavephasepower32620 2_prefrontalLFP80.mat';
preFileName{5}='spm_LFP_pz1ethyllwavephasepower0213020_prefrontalLFP80.mat';
preFileName{6}='spm_LFP_pz1propylwavephasepower013020_prefrontalLFP80.mat';
preFileName{7}='spm_LFP_pzz1ethyllwavephasepower043020_prefrontalLFP80.mat';
preFileName{8}='spm_LFP_pzz1propylwavephasepower071220_prefrontalLFP80.mat';


%Now process the hippocampus

%Load data
all_hippo=[];

for ii=1:length(hippFileName)
    load([hippPathName hippFileName{ii}])
    all_hippo(ii).handles_out=handles_out;
end

figNo=0;

%Now plot the average PRP for each electrode calculated per mouse
%(including all sessions for each mouse)
edges=[-15:0.1:10];
rand_offset=0.5;



for peak=0:1
    for pacii=[1 3]    %for amplitude bandwidths (beta, low gamma, high gamma)
        
        glm_PRP=[];
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
        
        set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
        hold on
        
        bar_lab_loc=[];
        no_ev_labels=0;
        ii_gr_included=0;
        bar_offset = 0;
        
       
        
        for evNo=1:2
            
            for per_ii=2:-1:1
                
                for grNo=1:3
                    
                    bar_offset = bar_offset + 1;
                    
                    %Get these MI values
                    these_PRP=[];
                    ii_PRP=0;
                    for ii=1:length(hippFileName)
                        this_jj=[];
                        for jj=1:all_hippo(ii).handles_out.PRP_ii
                            %                             fprintf(1, ['peak = %d\n'],all_hippo(ii).handles_out.PRP_values(jj).peak)
                            if all_hippo(ii).handles_out.PRP_values(jj).pacii==pacii
                                if all_hippo(ii).handles_out.PRP_values(jj).evNo==evNo
                                    if all_hippo(ii).handles_out.PRP_values(jj).per_ii==per_ii
                                        if all_hippo(ii).handles_out.PRP_values(jj).peak==peak
                                            
                                            if all_hippo(ii).handles_out.PRP_values(jj).groupNo==grNo
                                                this_jj=jj;
                                            end
                                        end
                                    end
                                    
                                end
                            end
                        end
                        if ~isempty(this_jj)
                            if mouse_op==1
                                these_PRP(ii_PRP+1:ii_PRP+length(all_hippo(ii).handles_out.PRP_values(this_jj).PRP_per_mouse))=all_hippo(ii).handles_out.PRP_values(this_jj).PRP_per_mouse;
                                ii_PRP=ii_PRP+length(all_hippo(ii).handles_out.PRP_values(this_jj).PRP_per_mouse);
                            else
                                ii_PRP=ii_PRP+1;
                                these_PRP(ii_PRP)=all_hippo(ii).handles_out.PRP_values(this_jj).PRP;
                            end
                        end
                    end
                    
                    
                    switch grNo
                        case 1
                            bar(bar_offset,mean(these_PRP),'g','LineWidth', 3,'EdgeColor','none')
                        case 2
                            bar(bar_offset,mean(these_PRP),'b','LineWidth', 3,'EdgeColor','none')
                        case 3
                            bar(bar_offset,mean(these_PRP),'y','LineWidth', 3,'EdgeColor','none')
                    end
                    
                    
                    %Violin plot
                    
                    [mean_out, CIout]=drgViolinPoint(these_PRP,edges,bar_offset,rand_offset,'k','k',3);
                    
                    %                     CI = bootci(1000, {@mean, these_PRP},'type','cper');
                    %                     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                    %                     plot(bar_offset*ones(1,length(these_PRP)),these_PRP,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                    %
                    
                    %                                 %Save data for glm and ranksum
                    
                    glm_PRP.data(glm_ii+1:glm_ii+length(these_PRP))=these_PRP;
                    glm_PRP.group(glm_ii+1:glm_ii+length(these_PRP))=grNo*ones(1,length(these_PRP));
                    glm_PRP.perCorr(glm_ii+1:glm_ii+length(these_PRP))=per_ii*ones(1,length(these_PRP));
                    glm_PRP.event(glm_ii+1:glm_ii+length(these_PRP))=evNo*ones(1,length(these_PRP));
                    glm_ii=glm_ii+length(these_PRP);
                    
                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_PRP;
                    input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                    
                    
                end
                bar_offset = bar_offset + 1;
                
            end
            bar_offset = bar_offset + 1;
            
        end
        
        title(['Average PRP calculated per mouse per odor pair for ' peak_label{peak+1} ' theta/' PACnames{pacii} ' hippocampus'])
        
        
        xticks([1 2 3 5 6 7 10 11 12 14 15 16])
        xticklabels({'nwtS+', 'nHS+', 'nKOS+', 'pwtS+', 'pHS+', 'pKOS+', 'nwtS-', 'nHS-', 'nKOS-', 'pwtS-', 'pHS-', 'pKOS-'})
        
        ylim([-15 10])
        
        ylabel('PRP')
        
        
        %Perform the glm
        fprintf(1, ['glm for average PRP calculated per mouse per odor pair for ' peak_label{peak+1} ' theta' PACnames{pacii} ' hippocampus\n'])
        
        fprintf(1, ['\n\nglm for PRP for Theta/' PACnames{pacii} '\n'])
        tbl = table(glm_PRP.data',glm_PRP.group',glm_PRP.perCorr',glm_PRP.event',...
            'VariableNames',{'MI','group','naive_vs_proficient','sp_vs_sm'});
        mdl = fitglm(tbl,'MI~group+naive_vs_proficient+sp_vs_sm+naive_vs_proficient*group*sp_vs_sm'...
            ,'CategoricalVars',[2,3,4])
        
        
        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for average PRP calculated per mouse per odor pair for ' peak_label{peak+1} ' theta' PACnames{pacii} ' hippocampus\n'])
        [output_data] = drgMutiRanksumorTtest(input_data);
        
        
    end
end

%Now plot the average AUC calculated per mouse for hippocampus
%(including all sessions for each mouse)
edges=[0:100:3500];
rand_offset=0.5;

for peak=0:1
    for pacii=[1 3]    %for amplitude bandwidths (beta, low gamma, high gamma)
        
        glm_AUC=[];
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
        
        %             try
        %                 close(figNo+pacii)
        %             catch
        %             end
        %             hFig=figure(figNo+pacii);
        
        set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
        hold on
        
        bar_lab_loc=[];
        no_ev_labels=0;
        ii_gr_included=0;
        bar_offset = 0;
        
        %     for evNo=1:2
        
        for per_ii=2:-1:1
            
            for grNo=1:3
                bar_offset = bar_offset +1;
                
                %                         if sum(eventType==3)>0
                %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
                %                         else
                %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                %                         end
                %
                %                         these_offsets(per_ii)=bar_offset;
                bar_offset = bar_offset + 1;
                
                %Get these MI values
                these_AUC=[];
                ii_AUC=0;
                for ii=1:length(hippFileName)
                    this_jj=[];
                    for jj=1:all_hippo(ii).handles_out.AUC_ii
                        if all_hippo(ii).handles_out.AUC_values(jj).pacii==pacii
                            if all_hippo(ii).handles_out.AUC_values(jj).per_ii==per_ii
                                if all_hippo(ii).handles_out.AUC_values(jj).peak==peak
                                    if all_hippo(ii).handles_out.AUC_values(jj).groupNo==grNo
                                        this_jj=jj;
                                    end
                                end
                            end
                        end
                        
                    end
                    if ~isempty(this_jj)
                        ii_AUC=ii_AUC+1;
                        these_AUC(ii_AUC)=all_hippo(ii).handles_out.AUC_values(this_jj).AUC;
                    end
                end
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(these_AUC),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(these_AUC),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(these_AUC),'y','LineWidth', 3,'EdgeColor','none')
                end
                
                
                %Violin plot
                
                %                 [mean_out, CIout]=drgViolinPoint(these_AUC,edges,bar_offset,rand_offset,'k','k',1);
                CI = bootci(1000, {@mean, these_AUC},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(these_AUC)),these_AUC,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                %                                 %Save data for glm and ranksum
                
                glm_AUC.data(glm_ii+1:glm_ii+length(these_AUC))=these_AUC;
                glm_AUC.group(glm_ii+1:glm_ii+length(these_AUC))=grNo*ones(1,length(these_AUC));
                glm_AUC.perCorr(glm_ii+1:glm_ii+length(these_AUC))=per_ii*ones(1,length(these_AUC));
                glm_AUC.event(glm_ii+1:glm_ii+length(these_AUC))=evNo*ones(1,length(these_AUC));
                glm_ii=glm_ii+length(these_AUC);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_AUC;
                input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                
                
            end
            bar_offset = bar_offset + 2;
            
        end
        bar_offset = bar_offset + 3;
        
        %     end
        
        title(['Average PRP auROC for each electrode calculated per mouse for ' peak_label{peak+1} ' theta/' PACnames{pacii} ' hippocampus'])
        
        
        %Annotations identifying groups
        x_interval=0.8/ii_gr_included;
        for ii=1:ii_gr_included
            annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
        end
        
        %Proficient/Naive annotations
        annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
        annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
        
        
        xticks([2 4 6 10 12 14 21 23 25 29 31 33])
        xticklabels({'nwtS+', 'nHS+', 'nKOS+', 'pwtS+', 'pHS+', 'pKOS+', 'nwtS-', 'nHS-', 'nKOS-', 'pwtS-', 'pHS-', 'pKOS-'})
        
        ylim([0 0.5])
        ylabel('auROC')
        
        %Perform the glm
        fprintf(1, ['glm for average PRP auROC for each electrode calculated per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' hippocampus\n'])
        
        fprintf(1, ['\n\nglm for auROC for Theta/' PACnames{pacii} '\n'])
        tbl = table(glm_AUC.data',glm_AUC.group',glm_AUC.perCorr',glm_AUC.event',...
            'VariableNames',{'MI','group','perCorr','event'});
        mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
            ,'CategoricalVars',[2,3,4])
        
        
        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for average PRP auROC for each electrode calculated per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' hippocampus\n'])
        [output_data] = drgMutiRanksumorTtest(input_data);
        
        
    end
end


%Now process the prefrontal



%Load data
all_pre=[];

for ii=1:length(preFileName)
    load([prePathName preFileName{ii}])
    all_pre(ii).handles_out=handles_out;
end


%Now plot the average MI for each electrode calculated per mouse
%(including all sessions for each mouse)
edges=[-15:0.1:10];
rand_offset=0.5;

for peak=0:1
    for pacii=[1 3]    %for amplitude bandwidths (beta, low gamma, high gamma)
        
        glm_PRP=[];
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
        
        set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
        hold on
        
        bar_lab_loc=[];
        no_ev_labels=0;
        ii_gr_included=0;
        bar_offset = 0;
        
        
        for evNo=1:2
            
            for per_ii=2:-1:1
                
                for grNo=1:3
                    bar_offset = bar_offset +1;
                    
                    
                    %Get these MI values
                    these_PRP=[];
                    ii_PRP=0;
                    for ii=1:length(hippFileName)
                        this_jj=[];
                        for jj=1:all_pre(ii).handles_out.PRP_ii
                            if all_pre(ii).handles_out.PRP_values(jj).pacii==pacii
                                if all_pre(ii).handles_out.PRP_values(jj).evNo==evNo
                                    if all_pre(ii).handles_out.PRP_values(jj).per_ii==per_ii
                                        if all_pre(ii).handles_out.PRP_values(jj).peak==peak
                                            if all_pre(ii).handles_out.PRP_values(jj).groupNo==grNo
                                                this_jj=jj;
                                            end
                                        end
                                    end
                                    
                                end
                            end
                            
                        end
                        if ~isempty(this_jj)
                            if mouse_op==1
                                these_PRP(ii_PRP+1:ii_PRP+length(all_pre(ii).handles_out.PRP_values(this_jj).PRP_per_mouse))=all_pre(ii).handles_out.PRP_values(this_jj).PRP_per_mouse;
                                ii_PRP=ii_PRP+length(all_pre(ii).handles_out.PRP_values(this_jj).PRP_per_mouse);
                            else
                                ii_PRP=ii_PRP+1;
                                these_PRP(ii_PRP)=all_pre(ii).handles_out.PRP_values(this_jj).PRP;
                            end
                        end
                    end
                    
        
                    switch grNo
                        case 1
                            bar(bar_offset,mean(these_PRP),'g','LineWidth', 3,'EdgeColor','none')
                        case 2
                            bar(bar_offset,mean(these_PRP),'b','LineWidth', 3,'EdgeColor','none')
                        case 3
                            bar(bar_offset,mean(these_PRP),'y','LineWidth', 3,'EdgeColor','none')
                    end
                    
                    
                    %Violin plot
                    
                    [mean_out, CIout]=drgViolinPoint(these_PRP,edges,bar_offset,rand_offset,'k','k',3);
                    
                    %                     CI = bootci(1000, {@mean, these_PRP},'type','cper');
                    %                     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                    %                     plot(bar_offset*ones(1,length(these_PRP)),these_PRP,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                    %
                    
                    %                                 %Save data for glm and ranksum
                    
                    glm_PRP.data(glm_ii+1:glm_ii+length(these_PRP))=these_PRP;
                    glm_PRP.group(glm_ii+1:glm_ii+length(these_PRP))=grNo*ones(1,length(these_PRP));
                    glm_PRP.perCorr(glm_ii+1:glm_ii+length(these_PRP))=per_ii*ones(1,length(these_PRP));
                    glm_PRP.event(glm_ii+1:glm_ii+length(these_PRP))=evNo*ones(1,length(these_PRP));
                    glm_ii=glm_ii+length(these_PRP);
                    
                    id_ii=id_ii+1;
                    input_data(id_ii).data=these_PRP;
                    input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                    
                    
                end
                bar_offset = bar_offset + 1;
                
            end
            bar_offset = bar_offset + 1;
            
        end
        
        title(['Average PRP calculated per mouse per odor pair for ' peak_label{peak+1} ' theta/' PACnames{pacii} ' prefrontal'])
        
        
        
        
        xticks([1 2 3 5 6 7 10 11 12 14 15 16])
        xticklabels({'nwtS+', 'nHS+', 'nKOS+', 'pwtS+', 'pHS+', 'pKOS+', 'nwtS-', 'nHS-', 'nKOS-', 'pwtS-', 'pHS-', 'pKOS-'})
        
        ylim([-15 10])
        
        
        ylabel('PRP')
        ylim([-15 10])
        
        %Perform the glm
        fprintf(1, ['glm for average PRP calculated per mouse per odor pair for ' peak_label{peak+1} ' theta' PACnames{pacii} ' prefrontal\n'])
        
        fprintf(1, ['\n\nglm for PRP for Theta/' PACnames{pacii} '\n'])
        tbl = table(glm_PRP.data',glm_PRP.group',glm_PRP.perCorr',glm_PRP.event',...
            'VariableNames',{'MI','group','perCorr','event'});
        mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
            ,'CategoricalVars',[2,3,4])
        
        
        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for average PRP calculated per mouse per odor pair for ' peak_label{peak+1} ' PAC theta' PACnames{pacii} ' prefrontal\n'])
        [output_data] = drgMutiRanksumorTtest(input_data);
        
        
    end
end

%Now plot the average PA variance for each electrode calculated per mouse
%(including all sessions for each mouse)
edges=[0:100:3500];
rand_offset=0.5;

for peak=0:1
    for pacii=[1 3]    %for amplitude bandwidths (beta, low gamma, high gamma)
        
        glm_AUC=[];
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
        
        %             try
        %                 close(figNo+pacii)
        %             catch
        %             end
        %             hFig=figure(figNo+pacii);
        
        set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
        hold on
        
        bar_lab_loc=[];
        no_ev_labels=0;
        ii_gr_included=0;
        bar_offset = 0;
        
        %     for evNo=1:2
        
        for per_ii=2:-1:1
            
            for grNo=1:3
                bar_offset = bar_offset +1;
                
                %                         if sum(eventType==3)>0
                %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
                %                         else
                %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                %                         end
                %
                %                         these_offsets(per_ii)=bar_offset;
                bar_offset = bar_offset + 1;
                
                %Get these MI values
                these_AUC=[];
                ii_AUC=0;
                for ii=1:length(hippFileName)
                    this_jj=[];
                    for jj=1:all_pre(ii).handles_out.AUC_ii
                        if all_pre(ii).handles_out.AUC_values(jj).pacii==pacii
                            %                             if all_pre(ii).handles_out.AUC_values(jj).evNo==evNo
                            if all_pre(ii).handles_out.AUC_values(jj).per_ii==per_ii
                                if all_pre(ii).handles_out.AUC_values(jj).peak==peak
                                    if all_pre(ii).handles_out.AUC_values(jj).groupNo==grNo
                                        this_jj=jj;
                                    end
                                end
                            end
                            
                            %                             end
                        end
                        
                    end
                    if ~isempty(this_jj)
                        ii_AUC=ii_AUC+1;
                        these_AUC(ii_AUC)=all_pre(ii).handles_out.AUC_values(this_jj).AUC;
                    end
                end
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(these_AUC),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(these_AUC),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(these_AUC),'y','LineWidth', 3,'EdgeColor','none')
                end
                
                
                %Violin plot
                
                %                 [mean_out, CIout]=drgViolinPoint(these_AUC,edges,bar_offset,rand_offset,'k','k',1);
                CI = bootci(1000, {@mean, these_AUC},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(these_AUC)),these_AUC,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                %                                 %Save data for glm and ranksum
                
                glm_AUC.data(glm_ii+1:glm_ii+length(these_AUC))=these_AUC;
                glm_AUC.group(glm_ii+1:glm_ii+length(these_AUC))=grNo*ones(1,length(these_AUC));
                glm_AUC.perCorr(glm_ii+1:glm_ii+length(these_AUC))=per_ii*ones(1,length(these_AUC));
                glm_AUC.event(glm_ii+1:glm_ii+length(these_AUC))=evNo*ones(1,length(these_AUC));
                glm_ii=glm_ii+length(these_AUC);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_AUC;
                input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                
                
            end
            bar_offset = bar_offset + 2;
            
        end
        bar_offset = bar_offset + 3;
        
        %     end
        
        title(['Average PRP auROC for each electrode calculated per mouse for ' peak_label{peak+1} ' theta/' PACnames{pacii} ' prefrontal'])
        
        
        %Annotations identifying groups
        x_interval=0.8/ii_gr_included;
        for ii=1:ii_gr_included
            annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
        end
        
        %Proficient/Naive annotations
        annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
        annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
        
        
        xticks([2 4 6 10 12 14 21 23 25 29 31 33])
        xticklabels({'nwtS+', 'nHS+', 'nKOS+', 'pwtS+', 'pHS+', 'pKOS+', 'nwtS-', 'nHS-', 'nKOS-', 'pwtS-', 'pHS-', 'pKOS-'})
        
        ylim([0 0.5])
        
        ylabel('auROC')
        
        %Perform the glm
        fprintf(1, ['glm for average PRP auROC for each electrode calculated per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' prefrontal\n'])
        
        fprintf(1, ['\n\nglm for PRP auROC for Theta/' PACnames{pacii} '\n'])
        tbl = table(glm_AUC.data',glm_AUC.group',glm_AUC.perCorr',glm_AUC.event',...
            'VariableNames',{'MI','group','perCorr','event'});
        mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
            ,'CategoricalVars',[2,3,4])
        
        
        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for average PA variance for each electrode calculated per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' prefrontal\n'])
        [output_data] = drgMutiRanksumorTtest(input_data);
        
        
    end
end

%Get the hipp/pre delta PRPs
all_hipp_PRP=[];
all_pre_PRP=[];
all_slopes=[];

for peak=0:1
    for pacii=[1 3]    %for amplitude bandwidths (beta, low gamma, high gamma)
        
        
        for grNo=1:3
            all_hipp_PRP.peak(peak+1).pacii(pacii).group(grNo).PRP=[];
            all_pre_PRP.peak(peak+1).pacii(pacii).group(grNo).PRP=[];
        end
        
        for evNo=1:2
            
            for per_ii=2:-1:1

                for grNo=1:3
                    
                     all_slopes.peak(peak+1).pacii(pacii).group(grNo).slope=[];
                     all_slopes.peak(peak+1).pacii(pacii).group(grNo).rho=[];
                    
                    %Get these MI values
                    these_PRP_hipp=[];
                    ii_PRP_hipp=0;
                    these_PRP_pre=[];
                    ii_PRP_pre=0;
                    for ii=1:length(hippFileName)
                        this_jj_hipp=[];
                        for jj=1:all_hippo(ii).handles_out.PRP_ii
                            %                             fprintf(1, ['peak = %d\n'],all_hippo(ii).handles_out.PRP_values(jj).peak)
                            if all_hippo(ii).handles_out.PRP_values(jj).pacii==pacii
                                if all_hippo(ii).handles_out.PRP_values(jj).evNo==evNo
                                    if all_hippo(ii).handles_out.PRP_values(jj).per_ii==per_ii
                                        if all_hippo(ii).handles_out.PRP_values(jj).peak==peak
                                            
                                            if all_hippo(ii).handles_out.PRP_values(jj).groupNo==grNo
                                                this_jj_hipp=jj;
                                            end
                                        end
                                    end
                                    
                                end
                            end
                        end
                        
                        this_jj_pre=[];
                        for jj=1:all_pre(ii).handles_out.PRP_ii
                            if all_pre(ii).handles_out.PRP_values(jj).pacii==pacii
                                if all_pre(ii).handles_out.PRP_values(jj).evNo==evNo
                                    if all_pre(ii).handles_out.PRP_values(jj).per_ii==per_ii
                                        if all_pre(ii).handles_out.PRP_values(jj).peak==peak
                                            if all_pre(ii).handles_out.PRP_values(jj).groupNo==grNo
                                                this_jj_pre=jj;
                                            end
                                        end
                                    end
                                    
                                end
                            end
                            
                        end
                        
                        if (~isempty(this_jj_hipp))&(~isempty(this_jj_pre))
                            if length(all_hippo(ii).handles_out.PRP_values(this_jj_hipp).PRP_per_mouse)==length(all_pre(ii).handles_out.PRP_values(this_jj_pre).PRP_per_mouse)
                                these_PRP_hipp(ii_PRP_hipp+1:ii_PRP_hipp+length(all_hippo(ii).handles_out.PRP_values(this_jj_hipp).PRP_per_mouse))=all_hippo(ii).handles_out.PRP_values(this_jj_hipp).PRP_per_mouse;
                                ii_PRP_hipp=ii_PRP_hipp+length(all_hippo(ii).handles_out.PRP_values(this_jj_hipp).PRP_per_mouse);
                                these_PRP_pre(ii_PRP_pre+1:ii_PRP_pre+length(all_pre(ii).handles_out.PRP_values(this_jj_pre).PRP_per_mouse))=all_pre(ii).handles_out.PRP_values(this_jj_pre).PRP_per_mouse;
                                ii_PRP_pre=ii_PRP_pre+length(all_pre(ii).handles_out.PRP_values(this_jj_pre).PRP_per_mouse);
                            else
                                pffft=1;
                            end
                        end
                        
                        x=these_PRP_hipp';
                        y=these_PRP_pre';
                        
                        %Fit a line
                        c = polyfit(x,y,1);
                        
                        all_slopes.peak(peak+1).pacii(pacii).group(grNo).slope = [all_slopes.peak(peak+1).pacii(pacii).group(grNo).slope c(1)];
                        
                        [rho,pval]=corr(x,y);
                        
                        all_slopes.peak(peak+1).pacii(pacii).group(grNo).rho = [all_slopes.peak(peak+1).pacii(pacii).group(grNo).rho rho];
                        
                    end
                    
                    all_hipp_PRP.peak(peak+1).pacii(pacii).group(grNo).PRP=[all_hipp_PRP.peak(peak+1).pacii(pacii).group(grNo).PRP these_PRP_hipp];
                    all_pre_PRP.peak(peak+1).pacii(pacii).group(grNo).PRP=[all_pre_PRP.peak(peak+1).pacii(pacii).group(grNo).PRP these_PRP_pre];
                    
                    

                end
                
                
            end
            
            
        end
        
        
        
    end
end

%Plot the correlation between delta PRP of hippocampus vs prefrontal
for peak=0:1
    for pacii=[1 3]    %for amplitude bandwidths (beta, low gamma, high gamma)
        
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
        
        for grNo=1:3
            
            x=all_hipp_PRP.peak(peak+1).pacii(pacii).group(grNo).PRP';
            y=all_pre_PRP.peak(peak+1).pacii(pacii).group(grNo).PRP';
            [rho,pval]=corr(x,y);
            
            %Fit a line
            c = polyfit(x,y,1);
            % Display evaluated equation y = m*x + b
            disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
            % Evaluate fit equation using polyval
            y_est = polyval(c,[-15 10]);
           
            
            fprintf(1, ['rho = %d, p value = %d for ' peak_label{peak+1} ' theta' PACnames{pacii} ' ' group_legend{grNo} '\n\n'],rho,pval)
   
            switch grNo
                case 1
                    plot(all_hipp_PRP.peak(peak+1).pacii(pacii).group(grNo).PRP,all_pre_PRP.peak(peak+1).pacii(pacii).group(grNo).PRP,'og')
                     % Add trend line to plot
                    plot([-15 10],y_est,'g-','LineWidth',2)
                case 2
                    plot(all_hipp_PRP.peak(peak+1).pacii(pacii).group(grNo).PRP,all_pre_PRP.peak(peak+1).pacii(pacii).group(grNo).PRP,'ob')
                     % Add trend line to plot
                    plot([-15 10],y_est,'b-','LineWidth',2)
                case 3
                    plot(all_hipp_PRP.peak(peak+1).pacii(pacii).group(grNo).PRP,all_pre_PRP.peak(peak+1).pacii(pacii).group(grNo).PRP,'oy')
                     % Add trend line to plot
                    plot([-15 10],y_est,'y-','LineWidth',2)
            end
            
          
        end
        
        title(['delta PRP calculated per mouse per odor pair for ' peak_label{peak+1} ' theta/' PACnames{pacii}])
        
        
        
        xlim([-15 10])
        ylim([-15 10])
        
        ylabel('Prefrontal delta PRP')
        ylabel('Hippocampal delta PRP')
    end
end


%Plot the slope of the relationship between delta PRP of hippocampus vs prefrontal
for peak=0:1
    for pacii=[1 3]    %for amplitude bandwidths (beta, low gamma, high gamma)
        
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
        bar_offset=0;
        
        for grNo=1:3
            
            these_slopes=all_slopes.peak(peak+1).pacii(pacii).group(grNo).slope;
            
            switch grNo
                case 1
                    bar(bar_offset,mean(these_slopes),'g','LineWidth', 3,'EdgeColor','none')
                case 2
                    bar(bar_offset,mean(these_slopes),'b','LineWidth', 3,'EdgeColor','none')
                case 3
                    bar(bar_offset,mean(these_slopes),'y','LineWidth', 3,'EdgeColor','none')
            end
            
            CI = bootci(1000, {@mean, these_slopes'},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_slopes)),these_slopes,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            bar_offset=bar_offset+1;
        end
        
        title(['Slope for ' peak_label{peak+1} ' theta/' PACnames{pacii}])
        
        ylabel('Slope')
        
    end
end
% 
% 
% %Plot the rho of the relationship between delta PRP of hippocampus vs prefrontal
% for peak=0:1
%     for pacii=[1 3]    %for amplitude bandwidths (beta, low gamma, high gamma)
%         
%         %Plot the average
%         figNo = figNo +1;
%         
%         try
%             close(figNo)
%         catch
%         end
%         hFig=figure(figNo);
%         
%         %             try
%         %                 close(figNo+pacii)
%         %             catch
%         %             end
%         %             hFig=figure(figNo+pacii);
%         
%         set(hFig, 'units','normalized','position',[.1 .5 .3 .4])
%         hold on
%         bar_offset=0;
%         
%         for grNo=1:3
%             
%             these_rhos=all_slopes.peak(peak+1).pacii(pacii).group(grNo).rho;
%             
%             switch grNo
%                 case 1
%                     bar(bar_offset,mean(these_rhos),'g','LineWidth', 3,'EdgeColor','none')
%                 case 2
%                     bar(bar_offset,mean(these_rhos),'b','LineWidth', 3,'EdgeColor','none')
%                 case 3
%                     bar(bar_offset,mean(these_rhos),'y','LineWidth', 3,'EdgeColor','none')
%             end
%             
%             CI = bootci(1000, {@mean, these_rhos'},'type','cper');
%             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%             plot(bar_offset*ones(1,length(these_rhos)),these_rhos,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%             
%             bar_offset=bar_offset+1;
%         end
%         
%         title(['Rho for ' peak_label{peak+1} ' theta/' PACnames{pacii}])
%         
%         ylabel('Slope')
%         
%     end
% end
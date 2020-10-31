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

%Location of files
hippPathName='/Users/restrepd/OneDrive - The University of Colorado Denver/CaMKII Paper/PRP Summary/Hippo/';

%Hippocampus
hippFileName{1}='spm_LFP_acetowavephasepower32620_hippoPRPnew.mat';
hippFileName{2}='spm_LFP_acetowavephasepower32620_hippoPRPnewpropyl.mat';
hippFileName{3}='spm_LFP_ethylbenwavephasepower41420_hippoPRPnew.mat';

prePathName='/Users/restrepd/OneDrive - The University of Colorado Denver/CaMKII Paper/PRP Summary/prefront/';

%Prefrontal
preFileName{1}='spm_LFP_acetowavephasepower32620_prefrontPRPnew.mat';
preFileName{2}='spm_LFP_acetowavephasepower32620_prefrontPRPnewpropyl.mat';
preFileName{3}='spm_LFP_ethylbenwavephasepower41420_prefrontPRPnew.mat';


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
edges=[0:0.001:0.02];
rand_offset=0.8;

for peak=0:1
    for pacii=1:3    %for amplitude bandwidths (beta, low gamma, high gamma)
        
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
        
        for evNo=1:2
            
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
                            ii_PRP=ii_PRP+1;
                            these_PRP(ii_PRP)=all_hippo(ii).handles_out.PRP_values(this_jj).PRP;
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
                     
                    %[mean_out, CIout]=drgViolinPoint(these_PRP,edges,bar_offset,rand_offset,'k','k',1);
                    CI = bootci(1000, {@mean, these_PRP},'type','cper');
                    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                    plot(bar_offset*ones(1,length(these_PRP)),these_PRP,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                    
                    
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
                bar_offset = bar_offset + 2;
                
            end
            bar_offset = bar_offset + 3;
            
        end
        
        title(['Average PRP for each electrode calculated per mouse for ' peak_label{peak+1} ' theta/' PACnames{pacii} ' hippocampus'])
        
        
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
        
        ylim([-15 10])
        
        ylabel('PRP')
        
        
        %Perform the glm
        fprintf(1, ['glm for average PRP for each electrode calculated per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' hippocampus\n'])
        
        fprintf(1, ['\n\nglm for PRP for Theta/' PACnames{pacii} '\n'])
        tbl = table(glm_PRP.data',glm_PRP.group',glm_PRP.perCorr',glm_PRP.event',...
            'VariableNames',{'MI','group','perCorr','event'});
        mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
            ,'CategoricalVars',[2,3,4])
        
        
        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for average PRP for each electrode calculated per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' hippocampus\n'])
        [output_data] = drgMutiRanksumorTtest(input_data);
        
        
    end
end
%Now plot the average AUC variance for each electrode calculated per mouse
%(including all sessions for each mouse)
edges=[0:100:3500];
rand_offset=0.8;

for peak=0:1
    for pacii=1:3    %for amplitude bandwidths (beta, low gamma, high gamma)
        
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
edges=[0:0.001:0.02];
rand_offset=0.8;

for peak=0:1
    for pacii=1:3    %for amplitude bandwidths (beta, low gamma, high gamma)
        
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
        
        for evNo=1:2
            
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
                            ii_PRP=ii_PRP+1;
                            these_PRP(ii_PRP)=all_pre(ii).handles_out.PRP_values(this_jj).PRP;
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
                    
                    %[mean_out, CIout]=drgViolinPoint(these_PRP,edges,bar_offset,rand_offset,'k','k',1);
                    CI = bootci(1000, {@mean, these_PRP},'type','cper');
                    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                    plot(bar_offset*ones(1,length(these_PRP)),these_PRP,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                    
                    
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
                bar_offset = bar_offset + 2;
                
            end
            bar_offset = bar_offset + 3;
            
        end
        
        title(['Average PRP for each electrode calculated per mouse for ' peak_label{peak+1} ' theta/' PACnames{pacii} ' prefrontal'])
        
        
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
        
        
        ylabel('PRP')
        ylim([-15 10])
        
        %Perform the glm
        fprintf(1, ['glm for average PRP for each electrode calculated per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' prefrontal\n'])
        
        fprintf(1, ['\n\nglm for PRP for Theta/' PACnames{pacii} '\n'])
        tbl = table(glm_PRP.data',glm_PRP.group',glm_PRP.perCorr',glm_PRP.event',...
            'VariableNames',{'MI','group','perCorr','event'});
        mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
            ,'CategoricalVars',[2,3,4])
        
        
        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for average PRP for each electrode calculated per mouse for ' peak_label{peak+1} ' PAC theta' PACnames{pacii} ' prefrontal\n'])
        [output_data] = drgMutiRanksumorTtest(input_data);
        
        
    end
end

%Now plot the average PA variance for each electrode calculated per mouse
%(including all sessions for each mouse)
edges=[0:100:3500];
rand_offset=0.8;

for peak=0:1
    for pacii=1:3    %for amplitude bandwidths (beta, low gamma, high gamma)
        
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
pffft=1;
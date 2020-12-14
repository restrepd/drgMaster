function drgSummaryBatchPACCaMKII
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

%Location of files
PathName='D:\CaMKIIpaper\datos sumarry\';

%Hippocampus
hippFileName{1}='CaMKIIacetoLFPPACall121319_hippoPAC6figures.mat';
hippFileName{2}='CaMKIIethylbenLFPPACal222019_hippoPAC6figuresonly.mat';
hippFileName{3}='CaMKIIethylaceLFPPACall12219_hippoPAC6figuresonly.mat';
hippFileName{4}='CaMKIIpropylaceLFPPACall1220_hippoPAConlyfigures.mat';
hippFileName{5}='CaMKIIpointzero1ethylacepropylaceLFPPACall122019_hippoPAConlyfigures.mat';
hippFileName{6}='CaMKIIpointzero1propylaceLFPPACall11620_hippoPAConlyfigures.mat';
hippFileName{7}='CaMKIIpzz1ethylaceLFPPACall121119_hippoPAConlyfigures.mat';
hippFileName{8}='CaMKIIpzz1propylaceLFPPACall121919_hippoPAConly6figures.mat';

%Prefrontal
preFileName{1}='CaMKIIacetoLFPPACall121319_prefrontPAC6figures.mat';
preFileName{2}='CaMKIIethylbenLFPPACal222019_prefrontPAC6figuresonly.mat';
preFileName{3}='CaMKIIethylaceLFPPACall12219_prefrontPAC6figuresonly.mat';
preFileName{4}='CaMKIIpropylaceLFPPACall1220_prfrontPAConlyfigures.mat';
preFileName{5}='CaMKIIpointzero1ethylacepropylaceLFPPACall122019_prefrontPAConlyfigures.mat';
preFileName{6}='CaMKIIpointzero1propylaceLFPPACall11620_prefrontPAConlyfigures.mat';
preFileName{7}='CaMKIIpzz1ethylaceLFPPACall121119_prefrontPAConlyfigures.mat';
preFileName{8}='CaMKIIpzz1propylaceLFPPACall121919_prefrontPAConly6figures.mat';



%Now process the hippocampus

%Load data
all_hippo=[];

for ii=1:length(hippFileName)
   load([PathName hippFileName{ii}]) 
   all_hippo(ii).handles_out=handles_out;
end

figNo=0;

%Now plot the average MI for each electrode calculated per mouse
%(including all sessions for each mouse)
edges=[0:0.001:0.02];
rand_offset=0.8;

for pacii=1:3    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_mi=[];
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
                        ii_mi=ii_mi+1;
                        these_mi(ii_mi)=all_hippo(ii).handles_out.mi_values(this_jj).MI;
                    end
                end
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(these_mi),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(these_mi),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(these_mi),'y','LineWidth', 3,'EdgeColor','none')
                end
                

                %Violin plot
                
                %[mean_out, CIout]=drgViolinPoint(these_mi,edges,bar_offset,rand_offset,'k','k',1);
                CI = bootci(1000, {@mean, these_mi},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(these_mi)),these_mi,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                
                %                                 %Save data for glm and ranksum
                
                glm_mi.data(glm_ii+1:glm_ii+length(these_mi))=these_mi;
                glm_mi.group(glm_ii+1:glm_ii+length(these_mi))=grNo*ones(1,length(these_mi));
                glm_mi.perCorr(glm_ii+1:glm_ii+length(these_mi))=per_ii*ones(1,length(these_mi));
                glm_mi.event(glm_ii+1:glm_ii+length(these_mi))=evNo*ones(1,length(these_mi));
                glm_ii=glm_ii+length(these_mi);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_mi;
                input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                

            end
            bar_offset = bar_offset + 2;

        end
        bar_offset = bar_offset + 3;
  
    end
    
    title(['Average MI for each electrode calculated per mouse for PAC theta/' PACnames{pacii} ' hippocampus'])
    
   
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
    
    
    ylabel('Modulation Index')
    
    %Perform the glm
    fprintf(1, ['glm for average MI for each electrode calculated per mouse for PAC theta' PACnames{pacii} ' hippocampus\n'])
    
    fprintf(1, ['\n\nglm for MI for Theta/' PACnames{pacii} '\n'])
    tbl = table(glm_mi.data',glm_mi.group',glm_mi.perCorr',glm_mi.event',...
        'VariableNames',{'MI','group','perCorr','event'});
    mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
        ,'CategoricalVars',[2,3,4])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for average MI for each electrode calculated per mouse for PAC theta' PACnames{pacii} ' hippocampus\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end

%Now plot the average PA variance for each electrode calculated per mouse
%(including all sessions for each mouse)
edges=[0:100:3500];
rand_offset=0.8;

for pacii=1:3    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_PA=[];
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
                        ii_PA=ii_PA+1;
                        these_PA(ii_PA)=all_hippo(ii).handles_out.PA_values(this_jj).MI;
                    end
                end
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(these_PA),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(these_PA),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(these_PA),'y','LineWidth', 3,'EdgeColor','none')
                end
                

                %Violin plot
                
                %                 [mean_out, CIout]=drgViolinPoint(these_PA,edges,bar_offset,rand_offset,'k','k',1);
                CI = bootci(1000, {@mean, these_PA},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(these_PA)),these_PA,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                %                                 %Save data for glm and ranksum
                
                glm_PA.data(glm_ii+1:glm_ii+length(these_PA))=these_PA;
                glm_PA.group(glm_ii+1:glm_ii+length(these_PA))=grNo*ones(1,length(these_PA));
                glm_PA.perCorr(glm_ii+1:glm_ii+length(these_PA))=per_ii*ones(1,length(these_PA));
                glm_PA.event(glm_ii+1:glm_ii+length(these_PA))=evNo*ones(1,length(these_PA));
                glm_ii=glm_ii+length(these_PA);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_PA;
                input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                

            end
            bar_offset = bar_offset + 2;

        end
        bar_offset = bar_offset + 3;
  
    end
    
    title(['Average PA variance for each electrode calculated per mouse for PAC theta/' PACnames{pacii} ' hippocampus'])
    
   
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
    
    
    ylabel('Modulation Index')
    
    %Perform the glm
    fprintf(1, ['glm for average PA variance for each electrode calculated per mouse for PAC theta' PACnames{pacii} ' hippocampus\n'])
    
    fprintf(1, ['\n\nglm for PA variance for Theta/' PACnames{pacii} '\n'])
    tbl = table(glm_PA.data',glm_PA.group',glm_PA.perCorr',glm_PA.event',...
        'VariableNames',{'MI','group','perCorr','event'});
    mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
        ,'CategoricalVars',[2,3,4])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for average PA variance for each electrode calculated per mouse for PAC theta' PACnames{pacii} ' hippocampus\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end



%Now process the prefrontal

%Load data
all_pre=[];

for ii=1:length(preFileName)
   load([PathName preFileName{ii}]) 
   all_pre(ii).handles_out=handles_out;
end


%Now plot the average MI for each electrode calculated per mouse
%(including all sessions for each mouse)
edges=[0:0.001:0.02];
rand_offset=0.8;

for pacii=1:3    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_mi=[];
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
                        ii_mi=ii_mi+1;
                        these_mi(ii_mi)=all_pre(ii).handles_out.mi_values(this_jj).MI;
                    end
                end
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(these_mi),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(these_mi),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(these_mi),'y','LineWidth', 3,'EdgeColor','none')
                end
                

                %Violin plot
                
                %[mean_out, CIout]=drgViolinPoint(these_mi,edges,bar_offset,rand_offset,'k','k',1);
                CI = bootci(1000, {@mean, these_mi},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(these_mi)),these_mi,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                
                %                                 %Save data for glm and ranksum
                
                glm_mi.data(glm_ii+1:glm_ii+length(these_mi))=these_mi;
                glm_mi.group(glm_ii+1:glm_ii+length(these_mi))=grNo*ones(1,length(these_mi));
                glm_mi.perCorr(glm_ii+1:glm_ii+length(these_mi))=per_ii*ones(1,length(these_mi));
                glm_mi.event(glm_ii+1:glm_ii+length(these_mi))=evNo*ones(1,length(these_mi));
                glm_ii=glm_ii+length(these_mi);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_mi;
                input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                

            end
            bar_offset = bar_offset + 2;

        end
        bar_offset = bar_offset + 3;
  
    end
    
    title(['Average MI for each electrode calculated per mouse for PAC theta/' PACnames{pacii} ' prefrontal'])
    
   
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
    
    
    ylabel('Modulation Index')
    
    %Perform the glm
    fprintf(1, ['glm for average MI for each electrode calculated per mouse for PAC theta' PACnames{pacii} ' prefrontal\n'])
    
    fprintf(1, ['\n\nglm for MI for Theta/' PACnames{pacii} '\n'])
    tbl = table(glm_mi.data',glm_mi.group',glm_mi.perCorr',glm_mi.event',...
        'VariableNames',{'MI','group','perCorr','event'});
    mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
        ,'CategoricalVars',[2,3,4])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for average MI for each electrode calculated per mouse for PAC theta' PACnames{pacii} ' prefrontal\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end

%Now plot the average PA variance for each electrode calculated per mouse
%(including all sessions for each mouse)
edges=[0:100:3500];
rand_offset=0.8;

for pacii=1:3    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_PA=[];
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
                        ii_PA=ii_PA+1;
                        these_PA(ii_PA)=all_pre(ii).handles_out.PA_values(this_jj).MI;
                    end
                end
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(these_PA),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(these_PA),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(these_PA),'y','LineWidth', 3,'EdgeColor','none')
                end
                

                %Violin plot
                
                %                 [mean_out, CIout]=drgViolinPoint(these_PA,edges,bar_offset,rand_offset,'k','k',1);
                CI = bootci(1000, {@mean, these_PA},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(these_PA)),these_PA,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                %                                 %Save data for glm and ranksum
                
                glm_PA.data(glm_ii+1:glm_ii+length(these_PA))=these_PA;
                glm_PA.group(glm_ii+1:glm_ii+length(these_PA))=grNo*ones(1,length(these_PA));
                glm_PA.perCorr(glm_ii+1:glm_ii+length(these_PA))=per_ii*ones(1,length(these_PA));
                glm_PA.event(glm_ii+1:glm_ii+length(these_PA))=evNo*ones(1,length(these_PA));
                glm_ii=glm_ii+length(these_PA);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_PA;
                input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                

            end
            bar_offset = bar_offset + 2;

        end
        bar_offset = bar_offset + 3;
  
    end
    
    title(['Average PA variance for each electrode calculated per mouse for PAC theta/' PACnames{pacii} ' prefrontal'])
    
   
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
    
    
    ylabel('Modulation Index')
    
    %Perform the glm
    fprintf(1, ['glm for average PA variance for each electrode calculated per mouse for PAC theta' PACnames{pacii} ' prefrontal\n'])
    
    fprintf(1, ['\n\nglm for PA variance for Theta/' PACnames{pacii} '\n'])
    tbl = table(glm_PA.data',glm_PA.group',glm_PA.perCorr',glm_PA.event',...
        'VariableNames',{'MI','group','perCorr','event'});
    mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
        ,'CategoricalVars',[2,3,4])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for average PA variance for each electrode calculated per mouse for PAC theta' PACnames{pacii} ' prefrontal\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end
pffft=1;
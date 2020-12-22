function drgSummaryBatchCohCaMKIIWT
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'drgbDiscPar' file listing 'Discriminant_*.mat' output files
%from drgLFPDiscriminantBatch
%
%Performs summary analises for LDA and  PAC analysis performed with case 19
%of drgAnalysisBatchLFP

warning('off')

close all
clear all


bandwidth_names{1}='Theta';
bandwidth_names{2}='Beta';
bandwidth_names{3}='Low gamma';
bandwidth_names{4}='High gamma';

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
% hippPathName='E:\CaMKIIpaper\datos sumarry\coherence\';
hippPathName='F:\Datos summary CaMKII111720\Coherence summary\';

%Files
FileName{1}='CaMKIIacetocoh122119_coherenew.mat';
FileName{2}='CaMKIIethylbencohe1292019_coherenew.mat';
FileName{3}='CaMKIIpropylacecohe12302019_coherenew.mat';
FileName{4}='CaMKIIethylacepropylacecohe1272019_coherenew.mat';
FileName{5}='CaMKIIpointzero1ethylacepropylacecohe12282019_coherenew.mat';
FileName{6}='CaMKIIpzpropylacecohe1132020_coherenew.mat';
FileName{7}='CaMKIIpzz1ethylacepropylacecohe10212019_coherenew.mat';
FileName{8}='CaMKIIpzz1propylacecohe12262019_coherenew.mat';


%Load data
all_files=[];

for ii=1:length(FileName)
    load([hippPathName FileName{ii}])
    all_files(ii).handles_out=handles_out;
end

figNo=0;

%Now plot the average PRP for each electrode calculated per mouse
%(including all sessions for each mouse)
edges=[0:0.001:0.02];
rand_offset=0.8;


for bwii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_coh=[];
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
    
    set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
    hold on
    
    bar_lab_loc=[];
    no_ev_labels=0;
    ii_gr_included=0;
    bar_offset = 0;
    
    for evNo=1:2
        
        for per_ii=2:-1:1
            
            grNo=1;
            
            bar_offset = bar_offset + 1;
            
            %Get these MI values
            these_coh=[];
            ii_coh=0;
            for ii=1:length(FileName)
                this_jj=[];
                for jj=1:all_files(ii).handles_out.dcoh_ii
                    
                    if all_files(ii).handles_out.dcoh_values(jj).pacii==bwii
                        if all_files(ii).handles_out.dcoh_values(jj).evNo==evNo
                            if all_files(ii).handles_out.dcoh_values(jj).per_ii==per_ii
                                
                                
                                if all_files(ii).handles_out.dcoh_values(jj).groupNo==grNo
                                    this_jj=jj;
                                end
                                
                            end
                            
                        end
                    end
                end
                if ~isempty(this_jj)
                    ii_coh=ii_coh+1;
                    these_coh(ii_coh)=all_files(ii).handles_out.dcoh_values(this_jj).dcoh;
                end
            end
            
            
            if evNo==2
                if per_ii==1
                    %S+ Proficient
                    bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                else
                    %S- Naive
                    bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                end
            else
                if per_ii==1
                    %S- Proficient
                    bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                else
                    %S- naive
                    bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                end
            end
            
            
            %Violin plot
            
            %[mean_out, CIout]=drgViolinPoint(these_coh,edges,bar_offset,rand_offset,'k','k',1);
            CI = bootci(1000, {@mean, these_coh},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_coh)),these_coh,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            
            
            glm_coh.data(glm_ii+1:glm_ii+length(these_coh))=these_coh;
            %                 glm_coh.group(glm_ii+1:glm_ii+length(these_coh))=grNo*ones(1,length(these_coh));
            glm_coh.perCorr(glm_ii+1:glm_ii+length(these_coh))=per_ii*ones(1,length(these_coh));
            glm_coh.event(glm_ii+1:glm_ii+length(these_coh))=evNo*ones(1,length(these_coh));
            glm_ii=glm_ii+length(these_coh);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=these_coh;
            input_data(id_ii).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
            
        end
        bar_offset = bar_offset + 1;
        
    end
    
    title(['Average coherence for each electrode calculated per mouse for ' bandwidth_names{bwii}])
    
    
    %     %Annotations identifying groups
    %     x_interval=0.8/ii_gr_included;
    %     for ii=1:ii_gr_included
    %         annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
    %     end
    %
    %     %Proficient/Naive annotations
    %     annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
    %     annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
    
    
    xticks([1 2 4 5])
    xticklabels({'nS+', 'pS+','nS-', 'pS-'})
    
    ylabel('PRP')
    
    
    %Perform the glm
    fprintf(1, ['glm for average coherence for each electrode calculated per mouse for '  bandwidth_names{bwii} '\n'])
    
    fprintf(1, ['\n\nglm for PRP for' bandwidth_names{bwii} '\n'])
    tbl = table(glm_coh.data',glm_coh.perCorr',glm_coh.event',...
        'VariableNames',{'MI','perCorr','event'});
    mdl = fitglm(tbl,'MI~perCorr+event+perCorr*event'...
        ,'CategoricalVars',[2,3])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for average PRP for each electrode calculated per mouse for ' bandwidth_names{bwii} ' hippocampus\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end



% %Now plot the average AUC variance for each electrode calculated per mouse
% %(including all sessions for each mouse)
% edges=[0:100:3500];
% rand_offset=0.8;
%
%
% for bwii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
%
%     glm_AUC=[];
%     glm_ii=0;
%
%     id_ii=0;
%     input_data=[];
%
%     %Plot the average
%     figNo = figNo +1;
%
%     try
%         close(figNo)
%     catch
%     end
%     hFig=figure(figNo);
%
%     %             try
%     %                 close(figNo+bwii)
%     %             catch
%     %             end
%     %             hFig=figure(figNo+bwii);
%
%     set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
%     hold on
%
%     bar_lab_loc=[];
%     no_ev_labels=0;
%     ii_gr_included=0;
%     bar_offset = 0;
%
%     %     for evNo=1:2
%
%     for per_ii=2:-1:1
%
%         for grNo=1:3
%             bar_offset = bar_offset +1;
%
%             %                         if sum(eventType==3)>0
%             %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
%             %                         else
%             %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
%             %                         end
%             %
%             %                         these_offsets(per_ii)=bar_offset;
%             bar_offset = bar_offset + 1;
%
%             %Get these MI values
%             these_AUC=[];
%             ii_AUC=0;
%             for ii=1:length(FileName)
%                 this_jj=[];
%                 for jj=1:all_files(ii).handles_out.auc_ii
%                     if all_files(ii).handles_out.auc_values(jj).pacii==bwii
%                         if all_files(ii).handles_out.auc_values(jj).per_ii==per_ii
%
%                             if all_files(ii).handles_out.auc_values(jj).groupNo==grNo
%                                 this_jj=jj;
%                             end
%
%                         end
%                     end
%
%                 end
%                 if ~isempty(this_jj)
%                     ii_AUC=ii_AUC+1;
%                     these_AUC(ii_AUC)=all_files(ii).handles_out.auc_values(this_jj).auc_coh;
%                 end
%             end
%
%             switch grNo
%                 case 1
%                     bar(bar_offset,mean(these_AUC),'g','LineWidth', 3,'EdgeColor','none')
%                 case 2
%                     bar(bar_offset,mean(these_AUC),'b','LineWidth', 3,'EdgeColor','none')
%                 case 3
%                     bar(bar_offset,mean(these_AUC),'y','LineWidth', 3,'EdgeColor','none')
%             end
%
%
%             %Violin plot
%
%             %                 [mean_out, CIout]=drgViolinPoint(these_AUC,edges,bar_offset,rand_offset,'k','k',1);
%             CI = bootci(1000, {@mean, these_AUC},'type','cper');
%             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%             plot(bar_offset*ones(1,length(these_AUC)),these_AUC,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%
%             %                                 %Save data for glm and ranksum
%
%             glm_AUC.data(glm_ii+1:glm_ii+length(these_AUC))=these_AUC;
%             glm_AUC.group(glm_ii+1:glm_ii+length(these_AUC))=grNo*ones(1,length(these_AUC));
%             glm_AUC.perCorr(glm_ii+1:glm_ii+length(these_AUC))=per_ii*ones(1,length(these_AUC));
%             glm_AUC.event(glm_ii+1:glm_ii+length(these_AUC))=evNo*ones(1,length(these_AUC));
%             glm_ii=glm_ii+length(these_AUC);
%
%             id_ii=id_ii+1;
%             input_data(id_ii).data=these_AUC;
%             input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
%
%
%         end
%         bar_offset = bar_offset + 2;
%
%     end
%     bar_offset = bar_offset + 3;
%
%     %     end
%
%     title(['Average coherence auROC for each electrode calculated per mouse for ' bandwidth_names{bwii}])
%
%
%     %Annotations identifying groups
%     x_interval=0.8/ii_gr_included;
%     for ii=1:ii_gr_included
%         annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
%     end
%
%     %Proficient/Naive annotations
%     annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
%     annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
%
%
%     xticks([2 4 6 10 12 14 21 23 25 29 31 33])
%     xticklabels({'nwtS+', 'nHS+', 'nKOS+', 'pwtS+', 'pHS+', 'pKOS+', 'nwtS-', 'nHS-', 'nKOS-', 'pwtS-', 'pHS-', 'pKOS-'})
%
%     ylim([0 0.5])
%     ylabel('auROC')
%
%     %Perform the glm
%     fprintf(1, ['glm for average coherence auROC for each electrode calculated per mouse for '  bandwidth_names{bwii} '\n'])
%
%     fprintf(1, ['\n\nglm for auROC for Theta/' bandwidth_names{bwii} '\n'])
%     tbl = table(glm_AUC.data',glm_AUC.group',glm_AUC.perCorr',glm_AUC.event',...
%         'VariableNames',{'MI','group','perCorr','event'});
%     mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
%         ,'CategoricalVars',[2,3,4])
%
%
%     %Do the ranksum/t-test
%     fprintf(1, ['\n\nRanksum or t-test p values for average coherence auROC for each electrode calculated per mouse for ' bandwidth_names{bwii} '\n'])
%     [output_data] = drgMutiRanksumorTtest(input_data);
%
%
% end


%Plot the bounded lines
maxlP=-200000;
minlP=200000;

frequency=all_files(1).handles_out.frequency;

figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);

set(hFig, 'units','normalized','position',[.1 .5 .4 .4])

set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
hold on

for evNo=1:2
    
    
    for per_ii=2:-1:1      %performance bins. blue = naive, red = proficient
        

        
        
        grNo=1
        
        %Get the deltaCxy_af
        %Get these MI values
        these_deltaCxy_af=[];
        ii_coh=0;
        for ii=1:length(FileName)
            this_jj=[];
            for jj=1:all_files(ii).handles_out.dcohaf_ii
                
                if all_files(ii).handles_out.dcohaf_values(jj).bwii==bwii
                    if all_files(ii).handles_out.dcohaf_values(jj).evNo==evNo
                        if all_files(ii).handles_out.dcohaf_values(jj).per_ii==per_ii
                            if all_files(ii).handles_out.dcohaf_values(jj).groupNo==grNo
                                this_jj=jj;
                            end
                        end
                    end
                end
            end
            if (~isempty(this_jj))&(~isempty(all_files(ii).handles_out.dcohaf_values(this_jj).dcohaf))
                ii_coh=ii_coh+1;
                these_deltaCxy_af(ii_coh,:)=all_files(ii).handles_out.dcohaf_values(this_jj).dcohaf;
            end
        end
        
        
        
        
        
        mean_deltaCxy=[];
        mean_deltaCxy=mean(these_deltaCxy_af,1);
        
        CI=[];
        CI = bootci(1000, {@mean, these_deltaCxy_af})';
        maxlP=max([maxlP max(CI(:))]);
        minlP=min([minlP min(CI(:))]);
        CI(:,1)= mean_deltaCxy'-CI(:,1);
        CI(:,2)=CI(:,2)- mean_deltaCxy';
        
        
        
        
        if evNo==2
            if per_ii==1
                %S+ Proficient
                [hlCR, hpCR] = boundedline(frequency',mean_deltaCxy', CI, 'cmap',[158/255 31/255 99/255]);
                %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
            else
                %S- Naive
                [hlCR, hpCR] = boundedline(frequency',mean_deltaCxy', CI, 'cmap',[238/255 111/255 179/255]);
                %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
            end
        else
            if per_ii==1
                %S- Proficient
                [hlCR, hpCR] = boundedline(frequency',mean_deltaCxy', CI, 'cmap',[0 114/255 178/255]);
                %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
            else
                %S- naive
                [hlCR, hpCR] = boundedline(frequency',mean_deltaCxy', CI, 'cmap',[80/255 194/255 255/255]);
                %                     bar(bar_offset,mean(these_coh),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
            end
        end
        
        pffft=1;
        
      
    end
    
end

title(['delta coherence '])
xlabel('Frequency (Hz')
ylabel('delta coherence')



pffft=1;




pffft=1;
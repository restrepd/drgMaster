function drgSummaryBatchPLVCaMKII
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
% hippPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/PLV/';
% 
% %Files
% FileName{1}='CaMKIIPLVCamKIAPEB0301202_out.mat';
% FileName{2}='CaMKIIEBAPPLV02282021_out.mat';
% FileName{3}='CaMKIIPAEAPLV03062021_out.mat';
% FileName{4}='CaMKIIEAPAPLV03112021_out.mat';
% FileName{5}='CaMKIIpz1EAPAPLV03052021_out.mat';
% FileName{6}='CaMKIIpz1PAEAPLV03042021_out.mat';
% FileName{7}='CaMKIIpzz1EAPAPLV02262021_out.mat';
% FileName{8}='CaMKIIpzz1propylacecPLV02222021_out.mat';

%Location of files
% hippPathName='E:\CaMKIIpaper\datos sumarry\coherence\';
hippPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/PLV80/';

%Files
FileName{1}='CaMKIIPLVCamKIAPEB8005012021_out.mat';
FileName{2}='CaMKIIEBAPPLV04292021_out.mat';
FileName{3}='CaMKIIPAEAPLV8004282021_out.mat';
FileName{4}='CaMKIIEAPAPLV8004292021_out.mat';
FileName{5}='CaMKIIpz1EAPAPLV8004272021_out.mat';
FileName{6}='CaMKIIpz1PAEAPLV8004262021_out.mat';
FileName{7}='CaMKII80pzz1EAPAPLV04242021_out.mat';
FileName{8}='CaMKIIpzz1PAEA80PLV04242021_out.mat';



%Load data
all_files=[];

for ii=1:length(FileName)
    load([hippPathName FileName{ii}])
    all_files(ii).delta_PLV_odor_minus_ref=delta_PLV_odor_minus_ref;
    all_files(ii).delta_phase_odor=delta_phase_odor;
    all_files(ii).delta_PLV_odor_minus_ref_per_mouse=delta_PLV_odor_minus_ref_per_mouse;
    all_files(ii).delta_phase_odor_per_mouse=delta_phase_odor_per_mouse;
    all_files(ii).group_no_per_mouse=group_no_per_mouse;
end

figNo=0;

%Now plot the average delta PLV for each odor pair
%(including all sessions for each mouse)
edges=[0:0.001:0.02];
rand_offset=0.8;


for bwii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_PLV=[];
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
                
                %Get these PLV values
                if (grNo==3)&(evNo==2)&(per_ii==1)
                    pffft=1;
                end
                
                these_PLV=[];
                ii_PLV=0;
                for ii=1:length(FileName)
                    ii_PLV=ii_PLV+1;
                    these_PLV(ii_PLV)=all_files(ii).delta_PLV_odor_minus_ref(grNo,bwii,per_ii,evNo);
                end
                
                
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(these_PLV),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(these_PLV),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(these_PLV),'y','LineWidth', 3,'EdgeColor','none')
                end
                
                
                %Violin plot
                
                %[mean_out, CIout]=drgViolinPoint(these_PLV,edges,bar_offset,rand_offset,'k','k',1);
                CI = bootci(1000, {@mean, these_PLV},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(these_PLV)),these_PLV,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                
                %                                 %Save data for glm and ranksum
                
                glm_PLV.data(glm_ii+1:glm_ii+length(these_PLV))=these_PLV;
                glm_PLV.group(glm_ii+1:glm_ii+length(these_PLV))=grNo*ones(1,length(these_PLV));
                glm_PLV.perCorr(glm_ii+1:glm_ii+length(these_PLV))=per_ii*ones(1,length(these_PLV));
                glm_PLV.event(glm_ii+1:glm_ii+length(these_PLV))=evNo*ones(1,length(these_PLV));
                glm_ii=glm_ii+length(these_PLV);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_PLV;
                input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                
                
            end
            bar_offset = bar_offset + 2;
            
        end
        bar_offset = bar_offset + 3;
        
    end
    
    title(['Average delta PLV for each odor_pair for ' bandwidth_names{bwii}])
    
    
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
    
    ylabel('delta PLV')
    
    
    %Perform the glm
    fprintf(1, ['glm for delta PLV for each odor pair for '  bandwidth_names{bwii} '\n'])
    
    fprintf(1, ['\n\nglm for PRP for' bandwidth_names{bwii} '\n'])
    tbl = table(glm_PLV.data',glm_PLV.group',glm_PLV.perCorr',glm_PLV.event',...
        'VariableNames',{'MI','group','perCorr','event'});
    mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
        ,'CategoricalVars',[2,3,4])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for delta PLV for each odor pair for ' bandwidth_names{bwii} ' hippocampus\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end




%Now plot the average delta phase for each odor pair
%(including all sessions for each mouse)
edges=[0:0.001:0.02];
rand_offset=0.8;


for bwii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_delta_phase=[];
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
                
                %Get these PLV values
                
                %Get these MI values
            these_delta_phase=[];
            ii_delta_phase=0;
            for ii=1:length(FileName)
                    ii_delta_phase=ii_delta_phase+1;
                    these_delta_phase(ii_delta_phase)=all_files(ii).delta_phase_odor(grNo,bwii,per_ii,evNo);
            end
                
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(these_delta_phase),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(these_delta_phase),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(these_delta_phase),'y','LineWidth', 3,'EdgeColor','none')
                end
                
                
                %Violin plot
                
                %[mean_out, CIout]=drgViolinPoint(these_delta_phase,edges,bar_offset,rand_offset,'k','k',1);
                CI = bootci(1000, {@mean, these_delta_phase},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(these_delta_phase)),these_delta_phase,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                
                %                                 %Save data for glm and ranksum
                
                glm_delta_phase.data(glm_ii+1:glm_ii+length(these_delta_phase))=these_delta_phase;
                glm_delta_phase.group(glm_ii+1:glm_ii+length(these_delta_phase))=grNo*ones(1,length(these_delta_phase));
                glm_delta_phase.perCorr(glm_ii+1:glm_ii+length(these_delta_phase))=per_ii*ones(1,length(these_delta_phase));
                glm_delta_phase.event(glm_ii+1:glm_ii+length(these_delta_phase))=evNo*ones(1,length(these_delta_phase));
                glm_ii=glm_ii+length(these_delta_phase);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_delta_phase;
                input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                
                
            end
            bar_offset = bar_offset + 2;
            
        end
        bar_offset = bar_offset + 3;
        
    end
    
    title(['Average delta phase for each odor_pair for ' bandwidth_names{bwii}])
    
    
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
    
    ylabel('delta phase')
    
    
    %Perform the glm
    fprintf(1, ['glm for delta PLV for each odor pair for '  bandwidth_names{bwii} '\n'])
    
    fprintf(1, ['\n\nglm for delta phase for' bandwidth_names{bwii} '\n'])
    tbl = table(glm_delta_phase.data',glm_delta_phase.group',glm_delta_phase.perCorr',glm_delta_phase.event',...
        'VariableNames',{'MI','group','perCorr','event'});
    mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
        ,'CategoricalVars',[2,3,4])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for delta phase for each odor pair for ' bandwidth_names{bwii} ' hippocampus\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end


%Now plot the  delta PLV for each odor pair/mouse
%(including all sessions for each mouse)
edges=[-0.6:0.05:0.6];
rand_offset=0.7;

for bwii=[1 2 4]    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_PLV=[];
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
      
                
                %Get these PLV values
                these_PLV=[];
                ii_PLV=0;
                for ii=1:length(FileName)
                    these_PLVs=zeros(1,sum(group_no_per_mouse==grNo));
                    these_PLVs(1,:)=all_files(ii).delta_PLV_odor_minus_ref_per_mouse(group_no_per_mouse==grNo,bwii,per_ii,evNo);
                    these_PLV(ii_PLV+1:ii_PLV+length(these_PLVs))=these_PLVs;
                    ii_PLV=ii_PLV+length(these_PLVs);
                end
                
                
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(these_PLV),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(these_PLV),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(these_PLV),'y','LineWidth', 3,'EdgeColor','none')
                end
                
                
                %Violin plot
                
                [mean_out, CIout]=drgViolinPoint(these_PLV,edges,bar_offset,rand_offset,'k','k',3);
%                 CI = bootci(1000, {@mean, these_PLV},'type','cper');
%                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%                 plot(bar_offset*ones(1,length(these_PLV)),these_PLV,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%                 
                
                %                                 %Save data for glm and ranksum
                
                glm_PLV.data(glm_ii+1:glm_ii+length(these_PLV))=these_PLV;
                glm_PLV.group(glm_ii+1:glm_ii+length(these_PLV))=grNo*ones(1,length(these_PLV));
                glm_PLV.perCorr(glm_ii+1:glm_ii+length(these_PLV))=per_ii*ones(1,length(these_PLV));
                glm_PLV.event(glm_ii+1:glm_ii+length(these_PLV))=evNo*ones(1,length(these_PLV));
                glm_ii=glm_ii+length(these_PLV);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_PLV;
                input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                
                
            end
            bar_offset = bar_offset + 1;
            
        end
        bar_offset = bar_offset + 1;
        
    end
    
    title(['Average delta PLV for each odor_pair/mouse for ' bandwidth_names{bwii}])
    
    
    xticks([1 2 3 5 6 7 10 11 12 14 15 16])
    xticklabels({'nwtS+', 'nHS+', 'nKOS+', 'pwtS+', 'pHS+', 'pKOS+', 'nwtS-', 'nHS-', 'nKOS-', 'pwtS-', 'pHS-', 'pKOS-'})
    
    ylabel('delta PLV')
    ylim([-0.6 0.6])
    
    %Perform the glm
    fprintf(1, ['glm for delta PLV for each odor pair/mouse for '  bandwidth_names{bwii} '\n'])
    
    fprintf(1, ['\n\nglm for delta PLV for each odor pair/mouse  for' bandwidth_names{bwii} '\n'])
    tbl = table(glm_PLV.data',glm_PLV.group',glm_PLV.perCorr',glm_PLV.event',...
        'VariableNames',{'MI','group','perCorr','event'});
    mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
        ,'CategoricalVars',[2,3,4])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for delta PLV for each odor pair for ' bandwidth_names{bwii} ' hippocampus\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end


 

%Now plot the average delta phase for each odor pair/mouse
%(including all sessions for each mouse)
edges=[-3:0.2:3];
rand_offset=0.7;

for bwii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    glm_delta_phase=[];
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
                
                %Get these PLV values
                
                %Get these MI values
             these_delta_phase=[];
            ii_delta_phase=0;
            for ii=1:length(FileName)
%                     ii_delta_phase=ii_delta_phase+1;
%                     these_delta_phase(ii_delta_phase)=all_files(ii).delta_phase_odor(grNo,bwii,per_ii,evNo);
                    these_delta_phases=zeros(1,sum(group_no_per_mouse==grNo));
                    these_delta_phases(1,:)=all_files(ii).delta_phase_odor_per_mouse(group_no_per_mouse==grNo,bwii,per_ii,evNo);
                    these_delta_phase(ii_delta_phase+1:ii_delta_phase+length(these_delta_phases))=these_delta_phases;
                    ii_delta_phase=ii_delta_phase+length(these_delta_phases);
            end
            
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(these_delta_phase),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(these_delta_phase),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(these_delta_phase),'y','LineWidth', 3,'EdgeColor','none')
                end
                
                
                %Violin plot
                
                [mean_out, CIout]=drgViolinPoint(these_delta_phase,edges,bar_offset,rand_offset,'k','k',3);
%                 CI = bootci(1000, {@mean, these_delta_phase},'type','cper');
%                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%                 plot(bar_offset*ones(1,length(these_delta_phase)),these_delta_phase,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%                 
                
                %                                 %Save data for glm and ranksum
                
                glm_delta_phase.data(glm_ii+1:glm_ii+length(these_delta_phase))=these_delta_phase;
                glm_delta_phase.group(glm_ii+1:glm_ii+length(these_delta_phase))=grNo*ones(1,length(these_delta_phase));
                glm_delta_phase.perCorr(glm_ii+1:glm_ii+length(these_delta_phase))=per_ii*ones(1,length(these_delta_phase));
                glm_delta_phase.event(glm_ii+1:glm_ii+length(these_delta_phase))=evNo*ones(1,length(these_delta_phase));
                glm_ii=glm_ii+length(these_delta_phase);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_delta_phase;
                input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                
                
            end
            bar_offset = bar_offset + 2;
            
        end
        bar_offset = bar_offset + 3;
        
    end
    
    title(['Average delta phase for each odor_pair/mouse for ' bandwidth_names{bwii}])
    
    
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
    
    ylabel('delta phase')
    
    
    %Perform the glm
    fprintf(1, ['glm for delta PLV for each odor pair/mouse for '  bandwidth_names{bwii} '\n'])
    
    fprintf(1, ['\n\nglm for delta phase for' bandwidth_names{bwii} '\n'])
    tbl = table(glm_delta_phase.data',glm_delta_phase.group',glm_delta_phase.perCorr',glm_delta_phase.event',...
        'VariableNames',{'MI','group','perCorr','event'});
    mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
        ,'CategoricalVars',[2,3,4])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for delta phase for each odor pair/mouse for ' bandwidth_names{bwii} ' hippocampus\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end



pffft=1;
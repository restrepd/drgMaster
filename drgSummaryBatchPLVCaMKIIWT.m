function drgSummaryBatchPLVCaMKIIWT
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

% %Location of files
% % hippPathName='E:\CaMKIIpaper\datos sumarry\coherence\';
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
FileName{2}='CaMKIIEAPAPLV8004292021_out.mat';
FileName{3}='CaMKIIEBAPPLV04292021_out.mat';
FileName{4}='CaMKIIPAEAPLV8004282021_out.mat';
FileName{5}='CaMKIIpz1EAPAPLV8004272021_out.mat';
FileName{6}='CaMKIIpz1PAEAPLV8004262021_out.mat';
FileName{7}='CaMKII80pzz1EAPAPLV04242021_out.mat';
FileName{8}='CaMKIIpzz1PAEA80PLV04242021_out.mat';


%Load the table of mouse numbers
mouse_no_table='/Users/restrepd/Documents/Projects/CaMKII_analysis/Reply_to_reviewers/camkii_mice_per_odor_pair_for_PLV.xlsx';
T_mouse_no = readtable(mouse_no_table);


%Text file for statistical output
fileID = fopen([hippPathName 'drgSummaryBatchPLVCaMKIIWT.txt'],'w');


%Load data
all_files=[];

for ii=1:length(FileName)
    load([hippPathName FileName{ii}])
    all_files(ii).delta_PLV_odor_minus_ref=delta_PLV_odor_minus_ref;
    all_files(ii).delta_phase_odor=delta_phase_odor;
    all_files(ii).delta_PLV_odor_minus_ref_per_mouse=delta_PLV_odor_minus_ref_per_mouse;
    all_files(ii).delta_phase_odor_per_mouse=delta_phase_odor_per_mouse;
    all_files(ii).group_no_per_mouse=group_no_per_mouse;
    all_files(ii).which_mice=which_mice;
end


figNo=0;

%Now plot the average PLV calculated per odor pair
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
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
    hold on
    

    bar_offset = 0;
    
    for evNo=1:2
        
        for per_ii=2:-1:1
            
            grNo=1;
            
            bar_offset = bar_offset + 1;
            
            %Get these MI values
            these_PLV=[];
            ii_PLV=0;
            for ii=1:length(FileName)
                    ii_PLV=ii_PLV+1;
                    these_PLV(ii_PLV)=all_files(ii).delta_PLV_odor_minus_ref(grNo,bwii,per_ii,evNo);
            end
            
            
            if evNo==2
                if per_ii==1
                    %S+ Proficient
                    bar(bar_offset,mean(these_PLV),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                else
                    %S+ Naive
                    bar(bar_offset,mean(these_PLV),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                end
            else
                if per_ii==1
                    %S- Proficient
                    bar(bar_offset,mean(these_PLV),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                else
                    %S- naive
                    bar(bar_offset,mean(these_PLV),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                end
            end
            
            
            %Violin plot
            
            %[mean_out, CIout]=drgViolinPoint(these_PLV,edges,bar_offset,rand_offset,'k','k',1);
            CI = bootci(1000, {@mean, these_PLV},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_PLV)),these_PLV,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            
            
            glm_PLV.data(glm_ii+1:glm_ii+length(these_PLV))=these_PLV;
            %                 glm_PLV.group(glm_ii+1:glm_ii+length(these_PLV))=grNo*ones(1,length(these_PLV));
            glm_PLV.perCorr(glm_ii+1:glm_ii+length(these_PLV))=per_ii*ones(1,length(these_PLV));
            glm_PLV.event(glm_ii+1:glm_ii+length(these_PLV))=evNo*ones(1,length(these_PLV));
            glm_ii=glm_ii+length(these_PLV);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=these_PLV;
            input_data(id_ii).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
            
        end
        bar_offset = bar_offset + 1;
        
    end
    
    title(['Average delta PLV for each electrode calculated per odor pair for ' bandwidth_names{bwii}])
    
    
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
    
    ylabel('delta PLV')
    
    
    %Perform the glm
    fprintf(1, ['glm for average delta PLV for each electrode calculated per odor pair for '  bandwidth_names{bwii} '\n'])
    
    fprintf(1, ['\n\nglm for PLV for' bandwidth_names{bwii} '\n'])
    tbl = table(glm_PLV.data',glm_PLV.perCorr',glm_PLV.event',...
        'VariableNames',{'dPLV','perCorr','event'});
    mdl = fitglm(tbl,'dPLV~perCorr+event+perCorr*event'...
        ,'CategoricalVars',[2,3])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for average delta PLV for each electrode calculated per odor pair for ' bandwidth_names{bwii} ' hippocampus\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end

 

%Now plot the average delta phase calculated per odor pair
for bwii=[1 2 4]    %for amplitude bandwidths (beta, low gamma, high gamma)
    
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
    
     ax=gca;ax.LineWidth=3;
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
    hold on
    
    
    bar_offset = 0;
    
    for evNo=1:2
        
        for per_ii=2:-1:1
            
            grNo=1;
            
            bar_offset = bar_offset + 1;
            
            %Get these MI values
            these_delta_phase=[];
            ii_delta_phase=0;
            for ii=1:length(FileName)
                    ii_delta_phase=ii_delta_phase+1;
                    these_delta_phase(ii_delta_phase)=all_files(ii).delta_phase_odor(grNo,bwii,per_ii,evNo);
            end
            
            
            if evNo==2
                if per_ii==1
                    %S+ Proficient
                    bar(bar_offset,circ_mean(these_delta_phase'),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                else
                    %S+ Naive
                    bar(bar_offset,circ_mean(these_delta_phase'),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                end
            else
                if per_ii==1
                    %S- Proficient
                    bar(bar_offset,circ_mean(these_delta_phase'),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                else
                    %S- naive
                    bar(bar_offset,circ_mean(these_delta_phase'),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                end
            end
             
            
            %Violin plot
            
            %[mean_out, CIout]=drgViolinPoint(these_delta_phase,edges,bar_offset,rand_offset,'k','k',1);
            CI = bootci(1000, {@circ_mean, these_delta_phase'},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            plot(bar_offset*ones(1,length(these_delta_phase')),these_delta_phase','o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            
            
            
            glm_delta_phase.data(glm_ii+1:glm_ii+length(these_delta_phase))=these_delta_phase;
            %                 glm_delta_phase.group(glm_ii+1:glm_ii+length(these_delta_phase))=grNo*ones(1,length(these_delta_phase));
            glm_delta_phase.perCorr(glm_ii+1:glm_ii+length(these_delta_phase))=per_ii*ones(1,length(these_delta_phase));
            glm_delta_phase.event(glm_ii+1:glm_ii+length(these_delta_phase))=evNo*ones(1,length(these_delta_phase));
            glm_ii=glm_ii+length(these_delta_phase);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=these_delta_phase;
            input_data(id_ii).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
            
        end
        bar_offset = bar_offset + 1;
        
    end
    
    title(['Average delta phase for each electrode calculated per mouse for ' bandwidth_names{bwii}])
    
    
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
    
    ylabel('delta PLV')
    
    
    %Perform the glm
    fprintf(1, ['glm for average delta phase for each electrode calculated per mouse for '  bandwidth_names{bwii} '\n'])
    
    fprintf(1, ['\n\nglm for phase for' bandwidth_names{bwii} '\n'])
    tbl = table(glm_delta_phase.data',glm_delta_phase.perCorr',glm_delta_phase.event',...
        'VariableNames',{'delta_phase','perCorr','event'});
    mdl = fitglm(tbl,'delta_phase~perCorr+event+perCorr*event'...
        ,'CategoricalVars',[2,3])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for average delta phase for each electrode calculated per mouse for ' bandwidth_names{bwii} ' hippocampus\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end


%Now plot the per odor pair/mouse in violin plots

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
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
    hold on
    

    bar_offset = 0;
    
    for evNo=1:2
        
        for per_ii=2:-1:1
            
            grNo=1;
            
            bar_offset = bar_offset + 1;
            
            %Get these PLV values
            these_PLV=[];
            ii_PLV=0;
            these_mice=[];
            for ii=1:length(FileName)
                
                group_no_per_mouse=all_files(ii).group_no_per_mouse;
                these_PLVs=zeros(1,sum(group_no_per_mouse==grNo));
                these_PLVs(1,:)=all_files(ii).delta_PLV_odor_minus_ref_per_mouse(group_no_per_mouse==grNo,bwii,per_ii,evNo);
                these_PLV(ii_PLV+1:ii_PLV+length(these_PLVs))=these_PLVs;
                ii_PLV=ii_PLV+length(these_PLVs);
                
                
                %Assign to overall mouse numbers
                these_mouse_no_per_op=T_mouse_no.mouse_no_per_op(T_mouse_no.odor_pair_no==ii);
                these_mouse_no=T_mouse_no.mouse_no(T_mouse_no.odor_pair_no==ii);
                these_mouse_nos=all_files(ii).which_mice(group_no_per_mouse==grNo);
                these_mouse_no_overall=[];
                for ii_mouse=1:length(these_mouse_nos)
                    these_mouse_no_overall(ii_mouse)=these_mouse_no(find(these_mouse_nos(ii_mouse)==these_mouse_no_per_op));
                end
                these_mice=[these_mice these_mouse_no_overall];
            end
            
            mean_PLV_mm=[];
            if per_ii==1
                mean_PLV_mm_per_ii1=[];
            else
                mean_PLV_mm_per_ii2=[];
            end
            
            for mouseNo=unique(these_mice)
                if per_ii==1
                    mean_PLV_mm_per_ii1=[mean_PLV_mm_per_ii1 mean(these_PLV(these_mice==mouseNo))];
                else
                    mean_PLV_mm_per_ii2=[mean_PLV_mm_per_ii2 mean(these_PLV(these_mice==mouseNo))];
                end
            end
            
            if evNo==2
                if per_ii==1
                    %S+ Proficient
                    bar(bar_offset,mean(these_PLV),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                else
                    %S+ Naive
                    bar(bar_offset,mean(these_PLV),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                end
            else
                if per_ii==1
                    %S- Proficient
                    bar(bar_offset,mean(these_PLV),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                else
                    %S- naive
                    bar(bar_offset,mean(these_PLV),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                end
            end
            
            if per_ii==1
                for mouseNo=1:length(mean_PLV_mm_per_ii2)
                    plot([bar_offset-1 bar_offset],[mean_PLV_mm_per_ii2(mouseNo) mean_PLV_mm_per_ii1(mouseNo)],'-ok','MarkerSize',6,'LineWidth',2,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor',[0.6 0.6 0.6],'Color',[0.6 0.6 0.6])
                end
            end
            
            %Violin plot
            
            [mean_out, CIout]=drgViolinPoint(these_PLV,edges,bar_offset,rand_offset,'k','k',3);
%             CI = bootci(1000, {@mean, these_PLV},'type','cper');
%             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%             plot(bar_offset*ones(1,length(these_PLV)),these_PLV,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%             
            
            
            glm_PLV.data(glm_ii+1:glm_ii+length(these_PLV))=these_PLV;
            %                 glm_PLV.group(glm_ii+1:glm_ii+length(these_PLV))=grNo*ones(1,length(these_PLV));
            glm_PLV.perCorr(glm_ii+1:glm_ii+length(these_PLV))=per_ii*ones(1,length(these_PLV));
            glm_PLV.event(glm_ii+1:glm_ii+length(these_PLV))=evNo*ones(1,length(these_PLV));
            glm_PLV.mouse_no(glm_ii+1:glm_ii+length(these_PLV))=these_mice;
            glm_ii=glm_ii+length(these_PLV);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=these_PLV;
            input_data(id_ii).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
            
        end
        bar_offset = bar_offset + 1;
        
    end
    
    title(['Average delta PLV for each odor pair/mouse for ' bandwidth_names{bwii}])
    
    
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
    
    ylabel('delta PLV')
    ylim([-0.6 0.6])
    
    %Perform the glm
    fprintf(1, ['glm for average delta PLV per mouse per odor pair for '  bandwidth_names{bwii} '\n'])
    fprintf(fileID, ['glm for average delta PLV per mouse per odor pair for '  bandwidth_names{bwii} '\n']);
    
    fprintf(1, ['\n\nglm for PLV for' bandwidth_names{bwii} '\n'])
    tbl = table(glm_PLV.data',glm_PLV.perCorr',glm_PLV.event',...
        'VariableNames',{'dPLV','naive_vs_proficient','sp_vs_sm'});
    mdl = fitglm(tbl,'dPLV~naive_vs_proficient+sp_vs_sm+naive_vs_proficient*sp_vs_sm'...
        ,'CategoricalVars',[2,3])
    
       txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
    fprintf(fileID,'%s\n', txt);
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for average delta PLV  per mouse per odor pair for ' bandwidth_names{bwii} '\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for average delta PLV  per mouse per odor pair for ' bandwidth_names{bwii} '\n']);
    
    [output_data] = drgMutiRanksumorTtest(input_data, fileID);
    
    
    %Nested ANOVAN
    %https://www.mathworks.com/matlabcentral/answers/491365-within-between-subjects-in-anovan
    nesting=[0 0 0; ... % This line indicates that group factor is not nested in any other factor.
         0 0 0; ... % This line indicates that perCorr is not nested in any other factor.
         1 1 0];    % This line indicates that mouse_no (the last factor) is nested under group, perCorr and event
                    % (the 1 in position 1 on the line indicates nesting under the first factor).
    figNo=figNo+1;
                     
    [p anovanTbl stats]=anovan(glm_PLV.data,{glm_PLV.perCorr glm_PLV.event glm_PLV.mouse_no},...
    'model','interaction',...
    'nested',nesting,...
    'varnames',{'naive_vs_proficient', 'sp_vs_sm','mouse_no'});
  
    fprintf(fileID, ['\n\nNested ANOVAN for average delta PLV  per mouse per odor pair for '  bandwidth_names{bwii} '\n'])
    drgWriteANOVANtbl(anovanTbl,fileID);
    fprintf(fileID, '\n\n')
end

 

%Now plot the average delta phase calculated per odor pair

edges=[-3:0.2:3];
rand_offset=0.7;
all_theta_phases=[];

for bwii=[1 2 4]    %for amplitude bandwidths (beta, low gamma, high gamma)
    
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
    
    ax=gca;ax.LineWidth=3;
    
    
    set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
    hold on
    
    
    bar_offset = 0;
    
    for evNo=1:2
        
        for per_ii=2:-1:1
            
            grNo=1;
            
            bar_offset = bar_offset + 1;
            
            %Get these delta phase values

            
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
            
            if bwii==1
                all_theta_phases=[all_theta_phases these_delta_phase];
            end
            
            if evNo==2
                if per_ii==1
                    %S- Proficient
                    bar(bar_offset,circ_mean(these_delta_phase'),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                else
                    %S- Naive
                    bar(bar_offset,circ_mean(these_delta_phase'),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                end
            else
                if per_ii==1
                    %S+ Proficient
                    bar(bar_offset,circ_mean(these_delta_phase'),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                else
                    %S+ naive
                    bar(bar_offset,circ_mean(these_delta_phase'),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                end
            end
             
            
            %Violin plot
            
            [mean_out, CIout]=drgViolinPoint(these_delta_phase,edges,bar_offset,rand_offset,'k','k',3);
%             CI = bootci(1000, {@circ_mean, these_delta_phase'},'type','cper');
%             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%             plot(bar_offset*ones(1,length(these_delta_phase')),these_delta_phase','o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%             
            
            
            glm_delta_phase.data(glm_ii+1:glm_ii+length(these_delta_phase))=these_delta_phase;
            %                 glm_delta_phase.group(glm_ii+1:glm_ii+length(these_delta_phase))=grNo*ones(1,length(these_delta_phase));
            glm_delta_phase.perCorr(glm_ii+1:glm_ii+length(these_delta_phase))=per_ii*ones(1,length(these_delta_phase));
            glm_delta_phase.event(glm_ii+1:glm_ii+length(these_delta_phase))=evNo*ones(1,length(these_delta_phase));
            glm_ii=glm_ii+length(these_delta_phase);
            
            id_ii=id_ii+1;
            input_data(id_ii).data=these_delta_phase;
            input_data(id_ii).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
            
        end
        bar_offset = bar_offset + 1;
        
    end
    
    title(['Average delta phase for odor pair/mouse for ' bandwidth_names{bwii}])
    
    
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
    
    ylabel('delta phase')
    ylim([-1 1])
    
    
    %Perform the glm
    fprintf(1, ['glm for average delta phase for odor pair/mouse for '  bandwidth_names{bwii} '\n'])
    
    fprintf(1, ['\n\nglm for phase for' bandwidth_names{bwii} '\n'])
    tbl = table(glm_delta_phase.data',glm_delta_phase.perCorr',glm_delta_phase.event',...
        'VariableNames',{'delta_phase','perCorr','event'});
    mdl = fitglm(tbl,'delta_phase~perCorr+event+perCorr*event'...
        ,'CategoricalVars',[2,3])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for odor pair/mouse mouse for ' bandwidth_names{bwii} ' hippocampus\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end

edges=[-1:0.1:1]

%Plot the average
figNo = figNo +1;

try
    close(figNo)
catch
end
hFig=figure(figNo);

ax=gca;ax.LineWidth=3;


set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
hold on

histogram(all_theta_phases,edges)

[h,p]=ttest(all_theta_phases);

fprintf(1, ['mean delta theta phase %d radians, %d msec\n\n'],mean(all_theta_phases),100*mean(all_theta_phases)/(2*pi()))
fprintf(1, ['p value = %d for theta delta phase vs zero\n\n'],p)
    
ylabel('Number of observations')
xlabel('delta theta phase hipp - pre')

fclose(fileID)

pffft=1;
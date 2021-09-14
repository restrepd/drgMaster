function drgSummaryBatchPLVvsCohCaMKIIWT
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

ii_pval=0;
pval=[];

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
% PLVPathName='E:\CaMKIIpaper\datos sumarry\coherence\';
PLVPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/PLV80/';

%Files
FileNamePLV{1}='CaMKIIPLVCamKIAPEB8005012021_out.mat';
FileNamePLV{2}='CaMKIIEBAPPLV04292021_out.mat';
FileNamePLV{3}='CaMKIIPAEAPLV8004282021_out.mat';
FileNamePLV{4}='CaMKIIEAPAPLV8004292021_out.mat';
FileNamePLV{5}='CaMKIIpz1EAPAPLV8004272021_out.mat';
FileNamePLV{6}='CaMKIIpz1PAEAPLV8004262021_out.mat';
FileNamePLV{7}='CaMKII80pzz1EAPAPLV04242021_out.mat';
FileNamePLV{8}='CaMKIIpzz1PAEA80PLV04242021_out.mat';

%Location of files
% hippPathName='E:\CaMKIIpaper\datos sumarry\coherence\';
cohPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Coherence new/';

cohFileName{1}='CaMKIIacetocohe02012021_out80.mat';
cohFileName{2}='CaMKIIethylbenacetocohe2262021_out80.mat';
cohFileName{3}='CaMKIIPAEAcohe02082021_out80.mat';
cohFileName{4}='CaMKIIEAPAcohe2262021_out80.mat';
cohFileName{5}='CaMKIIpz1EAcohe02142021_out80.mat';
cohFileName{6}='CaMKIIPZ1PAEAcohe202102021_out80.mat';
cohFileName{7}='CaMKIIpzz1EAPAcohe02112021_out80.mat';
cohFileName{8}='CaMKIIpzz1propylacecohe02092021_out80.mat';


%Load PLV data
all_filesPLV=[];

for ii=1:length(FileNamePLV)
    load([PLVPathName FileNamePLV{ii}])
    all_filesPLV(ii).delta_PLV_odor_minus_ref=delta_PLV_odor_minus_ref;
    all_filesPLV(ii).delta_phase_odor=delta_phase_odor;
    all_filesPLV(ii).delta_PLV_odor_minus_ref_per_mouse=delta_PLV_odor_minus_ref_per_mouse;
    all_filesPLV(ii).delta_phase_odor_per_mouse=delta_phase_odor_per_mouse;
    all_filesPLV(ii).group_no_per_mouse=group_no_per_mouse;
    all_filesPLV(ii).which_mice=which_mice;
end

%Load coherence data
all_files=[];

for ii=1:length(cohFileName)
    load([cohPathName cohFileName{ii}])
    all_filescoh(ii).handles_out=handles_out;
end

figNo=0;

edges=[-0.6:0.05:0.6];
rand_offset=0.7;

%
% %Now plot the average PLV calculated per odor pair
% edges=[-0.6:0.05:0.6];
% rand_offset=0.7;
%
%
% for bwii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
%
%     glm_PLV=[];
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
%      ax=gca;ax.LineWidth=3;
%
%
%     set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
%     hold on
%
%
%     bar_offset = 0;
%
%     for evNo=1:2
%
%         for per_ii=2:-1:1
%
%             grNo=1;
%
%             bar_offset = bar_offset + 1;
%
%             %Get these MI values
%             these_PLV=[];
%             ii_PLV=0;
%             for ii=1:length(FileName)
%                     ii_PLV=ii_PLV+1;
%                     these_PLV(ii_PLV)=all_filesPLV(ii).delta_PLV_odor_minus_ref(grNo,bwii,per_ii,evNo);
%             end
%
%
%             if evNo==2
%                 if per_ii==1
%                     %S+ Proficient
%                     bar(bar_offset,mean(these_PLV),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
%                 else
%                     %S+ Naive
%                     bar(bar_offset,mean(these_PLV),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
%                 end
%             else
%                 if per_ii==1
%                     %S- Proficient
%                     bar(bar_offset,mean(these_PLV),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
%                 else
%                     %S- naive
%                     bar(bar_offset,mean(these_PLV),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
%                 end
%             end
%
%
%             %Violin plot
%
%             %[mean_out, CIout]=drgViolinPoint(these_PLV,edges,bar_offset,rand_offset,'k','k',1);
%             CI = bootci(1000, {@mean, these_PLV},'type','cper');
%             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%             plot(bar_offset*ones(1,length(these_PLV)),these_PLV,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%
%
%
%             glm_PLV.data(glm_ii+1:glm_ii+length(these_PLV))=these_PLV;
%             %                 glm_PLV.group(glm_ii+1:glm_ii+length(these_PLV))=grNo*ones(1,length(these_PLV));
%             glm_PLV.perCorr(glm_ii+1:glm_ii+length(these_PLV))=per_ii*ones(1,length(these_PLV));
%             glm_PLV.event(glm_ii+1:glm_ii+length(these_PLV))=evNo*ones(1,length(these_PLV));
%             glm_ii=glm_ii+length(these_PLV);
%
%             id_ii=id_ii+1;
%             input_data(id_ii).data=these_PLV;
%             input_data(id_ii).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
%
%         end
%         bar_offset = bar_offset + 1;
%
%     end
%
%     title(['Average delta PLV for each electrode calculated per mouse for ' bandwidth_names{bwii}])
%
%
%     %     %Annotations identifying groups
%     %     x_interval=0.8/ii_gr_included;
%     %     for ii=1:ii_gr_included
%     %         annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
%     %     end
%     %
%     %     %Proficient/Naive annotations
%     %     annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
%     %     annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
%
%
%     xticks([1 2 4 5])
%     xticklabels({'nS+', 'pS+','nS-', 'pS-'})
%
%     ylabel('delta PLV')
%
%
%     %Perform the glm
%     fprintf(1, ['glm for average delta PLV for each electrode calculated per mouse for '  bandwidth_names{bwii} '\n'])
%
%     fprintf(1, ['\n\nglm for PLV for' bandwidth_names{bwii} '\n'])
%     tbl = table(glm_PLV.data',glm_PLV.perCorr',glm_PLV.event',...
%         'VariableNames',{'dPLV','perCorr','event'});
%     mdl = fitglm(tbl,'dPLV~perCorr+event+perCorr*event'...
%         ,'CategoricalVars',[2,3])
%
%
%     %Do the ranksum/t-test
%     fprintf(1, ['\n\nRanksum or t-test p values for average delta PLV for each electrode calculated per mouse for ' bandwidth_names{bwii} ' hippocampus\n'])
%     [output_data] = drgMutiRanksumorTtest(input_data);
%
%
% end
%
%
%
% %Now plot the average delta phase calculated per odor pair
% for bwii=1:4    %for amplitude bandwidths (beta, low gamma, high gamma)
%
%     glm_delta_phase=[];
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
%      ax=gca;ax.LineWidth=3;
%
%
%     set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
%     hold on
%
%
%     bar_offset = 0;
%
%     for evNo=1:2
%
%         for per_ii=2:-1:1
%
%             grNo=1;
%
%             bar_offset = bar_offset + 1;
%
%             %Get these MI values
%             these_delta_phase=[];
%             ii_delta_phase=0;
%             for ii=1:length(FileName)
%                     ii_delta_phase=ii_delta_phase+1;
%                     these_delta_phase(ii_delta_phase)=all_filesPLV(ii).delta_phase_odor(grNo,bwii,per_ii,evNo);
%             end
%
%
%             if evNo==2
%                 if per_ii==1
%                     %S+ Proficient
%                     bar(bar_offset,circ_mean(these_delta_phase'),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
%                 else
%                     %S+ Naive
%                     bar(bar_offset,circ_mean(these_delta_phase'),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
%                 end
%             else
%                 if per_ii==1
%                     %S- Proficient
%                     bar(bar_offset,circ_mean(these_delta_phase'),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
%                 else
%                     %S- naive
%                     bar(bar_offset,circ_mean(these_delta_phase'),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
%                 end
%             end
%
%
%             %Violin plot
%
%             %[mean_out, CIout]=drgViolinPoint(these_delta_phase,edges,bar_offset,rand_offset,'k','k',1);
%             CI = bootci(1000, {@circ_mean, these_delta_phase'},'type','cper');
%             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%             plot(bar_offset*ones(1,length(these_delta_phase')),these_delta_phase','o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%
%
%
%             glm_delta_phase.data(glm_ii+1:glm_ii+length(these_delta_phase))=these_delta_phase;
%             %                 glm_delta_phase.group(glm_ii+1:glm_ii+length(these_delta_phase))=grNo*ones(1,length(these_delta_phase));
%             glm_delta_phase.perCorr(glm_ii+1:glm_ii+length(these_delta_phase))=per_ii*ones(1,length(these_delta_phase));
%             glm_delta_phase.event(glm_ii+1:glm_ii+length(these_delta_phase))=evNo*ones(1,length(these_delta_phase));
%             glm_ii=glm_ii+length(these_delta_phase);
%
%             id_ii=id_ii+1;
%             input_data(id_ii).data=these_delta_phase;
%             input_data(id_ii).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
%
%         end
%         bar_offset = bar_offset + 1;
%
%     end
%
%     title(['Average delta phase for each electrode calculated per mouse for ' bandwidth_names{bwii}])
%
%
%     %     %Annotations identifying groups
%     %     x_interval=0.8/ii_gr_included;
%     %     for ii=1:ii_gr_included
%     %         annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
%     %     end
%     %
%     %     %Proficient/Naive annotations
%     %     annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
%     %     annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
%
%
%     xticks([1 2 4 5])
%     xticklabels({'nS+', 'pS+','nS-', 'pS-'})
%
%     ylabel('delta PLV')
%
%
%     %Perform the glm
%     fprintf(1, ['glm for average delta phase for each electrode calculated per mouse for '  bandwidth_names{bwii} '\n'])
%
%     fprintf(1, ['\n\nglm for phase for' bandwidth_names{bwii} '\n'])
%     tbl = table(glm_delta_phase.data',glm_delta_phase.perCorr',glm_delta_phase.event',...
%         'VariableNames',{'delta_phase','perCorr','event'});
%     mdl = fitglm(tbl,'delta_phase~perCorr+event+perCorr*event'...
%         ,'CategoricalVars',[2,3])
%
%
%     %Do the ranksum/t-test
%     fprintf(1, ['\n\nRanksum or t-test p values for average delta phase for each electrode calculated per mouse for ' bandwidth_names{bwii} ' hippocampus\n'])
%     [output_data] = drgMutiRanksumorTtest(input_data);
%
%
% end


%Now plot the per odor pair/mouse correlation of PLV and coherence
for grNo=1:3
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
        
        all_cohs=[];
        all_PLVs=[];
        
        for evNo=1:2
            
            for per_ii=2:-1:1
                
                
                
                bar_offset = bar_offset + 1;
                
                %Get these PLV and coherence values and sort them per mouse
                these_PLV=[];
                ii_PLV=0;
                these_coh=[];
                ii_coh=0;
                
                for ii=1:length(FileNamePLV)
                    these_PLVs=zeros(1,sum(group_no_per_mouse==grNo));
                    these_PLV_mice=zeros(1,sum(group_no_per_mouse==grNo));
                    these_PLVs(1,:)=all_filesPLV(ii).delta_PLV_odor_minus_ref_per_mouse(group_no_per_mouse==grNo,bwii,per_ii,evNo);
                    these_PLV_mice(1,:)=all_filesPLV(ii).which_mice(group_no_per_mouse==grNo);
                    
                    
                    
                    this_jj=[];
                    for jj=1:all_filescoh(ii).handles_out.dcoh_ii
                        
                        if all_filescoh(ii).handles_out.dcoh_values(jj).pacii==bwii
                            if all_filescoh(ii).handles_out.dcoh_values(jj).evNo==evNo
                                if all_filescoh(ii).handles_out.dcoh_values(jj).per_ii==per_ii
                                    
                                    
                                    if all_filescoh(ii).handles_out.dcoh_values(jj).groupNo==grNo
                                        this_jj=jj;
                                    end
                                    
                                end
                                
                            end
                        end
                    end
                    
                    
                    these_cohs=all_filescoh(ii).handles_out.dcoh_values(this_jj).dcoh_per_mouse;
                    these_coh_mice=all_filescoh(ii).handles_out.dcoh_values(this_jj).mouseNo;
                    
                    for kk=1:length(these_coh_mice)
                        this_PLVmouse_ii=find(these_PLV_mice==these_coh_mice(kk));
                        if ~isempty(this_PLVmouse_ii)
                            these_PLV(ii_PLV+1)=these_PLVs(this_PLVmouse_ii);
                            ii_PLV=ii_PLV+1;
                            these_coh(ii_coh+1)=these_cohs(kk);
                            ii_coh=ii_coh+1;
                        end
                    end
                    
                end
                
                all_cohs=[all_cohs these_coh];
                all_PLVs=[all_PLVs these_PLV];
                
                if evNo==2
                    if per_ii==1
                        %S+ Proficient
                        plot(these_coh,these_PLV,'o','MarkerEdgeColor','none','MarkerFaceColor',[158/255 31/255 99/255])
                    else
                        %S+ Naive
                        plot(these_coh,these_PLV,'o','MarkerEdgeColor','none','MarkerFaceColor',[238/255 111/255 179/255])
                    end
                else
                    if per_ii==1
                        %S- Proficient
                        plot(these_coh,these_PLV,'o','MarkerEdgeColor','none','MarkerFaceColor',[0 114/255 178/255])
                    else
                        %S- naive
                        plot(these_coh,these_PLV,'o','MarkerEdgeColor','none','MarkerFaceColor',[80/255 194/255 255/255])
                    end
                end
                
                
                
            end
            
            
        end
        
        title(['delta PLV vs. delta coherence for each odor pair/mouse for ' bandwidth_names{bwii} ' and ' group_legend{grNo}])
        
        
        %     %Annotations identifying groups
        %     x_interval=0.8/ii_gr_included;
        %     for ii=1:ii_gr_included
        %         annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
        %     end
        %
        %     %Proficient/Naive annotations
        %     annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
        %     annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
        
        
        
        
        ylabel('delta PLV')
        xlabel('delta coherence')
        xlim([-0.4 0.3])
        ylim([-0.6 0.5])
        
        ii_pval=ii_pval+1;
        x=all_cohs';
        y=all_PLVs';
        [rho,pval(ii_pval)]=corr(x,y);
        
        %Fit a line
        c = polyfit(x,y,1);
        % Display evaluated equation y = m*x + b
        disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
        % Evaluate fit equation using polyval
        y_est = polyval(c,[-0.4 0.3]);
        % Add trend line to plot
        plot([-0.4 0.3],y_est,'k-','LineWidth',2)
        
        fprintf(1, ['rho = %d, p value = %d for delta PLV vs delta coherence for '  bandwidth_names{bwii} ' and ' group_legend{grNo} '\n\n'],rho,pval(ii_pval))
        
        
    end
end


%Now plot the average delta phase calculated per odor pair

ii_pval=0;
pval=[];

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
    
    
    all_PLVs=[];
    all_phases=[];
    
    for evNo=1:2
        
        for per_ii=2:-1:1
            
            grNo=1;
            
            bar_offset = bar_offset + 1;
            
            %Get these delta phase values
            
            
            these_delta_phase=[];
            ii_delta_phase=0;
            these_PLV=[];
            ii_PLV=0;
            for ii=1:length(FileNamePLV)
                %                     ii_delta_phase=ii_delta_phase+1;
                %                     these_delta_phase(ii_delta_phase)=all_filesPLV(ii).delta_phase_odor(grNo,bwii,per_ii,evNo);
                these_delta_phases=zeros(1,sum(group_no_per_mouse==grNo));
                these_delta_phases(1,:)=all_filesPLV(ii).delta_phase_odor_per_mouse(group_no_per_mouse==grNo,bwii,per_ii,evNo);
                these_delta_phase(ii_delta_phase+1:ii_delta_phase+length(these_delta_phases))=these_delta_phases;
                ii_delta_phase=ii_delta_phase+length(these_delta_phases);
                
                these_PLVs=zeros(1,sum(group_no_per_mouse==grNo));
                these_PLVs(1,:)=all_filesPLV(ii).delta_PLV_odor_minus_ref_per_mouse(group_no_per_mouse==grNo,bwii,per_ii,evNo);
                these_PLV(ii_PLV+1:ii_PLV+length(these_PLVs))=these_PLVs;
                ii_PLV=ii_PLV+length(these_PLVs);
            end
            
            all_PLVs=[all_PLVs these_PLV];
            all_phases=[all_phases these_delta_phase];
            
            if evNo==2
                if per_ii==1
                    %S+ Proficient
                    plot(these_PLV,these_delta_phase,'o','MarkerEdgeColor','none','MarkerFaceColor',[158/255 31/255 99/255])
                else
                    %S+ Naive
                    plot(these_PLV,these_delta_phase,'o','MarkerEdgeColor','none','MarkerFaceColor',[238/255 111/255 179/255])
                end
            else
                if per_ii==1
                    %S- Proficient
                    plot(these_PLV,these_delta_phase,'o','MarkerEdgeColor','none','MarkerFaceColor',[0 114/255 178/255])
                else
                    %S- naive
                    plot(these_PLV,these_delta_phase,'o','MarkerEdgeColor','none','MarkerFaceColor',[80/255 194/255 255/255])
                end
            end
            
            
            
            
            
        end
        
        
    end
    
    title(['Average delta phase vs delta PLV for odor pair/mouse for ' bandwidth_names{bwii}])
    
    
    %     %Annotations identifying groups
    %     x_interval=0.8/ii_gr_included;
    %     for ii=1:ii_gr_included
    %         annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
    %     end
    %
    %     %Proficient/Naive annotations
    %     annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
    %     annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
    
    
    ylabel('delta phase')
    xlabel('delta PLV')
    xlim([-0.4 0.3])
    ylim([-1 1.5])
    
    ii_pval=ii_pval+1;
    x=all_phases';
    y=all_PLVs';
    [rho,pval(ii_pval)]=corr(x,y);
    
    %Fit a line
    c = polyfit(x,y,1);
    % Display evaluated equation y = m*x + b
    disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
    % Evaluate fit equation using polyval
    y_est = polyval(c,[-0.4 0.3]);
    % Add trend line to plot
    plot([-0.4 0.3],y_est,'k-','LineWidth',2)
    
    fprintf(1, ['rho = %d, p value = %d for delta phase vs delta PLV for '  bandwidth_names{bwii} '\n\n'],rho,pval(ii_pval))
    
    
    
end




pffft=1;
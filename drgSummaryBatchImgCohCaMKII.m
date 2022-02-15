function drgSummaryBatchImgCohCaMKII
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
hippPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Imaginary coherence/';
%
% %Files
% FileName{1}='CaMKIIacetocohe02012021_out.mat';
% FileName{2}='CaMKIIethylbenacetocohe2262021_out.mat';
% FileName{3}='CaMKIIPAEAcohe02082021_out.mat';
% FileName{4}='CaMKIIEAPAcohe2262021_out.mat';
% FileName{5}='CaMKIIpz1EAcohe02142021_out.mat';
% FileName{6}='CaMKIIPZ1PAEAcohe202102021_out.mat';
% FileName{7}='CaMKIIpzz1EAPAcohe02112021_out.mat';
% FileName{8}='CaMKIIpzz1propylacecohe02092021_out.mat';


FileName{1}='CaMKIIAPEBimgcohe01222022_out80.mat';
FileName{2}='spm_LFP_EBAP_imgcoh01282022_out80.mat';
FileName{3}='CaMKIIPAEAimagcohe01112022_out80';
FileName{4}='CaMKIIEAPAimagcohe01072022_out80.mat';
FileName{5}='CaMKIIpz1EAPAimagcohe01172022_out80.mat';
FileName{6}='CaMKIIpz1PAEAimagcohe01122022_out80.mat';
FileName{7}='spm_LFP_pzz1EAPA_imgcoh01252022_out80.mat';
FileName{8}='CaMKIIpzz1PAEAimgcoh01222021_out80.mat';

%Load the table of mouse numbers
%Note: This may need to be revised for PRP
mouse_no_table='/Users/restrepd/Documents/Projects/CaMKII_analysis/Reply_to_reviewers/camkii_mice_per_odor_pair_for_img_coh.xlsx';
T_mouse_no = readtable(mouse_no_table);


%Text file for statistical output
fileID = fopen([hippPathName 'drgSummaryBatchImgCohCaMKII.txt'],'w');

%Load data
all_files=[];

for ii=1:length(FileName)
    load([hippPathName FileName{ii}])
    all_files(ii).handles_out=handles_out;
end

handles_out=[];

figNo=0;

%Now plot the average PRP for each electrode calculated per mouse
%(including all sessions for each mouse)
edges=[-0.3:0.02:0.3];
rand_offset=0.5;


for bwii=[1 2 4]    %for amplitude bandwidths (beta, low gamma, high gamma)
    
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
                        if mouse_op==1
                            these_coh(ii_coh+1:ii_coh+length(all_files(ii).handles_out.dcoh_values(this_jj).dcoh_per_mouse))=all_files(ii).handles_out.dcoh_values(this_jj).dcoh_per_mouse;
                            ii_coh=ii_coh+length(all_files(ii).handles_out.dcoh_values(this_jj).dcoh_per_mouse);
                        else
                            ii_coh=ii_coh+1;
                            these_coh(ii_coh)=all_files(ii).handles_out.dcoh_values(this_jj).dcoh;
                        end
                    end
                end
                
                switch grNo
                    case 1
                        bar(bar_offset,mean(these_coh),'g','LineWidth', 3,'EdgeColor','none')
                    case 2
                        bar(bar_offset,mean(these_coh),'b','LineWidth', 3,'EdgeColor','none')
                    case 3
                        bar(bar_offset,mean(these_coh),'y','LineWidth', 3,'EdgeColor','none')
                end
                
                
                %Violin plot
                
                [mean_out, CIout]=drgViolinPoint(these_coh,edges,bar_offset,rand_offset,'k','k',3);
%                 CI = bootci(1000, {@mean, these_coh},'type','cper');
%                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%                 plot(bar_offset*ones(1,length(these_coh)),these_coh,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%                 
                
                %                                 %Save data for glm and ranksum
                
                glm_coh.data(glm_ii+1:glm_ii+length(these_coh))=these_coh;
                glm_coh.group(glm_ii+1:glm_ii+length(these_coh))=grNo*ones(1,length(these_coh));
                glm_coh.perCorr(glm_ii+1:glm_ii+length(these_coh))=per_ii*ones(1,length(these_coh));
                glm_coh.event(glm_ii+1:glm_ii+length(these_coh))=evNo*ones(1,length(these_coh));
                glm_ii=glm_ii+length(these_coh);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_coh;
                input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii} ' ' group_legend{grNo}];
                
                
            end
            bar_offset = bar_offset + 1;
            
        end
        bar_offset = bar_offset + 1;
        
    end
    
    title(['Average delta coherence per odor pair per mouse for ' bandwidth_names{bwii}])
    
    
    
    xticks([1 2 3 5 6 7 10 11 12 14 15 16])
    xticklabels({'nwtS+', 'nHS+', 'nKOS+', 'pwtS+', 'pHS+', 'pKOS+', 'nwtS-', 'nHS-', 'nKOS-', 'pwtS-', 'pHS-', 'pKOS-'})
    
    ylim([-0.15 0.15])
    ylabel('delta coherence')
    
    
    %Perform the glm
    fprintf(1, ['glm for average coherence per odor pair per mouse for '  bandwidth_names{bwii} '\n'])
    
    fprintf(1, ['\n\nglm for PRP for' bandwidth_names{bwii} '\n'])
    tbl = table(glm_coh.data',glm_coh.group',glm_coh.perCorr',glm_coh.event',...
        'VariableNames',{'MI','group','perCorr','event'});
    mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
        ,'CategoricalVars',[2,3,4])
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values per odor pair per mouse for ' bandwidth_names{bwii} ' hippocampus\n'])
    [output_data] = drgMutiRanksumorTtest(input_data);
    
    
end


%Now plot the percent significant delta Cxy
edges=[0:5:100];
rand_offset=0.5;


for bwii=[1 2 4]    %for amplitude bandwidths (beta, low gamma, high gamma)
    
    for grNo=1:3
        glm_coh=[];
        glm_ii=0;
        
        id_ii=0;
        input_data=[];
        
        glm_coh_mm=[];
        glm_ii_mm=0;
        
        id_ii_mm=0;
        input_data_mm=[];
        
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
        
        ax=gca;ax.LineWidth=3;
        
        bar_lab_loc=[];
        no_ev_labels=0;
        ii_gr_included=0;
        bar_offset = 0;
        
        for evNo=1:2
            
            %         these_coh_per_Ev=[];
            %         ii_coh=0;
            per_sig_thisEv=[];
            mouse_no_thisEv=[];
            per_ii_thisEv=[];
            
            for per_ii=2:-1:1
                
                
                
                bar_offset = bar_offset + 1;
                
                %Get these coherence values
                these_per_sig=[];
                ii_these_per_sig=0;
                for ii=1:length(FileName)
                    these_mouse_no_per_op=T_mouse_no.mouse_no_per_op(T_mouse_no.odor_pair_no==ii);
                    these_mouse_no=T_mouse_no.mouse_no(T_mouse_no.odor_pair_no==ii);
                    this_jj=[];
                    for jj=1:all_files(ii).handles_out.dcoh_ii
                        
                        if all_files(ii).handles_out.per_sig_values(jj).bwii==bwii
                            if all_files(ii).handles_out.per_sig_values(jj).evNo==evNo
                                if all_files(ii).handles_out.per_sig_values(jj).per_ii==per_ii
                                    
                                    
                                    if all_files(ii).handles_out.per_sig_values(jj).grNo==grNo
                                        this_jj=jj;
                                    end
                                    
                                end
                                
                            end
                        end
                    end
                    if ~isempty(this_jj)
                        if mouse_op==1
                            these_per_sig(ii_these_per_sig+1:ii_these_per_sig+length(all_files(ii).handles_out.per_sig_values(this_jj).percent_sig_per_mouse))=all_files(ii).handles_out.per_sig_values(this_jj).percent_sig_per_mouse;
                            ii_these_per_sig=ii_these_per_sig+length(all_files(ii).handles_out.per_sig_values(this_jj).percent_sig_per_mouse);
                        else
                            ii_these_per_sig=ii_these_per_sig+1;
                            these_per_sig(ii_these_per_sig)=all_files(ii).handles_out.per_sig_values(this_jj).dcoh;
                        end
                        
                        %Assign to overall mouse numbers
                        these_percent_sig=all_files(ii).handles_out.per_sig_values(this_jj).percent_sig_per_mouse;
                        these_mouse_nos=all_files(ii).handles_out.per_sig_values(this_jj).mouseNo;
                        these_mouse_no_overall=[];
                        for ii_mouse=1:length(these_mouse_nos)
                            these_mouse_no_overall(ii_mouse)=these_mouse_no(find(these_mouse_nos(ii_mouse)==these_mouse_no_per_op));
                        end
                        per_sig_thisEv=[per_sig_thisEv these_percent_sig];
                        mouse_no_thisEv=[mouse_no_thisEv these_mouse_no_overall];
                        per_ii_thisEv=[per_ii_thisEv per_ii*ones(1,length(these_percent_sig))];
                        
                    end
                end
                
                
                if evNo==2
                    if per_ii==1
                        %S- Proficient
                        bar(bar_offset,mean(these_per_sig),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                        %                     drgViolinPlot([bar_offset-1 bar_offset],[last_these_per_sig;these_per_sig],edges,rand_offset,'k','k',4,1);
                    else
                        %S- Naive
                        bar(bar_offset,mean(these_per_sig),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                        last_these_per_sig=these_per_sig;
                    end
                else
                    if per_ii==1
                        %S+ Proficient
                        bar(bar_offset,mean(these_per_sig),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                        %                     drgViolinPlot([bar_offset-1 bar_offset],[last_these_per_sig;these_per_sig],edges,rand_offset,'k','k',4,1);
                    else
                        %S+ naive
                        bar(bar_offset,mean(these_per_sig),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                        last_these_per_sig=these_per_sig;
                    end
                end
                
                %             ii_these_per_sig=ii_these_per_sig+1;
                %             these_per_sig_per_Ev(ii_these_per_sig,:)=these_per_sig;
                
                %Violin plot
                
                %             [mean_out, CIout]=drgViolinPoint(these_per_sig,edges,bar_offset,rand_offset,'k','k',3);
                %             CI = bootci(1000, {@mean, these_per_sig},'type','cper');
                %             plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                %             plot(bar_offset*ones(1,length(these_per_sig)),these_per_sig,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                %
                
                
                glm_coh.data(glm_ii+1:glm_ii+length(these_per_sig))=these_per_sig;
                %                 glm_coh.group(glm_ii+1:glm_ii+length(these_per_sig))=grNo*ones(1,length(these_per_sig));
                glm_coh.perCorr(glm_ii+1:glm_ii+length(these_per_sig))=per_ii*ones(1,length(these_per_sig));
                glm_coh.event(glm_ii+1:glm_ii+length(these_per_sig))=evNo*ones(1,length(these_per_sig));
                glm_coh.group(glm_ii+1:glm_ii+length(these_per_sig))=grNo*ones(1,length(these_per_sig));
                glm_ii=glm_ii+length(these_per_sig);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_per_sig;
                input_data(id_ii).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{per_ii} ' ' group_legend{grNo}];
                
            end
            %         drgViolinPlot([bar_offset-1 bar_offset],these_per_sig_per_Ev,edges,rand_offset,'k','k',4,1);
            %Calculate the per mouse means
            
            per_ii2_mean_per_sig_mm=[];
            per_ii1_mean_per_sig_mm=[];
            for mouseNo=1:max(mouse_no_thisEv)
                per_ii1_mean_per_sig_mm=[per_ii1_mean_per_sig_mm mean(per_sig_thisEv((mouse_no_thisEv==mouseNo)&(per_ii_thisEv==1)))];
                per_ii2_mean_per_sig_mm=[per_ii2_mean_per_sig_mm mean(per_sig_thisEv((mouse_no_thisEv==mouseNo)&(per_ii_thisEv==2)))];
            end
            
            %per_ii=1
            glm_coh_mm.data(glm_ii_mm+1:glm_ii_mm+length(per_ii1_mean_per_sig_mm))=per_ii1_mean_per_sig_mm;
            %                 glm_coh_mm.group(glm_ii_mm+1:glm_ii_mm+length(these_per_sig))=grNo*ones(1,length(per_ii1_mean_per_sig_mm));
            glm_coh_mm.perCorr(glm_ii_mm+1:glm_ii_mm+length(per_ii1_mean_per_sig_mm))=ones(1,length(per_ii1_mean_per_sig_mm));
            glm_coh_mm.group(glm_ii_mm+1:glm_ii_mm+length(per_ii1_mean_per_sig_mm))=grNo*ones(1,length(per_ii1_mean_per_sig_mm));
            glm_coh_mm.event(glm_ii_mm+1:glm_ii_mm+length(per_ii1_mean_per_sig_mm))=evNo*ones(1,length(per_ii1_mean_per_sig_mm));
            glm_ii_mm=glm_ii_mm+length(per_ii1_mean_per_sig_mm);
            
            id_ii_mm=id_ii_mm+1;
            input_data_mm(id_ii_mm).data=per_ii1_mean_per_sig_mm;
            input_data_mm(id_ii_mm).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{1} ' ' group_legend{grNo}];
            
            %per_ii=1
            glm_coh_mm.data(glm_ii_mm+1:glm_ii_mm+length(per_ii2_mean_per_sig_mm))=per_ii2_mean_per_sig_mm;
            %                 glm_coh_mm.group(glm_ii_mm+1:glm_ii_mm+length(these_per_sig))=grNo*ones(1,length(per_ii2_mean_per_sig_mm));
            glm_coh_mm.perCorr(glm_ii_mm+1:glm_ii_mm+length(per_ii2_mean_per_sig_mm))=2*ones(1,length(per_ii2_mean_per_sig_mm));
            glm_coh_mm.group(glm_ii_mm+1:glm_ii_mm+length(per_ii2_mean_per_sig_mm))=grNo*ones(1,length(per_ii2_mean_per_sig_mm));
            glm_coh_mm.event(glm_ii_mm+1:glm_ii_mm+length(per_ii2_mean_per_sig_mm))=evNo*ones(1,length(per_ii2_mean_per_sig_mm));
            glm_ii_mm=glm_ii_mm+length(per_ii2_mean_per_sig_mm);
            
            id_ii_mm=id_ii_mm+1;
            input_data_mm(id_ii_mm).data=per_ii2_mean_per_sig_mm;
            input_data_mm(id_ii_mm).description=[ evTypeLabels{evNo} ' ' prof_naive_leg{2} ' ' group_legend{grNo}];
            
            for mouseNo=1:length(per_ii2_mean_per_sig_mm)
                plot([bar_offset-1 bar_offset],[per_ii2_mean_per_sig_mm(mouseNo) per_ii1_mean_per_sig_mm(mouseNo)],'-ok','MarkerSize',6,'LineWidth',1,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor',[0.6 0.6 0.6],'Color',[0.6 0.6 0.6])
            end
            
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_per_sig...
                ,edges,bar_offset,rand_offset,'k','k',3);
            
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(last_these_per_sig...
                ,edges,bar_offset-1,rand_offset,'k','k',3);
            
            bar_offset = bar_offset + 1;
            
        end
        
        title(['Percent significant delta coherence for ' bandwidth_names{bwii} ' ' group_legend{grNo}])
        ylim([0 100])
        
        %     %Annotations identifying groups
        %     x_interval=0.8/ii_gr_included;
        %     for ii=1:ii_gr_included
        %         annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
        %     end
        %
        %     %Proficient/Naive annotations
        %     annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
        %     annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
        
        plot([0 6],[0 0],'-k','LineWidth',2)
        xlim([0 6])
        
        xticks([1 2 4 5])
        xticklabels({'nS+', 'pS+','nS-', 'pS-'})
        
        ylabel('Percent significant')
        
    end
    %Perform the glm per mouse per odor pair
    fprintf(1, ['glm for percent significant delta imaginary coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n'])
    fprintf(fileID, ['glm for percent significant delta imaginary coherence per mouse per odor pair for '  bandwidth_names{bwii} '\n']);
    
    tbl = table(glm_coh.data',glm_coh.perCorr',glm_coh.event',glm_coh.group',...
        'VariableNames',{'per_sig','naive_vs_proficient','sp_vs_sm','genotype'});
    mdl = fitglm(tbl,'per_sig~naive_vs_proficient+sp_vs_sm+genotype+naive_vs_proficient*sp_vs_sm*genotype'...
        ,'CategoricalVars',[2,3,4])
    
    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
    fprintf(fileID,'%s\n', txt);
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for percent significant delta imaginary coherence per mouse per odor pair for ' bandwidth_names{bwii} ' hippocampus\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for percent significant delta imaginary coherence per mouse per odor pair for ' bandwidth_names{bwii} ' hippocampus\n']);
    [output_data] = drgMutiRanksumorTtest(input_data, fileID);
    
    %Perform the glm per mouse
    fprintf(1, ['glm for percent significant delta imaginary coherence per mouse for '  bandwidth_names{bwii} '\n'])
    fprintf(fileID, ['glm for percent significant delta imaginary coherence per mouse for '  bandwidth_names{bwii} '\n']);
    
    tbl = table(glm_coh_mm.data',glm_coh_mm.perCorr',glm_coh_mm.event',glm_coh_mm.group',...
        'VariableNames',{'per_sig','naive_vs_proficient','sp_vs_sm','genotype'});
    mdl = fitglm(tbl,'per_sig~naive_vs_proficient+sp_vs_sm+genotype+naive_vs_proficient*sp_vs_sm*genotype'...
        ,'CategoricalVars',[2,3,4])
    
    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
    fprintf(fileID,'%s\n', txt);
    
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for percent significant delta imaginary coherence per mouse for ' bandwidth_names{bwii} ' hippocampus\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for percent significant delta imaginary coherence per mouse for ' bandwidth_names{bwii} ' hippocampus\n']);
    [output_data] = drgMutiRanksumorTtest(input_data_mm, fileID);
    
end

fclose(fileID);


pffft=1;
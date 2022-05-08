function drgSummaryBatchPRPCaMKIIWT
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'drgbDiscPar' file listing 'Discriminant_*.mat' output files
%from drgLFPDiscriminantBatch
%
%Performs summary analises for LDA and  PAC analysis performed with case 19
%of drgAnalysisBatchLFP

warning('off')

close all hidden
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
% hippPathName='F:\Datos summary CaMKII111720\PRP drgAnalysisBatchLFPCaMKII case 24 output for summary\';
hippPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/PRP drgAnalysisBatchLFPCaMKII case 24 output for 80/';


%Text file for statistical output
fileID = fopen([hippPathName 'drgSummaryBatchPRPCaMKIIWTstats.txt'],'w');

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

%Load the table of mouse numbers
mouse_no_table='/Users/restrepd/Documents/Projects/CaMKII_analysis/Reply_to_reviewers/camkii_mice_per_odor_pair_for_PRP.xlsx';
T_mouse_no = readtable(mouse_no_table);

%Now process the hippocampus

%Load data
all_hippo=[];

for ii=1:length(hippFileName)
    load([hippPathName hippFileName{ii}])
    all_hippo(ii).handles_out=handles_out;
end

all_pre=[];

for ii=1:length(preFileName)
    load([prePathName preFileName{ii}])
    all_pre(ii).handles_out=handles_out;
end

figNo=0;

%Now plot the average PRP for each electrode calculated per mouse
%(including all sessions for each mouse)
edges=[-15:0.5:15];
rand_offset=0.5;

for pacii=[1 3]    %for amplitude bandwidths (beta, low gamma, high gamma)
    glm_PRP_hipp=[];
    glm_ii_hipp=0;
    
    glm_PRP_both=[];
    glm_PRP_ii_both=0;
    
    glm_PRP_pre=[];
    glm_ii_pre=0;
    
    for peak=0:1
        
        
        
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
        
        ax=gca;ax.LineWidth=3;
        
        bar_lab_loc=[];
        no_ev_labels=0;
        ii_gr_included=0;
        bar_offset = 0;
        
        for evNo=1:2
            per_ii_PRP=[];
            unique_mouse_nos=[];
            for per_ii=2:-1:1
                
                grNo=1;
                bar_offset = bar_offset +1;
                
                
                
                %Get these PRP values
                these_PRP=[];
                ii_PRP=0;
                for ii=1:length(hippFileName)
                    this_jj=[];
                    these_mouse_no_per_op=T_mouse_no.mouse_no_per_op(T_mouse_no.odor_pair_no==ii);
                    these_mouse_no=T_mouse_no.mouse_no(T_mouse_no.odor_pair_no==ii);
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
                        these_PRP(ii_PRP+1:ii_PRP+length(all_hippo(ii).handles_out.PRP_values(this_jj).PRP_per_mouse))=all_hippo(ii).handles_out.PRP_values(this_jj).PRP_per_mouse;
                        
                        for ww=1:length(all_hippo(ii).handles_out.PRP_values(this_jj).PRP_per_mouse)
                            this_nn=find(all_hippo(ii).handles_out.PRP_values(this_jj).mouseNo(ww)==these_mouse_no_per_op);
                            these_mouse_nos(ii_PRP+ww)=these_mouse_no(this_nn);
                        end
                        
                        ii_PRP=ii_PRP+length(all_hippo(ii).handles_out.PRP_values(this_jj).PRP_per_mouse);
                    end
                end
                
                per_ii_PRP(per_ii).these_PRP=these_PRP;
                per_ii_PRP(per_ii).these_mouse_nos=these_mouse_nos;
                unique_mouse_nos=unique([unique_mouse_nos these_mouse_nos]);
                
                if evNo==2
                    if per_ii==1
                        %S+ Proficient
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    else
                        %S+ Naive
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                    end
                else
                    if per_ii==1
                        %S- Proficient
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                    else
                        %S- naive
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                    end
                end
                
                %Mean per mouse
                if per_ii==1
                    for msNo=unique_mouse_nos
                        
                        this_per_ii=1;
                        this_mouse_mean1=mean(per_ii_PRP(this_per_ii).these_PRP(per_ii_PRP(this_per_ii).these_mouse_nos==msNo));
                        
                        this_per_ii=2;
                        this_mouse_mean2=mean(per_ii_PRP(this_per_ii).these_PRP(per_ii_PRP(this_per_ii).these_mouse_nos==msNo));
                        
                        plot([bar_offset-1 bar_offset],[this_mouse_mean2 this_mouse_mean1],'-ok','MarkerSize',6,'LineWidth',1,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor',[0.6 0.6 0.6],'Color',[0.6 0.6 0.6])
                        
                    end
                end
                
                
                %Violin plot
                
                [mean_out, CIout]=drgViolinPoint(these_PRP,edges,bar_offset,rand_offset,'k','k',3);
                %                 CI = bootci(1000, {@mean, these_PRP},'type','cper');
                %                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                %                 plot(bar_offset*ones(1,length(these_PRP)),these_PRP,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                %
                
                %                                 %Save data for glm and ranksum
                
                glm_PRP_hipp.data(glm_ii_hipp+1:glm_ii_hipp+length(these_PRP))=these_PRP;
                glm_PRP_hipp.perCorr(glm_ii_hipp+1:glm_ii_hipp+length(these_PRP))=per_ii*ones(1,length(these_PRP));
                glm_PRP_hipp.event(glm_ii_hipp+1:glm_ii_hipp+length(these_PRP))=evNo*ones(1,length(these_PRP));
                glm_PRP_hipp.peak(glm_ii_hipp+1:glm_ii_hipp+length(these_PRP))=peak*ones(1,length(these_PRP));
                glm_PRP_hipp.mouseNo(glm_ii_hipp+1:glm_ii_hipp+length(these_PRP))=these_mouse_nos;
                glm_ii_hipp=glm_ii_hipp+length(these_PRP);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_PRP;
                input_data(id_ii).description=[evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                
                
                glm_PRP_both.data(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=these_PRP;
                glm_PRP_both.brain_region(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=ones(1,length(these_PRP));
                glm_PRP_both.perCorr(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=per_ii*ones(1,length(these_PRP));
                glm_PRP_both.event(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=evNo*ones(1,length(these_PRP));
                glm_PRP_both.peak(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=peak*ones(1,length(these_PRP));
                glm_PRP_ii_both=glm_PRP_ii_both+length(these_PRP);
                
                
                
                
            end
            bar_offset = bar_offset + 1;
            
        end
        
        title(['Average PRP per odor pair per mouse for ' peak_label{peak+1} ' theta/' PACnames{pacii} ' hippocampus'])
        
        
        %         %Annotations identifying groups
        %         x_interval=0.8/ii_gr_included;
        %         for ii=1:ii_gr_included
        %             annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
        %         end
        
        %         %Proficient/Naive annotations
        %         annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
        %         annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
        %
        %
        
        
        plot([0 6],[0 0],'-k','LineWidth',2)
        xlim([0 6])
        
        xticks([1 2 4 5])
        xticklabels({'nS+', 'pS+','nS-', 'pS-'})
        
        ylim([-15 15])
        
        ylabel('PRP')
        
        
        
        
        
        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for average tPRP for each electrode calculated per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' hippocampus\n'])
        fprintf(fileID, ['\n\nRanksum or t-test p values for average tPRP for each electrode calculated per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' hippocampus\n']);
        [output_data] = drgMutiRanksumorTtest(input_data,fileID);
        
        
    end
    
    %Perform the glm
    fprintf(1, ['glm for average tPRP per odor pair per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' hippocampus\n'])
    fprintf(fileID, ['glm for average tPRP per odor pair per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' hippocampus\n']);
    
    tbl = table(glm_PRP_hipp.data',glm_PRP_hipp.perCorr',glm_PRP_hipp.event',glm_PRP_hipp.peak',...
        'VariableNames',{'tPRP','naive_vs_proficient','rewarded_vs_unrewarded','peak_trough'});
    mdl = fitglm(tbl,'tPRP~naive_vs_proficient+rewarded_vs_unrewarded+peak_trough+naive_vs_proficient*rewarded_vs_unrewarded*peak_trough'...
        ,'CategoricalVars',[2,3,4])
    
    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
    fprintf(fileID,'%s\n', txt);
    
     %Nested ANOVAN
    %https://www.mathworks.com/matlabcentral/answers/491365-within-between-subjects-in-anovan
    nesting=[0 0 0 0; ... % This line indicates that group factor is not nested in any other factor.
         0 0 0 0; ... % This line indicates that perCorr is not nested in any other factor.
         0 0 0 0; ... % This line indicates that event is not nested in any other factor.
         1 1 1 0];    % This line indicates that mouse_no (the last factor) is nested under group, perCorr and event
                    % (the 1 in position 1 on the line indicates nesting under the first factor).
    figNo=figNo+1;
                     
    [p anovanTbl stats]=anovan(glm_PRP_hipp.data,{glm_PRP_hipp.perCorr glm_PRP_hipp.event glm_PRP_hipp.peak glm_PRP_hipp.mouseNo},...
    'model','interaction',...
    'nested',nesting,...
    'varnames',{'naive_vs_proficient', 'rewarded_vs_unrewarded','peak_vs_trough','mouse_no'});
  
    fprintf(fileID, ['\n\nNested ANOVAN for tPRP during odor per mouse per odor pair for '  PACnames{pacii} ' for hippocampus\n']);
    drgWriteANOVANtbl(anovanTbl,fileID);
    fprintf(fileID, '\n\n');
    
    
    %Now process the prefrontal
    for peak=0:1
        
        
        
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
        
        ax=gca;ax.LineWidth=3;
        
        bar_lab_loc=[];
        no_ev_labels=0;
        ii_gr_included=0;
        bar_offset = 0;
        
        for evNo=1:2
             per_ii_PRP=[];
            unique_mouse_nos=[]
            for per_ii=2:-1:1
                
                grNo=1;
                bar_offset = bar_offset +1;
                
                
                
                %Get these PRP values
                these_PRP=[];
                ii_PRP=0;
                for ii=1:length(hippFileName)
                    this_jj=[];
                     these_mouse_no_per_op=T_mouse_no.mouse_no_per_op(T_mouse_no.odor_pair_no==ii);
                    these_mouse_no=T_mouse_no.mouse_no(T_mouse_no.odor_pair_no==ii);
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
                        
                        these_PRP(ii_PRP+1:ii_PRP+length(all_pre(ii).handles_out.PRP_values(this_jj).PRP_per_mouse))=all_pre(ii).handles_out.PRP_values(this_jj).PRP_per_mouse;
                        
                        for ww=1:length(all_hippo(ii).handles_out.PRP_values(this_jj).PRP_per_mouse)
                            this_nn=find(all_hippo(ii).handles_out.PRP_values(this_jj).mouseNo(ww)==these_mouse_no_per_op);
                            these_mouse_nos(ii_PRP+ww)=these_mouse_no(this_nn);
                        end
                        
                        ii_PRP=ii_PRP+length(all_pre(ii).handles_out.PRP_values(this_jj).PRP_per_mouse);
                        
                    end
                end
                
                per_ii_PRP(per_ii).these_PRP=these_PRP;
                per_ii_PRP(per_ii).these_mouse_nos=these_mouse_nos;
                unique_mouse_nos=unique([unique_mouse_nos these_mouse_nos]);
                
                if evNo==2
                    if per_ii==1
                        %S+ Proficient
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    else
                        %S- Naive
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                    end
                else
                    if per_ii==1
                        %S- Proficient
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                    else
                        %S- naive
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                    end
                end
                
                %Mean per mouse
                if per_ii==1
                    for msNo=unique_mouse_nos
                        
                        this_per_ii=1;
                        this_mouse_mean1=mean(per_ii_PRP(this_per_ii).these_PRP(per_ii_PRP(this_per_ii).these_mouse_nos==msNo));
                        
                        this_per_ii=2;
                        this_mouse_mean2=mean(per_ii_PRP(this_per_ii).these_PRP(per_ii_PRP(this_per_ii).these_mouse_nos==msNo));
                        
                        plot([bar_offset-1 bar_offset],[this_mouse_mean2 this_mouse_mean1],'-ok','MarkerSize',6,'LineWidth',1,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor',[0.6 0.6 0.6],'Color',[0.6 0.6 0.6])
                        
                    end
                end
                
                %Violin plot
                
                [mean_out, CIout]=drgViolinPoint(these_PRP,edges,bar_offset,rand_offset,'k','k',3);
%                 CI = bootci(1000, {@mean, these_PRP},'type','cper');
%                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%                 plot(bar_offset*ones(1,length(these_PRP)),these_PRP,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
%                 
                
                %                                 %Save data for glm and ranksum
                
                glm_PRP_pre.data(glm_ii_pre+1:glm_ii_pre+length(these_PRP))=these_PRP;
                glm_PRP_pre.perCorr(glm_ii_pre+1:glm_ii_pre+length(these_PRP))=per_ii*ones(1,length(these_PRP));
                glm_PRP_pre.event(glm_ii_pre+1:glm_ii_pre+length(these_PRP))=evNo*ones(1,length(these_PRP));
                glm_PRP_pre.peak(glm_ii_pre+1:glm_ii_pre+length(these_PRP))=peak*ones(1,length(these_PRP));
                glm_PRP_pre.mouseNo(glm_ii_pre+1:glm_ii_pre+length(these_PRP))=these_mouse_nos;
                glm_ii_pre=glm_ii_pre+length(these_PRP);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_PRP;
                input_data(id_ii).description=[evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                
                
                glm_PRP_both.data(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=these_PRP;
                glm_PRP_both.brain_region(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=2*ones(1,length(these_PRP));
                glm_PRP_both.perCorr(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=per_ii*ones(1,length(these_PRP));
                glm_PRP_both.event(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=evNo*ones(1,length(these_PRP));
                glm_PRP_both.peak(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=peak*ones(1,length(these_PRP));
                glm_PRP_ii_both=glm_PRP_ii_both+length(these_PRP);
                
            end
            bar_offset = bar_offset + 1;
            
        end
        
        title(['Average PRP per odor pair per mouse for ' peak_label{peak+1} ' theta/' PACnames{pacii} ' prefrontal'])
        
        
        %         %Annotations identifying groups
        %         x_interval=0.8/ii_gr_included;
        %         for ii=1:ii_gr_included
        %             annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
        %         end
        %
        %         %Proficient/Naive annotations
        %         annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
        %         annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
        
        plot([0 6],[0 0],'-k','LineWidth',2)
        xlim([0 6])
    
        xticks([1 2 4 5])
        xticklabels({'nS+', 'pS+','nS-', 'pS-'})
        
        ylabel('PRP')
        ylim([-15 15])

        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for average PRP for each electrode calculated per mouse for ' peak_label{peak+1} ' PAC theta' PACnames{pacii} ' prefrontal\n'])
        fprintf(fileID, ['\n\nRanksum or t-test p values for average PRP for each electrode calculated per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' hippocampus\n']);
        [output_data] = drgMutiRanksumorTtest(input_data,fileID);

    end
    
    %Perform the glm
    fprintf(1, ['glm for average PRP per odor pair per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' prefrontal\n'])
    fprintf(fileID, ['glm for average PRP per odor pair per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' prefrontal\n']);
    
    fprintf(1, ['\n\nglm for PRP for Theta/' PACnames{pacii} '\n'])
    tbl = table(glm_PRP_pre.data',glm_PRP_pre.perCorr',glm_PRP_pre.event',glm_PRP_pre.peak',...
        'VariableNames',{'tPRP','naive_vs_proficient','rewarded_vs_unrewarded','peak_trough'});
    mdl = fitglm(tbl,'tPRP~naive_vs_proficient+rewarded_vs_unrewarded+peak_trough+naive_vs_proficient*rewarded_vs_unrewarded*peak_trough'...
        ,'CategoricalVars',[2,3,4])
    
    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
    fprintf(fileID,'%s\n', txt);
    
    %Nested ANOVAN
    %https://www.mathworks.com/matlabcentral/answers/491365-within-between-subjects-in-anovan
    nesting=[0 0 0 0; ... % This line indicates that group factor is not nested in any other factor.
         0 0 0 0; ... % This line indicates that perCorr is not nested in any other factor.
         0 0 0 0; ... % This line indicates that event is not nested in any other factor.
         1 1 1 0];    % This line indicates that mouse_no (the last factor) is nested under group, perCorr and event
                    % (the 1 in position 1 on the line indicates nesting under the first factor).
    figNo=figNo+1;
                     
    [p anovanTbl stats]=anovan(glm_PRP_pre.data,{glm_PRP_pre.perCorr glm_PRP_pre.event glm_PRP_pre.peak glm_PRP_pre.mouseNo},...
    'model','interaction',...
    'nested',nesting,...
    'varnames',{'naive_vs_proficient', 'rewarded_vs_unrewarded','peak_vs_trough','mouse_no'});
  
    fprintf(fileID, ['\n\nNested ANOVAN for tPRP during odor per mouse per odor pair for '  PACnames{pacii} ' for prefrontal\n']);
    drgWriteANOVANtbl(anovanTbl,fileID);
    fprintf(fileID, '\n\n');
    
    
    %Perform the glm for both brain regions
    fprintf(1, ['glm for average PRP per odor pair per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' prefrontal and hippocampus\n'])
    
    fprintf(1, ['\n\nglm for PRP for Theta/' PACnames{pacii} '\n'])
    tbl = table(glm_PRP_both.data',glm_PRP_both.brain_region',glm_PRP_both.perCorr',glm_PRP_both.event',glm_PRP_both.peak',...
        'VariableNames',{'PRP','brain_region','naive_vs_proficient','sp_vs_sm','peak_trough'});
    mdl = fitglm(tbl,'PRP~brain_region+naive_vs_proficient+sp_vs_sm+peak_trough+naive_vs_proficient*sp_vs_sm*peak_trough*brain_region'...
        ,'CategoricalVars',[2,3,4,5])
    
    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
%     fprintf(fileID,'%s\n', txt);

end

fclose(fileID);

pffft=1;
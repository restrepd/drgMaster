function drgSummaryBatchPRPCaMKIIWTrev
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'drgbDiscPar' file listing 'Discriminant_*.mat' output files
%from drgLFPDiscriminantBatch
%
%Performs summary analises for LDA and  PAC analysis performed with case 19
%of drgAnalysisBatchLFP

warning('off')

close all
clear all

grNo=1;

PACnames{1}='Beta';
PACnames{2}='Low gamma';
PACnames{3}='High gamma';

prof_naive_leg{1}='Proficient';
prof_naive_leg{2}='Naive';

fwd_rev_leg{1}='Forward';
fwd_rev_leg{2}='Reverse';

group_legend{1}='WT';
group_legend{2}='Het';
group_legend{3}='KO';

evTypeLabels{1}='S+';
evTypeLabels{2}='S-';

peak_label{1}='Trough';
peak_label{2}='Peak';

%Location of files
% hippPathName='F:\Datos summary CaMKII111720\PRP drgAnalysisBatchLFPCaMKII case 24 output for summary\';
hippPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/PRP drgAnalysisBatchLFPCaMKII case 24 output for summary/';
fwd_rev_hippo=[1 2 2 1 1 2 1 2]; %1=forward, 2=reverse
    
%Hippocampus
hippFileName{1}='spm_LFP_acetowavephasepower32620_hippoPRPnew.mat';
hippFileName{2}='spm_LFP_acetowavephasepower32620_hippoPRPnewpropyl.mat';
hippFileName{3}='spm_LFP_ethylbenwavephasepower41420_hippoPRPnew.mat';
hippFileName{4}='spm_LFP_ethylwavephasepower3520_hippoPRPnew.mat';
hippFileName{5}='spm_LFP_pz1ethyllwavephasepower0213020_hippoPRPnew.mat';
hippFileName{6}='spm_LFP_pz1propylwavephasepower013020_hippoPRPnew.mat';
hippFileName{7}='spm_LFP_pzz1ethyllwavephasepower043020_hippoPRPnew.mat';
hippFileName{8}='spm_LFP_pzz1propylwavephasepower071220_hippoPRPnew.mat';

% prePathName='F:\Datos summary CaMKII111720\PRP drgAnalysisBatchLFPCaMKII case 24 output for summary\';
prePathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/PRP drgAnalysisBatchLFPCaMKII case 24 output for summary/';
fwd_rev_pre=[1 2 2 1 1 2 1 2]; %1=forward, 2=reverse

%Prefrontal
preFileName{1}='spm_LFP_acetowavephasepower32620_prefrontPRPnew.mat';
preFileName{2}='spm_LFP_acetowavephasepower32620_prefrontPRPnewpropyl.mat';
preFileName{3}='spm_LFP_ethylbenwavephasepower41420_prefrontPRPnew.mat';
preFileName{4}='spm_LFP_ethylwavephasepower3520_prefrontPRPnew.mat';
preFileName{5}='spm_LFP_pz1ethyllwavephasepower0213020_prefrontPRPnew.mat';
preFileName{6}='spm_LFP_pz1propylwavephasepower013020_prefrontPRPnew.mat';
preFileName{7}='spm_LFP_pzz1ethyllwavephasepower043020_prefrontPRPnew.mat';
preFileName{8}='spm_LFP_pzz1propylwavephasepower071220_prefrontPRPnew.mat';


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
edges=[0:0.001:0.02];
rand_offset=0.8;


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
            
            for fwd_rev=1:2
                per_ii=1;
                
                
                bar_offset = bar_offset +1;
                
                
                
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
                                            if fwd_rev_hippo(ii)==fwd_rev
                                                this_jj=jj;
                                            end
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
                
                
                if evNo==2
                    if fwd_rev==1
                        %S+ Forward
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    else
                        %S+ Reverse
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                    end
                else
                    if fwd_rev==1
                        %S- Forward
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                    else
                        %S- Reverse
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                    end
                end
                
                
                %Violin plot
                
                %[mean_out, CIout]=drgViolinPoint(these_PRP,edges,bar_offset,rand_offset,'k','k',1);
                CI = bootci(1000, {@mean, these_PRP},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(these_PRP)),these_PRP,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                
                %                                 %Save data for glm and ranksum
                
                glm_PRP_hipp.data(glm_ii_hipp+1:glm_ii_hipp+length(these_PRP))=these_PRP;
                glm_PRP_hipp.fwd_rev(glm_ii_hipp+1:glm_ii_hipp+length(these_PRP))=fwd_rev*ones(1,length(these_PRP));
                glm_PRP_hipp.event(glm_ii_hipp+1:glm_ii_hipp+length(these_PRP))=evNo*ones(1,length(these_PRP));
                glm_PRP_hipp.peak(glm_ii_hipp+1:glm_ii_hipp+length(these_PRP))=peak*ones(1,length(these_PRP));
                glm_ii_hipp=glm_ii_hipp+length(these_PRP);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_PRP;
                input_data(id_ii).description=[evTypeLabels{evNo} ' ' fwd_rev_leg{per_ii}];
                
                
                glm_PRP_both.data(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=these_PRP;
                glm_PRP_both.brain_region(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=ones(1,length(these_PRP));
                glm_PRP_both.fwd_rev(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=fwd_rev*ones(1,length(these_PRP));
                glm_PRP_both.event(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=evNo*ones(1,length(these_PRP));
                glm_PRP_both.peak(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=peak*ones(1,length(these_PRP));
                glm_PRP_ii_both=glm_PRP_ii_both+length(these_PRP);
                
                
                
                
            end
            bar_offset = bar_offset + 1;
            
        end
        
        title(['Average PRP for each electrode calculated per mouse for ' peak_label{peak+1} ' theta/' PACnames{pacii} ' hippocampus'])
        
        
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
        xticks([1 2 4 5])
        xticklabels({'FS+', 'RS+','FS-', 'RS-'})
        
        ylim([-13 10])
        
        ylabel('PRP')
        
        
        
        
        
        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for average PRP for each electrode calculated per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' hippocampus\n'])
        [output_data] = drgMutiRanksumorTtest(input_data);
        
        
    end
    
    %Perform the glm
    fprintf(1, ['glm for average PRP for each electrode calculated per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' hippocampus\n'])
    
    fprintf(1, ['\n\nglm for PRP for Theta/' PACnames{pacii} '\n'])
    tbl = table(glm_PRP_hipp.data',glm_PRP_hipp.fwd_rev',glm_PRP_hipp.event',glm_PRP_hipp.peak',...
        'VariableNames',{'PRP','fwd_rev','event','peak_trough'});
    mdl = fitglm(tbl,'PRP~fwd_rev+event+peak_trough+fwd_rev*event*peak_trough'...
        ,'CategoricalVars',[2,3,4])
    
    
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
            
            for fwd_rev=1:2
                per_ii=1;
                
%                 grNo=1;
                bar_offset = bar_offset +1;
                
                
                
                %Get these PRP values
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
                                            if fwd_rev_hippo(ii)==fwd_rev
                                                this_jj=jj;
                                            end
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
                
                
                if evNo==2
                    if fwd_rev==1
                        %S+ Fwd
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                    else
                        %S+ Rev
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                    end
                else
                    if fwd_rev==1
                        %S- Fwd
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                    else
                        %S- Rev
                        bar(bar_offset,mean(these_PRP),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                    end
                end
                
                
                %Violin plot
                
                %[mean_out, CIout]=drgViolinPoint(these_PRP,edges,bar_offset,rand_offset,'k','k',1);
                CI = bootci(1000, {@mean, these_PRP},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                plot(bar_offset*ones(1,length(these_PRP)),these_PRP,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                
                %                                 %Save data for glm and ranksum
                
                glm_PRP_pre.data(glm_ii_pre+1:glm_ii_pre+length(these_PRP))=these_PRP;
                glm_PRP_pre.fwd_rev(glm_ii_pre+1:glm_ii_pre+length(these_PRP))=fwd_rev*ones(1,length(these_PRP));
                glm_PRP_pre.event(glm_ii_pre+1:glm_ii_pre+length(these_PRP))=evNo*ones(1,length(these_PRP));
                glm_PRP_pre.peak(glm_ii_pre+1:glm_ii_pre+length(these_PRP))=peak*ones(1,length(these_PRP));
                glm_ii_pre=glm_ii_pre+length(these_PRP);
                
                id_ii=id_ii+1;
                input_data(id_ii).data=these_PRP;
                input_data(id_ii).description=[evTypeLabels{evNo} ' ' fwd_rev_leg{per_ii}];
                
                
                glm_PRP_both.data(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=these_PRP;
                glm_PRP_both.brain_region(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=2*ones(1,length(these_PRP));
                glm_PRP_both.fwd_rev(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=fwd_rev*ones(1,length(these_PRP));
                glm_PRP_both.event(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=evNo*ones(1,length(these_PRP));
                glm_PRP_both.peak(glm_PRP_ii_both+1:glm_PRP_ii_both+length(these_PRP))=peak*ones(1,length(these_PRP));
                glm_PRP_ii_both=glm_PRP_ii_both+length(these_PRP);
                
            end
            bar_offset = bar_offset + 1;
            
        end
        
        title(['Average PRP for each electrode calculated per mouse for ' peak_label{peak+1} ' theta/' PACnames{pacii} ' prefrontal'])
        
        
        %         %Annotations identifying groups
        %         x_interval=0.8/ii_gr_included;
        %         for ii=1:ii_gr_included
        %             annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
        %         end
        %
        %         %Proficient/Naive annotations
        %         annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
        %         annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
        
        
        xticks([1 2 4 5])
        xticklabels({'FS+', 'RS+','FS-', 'RS-'})
        
        ylabel('PRP')
        ylim([-13 10])
        
        
        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for average PRP for each electrode calculated per mouse for ' peak_label{peak+1} ' PAC theta' PACnames{pacii} ' prefrontal\n'])
        [output_data] = drgMutiRanksumorTtest(input_data);
        
        
    end
    
    %Perform the glm
    fprintf(1, ['glm for average PRP for each electrode calculated per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' prefrontal\n'])
    
    fprintf(1, ['\n\nglm for PRP for Theta/' PACnames{pacii} '\n'])
    tbl = table(glm_PRP_pre.data',glm_PRP_pre.fwd_rev',glm_PRP_pre.event',glm_PRP_pre.peak',...
        'VariableNames',{'PRP','fwd_rev','event','peak_trough'});
    mdl = fitglm(tbl,'PRP~fwd_rev+event+peak_trough+fwd_rev*event*peak_trough'...
        ,'CategoricalVars',[2,3,4])
    
    %Perform the glm for both brain regions
    fprintf(1, ['glm for average PRP for each electrode calculated per mouse for ' peak_label{peak+1} ' theta' PACnames{pacii} ' prefrontal and hippocampus\n'])
    
    fprintf(1, ['\n\nglm for PRP for Theta/' PACnames{pacii} '\n'])
    tbl = table(glm_PRP_both.data',glm_PRP_both.brain_region',glm_PRP_both.fwd_rev',glm_PRP_both.event',glm_PRP_both.peak',...
        'VariableNames',{'PRP','brain_region','fwd_rev','event','peak_trough'});
    mdl = fitglm(tbl,'PRP~brain_region+fwd_rev+event+peak_trough+fwd_rev*event*peak_trough*brain_region'...
        ,'CategoricalVars',[2,3,4,5])
    
end

pffft=1;
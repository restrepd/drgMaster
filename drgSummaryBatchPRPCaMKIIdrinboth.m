function drgSummaryBatchPRPCaMKIIdrinboth
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
PACnames{2}='High gamma';
PACnames{3}='SWR';

prof_naive_leg{1}='Proficient';
prof_naive_leg{2}='Naive';

group_legend{1}='tatCN21_exp';
group_legend{2}='tatCN21_cont';

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
hippPathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/Drug infusion/';
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
hippFileName{1}='spm_LFP_druginfusionaceto_one_window_31522_case1__hippoinfLFP80.mat';
hippFileName{2}='spm_LFP_druginfusionpropylace_one_window_32222_case1__hippoinfLFP80.mat';
hippFileName{3}='spm_LFP_druginfusionaceto_one_window_31522_case1__hippoinfLFP80.mat';
hippFileName{4}='spm_LFP_druginfusionethylben_one_window_31522_case1__hippoinfLFP80.mat';
% hippFileName{3}='spm_LFP_ethylwavephasepower3520_hippocampusLFP80.mat';
% hippFileName{4}='spm_LFP_acetowavephasepower32620 2_hippocampusLFP80.mat';
% hippFileName{5}='spm_LFP_pz1ethyllwavephasepower0213020_hippocampusLFP80.mat';
% hippFileName{6}='spm_LFP_pz1propylwavephasepower013020_hippocampusLFP80.mat';
% hippFileName{7}='spm_LFP_pzz1ethyllwavephasepower043020_hippocampusLFP80.mat';
% hippFileName{8}='spm_LFP_pzz1propylwavephasepower071220_hippocampusLFP80.mat';

% prePathName='F:\Datos summary CaMKII111720\PRP drgAnalysisBatchLFPCaMKII case 24 output for summary\';
% prePathName='/Users/restrepd/Documents/Projects/CaMKII_analysis/PRP drgAnalysisBatchLFPCaMKII case 24 output for 80/';

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
% preFileName{1}='spm_LFP_acetowavephasepower32620_prefrontalLFP80.mat';
% preFileName{2}='spm_LFP_ethylbenwavephasepower41420_prefrontalLFP80.mat';
% preFileName{3}='spm_LFP_ethylwavephasepower3520_prefrontalLFP80.mat';
% preFileName{4}='spm_LFP_acetowavephasepower32620 2_prefrontalLFP80.mat';
% preFileName{5}='spm_LFP_pz1ethyllwavephasepower0213020_prefrontalLFP80.mat';
% preFileName{6}='spm_LFP_pz1propylwavephasepower013020_prefrontalLFP80.mat';
% preFileName{7}='spm_LFP_pzz1ethyllwavephasepower043020_prefrontalLFP80.mat';
% preFileName{8}='spm_LFP_pzz1propylwavephasepower071220_prefrontalLFP80.mat';


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
        
        ax=gca;ax.LineWidth=3;
        
        set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
        hold on
        
        bar_lab_loc=[];
        no_ev_labels=0;
        ii_gr_included=0;
        bar_offset = 0;
        
       
        
        for evNo=1:2
            
            for per_ii=2:-1:1
                
                for grNo=1:2
                    
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
            'VariableNames',{'MI','group','perCorr','event'});
        mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
            ,'CategoricalVars',[2,3,4])
        
        
        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for average PRP calculated per mouse per odor pair for ' peak_label{peak+1} ' theta' PACnames{pacii} ' hippocampus\n'])
        [output_data] = drgMutiRanksumorTtest(input_data);
        
        
    end
end

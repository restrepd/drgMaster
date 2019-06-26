function drgSummaryBatchPAC
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'drgbDiscPar' file listing 'Discriminant_*.mat' output files
%from drgLFPDiscriminantBatch
%
%Performs summary analises for LDA and  PAC analysis performed with case 19
%of drgAnalysisBatchLFP


close all
clear all

evTypeLabels={'S+','S-',};
per_ii_labels={'Proficient','Naive'}

PACnames{1}='Beta';
PACnames{2}='Low gamma';
PACnames{3}='High gamma';


PathName{1}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 3 Batch PAC/Daniel_data/';
FileName{1}='Olfactorypaper04202019_APEB.mat';
odorPairName{1}='APEBloc1';

PathName{2}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 3 Batch PAC/Daniel_data/';
FileName{2}='Olfactorypaper04202019_EAPA.mat';
odorPairName{2}='EAPA1loc1';

PathName{3}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 3 Batch PAC/Daniel_data/';
FileName{3}='Olfactorypaper04202019_IAMO.mat';
odorPairName{3}='IAMOloc1';

PathName{4}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 3 Batch PAC/Justin_data/';
FileName{4}='spm_LFP_PACpower04172019_EAPA.mat';
odorPairName{4}='EAPA2loc2';

PathName{5}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 3 Batch PAC/Justin_data/';
FileName{5}='spm_LFP_PACpower04172019_IAAP.mat';
odorPairName{5}='IAAPloc2';

PathName{6}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 3 Batch PAC/Justin_data/';
FileName{6}='spm_LFP_PACpower04172019_IAMO.mat';
odorPairName{6}='IAMOloc2';

figNo=0;



for PACii=1:3
    
    %Plot the bar graph for MI
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.3 .3 .6 .4])
    hold on
    
    glm_pav_ii=0;
    glm_pav=[];
    glm_mi_ii=0;
    glm_mi=[];
    ii_stats_mi=0;
    ii_stats_pav=0;
    p_pav=[];
    p_mi_stats=[];
    
    bar_ii=1;
    for fileNo=1:length(PathName)
        
        load([PathName{fileNo} FileName{fileNo}])
        
        for ii_rank=1:length(out_mi_rank(PACii).mi_rank)
            per_ii=out_mi_rank(PACii).mi_rank(ii_rank).per_ii;
            evNo=out_mi_rank(PACii).mi_rank(ii_rank).evNo;
            if evNo==1
                if per_ii==1
                    %S+ Proficient
                    this_bar_ii=bar_ii+4;
                    mean_mi=mean(out_mi_rank(PACii).mi_rank(ii_rank).mi);
                    miCI = bootci(1000, {@mean, out_mi_rank(PACii).mi_rank(ii_rank).mi},'type','cper');
                    hp1=bar(this_bar_ii,mean_mi,'r');
                    plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
                    
                    these_mi=out_mi_rank(PACii).mi_rank(ii_rank).mi;
                    glm_mi.data(glm_mi_ii+1:glm_mi_ii+length(these_mi))=these_mi;
                    glm_mi.perCorr(glm_mi_ii+1:glm_mi_ii+length(these_mi))=per_ii*ones(1,length(these_mi));
                    glm_mi.event(glm_mi_ii+1:glm_mi_ii+length(these_mi))=1*ones(1,length(these_mi));
                    if fileNo<=3
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=1*ones(1,length(these_mi));
                    else
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=2*ones(1,length(these_mi));
                    end
                    glm_mi_ii=glm_mi_ii+length(these_mi);
                    
                    ii_stats_mi=ii_stats_mi+1;
                    p_mi_stats(ii_stats_mi).data=these_mi;
                    if per_ii==1
                        p_mi_stats(ii_stats_mi).description=['Proficient CS+ ' odorPairName{fileNo}];
                    else
                        p_mi_stats(ii_stats_mi).description=['Naive CS+ ' odorPairName{fileNo}];
                    end
                    
                else
                    %S+ Naive
                    this_bar_ii=bar_ii+3;
                    mean_mi=mean(out_mi_rank(PACii).mi_rank(ii_rank).mi);
                    miCI = bootci(1000, {@mean, out_mi_rank(PACii).mi_rank(ii_rank).mi},'type','cper');
                    hp2=bar(this_bar_ii,mean_mi,'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7]);
                    plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
                    
                    these_mi=out_mi_rank(PACii).mi_rank(ii_rank).mi;
                    glm_mi.data(glm_mi_ii+1:glm_mi_ii+length(these_mi))=these_mi;
                    glm_mi.perCorr(glm_mi_ii+1:glm_mi_ii+length(these_mi))=per_ii*ones(1,length(these_mi));
                    glm_mi.event(glm_mi_ii+1:glm_mi_ii+length(these_mi))=1*ones(1,length(these_mi));
                    if fileNo<=3
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=1*ones(1,length(these_mi));
                    else
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=2*ones(1,length(these_mi));
                    end
                    glm_mi_ii=glm_mi_ii+length(these_mi);
                    
                    ii_stats_mi=ii_stats_mi+1;
                    p_mi_stats(ii_stats_mi).data=these_mi;
                    if per_ii==1
                        p_mi_stats(ii_stats_mi).description=['Proficient CS+ ' odorPairName{fileNo}];
                    else
                        p_mi_stats(ii_stats_mi).description=['Naive CS+ ' odorPairName{fileNo}];
                    end
                end
            else
                if per_ii==1
                    %S- Proficient
                    this_bar_ii=bar_ii+1;
                    mean_mi=mean(out_mi_rank(PACii).mi_rank(ii_rank).mi);
                    miCI = bootci(1000, {@mean, out_mi_rank(PACii).mi_rank(ii_rank).mi},'type','cper');
                    hp3=bar(this_bar_ii,mean_mi,'b');
                    plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
                    
                    these_mi=out_mi_rank(PACii).mi_rank(ii_rank).mi;
                    glm_mi.data(glm_mi_ii+1:glm_mi_ii+length(these_mi))=these_mi;
                    glm_mi.perCorr(glm_mi_ii+1:glm_mi_ii+length(these_mi))=per_ii*ones(1,length(these_mi));
                    glm_mi.event(glm_mi_ii+1:glm_mi_ii+length(these_mi))=2*ones(1,length(these_mi));
                    if fileNo<=3
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=1*ones(1,length(these_mi));
                    else
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=2*ones(1,length(these_mi));
                    end
                    glm_mi_ii=glm_mi_ii+length(these_mi);
                    
                    ii_stats_mi=ii_stats_mi+1;
                    p_mi_stats(ii_stats_mi).data=these_mi;
                    if per_ii==1
                        p_mi_stats(ii_stats_mi).description=['Proficient CS- ' odorPairName{fileNo}];
                    else
                        p_mi_stats(ii_stats_mi).description=['Naive CS- ' odorPairName{fileNo}];
                    end
                else
                    %S- Naive
                    this_bar_ii=bar_ii;
                    mean_mi=mean(out_mi_rank(PACii).mi_rank(ii_rank).mi);
                    miCI = bootci(1000, {@mean, out_mi_rank(PACii).mi_rank(ii_rank).mi},'type','cper');
                    hp4=bar(this_bar_ii,mean_mi,'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1]);
                    plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
                    
                    these_mi=out_mi_rank(PACii).mi_rank(ii_rank).mi;
                    glm_mi.data(glm_mi_ii+1:glm_mi_ii+length(these_mi))=these_mi;
                    glm_mi.perCorr(glm_mi_ii+1:glm_mi_ii+length(these_mi))=per_ii*ones(1,length(these_mi));
                    glm_mi.event(glm_mi_ii+1:glm_mi_ii+length(these_mi))=2*ones(1,length(these_mi));
                    if fileNo<=3
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=1*ones(1,length(these_mi));
                    else
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=2*ones(1,length(these_mi));
                    end
                    glm_mi_ii=glm_mi_ii+length(these_mi);
                    
                    ii_stats_mi=ii_stats_mi+1;
                    p_mi_stats(ii_stats_mi).data=these_mi;
                    if per_ii==1
                        p_mi_stats(ii_stats_mi).description=['Proficient CS- ' odorPairName{fileNo}];
                    else
                        p_mi_stats(ii_stats_mi).description=['Naive CS- ' odorPairName{fileNo}];
                    end
                end
                
            end
            
            
        end
        bar_ii=bar_ii+7;
    end
    
    legend([hp1 hp2 hp3 hp4],{'S+ Proficient','S+ Naive','S- Proficient','S- Naive'})
    xticks([3.5 10.5 17.5 24.5 31.5 38.5])
    xticklabels({'APEB','EAPA1loc1','EAPA2loc2','IAAPloc2','IAMOloc1','IAMOloc2'})
    ylabel('MI')
    
    switch PACii
        case 1
            title('Modulation index Theta/Beta')
        case 2
            title('Modulation index Theta/Low Gamma')
        case 3
            title('Modulation index Theta/High Gamma')
    end
    
    %Perform the glm for mi
    fprintf(1, ['\n\nglm for MI for Theta/' PACnames{PACii} '\n'])
    tbl = table(glm_mi.data',glm_mi.perCorr',glm_mi.event',glm_mi.location',...
        'VariableNames',{'MI','proficiency','event','location'});
    mdl = fitglm(tbl,'MI~proficiency+location+event+proficiency*event*location'...
        ,'CategoricalVars',[2,3 4])
     
    %Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for MI per mouse for Theta/' PACnames{PACii} '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(p_mi_stats);
        fprintf(1, '\n\n')
    catch
    end
    
    %Plot the bar graph for MI
    figNo=figNo+1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.3 .3 .6 .4])
    hold on
    
    
    bar_ii=1;
    for fileNo=1:length(PathName)
        
        load([PathName{fileNo} FileName{fileNo}])
        
        for ii_rank=1:length(out_PAvar_rank(PACii).PAvar_rank)
            per_ii=out_PAvar_rank(PACii).PAvar_rank(ii_rank).per_ii;
            evNo=out_PAvar_rank(PACii).PAvar_rank(ii_rank).evNo;
            if evNo==1
                if per_ii==1
                    this_bar_ii=bar_ii+4;
                    mean_pav=mean(out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar);
                    miCI = bootci(1000, {@mean, out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar},'type','cper');
                    hp1=bar(this_bar_ii,mean_pav,'r');
                    plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
                    
                    these_pav=out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar;
                    glm_pav.data(glm_pav_ii+1:glm_pav_ii+length(these_pav))=these_pav;
                    glm_pav.perCorr(glm_pav_ii+1:glm_pav_ii+length(these_pav))=per_ii*ones(1,length(these_pav));
                    glm_pav.event(glm_pav_ii+1:glm_pav_ii+length(these_pav))=1*ones(1,length(these_pav));
                    if fileNo<=3
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=1*ones(1,length(these_pav));
                    else
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=2*ones(1,length(these_pav));
                    end
                    glm_pav_ii=glm_pav_ii+length(these_pav);
                    
                    ii_stats_pav=ii_stats_pav+1;
                    p_pav_stats(ii_stats_pav).data=these_pav;
                    if per_ii==1
                        p_pav_stats(ii_stats_pav).description=['Proficient CS+ ' odorPairName{fileNo}];
                    else
                        p_pav_stats(ii_stats_pav).description=['Naive CS+ ' odorPairName{fileNo}];
                    end
                else
                    this_bar_ii=bar_ii+3;
                    mean_pav=mean(out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar);
                    miCI = bootci(1000, {@mean, out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar},'type','cper');
                    hp2=bar(this_bar_ii,mean_pav,'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7]);
                    plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
                    
                    these_pav=out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar;
                    glm_pav.data(glm_pav_ii+1:glm_pav_ii+length(these_pav))=these_pav;
                    glm_pav.perCorr(glm_pav_ii+1:glm_pav_ii+length(these_pav))=per_ii*ones(1,length(these_pav));
                    glm_pav.event(glm_pav_ii+1:glm_pav_ii+length(these_pav))=1*ones(1,length(these_pav));
                    if fileNo<=3
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=1*ones(1,length(these_pav));
                    else
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=2*ones(1,length(these_pav));
                    end
                    glm_pav_ii=glm_pav_ii+length(these_pav);
                    
                    ii_stats_pav=ii_stats_pav+1;
                    p_pav_stats(ii_stats_pav).data=these_pav;
                    if per_ii==1
                        p_pav_stats(ii_stats_pav).description=['Proficient CS+ ' odorPairName{fileNo}];
                    else
                        p_pav_stats(ii_stats_pav).description=['Naive CS+ ' odorPairName{fileNo}];
                    end
                end
            else
                if per_ii==1
                    this_bar_ii=bar_ii+1;
                    mean_pav=mean(out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar);
                    miCI = bootci(1000, {@mean, out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar},'type','cper');
                    hp3=bar(this_bar_ii,mean_pav,'b');
                    plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
                    
                    these_pav=out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar;
                    glm_pav.data(glm_pav_ii+1:glm_pav_ii+length(these_pav))=these_pav;
                    glm_pav.perCorr(glm_pav_ii+1:glm_pav_ii+length(these_pav))=per_ii*ones(1,length(these_pav));
                    glm_pav.event(glm_pav_ii+1:glm_pav_ii+length(these_pav))=2*ones(1,length(these_pav));
                    if fileNo<=3
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=1*ones(1,length(these_pav));
                    else
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=2*ones(1,length(these_pav));
                    end
                    glm_pav_ii=glm_pav_ii+length(these_pav);
                    
                    ii_stats_pav=ii_stats_pav+1;
                    p_pav_stats(ii_stats_pav).data=these_pav;
                    if per_ii==1
                        p_pav_stats(ii_stats_pav).description=['Proficient CS+ ' odorPairName{fileNo}];
                    else
                        p_pav_stats(ii_stats_pav).description=['Naive CS+ ' odorPairName{fileNo}];
                    end
                else
                    this_bar_ii=bar_ii;
                    mean_pav=mean(out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar);
                    miCI = bootci(1000, {@mean, out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar},'type','cper');
                    hp4=bar(this_bar_ii,mean_pav,'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1]);
                    plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
                    
                    these_pav=out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar;
                    glm_pav.data(glm_pav_ii+1:glm_pav_ii+length(these_pav))=these_pav;
                    glm_pav.perCorr(glm_pav_ii+1:glm_pav_ii+length(these_pav))=per_ii*ones(1,length(these_pav));
                    glm_pav.event(glm_pav_ii+1:glm_pav_ii+length(these_pav))=2*ones(1,length(these_pav));
                    if fileNo<=3
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=1*ones(1,length(these_pav));
                    else
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=2*ones(1,length(these_pav));
                    end
                    glm_pav_ii=glm_pav_ii+length(these_pav);
                    
                    ii_stats_pav=ii_stats_pav+1;
                    p_pav_stats(ii_stats_pav).data=these_pav;
                    if per_ii==1
                        p_pav_stats(ii_stats_pav).description=['Proficient CS+ ' odorPairName{fileNo}];
                    else
                        p_pav_stats(ii_stats_pav).description=['Naive CS+ ' odorPairName{fileNo}];
                    end
                end
                
            end
            
            
        end
        bar_ii=bar_ii+7;
    end
    
    legend([hp1 hp2 hp3 hp4],{'S+ Proficient','S+ Naive','S- Proficient','S- Naive'})
    xticks([3.5 10.5 17.5 24.5 31.5 38.5])
    xticklabels({'APEB','EAPA1loc1','EAPA2loc2','IAAPloc2','IAMOloc1','IAMOloc2'})
    ylabel('PA variance degrees^2')
    switch PACii
        case 1
            title('Peak angle variance Theta/Beta')
        case 2
            title('Peak angle variance Theta/Low Gamma')
        case 3
            title('Peak angle variance Theta/High Gamma')
    end
    
    ylim([0 3500])
    
    %Perform the glm for pav
    fprintf(1, ['\n\nglm for peak angle variance for Theta/' PACnames{PACii} '\n'])
    tbl = table(glm_pav.data',glm_pav.perCorr',glm_pav.event',glm_pav.location',...
        'VariableNames',{'pav','proficiency','event','location'});
    mdl = fitglm(tbl,'pav~proficiency+location+event+proficiency*event*location'...
        ,'CategoricalVars',[2,3 4])
    
    %Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for peak angle variance for Theta/' PACnames{PACii} '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(p_pav_stats);
        fprintf(1, '\n\n')
    catch
    end
    
end
pffft=1;
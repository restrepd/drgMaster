function drgSummaryBatchPACfig2h
%Performs summary analysis for PAC analysis performed with case 19
%of drgAnalysisBatchLFP

warning('off')

close all
clear all

evTypeLabels={'S+','S-',};
per_ii_labels={'Proficient','Naive'}

PACnames{1}='Beta';
PACnames{2}='Low gamma';
PACnames{3}='High gamma';

try
    PathName{1}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 2 ext Batch PAC/Daniel_data/Final Data/';
    cd(PathName{1})
catch
    PathName{1}='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 2 ext Batch PAC\Daniel_data\Final Data\';
end
FileName{1}='Olfactorypaper04202019_APEB.mat';
odorPairName{1}='APEBloc1';

try
    PathName{2}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 2 ext Batch PAC/Daniel_data/Final Data/';
    cd(PathName{2})
catch
    PathName{2}='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 2 ext Batch PAC\Daniel_data\Final Data\';
end
FileName{2}='Olfactorypaper04202019_EAPA.mat';
odorPairName{2}='EAPA1loc1';

try
    PathName{3}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 2 ext Batch PAC/Daniel_data/Final Data/';
    cd(PathName{3})
catch
    PathName{3}='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 2 ext Batch PAC\Daniel_data\Final Data\';
end
FileName{3}='Olfactorypaper04202019_IAMO.mat';
odorPairName{3}='IAMOloc1';

try
    PathName{4}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 2 ext Batch PAC/Justin_data/Final Data/';
    cd(PathName{4})
catch
    PathName{4}='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 2 ext Batch PAC\Justin_data\Final Data\';
end
FileName{4}='spm_LFP_PACpower04172019_EAPA.mat';
odorPairName{4}='EAPA2loc2';

try
    PathName{5}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 2 ext Batch PAC/Justin_data/Final Data/';
    cd(PathName{5})
catch
    PathName{5}='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 2 ext Batch PAC\Justin_data\Final Data\';
end
FileName{5}='spm_LFP_PACpower04172019_IAAP.mat';
odorPairName{5}='IAAPloc2';

try
    PathName{6}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 2 ext Batch PAC/Justin_data/Final Data/';
    cd(PathName{6})
catch
    PathName{6}='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 2 ext Batch PAC\Justin_data\Final Data\';
end
FileName{6}='spm_LFP_PACpower04172019_IAMO.mat';
odorPairName{6}='IAMOloc2';

figNo=0;



for PACii=[1 3]
    
    %Plot the bar graph for MI for each odor pair
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
    
    MI_mean_ii=0;
    MI_mean=[];
    MI_mean_evNo=[];
    MI_mean_per_ii=[];
    MI_mean_exp=[];
    
    for fileNo=1:length(PathName)
        
        out_mi_rank=[];
        out_PAvar_rank=[];
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
                    
                    %Save the mean values
                    MI_mean_ii=MI_mean_ii+1;
                    MI_mean(MI_mean_ii)=mean_mi;
                    MI_mean_evNo(MI_mean_ii)=evNo;
                    MI_mean_per_ii(MI_mean_ii)=per_ii;
                    
                    these_mi=out_mi_rank(PACii).mi_rank(ii_rank).mi;
                    glm_mi.data(glm_mi_ii+1:glm_mi_ii+length(these_mi))=these_mi;
                    glm_mi.perCorr(glm_mi_ii+1:glm_mi_ii+length(these_mi))=per_ii*ones(1,length(these_mi));
                    glm_mi.event(glm_mi_ii+1:glm_mi_ii+length(these_mi))=1*ones(1,length(these_mi));
                    if fileNo<=3
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=1*ones(1,length(these_mi));
                        MI_mean_exp(MI_mean_ii)=1;
                    else
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=2*ones(1,length(these_mi));
                        MI_mean_exp(MI_mean_ii)=2;
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
                    
                    %Save the mean values
                    MI_mean_ii=MI_mean_ii+1;
                    MI_mean(MI_mean_ii)=mean_mi;
                    MI_mean_evNo(MI_mean_ii)=evNo;
                    MI_mean_per_ii(MI_mean_ii)=per_ii;
                    
                    these_mi=out_mi_rank(PACii).mi_rank(ii_rank).mi;
                    glm_mi.data(glm_mi_ii+1:glm_mi_ii+length(these_mi))=these_mi;
                    glm_mi.perCorr(glm_mi_ii+1:glm_mi_ii+length(these_mi))=per_ii*ones(1,length(these_mi));
                    glm_mi.event(glm_mi_ii+1:glm_mi_ii+length(these_mi))=1*ones(1,length(these_mi));
                    if fileNo<=3
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=1*ones(1,length(these_mi));
                        MI_mean_exp(MI_mean_ii)=1;
                    else
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=2*ones(1,length(these_mi));
                        MI_mean_exp(MI_mean_ii)=2;
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
                    
                    %Save the mean values
                    MI_mean_ii=MI_mean_ii+1;
                    MI_mean(MI_mean_ii)=mean_mi;
                    MI_mean_evNo(MI_mean_ii)=evNo;
                    MI_mean_per_ii(MI_mean_ii)=per_ii;
                    
                    these_mi=out_mi_rank(PACii).mi_rank(ii_rank).mi;
                    glm_mi.data(glm_mi_ii+1:glm_mi_ii+length(these_mi))=these_mi;
                    glm_mi.perCorr(glm_mi_ii+1:glm_mi_ii+length(these_mi))=per_ii*ones(1,length(these_mi));
                    glm_mi.event(glm_mi_ii+1:glm_mi_ii+length(these_mi))=2*ones(1,length(these_mi));
                    if fileNo<=3
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=1*ones(1,length(these_mi));
                        MI_mean_exp(MI_mean_ii)=1;
                    else
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=2*ones(1,length(these_mi));
                        MI_mean_exp(MI_mean_ii)=2;
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
                    
                    %Save the mean values
                    MI_mean_ii=MI_mean_ii+1;
                    MI_mean(MI_mean_ii)=mean_mi;
                    MI_mean_evNo(MI_mean_ii)=evNo;
                    MI_mean_per_ii(MI_mean_ii)=per_ii;
                    
                    these_mi=out_mi_rank(PACii).mi_rank(ii_rank).mi;
                    glm_mi.data(glm_mi_ii+1:glm_mi_ii+length(these_mi))=these_mi;
                    glm_mi.perCorr(glm_mi_ii+1:glm_mi_ii+length(these_mi))=per_ii*ones(1,length(these_mi));
                    glm_mi.event(glm_mi_ii+1:glm_mi_ii+length(these_mi))=2*ones(1,length(these_mi));
                    if fileNo<=3
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=1*ones(1,length(these_mi));
                        MI_mean_exp(MI_mean_ii)=1;
                    else
                        glm_mi.location(glm_mi_ii+1:glm_mi_ii+length(these_mi))=2*ones(1,length(these_mi));
                        MI_mean_exp(MI_mean_ii)=2;
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
    xticklabels({'APEBloc1','EAPA1loc1','IAMOloc1','EAPA2loc2','IAAPloc2','IAMOloc2'})
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
    
    %Now plot the MI mean per experiment (location)
    figNo=figNo+1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.3 .3 .6 .4])
    hold on
    
    glm_MI_mean=[];
    glm_MI_mii=0;
    
    ii_stats_m_MI=0;
    p_MI_mean_stats=[]
    
    this_bar_ii=0;
    
    for exp=1:2
        
        %S- naive
        evNo=2;
        per_ii=2;
        this_bar_ii=this_bar_ii+1;
        these_ii=(MI_mean_evNo==evNo)&(MI_mean_exp==exp)&(MI_mean_per_ii==per_ii);
        mean_MI_m=mean(MI_mean(these_ii));
        miCI = bootci(1000, {@mean, MI_mean(these_ii)},'type','cper');
        hp4=bar(this_bar_ii,mean_MI_m,'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1]);
        plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
        
        glm_MI_mean.data(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=MI_mean(these_ii);
        glm_MI_mean.perCorr(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=per_ii*ones(1,sum(these_ii));
        glm_MI_mean.event(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=evNo*ones(1,sum(these_ii));
        glm_MI_mean.exp(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=exp*ones(1,sum(these_ii));
        glm_MI_mii=glm_MI_mii+sum(these_ii);
        
        ii_stats_m_MI=ii_stats_m_MI+1;
        p_MI_mean_stats(ii_stats_m_MI).data=MI_mean(these_ii);
        p_MI_mean_stats(ii_stats_m_MI).description=['Naive S- Exp' num2str(exp)];
        
        %S- proficient
        evNo=2;
        per_ii=1;
        this_bar_ii=this_bar_ii+1;
        these_ii=(MI_mean_evNo==evNo)&(MI_mean_exp==exp)&(MI_mean_per_ii==per_ii);
        mean_MI_m=mean(MI_mean(these_ii));
        miCI = bootci(1000, {@mean, MI_mean(these_ii)},'type','cper');
        hp3=bar(this_bar_ii,mean_MI_m,'b');
        plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
        
        glm_MI_mean.data(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=MI_mean(these_ii);
        glm_MI_mean.perCorr(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=per_ii*ones(1,sum(these_ii));
        glm_MI_mean.event(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=evNo*ones(1,sum(these_ii));
        glm_MI_mean.exp(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=exp*ones(1,sum(these_ii));
        glm_MI_mii=glm_MI_mii+sum(these_ii);
        
        ii_stats_m_MI=ii_stats_m_MI+1;
        p_MI_mean_stats(ii_stats_m_MI).data=MI_mean(these_ii);
        p_MI_mean_stats(ii_stats_m_MI).description=['Proficient S- Exp' num2str(exp)];
        
        %S+ naive
        evNo=1;
        per_ii=2;
        this_bar_ii=this_bar_ii+2;
        these_ii=(MI_mean_evNo==evNo)&(MI_mean_exp==exp)&(MI_mean_per_ii==per_ii);
        mean_MI_m=mean(MI_mean(these_ii));
        miCI = bootci(1000, {@mean, MI_mean(these_ii)},'type','cper');
        hp2=bar(this_bar_ii,mean_MI_m,'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7]);
        plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
        
        glm_MI_mean.data(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=MI_mean(these_ii);
        glm_MI_mean.perCorr(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=per_ii*ones(1,sum(these_ii));
        glm_MI_mean.event(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=evNo*ones(1,sum(these_ii));
        glm_MI_mean.exp(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=exp*ones(1,sum(these_ii));
        glm_MI_mii=glm_MI_mii+sum(these_ii);
        
        ii_stats_m_MI=ii_stats_m_MI+1;
        p_MI_mean_stats(ii_stats_m_MI).data=MI_mean(these_ii);
        p_MI_mean_stats(ii_stats_m_MI).description=['Naive S+  Exp' num2str(exp)];
        
        %S+ proficient
        evNo=1;
        per_ii=1;
        this_bar_ii=this_bar_ii+1;
        these_ii=(MI_mean_evNo==evNo)&(MI_mean_exp==exp)&(MI_mean_per_ii==per_ii);
        mean_MI_m=mean(MI_mean(these_ii));
        miCI = bootci(1000, {@mean, MI_mean(these_ii)},'type','cper');
        hp1=bar(this_bar_ii,mean_MI_m,'r');
        plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
        
        glm_MI_mean.data(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=MI_mean(these_ii);
        glm_MI_mean.perCorr(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=per_ii*ones(1,sum(these_ii));
        glm_MI_mean.event(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=evNo*ones(1,sum(these_ii));
        glm_MI_mean.exp(glm_MI_mii+1:glm_MI_mii+sum(these_ii))=exp*ones(1,sum(these_ii));
        glm_MI_mii=glm_MI_mii+sum(these_ii);
        
        ii_stats_m_MI=ii_stats_m_MI+1;
        p_MI_mean_stats(ii_stats_m_MI).data=MI_mean(these_ii);
        p_MI_mean_stats(ii_stats_m_MI).description=['Proficient S+ Exp' num2str(exp)];
        
        this_bar_ii=this_bar_ii+3;
    end
    
    legend([hp1 hp2 hp3 hp4],{'S+ Proficient','S+ Naive','S- Proficient','S- Naive'})
    xticks([3 11])
    xticklabels({'Exp1','Exp2'})
    ylabel('MI')
    switch PACii
        case 1
            title('Mean MI per experiment Theta/Beta')
        case 2
            title('Mena MI per experiment Theta/Low Gamma')
        case 3
            title('Mean MI per experiment Theta/High Gamma')
    end
    
    ylim([0 0.02])
    
    %Perform the glm for MI_mean
    fprintf(1, ['\n\nglm for per experiment mean MI for Theta/' PACnames{PACii} '\n'])
    tbl = table(glm_MI_mean.data',glm_MI_mean.perCorr',glm_MI_mean.event',glm_MI_mean.exp',...
        'VariableNames',{'pav','proficiency','event','experiment'});
    mdl = fitglm(tbl,'pav~proficiency+experiment+event+proficiency*event*experiment'...
        ,'CategoricalVars',[2,3,4])
    
    %Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for per experiment mean MI for Theta/' PACnames{PACii} '\n'])
    try
        [output_data] = drgMutiRanksumorTtest( p_MI_mean_stats);
        fprintf(1, '\n\n')
    catch
    end
    
    %Plot the bar graph for peak variance for each odor pair
    figNo=figNo+1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.3 .3 .6 .4])
    hold on
    
    %Save the mean values
    pavar_mean=[];
    pavar_mean_evNo=[];
    pavar_mean_per_ii=[];
    pavar_mean_exp=[];
    pavar_mean_ii=0;
    
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
                    
                    %Save the mean values
                    pavar_mean_ii=pavar_mean_ii+1;
                    pavar_mean(pavar_mean_ii)=mean_pav;
                    pavar_mean_evNo(pavar_mean_ii)=evNo;
                    pavar_mean_per_ii(pavar_mean_ii)=per_ii;
                    
                    these_pav=out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar;
                    glm_pav.data(glm_pav_ii+1:glm_pav_ii+length(these_pav))=these_pav;
                    glm_pav.perCorr(glm_pav_ii+1:glm_pav_ii+length(these_pav))=per_ii*ones(1,length(these_pav));
                    glm_pav.event(glm_pav_ii+1:glm_pav_ii+length(these_pav))=1*ones(1,length(these_pav));
                    if fileNo<=3
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=1*ones(1,length(these_pav));
                        pavar_mean_exp(pavar_mean_ii)=1;
                    else
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=2*ones(1,length(these_pav));
                        pavar_mean_exp(pavar_mean_ii)=2;
                    end
                    glm_pav_ii=glm_pav_ii+length(these_pav);
                    
                    ii_stats_pav=ii_stats_pav+1;
                    p_pav_stats(ii_stats_pav).data=these_pav;
                    
                    p_pav_stats(ii_stats_pav).description=['Proficient CS+ ' odorPairName{fileNo}];
                    
                else
                    this_bar_ii=bar_ii+3;
                    mean_pav=mean(out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar);
                    miCI = bootci(1000, {@mean, out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar},'type','cper');
                    hp2=bar(this_bar_ii,mean_pav,'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7]);
                    plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
                    
                    %Save the mean values
                    pavar_mean_ii=pavar_mean_ii+1;
                    pavar_mean(pavar_mean_ii)=mean_pav;
                    pavar_mean_evNo(pavar_mean_ii)=evNo;
                    pavar_mean_per_ii(pavar_mean_ii)=per_ii;
                    
                    these_pav=out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar;
                    glm_pav.data(glm_pav_ii+1:glm_pav_ii+length(these_pav))=these_pav;
                    glm_pav.perCorr(glm_pav_ii+1:glm_pav_ii+length(these_pav))=per_ii*ones(1,length(these_pav));
                    glm_pav.event(glm_pav_ii+1:glm_pav_ii+length(these_pav))=1*ones(1,length(these_pav));
                    if fileNo<=3
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=1*ones(1,length(these_pav));
                        pavar_mean_exp(pavar_mean_ii)=1;
                        
                    else
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=2*ones(1,length(these_pav));
                        pavar_mean_exp(pavar_mean_ii)=2;
                    end
                    glm_pav_ii=glm_pav_ii+length(these_pav);
                    
                    ii_stats_pav=ii_stats_pav+1;
                    p_pav_stats(ii_stats_pav).data=these_pav;
                    
                    p_pav_stats(ii_stats_pav).description=['Naive CS+ ' odorPairName{fileNo}];
                    
                end
            else
                if per_ii==1
                    this_bar_ii=bar_ii+1;
                    mean_pav=mean(out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar);
                    miCI = bootci(1000, {@mean, out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar},'type','cper');
                    hp3=bar(this_bar_ii,mean_pav,'b');
                    plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
                    
                    %Save the mean values
                    pavar_mean_ii=pavar_mean_ii+1;
                    pavar_mean(pavar_mean_ii)=mean_pav;
                    pavar_mean_evNo(pavar_mean_ii)=evNo;
                    pavar_mean_per_ii(pavar_mean_ii)=per_ii;
                    
                    these_pav=out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar;
                    glm_pav.data(glm_pav_ii+1:glm_pav_ii+length(these_pav))=these_pav;
                    glm_pav.perCorr(glm_pav_ii+1:glm_pav_ii+length(these_pav))=per_ii*ones(1,length(these_pav));
                    glm_pav.event(glm_pav_ii+1:glm_pav_ii+length(these_pav))=2*ones(1,length(these_pav));
                    if fileNo<=3
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=1*ones(1,length(these_pav));
                        pavar_mean_exp(pavar_mean_ii)=1;
                    else
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=2*ones(1,length(these_pav));
                        pavar_mean_exp(pavar_mean_ii)=2;
                    end
                    glm_pav_ii=glm_pav_ii+length(these_pav);
                    
                    ii_stats_pav=ii_stats_pav+1;
                    p_pav_stats(ii_stats_pav).data=these_pav;
                    
                    p_pav_stats(ii_stats_pav).description=['Proficient CS- ' odorPairName{fileNo}];
                    
                else
                    this_bar_ii=bar_ii;
                    mean_pav=mean(out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar);
                    miCI = bootci(1000, {@mean, out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar},'type','cper');
                    hp4=bar(this_bar_ii,mean_pav,'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1]);
                    plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
                    
                    %Save the mean values
                    pavar_mean_ii=pavar_mean_ii+1;
                    pavar_mean(pavar_mean_ii)=mean_pav;
                    pavar_mean_evNo(pavar_mean_ii)=evNo;
                    pavar_mean_per_ii(pavar_mean_ii)=per_ii;
                    
                    these_pav=out_PAvar_rank(PACii).PAvar_rank(ii_rank).meanPAvar;
                    glm_pav.data(glm_pav_ii+1:glm_pav_ii+length(these_pav))=these_pav;
                    glm_pav.perCorr(glm_pav_ii+1:glm_pav_ii+length(these_pav))=per_ii*ones(1,length(these_pav));
                    glm_pav.event(glm_pav_ii+1:glm_pav_ii+length(these_pav))=2*ones(1,length(these_pav));
                    if fileNo<=3
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=1*ones(1,length(these_pav));
                        pavar_mean_exp(pavar_mean_ii)=1;
                    else
                        glm_pav.location(glm_pav_ii+1:glm_pav_ii+length(these_pav))=2*ones(1,length(these_pav));
                        pavar_mean_exp(pavar_mean_ii)=2;
                    end
                    glm_pav_ii=glm_pav_ii+length(these_pav);
                    
                    ii_stats_pav=ii_stats_pav+1;
                    p_pav_stats(ii_stats_pav).data=these_pav;
                    
                    p_pav_stats(ii_stats_pav).description=['Naive CS- ' odorPairName{fileNo}];
                    
                end
                
            end
            
            
        end
        bar_ii=bar_ii+7;
    end
    
    legend([hp1 hp2 hp3 hp4],{'S+ Proficient','S+ Naive','S- Proficient','S- Naive'})
    xticks([3.5 10.5 17.5 24.5 31.5 38.5])
    xticklabels({'APEBloc1','EAPA1loc1','IAMOloc1','EAPA2loc2','IAAPloc2','IAMOloc2'})
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
    
   
    %Now plot the peak variance mean per experiment (location)
    figNo=figNo+1;
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.3 .3 .6 .4])
    hold on
    
    glm_pav_mean=[];
    glm_pav_mii=0;
    
    p_pav_mean_stats=[];
    ii_stats_m_pav=0;
    
    this_bar_ii=0;
    
    for exp=1:2
        
        %S- naive
        evNo=2;
        per_ii=2;
        this_bar_ii=this_bar_ii+1;
        these_ii=(pavar_mean_evNo==evNo)&(pavar_mean_exp==exp)&(pavar_mean_per_ii==per_ii);
        mean_pav_m=mean(pavar_mean(these_ii));
        miCI = bootci(1000, {@mean, pavar_mean(these_ii)},'type','cper');
        hp4=bar(this_bar_ii,mean_pav_m,'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1]);
        plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
        
        glm_pav_mean.data(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=pavar_mean(these_ii);
        glm_pav_mean.perCorr(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=per_ii*ones(1,sum(these_ii));
        glm_pav_mean.event(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=evNo*ones(1,sum(these_ii));
        glm_pav_mean.exp(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=exp*ones(1,sum(these_ii));
        glm_pav_mii=glm_pav_mii+sum(these_ii);
        
        ii_stats_m_pav=ii_stats_m_pav+1;
        p_pav_mean_stats(ii_stats_m_pav).data=pavar_mean(these_ii);
        p_pav_mean_stats(ii_stats_m_pav).description=['Naive S- Exp' num2str(exp)];
        
        %S- proficient
        evNo=2;
        per_ii=1;
        this_bar_ii=this_bar_ii+1;
        these_ii=(pavar_mean_evNo==evNo)&(pavar_mean_exp==exp)&(pavar_mean_per_ii==per_ii);
        mean_pav_m=mean(pavar_mean(these_ii));
        miCI = bootci(1000, {@mean, pavar_mean(these_ii)},'type','cper');
        hp3=bar(this_bar_ii,mean_pav_m,'b');
        plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
        
        glm_pav_mean.data(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=pavar_mean(these_ii);
        glm_pav_mean.perCorr(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=per_ii*ones(1,sum(these_ii));
        glm_pav_mean.event(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=evNo*ones(1,sum(these_ii));
        glm_pav_mean.exp(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=exp*ones(1,sum(these_ii));
        glm_pav_mii=glm_pav_mii+sum(these_ii);
        
        ii_stats_m_pav=ii_stats_m_pav+1;
        p_pav_mean_stats(ii_stats_m_pav).data=pavar_mean(these_ii);
        p_pav_mean_stats(ii_stats_m_pav).description=['Proficient S- Exp' num2str(exp)];
        
        %S+ naive
        evNo=1;
        per_ii=2;
        this_bar_ii=this_bar_ii+2;
        these_ii=(pavar_mean_evNo==evNo)&(pavar_mean_exp==exp)&(pavar_mean_per_ii==per_ii);
        mean_pav_m=mean(pavar_mean(these_ii));
        miCI = bootci(1000, {@mean, pavar_mean(these_ii)},'type','cper');
        hp2=bar(this_bar_ii,mean_pav_m,'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7]);
        plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
        
        glm_pav_mean.data(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=pavar_mean(these_ii);
        glm_pav_mean.perCorr(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=per_ii*ones(1,sum(these_ii));
        glm_pav_mean.event(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=evNo*ones(1,sum(these_ii));
        glm_pav_mean.exp(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=exp*ones(1,sum(these_ii));
        glm_pav_mii=glm_pav_mii+sum(these_ii);
        
        ii_stats_m_pav=ii_stats_m_pav+1;
        p_pav_mean_stats(ii_stats_m_pav).data=pavar_mean(these_ii);
        p_pav_mean_stats(ii_stats_m_pav).description=['Naive S+ Exp' num2str(exp)];
        
        %S+ proficient
        evNo=1;
        per_ii=1;
        this_bar_ii=this_bar_ii+1;
        these_ii=(pavar_mean_evNo==evNo)&(pavar_mean_exp==exp)&(pavar_mean_per_ii==per_ii);
        mean_pav_m=mean(pavar_mean(these_ii));
        miCI = bootci(1000, {@mean, pavar_mean(these_ii)},'type','cper');
        hp1=bar(this_bar_ii,mean_pav_m,'r');
        plot([this_bar_ii this_bar_ii],miCI,'-k','LineWidth',3);
        
        glm_pav_mean.data(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=pavar_mean(these_ii);
        glm_pav_mean.perCorr(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=per_ii*ones(1,sum(these_ii));
        glm_pav_mean.event(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=evNo*ones(1,sum(these_ii));
        glm_pav_mean.exp(glm_pav_mii+1:glm_pav_mii+sum(these_ii))=exp*ones(1,sum(these_ii));
        glm_pav_mii=glm_pav_mii+sum(these_ii);
        
        ii_stats_m_pav=ii_stats_m_pav+1;
        p_pav_mean_stats(ii_stats_m_pav).data=pavar_mean(these_ii);
        p_pav_mean_stats(ii_stats_m_pav).description=['Proficient S+ Exp' num2str(exp)];
        
        this_bar_ii=this_bar_ii+3;
    end
    
    legend([hp1 hp2 hp3 hp4],{'S+ Proficient','S+ Naive','S- Proficient','S- Naive'})
    xticks([3 11])
    xticklabels({'Exp1','Exp2'})
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
    
    %Perform the glm for pav_mean
    fprintf(1, ['\n\nglm for per experiment mean peak angle variance for Theta/' PACnames{PACii} '\n'])
    tbl = table(glm_pav_mean.data',glm_pav_mean.perCorr',glm_pav_mean.event',glm_pav_mean.exp',...
        'VariableNames',{'pav','proficiency','event','experiment'});
    mdl = fitglm(tbl,'pav~proficiency+experiment+event+proficiency*event*experiment'...
        ,'CategoricalVars',[2,3,4])
    
    %Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for per experiment mean peak angle variance for Theta/' PACnames{PACii} '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(p_pav_mean_stats);
        fprintf(1, '\n\n')
    catch
    end
    
    if PACii==3
       pffft=1; 
    end
end
pffft=1;
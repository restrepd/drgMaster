function drgSummaryBatchWavPwr
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

PathName{1}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 wavelet peak power/Daniel Batch analysis/';
FileName{1}='spm_LFP_wavephasepower04262019_APEB.mat';
odorPairName{1}='APEBloc1';

PathName{2}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 wavelet peak power/Daniel Batch analysis/';
FileName{2}='spm_LFP_wavephasepower04262019_EAPA.mat';
odorPairName{2}='EAPA1loc1';

PathName{3}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 wavelet peak power/Daniel Batch analysis/';
FileName{3}='spm_LFP_wavephasepower04262019_IAMO.mat';
odorPairName{3}='IAMOloc1';

PathName{4}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 wavelet peak power/Justin Batch analysis/';
FileName{4}='spm_LFP_wavephasepower04202019_EAPA.mat';
odorPairName{4}='EAPA2loc2';

PathName{5}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 wavelet peak power/Justin Batch analysis/';
FileName{5}='spm_LFP_wavephasepower04202019_IAAP.mat';
odorPairName{5}='IAAPloc2';

PathName{6}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 wavelet peak power/Justin Batch analysis/';
FileName{6}='spm_LFP_wavephasepower04202019_IAMO.mat';
odorPairName{6}='IAMOloc2';

figNo=1;

for PACii=1:3
    
    %Plot the bar graph for MI
    
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.3 .3 .6 .4])
    hold on
    
    try
        close(figNo+1)
    catch
    end
    hFig=figure(figNo+1);
    set(hFig, 'units','normalized','position',[.3 .3 .6 .4])
    hold on
    
    glm_prp_ii=0;
    glm_prp=[];
    ii_stats_prp=0;
    p_prp_stats=[];
    
    bar_ii=1;
    for fileNo=1:length(PathName)
        
        load([PathName{fileNo} FileName{fileNo}])
        
        
        for ii_rank=1:length(wave_power(PACii).input_datapt)
            
            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                figure(figNo)
            else
                figure(figNo+1)
            end
            per_ii=wave_power(PACii).input_datapt(ii_rank).per_ii;
            evNo=wave_power(PACii).input_datapt(ii_rank).evNo;
            if evNo==1
                if per_ii==1
                    %S+ Proficient
                    this_bar_ii=bar_ii+4;
                    mean_power=mean(wave_power(PACii).input_datapt(ii_rank).data);
                    mean_powerCI = bootci(1000, {@mean, wave_power(PACii).input_datapt(ii_rank).data},'type','cper');
                    hp1=bar(this_bar_ii,mean_power,'r');
                    plot([this_bar_ii this_bar_ii],mean_powerCI,'-k','LineWidth',3);
                    
                    these_prp=wave_power(PACii).input_datapt(ii_rank).data;
                    glm_prp.data(glm_prp_ii+1:glm_prp_ii+length(these_prp))=these_prp;
                    glm_prp.perCorr(glm_prp_ii+1:glm_prp_ii+length(these_prp))=1*ones(1,length(these_prp));
                    glm_prp.event(glm_prp_ii+1:glm_prp_ii+length(these_prp))=1*ones(1,length(these_prp));
                    if fileNo<=3
                        glm_prp.location(glm_prp_ii+1:glm_prp_ii+length(these_prp))=1*ones(1,length(these_prp));
                    else
                        glm_prp.location(glm_prp_ii+1:glm_prp_ii+length(these_prp))=2*ones(1,length(these_prp));
                    end
                    glm_prp.peak_trough(glm_prp_ii+1:glm_prp_ii+length(these_prp))=wave_power(PACii).input_datapt(ii_rank).peak_trough*ones(1,length(these_prp));
                    glm_prp_ii=glm_prp_ii+length(these_prp);
                    
                    ii_stats_prp=ii_stats_prp+1;
                    p_prp_stats(ii_stats_prp).data=these_prp;
                    if evNo==1
                        if per_ii==1
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Proficient CS+ peak ' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Proficient CS+ trough ' odorPairName{fileNo}];
                            end
                        else
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ peak ' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ trough' odorPairName{fileNo}];
                            end
                        end
                    else
                        if per_ii==1
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Proficient CS- peak' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Proficient CS- trough' odorPairName{fileNo}];
                            end
                        else
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ peak' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ torugh' odorPairName{fileNo}];
                            end
                        end
                    end
                    
                else
                    %S+ Naive
                    this_bar_ii=bar_ii+3;
                    mean_power=mean(wave_power(PACii).input_datapt(ii_rank).data);
                    mean_powerCI = bootci(1000, {@mean, wave_power(PACii).input_datapt(ii_rank).data},'type','cper');
                    hp2=bar(this_bar_ii,mean_power,'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7]);
                    plot([this_bar_ii this_bar_ii],mean_powerCI,'-k','LineWidth',3);
                    
                    these_prp=wave_power(PACii).input_datapt(ii_rank).data;
                    glm_prp.data(glm_prp_ii+1:glm_prp_ii+length(these_prp))=these_prp;
                    glm_prp.perCorr(glm_prp_ii+1:glm_prp_ii+length(these_prp))=2*ones(1,length(these_prp));
                    glm_prp.event(glm_prp_ii+1:glm_prp_ii+length(these_prp))=1*ones(1,length(these_prp));
                    if fileNo<=3
                        glm_prp.location(glm_prp_ii+1:glm_prp_ii+length(these_prp))=1*ones(1,length(these_prp));
                    else
                        glm_prp.location(glm_prp_ii+1:glm_prp_ii+length(these_prp))=2*ones(1,length(these_prp));
                    end
                    glm_prp.peak_trough(glm_prp_ii+1:glm_prp_ii+length(these_prp))=wave_power(PACii).input_datapt(ii_rank).peak_trough*ones(1,length(these_prp));
                    glm_prp_ii=glm_prp_ii+length(these_prp);
                    
                    ii_stats_prp=ii_stats_prp+1;
                    p_prp_stats(ii_stats_prp).data=these_prp;
                    if evNo==1
                        if per_ii==1
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Proficient CS+ peak ' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Proficient CS+ trough ' odorPairName{fileNo}];
                            end
                        else
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ peak ' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ trough' odorPairName{fileNo}];
                            end
                        end
                    else
                        if per_ii==1
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Proficient CS- peak' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Proficient CS- trough' odorPairName{fileNo}];
                            end
                        else
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ peak' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ torugh' odorPairName{fileNo}];
                            end
                        end
                    end
                end
            else
                if per_ii==1
                    %S- Proficient
                    this_bar_ii=bar_ii+1;
                    mean_power=mean(wave_power(PACii).input_datapt(ii_rank).data);
                    mean_powerCI = bootci(1000, {@mean, wave_power(PACii).input_datapt(ii_rank).data},'type','cper');
                    hp3=bar(this_bar_ii,mean_power,'b');
                    plot([this_bar_ii this_bar_ii],mean_powerCI,'-k','LineWidth',3);
                    
                    these_prp=wave_power(PACii).input_datapt(ii_rank).data;
                    glm_prp.data(glm_prp_ii+1:glm_prp_ii+length(these_prp))=these_prp;
                    glm_prp.perCorr(glm_prp_ii+1:glm_prp_ii+length(these_prp))=1*ones(1,length(these_prp));
                    glm_prp.event(glm_prp_ii+1:glm_prp_ii+length(these_prp))=2*ones(1,length(these_prp));
                    if fileNo<=3
                        glm_prp.location(glm_prp_ii+1:glm_prp_ii+length(these_prp))=1*ones(1,length(these_prp));
                    else
                        glm_prp.location(glm_prp_ii+1:glm_prp_ii+length(these_prp))=2*ones(1,length(these_prp));
                    end
                    glm_prp.peak_trough(glm_prp_ii+1:glm_prp_ii+length(these_prp))=wave_power(PACii).input_datapt(ii_rank).peak_trough*ones(1,length(these_prp));
                    glm_prp_ii=glm_prp_ii+length(these_prp);
                    
                    ii_stats_prp=ii_stats_prp+1;
                    p_prp_stats(ii_stats_prp).data=these_prp;
                    if evNo==1
                        if per_ii==1
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Proficient CS+ peak ' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Proficient CS+ trough ' odorPairName{fileNo}];
                            end
                        else
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ peak ' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ trough' odorPairName{fileNo}];
                            end
                        end
                    else
                        if per_ii==1
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Proficient CS- peak' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Proficient CS- trough' odorPairName{fileNo}];
                            end
                        else
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ peak' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ torugh' odorPairName{fileNo}];
                            end
                        end
                    end
                else
                    %S- Naive
                    this_bar_ii=bar_ii;
                    mean_power=mean(wave_power(PACii).input_datapt(ii_rank).data);
                    mean_powerCI = bootci(1000, {@mean, wave_power(PACii).input_datapt(ii_rank).data},'type','cper');
                    hp4=bar(this_bar_ii,mean_power,'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1]);
                    plot([this_bar_ii this_bar_ii],mean_powerCI,'-k','LineWidth',3);
                    
                    these_prp=wave_power(PACii).input_datapt(ii_rank).data;
                    glm_prp.data(glm_prp_ii+1:glm_prp_ii+length(these_prp))=these_prp;
                    glm_prp.perCorr(glm_prp_ii+1:glm_prp_ii+length(these_prp))=1*ones(1,length(these_prp));
                    glm_prp.event(glm_prp_ii+1:glm_prp_ii+length(these_prp))=1*ones(1,length(these_prp));
                    if fileNo<=3
                        glm_prp.location(glm_prp_ii+1:glm_prp_ii+length(these_prp))=2*ones(1,length(these_prp));
                    else
                        glm_prp.location(glm_prp_ii+1:glm_prp_ii+length(these_prp))=2*ones(1,length(these_prp));
                    end
                    glm_prp.peak_trough(glm_prp_ii+1:glm_prp_ii+length(these_prp))=wave_power(PACii).input_datapt(ii_rank).peak_trough*ones(1,length(these_prp));
                    glm_prp_ii=glm_prp_ii+length(these_prp);
                    
                    ii_stats_prp=ii_stats_prp+1;
                    p_prp_stats(ii_stats_prp).data=these_prp;
                    if evNo==1
                        if per_ii==1
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Proficient CS+ peak ' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Proficient CS+ trough ' odorPairName{fileNo}];
                            end
                        else
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ peak ' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ trough' odorPairName{fileNo}];
                            end
                        end
                    else
                        if per_ii==1
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Proficient CS- peak' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Proficient CS- trough' odorPairName{fileNo}];
                            end
                        else
                            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ peak' odorPairName{fileNo}];
                            else
                                p_prp_stats(ii_stats_prp).description=['Naive CS+ torugh' odorPairName{fileNo}];
                            end
                        end
                    end
                end
                
            end
            
            
        end
        
        bar_ii=bar_ii+7;
        
    end
    
    figure(figNo)
    ylim1=ylim;
    figure(figNo+1)
    ylim2=ylim;
    
    newylim=[min([ylim1(1) ylim2(1)]) max([ylim1(2) ylim2(2)])];
    
    for fNo=figNo:figNo+1
        figure(fNo)
        ylim(newylim)
        legend([hp1 hp2 hp3 hp4],{'S+ Proficient','S+ Naive','S- Proficient','S- Naive'})
        xticks([3.5 10.5 17.5 24.5 31.5 38.5])
        xticklabels({'APEB','EAPA1loc1','EAPA2loc2','IAAPloc2','IAMOloc1','IAMOloc2'})
        ylabel('Wavelet power')
        
        if fNo==figNo
            
            switch PACii
                case 1
                    title('Peak wavelet power Theta/Beta')
                case 2
                    title('Peak wavelet power Theta/Low Gamma')
                case 3
                    title('Peak wavelet power Theta/High Gamma')
            end
        else
            switch PACii
                case 1
                    title('Trough wavelet power Theta/Beta')
                case 2
                    title('Trough wavelet power Theta/Low Gamma')
                case 3
                    title('Trough wavelet power Theta/High Gamma')
            end
        end
    end
    
   figNo=figNo+2;
   
   %Perform the glm for mi
    fprintf(1, ['\n\nglm for phase-referenced power for Theta/' PACnames{PACii} '\n'])
    tbl = table(glm_prp.data',glm_prp.perCorr',glm_prp.event',glm_prp.location',glm_prp.peak_trough',...
        'VariableNames',{'prp','proficiency','event','location','peak_trough'});
    mdl = fitglm(tbl,'prp~proficiency+location+event+peak_trough+proficiency*event*location*peak_trough'...
        ,'CategoricalVars',[2,3 4,5])
     
    %Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for phase-referenced power per mouse for Theta/' PACnames{PACii} '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(p_prp_stats);
        fprintf(1, '\n\n')
    catch
    end
end




pffft=1;
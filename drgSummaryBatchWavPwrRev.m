function drgSummaryBatchWavPwrRev
%Performs the reversal analysis for Figure 4 of Losacco, Ramirez-Gordillo
%et al. eLife 2020


close all
clear all

% evTypeLabels={'S+','S-',};
% per_ii_labels={'Proficient','Naive'}

evTypeLabels={'Odor A','Odor B',};




PACnames{1}='Beta';
PACnames{2}='Low gamma';
PACnames{3}='High gamma';

%Dnaiel's reversal
PathName{1}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 Reversal/Daniel rev Batch analysis/';
FileName{1}='spm_LFP_wavephasepower04262019_EAPAfr.mat';
odorPairName{1}='EAPAexp1 fwd';
fwd_rev(1)=1;

PathName{2}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 Reversal/Daniel rev Batch analysis/';
FileName{2}='spm_LFP_wavephasepower04262019_EAPAr.mat';
odorPairName{2}='EAPAexp1 rev';
fwd_rev(2)=2;

%Justin's reversals
PathName{3}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 Reversal/Justin rev Batch analysis/';
FileName{3}='spm_LFP_wavephasepower04202019_EAPAfr.mat';
odorPairName{3}='EAPAexp2 fwd';
fwd_rev(3)=1;

PathName{4}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 Reversal/Justin rev Batch analysis/';
FileName{4}='spm_LFP_wavephasepower04202019_EAPAr.mat';
odorPairName{4}='EAPAexp2 rev';
fwd_rev(4)=2;

PathName{5}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 Reversal/Justin rev Batch analysis/';
FileName{5}='spm_LFP_wavephasepower04202019_IAAPfr.mat';
odorPairName{5}='IAAPexp2 fwd';
fwd_rev(5)=1;

PathName{6}='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 Reversal/Justin rev Batch analysis/';
FileName{6}='spm_LFP_wavephasepower04202019_IAAPr.mat';
odorPairName{6}='IAAPexp2 rev';
fwd_rev(6)=2;



figNo=1;

for PACii=[1 3]
    
    %Plot the bar graph for peak forward vs reversed
 
    %Forward
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.3 .3 .6 .4])
    hold on
    
    %Reversed
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
        
        odorPNo=floor((fileNo+1)/2);
        bar_ii=(odorPNo-1)*3;
        
        for ii_rank=1:length(wave_power(PACii).input_datapt)
            
            %Plot peak only
            if wave_power(PACii).input_datapt(ii_rank).peak_trough==1
                
                if fwd_rev(fileNo)==1
                    figure(figNo)
                else
                    figure(figNo+1)
                end
                
                
                
                per_ii=wave_power(PACii).input_datapt(ii_rank).per_ii;
                evNo=wave_power(PACii).input_datapt(ii_rank).evNo;
                
                if fwd_rev(fileNo)==1
                    %These are forward, odor A is S+
                    if evNo==1
                        if per_ii==1
                            %S+ Proficient Odor A
                            this_bar_ii=bar_ii+2;
                            mean_power=mean(wave_power(PACii).input_datapt(ii_rank).data);
                            mean_powerCI = bootci(1000, {@mean, wave_power(PACii).input_datapt(ii_rank).data},'type','cper');
                            hp1=bar(this_bar_ii,mean_power,'r');
                            plot([this_bar_ii this_bar_ii],mean_powerCI,'-k','LineWidth',3);
                            
                            these_prp=wave_power(PACii).input_datapt(ii_rank).data;
                            glm_prp.data(glm_prp_ii+1:glm_prp_ii+length(these_prp))=these_prp;
                            glm_prp.odorant(glm_prp_ii+1:glm_prp_ii+length(these_prp))=1*ones(1,length(these_prp));
                            glm_prp.fwdvsrev(glm_prp_ii+1:glm_prp_ii+length(these_prp))=1*ones(1,length(these_prp));
                            glm_prp_ii=glm_prp_ii+length(these_prp);
                            
                            ii_stats_prp=ii_stats_prp+1;
                            p_prp_stats(ii_stats_prp).data=these_prp;
                            p_prp_stats(ii_stats_prp).description=['Forward odor A S+ ' odorPairName{fileNo}];
                            
                        end
                    else
                        if per_ii==1
                            %S- Proficient, Odor B 
                            this_bar_ii=bar_ii+1;
                            mean_power=mean(wave_power(PACii).input_datapt(ii_rank).data);
                            mean_powerCI = bootci(1000, {@mean, wave_power(PACii).input_datapt(ii_rank).data},'type','cper');
                            hp2=bar(this_bar_ii,mean_power,'b');
                            plot([this_bar_ii this_bar_ii],mean_powerCI,'-k','LineWidth',3);
                            
                            these_prp=wave_power(PACii).input_datapt(ii_rank).data;
                            glm_prp.data(glm_prp_ii+1:glm_prp_ii+length(these_prp))=these_prp;
                            glm_prp.odorant(glm_prp_ii+1:glm_prp_ii+length(these_prp))=2*ones(1,length(these_prp));
                            glm_prp.fwdvsrev(glm_prp_ii+1:glm_prp_ii+length(these_prp))=1*ones(1,length(these_prp));
                            glm_prp_ii=glm_prp_ii+length(these_prp);
                            
                            ii_stats_prp=ii_stats_prp+1;
                            p_prp_stats(ii_stats_prp).data=these_prp;
                            p_prp_stats(ii_stats_prp).description=['Forward odor B S- ' odorPairName{fileNo}];
                            
                        end
                        
                    end
                    
                else
                    %Reversed
                    if evNo==1
                        if per_ii==1
                            %S+ Proficient, Odor B
                            this_bar_ii=bar_ii+1;
                            mean_power=mean(wave_power(PACii).input_datapt(ii_rank).data);
                            mean_powerCI = bootci(1000, {@mean, wave_power(PACii).input_datapt(ii_rank).data},'type','cper');
                            hp3=bar(this_bar_ii,mean_power,'b');
                            plot([this_bar_ii this_bar_ii],mean_powerCI,'-k','LineWidth',3);
                            
                            these_prp=wave_power(PACii).input_datapt(ii_rank).data;
                            glm_prp.data(glm_prp_ii+1:glm_prp_ii+length(these_prp))=these_prp;
                            glm_prp.odorant(glm_prp_ii+1:glm_prp_ii+length(these_prp))=2*ones(1,length(these_prp));
                            glm_prp.fwdvsrev(glm_prp_ii+1:glm_prp_ii+length(these_prp))=2*ones(1,length(these_prp));
                            glm_prp_ii=glm_prp_ii+length(these_prp);
                            
                            ii_stats_prp=ii_stats_prp+1;
                            p_prp_stats(ii_stats_prp).data=these_prp;
                            p_prp_stats(ii_stats_prp).description=['Reversed odor B S+ ' odorPairName{fileNo}];
                            
                            

                        end
                    else
                        if per_ii==1
                            %S- Proficient Odor A
                            this_bar_ii=bar_ii+2;
                            mean_power=mean(wave_power(PACii).input_datapt(ii_rank).data);
                            mean_powerCI = bootci(1000, {@mean, wave_power(PACii).input_datapt(ii_rank).data},'type','cper');
                            hp4=bar(this_bar_ii,mean_power,'r');
                            plot([this_bar_ii this_bar_ii],mean_powerCI,'-k','LineWidth',3);
                            
                            these_prp=wave_power(PACii).input_datapt(ii_rank).data;
                            glm_prp.data(glm_prp_ii+1:glm_prp_ii+length(these_prp))=these_prp;
                            glm_prp.odorant(glm_prp_ii+1:glm_prp_ii+length(these_prp))=1*ones(1,length(these_prp));
                            glm_prp.fwdvsrev(glm_prp_ii+1:glm_prp_ii+length(these_prp))=2*ones(1,length(these_prp));
                            glm_prp_ii=glm_prp_ii+length(these_prp);
                            
                            ii_stats_prp=ii_stats_prp+1;
                            p_prp_stats(ii_stats_prp).data=these_prp;
                            p_prp_stats(ii_stats_prp).description=['Reversed odor A S- ' odorPairName{fileNo}];
                            

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
        if fNo==figNo
            legend([hp1 hp2],{'Odor A','Odor B'})
        end
        xticks([1.5 4.5 7.5])
        xticklabels({'EAPAexp1','EAPAexp2','IAAPexp2'})
        ylabel('Wavelet power')
        
        if fNo==figNo
            
            switch PACii
                case 1
                    title('Peak wavelet power Theta/Beta, forward')
                case 2
                    title('Peak wavelet power Theta/Low Gamma, forward')
                case 3
                    title('Peak wavelet power Theta/High Gamma, forward')
            end
        else
            switch PACii
                   case 1
                    title('Peak wavelet power Theta/Beta, reverse')
                case 2
                    title('Peak wavelet power Theta/Low Gamma, reverse')
                case 3
                    title('Peak wavelet power Theta/High Gamma, reverse')
            end
        end
    end
    
    figNo=figNo+2;
    
    %Perform the glm for mi
    fprintf(1, ['\n\nglm for phase-referenced power for Theta/' PACnames{PACii} '\n'])
    tbl = table(glm_prp.data',glm_prp.fwdvsrev',glm_prp.odorant',...
        'VariableNames',{'prp','fwd_rev','odorant'});
    mdl = fitglm(tbl,'prp~fwd_rev+odorant+fwd_rev*odorant'...
        ,'CategoricalVars',[2,3])
    
    %Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for phase-referenced power per mouse for Theta/' PACnames{PACii} '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(p_prp_stats);
        fprintf(1, '\n\n')
    catch
    end
end




pffft=1;
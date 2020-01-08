close all
clear all

no_odor_pairs=6;
 
odor_pair_no=0;

%Daniel's
outFileName='dim_Discriminant_spmc_discriminantolfac_all_PCA_LFP_aceto42919.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 LDA/APEBDR/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 LDA\APEBDR\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out2=handles_out2;
odor_pair_label{1}='APEBloc1';
location(1)=1;
 
outFileName='dim_Discriminant_spmc_discriminantolfac_all_PCA_LFP_EAPA42919.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 LDA/EAPADR/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 LDA\EAPADR\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out2=handles_out2;
odor_pair_label{2}='EAPAloc1';
location(2)=1;

outFileName='dim_Discriminant_spmc_discriminantolfac_all_PCA_LFP_iso42919.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 LDA/IAMODR/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 LDA\IAMODR\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out2=handles_out2;
odor_pair_label{3}='IAMOloc1';
location(3)=1;

%Justin's
outFileName='dim_Discriminant_spm_discriminant_LFP_wavephase_05172019_IAAP.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 LDA/IAAPJL/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 LDA\IAAPJL\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out2=handles_out2;
odor_pair_label{4}='IAAPloc2';
location(4)=2;

outFileName='dim_Discriminant_spm_discriminant_LFP_wavephase_04252019_EAPA.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 LDA/EAPAJL/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 LDA\EAPAJL\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out2=handles_out2;
odor_pair_label{5}='EAPAloc2';
location(5)=2;

outFileName='dim_Discriminant_spm_discriminant_LFP_wavephase_04252019_IsoAA_mo.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 5 LDA/IAMOJL/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 5 LDA\IAMOJL\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_dim_per_odor_pair(odor_pair_no).handles_out2=handles_out2;
odor_pair_label{6}='IAMOloc2';
location(6)=2;

PACnames{1}='Beta';
PACnames{2}='Low gamma';
PACnames{3}='High gamma';

window_label{1}='Pre-odor';
window_label{2}='Odor';


%Plot the relationship of the detection times
figNo=0;

%Plot decision time reationship using mean decision lick times per mouse
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

bar_ii=1;

figNo=0;


for PACii=1:3
    glm_ii=0;
    ii_stats=0;
    glm_dtime=[];
    p_dtime_stats=[];
    for winNo=2:-1:1
        
        %Plot dimensionality for odor window for proficient mice
        
        groupNo=1;
        
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        set(hFig, 'units','normalized','position',[.1 .4 .75 .47])
        hold on
        
        groupNo=1;
        bar_ii=0;
        
        for ii_op=1:odor_pair_no
            
            
            %Naive trough
            percent_correct_ii=2;
            peak_trough=0;
            
            these_groups=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.group;
            these_windows=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.odor_window;
            these_perCorr=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.perCorr;
            these_peaks=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.peak;
            
            these_dims=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.data(...
                (these_groups==groupNo)&(these_windows==winNo)&(these_perCorr==percent_correct_ii)...
                &(these_peaks==peak_trough));
            trough_dims=these_dims;
            
            mean_dims=mean(these_dims)';
            CIdims = bootci(1000, {@mean, these_dims})';
            p1=bar(bar_ii,mean_dims,'EdgeColor',[0.7 0.7 1],'FaceColor',[0.7 0.7 1]);
            plot(bar_ii,mean_dims,'ok','MarkerFaceColor','k','MarkerSize',5)
            plot(bar_ii*ones(1,length(these_dims)),these_dims,'ok')
            plot([bar_ii bar_ii],CIdims,'-k','LineWidth',3)
            
            glm_dims.data(glm_ii+1:glm_ii+length(these_dims))=these_dims;
            glm_dims.perCorr(glm_ii+1:glm_ii+length(these_dims))=percent_correct_ii;
            glm_dims.trough_peak(glm_ii+1:glm_ii+length(these_dims))=peak_trough;
            glm_dims.location(glm_ii+1:glm_ii+length(these_dims))=location(ii_op);
            glm_dims.window(glm_ii+1:glm_ii+length(these_dims))=winNo;
            glm_ii=glm_ii+length(these_dims);
            
            ii_stats=ii_stats+1;
            p_dims_stats(ii_stats).data=these_dims;
            p_dims_stats(ii_stats).description=['Naive trough ' odor_pair_label{ii_op}];
            p_dims_stats(ii_stats).perCorr=percent_correct_ii;
            p_dims_stats(ii_stats).trough_peak=peak_trough;
            
            %Naive peak
            bar_ii=bar_ii+1;
            peak_trough=1;
            
            these_groups=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.group;
            these_windows=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.odor_window;
            these_perCorr=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.perCorr;
            these_peaks=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.peak;
            
            these_dims=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.data(...
                (these_groups==groupNo)&(these_windows==winNo)&(these_perCorr==percent_correct_ii)...
                &(these_peaks==peak_trough));
            peak_dims=these_dims;
            
            mean_dims=mean(these_dims)';
            CIdims = bootci(1000, {@mean, these_dims})';
            p2=bar(bar_ii,mean_dims,'EdgeColor',[1 0.7 0.7],'FaceColor',[1 0.7 0.7]);
            plot(bar_ii,mean_dims,'ok','MarkerFaceColor','k','MarkerSize',5)
            plot(bar_ii*ones(1,length(these_dims)),these_dims,'ok')
            plot([bar_ii bar_ii],CIdims,'-k','LineWidth',3)
            
            glm_dims.data(glm_ii+1:glm_ii+length(these_dims))=these_dims;
            glm_dims.perCorr(glm_ii+1:glm_ii+length(these_dims))=percent_correct_ii;
            glm_dims.trough_peak(glm_ii+1:glm_ii+length(these_dims))=peak_trough;
            glm_dims.location(glm_ii+1:glm_ii+length(these_dims))=location(ii_op);
            glm_dims.window(glm_ii+1:glm_ii+length(these_dims))=winNo;
            glm_ii=glm_ii+length(these_dims);
            
            ii_stats=ii_stats+1;
            p_dims_stats(ii_stats).data=these_dims;
            p_dims_stats(ii_stats).description=['Naive peak ' odor_pair_label{ii_op}];
            p_dims_stats(ii_stats).perCorr=percent_correct_ii;
            p_dims_stats(ii_stats).trough_peak=peak_trough;
            
            for ii_dim=1:length(these_dims)
                plot([bar_ii-1 bar_ii], [trough_dims(ii_dim) peak_dims(ii_dim)],'-k')
            end
            
            bar_ii=bar_ii+2;
            %Proficient trough
            percent_correct_ii=1;
            peak_trough=0;
            
            these_groups=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.group;
            these_windows=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.odor_window;
            these_perCorr=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.perCorr;
            these_peaks=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.peak;
            
            these_dims=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.data(...
                (these_groups==groupNo)&(these_windows==winNo)&(these_perCorr==percent_correct_ii)...
                &(these_peaks==peak_trough));
            trough_dims=these_dims;
            
            mean_dims=mean(these_dims)';
            CIdims = bootci(1000, {@mean, these_dims})';
            p3=bar(bar_ii,mean_dims,'EdgeColor','b','FaceColor','b');
            plot(bar_ii,mean_dims,'ok','MarkerFaceColor','k','MarkerSize',5)
            plot(bar_ii*ones(1,length(these_dims)),these_dims,'ok')
            plot([bar_ii bar_ii],CIdims,'-k','LineWidth',3)
            
            glm_dims.data(glm_ii+1:glm_ii+length(these_dims))=these_dims;
            glm_dims.perCorr(glm_ii+1:glm_ii+length(these_dims))=percent_correct_ii;
            glm_dims.trough_peak(glm_ii+1:glm_ii+length(these_dims))=peak_trough;
            glm_dims.location(glm_ii+1:glm_ii+length(these_dims))=location(ii_op);
            glm_dims.window(glm_ii+1:glm_ii+length(these_dims))=winNo;
            glm_ii=glm_ii+length(these_dims);
            
            ii_stats=ii_stats+1;
            p_dims_stats(ii_stats).data=these_dims;
            p_dims_stats(ii_stats).description=['Proficient trough ' odor_pair_label{ii_op}];
            p_dims_stats(ii_stats).perCorr=percent_correct_ii;
            p_dims_stats(ii_stats).trough_peak=peak_trough;
            
            %Proficient peak
            bar_ii=bar_ii+1;
            peak_trough=1;
            
            these_groups=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.group;
            these_windows=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.odor_window;
            these_perCorr=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.perCorr;
            these_peaks=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.peak;
            
            these_dims=handles_dim_per_odor_pair(ii_op).handles_out2.PACii(PACii).glm_dim.data(...
                (these_groups==groupNo)&(these_windows==winNo)&(these_perCorr==percent_correct_ii)...
                &(these_peaks==peak_trough));
            peak_dims=these_dims;
            
            mean_dims=mean(these_dims)';
            CIdims = bootci(1000, {@mean, these_dims})';
            p4=bar(bar_ii,mean_dims,'EdgeColor','r','FaceColor','r');
            plot(bar_ii,mean_dims,'ok','MarkerFaceColor','k','MarkerSize',5)
            plot(bar_ii*ones(1,length(these_dims)),these_dims,'ok')
            plot([bar_ii bar_ii],CIdims,'-k','LineWidth',3)
            
            glm_dims.data(glm_ii+1:glm_ii+length(these_dims))=these_dims;
            glm_dims.perCorr(glm_ii+1:glm_ii+length(these_dims))=percent_correct_ii;
            glm_dims.trough_peak(glm_ii+1:glm_ii+length(these_dims))=peak_trough;
            glm_dims.location(glm_ii+1:glm_ii+length(these_dims))=location(ii_op);
            glm_dims.window(glm_ii+1:glm_ii+length(these_dims))=winNo;
            glm_ii=glm_ii+length(these_dims);
            
            ii_stats=ii_stats+1;
            p_dims_stats(ii_stats).data=these_dims;
            p_dims_stats(ii_stats).description=['Proficient peak ' odor_pair_label{ii_op}];
            p_dims_stats(ii_stats).perCorr=percent_correct_ii;
            p_dims_stats(ii_stats).trough_peak=peak_trough;
            
            for ii_dim=1:length(these_dims)
                plot([bar_ii-1 bar_ii], [trough_dims(ii_dim) peak_dims(ii_dim)],'-k')
            end
            
            bar_ii=bar_ii+4;
        end
        
        title(['[Dimensionality per mouse for theta/' PACnames{PACii} ', ' window_label{winNo} ' window'])
        ylabel('Dimensionality')
        legend([p1 p2 p3 p4],{'Naive, Trough','Naive, Peak','Proficient, Trough','Proficient, Peak'})
        xticks([2 10 18 26 34 42])
        xticklabels({'APEBloc1','EAPAloc1','IAMOloc1','EAPAloc2','IAAPloc2','IAMOloc2'})
        ylim([0 6])
        
    end
    
    %Perform the glm
    fprintf(1, ['\n\nglm for dimensionality per mouse for Theta/' PACnames{PACii} '\n'])
    tbl = table(glm_dims.data',glm_dims.perCorr',glm_dims.trough_peak',glm_dims.location',glm_dims.window',...
        'VariableNames',{'dimensionality','proficiency','peak_trough','location','window'});
    mdl = fitglm(tbl,'dimensionality~proficiency+location+peak_trough+window'...
        ,'CategoricalVars',[2,3 4])
    
    %Do ranksum/t test
    fprintf(1, ['\n\nRanksum or t-test p values for dimensionality per mouse for Theta/' PACnames{PACii} '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(p_dims_stats);
        fprintf(1, '\n\n')
    catch
    end
end






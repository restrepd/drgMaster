close all
clear all

no_odor_pairs=6;
 
odor_pair_no=0;

%Justin's
outFileName='LDA_all_m_spm_discriminant_LFP_wavephase_05172019_IAAP.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 LDA/IAAPJL/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 4 LDA\IAAPJL\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_lda_per_odor_pair(odor_pair_no).handles_outp=handles_outp;
odorPairName{1}='IAAPloc2';
location(1)=2;

outFileName='LDA_all_m_spm_discriminant_LFP_wavephase_04252019_EAPA.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 LDA/EAPAJL/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 4 LDA\EAPAJL\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_lda_per_odor_pair(odor_pair_no).handles_outp=handles_outp;
odorPairName{2}='EAPAloc2';
location(2)=2;


%Daniel's
outFileName='LDA_all_m_spmc_discriminantolfac_all_PCA_LFP_aceto42919.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 LDA/APEBDR/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 4 LDA\APEBDR\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_lda_per_odor_pair(odor_pair_no).handles_outp=handles_outp;
odorPairName{3}='APEBloc1';
location(3)=1;

outFileName='LDA_all_m_spmc_discriminantolfac_all_PCA_LFP_EAPA42919.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 LDA/EAPADR/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 4 LDA\EAPADR\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_lda_per_odor_pair(odor_pair_no).handles_outp=handles_outp;
odorPairName{4}='EAPAloc1';
location(4)=1;

outFileName='LDA_all_m_spmc_discriminantolfac_all_PCA_LFP_iso42919.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 LDA/IAMODR/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 4 LDA\IAMODR\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_lda_per_odor_pair(odor_pair_no).handles_outp=handles_outp;
odorPairName{5}='IAMOloc1';
location(5)=1;

outFileName='LDA_all_m_spm_discriminant_LFP_wavephase_04252019_IsoAA_mo.mat';
try
    outPathName='/Users/restrepd/Dropbox/SPM_manuscript_Justin_Daniel/Fig 4 LDA/IAMOJL/';
    load([outPathName outFileName])
catch
    outPathName='C:\Users\Justin Losacco\Dropbox\Restrepo_Lab\Papers\JL paper\Identity paper\SPM manuscript  Justin-Daniel\Fig 4 LDA\IAMOJL\';
    load([outPathName outFileName])
end
odor_pair_no=odor_pair_no+1;
handles_lda_per_odor_pair(odor_pair_no).handles_outp=handles_outp;
odorPairName{6}='IAMOloc2';
location(6)=2;

PACnames{1}='Beta';
PACnames{2}='Low gamma';
PACnames{3}='High gamma';

%Plot the relationship of the detection times
figNo=0;

%Do decision time reationship for pooled mice proficient 
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

bar_ii=1;
glm_ii=0;
ii_stats=0;
glm_dtime=[];
p_dtime_stats=[];

%mean_lick_t_detect is the lick decision time computed per mouse
%t_detect_licks_all is the lick decision time computed with all mice pooled

percent_correct_ii=1;

groupNo=1;

for PACii=[1 3]
    
    
    
    %Plot decision time reationship between per mouse and all mouse decision
    %times
    
    lick_t_detect=[];
    t_detect_licks_all=[];
    
    %Licks
    for ii_odor_pair=1:no_odor_pairs
        lick_t_detect(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_lick_t_detect;
    end
    for ii_odor_pair=1:no_odor_pairs
        t_detect_licks_all(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).t_detect_licks_all;
    end
    
    mean_lick_t_detect=mean(lick_t_detect)';
    CIlicks = bootci(1000, {@mean, lick_t_detect})';
    p1=bar(bar_ii,mean_lick_t_detect,'EdgeColor','k','FaceColor',[0.7 1 0.7]);
    plot(bar_ii,mean_lick_t_detect,'ok','MarkerFaceColor','k','MarkerSize',5)
    plot(bar_ii*ones(1,length(lick_t_detect)),lick_t_detect,'ok')
    plot([bar_ii bar_ii],CIlicks,'-k','LineWidth',3)
    
    glm_dtime.data(glm_ii+1:glm_ii+length(lick_t_detect))=lick_t_detect;
    glm_dtime.PACii(glm_ii+1:glm_ii+length(lick_t_detect))=PACii;
    glm_dtime.per_vs_all(glm_ii+1:glm_ii+length(lick_t_detect))=1;
    glm_dtime.location(glm_ii+1:glm_ii+length(lick_t_detect))=location;
    glm_ii=glm_ii+length(lick_t_detect);
    
    ii_stats=ii_stats+1;
    p_dtime_stats(ii_stats).data=lick_t_detect;
    p_dtime_stats(ii_stats).description=[PACnames{PACii} ' lick mean per mouse'];
    p_dtime_stats(ii_stats).PACii=PACii;
    p_dtime_stats(ii_stats).lick_peak=1;
    
    bar_ii=bar_ii+1;
    mean_t_detect_licks_all=mean(t_detect_licks_all)';
    CIlicks_all = bootci(1000, {@mean, t_detect_licks_all})';
    p2=bar(bar_ii,mean_t_detect_licks_all,'EdgeColor','g','FaceColor','g');
    plot(bar_ii,mean_t_detect_licks_all,'ok','MarkerFaceColor','k','MarkerSize',5)
    plot(bar_ii*ones(1,length(t_detect_licks_all)),t_detect_licks_all,'ok')
    plot([bar_ii bar_ii],CIlicks_all,'-k','LineWidth',3)
    
    glm_dtime.data(glm_ii+1:glm_ii+length(t_detect_licks_all))=t_detect_licks_all;
    glm_dtime.PACii(glm_ii+1:glm_ii+length(t_detect_licks_all))=PACii;
    glm_dtime.per_vs_all(glm_ii+1:glm_ii+length(t_detect_licks_all))=2;
    glm_dtime.location(glm_ii+1:glm_ii+length(lick_t_detect))=location;
    glm_ii=glm_ii+length(t_detect_licks_all);
    
    ii_stats=ii_stats+1;
    p_dtime_stats(ii_stats).data=t_detect_licks_all;
    p_dtime_stats(ii_stats).description=[PACnames{PACii} ' lick pooled mice'];
    p_dtime_stats(ii_stats).PACii=PACii;
    p_dtime_stats(ii_stats).lick_peak=1;
    
    
    
    for ii_odor_pair=1:no_odor_pairs
        plot([bar_ii-1 bar_ii], [lick_t_detect(ii_odor_pair) t_detect_licks_all(ii_odor_pair)],'-k')
    end
    
    bar_ii=bar_ii+2;
end

title('Lick decision time for proficient mice')
ylabel('Decision time (sec)')
legend([p1 p2],{'Mean per mouse','Pooled mice'})
xticks([1.5 4.5])
xticklabels({'Beta','High gamma'})
ylim([0 5])

%Perform the glm
fprintf(1, ['\n\nglm for lick decision time for pooled mice vs per mouse (proficient) \n'])
tbl = table(glm_dtime.data',glm_dtime.PACii',glm_dtime.per_vs_all',glm_dtime.location',...
    'VariableNames',{'detection_time','PAC','per_vs_all','location'});
mdl = fitglm(tbl,'detection_time~PAC+per_vs_all+location'...
    ,'CategoricalVars',[2,3])

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for lick decision times for pooled mice vs per mouse (proficient) \n'])
try
    [output_data] = drgMutiRanksumorTtest(p_dtime_stats);
    fprintf(1, '\n\n')
catch
end

figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

bar_ii=1;
glm_ii=0;
ii_stats=0;
glm_dtime=[];
p_dtime_stats=[];

for PACii=[1 3]
    %Peak
    bar_ii=bar_ii+1;
    for ii_odor_pair=1:no_odor_pairs
        peak_t_detect_all(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).t_detect_peak_all;
    end
    for ii_odor_pair=1:no_odor_pairs
        peak_t_detect(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_peak_t_detect;
    end
    
    mean_peak_t_detect=mean(peak_t_detect)';
    handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_peak_t_detect=mean_peak_t_detect;
    CIpeak = bootci(1000, {@mean, peak_t_detect})';
    p1=bar(bar_ii,mean_peak_t_detect,'EdgeColor','k','FaceColor',[1 0.7 0.7]);
    plot(bar_ii,mean_peak_t_detect,'ok','MarkerFaceColor','k','MarkerSize',5)
    plot(bar_ii*ones(1,length(peak_t_detect)),peak_t_detect,'ok')
    plot([bar_ii bar_ii],CIpeak,'-k','LineWidth',3)
    
    glm_dtime.data(glm_ii+1:glm_ii+length(peak_t_detect))=peak_t_detect;
    glm_dtime.PACii(glm_ii+1:glm_ii+length(peak_t_detect))=PACii;
    glm_dtime.per_vs_all(glm_ii+1:glm_ii+length(t_detect_licks_all))=1;
    glm_dtime.location(glm_ii+1:glm_ii+length(lick_t_detect))=location;
    glm_ii=glm_ii+length(peak_t_detect);
   
    ii_stats=ii_stats+1;
    p_dtime_stats(ii_stats).data=peak_t_detect;
    p_dtime_stats(ii_stats).description=[PACnames{PACii} ' peak mean per mouse'];
    p_dtime_stats(ii_stats).PACii=PACii;
    p_dtime_stats(ii_stats).lick_peak=2;
    
    bar_ii=bar_ii+1;
    mean_peak_t_detect_all=mean(peak_t_detect_all)';
    handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_peak_t_detect_all=mean_peak_t_detect_all;
    CIpeak = bootci(1000, {@mean, peak_t_detect_all})';
    p3=bar(bar_ii,mean_peak_t_detect_all,'r');
    plot(bar_ii,mean_peak_t_detect_all,'ok','MarkerFaceColor','k','MarkerSize',5)
    plot(bar_ii*ones(1,length(peak_t_detect_all)),peak_t_detect_all,'ok')
    plot([bar_ii bar_ii],CIpeak,'-k','LineWidth',3)
    
    glm_dtime.data(glm_ii+1:glm_ii+length(peak_t_detect_all))=peak_t_detect_all;
    glm_dtime.PACii(glm_ii+1:glm_ii+length(peak_t_detect_all))=PACii;
    glm_dtime.per_vs_all(glm_ii+1:glm_ii+length(t_detect_licks_all))=2;
    glm_dtime.location(glm_ii+1:glm_ii+length(lick_t_detect))=location;
    glm_ii=glm_ii+length(peak_t_detect_all);
    
    ii_stats=ii_stats+1;
    p_dtime_stats(ii_stats).data=peak_t_detect_all;
    p_dtime_stats(ii_stats).description=[PACnames{PACii} ' peak pooled mice'];
    p_dtime_stats(ii_stats).PACii=PACii;
    p_dtime_stats(ii_stats).lick_peak=2;
    
    for ii_odor_pair=1:no_odor_pairs
        plot([bar_ii-1 bar_ii], [peak_t_detect(ii_odor_pair) peak_t_detect_all(ii_odor_pair)],'-k')
    end
    
    bar_ii=bar_ii+2;
    
end

title('Peak decision making time for proficient mice')
ylabel('Decision time (sec)')
legend([p1 p3],{'Mean per mouse','Pooled mice'})
xticks([1.5 4.5])
xticklabels({'Beta','High gamma'})
ylim([0 5])

%Perform the glm
fprintf(1, ['\n\nglm for peak decision time for pooled mice vs per mouse mice (proficient) \n'])
tbl = table(glm_dtime.data',glm_dtime.PACii',glm_dtime.per_vs_all',glm_dtime.location',...
    'VariableNames',{'detection_time','PAC','per_vs_all','location'});
mdl = fitglm(tbl,'detection_time~PAC+per_vs_all+location'...
    ,'CategoricalVars',[2,3])

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for peak decision times for pooled mice vs per mouse (proficient) \n'])
try
    [output_data] = drgMutiRanksumorTtest(p_dtime_stats);
    fprintf(1, '\n\n')
catch
end


%Now do Naive
%Plot decision time reationship using mean decision lick times per mouse
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

bar_ii=1;
glm_ii=0;
ii_stats=0;
glm_dtime=[];
p_dtime_stats=[];

%mean_lick_t_detect is the lick decision time computed per mouse
%t_detect_licks_all is the lick decision time computed with all mice pooled

percent_correct_ii=2;

groupNo=1;

for PACii=[1 3]
    
    
    
    %Plot decision time reationship between per mouse and all mouse decision
    %times
    
    %Licks
    for ii_odor_pair=1:no_odor_pairs
        lick_t_detect(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_lick_t_detect;
    end
    for ii_odor_pair=1:no_odor_pairs
        t_detect_licks_all(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).t_detect_licks_all;
    end
    
    mean_lick_t_detect=mean(lick_t_detect)';
    CIlicks = bootci(1000, {@mean, lick_t_detect})';
    p1=bar(bar_ii,mean_lick_t_detect,'EdgeColor','k','FaceColor',[0.7 1 0.7]);
    plot(bar_ii,mean_lick_t_detect,'ok','MarkerFaceColor','k','MarkerSize',5)
    plot(bar_ii*ones(1,length(lick_t_detect)),lick_t_detect,'ok')
    plot([bar_ii bar_ii],CIlicks,'-k','LineWidth',3)
    
    glm_dtime.data(glm_ii+1:glm_ii+length(lick_t_detect))=lick_t_detect;
    glm_dtime.PACii(glm_ii+1:glm_ii+length(lick_t_detect))=PACii;
    glm_dtime.per_vs_all(glm_ii+1:glm_ii+length(lick_t_detect))=1;
    glm_dtime.location(glm_ii+1:glm_ii+length(lick_t_detect))=location;
    glm_ii=glm_ii+length(lick_t_detect);
    
    ii_stats=ii_stats+1;
    p_dtime_stats(ii_stats).data=lick_t_detect;
    p_dtime_stats(ii_stats).description=[PACnames{PACii} ' lick mean per mouse'];
    p_dtime_stats(ii_stats).PACii=PACii;
    p_dtime_stats(ii_stats).lick_peak=1;
    
    bar_ii=bar_ii+1;
    mean_t_detect_licks_all=mean(t_detect_licks_all)';
    CIlicks_all = bootci(1000, {@mean, t_detect_licks_all})';
    p2=bar(bar_ii,mean_t_detect_licks_all,'EdgeColor','g','FaceColor','g');
    plot(bar_ii,mean_t_detect_licks_all,'ok','MarkerFaceColor','k','MarkerSize',5)
    plot(bar_ii*ones(1,length(t_detect_licks_all)),t_detect_licks_all,'ok')
    plot([bar_ii bar_ii],CIlicks_all,'-k','LineWidth',3)
    
    glm_dtime.data(glm_ii+1:glm_ii+length(t_detect_licks_all))=t_detect_licks_all;
    glm_dtime.PACii(glm_ii+1:glm_ii+length(t_detect_licks_all))=PACii;
    glm_dtime.per_vs_all(glm_ii+1:glm_ii+length(t_detect_licks_all))=2;
    glm_dtime.location(glm_ii+1:glm_ii+length(lick_t_detect))=location;
    glm_ii=glm_ii+length(t_detect_licks_all);
    
    ii_stats=ii_stats+1;
    p_dtime_stats(ii_stats).data=t_detect_licks_all;
    p_dtime_stats(ii_stats).description=[PACnames{PACii} ' lick pooled mice'];
    p_dtime_stats(ii_stats).PACii=PACii;
    p_dtime_stats(ii_stats).lick_peak=1;
    
    
    for ii_odor_pair=1:no_odor_pairs
        plot([bar_ii-1 bar_ii], [lick_t_detect(ii_odor_pair) t_detect_licks_all(ii_odor_pair)],'-k')
    end
    
    bar_ii=bar_ii+2;
end

title('Lick decision time for naive mice')
ylabel('Decision time (sec)')
legend([p1 p2],{'Mean per mouse','Pooled mice'})
xticks([1.5 4.5])
xticklabels({'Beta','High gamma'})
ylim([0 5])

%Perform the glm
fprintf(1, ['\n\nglm for lick decision time for pooled mice vs per mouse  (naive) \n'])
tbl = table(glm_dtime.data',glm_dtime.PACii',glm_dtime.per_vs_all',glm_dtime.location',...
    'VariableNames',{'detection_time','PAC','per_vs_all','location'});
mdl = fitglm(tbl,'detection_time~PAC+per_vs_all+location'...
    ,'CategoricalVars',[2,3])

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for for pooled mice vs per mouse lick decision times (naive) \n'])
try
    [output_data] = drgMutiRanksumorTtest(p_dtime_stats);
    fprintf(1, '\n\n')
catch
end

figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

bar_ii=1;
glm_ii=0;
ii_stats=0;
glm_dtime=[];
p_dtime_stats=[];

for PACii=[1 3]
    %Peak
    bar_ii=bar_ii+1;
    for ii_odor_pair=1:no_odor_pairs
        peak_t_detect_all(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).t_detect_peak_all;
    end
    for ii_odor_pair=1:no_odor_pairs
        peak_t_detect(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_peak_t_detect;
    end
    
    mean_peak_t_detect=mean(peak_t_detect)';
    handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_peak_t_detect=mean_peak_t_detect;
    CIpeak = bootci(1000, {@mean, peak_t_detect})';
    p1=bar(bar_ii,mean_peak_t_detect,'EdgeColor','k','FaceColor',[1 0.7 0.7]);
    plot(bar_ii,mean_peak_t_detect,'ok','MarkerFaceColor','k','MarkerSize',5)
    plot(bar_ii*ones(1,length(peak_t_detect)),peak_t_detect,'ok')
    plot([bar_ii bar_ii],CIpeak,'-k','LineWidth',3)
    
    glm_dtime.data(glm_ii+1:glm_ii+length(peak_t_detect))=peak_t_detect;
    glm_dtime.PACii(glm_ii+1:glm_ii+length(peak_t_detect))=PACii;
    glm_dtime.per_vs_all(glm_ii+1:glm_ii+length(t_detect_licks_all))=1;
    glm_dtime.location(glm_ii+1:glm_ii+length(lick_t_detect))=location;
    glm_ii=glm_ii+length(peak_t_detect);
   
    ii_stats=ii_stats+1;
    p_dtime_stats(ii_stats).data=peak_t_detect;
    p_dtime_stats(ii_stats).description=[PACnames{PACii} ' peak mean per mouse'];
    p_dtime_stats(ii_stats).PACii=PACii;
    p_dtime_stats(ii_stats).lick_peak=2;
    
    bar_ii=bar_ii+1;
    mean_peak_t_detect_all=mean(peak_t_detect_all)';
    handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_peak_t_detect_all=mean_peak_t_detect_all;
    CIpeak = bootci(1000, {@mean, peak_t_detect_all})';
    p3=bar(bar_ii,mean_peak_t_detect_all,'r');
    plot(bar_ii,mean_peak_t_detect_all,'ok','MarkerFaceColor','k','MarkerSize',5)
    plot(bar_ii*ones(1,length(peak_t_detect_all)),peak_t_detect_all,'ok')
    plot([bar_ii bar_ii],CIpeak,'-k','LineWidth',3)
    
    glm_dtime.data(glm_ii+1:glm_ii+length(peak_t_detect_all))=peak_t_detect_all;
    glm_dtime.PACii(glm_ii+1:glm_ii+length(peak_t_detect_all))=PACii;
    glm_dtime.per_vs_all(glm_ii+1:glm_ii+length(t_detect_licks_all))=2;
    glm_dtime.location(glm_ii+1:glm_ii+length(lick_t_detect))=location;
    glm_ii=glm_ii+length(peak_t_detect_all);
    
    ii_stats=ii_stats+1;
    p_dtime_stats(ii_stats).data=peak_t_detect_all;
    p_dtime_stats(ii_stats).description=[PACnames{PACii} ' peak pooled mice'];
    p_dtime_stats(ii_stats).PACii=PACii;
    p_dtime_stats(ii_stats).lick_peak=2;
    
    for ii_odor_pair=1:no_odor_pairs
        plot([bar_ii-1 bar_ii], [peak_t_detect(ii_odor_pair) peak_t_detect_all(ii_odor_pair)],'-k')
    end
    
    bar_ii=bar_ii+2;
    
end

title('Peak decision making time for naive mice')
ylabel('Decision time (sec)')
legend([p1 p3],{'Mean per mouse','Pooled mice'})
xticks([1.5 4.5])
xticklabels({'Beta','High gamma'})
ylim([0 5])

%Perform the glm
fprintf(1, ['\n\nglm for peak decision time for for pooled mice vs per mouse  (naive) \n'])
tbl = table(glm_dtime.data',glm_dtime.PACii',glm_dtime.per_vs_all',glm_dtime.location',...
    'VariableNames',{'detection_time','PAC','per_vs_all','location'});
mdl = fitglm(tbl,'detection_time~PAC+per_vs_all+location'...
    ,'CategoricalVars',[2,3])

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for peak decision times for pooled mice vs per mouse (naive) \n'])
try
    [output_data] = drgMutiRanksumorTtest(p_dtime_stats);
    fprintf(1, '\n\n')
catch
end


%Plot decision time reationship using all mouse decision lick times per mouse
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

bar_ii=1;
glma_ii=0;
ii_stats=0;
glma_dtime=[];
p_dtime_statsa=[];

glma_pn_ii=0;
ii_stats_pn=0;
glma_pn_dtime=[];
p_dtime_statsa_pn=[];

for PACii=[1 3]
    
    percent_correct_ii=1
        
        groupNo=1;
        
              %Plot decision time reationship
          
            
            %Licks
            for ii_odor_pair=1:no_odor_pairs
                lick_t_detect(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_lick_t_detect;
            end
            for ii_odor_pair=1:no_odor_pairs
                t_detect_licks_all(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).t_detect_licks_all;
            end
            mean_t_detect_licks_all=mean(t_detect_licks_all)';
            CIlicks = bootci(1000, {@mean, t_detect_licks_all})';
            p1=bar(bar_ii,mean_t_detect_licks_all,'g','LineWidth',2);
            plot(bar_ii,mean_t_detect_licks_all,'ok','MarkerFaceColor','k','MarkerSize',8)
            plot(bar_ii*ones(1,length(t_detect_licks_all)),t_detect_licks_all,'ok','MarkerSize',8)
            plot([bar_ii bar_ii],CIlicks,'-k','LineWidth',2)
            
            glma_dtime.data(glma_ii+1:glma_ii+length(t_detect_licks_all))=t_detect_licks_all;
            glma_dtime.PACii(glma_ii+1:glma_ii+length(t_detect_licks_all))=PACii;
            glma_dtime.lick_peak(glma_ii+1:glma_ii+length(t_detect_licks_all))=1;
            glma_dtime.proficient(glma_ii+1:glma_ii+length(t_detect_licks_all))=1;
            glma_ii=glma_ii+length(t_detect_licks_all);
            
            ii_stats=ii_stats+1;
            p_dtime_statsa(ii_stats).data=t_detect_licks_all;
            p_dtime_statsa(ii_stats).description=[PACnames{PACii} ' lick'];
            p_dtime_statsa(ii_stats).PACii=PACii;
            p_dtime_statsa(ii_stats).lick_peak=1;
            
            glma_pn_dtime.data(glma_pn_ii+1:glma_pn_ii+length(t_detect_licks_all))=t_detect_licks_all;
            glma_pn_dtime.PACii(glma_pn_ii+1:glma_pn_ii+length(t_detect_licks_all))=PACii;
            glma_pn_dtime.lick_peak(glma_pn_ii+1:glma_pn_ii+length(t_detect_licks_all))=1;
            glma_pn_dtime.proficient(glma_pn_ii+1:glma_pn_ii+length(t_detect_licks_all))=1;
            glma_pn_ii=glma_pn_ii+length(t_detect_licks_all);
            
            ii_stats_pn=ii_stats_pn+1;
            p_dtime_statsa_pn(ii_stats_pn).data=t_detect_licks_all;
            p_dtime_statsa_pn(ii_stats_pn).description=[PACnames{PACii} ' lick proficient'];
            p_dtime_statsa_pn(ii_stats_pn).PACii=PACii;
            p_dtime_statsa_pn(ii_stats_pn).lick_peak=1;
            
            %Peak
            bar_ii=bar_ii+1;
            for ii_odor_pair=1:no_odor_pairs
                peak_t_detect(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).t_detect_peak_all;
            end
            mean_peak_t_detect=mean(peak_t_detect)';
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_peak_t_detect=mean_peak_t_detect;
            CIpeak = bootci(1000, {@mean, peak_t_detect})';
            p3=bar(bar_ii,mean_peak_t_detect,'r','LineWidth',2);
            plot(bar_ii,mean_peak_t_detect,'ok','MarkerFaceColor','k','MarkerSize',8)
            plot(bar_ii*ones(1,length(peak_t_detect)),peak_t_detect,'ok','MarkerSize',8)
            plot([bar_ii bar_ii],CIpeak,'-k','LineWidth',2)
            
            glma_dtime.data(glma_ii+1:glma_ii+length(peak_t_detect))=peak_t_detect;
            glma_dtime.PACii(glma_ii+1:glma_ii+length(peak_t_detect))=PACii;
            glma_dtime.lick_peak(glma_ii+1:glma_ii+length(peak_t_detect))=2;
            glma_dtime.proficient(glma_ii+1:glma_ii+length(peak_t_detect))=1;
            glma_ii=glma_ii+length(peak_t_detect);
            
            ii_stats=ii_stats+1;
            p_dtime_statsa(ii_stats).data=peak_t_detect;
            p_dtime_statsa(ii_stats).description=[PACnames{PACii} ' peak'];
            p_dtime_statsa(ii_stats).PACii=PACii;
            p_dtime_statsa(ii_stats).lick_peak=2;
            
            glma_pn_dtime.data(glma_pn_ii+1:glma_pn_ii+length(peak_t_detect))=peak_t_detect;
            glma_pn_dtime.PACii(glma_pn_ii+1:glma_pn_ii+length(peak_t_detect))=PACii;
            glma_pn_dtime.lick_peak(glma_pn_ii+1:glma_pn_ii+length(peak_t_detect))=2;
            glma_pn_dtime.proficient(glma_pn_ii+1:glma_pn_ii+length(peak_t_detect))=1;
            glma_pn_ii=glma_pn_ii+length(peak_t_detect);
            
            ii_stats_pn=ii_stats_pn+1;
            p_dtime_statsa_pn(ii_stats_pn).data=peak_t_detect;
            p_dtime_statsa_pn(ii_stats_pn).description=[PACnames{PACii} ' peak proficient'];
            p_dtime_statsa_pn(ii_stats_pn).PACii=PACii;
            p_dtime_statsa_pn(ii_stats_pn).lick_peak=2;
            
            for ii_odor_pair=1:no_odor_pairs
                plot([bar_ii-1 bar_ii], [t_detect_licks_all(ii_odor_pair) peak_t_detect(ii_odor_pair)],'-k','LineWidth',2)
            end
             
            bar_ii=bar_ii+2;

end

title('Decision making time for proficient mice: pooled mouse LDA, Figure 5D')
ylabel('Decision time (sec)')
legend([p1 p3],{'Lick','Peak'})
xticks([1.5 4.5])
xticklabels({'Beta','High gamma'})
ylim([0 5])

%Perform the glma
fprintf(1, ['\n\nglm for LDA and lick decision time for pooled mice (proficient) Theta/' PACnames{PACii} ' Fig 5D\n'])
tbl = table(glma_dtime.data',glma_dtime.PACii',glma_dtime.lick_peak',...
    'VariableNames',{'detection_time','PAC','lick_peak'});
mdl = fitglm(tbl,'detection_time~PAC+lick_peak'...
    ,'CategoricalVars',[2,3])

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for proficient mice area all mouse decision time for Theta/' PACnames{PACii} '\n'])
try
    [output_data] = drgMutiRanksumorTtest(p_dtime_statsa);
    fprintf(1, '\n\n')
catch
end

%Plot decision times for the two locations
%Now do Naive
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

%First do licks
bar_ii=0;

%Peak
ii_loc1=0;
ii_loc2=0;
for ii_odor_pair=1:no_odor_pairs
    if location(ii_odor_pair)==1
        ii_loc1=ii_loc1+1;
        peak_t_detect_loc1(ii_loc1)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).t_detect_peak_all;
    else
        ii_loc2=ii_loc2+1;
        peak_t_detect_loc2(ii_loc2)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).t_detect_peak_all;
    end
end

bar_ii=bar_ii+1;
mean_peak_t_detect_loc1=mean(peak_t_detect_loc1)';
CIpeak_loc1 = bootci(1000, {@mean, peak_t_detect_loc1})';
p1=bar(bar_ii,mean_peak_t_detect_loc1,'r');
plot(bar_ii,mean_peak_t_detect_loc1,'ok','MarkerFaceColor','k','MarkerSize',5)
plot(bar_ii*ones(1,length(peak_t_detect_loc1)),peak_t_detect_loc1,'ok')
plot([bar_ii bar_ii],CIpeak_loc1,'-k','LineWidth',3)

bar_ii=bar_ii+1;
mean_peak_t_detect_loc2=mean(peak_t_detect_loc2)';
CIpeak_loc2 = bootci(1000, {@mean, peak_t_detect_loc2})';
bar(bar_ii,mean_peak_t_detect_loc2,'r');
plot(bar_ii,mean_peak_t_detect_loc2,'ok','MarkerFaceColor','k','MarkerSize',5)
plot(bar_ii*ones(1,length(peak_t_detect_loc2)),peak_t_detect_loc2,'ok')
plot([bar_ii bar_ii],CIpeak_loc2,'-k','LineWidth',3)

bar_ii=bar_ii+2;

%Now do licks
ii_loc1=0;
ii_loc2=0;
for ii_odor_pair=1:no_odor_pairs
    if location(ii_odor_pair)==1
        ii_loc1=ii_loc1+1;
    t_detect_licks_all_loc1(ii_loc1)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).t_detect_licks_all;
    else
       ii_loc2=ii_loc2+1;
    t_detect_licks_all_loc2(ii_loc2)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).t_detect_licks_all;  
    end
end

bar_ii=bar_ii+1;
mean_t_detect_licks_all_loc1=mean(t_detect_licks_all_loc1)';
CIlicks_loc1 = bootci(1000, {@mean, t_detect_licks_all_loc1})';
p3=bar(bar_ii,mean_t_detect_licks_all_loc1,'EdgeColor','k','FaceColor',[0.7 0.7 0.7]);
plot(bar_ii,mean_t_detect_licks_all_loc1,'ok','MarkerFaceColor','k','MarkerSize',5)
plot(bar_ii*ones(1,length(t_detect_licks_all_loc1)),t_detect_licks_all_loc1,'ok')
plot([bar_ii bar_ii],CIlicks_loc1,'-k','LineWidth',3)

bar_ii=bar_ii+1;
mean_t_detect_licks_all_loc2=mean(t_detect_licks_all_loc2)';
CIlicks_loc2 = bootci(1000, {@mean, t_detect_licks_all_loc2})';
bar(bar_ii,mean_t_detect_licks_all_loc2,'EdgeColor','k','FaceColor',[0.7 0.7 0.7]);
plot(bar_ii,mean_t_detect_licks_all_loc2,'ok','MarkerFaceColor','k','MarkerSize',5)
plot(bar_ii*ones(1,length(t_detect_licks_all_loc2)),t_detect_licks_all_loc2,'ok')
plot([bar_ii bar_ii],CIlicks_loc2,'-k','LineWidth',3)

title('Decision making time for proficient mice: pooled mouse, separated by location')
ylabel('Decision time (sec)')
legend([p1 p3],{'Peak','Lick'})
xticks([1 2 5 6])
xticklabels({'loc1','loc2','loc1','loc2'})
ylim([0 5])

%Now do Naive
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

bar_ii=1;
glma_ii=0;
ii_stats=0;
glma_dtime=[];
p_dtime_statsa=[];

for PACii=[1 3]
    
    percent_correct_ii=2;
        
        groupNo=1;
        
              %Plot decision time reationship
          
            
            %Licks
            for ii_odor_pair=1:no_odor_pairs
                lick_t_detect(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_lick_t_detect;
            end
            for ii_odor_pair=1:no_odor_pairs
                t_detect_licks_all(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).t_detect_licks_all;
            end
            mean_t_detect_licks_all=mean(t_detect_licks_all)';
            CIlicks = bootci(1000, {@mean, t_detect_licks_all})';
            p1=bar(bar_ii,mean_t_detect_licks_all,'g','LineWidth',2);
            plot(bar_ii,mean_t_detect_licks_all,'ok','MarkerFaceColor','k','MarkerSize',8)
            plot(bar_ii*ones(1,length(t_detect_licks_all)),t_detect_licks_all,'ok','MarkerSize',8)
            plot([bar_ii bar_ii],CIlicks,'-k','LineWidth',2)
            
            glma_dtime.data(glma_ii+1:glma_ii+length(t_detect_licks_all))=t_detect_licks_all;
            glma_dtime.PACii(glma_ii+1:glma_ii+length(t_detect_licks_all))=PACii;
            glma_dtime.lick_peak(glma_ii+1:glma_ii+length(t_detect_licks_all))=1;
            glma_dtime.location(glma_ii+1:glma_ii+length(t_detect_licks_all))=location;
            glma_ii=glma_ii+length(t_detect_licks_all);
            
            ii_stats=ii_stats+1;
            p_dtime_statsa(ii_stats).data=t_detect_licks_all;
            p_dtime_statsa(ii_stats).description=[PACnames{PACii} ' lick'];
            p_dtime_statsa(ii_stats).PACii=PACii;
            p_dtime_statsa(ii_stats).lick_peak=1;
            
            glma_pn_dtime.data(glma_pn_ii+1:glma_pn_ii+length(t_detect_licks_all))=t_detect_licks_all;
            glma_pn_dtime.PACii(glma_pn_ii+1:glma_pn_ii+length(t_detect_licks_all))=PACii;
            glma_pn_dtime.lick_peak(glma_pn_ii+1:glma_pn_ii+length(t_detect_licks_all))=1;
            glma_pn_dtime.proficient(glma_pn_ii+1:glma_pn_ii+length(t_detect_licks_all))=0;
            glma_pn_ii=glma_pn_ii+length(t_detect_licks_all);
            
            ii_stats_pn=ii_stats_pn+1;
            p_dtime_statsa_pn(ii_stats_pn).data=t_detect_licks_all;
            p_dtime_statsa_pn(ii_stats_pn).description=[PACnames{PACii} ' lick'];
            p_dtime_statsa_pn(ii_stats_pn).PACii=PACii;
            p_dtime_statsa_pn(ii_stats_pn).lick_peak=1;
            
            %Peak
            bar_ii=bar_ii+1;
            for ii_odor_pair=1:no_odor_pairs
                peak_t_detect(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).t_detect_peak_all;
            end
            mean_peak_t_detect=mean(peak_t_detect)';
            handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_peak_t_detect=mean_peak_t_detect;
            CIpeak = bootci(1000, {@mean, peak_t_detect})';
            p3=bar(bar_ii,mean_peak_t_detect,'r','LineWidth',2);
            plot(bar_ii,mean_peak_t_detect,'ok','MarkerFaceColor','k','MarkerSize',8)
            plot(bar_ii*ones(1,length(peak_t_detect)),peak_t_detect,'ok','MarkerSize',8)
            plot([bar_ii bar_ii],CIpeak,'-k','LineWidth',2)
            
            glma_dtime.data(glma_ii+1:glma_ii+length(peak_t_detect))=peak_t_detect;
            glma_dtime.PACii(glma_ii+1:glma_ii+length(peak_t_detect))=PACii;
            glma_dtime.lick_peak(glma_ii+1:glma_ii+length(peak_t_detect))=2;
            glma_dtime.location(glma_ii+1:glma_ii+length(peak_t_detect))=location;
            glma_ii=glma_ii+length(peak_t_detect);
            
            ii_stats=ii_stats+1;
            p_dtime_statsa(ii_stats).data=peak_t_detect;
            p_dtime_statsa(ii_stats).description=[PACnames{PACii} ' peak'];
            p_dtime_statsa(ii_stats).PACii=PACii;
            p_dtime_statsa(ii_stats).lick_peak=2;
            
            glma_pn_dtime.data(glma_pn_ii+1:glma_pn_ii+length(peak_t_detect))=peak_t_detect;
            glma_pn_dtime.PACii(glma_pn_ii+1:glma_pn_ii+length(peak_t_detect))=PACii;
            glma_pn_dtime.lick_peak(glma_pn_ii+1:glma_pn_ii+length(peak_t_detect))=2;
            glma_pn_dtime.proficient(glma_pn_ii+1:glma_pn_ii+length(peak_t_detect))=0;
            glma_pn_ii=glma_pn_ii+length(peak_t_detect);
            
            ii_stats_pn=ii_stats_pn+1;
            p_dtime_statsa_pn(ii_stats_pn).data=peak_t_detect;
            p_dtime_statsa_pn(ii_stats_pn).description=[PACnames{PACii} ' peak naive'];
            p_dtime_statsa_pn(ii_stats_pn).PACii=PACii;
            p_dtime_statsa_pn(ii_stats_pn).lick_peak=2;
            
            for ii_odor_pair=1:no_odor_pairs
                plot([bar_ii-1 bar_ii], [t_detect_licks_all(ii_odor_pair) peak_t_detect(ii_odor_pair)],'-k','LineWidth',2)
            end
            
            bar_ii=bar_ii+2;

end

title('Decision making time for naive mice: pooled mouse LDA, Figure 5D')
ylabel('Decision time (sec)')
legend([p1 p3],{'Lick','Peak'})
xticks([1.5 4.5])
xticklabels({'Beta','High gamma'})
ylim([0 5])



%Perform the glma
fprintf(1, ['\n\nglm for LDA and lick decision time for pooled mice for naive mice, Figure 5D Theta/' PACnames{PACii} '\n'])
tbl = table(glma_dtime.data',glma_dtime.PACii',glma_dtime.lick_peak',...
    'VariableNames',{'detection_time','PAC','lick_peak'});
mdl = fitglm(tbl,'detection_time~PAC+lick_peak'...
    ,'CategoricalVars',[2,3])

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for area under the curve for naive mice Theta/' PACnames{PACii} ' Figure 5D\n'])
try
    [output_data] = drgMutiRanksumorTtest(p_dtime_statsa);
    fprintf(1, '\n\n')
catch
end

%Perform the glma_pn
fprintf(1, ['\n\nglm for LDA and lick decision time for pooled mice for proficient and naive mice, Figure 5D Theta/' PACnames{PACii} '\n'])
tbl = table(glma_pn_dtime.data',glma_pn_dtime.PACii',glma_pn_dtime.lick_peak',glma_pn_dtime.proficient',...
    'VariableNames',{'detection_time','PAC','lick_peak','proficient'});
mdl = fitglm(tbl,'detection_time~PAC+lick_peak+proficient'...
    ,'CategoricalVars',[2,3,4])

%Now do the analysis for AUC
%Proficient
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

bar_ii=1;
glm_ii=0;
ii_stats=0;
glm_dtime=[];
p_dtime_stats=[];

percent_correct_ii=1;

for PACii=[1 3]
    %Peak
    for ii_odor_pair=1:no_odor_pairs
        peak_AUC_all(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auc_peak_all;
    end
    for ii_odor_pair=1:no_odor_pairs
        peak_AUC(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_peak_auc;
    end
    
    mean_peak_AUC=mean(peak_AUC)';
    handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_peak_AUC=mean_peak_AUC;
    CIpeak = bootci(1000, {@mean, peak_AUC})';
    p1=bar(bar_ii,mean_peak_AUC,'EdgeColor','k','FaceColor',[1 0.7 0.7]);
    plot(bar_ii,mean_peak_AUC,'ok','MarkerFaceColor','k','MarkerSize',5)
    plot(bar_ii*ones(1,length(peak_AUC)),peak_AUC,'ok')
    plot([bar_ii bar_ii],CIpeak,'-k','LineWidth',3)
    
    glm_dtime.data(glm_ii+1:glm_ii+length(peak_AUC))=peak_AUC;
    glm_dtime.PACii(glm_ii+1:glm_ii+length(peak_AUC))=PACii;
    glm_dtime.per_vs_all(glm_ii+1:glm_ii+length(peak_AUC))=1;
    glm_dtime.location(glm_ii+1:glm_ii+length(peak_AUC))=location;
    glm_ii=glm_ii+length(peak_AUC);
   
    ii_stats=ii_stats+1;
    p_dtime_stats(ii_stats).data=peak_AUC;
    p_dtime_stats(ii_stats).description=[PACnames{PACii} ' peak AUC mean per mouse'];
    p_dtime_stats(ii_stats).PACii=PACii;
    p_dtime_stats(ii_stats).lick_peak=2;
    
    bar_ii=bar_ii+1;
    mean_peak_AUC_all=mean(peak_AUC_all)';
    handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_peak_AUC_all=mean_peak_AUC_all;
    CIpeak = bootci(1000, {@mean, peak_AUC_all})';
    p3=bar(bar_ii,mean_peak_AUC_all,'r');
    plot(bar_ii,mean_peak_AUC_all,'ok','MarkerFaceColor','k','MarkerSize',5)
    plot(bar_ii*ones(1,length(peak_AUC_all)),peak_AUC_all,'ok')
    plot([bar_ii bar_ii],CIpeak,'-k','LineWidth',3)
    
    glm_dtime.data(glm_ii+1:glm_ii+length(peak_AUC_all))=peak_AUC_all;
    glm_dtime.PACii(glm_ii+1:glm_ii+length(peak_AUC_all))=PACii;
    glm_dtime.per_vs_all(glm_ii+1:glm_ii+length(peak_AUC_all))=2;
    glm_dtime.location(glm_ii+1:glm_ii+length(peak_AUC_all))=location;
    glm_ii=glm_ii+length(peak_AUC_all);
    
    ii_stats=ii_stats+1;
    p_dtime_stats(ii_stats).data=peak_AUC_all;
    p_dtime_stats(ii_stats).description=[PACnames{PACii} ' peak AUC pooled mice'];
    p_dtime_stats(ii_stats).PACii=PACii;
    p_dtime_stats(ii_stats).lick_peak=2;
    
    for ii_odor_pair=1:no_odor_pairs
        plot([bar_ii-1 bar_ii], [peak_AUC(ii_odor_pair) peak_AUC_all(ii_odor_pair)],'-k')
    end
    
    bar_ii=bar_ii+2;
    
end

title('Peak AUC for proficient mice')
ylabel('AUC')
legend([p1 p3],{'Mean per mouse','Pooled mice'})
xticks([1.5 4.5])
xticklabels({'Beta','High gamma'})
ylim([-0.1 1])

%Perform the glm
fprintf(1, ['\n\nglm for peak AUC for pooled mice (proficient) \n'])
tbl = table(glm_dtime.data',glm_dtime.PACii',glm_dtime.per_vs_all',glm_dtime.location',...
    'VariableNames',{'detection_time','PAC','per_vs_all','location'});
mdl = fitglm(tbl,'detection_time~PAC+per_vs_all+location'...
    ,'CategoricalVars',[2,3])

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for peak AUC (proficient) \n'])
try
    [output_data] = drgMutiRanksumorTtest(p_dtime_stats);
    fprintf(1, '\n\n')
catch
end


%Now do Naive
%Plot decision time reationship using mean decision lick times per mouse
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

bar_ii=1;
glm_ii=0;
ii_stats=0;
glm_dtime=[];
p_dtime_stats=[];

%mean_peak_AUC is the lick decision time computed per mouse
%AUC_peak_all is the lick decision time computed with all mice pooled

percent_correct_ii=2;

groupNo=1;

for PACii=[1 3]
    
    
    
    %Plot decision time reationship between per mouse and all mouse decision
    %times
    
    %Licks
    for ii_odor_pair=1:no_odor_pairs
        peak_AUC(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auc_peak_all;
    end
    for ii_odor_pair=1:no_odor_pairs
        AUC_peak_all(ii_odor_pair)=handles_lda_per_odor_pair(ii_odor_pair).handles_outp.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).mean_peak_auc;
    end
    
    mean_peak_AUC=mean(peak_AUC)';
    CIlicks = bootci(1000, {@mean, peak_AUC})';
    p1=bar(bar_ii,mean_peak_AUC,'EdgeColor','k','FaceColor',[1 0.7 0.7]);
    plot(bar_ii,mean_peak_AUC,'ok','MarkerFaceColor','k','MarkerSize',5)
    plot(bar_ii*ones(1,length(peak_AUC)),peak_AUC,'ok')
    plot([bar_ii bar_ii],CIlicks,'-k','LineWidth',3)
    
    glm_dtime.data(glm_ii+1:glm_ii+length(peak_AUC))=peak_AUC;
    glm_dtime.PACii(glm_ii+1:glm_ii+length(peak_AUC))=PACii;
    glm_dtime.per_vs_all(glm_ii+1:glm_ii+length(peak_AUC))=1;
    glm_dtime.location(glm_ii+1:glm_ii+length(peak_AUC))=location;
    glm_ii=glm_ii+length(peak_AUC);
    
    ii_stats=ii_stats+1;
    p_dtime_stats(ii_stats).data=peak_AUC;
    p_dtime_stats(ii_stats).description=[PACnames{PACii} ' AUC mean per mouse'];
    p_dtime_stats(ii_stats).PACii=PACii;
    p_dtime_stats(ii_stats).lick_peak=1;
    
    bar_ii=bar_ii+1;
    mean_AUC_peak_all=mean(AUC_peak_all)';
    CIlicks_all = bootci(1000, {@mean, AUC_peak_all})';
    p2=bar(bar_ii,mean_AUC_peak_all,'EdgeColor','r','FaceColor','r');
    plot(bar_ii,mean_AUC_peak_all,'ok','MarkerFaceColor','k','MarkerSize',5)
    plot(bar_ii*ones(1,length(AUC_peak_all)),AUC_peak_all,'ok')
    plot([bar_ii bar_ii],CIlicks_all,'-k','LineWidth',3)
    
    glm_dtime.data(glm_ii+1:glm_ii+length(AUC_peak_all))=AUC_peak_all;
    glm_dtime.PACii(glm_ii+1:glm_ii+length(AUC_peak_all))=PACii;
    glm_dtime.per_vs_all(glm_ii+1:glm_ii+length(AUC_peak_all))=2;
    glm_dtime.location(glm_ii+1:glm_ii+length(AUC_peak_all))=location;
    glm_ii=glm_ii+length(AUC_peak_all);
    
    ii_stats=ii_stats+1;
    p_dtime_stats(ii_stats).data=AUC_peak_all;
    p_dtime_stats(ii_stats).description=[PACnames{PACii} ' AUC pooled mice'];
    p_dtime_stats(ii_stats).PACii=PACii;
    p_dtime_stats(ii_stats).lick_peak=1;
    
    
    for ii_odor_pair=1:no_odor_pairs
        plot([bar_ii-1 bar_ii], [peak_AUC(ii_odor_pair) AUC_peak_all(ii_odor_pair)],'-k')
    end
    
    bar_ii=bar_ii+2;
end

title('Peak AUC for naive mice')
ylabel('AUC')
legend([p1 p2],{'Mean per mouse','Pooled mice'})
xticks([1.5 4.5])
xticklabels({'Beta','High gamma'})
ylim([-0.1 1])

%Perform the glm
fprintf(1, ['\n\nglm for peak AUC for pooled mice (naive) \n'])
tbl = table(glm_dtime.data',glm_dtime.PACii',glm_dtime.per_vs_all',glm_dtime.location',...
    'VariableNames',{'detection_time','PAC','per_vs_all','location'});
mdl = fitglm(tbl,'detection_time~PAC+per_vs_all+location'...
    ,'CategoricalVars',[2,3])

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for peak AUC (naive) \n'])
try
    [output_data] = drgMutiRanksumorTtest(p_dtime_stats);
    fprintf(1, '\n\n')
catch
end
pffft=1


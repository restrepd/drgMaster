%drgAnalysisBatchLFPCaMKIIv2_multi_file
%Runs drgAnalysisBatchLFPCaMKIIv2 for multiple files

pathName_pars_choice='/Users/restrepd/Documents/Projects/CaMKII_analysis/Lick-theta analysis/';
fileName_pars_choice='drgLFPBatchAnalPars_CaMKII_spm_dr_case2PRPtime80.m';

pathName_data='/Users/restrepd/Documents/Projects/CaMKII_analysis/Lick-theta analysis/';

fileName{1}='spm_LFP_APEB_one_window_12102021.mat';
fileName{2}='spm_LFP_ethyl_one_win12182021.mat';
fileName{3}='spm_LFP_EBAP_one_window_12192021.mat';
fileName{4}='spm_LFP_PAEA_one_window122121.mat';
fileName{5}='spm_LFP_pz1ethyl_one_window_122321.mat';
fileName{6}='spm_LFP_pz1PAEA_one_win12302021.mat';
fileName{7}='spm_LFP_pzz1EAPA_one_window01012022.mat';
fileName{8}='spm_LFP_pzz1PAEA_one_window12302021.mat';

first_file=1;
show_figures=0;

%Parallel batch processing for each file
all_files_present=1;

if exist([pathName_pars_choice fileName_pars_choice])==0
    fprintf(1, ['Program will be terminated because choice file ' fileName_pars_choice ' does not exist\n'],filNum);
    all_files_present=0;
end

for filNum=first_file:length(fileName)
    %Make sure that all the files exist
    if exist([pathName_data fileName{filNum}])==0
        fprintf(1, ['Program will be terminated because file No %d, ' fileName{filNum} ' does not exist\n'],filNum);
        all_files_present=0;
    end
end

if  all_files_present==1
    for filNum=first_file:length(fileName)
        fprintf(1, ['Processing file No %d, ' fileName{filNum} '\n'],filNum);
        drgAnalysisBatchLFPCaMKIIv2(fileName_pars_choice,pathName_pars_choice,fileName{filNum},pathName_data,show_figures)
    end
end
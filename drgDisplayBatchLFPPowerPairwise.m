function drgDisplayBatchLFPPowerPairwise(handles)

%This function displays the LFP power spectrum for drgRunBatch-generated
%data

%Which analysis is performed is determined by the value enterd in the
%variable which_display:
%
%
% 1 ERP analysis compare auROC in the last few trials of pre with first few trials of post
%
% 2 Generate delta LFP power and auROC for reversal figure for Daniel's paper
%
% 3 For a subset of first files for events 1 and 2 plot the LFP bandwide spectrum,
%   LFP power histograms for different bandwidths for each electrode and LFP auROCs.
%   To generate Fig. 2 for Daniels' LFP power paper enter the proficient files
%
% 4 Generate LFP power auROC for Fig. 3 for Daniel's paper. first vs last.
%
%
% 5 Compare auROC in the last few trials of pre with first few trials of post
%    Used for old Fig. 4 of Daniel's paper with acetoethylben_electrode9202017.mat
%
% 6 For a subset of first files for events 1 and 2 plot the ERP LFP bandwide spectrum,
%   ERP LFP power histograms for different bandwidths for each electrode and ERP LFP auROCs.
%
% 7 Generate ERP LFP power auROC. first vs last.
%
% 8 Compare auROC for ERP LFP in the last few trials of pre with first few trials of post
%   Used for New Fig. 7 of Daniel's paper
%
% 9 Compare auROC for power LFP in two percent windows for all of the files 
%
% 10 Compare auROC for power LFP for two groups (e.g. NRG1 vs control)
% within one precent window

%% First enter the choices for what will be analyzed.
% THESE VALUES ARE IMPORTANT and differ for each user

% 
% 
% % For Daniel's acetoethylben_nlvsl_02042018.mat
% % drgbChoicesDanielAPEBEnlvsl02042018
% %DO NOT USE 159867
% % For Fig. 6 of Daniel's paper run with which_display=5
% % NOTE: I am not using file 6 because the animal was not proficient
% winNo=2;
% refWin=1;
% which_display=5;
% 
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% 
% %Experiment pairs
% %Important: The first file must be the experiment performed first
% %For example in acetophenone ethyl benzoate no laser is first, laser is
% %second
% file_pairs=[
%     7 1;
%     8 2;
%     9 3;
%     10 4;
%     11 5;
%     12 6;
%     17 13;
%     18 14;
%     19 15;
%     20 16];
% no_file_pairs=10;
% 
% 
% trials_to_process=20;
% min_trials_per_event=4;
% 
% grpre=[1 3];
% grpost=[2 4];

% % For Daniel's isoamylmintstart_electrodeChR281217.mat
% % drgbChoicesDanielIAMOtstartElectrodes
% % For old Fig. 4 of Daniel's paper run with which_display=5
% 
% winNo=2;
% refWin=1;
% which_display=5;
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% 
% 
% %Experiment pairs
% %Important: The first file must be the experiment performed first
% %For example in acetophenone ethyl benzoate no laser is first, laser is
% %second
% file_pairs=[
%     2 11;
%     3 12;
%     1 13;
%     4 14;
%     7 17;
%     8 18;
%     22 19;
%     9 20;
%     10 21];
% no_file_pairs=9;
% 
% trials_to_process=20;
% min_trials_per_event=4;
% 
% grpre=[1 3];
% grpost=[2 4];



% % For Daniel's heptaoctanol11317.mat
% % For old Fig. 5 of Daniel's paper run with which_display=5
% winNo=2;
% refWin=1;
% which_display=5;
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
%
%
% %Experiment pairs
% %Important: The first file must be the experiment performed first
% %For example in acetophenone ethyl benzoate no laser is first, laser is
% %second
% file_pairs=[
%     9 1;
%     10 2;
%     11 3;
%     12 4;
%     13 5;
%     14 6;
%     15 7
%     16 8
%     ];
% no_file_pairs=8;
%
% comp_window=10; %Note: The ancova is not significant for all bandwidths when the comp_window is increaesed to 15
%
% grpre=[1 3];
% grpost=[2 4];

% 
% % % For Daniel's acetoethylben_firstandlast91117
% % drgbChoicesDanielAPEBEfirstandlast9617
% %New Fig. 2 G2 run with which_display=3;
% 
% 
% winNo=2;
% refWin=1;
% % which_display=3;
% which_display=13;
% 
% %
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % comp_window=10;
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% 
% % Enter the files to be processed
% files=[1 2 3 4 5 13 14 15];
% 
% trials_to_process=30;
% min_trials_per_event=4;
% output_suffix='_test.mat';

% 
% % % For Daniel's acetoethylben_firstandlast91117
% % drgbChoicesDanielAPEBEfirstandlast9617
% % Used to test comparing low to high percent correct
% 
% winNo=2;
% refWin=1;
% which_display=9;
% 
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % trials_to_process=20;
% % output_suffix='_spmouthit.mat';
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% output_suffix='pertest_spmout.mat';
% 
% 
% %Files to analyze
% files=[1 11 2 7 3 8 4 9 5 10 13 16 14 17 15 18];
% 
% percent_windows=[80 100;45 65];
% 
% file_label{1}='proficient';
% file_label{2}='naive';
% 
% min_trials_per_event=4;


<<<<<<< HEAD
% % % For Ming's Odorandlaser_firstandlast_Allfiles_fourgroups
% % drgbChoicesMingFirstLast022318_Allfiles_fourgroups
% % Used to test comparing low to high percent correct
% 
% winNo=2;
% refWin=1;
% which_display=9;
% 
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % trials_to_process=20;
% % output_suffix='_spmouthit.mat';
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% output_suffix='pernaivefirstlast_spmout.mat';
% 
% 
% %Files to analyze
% % files=[53:71];  %Naive odor
% % files=[39:52];  %ChR2
% % files=[1:16];  %Experienced odor
% files=[17:38];  %ChETA
% 
% percent_windows=[80 100;45 65];
% 
% file_label{1}='proficient';
% file_label{2}='naive';
% 
% min_trials_per_event=4;


=======
% % For Ming's Odorandlaser_firstandlast_Allfourgroups
% drgbChoicesMingFirstLast022218_Allfourgroups
% Used to test comparing low to high percent correct

winNo=2;
refWin=1;
which_display=9;

% eventType=[2 5];
% evTypeLabels={'Hit','CR'};
% trials_to_process=20;
% output_suffix='_spmouthit.mat';

eventType=[3 6];
evTypeLabels={'S+','S-'};
output_suffix='pernaivefirstlast_spmout.mat';


%Files to analyze
% files=[1:8];
files=[9:16];

percent_windows=[80 100;45 65];

file_label{1}='proficient';
file_label{2}='naive';

min_trials_per_event=4;
>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855


% 
% % % For Daniel's acetoethylben_firstandlast91117
% % drgbChoicesDanielAPEBEfirstandlast9617
% % New Fig. 3 run with which_display=4;
% 
% winNo=2;
% refWin=1;
% which_display=4;
% 
% 
% 
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % trials_to_process=20;
% % output_suffix='_spmouthit.mat';
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% trials_to_process=20;
% output_suffix='_spmout.mat';
% 
% 
% %Experiment pairs
% % Here the last session file is the first file and the first session the second file in the pair
% % file_pairs=[file_for_last_session file_for_first_session]
% % because of this front_mask=[0 1];
% file_pairs=[
%     1 11;
%     2 7;
%     3 8;
%     4 9;
%     5 10;
% %     6 12; %Exclude,159867 should NOT be used because the mouse was not
% %     proficient in the "last" session
%     13 16;
%     14 17;
%     15 18];
% no_file_pairs=8;
% 
% front_mask=[0 1];  %For each file pair enter whether you will analyze the first few trials (1) or the last(0)
% 
% file_label{1}='proficient';
% file_label{2}='learning';
% 
% min_trials_per_event=4;


% % % For Daniel's ethylacetatepropylacetatefirstandlast92617.mat
% %drgbChoicesDanielEAPAfirstandlast92617
% %New Fig. 3C run with which_display=4;
% winNo=2;
% refWin=1;
% which_display=4;
%
%
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % trials_to_process=20;
% % output_suffix='_spmouthit.mat';
%
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% trials_to_process=20;
% output_suffix='_spmout.mat';
%
%
%
% %Experiment pairs
% % Here the last session file is the first file and the first session the second file in the pair
% % file_pairs=[file_for_last_session file_for_first_session]
% % because of this front_mask=[0 1];
% file_pairs=[
%     1 6;
%     2 7;
%     3 8;
%     4 9;
%     5 10];
% no_file_pairs=5;
%
% front_mask=[0 1];  %For each file pair enter whether you will analyze the first few trials (1) or the last(0)
%
% file_label{1}='proficient';
% file_label{2}='learning';
%
%
% min_trials_per_event=4;





% % % For Daniel's ethylacetatepropylacetatefirstandlast92617.mat
% %drgbChoicesDanielEAPAfirstandlast92617
% %Fig. 1 G3 run with which_display=3;
%
% winNo=2;
% refWin=1;
% which_display=3;
%
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % comp_window=10; %This does not work with 15
%
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
%
% %Experiment pairs
% %Important: For first and last the first file must be the file for the last session
% % and the second the file for the first session
% % the program will analyze the last few trials in the last session
% % and the first few trials in the first session
% % file_pairs=[file_for_last_session file_for_first_session]
% files=[1 2 3 4 5];
% no_files=5;
%
% trials_to_process=20;
% min_trials_per_event=4;





%
% % drgbChoicesDanielIMMOfirstandlast91817
% % % For Daniel's New Fig. 3B run with isomin_firstandlastIAMO91817
% % % and which_display=4
% winNo=2;
% refWin=1;
% which_display=4;
%
%
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % trials_to_process=20;
% % output_suffix='_spmouthit.mat';
%
%
% eventType=[3 6];
% evTypeLabels={'S+','S-'}; %Note: This is currently not used
% trials_to_process=20;
% output_suffix='_spmout.mat';
%
% %Experiment pairs
% % Here the last session file is the first file and the first session the second file in the pair
% % file_pairs=[file_for_last_session file_for_first_session]
% % because of this front_mask=[0 1];
%
% file_pairs=[
% %     1 5; %Do not use this pair, the learning was deficient
%     2 4;
%     3 6;
%     7 11;
%     9 12;
% %     8 14; Taken out because of a small number of trials
%     10 13];
% no_file_pairs=5;
%
%
% front_mask=[0 1];  %For each file pair enter whether you will analyze the first few trials (1) or the last(0)
%
% file_label{1}='proficient';
% file_label{2}='learning';
%
%
% min_trials_per_event=4;



%
% % IMPORTANT: IAMO was used for Daniel's New Figs. 2A to G1
% % drgbChoicesDanielIMMOfirstandlast91817
% % % For Daniel's Fig. 1 run with isomin_firstandlastIAMO91817
% % % and which_display=3
%
% winNo=2;
% refWin=1;
% which_display=3;
%
%
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % comp_window=10; %works well with 8-12
%
%
% eventType=[3 6];
% evTypeLabels={'S+','S-'}; %Note: This is currently not used
%
%
% % Enter the files to be processed
% files=[2 3 7 8 9 10];
% no_files=6;
%
% trials_to_process=30;
% min_trials_per_event=4;


% 
% % For Daniel's ethylacetatepropylacetate122117
% % drgbChoicesDanielEAPA122117
% % For old Fig. 4 of Daniel's paper run with which_display=5
% winNo=2;
% refWin=1;
% which_display=5;
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% 
% 
% %Experiment pairs
% %Important: The first file must be the experiment performed first
% %For example in acetophenone ethyl benzoate no laser is first, laser is
% %second
% file_pairs=[
%     13 1;
%     14 2;
%     15 3;
%     16 4;
%     17 5;
%     18 6;
%     19 7;
%     20 8;
%     21 9;
%     22 10;
%     23 11;
%     24 12
%     ];
% no_file_pairs=12;
% 
% trials_to_process=20;
% min_trials_per_event=4;
% 
% grpre=[1 3];
% grpost=[2 4];


% % For Daniel's EAPA_ERP_1172018.mat
% % drgbChoicesDanielEAPA_ERP1172018
% % logP for ERP in was saved with no lag from the event
% % For New Fig. 6 ERP figure  of Daniel's paper run with which_display=1
% winNo=1;
% which_display=1;
% 
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% 
% 
% %Experiment pairs
% %Important: The first file must be the experiment performed first
% %For example in acetophenone ethyl benzoate no laser is first, laser is
% %second
% file_pairs=[
%     13 1;
%     14 2;
%     15 3;
%     16 4;
%     17 5;
%     18 6;
%     19 7;
%     20 8;
%     21 9;
%     22 10;
%     23 11;
%     24 12;
%     ];
% no_file_pairs=12;
% 
% trials_to_process=20;
% min_trials_per_event=4;
% 
% grpre=[1 3];
% grpost=[2 4];

% 
% % Odorandlaser_firstandlast_orderchange_Allin.mat
% % For  New Fig. 2 analysis using Ming's drgbChoicesMingFirstLast01102018_orderchange_allin
% % To compare naive with experienced LFP power
% winNo=2;
% refWin=1;
% which_display=3;
% 
% 
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% 
% 
% 
% % Enter the files to be processed
% files=[13 14 15 16]; %ChETA
% 
% 
% % Enter the files to be processed
% % files=[19 20]; %ChR2
% % no_files=2;
% 
% trials_to_process=30;
% min_trials_per_event=4;


% 
% % Odorandlaser_firstandlast_orderchange_Allin.mat
% % % For  New Fig. 3 analysis using Ming's drgbChoicesMingFirstLast01102018_orderchange_allin
% % % and which_display=4
% winNo=2;
% refWin=1;
% which_display=4;
% 
% 
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % trials_to_process=20;
% % output_suffix='_spmouthit.mat';
% 
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'}; %Note: This is currently not used
% trials_to_process=30;
% output_suffix='_spmout.mat';
% 
% %Experiment pairs
% % Here the last session file is the first file and the first session the second file in the pair
% % file_pairs=[file_for_last_session file_for_first_session]
% % because of this front_mask=[0 1];
% 
% % %Odor
% % file_pairs=[
% %     5 1;
% %     6 2;
% %     7 3;
% %     8 4;
% % ];
% % no_file_pairs=4;
% 
% %ChETA
% file_pairs=[
%     13 9;
%     14 10;
%     15 11;
%     16 12;
% ];
% 
% % file_pairs=[
% %     16 12
% % ];
% 
% 
% %ChR2
% % file_pairs=[
% %     19 17;
% %     20 18;
% % ];
% % no_file_pairs=2;
% 
% 
% 
% 
% front_mask=[0 1];  %For each file pair enter whether you will analyze the first few trials (1) or the last(0)
% 
% file_label{1}='proficient';
% file_label{2}='learning';
% 
% 
% min_trials_per_event=4;




% %Forward-reverse using Daniel's Fig. 1 analysis
% % drgbChoicesDanielEAPAlastandreversallast11118
% % % For Daniel's reversal EAPA_lastandreversallastEAPA11118
% % % and which_display=3
%
% winNo=2;
% refWin=1;
% which_display=3;
%
%
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % comp_window=10; %works well with 8-12
%
%
% eventType=[3 6];
% evTypeLabels={'S+','S-'}; %Note: This is currently not used
%
%
% % Enter the files to be processed
% % files=[1 2 3 4 5];  %Forward
% files=[6 7 8 9 10];  %Reverse
%
%
% no_files=5;
%
% trials_to_process=30;
% min_trials_per_event=4;



% % Forward-reverse analysis for Daniel's New Fig. 5
% % drgbChoicesDanielEAPAlastandreversallast11118
% % % For Daniel's reversal EAPA_lastandreversallastEAPA11118
% % % and which_display=3
%
% winNo=2;
% refWin=1;
% which_display=2;
%
%
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % comp_window=10; %works well with 8-12
%
%
% eventType=[3 6];
% evTypeLabels={'S+','S-'}; %Note: This is currently not used
%
% odorant{1}='EA'; %Note: odorant 1 is S+ in forward
% odorant{2}='PA';
%
% %Experiment pairs
% % Here the last session file is the first file and the first session the second file in the pair
% % file_pairs=[file_for_last_session file_for_first_session]
% % because of this front_mask=[0 1];
% file_pairs=[
%     6 1;
%     7 2;
%     8 3;
%     9 4;
%     10 5;];
% no_file_pairs=5;
%
% file_label{1}='reverse';
% file_label{2}='forward';
%
%
% front_mask=[0 0];  %For each file pair enter whether you will analyze the first few trials (1) or the last(0)
%
% trials_to_process=20;
% min_trials_per_event=4;
%
% output_suffix='_frout.mat';


% % % For Daniel's ethylacetatepropylacetateERPfirstandlast112018.mat
% % drgbChoicesDanielEAPA_ERPfirstandlast1102018
% % logP for ERP in was saved with no lag from the event
% % run with which_display=6;
%
%
% winNo=1;
% which_display=6;
%
% %
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % comp_window=10;
%
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
%
% % Enter the files to be processed
% %files=[1 2 3 4 5];
% files=[6 7 8 9 10];
% no_files=5;
%
%
% trials_to_process=30;
% min_trials_per_event=4;



% % For Daniel's EAPA_ERP_1172018.mat
% % drgbChoicesDanielEAPA_ERP1172018
% % logP for ERP in was saved with no lag from the event
% % run with which_display=6;
%
%
% winNo=1;
% which_display=6;
%
% %
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % comp_window=10;
%
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
%
% % Enter the files to be processed
% %files=[1 2 3 4 5];
% files=[21 22 23 24];
% no_files=4;
%
%
% trials_to_process=20;
% min_trials_per_event=4;
% 

% 
% %acetoethylben_ERP_shift_02032018
% %drgbChoicesDanielAPEB_ERP_shift2_1202018
% % logP for ERP in was saved with no lag from the event
% % For ERP figure  of Daniel's paper run with which_display=1
% winNo=1;
% which_display=1;
% % which_display=8;
% eventType=[2 5];
% evTypeLabels={'Hit','CR'};
% % eventType=[3 6];
% % evTypeLabels={'S+','S-'};
% 
% 
% %Experiment pairs
% %Important: The first file must be the experiment performed first
% %For example in acetophenone ethyl benzoate no laser is first, laser is
% %second
% file_pairs=[
%     7 1;
%     8 2;
%     9 3;
%     10 4;
%     11 5;
%     12 6;
%     17 13;
%     18 14;
%     19 15;
%     20 16];
% no_file_pairs=10;
% 
% trials_to_process=20;
% min_trials_per_event=3;
% shift_time=0.3;
% shift_from_event=floor(shift_time/0.025);
% 
% 
% 
% grpre=[1 3];
% grpost=[2 4];

% % drgbChoicesDanielAPEB_ERP_shift2_1202018
% % acetoethylben_ERP_shift_1232018.mat
% % logP for ERP in was saved with no lag from the event
% % run with which_display=6;
% winNo=1;
% which_display=6;
% 
% %
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % comp_window=10;
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% 
% % Enter the files to be processed
% 
% %Laser experimental
% % files=[1 2 3 4 5 6];
% 
% %No laser experimental
% files=[7 8];
% 
% trials_to_process=20;
% min_trials_per_event=3;
% shift_time=0;
% shift_from_event=floor(shift_time/0.025);
% 

<<<<<<< HEAD

=======
% 
>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855
% % drgbChoicesDaniel_shift2_IAMOfirstandlast_ERP252018.m
% % isoamylacetatemineraloilERshiftfirstandlast252018.mat
% % logP for ERP in was saved with no lag from the event
% % run with which_display=7;
% 
% winNo=1;
% which_display=7;
% 
% 
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % trials_to_process=20;
% % output_suffix='_spmouthit.mat';
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% trials_to_process=30;
% output_suffix='_spmout.mat';
% 
% 
% %Experiment pairs
% % Here the last session file is the first file and the first session the second file in the pair
% % file_pairs=[file_for_last_session file_for_first_session]
% % because of this front_mask=[0 1];
% file_pairs=[
%     1 5;
%     2 4;
%     3 6;
%     7 10;
%     8 12;
%     9 11;
%     13 14;
% ];
% 
% front_mask=[0 1];  %For each file pair enter whether you will analyze the first few trials (1) or the last(0)
% 
% file_label{1}='proficient';
% file_label{2}='learning';
% 
% min_trials_per_event=3;
% 
% shift_time=0.3;
% shift_from_event=floor(shift_time/0.025);

<<<<<<< HEAD

% drgbChoicesMingERPspikes02242018_orderchange_allin
% Odorandlaser_ERPspikes_orderchange_Allin02242018
% logP for ERP in was saved with no lag from the event
% run with which_display=7;

winNo=1;
which_display=11;


% eventType=[2 5];
% evTypeLabels={'Hit','CR'};
% trials_to_process=20;
% output_suffix='_spmouthit.mat';

eventType=[3 6];
evTypeLabels={'S+','S-'};
trials_to_process=30;
output_suffix='_spmout.mat';


%Files to analyze
% files=[53:71];  %Naive odor
% files=[39:52];  %ChR2
% files=[1:16];  %Experienced odor
files=[17:38];  %ChETA

file_label{1}='proficient';
file_label{2}='learning';

min_trials_per_event=4;

shift_time=0.3;
shift_from_event=floor(shift_time/0.025);

percent_windows=[80 100;45 65];

% no_bandwidths=4;
% low_freq=[500 1000 500 1000];
% high_freq=[1000 5000 2000 3000];
% freq_names={'500-1000','1000-5000','500-2000','1000-3000'};
    
=======
>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855
% % % For Daniel's drgbChoicesDanielAPEB_ERPshiftfirstandlast1282018
% % acetoethylbenERshiftfirstandlast1302018.mat
% % logP for ERP in was saved with no lag from the event
% % run with which_display=7;
% 
% winNo=1;
% which_display=7;
% 
% 
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % trials_to_process=20;
% % output_suffix='_spmouthit.mat';
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% trials_to_process=30;
% output_suffix='_spmout.mat';
% 
% 
% %Experiment pairs
% % Here the last session file is the first file and the first session the second file in the pair
% % file_pairs=[file_for_last_session file_for_first_session]
% % because of this front_mask=[0 1];
% file_pairs=[
%     1 13;
%     2 14;
%     3 15;
%     4 16;
%     5 18;
%     6 19;
%     7 10;
%     8 11;
%     9 12;
% ];
% 
% front_mask=[0 1];  %For each file pair enter whether you will analyze the first few trials (1) or the last(0)
% 
% file_label{1}='proficient';
% file_label{2}='learning';
% 
% min_trials_per_event=3;
% 
% shift_time=0.3;
% shift_from_event=floor(shift_time/0.025);


% 
% % % For Daniel's drgbChoicesDanielEAPA_ERPshiftfirstandlast1232018
% % ethylacetatepropylacetateERshiftfirstandlast1232018.mat
% % logP for ERP in was saved with no lag from the event
% % run with which_display=7;
% 
% winNo=1;
% which_display=7;
% 
% 
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % trials_to_process=20;
% % output_suffix='_spmouthit.mat';
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% trials_to_process=30;
% output_suffix='_spmout.mat';
% 
% 
% %Experiment pairs
% % Here the last session file is the first file and the first session the second file in the pair
% % file_pairs=[file_for_last_session file_for_first_session]
% % because of this front_mask=[0 1];
% file_pairs=[
%     1 6;
%     2 7;
%     3 8;
%     4 9;
%     5 10;
% ];
% 
% front_mask=[0 1];  %For each file pair enter whether you will analyze the first few trials (1) or the last(0)
% 
% file_label{1}='proficient';
% file_label{2}='learning';
% 
% min_trials_per_event=3;
% 
% shift_time=0.3;
% shift_from_event=floor(shift_time/0.025);


% % % For Daniel's drgbChoicesDanielEAPA_ERPshiftfirstandlast1232018
% % ethylacetatepropylacetateERshiftfirstandlast1232018.mat
% % logP for ERP in was saved with no lag from the event
% % run with which_display=6;
% winNo=1;
% which_display=6;
% 
% %
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % comp_window=10;
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% 
% % Enter the files to be processed
% 
% %Last session
% files=[1 2 3 4 5];
% 
% front_mask=0; 
% 
% trials_to_process=30;
% min_trials_per_event=3;
% shift_time=0;
% shift_from_event=floor(shift_time/0.025);



%
% % % For Justin's spmc_LFP_20180111.mat
% % drgbChoicesJustinLFP20180111
% %New Fig. 2 G2 run with which_display=3;
%
%
% winNo=2;
% refWin=1;
% which_display=3;
%
% %
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % comp_window=10;
%
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
%
% % Enter the files to be processed
% files=[3]; %[3 13 20 26 29];
% no_files=1;
%
%
% trials_to_process=30;
% min_trials_per_event=4;


% % % For Justin's spmc_LFP_20180111.mat
% % drgbChoicesJustinLFP20180111
% % logP for ERP in was saved with no lag from the event
% % run with which_display=6;
%
%
% winNo=1;
% which_display=6;
%
% %
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % comp_window=10;
%
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
%
% % Enter the files to be processed
% %files=[1 2 3 4 5];
% files=[3 13 20 26 29];
% no_files=5;
%
%
% trials_to_process=30;
% min_trials_per_event=4;


% % For Justin's spmc_LFP_20180111.mat
% %drgbChoicesJustinLFP20180111
% % New Fig. 3 run with which_display=4;
%
% winNo=2;
% refWin=1;
% which_display=4;
%
%
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % trials_to_process=20;
% % output_suffix='_spmouthit.mat';
%
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% trials_to_process=20;
% output_suffix='_spmout.mat';
%
%
% %Experiment pairs
% % Here the last session file is the first file and the first session the second file in the pair
% % file_pairs=[file_for_last_session file_for_first_session]
% % because of this front_mask=[0 1];
% file_pairs=[
%  1 3];
% no_file_pairs=1;
%
% front_mask=[0 1];  %For each file pair enter whether you will analyze the first few trials (1) or the last(0)
%
% file_label{1}='proficient';
% file_label{2}='learning';
%
% min_trials_per_event=4;

% 
% % % For Daniel's isomin_firstandlast12318
% % drgbChoicesDanielNRG1IAMOfirstandlast12318
% %New Fig. 2 G2 run with which_display=3;
% 
% 
% winNo=2;
% refWin=1;
% which_display=3;
% 
% %
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % comp_window=10;
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% 
% % Enter the files to be processed
% files=[1 2 3]; %last experimental
% output_suffix='_spmoutnrg1.mat';
% 
% % 
% % files=[7 8 9 10]; %last control
% % output_suffix='_spmoutcontrol.mat';
% 
% 
% which_electrodes=[1:4,13:16]; %Prefrontal
% % which_electrodes=[5:12]; %Hippocampus
%     
% trials_to_process=20;
% min_trials_per_event=4;



% % drgbChoicesDanielNRG1IAMOfirstandlast12318
% % % For Daniel's isomin_firstandlast12318
% %New Fig. 2 G2 run with which_display=3;
% 
% 
% winNo=2;
% refWin=1;
% which_display=3;
% 
% %
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % comp_window=10;
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% 
% % Enter the files to be processed
% % files=[1 2 3]; %last experimental
% % output_suffix='_spmoutnrg1.mat';
% 
% % 
% files=[7 8 9 10]; %last control
% output_suffix='_spmoutcontrol.mat';
% 
% 
% % which_electrodes=[1:4,13:16]; %Prefrontal
% which_electrodes=[5:12]; %Hippocampus
%     
% trials_to_process=20;
% min_trials_per_event=4;

% 
% %Compare NRG delta LFP power with controls percent>80
% % drgbChoicesDanielNRG1IAMOallfiles21218
% % % For Daniel's NRG1isominallfiles21218.mat
% winNo=2;
% refWin=1;
% which_display=10;
% % 
% % 
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% 
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% 
% % Enter the files to be processed
% files=[1:25];
% 
% groups=[ones(1,13), 2*ones(1,12)];
% group_names{1}='Experimental';
% group_names{2}='Control';
% 
% output_suffix='_spmcompnrg1.mat';
% 
% 
% which_electrodes=[1:4,13:16]; %Prefrontal
% % which_electrodes=[5:12]; %Hippocampus
% 
% min_trials_per_event=4;
% percent_windows=[80 100];


<<<<<<< HEAD

=======
% %Compare NRG ERP delta LFP power with controls percent>80
% % drgbChoicesDanielNRG1IAMOfirstandlast12318
% % % For Daniel's isomin_firstandlast12318
% 
% winNo=1;
% which_display=11;
% 
% eventType=[2 5];
% evTypeLabels={'Hit','CR'};
% % eventType=[3 6];
% % evTypeLabels={'S+','S-'};
% 
% % Enter the files to be processed
% files=[1:14];
% 
% groups=[ones(1,6), 2*ones(1,8)];
% group_names{1}='Experimental';
% group_names{2}='Control';
% 
% output_suffix='_spmcompnrg1.mat';
% 
% 
% % which_electrodes=[1:4,13:16]; %Prefrontal
% which_electrodes=[5:12]; %Hippocampus
% 
% min_trials_per_event=4;
% percent_windows=[80 100];
% 
% shift_time=0.3;
% shift_from_event=floor(shift_time/0.025);
>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855



% % drgbChoicesDanielNRG1IAMOfirstandlast12318
% % For Daniel's isomin_firstandlast12318
% % Compare low percent vs high percent
% winNo=2;
% refWin=1;
% which_display=9;
% 
% % 
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % comp_window=10;
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% 
% % Enter the files to be processed
% % files=[7:14]; %Control files
% files=[1:6]; %Experimental files
% 
% output_suffix='_spmpernrg1.mat';
% 
% 
% which_electrodes=[1:4,13:16]; %Prefrontal
% % which_electrodes=[5:12]; %Hippocampus
%     
% min_trials_per_event=4;
% percent_windows=[80 100;45 65];
% 
% file_label{1}='proficient';
% file_label{2}='learning';



% % % drgbChoicesDanielNRG1IAMOfirstandlast12318
% % % % For Daniel's isomin_firstandlast12318
% 
% winNo=2;
% refWin=1;
% which_display=4;
% 
% 
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % trials_to_process=20;
% % output_suffix='_spmouthit.mat';
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% trials_to_process=30;
% output_suffix='_spmout.mat';
% 
% 
% %Experiment pairs
% % Here the last session file is the first file and the first session the second file in the pair
% % file_pairs=[file_for_last_session file_for_first_session]
% % because of this front_mask=[0 1];
% file_pairs=[
%     1 11;
%     2 7;
%     3 8;
%     4 9;
%     5 10;
% %     6 12; %Exclude,159867 should NOT be used because the mouse was not
% %     proficient in the "last" session
%     13 16;
%     14 17;
%     15 18];
% no_file_pairs=8;
% 
% front_mask=[0 0];  %For each file pair enter whether you will analyze the first few trials (1) or the last(0)
% 
% file_label{1}='Control';
% file_label{2}='NRG1-IV';
% 
% which_electrodes=[1:4,13:16]; %Prefrontal
% % which_electrodes=[5:12]; %Hippocampus
% 
% min_trials_per_event=4;



% % For Daniel's EAPA_ERP_shift2_1252018.mat
% % drgbChoicesDaniel_shift2_EAPA_ERP1252018.m
% % logP for ERP in was saved with no lag from the event
% % For ERP figure  of Daniel's paper run with which_display=1
% winNo=1;
% % which_display=1;
% which_display=8;
% eventType=[2 5];
% evTypeLabels={'Hit','CR'};
% % eventType=[3 6];
% % evTypeLabels={'S+','S-'};
% 
% 
% %Experiment pairs
% %Important: The first file must be the experiment performed first
% %For example in acetophenone ethyl benzoate no laser is first, laser is
% %second
% file_pairs=[
%     13 1;
%     14 2;
%     15 4;
%     16 5;
%     17 6;
%     18 7;
%     19 8;
%     20 9;
%     21 10;
%     22 11;
%     23 12;
%     24 3];
% no_file_pairs=12;
% 
% trials_to_process=20;
% min_trials_per_event=3;
% shift_time=0.3;
% shift_from_event=floor(shift_time/0.025);
% 
% 
% grpre=[1 3];
% grpost=[2 4];

% 
% % % For Daniel's isoamylacetatemineraloilERshiftfirstandlast1272018.mat
% % drgbChoicesDanielIAMO_ERPshiftfirstandlast1272018
% % logP for ERP in was saved with no lag from the event
% % run with which_display=7;
% 
% winNo=1;
% which_display=7;
% 
% 
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% % trials_to_process=20;
% % output_suffix='_spmouthit.mat';
% 
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
% trials_to_process=30;
% output_suffix='_ERPspmout.mat';
% 
% 
% %Experiment pairs
% % Here the last session file is the first file and the first session the second file in the pair
% % file_pairs=[file_for_last_session file_for_first_session]
% % because of this front_mask=[0 1];
% file_pairs=[
%     1 4;
%     2 5;
%     3 6;
%     7 10;
%     8 11;
%     9 12;
%     13 14;
% ];
% no_file_pairs=7;
% 
% front_mask=[0 1];  %For each file pair enter whether you will analyze the first few trials (1) or the last(0)
% 
% file_label{1}='proficient';
% file_label{2}='learning';
% 
% min_trials_per_event=4;
% shift_time=0.3;
% shift_from_event=floor(shift_time/0.025);



% %drgbChoicesDaniel_shift2_IAMO_ERP1312018
% %IAMO_ERP_shift2_02052018.mat
% % logP for ERP in was saved with no lag from the event
% % For ERP figure  of Daniel's paper run with which_display=1
% winNo=1;
% % which_display=1;
% which_display=8;
% eventType=[2 5];
% evTypeLabels={'Hit','CR'};
% % eventType=[3 6];
% % evTypeLabels={'S+','S-'};
% 
% 
% %Experiment pairs
% %Important: The first file must be the experiment performed first
% %For example in acetophenone ethyl benzoate no laser is first, laser is
% %second
% file_pairs=[
%     10 1;
%     11 2;
%     12 3;
%     13 4;
%     14 5;
%     15 6;
%     16 7;
%     17 8;
%     18 9];
% 
% trials_to_process=20;
% min_trials_per_event=3;
% shift_time=0.3;
% shift_from_event=floor(shift_time/0.025);
% 
% 
% grpre=[1 3];
% grpost=[2 4];

%% The code processing pairwise batch LFP starts here

close all
warning('off')


no_event_types=length(eventType);

%Statistics
%stats_method = 1, FDR, 2, statscond
stat_method=1;

%Mode for percent ANOVA
mode_statcond='perm';
% mode_statcond='bootstrap';

%Which percent correct bins do you want to use?
percent_low=[45 65 80];
percent_high=[65 80 100];
percent_bin_legend={' learning';' between learning and proficient';' proficient'};
no_percent_bins=3;

%Bandwidths
<<<<<<< HEAD
if exist('no_bandwidths')==0
    no_bandwidths=4;
    low_freq=[6 15 35 65];
    high_freq=[12 30 55 95];
    freq_names={'Theta','Beta','Low gamma','High gamma'};
end
=======
no_bandwidths=4;
low_freq=[6 15 35 65];
high_freq=[12 30 55 95];
freq_names={'Theta','Beta','Low gamma','High gamma'};
>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855

event1=eventType(1);
event2=eventType(2);

%To calculate percent correct
perevent1=2;
perevent2=5;


%Ask user for the drgb output .mat file and load those data
[handles.drgb.outFileName,handles.PathName] = uigetfile('*.mat','Select the drgb output file');
load([handles.PathName handles.drgb.outFileName])


fprintf(1, ['\ndrgDisplayBatchLFPPowerPairwise run for ' handles.drgb.outFileName '\nwhich_display= = %d\n\n'],which_display);

switch which_display
<<<<<<< HEAD
    case {1,6,7,8,11}
=======
    case {1,6,7,8}
>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855
        frequency=handles_drgb.drgb.lfpevpair(1).fERP;
        max_events_per_sec=(handles_drgb.drgbchoices.timeEnd(winNo)-handles_drgb.drgbchoices.timeStart(winNo))*handles_drgb.max_events_per_sec;
    otherwise
        frequency=handles_drgb.drgb.freq_for_LFPpower;
end



figNo=0;

%These are the colors for the different lines

these_colors{1}='b';
these_colors{2}='r';
these_colors{3}='m';
these_colors{8}='g';
these_colors{5}='y';
these_colors{6}='k';
these_colors{7}='c';
these_colors{4}='k';

these_lines{1}='-b';
these_lines{2}='-r';
these_lines{3}='-m';
these_lines{8}='-g';
these_lines{5}='-y';
these_lines{6}='-k';
these_lines{7}='-c';
these_lines{4}='-k';

%Initialize the variables
%Get files and electrode numbers
for lfpodNo=1:handles_drgb.drgb.lfpevpair_no
    files_per_lfp(lfpodNo)=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
    elec_per_lfp(lfpodNo)=handles_drgb.drgb.lfpevpair(lfpodNo).elecNo;
    window_per_lfp(lfpodNo)=handles_drgb.drgb.lfpevpair(lfpodNo).timeWindow;
end

switch which_display
    case 1
        %Compare auROC for ERP LFP in the last few trials of pre with first few trials of post
        %Used for New Fig. 5 of Daniel's paper
        no_dBs=1;
        delta_dB_power_pre=[];
        no_ROCs=0;
        ROCoutpre=[];
        ROCoutpost=[];
        p_vals_ROC=[];
        delta_dB_powerpreHit=[];
        no_hits=0;
        perCorr_pre=[];
        perCorr_post=[];
        group_pre=[];
        group_post=[];
        shift_ii=floor(length(handles_drgb.drgb.lfpevpair(1).out_times)/2)+1+shift_from_event;
        
        
        fprintf(1, ['Pairwise auROC analysis for ERP LFP power\n\n'],'perCorr_pre','perCorr_post')
        p_vals=[];
        for fps=1:no_file_pairs
            for elec=1:16
                
                lfpodNopre=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                lfpodNopost=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                
%                 if elec==1
%                     %Find percent correct for pre
%                     trials_Sp=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(event1,:)==1);
%                     trials_Sm=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(event2,:)==1);
%                     trials_Hit=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(perevent1,:)==1);
%                     trials_CR=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(perevent2,:)==1);
%                     perCorr_pre(fps)=100*(trials_Hit+trials_CR)/(trials_Sp+trials_Sm);
%                     group_pre(fps)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopre).fileNo);
%                     
%                     %Find percent correct for post
%                     trials_Sp=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(event1,:)==1);
%                     trials_Sm=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(event2,:)==1);
%                     trials_Hit=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(perevent1,:)==1);
%                     trials_CR=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(perevent2,:)==1);
%                     perCorr_post(fps)=100*(trials_Hit+trials_CR)/(trials_Sp+trials_Sm);
%                     group_post(fps)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopost).fileNo);
%                     
%                     
%                     
%                     fprintf(1, '\nPercent correct for session pair %d pre= %d, post= %d\n',fps,perCorr_pre(fps),perCorr_post(fps));
%                     
%                 end
                
                if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre)))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost)))
                    
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre).log_P_tERP))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost).log_P_tERP))
                        
                        if (length(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(1,:))>=trials_to_process) &...
                                (length(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(1,:))>=trials_to_process)
                            
                            length_pre=length(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(1,:));
                            pre_mask=logical([zeros(1,length_pre-trials_to_process) ones(1,trials_to_process)]);
                            trials_in_event_preHit=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(event1,:)==1);
                            trials_in_event_preCR=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(event2,:)==1);
                            
                            trials_with_event_pre=(handles_drgb.drgb.lfpevpair(lfpodNopre).no_events_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNopre).no_events_per_trial<=max_events_per_sec)...
                                &(handles_drgb.drgb.lfpevpair(lfpodNopre).no_ref_evs_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNopre).no_ref_evs_per_trial<=max_events_per_sec);
                            
                            
                            length_post=length(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(1,:));
                            post_mask=logical([ones(1,trials_to_process) zeros(1,length_post-trials_to_process)]);
                            trials_in_event_postHit=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(event1,:)==1);
                            trials_in_event_postCR=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(event2,:)==1);
                            
                            trials_with_event_post=(handles_drgb.drgb.lfpevpair(lfpodNopost).no_events_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNopost).no_events_per_trial<=max_events_per_sec)...
                                &(handles_drgb.drgb.lfpevpair(lfpodNopost).no_ref_evs_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNopost).no_ref_evs_per_trial<=max_events_per_sec);
                            
                            
                            if (sum(trials_in_event_preHit&trials_with_event_pre&pre_mask)>=min_trials_per_event) & (sum(trials_in_event_preCR&trials_with_event_pre&pre_mask)>=min_trials_per_event) & ...
                                    (sum(trials_in_event_postHit&trials_with_event_post&post_mask)>=min_trials_per_event) & (sum(trials_in_event_postCR&trials_with_event_post&post_mask)>=min_trials_per_event)
                                
                                
                                %pre Hits
                                this_dB_powerpreHit=zeros(sum(trials_in_event_preHit&pre_mask),length(frequency));
                                this_dB_powerpreHit(:,:)=handles_drgb.drgb.lfpevpair(lfpodNopre).log_P_tERP(trials_in_event_preHit&pre_mask,:,shift_ii);
                                
                                
                                %pre CRs
                                this_dB_powerpreCR=zeros(sum(trials_in_event_preCR&pre_mask),length(frequency));
                                this_dB_powerpreCR(:,:)=handles_drgb.drgb.lfpevpair(lfpodNopre).log_P_tERP(trials_in_event_preCR&pre_mask,:,shift_ii);
                                
                                
                                %post Hits
                                this_dB_powerpostHit=zeros(sum(trials_in_event_postHit&post_mask),length(frequency));
                                this_dB_powerpostHit(:,:)=handles_drgb.drgb.lfpevpair(lfpodNopost).log_P_tERP(trials_in_event_postHit&post_mask,:,shift_ii);
                                
                                
                                %post CRs
                                this_dB_powerpostCR=zeros(sum(trials_in_event_postCR&post_mask),length(frequency));
                                this_dB_powerpostCR(:,:)=handles_drgb.drgb.lfpevpair(lfpodNopost).log_P_tERP(trials_in_event_postCR&post_mask,:,shift_ii);
                                
                                for bwii=1:no_bandwidths
                                    
                                    no_ROCs=no_ROCs+1;
                                    this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                    
                                    %Enter the pre Hits
                                    this_delta_dB_powerpreHit=zeros(sum(trials_in_event_preHit&pre_mask),1);
                                    this_delta_dB_powerpreHit=mean(this_dB_powerpreHit(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:sum(trials_in_event_preHit&pre_mask),1)=this_delta_dB_powerpreHit;
                                    roc_data(1:sum(trials_in_event_preHit&pre_mask),2)=zeros(sum(trials_in_event_preHit&pre_mask),1);
                                    
                                    %Enter pre CR
                                    total_trials=sum(trials_in_event_preHit&pre_mask)+sum(trials_in_event_preCR&pre_mask);
                                    this_delta_dB_powerpreCR=zeros(sum(trials_in_event_preCR&pre_mask),1);
                                    this_delta_dB_powerpreCR=mean(this_dB_powerpreCR(:,this_band),2);
                                    roc_data(sum(trials_in_event_preHit&pre_mask)+1:total_trials,1)=this_delta_dB_powerpreCR;
                                    roc_data(sum(trials_in_event_preHit&pre_mask)+1:total_trials,2)=ones(sum(trials_in_event_preCR&pre_mask),1);
                                    
                                    
                                    %Find pre ROC
                                    ROCoutpre(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCoutpre(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopre).fileNo;
                                    ROCgroupNopre(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopre).fileNo);
                                    ROCoutpre(no_ROCs).timeWindow=winNo;
                                    ROCbandwidthpre(no_ROCs)=bwii;
                                    auROCpre(no_ROCs)=ROCoutpre(no_ROCs).roc.AUC-0.5;
                                    p_valROCpre(no_ROCs)=ROCoutpre(no_ROCs).roc.p;
                                    
                                    p_vals_ROC=[p_vals_ROC ROCoutpre(no_ROCs).roc.p];
                                    
                                    %Enter the post Hits
                                    this_delta_dB_powerpostHit=zeros(sum(trials_in_event_postHit&post_mask),1);
                                    this_delta_dB_powerpostHit=mean(this_dB_powerpostHit(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:sum(trials_in_event_postHit&post_mask),1)=this_delta_dB_powerpostHit;
                                    roc_data(1:sum(trials_in_event_postHit&post_mask),2)=zeros(sum(trials_in_event_postHit&post_mask),1);
                                    
                                    %Enter post CR
                                    total_trials=sum(trials_in_event_postHit&post_mask)+sum(trials_in_event_postCR&post_mask);
                                    this_delta_dB_powerpostCR=zeros(sum(trials_in_event_postCR&post_mask),1);
                                    this_delta_dB_powerpostCR=mean(this_dB_powerpostCR(:,this_band),2);
                                    roc_data(sum(trials_in_event_postHit&post_mask)+1:total_trials,1)=this_delta_dB_powerpostCR;
                                    roc_data(sum(trials_in_event_postHit&post_mask)+1:total_trials,2)=ones(sum(trials_in_event_postCR&post_mask),1);
                                    
                                    
                                    %Find post ROC
                                    ROCoutpost(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCoutpost(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopost).fileNo;
                                    ROCgroupNopost(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopost).fileNo);
                                    ROCoutpost(no_ROCs).timeWindow=winNo;
                                    ROCbandwidthpost(no_ROCs)=bwii;
                                    auROCpost(no_ROCs)=ROCoutpost(no_ROCs).roc.AUC-0.5;
                                    p_valROCpost(no_ROCs)=ROCoutpost(no_ROCs).roc.p;
                                    
                                    p_vals_ROC=[p_vals_ROC ROCoutpost(no_ROCs).roc.p];
                                    
                                    if (auROCpost(no_ROCs)<0.3)&(auROCpre(no_ROCs)>0.4)&(ROCgroupNopre(no_ROCs)==1)&(ROCbandwidthpre(no_ROCs)==2)
                                        fprintf(1, ['Decrease in auROC for file No %d vs file No %d electrode %d bandwidth No: %d\n'],file_pairs(fps,1),file_pairs(fps,2),elec,bwii);
                                    end
                                    
                                    %Are the delta dB LFP's different?
                                    
                                    %Hit
                                    p_val(no_dBs,bwii)=ranksum(this_delta_dB_powerpreHit,this_delta_dB_powerpostHit);
                                    p_vals=[p_vals p_val(no_dBs,bwii)];
                                    groupNopre(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                    groupNopost(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                    events(no_dBs)=1;
                                    
                                    
                                    %CR
                                    p_val(no_dBs+1,bwii)=ranksum(this_delta_dB_powerpreCR,this_delta_dB_powerpostCR);
                                    p_vals=[p_vals p_val(no_dBs+1,bwii)];
                                    groupNopre(no_dBs+1)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                    groupNopost(no_dBs+1)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                    events(no_dBs+1)=2;
                                    
                                    if p_val(no_dBs,bwii)<0.05
                                        dB_power_changeHit(no_ROCs)=1;
                                    else
                                        dB_power_changeHit(no_ROCs)=0;
                                    end
                                    
                                    if p_val(no_dBs+1,bwii)<0.05
                                        dB_power_changeCR(no_ROCs)=1;
                                    else
                                        dB_power_changeCR(no_ROCs)=0;
                                    end
                                    
                                    %Plot the points and save the data
                                    if groupNopre(no_dBs)==1
                                        
                                        
                                        %Hit, all points
                                        delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                        delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                        
                                        
                                        %CR, all points
                                        delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                        delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                        
                                    else
                                        if groupNopre(no_dBs)==3
                                            
                                            %Hit, all points
                                            delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                            delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                            
                                            
                                            %CR, all points
                                            delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                            delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                            %                                             figure(bwii+4+12)
                                            %                                             hold on
                                            %                                             plot([3 4],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        end
                                    end
                                end
                                
                                no_dBs=no_dBs+2;
                                
                            else
                                
                                if (sum(trials_in_event_preHit&trials_with_event_pre&pre_mask)<min_trials_per_event)
                                    fprintf(1, ['%d trials with lick events in ' evTypeLabels{find(eventType==event1)} ' fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_preHit&trials_with_event_pre&pre_mask), min_trials_per_event,file_pairs(fps,1),elec);
                                end
                                
                                if (sum(trials_in_event_preCR&trials_with_event_pre&pre_mask)<min_trials_per_event)
                                    fprintf(1, ['%d trials with lick events in ' evTypeLabels{find(eventType==event2)} ' fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_preCR&trials_with_event_pre&pre_mask),event2, min_trials_per_event,file_pairs(fps,1),elec);
                                end
                                
                                if (sum(trials_in_event_postHit&trials_with_event_post&post_mask)<min_trials_per_event)
                                    fprintf(1, ['%d trials with lick events in ' evTypeLabels{find(eventType==event1)} ' fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_postHit&trials_with_event_post&post_mask),event1, min_trials_per_event,file_pairs(fps,2),elec);
                                end
                                
                                if (sum(trials_in_event_postCR&trials_with_event_post&post_mask)<min_trials_per_event)
                                    fprintf(1, ['%d trials with lick events in ' evTypeLabels{find(eventType==event2)} ' fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_postCR&trials_with_event_post&post_mask),event2, min_trials_per_event,file_pairs(fps,2),elec);
                                end
                                
                            end
                            
                        else
                            
                            if (length(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(1,:))<trials_to_process)
                                fprintf(1, ['%d trials fewer than %d trials to process for file No %d electrode %d\n'],length(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(1,:)),trials_to_process,file_pairs(fps,1),elec);
                            end
                            
                            if (length(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(1,:))<trials_to_process)
                                fprintf(1, ['%d trials fewer than %d trials to process for file No %d electrode %d\n'],length(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(1,:)),trials_to_process,file_pairs(fps,1),elec);
                            end
                            
                        end
                    else
                        
                        if isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).log_P_tERP)
                            fprintf(1, ['Empty log_P_tERP for file No %d electrode %d\n'],file_pairs(fps,1),elec);
                        end
                        
                        if isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).log_P_tERP)
                            fprintf(1, ['Empty log_P_tERP for file No %d electrode %d\n'],file_pairs(fps,2),elec);
                        end
                        
                    end
                  
                else
                    
                    if isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref))
                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],file_pairs(fps,1),elec);
                    end
                    
                    if isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref))
                        fprintf(1, ['Empty lfpevpairfor file No %d electrode %d\n'],file_pairs(fps,2),elec);
                    end
                    
                end
            end
            
        end
        fprintf(1, '\n\n')
        
        
        pFDRROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for significant difference of auROC p value from 0.5  = %d\n\n'],pFDRROC);
        
        
        %Now plot the bar graphs and do anovan for LFP power
        p_vals_anovan=[];
        pvals_ancova=[];
        pvals_auROCancova=[];
        for bwii=1:4
            
            
            %Do ancova for auROC auROCpre
            this_auROCpre=[];
            this_auROCpre=auROCpre((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))';
            this_auROCpre=[this_auROCpre; auROCpre((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))'];
            
            
            this_auROCpost=[];
            this_auROCpost=auROCpost((ROCgroupNopre==1)&(ROCbandwidthpost==bwii))';
            this_auROCpost=[this_auROCpost; auROCpost((ROCgroupNopre==3)&(ROCbandwidthpost==bwii))'];
            
            pre_post=[];
            pre_post=[zeros(sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)),1); ones(sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)),1)];
            
            
            [h,atab,ctab,stats] = aoctool(this_auROCpre,this_auROCpost,pre_post,0.05,'','','','off');
            
            
            pvals_auROCancova=[pvals_auROCancova atab{4,6}];
            fprintf(1, ['ancova auROC p value ' freq_names{bwii} ' = %d\n\n'],atab{4,6});
            
            %Do ancova figure for auROC
            figure(10+bwii)
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            h1=plot(auROCpre((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)),auROCpost((ROCgroupNopre==1)&(ROCbandwidthpost==bwii)),'or','MarkerFace','r');
            hold on
            h2=plot(auROCpre((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)),auROCpost((ROCgroupNopre==3)&(ROCbandwidthpost==bwii)),'ob','MarkerFace','b');
            
            slope_pre=ctab{5,2}+ctab{6,2};
            int_pre=ctab{2,2}+ctab{3,2};
            min_x=min([min(auROCpre((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) min(auROCpre((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
            max_x=max([max(auROCpre((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) max(auROCpre((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
            x=[-0.2 0.5];
            plot(x,slope_pre*x+int_pre,'-r','LineWidth',2)
            
            slope_post=ctab{5,2}+ctab{7,2};
            int_post=ctab{2,2}+ctab{4,2};
            x=[-0.2 0.5];
            plot(x,slope_post*x+int_post,'-b','LineWidth',2)
            
            plot([-0.2 0.5],[-0.2 0.5],'-k','LineWidth',2)
            
            title(['post vs pre auROC for ' freq_names{bwii} ])
            xlabel('pre auROC')
            ylabel('post auROC')
            legend([h1 h2],'halo','no halo')
            xlim([-0.2 0.5])
            ylim([-0.2 0.5])
            ax=gca;
            ax.LineWidth=3;
        end
        
        %         pFDRanovan=drsFDRpval(p_vals_anovan);
        %         fprintf(1, ['pFDR for anovan p value  = %d\n\n'],pFDRanovan);
        %
        %         pFDRancova=drsFDRpval(pvals_ancova);
        %         fprintf(1, ['pFDR for power dB ancova p value  = %d\n\n'], pFDRancova);
        
        pFDRauROCancova=drsFDRpval(pvals_auROCancova);
        fprintf(1, ['pFDR for auROC ancova p value  = %d\n\n'], pFDRauROCancova);
        
        fprintf(1, '\n\n')
        
        %         p_chi=[];
        %         for evTN1=1:length(eventType)
        %             fprintf(1, ['Significant changes in pairwise LFP power analysis for event: ' evTypeLabels{evTN1} '\n\n'])
        %             for bwii=1:4
        %                 for grs=grpre
        %                     num_sig(grs)=sum(p_val((events==evTN1)&(groupNopre==grs),bwii)<=0.05);
        %                     tot_num(grs)=sum((events==evTN1)&(grs==groupNopre));
        %                     fprintf(1, ['Number significant for ' freq_names{bwii} ' and ' handles_drgb.drgbchoices.group_no_names{grs} ' = %d of %d\n'],num_sig(grs),tot_num(grs));
        %                 end
        %                 [p, Q]= chi2test([num_sig(grpre(1)), tot_num(grpre(1))-num_sig(grpre(1)); num_sig(grpre(2)), tot_num(grpre(2))-num_sig(grpre(2))]);
        %                 fprintf(1, ['Chi squared p value  = %d\n\n'],p);
        %                 p_chi=[p_chi p];
        %             end
        %             fprintf(1, '\n\n\n')
        %         end
        %
        %         pFDRchi=drsFDRpval(p_chi);
        %         fprintf(1, ['pFDR for Chi squared p value  = %d\n\n'],pFDRchi);
        
        %Plot cumulative histos for auROCs
        dB_power_change=logical(dB_power_changeHit+dB_power_changeCR);
        figNo=0;
        p_val_ROC=[];
        pvals_auROCperm=[];
        
        
        for bwii=1:4
            n_cum=0;
            this_legend=[];
            data_auROC=[];
            pre_post_auROC=[];
            gr_auROC=[];
            for grs=1:2
                if grs==1
                    try
                        close(figNo+1)
                    catch
                    end
                    figure(figNo+1)
                else
                    try
                        close(figNo+2)
                    catch
                    end
                    figure(figNo+2)
                end
                hold on
                
                %Plot the histograms
                maxauROC=max([max(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))) max(auROCpost((ROCgroupNopost==grpost(grs))&(ROCbandwidthpost==bwii)))]);
                minauROC=min([min(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))) min(auROCpost((ROCgroupNopost==grpost(grs))&(ROCbandwidthpost==bwii)))]);
                edges=[-0.5:0.05:0.5];
                pos2=[0.1 0.1 0.6 0.8];
                subplot('Position',pos2)
                set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                hold on
                
                h2=histogram(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)),edges);
                h2.FaceColor='b';
                h1=histogram(auROCpost((ROCgroupNopost==grpost(grs))&(ROCbandwidthpost==bwii)),edges);
                h1.FaceColor='r';
                
                xlabel('auROC')
                ylabel('# of electrodes')
                legend('Pre','Laser')
                if grs==1
                    title(['auROC DBh Cre x halo for ' freq_names{bwii}])
                else
                    title(['auROC DBh Cre for ' freq_names{bwii}])
                end
                xlim([-0.3 0.6])
                ylim([0 40])
                ax=gca;
                ax.LineWidth=3;
                %                 if grs==1
                %                     ylim([0 30])
                %                 else
                %                     ylim([0 40])
                %                 end
                
                %Plot the single electrodes
                pos2=[0.8 0.1 0.1 0.8];
                subplot('Position',pos2)
                hold on
                for ii=1:length(auROCpre)
                    if (ROCgroupNopre(ii)==grpre(grs))&(ROCbandwidthpre(ii)==bwii)
                        plot([0 1],[auROCpre(ii) auROCpost(ii)],'-o', 'Color',[0.7 0.7 0.7])
                    end
                end
                
                
                plot([0 1],[mean(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))) mean(auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))],'-k','LineWidth', 3)
                CI = bootci(1000, @mean, auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                plot([0 0],CI,'-b','LineWidth',3)
                plot(0,mean(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))),'ob','MarkerSize', 10,'MarkerFace','b')
                CI = bootci(1000, @mean, auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                plot([1 1],CI,'-r','LineWidth',3)
                plot(1,mean(auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))),'or','MarkerSize', 10,'MarkerFace','r')
                ylabel('auROC')
                ylim([-0.2 0.5])
                ax=gca;
                ax.LineWidth=3;
                %Do the statistics for auROC differences
                %                 a={auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))' auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))'};
                %                 mode_statcond='perm';
                %                 [F df pval_auROCperm] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
                %
                pval_auROCperm=ranksum(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)), auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                
                if grs==1
                    fprintf(1, ['p value for premuted anovan for auROC DBH Cre x halo pre vs laser ' freq_names{bwii} '= %d\n'],  pval_auROCperm);
                else
                    fprintf(1, ['p value for premuted anovan for auROC DBH Cre pre vs laser ' freq_names{bwii} '= %d\n'],  pval_auROCperm);
                end
                pvals_auROCperm=[pvals_auROCperm pval_auROCperm];
                
                %Save the data for anovan interaction
                %Pre
                data_auROC=[data_auROC auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))];
                gr_auROC=[gr_auROC grs*ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
                pre_post_auROC=[pre_post_auROC ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
                
                %Post
                data_auROC=[data_auROC auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))];
                gr_auROC=[gr_auROC grs*ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
                pre_post_auROC=[pre_post_auROC 2*ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
            end
            figNo=figNo+2;
            x=x+3;
            
            %Calculate anovan for inteaction
            [p,tbl,stats]=anovan(data_auROC,{pre_post_auROC gr_auROC},'model','interaction','varnames',{'pre_vs_post','halo_vs_no_halo'},'display','off');
            fprintf(1, ['p value for anovan auROC interaction for ' freq_names{bwii} '= %d\n'],  p(3));
            p_aovan_int(bwii)=p(3);
            
        end
        
        pFDRauROC=drsFDRpval(pvals_auROCperm);
        fprintf(1, ['pFDR for auROC  = %d\n\n'],pFDRauROC);
        
        pFDRauROCint=drsFDRpval(p_aovan_int);
        fprintf(1, ['pFDR for auROC anovan interaction  = %d\n\n'],pFDRauROCint);
        

        save([handles.PathName handles.drgb.outFileName(1:end-4) '_out.mat'],'perCorr_pre','perCorr_post','group_pre', 'group_post');
        pfft=1;
        
    case 2
        %Compare auROC in the last few trials of the last session file with
        %first few trials of session
        % Generate figure 2 for Daniel's paper. first vs last.
        no_dBs=1;
        delta_dB_power_fp1=[];
        no_ROCs=0;
        ROCoutfp1=[];
        ROCoutfp2=[];
        p_vals_ROC=[];
        delta_dB_powerfp1Ev1=[];
        no_Ev1=0;
        pvals_auROCperm=[];
        pvals_dBperm=[];
        perCorr_fp1=[];
        perCorr_fp2=[];
        
        fprintf(1, ['Pairwise auROC analysis for ' evTypeLabels{1} ' and ' evTypeLabels{2} ' LFP power\n\n'])
        p_vals=[];
        for fps=1:no_file_pairs
            
            
            for elec=1:16
                
                lfpodNofp1_ref=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                lfpodNofp2_ref=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                
                if elec==1
                    %Find percent correct for fp1 block
                    perCorr_fp1(fps)=handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).perCorrLFPPower(end);
                    
                    %Find percent correct for fp2 block
                    perCorr_fp2(fps)=handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).perCorrLFPPower(1);
                    
                    fprintf(1, '\nPercent correct for session pair %d last= %d, first= %d\n',fps,perCorr_fp1(fps),perCorr_fp2(fps));
                    
                end
                
                if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref)))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref)))
                    
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).allPower))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).allPower))
                        
                        if (length(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).which_eventLFPPower(1,:))>=trials_to_process) &...
                                (length(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).which_eventLFPPower(1,:))>=trials_to_process)
                            
                            length_fp1=length(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).which_eventLFPPower(1,:));
                            if front_mask(1)==1
                                fp1_mask=logical([ones(1,trials_to_process) zeros(1,length_fp1-trials_to_process)]);
                            else
                                fp1_mask=logical([zeros(1,length_fp1-trials_to_process) ones(1,trials_to_process)]);
                            end
                            trials_in_event_fp1Ev1=(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).which_eventLFPPower(event1,:)==1);
                            trials_in_event_fp1Ev2=(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).which_eventLFPPower(event2,:)==1);
                            
                            length_fp2=length(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).which_eventLFPPower(1,:));
                            if front_mask(2)==1
                                fp2_mask=logical([ones(1,trials_to_process) zeros(1,length_fp2-trials_to_process)]);
                            else
                                fp2_mask=logical([zeros(1,length_fp2-trials_to_process) ones(1,trials_to_process)]);
                            end
                            trials_in_event_fp2Ev1=(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).which_eventLFPPower(event1,:)==1);
                            trials_in_event_fp2Ev2=(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).which_eventLFPPower(event2,:)==1);
                            
                            
                            if (sum(trials_in_event_fp1Ev1)>=min_trials_per_event) & (sum( trials_in_event_fp1Ev2)>=min_trials_per_event) & ...
                                    (sum(trials_in_event_fp2Ev1)>=min_trials_per_event) & (sum(trials_in_event_fp2Ev2)>=min_trials_per_event)
                                
                                lfpodNofp1=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                lfpodNofp2=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                
                                %fp1 Ev1
                                this_dB_powerfp1refEv1=zeros(sum(trials_in_event_fp1Ev1&fp1_mask),length(frequency));
                                this_dB_powerfp1refEv1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).allPower(trials_in_event_fp1Ev1&fp1_mask,:));
                                
                                this_dB_powerfp1Ev1=zeros(sum(trials_in_event_fp1Ev1&fp1_mask),length(frequency));
                                this_dB_powerfp1Ev1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp1).allPower(trials_in_event_fp1Ev1&fp1_mask,:));
                                
                                %fp1 Ev2
                                this_dB_powerfp1refEv2=zeros(sum(trials_in_event_fp1Ev2&fp1_mask),length(frequency));
                                this_dB_powerfp1refEv2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).allPower(trials_in_event_fp1Ev2&fp1_mask,:));
                                
                                this_dB_powerfp1Ev2=zeros(sum(trials_in_event_fp1Ev2&fp1_mask),length(frequency));
                                this_dB_powerfp1Ev2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp1).allPower(trials_in_event_fp1Ev2&fp1_mask,:));
                                
                                
                                %fp2 Ev1
                                this_dB_powerfp2refEv1=zeros(sum(trials_in_event_fp2Ev1&fp2_mask),length(frequency));
                                this_dB_powerfp2refEv1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).allPower(trials_in_event_fp2Ev1&fp2_mask,:));
                                
                                this_dB_powerfp2Ev1=zeros(sum(trials_in_event_fp2Ev1&fp2_mask),length(frequency));
                                this_dB_powerfp2Ev1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp2).allPower(trials_in_event_fp2Ev1&fp2_mask,:));
                                
                                %fp2 Ev2
                                this_dB_powerfp2refEv2=zeros(sum(trials_in_event_fp2Ev2&fp2_mask),length(frequency));
                                this_dB_powerfp2refEv2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).allPower(trials_in_event_fp2Ev2&fp2_mask,:));
                                
                                this_dB_powerfp2Ev2=zeros(sum(trials_in_event_fp2Ev2&fp2_mask),length(frequency));
                                this_dB_powerfp2Ev2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp2).allPower(trials_in_event_fp2Ev2&fp2_mask,:));
                                
                                for bwii=1:no_bandwidths
                                    
                                    no_ROCs=no_ROCs+1;
                                    this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                    
                                    %Enter the fp1 Ev1
                                    this_delta_dB_powerfp1Ev1=zeros(sum(trials_in_event_fp1Ev1&fp1_mask),1);
                                    this_delta_dB_powerfp1Ev1=mean(this_dB_powerfp1Ev1(:,this_band)-this_dB_powerfp1refEv1(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:sum(trials_in_event_fp1Ev1&fp1_mask),1)=this_delta_dB_powerfp1Ev1;
                                    roc_data(1:sum(trials_in_event_fp1Ev1&fp1_mask),2)=zeros(sum(trials_in_event_fp1Ev1&fp1_mask),1);
                                    
                                    %Enter fp1 Ev2
                                    total_trials=sum(trials_in_event_fp1Ev1&fp1_mask)+sum(trials_in_event_fp1Ev2&fp1_mask);
                                    this_delta_dB_powerfp1Ev2=zeros(sum(trials_in_event_fp1Ev2&fp1_mask),1);
                                    this_delta_dB_powerfp1Ev2=mean(this_dB_powerfp1Ev2(:,this_band)-this_dB_powerfp1refEv2(:,this_band),2);
                                    roc_data(sum(trials_in_event_fp1Ev1&fp1_mask)+1:total_trials,1)=this_delta_dB_powerfp1Ev2;
                                    roc_data(sum(trials_in_event_fp1Ev1&fp1_mask)+1:total_trials,2)=ones(sum(trials_in_event_fp1Ev2&fp1_mask),1);
                                    
                                    
                                    %Find fp1 ROC
                                    ROCoutfp1(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCoutfp1(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).fileNo;
                                    ROCgroupNofp1(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).fileNo);
                                    ROCoutfp1(no_ROCs).timeWindow=winNo;
                                    ROCbandwidthfp1(no_ROCs)=bwii;
                                    auROCfp1(no_ROCs)=ROCoutfp1(no_ROCs).roc.AUC-0.5;
                                    p_valROCfp1(no_ROCs)=ROCoutfp1(no_ROCs).roc.p;
                                    
                                    p_vals_ROC=[p_vals_ROC ROCoutfp1(no_ROCs).roc.p];
                                    
                                    %Enter the fp2 Ev1
                                    this_delta_dB_powerfp2Ev1=zeros(sum(trials_in_event_fp2Ev1&fp2_mask),1);
                                    this_delta_dB_powerfp2Ev1=mean(this_dB_powerfp2Ev1(:,this_band)-this_dB_powerfp2refEv1(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:sum(trials_in_event_fp2Ev1&fp2_mask),1)=this_delta_dB_powerfp2Ev1;
                                    roc_data(1:sum(trials_in_event_fp2Ev1&fp2_mask),2)=zeros(sum(trials_in_event_fp2Ev1&fp2_mask),1);
                                    
                                    %Enter fp2 Ev2
                                    total_trials=sum(trials_in_event_fp2Ev1&fp2_mask)+sum(trials_in_event_fp2Ev2&fp2_mask);
                                    this_delta_dB_powerfp2Ev2=zeros(sum(trials_in_event_fp2Ev2&fp2_mask),1);
                                    this_delta_dB_powerfp2Ev2=mean(this_dB_powerfp2Ev2(:,this_band)-this_dB_powerfp2refEv2(:,this_band),2);
                                    roc_data(sum(trials_in_event_fp2Ev1&fp2_mask)+1:total_trials,1)=this_delta_dB_powerfp2Ev2;
                                    roc_data(sum(trials_in_event_fp2Ev1&fp2_mask)+1:total_trials,2)=ones(sum(trials_in_event_fp2Ev2&fp2_mask),1);
                                    
                                    
                                    %Find fp2 ROC
                                    ROCoutfp2(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCoutfp2(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).fileNo;
                                    ROCgroupNofp2(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).fileNo);
                                    ROCoutfp2(no_ROCs).timeWindow=winNo;
                                    ROCbandwidthfp2(no_ROCs)=bwii;
                                    auROCfp2(no_ROCs)=ROCoutfp2(no_ROCs).roc.AUC-0.5;
                                    p_valROCfp2(no_ROCs)=ROCoutfp2(no_ROCs).roc.p;
                                    
                                    p_vals_ROC=[p_vals_ROC ROCoutfp2(no_ROCs).roc.p];
                                    
                                    
                                    %Are the delta dB LFP's different?
                                    
                                    %Ev1
                                    
                                    p_val(no_dBs,bwii)=ranksum(this_delta_dB_powerfp1Ev1,this_delta_dB_powerfp2Ev1);
                                    p_vals=[p_vals p_val(no_dBs,bwii)];
                                    groupNofp1(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                    groupNofp2(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                    events(no_dBs)=1;
                                    
                                    
                                    %Ev2
                                    p_val(no_dBs+1,bwii)=ranksum(this_delta_dB_powerfp1Ev2,this_delta_dB_powerfp2Ev2);
                                    p_vals=[p_vals p_val(no_dBs+1,bwii)];
                                    groupNofp1(no_dBs+1)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                    groupNofp2(no_dBs+1)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                    events(no_dBs+1)=2;
                                    
                                    if p_val(no_dBs,bwii)<0.05
                                        dB_power_changeEv1(no_ROCs)=1;
                                    else
                                        dB_power_changeEv1(no_ROCs)=0;
                                    end
                                    
                                    if p_val(no_dBs+1,bwii)<0.05
                                        dB_power_changeEv2(no_ROCs)=1;
                                    else
                                        dB_power_changeEv2(no_ROCs)=0;
                                    end
                                    
                                    %Save the data
                                    
                                    %Ev1, all points
                                    delta_dB_powerfp1Ev1(no_ROCs)=mean(this_delta_dB_powerfp1Ev1);
                                    delta_dB_powerfp2Ev1(no_ROCs)=mean(this_delta_dB_powerfp2Ev1);
                                    
                                    %Ev2, all points
                                    delta_dB_powerfp1Ev2(no_ROCs)=mean(this_delta_dB_powerfp1Ev2);
                                    delta_dB_powerfp2Ev2(no_ROCs)=mean(this_delta_dB_powerfp2Ev2);
                                    
                                    
                                    
                                    %Plot these points
                                    %Odorant 1 S+ on the left
                                    figure(2*(bwii-1)+1)
                                    pos2=[0.8 0.1 0.1 0.8];
                                    subplot('Position',pos2)
                                    hold on
                                    plot([0 1],[delta_dB_powerfp2Ev1(no_ROCs) delta_dB_powerfp1Ev2(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                    set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                                    
                                    %Odorant 2 S+ on the right
                                    figure(2*(bwii-1)+2)
                                    pos2=[0.8 0.1 0.1 0.8];
                                    subplot('Position',pos2)
                                    hold on
                                    plot([1 0],[delta_dB_powerfp1Ev1(no_ROCs) delta_dB_powerfp2Ev2(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                    set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                                end
                                
                                no_dBs=no_dBs+2;
                                
                            else
                                
                                if (sum(trials_in_event_fp1Ev1)<min_trials_per_event)
                                    fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_fp1Ev1),event1, min_trials_per_event,file_pairs(fps,1),elec);
                                end
                                
                                if (sum(trials_in_event_fp1Ev2)<min_trials_per_event)
                                    fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_fp1Ev2),event2, min_trials_per_event,file_pairs(fps,1),elec);
                                end
                                
                                if (sum(trials_in_event_fp2Ev1)<min_trials_per_event)
                                    fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_fp2Ev1),event1, min_trials_per_event,file_pairs(fps,2),elec);
                                end
                                
                                if (sum(trials_in_event_fp2Ev2)<min_trials_per_event)
                                    fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_fp2Ev2),event2, min_trials_per_event,file_pairs(fps,2),elec);
                                end
                                
                            end
                            
                        else
                            
                            if (length(handles_drgb.drgb.lfpevpair(lfpodNofp1).which_eventLFPPower(1,:))<trials_to_process)
                                fprintf(1, ['%d trials fewer than %d trials to process for file No %d electrode %d\n'],length(handles_drgb.drgb.lfpevpair(lfpodNofp1).which_eventLFPPower(1,:)),trials_to_process,file_pairs(fps,1),elec);
                            end
                            
                            if (length(handles_drgb.drgb.lfpevpair(lfpodNofp2).which_eventLFPPower(1,:))<trials_to_process)
                                fprintf(1, ['%d trials fewer than %d trials to process for file No %d electrode %d\n'],length(handles_drgb.drgb.lfpevpair(lfpodNofp2).which_eventLFPPower(1,:)),trials_to_process,file_pairs(fps,1),elec);
                            end
                            
                        end
                    else
                        
                        if isempty(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).allPower)
                            fprintf(1, ['Empty allPower for file No %d electrode %d\n'],file_pairs(fps,1),elec);
                        end
                        
                        if isempty(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).allPower)
                            fprintf(1, ['Empty allPower for file No %d electrode %d\n'],file_pairs(fps,2),elec);
                        end
                        
                    end
                    
                else
                    
                    if isempty(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref))
                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],file_pairs(fps,1),elec);
                    end
                    
                    if isempty(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref))
                        fprintf(1, ['Empty lfpevpairfor file No %d electrode %d\n'],file_pairs(fps,2),elec);
                    end
                    
                end
            end
            
        end
        fprintf(1, '\n\n')
        
        %Now plot the histograms and the average for LFP power
        num_pv_perm=0;
        for bwii=1:4
            
            
            %Plot the mean and CI for odorant 1 S+ on the left
            figure(2*(bwii-1)+1)
            pos2=[0.8 0.1 0.1 0.8];
            subplot('Position',pos2)
            
            hold on
            plot([0 1],[mean(delta_dB_powerfp2Ev1(ROCbandwidthfp2==bwii)) mean(delta_dB_powerfp1Ev2(ROCbandwidthfp1==bwii))],'-k','LineWidth', 3)
            CI = bootci(1000, @mean, delta_dB_powerfp2Ev1(ROCbandwidthfp2==bwii));
            plot([0 0],CI,'-r','LineWidth',3)
            plot(0,mean(delta_dB_powerfp2Ev1(ROCbandwidthfp2==bwii)),'or','MarkerSize', 10,'MarkerFace','r')
            CI = bootci(1000, @mean, delta_dB_powerfp1Ev2(ROCbandwidthfp1==bwii));
            plot([1 1],CI,'-b','LineWidth',3)
            plot(1,mean(delta_dB_powerfp1Ev2(ROCbandwidthfp1==bwii)),'ob','MarkerSize', 10,'MarkerFace','b')
            ylabel('delta Power (dB)')
            ylim([-10 10])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Plot the mean and CI for odorant 2 S+ on the right
            figure(2*(bwii-1)+2)
            pos2=[0.8 0.1 0.1 0.8];
            subplot('Position',pos2)
            
            hold on
            plot([1 0],[mean(delta_dB_powerfp1Ev1(ROCbandwidthfp1==bwii)) mean(delta_dB_powerfp2Ev2(ROCbandwidthfp2==bwii))],'-k','LineWidth', 3)
            CI = bootci(1000, @mean, delta_dB_powerfp1Ev1(ROCbandwidthfp1==bwii));
            plot([1 1],CI,'-r','LineWidth',3)
            plot(1,mean(delta_dB_powerfp1Ev1(ROCbandwidthfp1==bwii)),'or','MarkerSize', 10,'MarkerFace','r')
            CI = bootci(1000, @mean, delta_dB_powerfp2Ev2(ROCbandwidthfp2==bwii));
            plot([0 0],CI,'-b','LineWidth',3)
            plot(0,mean(delta_dB_powerfp2Ev2(ROCbandwidthfp2==bwii)),'ob','MarkerSize', 10,'MarkerFace','b')
            ylabel('delta Power (dB)')
            ylim([-10 10])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Plot the histograms for odorant 1 S+ red
            figure(2*(bwii-1)+1)
            edges=[-15:1:15];
            pos2=[0.1 0.1 0.6 0.8];
            subplot('Position',pos2)
            hold on
            
            h1=histogram(delta_dB_powerfp2Ev1(ROCbandwidthfp2==bwii),edges);
            h1.FaceColor='r';
            h2=histogram(delta_dB_powerfp1Ev2(ROCbandwidthfp1==bwii),edges);
            h2.FaceColor='b';
            xlabel('delta Power (dB)')
            ylabel('# of electrodes')
            legend([odorant{1} ' S+'],[odorant{1} ' S-'])
            xlim([-12 12])
            ylim([0 30])
            title(['delta LFP power (dB) for ' odorant{1} ' ' freq_names{bwii}])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Plot the histograms for odorant 2 S+ red
            figure(2*(bwii-1)+2)
            edges=[-15:1:15];
            pos2=[0.1 0.1 0.6 0.8];
            subplot('Position',pos2)
            hold on
            
            h1=histogram(delta_dB_powerfp1Ev1(ROCbandwidthfp1==bwii),edges);
            h1.FaceColor='r';
            h2=histogram(delta_dB_powerfp2Ev2(ROCbandwidthfp2==bwii),edges);
            h2.FaceColor='b';
            xlabel('delta Power (dB)')
            ylabel('# of electrodes')
            legend([odorant{2} ' S+'],[odorant{2} ' S-'])
            xlim([-12 12])
            ylim([0 30])
            title(['delta LFP power (dB) for ' odorant{2} ' ' freq_names{bwii}])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            
            %Odorant 1
            a={delta_dB_powerfp2Ev1(ROCbandwidthfp2==bwii)' delta_dB_powerfp1Ev2(ROCbandwidthfp1==bwii)'};
            mode_statcond='perm';
            num_pv_perm=num_pv_perm+1;
            [F df pvals_perm(num_pv_perm)] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for premuted anovan dB delta power S+ vs S- for odorant ' odorant{1} ' ' freq_names{bwii} '= %d\n\n'],  pvals_perm(bwii));
            
            
            %Odorant 2
            a={delta_dB_powerfp1Ev1(ROCbandwidthfp1==bwii)' delta_dB_powerfp2Ev2(ROCbandwidthfp2==bwii)'};
            mode_statcond='perm';
            num_pv_perm=num_pv_perm+1;
            [F df pvals_perm(num_pv_perm)] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for premuted anovan dB delta power S+ vs S- for odorant ' odorant{2} ' ' freq_names{bwii} '= %d\n\n'],  pvals_perm(bwii));
            
        end
        
        pFDRanovan=drsFDRpval(pvals_perm);
        fprintf(1, ['pFDR for premuted anovan p value  = %d\n\n'],pFDRanovan);
        
        
        
        fprintf(1, '\n\n')
        
        pFDRROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for significant difference of auROC p value from 0.5  = %d\n\n'],pFDRROC);
        
        
        fprintf(1, '\n\n')
        
        
        %Plot cumulative histos for auROCs
        
        figNo=8;
        p_val_ROC=[];
        
        
        x=0;
        
        for bwii=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            %Plot the histograms
            
            edges=[-0.5:0.05:0.5];
            pos2=[0.1 0.1 0.6 0.8];
            subplot('Position',pos2)
            hold on
            
            h2=histogram(auROCfp2(ROCbandwidthfp1==bwii),edges);
            h2.FaceColor='b';
            h1=histogram(auROCfp1(ROCbandwidthfp1==bwii),edges);
            h1.FaceColor='r';
            
            xlabel('auROC')
            ylabel('# of electrodes')
            legend(file_label{2},file_label{1})
            title(['auROC for ' freq_names{bwii}])
            xlim([-0.3 0.6])
            ylim([0 30])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Plot the single electrodes
            pos2=[0.8 0.1 0.1 0.8];
            subplot('Position',pos2)
            hold on
            for ii=1:length(auROCfp1)
                if ROCbandwidthfp1(ii)==bwii
                    plot([0 1],[auROCfp2(ii) auROCfp1(ii)],'-o', 'Color',[0.7 0.7 0.7])
                end
            end
            
            %PLot the mean and 95% CI
            plot([0 1],[mean(auROCfp2(ROCbandwidthfp1==bwii)) mean(auROCfp1(ROCbandwidthfp1==bwii))],'-k','LineWidth', 3)
            CI = bootci(1000, @mean, auROCfp2(ROCbandwidthfp1==bwii));
            plot([0 0],CI,'-b','LineWidth',3)
            plot(0,mean(auROCfp2(ROCbandwidthfp1==bwii)),'ob','MarkerSize', 10,'MarkerFace','b')
            CI = bootci(1000, @mean, auROCfp1(ROCbandwidthfp1==bwii));
            plot([1 1],CI,'-r','LineWidth',3)
            plot(1,mean(auROCfp1(ROCbandwidthfp1==bwii)),'or','MarkerSize', 10,'MarkerFace','r')
            ylabel('auROC')
            ylim([-0.2 0.5])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Do the statistics for auROC differences
            a={auROCfp2(ROCbandwidthfp1==bwii)' auROCfp1(ROCbandwidthfp1==bwii)'};
            mode_statcond='perm';
            [F df pval_auROCperm] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for permuted anovan for auROC S+ vs S- ' freq_names{bwii} '= %d\n\n'],  pval_auROCperm);
            pvals_auROCperm=[pvals_auROCperm pval_auROCperm];
            
            
            figure(13)
            hold on
            
            percent_auROCfp2=100*sum(p_valROCfp2(ROCbandwidthfp1==bwii)<=pFDRROC)/sum(ROCbandwidthfp1==bwii);
            bar(x,percent_auROCfp2,'b')
            
            learn_sig(bwii)=sum(p_valROCfp2(ROCbandwidthfp1==bwii)<=pFDRROC);
            learn_not_sig(bwii)=sum(ROCbandwidthfp1==bwii)-sum(p_valROCfp2(ROCbandwidthfp1==bwii)<=pFDRROC);
            
            percent_auROCfp1=100*sum(p_valROCfp1(ROCbandwidthfp1==bwii)<=pFDRROC)/sum(ROCbandwidthfp1==bwii);
            bar(x+1,percent_auROCfp1,'r')
            
            prof_sig(bwii)=sum(p_valROCfp1(ROCbandwidthfp1==bwii)<=pFDRROC);
            prof_not_sig(bwii)=sum(ROCbandwidthfp1==bwii)-sum(p_valROCfp1(ROCbandwidthfp1==bwii)<=pFDRROC);
            
            
            
            x=x+3;
            
        end
        
        figure(13)
        title('Percent significant auROC')
        legend(file_label{2},file_label{1})
        ylim([0 100])
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        
        pFDRanovanauROC=drsFDRpval(pval_auROCperm);
        fprintf(1, ['\npFDR for premuted anovan p value for difference between ' file_label{1} ' and ' file_label{2} ' for auROC = %d\n\n'],pFDRanovanauROC);
        
        
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix],'learn_sig','learn_not_sig','prof_sig','prof_not_sig');
        
        pffft=1;
        
    case 3
        %Generate Fig. 2  for Daniels' LFP power paper. For the proficient mice in the first and last sessions
        %plot the LFP spectrum for S+ vs S-, plot LFP power for S+ vs S- for each electrode and plot auROCs
        %NOTE: This does the analysis in all the files and DOES not distinguish between groups!!!
        no_dBs=1;
        delta_dB_power=[];
        no_ROCs=0;
        ROCout=[];
        p_vals_ROC=[];
        delta_dB_powerEv1=[];
        no_Ev1=0;
        noWB=0;
        delta_dB_powerEv1WB=[];
        delta_dB_powerEv2WB=[];
        
        fprintf(1, ['Pairwise auROC analysis for Fig 1 of Daniel''s paper\n\n'])
        p_vals=[];
        no_files=length(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        for fileNo=1:no_files
            for elec=1:16
                if sum(which_electrodes==elec)>0
                    lfpodNo_ref=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref)))
                        
                        
                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower))
                            
                            if (length(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(1,:))>=trials_to_process)
                                
                                trials=length(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(1,:));
                                mask=logical([zeros(1,trials-trials_to_process) ones(1,trials_to_process)]);
                                trials_in_eventEv1=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(event1,:)==1);
                                trials_in_eventEv2=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(event2,:)==1);
                                
                                if (sum(trials_in_eventEv1)>=min_trials_per_event) & (sum(trials_in_eventEv2)>=min_trials_per_event)
                                    
                                    lfpodNo=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                    
                                    % Ev1
                                    this_dB_powerrefEv1=zeros(sum(trials_in_eventEv1&mask),length(frequency));
                                    this_dB_powerrefEv1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower(trials_in_eventEv1&mask,:));
                                    
                                    
                                    this_dB_powerEv1=zeros(sum(trials_in_eventEv1&mask),length(frequency));
                                    this_dB_powerEv1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_eventEv1&mask,:));
                                    
                                    % Ev2
                                    this_dB_powerrefEv2=zeros(sum(trials_in_eventEv2&mask),length(frequency));
                                    this_dB_powerrefEv2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower(trials_in_eventEv2&mask,:));
                                    
                                    
                                    this_dB_powerEv2=zeros(sum(trials_in_eventEv2&mask),length(frequency));
                                    this_dB_powerEv2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_eventEv2&mask,:));
                                    
                                    
                                    %Wide band spectrum
                                    noWB=noWB+1;
                                    
                                    delta_dB_powerEv1WB(noWB,:)=mean(this_dB_powerEv1-this_dB_powerrefEv1,1);
                                    delta_dB_powerEv2WB(noWB,:)=mean(this_dB_powerEv2-this_dB_powerrefEv2,1);
                                    
                                    
                                    %Do per badwidth analysis
                                    for bwii=1:no_bandwidths
                                        
                                        no_ROCs=no_ROCs+1;
                                        this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                        
                                        %Enter the  Ev1
                                        this_delta_dB_powerEv1=zeros(sum(trials_in_eventEv1&mask),1);
                                        this_delta_dB_powerEv1=mean(this_dB_powerEv1(:,this_band)-this_dB_powerrefEv1(:,this_band),2);
                                        roc_data=[];
                                        roc_data(1:sum(trials_in_eventEv1&mask),1)=this_delta_dB_powerEv1;
                                        roc_data(1:sum(trials_in_eventEv1&mask),2)=zeros(sum(trials_in_eventEv1&mask),1);
                                        
                                        %Enter  Ev2
                                        total_trials=sum(trials_in_eventEv1&mask)+sum(trials_in_eventEv2&mask);
                                        this_delta_dB_powerEv2=zeros(sum(trials_in_eventEv2&mask),1);
                                        this_delta_dB_powerEv2=mean(this_dB_powerEv2(:,this_band)-this_dB_powerrefEv2(:,this_band),2);
                                        roc_data(sum(trials_in_eventEv1&mask)+1:total_trials,1)=this_delta_dB_powerEv2;
                                        roc_data(sum(trials_in_eventEv1&mask)+1:total_trials,2)=ones(sum(trials_in_eventEv2&mask),1);
                                        
                                        
                                        %Find  ROC
                                        ROCout(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                        ROCout(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo_ref).fileNo;
                                        ROCgroupNo(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNo_ref).fileNo);
                                        ROCout(no_ROCs).timeWindow=winNo;
                                        ROCbandwidth(no_ROCs)=bwii;
                                        auROC(no_ROCs)=ROCout(no_ROCs).roc.AUC-0.5;
                                        p_valROC(no_ROCs)=ROCout(no_ROCs).roc.p;
                                        
                                        p_vals_ROC=[p_vals_ROC ROCout(no_ROCs).roc.p];
                                        
                                        
                                        delta_dB_powerEv1(no_ROCs)=mean(this_delta_dB_powerEv1);
                                        delta_dB_powerEv2(no_ROCs)=mean(this_delta_dB_powerEv2);
                                        
                                        
                                        %Plot this point
                                        figure(bwii+1)
                                        pos2=[0.8 0.1 0.1 0.8];
                                        subplot('Position',pos2)
                                        hold on
                                        plot([1 0],[delta_dB_powerEv1(no_ROCs) delta_dB_powerEv2(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                                        
                                        
                                    end
                                    
                                    
                                    
                                else
                                    
                                    if (sum(trials_in_eventEv1)<min_trials_per_event)
                                        fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_eventEv1),event1, min_trials_per_event,files(fileNo),elec);
                                    end
                                    
                                    if (sum(trials_in_eventEv2)<min_trials_per_event)
                                        fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_eventEv2),event2, min_trials_per_event,files(fileNo),elec);
                                    end
                                    
                                end
                                
                            else
                                
                                fprintf(1, ['%d trials fewer than %d trials to process for file No %d electrode %d\n'],length(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(1,:)),trials_to_process,files(fileNo),elec);
                                
                            end
                        else
                            
                            fprintf(1, ['Empty allPower for file No %d electrode %d\n'],files(fileNo),elec);
                            
                        end
                        
                        
                    else
                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],files(fileNo),elec);
                        
                        
                    end
                end
            end
            
        end
        fprintf(1, '\n\n')
        
        
        %Now plot the bounded line for
        
        %Calculate the mean and 95% CI for Ev1
        dB_Ev1_ci=zeros(length(frequency),2);
        for ifreq=1:length(frequency)
            %             pd=fitdist(delta_dB_powerEv1WB(:,ifreq),'Normal');
            %             ci=paramci(pd);
            %             dB_Ev1_ci(ifreq)=pd.mu-ci(1,1);
            dB_Ev1_mean(ifreq)=mean(delta_dB_powerEv1WB(:,ifreq));
            CI = bootci(1000, @mean, delta_dB_powerEv1WB(:,ifreq));
            dB_Ev1_ci(ifreq,1)=CI(2)-dB_Ev1_mean(ifreq);
            dB_Ev1_ci(ifreq,2)=-(CI(1)-dB_Ev1_mean(ifreq));
        end
        
        figure(1)
        [hl1, hp1] = boundedline(frequency,dB_Ev1_mean, dB_Ev1_ci, 'r');
        
        %Calculate the mean and 95% CI for Ev2
        dB_Ev2_ci=zeros(length(frequency),2);
        for ifreq=1:length(frequency)
            dB_Ev2_mean(ifreq)=mean(delta_dB_powerEv2WB(:,ifreq));
            CI = bootci(1000, @mean, delta_dB_powerEv2WB(:,ifreq));
            dB_Ev2_ci(ifreq,1)=CI(2)-dB_Ev2_mean(ifreq);
            dB_Ev2_ci(ifreq,2)=-(CI(1)-dB_Ev2_mean(ifreq));
        end
        
        hold on
        [hl2, hp2] = boundedline(frequency,dB_Ev2_mean, dB_Ev2_ci, 'b');
        xlabel('Frequency (Hz)')
        ylabel('delta Power (dB)')
        legend([hl1 hl2],'S+','S-')
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        
        %Now plot the histograms and the average
        for bwii=1:4
            %Plot the average
            figure(bwii+1)
            pos2=[0.8 0.1 0.1 0.8];
            subplot('Position',pos2)
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            plot([1 0],[mean(delta_dB_powerEv1(ROCbandwidth==bwii)) mean(delta_dB_powerEv2(ROCbandwidth==bwii))],'-k','LineWidth', 3)
            CI = bootci(1000, @mean, delta_dB_powerEv1(ROCbandwidth==bwii));
            plot([1 1],CI,'-r','LineWidth',3)
            plot(1,mean(delta_dB_powerEv1(ROCbandwidth==bwii)),'or','MarkerSize', 10,'MarkerFace','r')
            CI = bootci(1000, @mean, delta_dB_powerEv2(ROCbandwidth==bwii));
            plot([0 0],CI,'-b','LineWidth',3)
            plot(0,mean(delta_dB_powerEv2(ROCbandwidth==bwii)),'ob','MarkerSize', 10,'MarkerFace','b')
            ylabel('delta Power (dB)')
            ylim([-10 15])
            
            %Plot the histograms
            
            maxdB=max([max(delta_dB_powerEv1(ROCbandwidth==bwii)) max(delta_dB_powerEv2(ROCbandwidth==bwii))]);
            mindB=min([min(delta_dB_powerEv1(ROCbandwidth==bwii)) min(delta_dB_powerEv2(ROCbandwidth==bwii))]);
            edges=[-15:1:15];
            pos2=[0.1 0.1 0.6 0.8];
            subplot('Position',pos2)
            hold on
            
            h1=histogram(delta_dB_powerEv2(ROCbandwidth==bwii),edges);
            h1.FaceColor='b';
            h2=histogram(delta_dB_powerEv1(ROCbandwidth==bwii),edges);
            h2.FaceColor='r';
            xlabel('delta Power (dB)')
            ylabel('# of electrodes')
            legend('S-','S+')
            xlim([-12 12])
            ylim([0 70])
            title(freq_names{bwii})
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            
            
            a={ delta_dB_powerEv1(ROCbandwidth==bwii)' delta_dB_powerEv2(ROCbandwidth==bwii)'};
            mode_statcond='perm';
            [F df pvals_perm(bwii)] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for premuted anovan dB delta power S+ vs S- ' freq_names{bwii} '= %d\n'],  pvals_perm(bwii));
            
        end
        
        pFDRanovan=drsFDRpval(pvals_perm);
        fprintf(1, ['pFDR for premuted anovan p value  = %d\n\n'],pFDRanovan);
        
        
        
        fprintf(1, '\n\n')
        
        
        pFDRauROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for auROC  = %d\n\n'],pFDRauROC);
        %Plot cumulative histos for auROCs
        
        figNo=5;
        p_val_ROC=[];
        edges=-0.5:0.05:0.5;
        
        for bwii=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            n_cum=0;
            this_legend=[];
            
            histogram(auROC(( p_valROC>pFDRauROC)&(ROCbandwidth==bwii)),edges)
            histogram(auROC(( p_valROC<=pFDRauROC)&(ROCbandwidth==bwii)),edges)
            legend('auROC not singificant','auROC significant')
            title(['Histogram for ' freq_names{bwii} ' auROC for LFPs'])
            xlim([-0.2 0.6])
            ylim([0 30])
        end
        
         
        
        %Plot percent significant ROC
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        figure(figNo)
        
        hold on
        
        for bwii=1:4
            bar(bwii,100*sum(( p_valROC<=pFDRauROC)&(ROCbandwidth==bwii))/sum((ROCbandwidth==bwii)))
            auROC_sig.sig(bwii)=sum(( p_valROC<=pFDRauROC)&(ROCbandwidth==bwii));
            auROC_sig.not_sig(bwii)=sum((ROCbandwidth==bwii))-sum(( p_valROC<=pFDRauROC)&(ROCbandwidth==bwii));
        end
        title('Percent auROC significantly different from zero')
        ylim([0 100])
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        pffft=1;
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix],'auROC_sig');
        
        
    case 4
        %Compare auROC in the last few trials of the last session file with
        %first few trials of the first session
        %Generates Fig. 3 for Daniel's paper. first vs last.
        no_dBs=1;
        delta_dB_power_fp1=[];
        no_ROCs=0;
        ROCoutfp1=[];
        ROCoutfp2=[];
        p_vals_ROC=[];
        delta_dB_powerfp1Ev1=[];
        no_Ev1=0;
        pvals_auROCperm=[];
        pvals_dBperm=[];
        perCorr_fp1=[];
        perCorr_fp2=[];
        
        fprintf(1, ['Pairwise auROC analysis for ' evTypeLabels{1} ' and ' evTypeLabels{2} ' LFP power\n\n'])
        p_vals=[];
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        sz_fp=size(file_pairs);
        no_file_pairs=sz_fp(1);
        
        for fps=1:no_file_pairs
            
            
            for elec=1:16
                if sum(which_electrodes==elec)>0
                    lfpodNofp1_ref=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                    lfpodNofp2_ref=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                    
                    if elec==1
                        %Find percent correct for fp1 block
                        perCorr_fp1(fps)=handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).perCorrLFPPower(end);
                        
                        %Find percent correct for fp2 block
                        perCorr_fp2(fps)=handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).perCorrLFPPower(1);
                        
                        fprintf(1, '\nPercent correct for session pair %d last= %d, first= %d\n',fps,perCorr_fp1(fps),perCorr_fp2(fps));
                        
                    end
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref)))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref)))
                        
                        
                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).allPower))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).allPower))
                            
                            if (length(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).which_eventLFPPower(1,:))>=trials_to_process) &...
                                    (length(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).which_eventLFPPower(1,:))>=trials_to_process)
                                
                                length_fp1=length(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).which_eventLFPPower(1,:));
                                if front_mask(1)==1
                                    fp1_mask=logical([ones(1,trials_to_process) zeros(1,length_fp1-trials_to_process)]);
                                else
                                    fp1_mask=logical([zeros(1,length_fp1-trials_to_process) ones(1,trials_to_process)]);
                                end
                                trials_in_event_fp1Ev1=(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).which_eventLFPPower(event1,:)==1);
                                trials_in_event_fp1Ev2=(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).which_eventLFPPower(event2,:)==1);
                                
                                length_fp2=length(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).which_eventLFPPower(1,:));
                                if front_mask(2)==1
                                    fp2_mask=logical([ones(1,trials_to_process) zeros(1,length_fp2-trials_to_process)]);
                                else
                                    fp2_mask=logical([zeros(1,length_fp2-trials_to_process) ones(1,trials_to_process)]);
                                end
                                trials_in_event_fp2Ev1=(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).which_eventLFPPower(event1,:)==1);
                                trials_in_event_fp2Ev2=(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).which_eventLFPPower(event2,:)==1);
                                
                                
                                if (sum(trials_in_event_fp1Ev1)>=min_trials_per_event) & (sum( trials_in_event_fp1Ev2)>=min_trials_per_event) & ...
                                        (sum(trials_in_event_fp2Ev1)>=min_trials_per_event) & (sum(trials_in_event_fp2Ev2)>=min_trials_per_event)
                                    
                                    lfpodNofp1=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                    lfpodNofp2=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                    
                                    %fp1 Ev1
                                    this_dB_powerfp1refEv1=zeros(sum(trials_in_event_fp1Ev1&fp1_mask),length(frequency));
                                    this_dB_powerfp1refEv1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).allPower(trials_in_event_fp1Ev1&fp1_mask,:));
                                    
                                    this_dB_powerfp1Ev1=zeros(sum(trials_in_event_fp1Ev1&fp1_mask),length(frequency));
                                    this_dB_powerfp1Ev1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp1).allPower(trials_in_event_fp1Ev1&fp1_mask,:));
                                    
                                    %fp1 Ev2
                                    this_dB_powerfp1refEv2=zeros(sum(trials_in_event_fp1Ev2&fp1_mask),length(frequency));
                                    this_dB_powerfp1refEv2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).allPower(trials_in_event_fp1Ev2&fp1_mask,:));
                                    
                                    this_dB_powerfp1Ev2=zeros(sum(trials_in_event_fp1Ev2&fp1_mask),length(frequency));
                                    this_dB_powerfp1Ev2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp1).allPower(trials_in_event_fp1Ev2&fp1_mask,:));
                                    
                                    
                                    %fp2 Ev1
                                    this_dB_powerfp2refEv1=zeros(sum(trials_in_event_fp2Ev1&fp2_mask),length(frequency));
                                    this_dB_powerfp2refEv1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).allPower(trials_in_event_fp2Ev1&fp2_mask,:));
                                    
                                    this_dB_powerfp2Ev1=zeros(sum(trials_in_event_fp2Ev1&fp2_mask),length(frequency));
                                    this_dB_powerfp2Ev1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp2).allPower(trials_in_event_fp2Ev1&fp2_mask,:));
                                    
                                    %fp2 Ev2
                                    this_dB_powerfp2refEv2=zeros(sum(trials_in_event_fp2Ev2&fp2_mask),length(frequency));
                                    this_dB_powerfp2refEv2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).allPower(trials_in_event_fp2Ev2&fp2_mask,:));
                                    
                                    this_dB_powerfp2Ev2=zeros(sum(trials_in_event_fp2Ev2&fp2_mask),length(frequency));
                                    this_dB_powerfp2Ev2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNofp2).allPower(trials_in_event_fp2Ev2&fp2_mask,:));
                                    
                                    for bwii=1:no_bandwidths
                                        
                                        no_ROCs=no_ROCs+1;
                                        this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                        
                                        %Enter the fp1 Ev1
                                        this_delta_dB_powerfp1Ev1=zeros(sum(trials_in_event_fp1Ev1&fp1_mask),1);
                                        this_delta_dB_powerfp1Ev1=mean(this_dB_powerfp1Ev1(:,this_band)-this_dB_powerfp1refEv1(:,this_band),2);
                                        roc_data=[];
                                        roc_data(1:sum(trials_in_event_fp1Ev1&fp1_mask),1)=this_delta_dB_powerfp1Ev1;
                                        roc_data(1:sum(trials_in_event_fp1Ev1&fp1_mask),2)=zeros(sum(trials_in_event_fp1Ev1&fp1_mask),1);
                                        
                                        %Enter fp1 Ev2
                                        total_trials=sum(trials_in_event_fp1Ev1&fp1_mask)+sum(trials_in_event_fp1Ev2&fp1_mask);
                                        this_delta_dB_powerfp1Ev2=zeros(sum(trials_in_event_fp1Ev2&fp1_mask),1);
                                        this_delta_dB_powerfp1Ev2=mean(this_dB_powerfp1Ev2(:,this_band)-this_dB_powerfp1refEv2(:,this_band),2);
                                        roc_data(sum(trials_in_event_fp1Ev1&fp1_mask)+1:total_trials,1)=this_delta_dB_powerfp1Ev2;
                                        roc_data(sum(trials_in_event_fp1Ev1&fp1_mask)+1:total_trials,2)=ones(sum(trials_in_event_fp1Ev2&fp1_mask),1);
                                        
                                        
                                        %Find fp1 ROC
                                        ROCoutfp1(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                        ROCoutfp1(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).fileNo;
                                        ROCgroupNofp1(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).fileNo);
                                        ROCoutfp1(no_ROCs).timeWindow=winNo;
                                        ROCbandwidthfp1(no_ROCs)=bwii;
                                        auROCfp1(no_ROCs)=ROCoutfp1(no_ROCs).roc.AUC-0.5;
                                        p_valROCfp1(no_ROCs)=ROCoutfp1(no_ROCs).roc.p;
                                        
                                        p_vals_ROC=[p_vals_ROC ROCoutfp1(no_ROCs).roc.p];
                                        
                                        %Enter the fp2 Ev1
                                        this_delta_dB_powerfp2Ev1=zeros(sum(trials_in_event_fp2Ev1&fp2_mask),1);
                                        this_delta_dB_powerfp2Ev1=mean(this_dB_powerfp2Ev1(:,this_band)-this_dB_powerfp2refEv1(:,this_band),2);
                                        roc_data=[];
                                        roc_data(1:sum(trials_in_event_fp2Ev1&fp2_mask),1)=this_delta_dB_powerfp2Ev1;
                                        roc_data(1:sum(trials_in_event_fp2Ev1&fp2_mask),2)=zeros(sum(trials_in_event_fp2Ev1&fp2_mask),1);
                                        
                                        %Enter fp2 Ev2
                                        total_trials=sum(trials_in_event_fp2Ev1&fp2_mask)+sum(trials_in_event_fp2Ev2&fp2_mask);
                                        this_delta_dB_powerfp2Ev2=zeros(sum(trials_in_event_fp2Ev2&fp2_mask),1);
                                        this_delta_dB_powerfp2Ev2=mean(this_dB_powerfp2Ev2(:,this_band)-this_dB_powerfp2refEv2(:,this_band),2);
                                        roc_data(sum(trials_in_event_fp2Ev1&fp2_mask)+1:total_trials,1)=this_delta_dB_powerfp2Ev2;
                                        roc_data(sum(trials_in_event_fp2Ev1&fp2_mask)+1:total_trials,2)=ones(sum(trials_in_event_fp2Ev2&fp2_mask),1);
                                        
                                        
                                        %Find fp2 ROC
                                        ROCoutfp2(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                        ROCoutfp2(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).fileNo;
                                        ROCgroupNofp2(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).fileNo);
                                        ROCoutfp2(no_ROCs).timeWindow=winNo;
                                        ROCbandwidthfp2(no_ROCs)=bwii;
                                        auROCfp2(no_ROCs)=ROCoutfp2(no_ROCs).roc.AUC-0.5;
                                        p_valROCfp2(no_ROCs)=ROCoutfp2(no_ROCs).roc.p;
                                        
                                        p_vals_ROC=[p_vals_ROC ROCoutfp2(no_ROCs).roc.p];
                                        
                                        
                                        %Are the delta dB LFP's different?
                                        
                                        %Ev1
                                        
                                        p_val(no_dBs,bwii)=ranksum(this_delta_dB_powerfp1Ev1,this_delta_dB_powerfp2Ev1);
                                        p_vals=[p_vals p_val(no_dBs,bwii)];
                                        groupNofp1(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                        groupNofp2(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                        events(no_dBs)=1;
                                        
                                        
                                        %Ev2
                                        p_val(no_dBs+1,bwii)=ranksum(this_delta_dB_powerfp1Ev2,this_delta_dB_powerfp2Ev2);
                                        p_vals=[p_vals p_val(no_dBs+1,bwii)];
                                        groupNofp1(no_dBs+1)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                        groupNofp2(no_dBs+1)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                        events(no_dBs+1)=2;
                                        
                                        if p_val(no_dBs,bwii)<0.05
                                            dB_power_changeEv1(no_ROCs)=1;
                                        else
                                            dB_power_changeEv1(no_ROCs)=0;
                                        end
                                        
                                        if p_val(no_dBs+1,bwii)<0.05
                                            dB_power_changeEv2(no_ROCs)=1;
                                        else
                                            dB_power_changeEv2(no_ROCs)=0;
                                        end
                                        
                                        
                                        
                                        %Ev1, all points
                                        delta_dB_powerfp1Ev1(no_ROCs)=mean(this_delta_dB_powerfp1Ev1);
                                        delta_dB_powerfp2Ev1(no_ROCs)=mean(this_delta_dB_powerfp2Ev1);
                                        
                                        %Ev2, all points
                                        delta_dB_powerfp1Ev2(no_ROCs)=mean(this_delta_dB_powerfp1Ev2);
                                        delta_dB_powerfp2Ev2(no_ROCs)=mean(this_delta_dB_powerfp2Ev2);
                                        
                                        
                                    end
                                    
                                    no_dBs=no_dBs+2;
                                    
                                else
                                    
                                    if (sum(trials_in_event_fp1Ev1)<min_trials_per_event)
                                        fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_fp1Ev1),event1, min_trials_per_event,file_pairs(fps,1),elec);
                                    end
                                    
                                    if (sum(trials_in_event_fp1Ev2)<min_trials_per_event)
                                        fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_fp1Ev2),event2, min_trials_per_event,file_pairs(fps,1),elec);
                                    end
                                    
                                    if (sum(trials_in_event_fp2Ev1)<min_trials_per_event)
                                        fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_fp2Ev1),event1, min_trials_per_event,file_pairs(fps,2),elec);
                                    end
                                    
                                    if (sum(trials_in_event_fp2Ev2)<min_trials_per_event)
                                        fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_fp2Ev2),event2, min_trials_per_event,file_pairs(fps,2),elec);
                                    end
                                    
                                end
                                
                            else
                                
                                if (length(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).which_eventLFPPower(1,:))<trials_to_process)
                                    fprintf(1, ['%d trials fewer than %d trials to process for file No %d electrode %d\n'],length(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).which_eventLFPPower(1,:)),trials_to_process,file_pairs(fps,1),elec);
                                end
                                
                                if (length(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).which_eventLFPPower(1,:))<trials_to_process)
                                    fprintf(1, ['%d trials fewer than %d trials to process for file No %d electrode %d\n'],length(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).which_eventLFPPower(1,:)),trials_to_process,file_pairs(fps,1),elec);
                                end
                                
                            end
                        else
                            
                            if isempty(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref).allPower)
                                fprintf(1, ['Empty allPower for file No %d electrode %d\n'],file_pairs(fps,1),elec);
                            end
                            
                            if isempty(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref).allPower)
                                fprintf(1, ['Empty allPower for file No %d electrode %d\n'],file_pairs(fps,2),elec);
                            end
                            
                        end
                        
                    else
                        
                        if isempty(handles_drgb.drgb.lfpevpair(lfpodNofp1_ref))
                            fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],file_pairs(fps,1),elec);
                        end
                        
                        if isempty(handles_drgb.drgb.lfpevpair(lfpodNofp2_ref))
                            fprintf(1, ['Empty lfpevpairfor file No %d electrode %d\n'],file_pairs(fps,2),elec);
                        end
                        
                    end
                end
            end
            
        end
        fprintf(1, '\n\n')
        
        pFDRROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for significant difference of auROC p value from 0.5  = %d\n\n'],pFDRROC);
        
        
        
        fprintf(1, '\n\n')
        
        
        %Plot cumulative histos for auROCs
        dB_power_change=logical(dB_power_changeEv1+dB_power_changeEv2);
        figNo=0;
        p_val_ROC=[];
        
        try
            close(5)
        catch
        end
        figure(5)
        hold on
        x=0;
        
        for bwii=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            %Plot the histograms
            edges=[-0.5:0.05:0.5];
            pos2=[0.1 0.1 0.6 0.8];
            subplot('Position',pos2)
            hold on
            
            h2=histogram(auROCfp2(ROCbandwidthfp1==bwii),edges);
            h2.FaceColor='b';
            h1=histogram(auROCfp1(ROCbandwidthfp1==bwii),edges);
            h1.FaceColor='r';
            
            xlabel('auROC')
            ylabel('# of electrodes')
            legend(file_label{2},file_label{1})
            title(['auROC for ' freq_names{bwii}])
            xlim([-0.3 0.6])
            ylim([0 30])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Plot the single electrodes
            pos2=[0.8 0.1 0.1 0.8];
            subplot('Position',pos2)
            hold on
            for ii=1:length(auROCfp1)
                if ROCbandwidthfp1(ii)==bwii
                    plot([0 1],[auROCfp2(ii) auROCfp1(ii)],'-o', 'Color',[0.7 0.7 0.7])
                end
            end
            
            %PLot the mean and 95% CI
            plot([0 1],[mean(auROCfp2(ROCbandwidthfp1==bwii)) mean(auROCfp1(ROCbandwidthfp1==bwii))],'-k','LineWidth', 3)
            CI = bootci(1000, @mean, auROCfp2(ROCbandwidthfp1==bwii));
            plot([0 0],CI,'-b','LineWidth',3)
            plot(0,mean(auROCfp2(ROCbandwidthfp1==bwii)),'ob','MarkerSize', 10,'MarkerFace','b')
            CI = bootci(1000, @mean, auROCfp1(ROCbandwidthfp1==bwii));
            plot([1 1],CI,'-r','LineWidth',3)
            plot(1,mean(auROCfp1(ROCbandwidthfp1==bwii)),'or','MarkerSize', 10,'MarkerFace','r')
            ylabel('auROC')
            ylim([-0.2 0.5])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Do the statistics for auROC differences
            a={auROCfp2(ROCbandwidthfp1==bwii)' auROCfp1(ROCbandwidthfp1==bwii)'};
            mode_statcond='perm';
            [F df pval_auROCperm] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for permuted anovan for auROC S+ vs S- ' freq_names{bwii} '= %d\n\n'],  pval_auROCperm);
            pvals_auROCperm=[pvals_auROCperm pval_auROCperm];
            
            %Figure 5
            figure(5)
            
            percent_auROCfp2=100*sum(p_valROCfp2(ROCbandwidthfp1==bwii)<=pFDRROC)/sum(ROCbandwidthfp1==bwii);
            bar(x,percent_auROCfp2,'b')
            
            learn_sig(bwii)=sum(p_valROCfp2(ROCbandwidthfp1==bwii)<=pFDRROC);
            learn_not_sig(bwii)=sum(ROCbandwidthfp1==bwii)-sum(p_valROCfp2(ROCbandwidthfp1==bwii)<=pFDRROC);
            
            percent_auROCfp1=100*sum(p_valROCfp1(ROCbandwidthfp1==bwii)<=pFDRROC)/sum(ROCbandwidthfp1==bwii);
            bar(x+1,percent_auROCfp1,'r')
            
            prof_sig(bwii)=sum(p_valROCfp1(ROCbandwidthfp1==bwii)<=pFDRROC);
            prof_not_sig(bwii)=sum(ROCbandwidthfp1==bwii)-sum(p_valROCfp1(ROCbandwidthfp1==bwii)<=pFDRROC);
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            x=x+3;
            
        end
        
        figure(5)
        title('Percent singificant auROC')
        legend(file_label{2},file_label{1})
        ylim([0 100])
        
        pFDRanovanauROC=drsFDRpval(pval_auROCperm);
        fprintf(1, ['\npFDR for premuted anovan p value for difference between ' file_label{1} ' and ' file_label{2} ' for auROC = %d\n\n'],pFDRanovanauROC);
        
        
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix],'learn_sig','learn_not_sig','prof_sig','prof_not_sig');
        
        pffft=1;
        
    case 5
        %Compare auROC in the last few trials of pre with first few trials of post
        %Used for Fig. 5 of Daniel's paper
        no_dBs=1;
        delta_dB_power_pre=[];
        no_ROCs=0;
        ROCoutpre=[];
        ROCoutpost=[];
        p_vals_ROC=[];
        delta_dB_powerpreHit=[];
        no_hits=0;
        perCorr_pre=[];
        perCorr_post=[];
        group_pre=[];
        group_post=[];
        
        fprintf(1, ['Pairwise auROC analysis for ', evTypeLabels{1} ' and ' evTypeLabels{2} ' LFP power\n\n'])
        p_vals=[];
        for fps=1:no_file_pairs
            for elec=1:16
                
                lfpodNopre_ref=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                lfpodNopost_ref=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                
             
                
                if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref)))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref)))
                    
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower))
                        
                        if (length(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).which_eventLFPPower(1,:))>=trials_to_process) &...
                                (length(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).which_eventLFPPower(1,:))>=trials_to_process)
                            
                            length_pre=length(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).which_eventLFPPower(1,:));
                            pre_mask=logical([zeros(1,length_pre-trials_to_process) ones(1,trials_to_process)]);
                            trials_in_event_preHit=(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).which_eventLFPPower(event1,:)==1);
                            trials_in_event_preCR=(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).which_eventLFPPower(event2,:)==1);
                            
                            length_post=length(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).which_eventLFPPower(1,:));
                            post_mask=logical([ones(1,trials_to_process) zeros(1,length_post-trials_to_process)]);
                            trials_in_event_postHit=(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).which_eventLFPPower(event1,:)==1);
                            trials_in_event_postCR=(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).which_eventLFPPower(event2,:)==1);
                            
                            if (sum(trials_in_event_preHit)>=min_trials_per_event) & (sum(trials_in_event_preCR)>=min_trials_per_event) & ...
                                    (sum(trials_in_event_postHit)>=min_trials_per_event) & (sum(trials_in_event_postCR)>=min_trials_per_event)
                                
                                
                                lfpodNopre=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                lfpodNopost=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                
                                %pre Hits
                                this_dB_powerprerefHit=zeros(sum(trials_in_event_preHit&pre_mask),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerprerefHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower(trials_in_event_preHit&pre_mask,:));
                                
                                this_dB_powerpreHit=zeros(sum(trials_in_event_preHit&pre_mask),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpreHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre).allPower(trials_in_event_preHit&pre_mask,:));
                                
                                
                                %pre CRs
                                this_dB_powerprerefCR=zeros(sum(trials_in_event_preCR&pre_mask),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerprerefCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower(trials_in_event_preCR&pre_mask,:));
                                
                                this_dB_powerpreCR=zeros(sum(trials_in_event_preCR&pre_mask),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpreCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre).allPower(trials_in_event_preCR&pre_mask,:));
                                
                                %post Hits
                                this_dB_powerpostrefHit=zeros(sum(trials_in_event_postHit&post_mask),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostrefHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower(trials_in_event_postHit&post_mask,:));
                                
                                this_dB_powerpostHit=zeros(sum(trials_in_event_postHit&post_mask),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost).allPower(trials_in_event_postHit&post_mask,:));
                                
                                
                                %post CRs
                                this_dB_powerpostrefCR=zeros(sum(trials_in_event_postCR&post_mask),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostrefCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower(trials_in_event_postCR&post_mask,:));
                                
                                this_dB_powerpostCR=zeros(sum(trials_in_event_postCR&post_mask),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost).allPower(trials_in_event_postCR&post_mask,:));
                                
                                for bwii=1:no_bandwidths
                                    
                                    no_ROCs=no_ROCs+1;
                                    this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                    
                                    %Enter the pre Hits
                                    this_delta_dB_powerpreHit=zeros(sum(trials_in_event_preHit&pre_mask),1);
                                    this_delta_dB_powerpreHit=mean(this_dB_powerpreHit(:,this_band)-this_dB_powerprerefHit(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:sum(trials_in_event_preHit&pre_mask),1)=this_delta_dB_powerpreHit;
                                    roc_data(1:sum(trials_in_event_preHit&pre_mask),2)=zeros(sum(trials_in_event_preHit&pre_mask),1);
                                    
                                    %Enter pre CR
                                    total_trials=sum(trials_in_event_preHit&pre_mask)+sum(trials_in_event_preCR&pre_mask);
                                    this_delta_dB_powerpreCR=zeros(sum(trials_in_event_preCR&pre_mask),1);
                                    this_delta_dB_powerpreCR=mean(this_dB_powerpreCR(:,this_band)-this_dB_powerprerefCR(:,this_band),2);
                                    roc_data(sum(trials_in_event_preHit&pre_mask)+1:total_trials,1)=this_delta_dB_powerpreCR;
                                    roc_data(sum(trials_in_event_preHit&pre_mask)+1:total_trials,2)=ones(sum(trials_in_event_preCR&pre_mask),1);
                                    
                                    
                                    %Find pre ROC
                                    ROCoutpre(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCoutpre(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopre_ref).fileNo;
                                    ROCgroupNopre(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).fileNo);
                                    ROCoutpre(no_ROCs).timeWindow=winNo;
                                    ROCbandwidthpre(no_ROCs)=bwii;
                                    auROCpre(no_ROCs)=ROCoutpre(no_ROCs).roc.AUC-0.5;
                                    p_valROCpre(no_ROCs)=ROCoutpre(no_ROCs).roc.p;
                                    
                                    p_vals_ROC=[p_vals_ROC ROCoutpre(no_ROCs).roc.p];
                                    
                                    %Enter the post Hits
                                    this_delta_dB_powerpostHit=zeros(sum(trials_in_event_postHit&post_mask),1);
                                    this_delta_dB_powerpostHit=mean(this_dB_powerpostHit(:,this_band)-this_dB_powerpostrefHit(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:sum(trials_in_event_postHit&post_mask),1)=this_delta_dB_powerpostHit;
                                    roc_data(1:sum(trials_in_event_postHit&post_mask),2)=zeros(sum(trials_in_event_postHit&post_mask),1);
                                    
                                    %Enter post CR
                                    total_trials=sum(trials_in_event_postHit&post_mask)+sum(trials_in_event_postCR&post_mask);
                                    this_delta_dB_powerpostCR=zeros(sum(trials_in_event_postCR&post_mask),1);
                                    this_delta_dB_powerpostCR=mean(this_dB_powerpostCR(:,this_band)-this_dB_powerpostrefCR(:,this_band),2);
                                    roc_data(sum(trials_in_event_postHit&post_mask)+1:total_trials,1)=this_delta_dB_powerpostCR;
                                    roc_data(sum(trials_in_event_postHit&post_mask)+1:total_trials,2)=ones(sum(trials_in_event_postCR&post_mask),1);
                                    
                                    
                                    %Find post ROC
                                    ROCoutpost(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCoutpost(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopost_ref).fileNo;
                                    ROCgroupNopost(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).fileNo);
                                    ROCoutpost(no_ROCs).timeWindow=winNo;
                                    ROCbandwidthpost(no_ROCs)=bwii;
                                    auROCpost(no_ROCs)=ROCoutpost(no_ROCs).roc.AUC-0.5;
                                    p_valROCpost(no_ROCs)=ROCoutpost(no_ROCs).roc.p;
                                    
                                    p_vals_ROC=[p_vals_ROC ROCoutpost(no_ROCs).roc.p];
                                    
                                    if (auROCpost(no_ROCs)<0.3)&(auROCpre(no_ROCs)>0.4)&(ROCgroupNopre(no_ROCs)==1)&(ROCbandwidthpre(no_ROCs)==2)
                                        fprintf(1, ['Decrease in auROC for file No %d vs file No %d electrode %d bandwidth No: %d\n'],file_pairs(fps,1),file_pairs(fps,2),elec,bwii);
                                    end
                                    
                                    %Are the delta dB LFP's different?
                                    
                                    %Hit
                                    
                                    p_val(no_dBs,bwii)=ranksum(this_delta_dB_powerpreHit,this_delta_dB_powerpostHit);
                                    p_vals=[p_vals p_val(no_dBs,bwii)];
                                    groupNopre(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                    groupNopost(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                    events(no_dBs)=1;
                                    
                                    
                                    %CR
                                    p_val(no_dBs+1,bwii)=ranksum(this_delta_dB_powerpreCR,this_delta_dB_powerpostCR);
                                    p_vals=[p_vals p_val(no_dBs+1,bwii)];
                                    groupNopre(no_dBs+1)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                    groupNopost(no_dBs+1)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                    events(no_dBs+1)=2;
                                    
                                    if p_val(no_dBs,bwii)<0.05
                                        dB_power_changeHit(no_ROCs)=1;
                                    else
                                        dB_power_changeHit(no_ROCs)=0;
                                    end
                                    
                                    if p_val(no_dBs+1,bwii)<0.05
                                        dB_power_changeCR(no_ROCs)=1;
                                    else
                                        dB_power_changeCR(no_ROCs)=0;
                                    end
                                    
                                    %Plot the points and save the data
                                    if groupNopre(no_dBs)==1
                                        
                                        
                                        %Hit, all points
                                        delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                        delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                        
                                        
                                        %CR, all points
                                        delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                        delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                        
                                    else
                                        if groupNopre(no_dBs)==3
                                            
                                            %Hit, all points
                                            delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                            delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                            
                                            
                                            %CR, all points
                                            delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                            delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                            %                                             figure(bwii+4+12)
                                            %                                             hold on
                                            %                                             plot([3 4],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        end
                                    end
                                end
                                
                                no_dBs=no_dBs+2;
                                
                            else
                                
                                if (sum(trials_in_event_preHit)<min_trials_per_event)
                                    fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_preHit),event1, min_trials_per_event,file_pairs(fps,1),elec);
                                end
                                
                                if (sum(trials_in_event_preCR)<min_trials_per_event)
                                    fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_preCR),event2, min_trials_per_event,file_pairs(fps,1),elec);
                                end
                                
                                if (sum(trials_in_event_postHit)<min_trials_per_event)
                                    fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_postHit),event1, min_trials_per_event,file_pairs(fps,2),elec);
                                end
                                
                                if (sum(trials_in_event_postCR)<min_trials_per_event)
                                    fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_postCR),event2, min_trials_per_event,file_pairs(fps,2),elec);
                                end
                                
                            end
                            
                        else
                            
                            if (length(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventLFPPower(1,:))<trials_to_process)
                                fprintf(1, ['%d trials fewer than %d trials to process for file No %d electrode %d\n'],length(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventLFPPower(1,:)),trials_to_process,file_pairs(fps,1),elec);
                            end
                            
                            if (length(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventLFPPower(1,:))<trials_to_process)
                                fprintf(1, ['%d trials fewer than %d trials to process for file No %d electrode %d\n'],length(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventLFPPower(1,:)),trials_to_process,file_pairs(fps,1),elec);
                            end
                            
                            
                        end
                    else
                        
                        if isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower)
                            fprintf(1, ['Empty allPower for file No %d electrode %d\n'],file_pairs(fps,1),elec);
                        end
                        
                        if isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower)
                            fprintf(1, ['Empty allPower for file No %d electrode %d\n'],file_pairs(fps,2),elec);
                        end
                        
                    end
                    
                else
                    
                    if isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref))
                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],file_pairs(fps,1),elec);
                    end
                    
                    if isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref))
                        fprintf(1, ['Empty lfpevpairfor file No %d electrode %d\n'],file_pairs(fps,2),elec);
                    end
                    
                end
            end
            
        end
        fprintf(1, '\n\n')
        
        
        pFDRROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for significant difference of auROC p value from 0.5  = %d\n\n'],pFDRROC);
        
        
        %Now plot the bar graphs and do anovan for LFP power
        p_vals_anovan=[];
        pvals_ancova=[];
        pvals_auROCancova=[];
        for bwii=1:4
            
            
            %Do ancova for auROC auROCpre
            this_auROCpre=[];
            this_auROCpre=auROCpre((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))';
            this_auROCpre=[this_auROCpre; auROCpre((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))'];
            
            
            this_auROCpost=[];
            this_auROCpost=auROCpost((ROCgroupNopre==1)&(ROCbandwidthpost==bwii))';
            this_auROCpost=[this_auROCpost; auROCpost((ROCgroupNopre==3)&(ROCbandwidthpost==bwii))'];
            
            pre_post=[];
            pre_post=[zeros(sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)),1); ones(sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)),1)];
            
            
            [h,atab,ctab,stats] = aoctool(this_auROCpre,this_auROCpost,pre_post,0.05,'','','','off');
            
            
            pvals_auROCancova=[pvals_auROCancova atab{4,6}];
            fprintf(1, ['ancova auROC p value ' freq_names{bwii} ' = %d\n\n'],atab{4,6});
            
            %Do ancova figure for auROC
            figure(10+bwii)
            h1=plot(auROCpre((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)),auROCpost((ROCgroupNopre==1)&(ROCbandwidthpost==bwii)),'or','MarkerFace','r');
            hold on
            h2=plot(auROCpre((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)),auROCpost((ROCgroupNopre==3)&(ROCbandwidthpost==bwii)),'ob','MarkerFace','b');
            
            slope_pre=ctab{5,2}+ctab{6,2};
            int_pre=ctab{2,2}+ctab{3,2};
            min_x=min([min(auROCpre((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) min(auROCpre((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
            max_x=max([max(auROCpre((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) max(auROCpre((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
            x=[-0.2 0.5];
            plot(x,slope_pre*x+int_pre,'-r','LineWidth',2)
            
            slope_post=ctab{5,2}+ctab{7,2};
            int_post=ctab{2,2}+ctab{4,2};
            x=[-0.2 0.5];
            plot(x,slope_post*x+int_post,'-b','LineWidth',2)
            
            plot([-0.2 0.5],[-0.2 0.5],'-k','LineWidth',2)
            
            title(['post vs pre auROC for ' freq_names{bwii} ])
            xlabel('pre auROC')
            ylabel('post auROC')
            legend([h1 h2],'halo','no halo')
            xlim([-0.2 0.5])
            ylim([-0.2 0.5])
            ax=gca;
            ax.LineWidth=3;
        end
        
        %         pFDRanovan=drsFDRpval(p_vals_anovan);
        %         fprintf(1, ['pFDR for anovan p value  = %d\n\n'],pFDRanovan);
        %
        %         pFDRancova=drsFDRpval(pvals_ancova);
        %         fprintf(1, ['pFDR for power dB ancova p value  = %d\n\n'], pFDRancova);
        
        pFDRauROCancova=drsFDRpval(pvals_auROCancova);
        fprintf(1, ['pFDR for auROC ancova p value  = %d\n\n'], pFDRauROCancova);
        
        fprintf(1, '\n\n')
        
        %         p_chi=[];
        %         for evTN1=1:length(eventType)
        %             fprintf(1, ['Significant changes in pairwise LFP power analysis for event: ' evTypeLabels{evTN1} '\n\n'])
        %             for bwii=1:4
        %                 for grs=grpre
        %                     num_sig(grs)=sum(p_val((events==evTN1)&(groupNopre==grs),bwii)<=0.05);
        %                     tot_num(grs)=sum((events==evTN1)&(grs==groupNopre));
        %                     fprintf(1, ['Number significant for ' freq_names{bwii} ' and ' handles_drgb.drgbchoices.group_no_names{grs} ' = %d of %d\n'],num_sig(grs),tot_num(grs));
        %                 end
        %                 [p, Q]= chi2test([num_sig(grpre(1)), tot_num(grpre(1))-num_sig(grpre(1)); num_sig(grpre(2)), tot_num(grpre(2))-num_sig(grpre(2))]);
        %                 fprintf(1, ['Chi squared p value  = %d\n\n'],p);
        %                 p_chi=[p_chi p];
        %             end
        %             fprintf(1, '\n\n\n')
        %         end
        %
        %         pFDRchi=drsFDRpval(p_chi);
        %         fprintf(1, ['pFDR for Chi squared p value  = %d\n\n'],pFDRchi);
        
        %Plot cumulative histos for auROCs
        dB_power_change=logical(dB_power_changeHit+dB_power_changeCR);
        figNo=0;
        p_val_ROC=[];
        pvals_auROCperm=[];
        
        x=0
        for bwii=1:4
            n_cum=0;
            this_legend=[];
            data_auROC=[];
            pre_post_auROC=[];
            gr_auROC=[];
            for grs=1:2
                if grs==1
                    try
                        close(figNo+1)
                    catch
                    end
                    figure(figNo+1)
                else
                    try
                        close(figNo+2)
                    catch
                    end
                    figure(figNo+2)
                end
                hold on
                
                %Plot the histograms
                maxauROC=max([max(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))) max(auROCpost((ROCgroupNopost==grpost(grs))&(ROCbandwidthpost==bwii)))]);
                minauROC=min([min(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))) min(auROCpost((ROCgroupNopost==grpost(grs))&(ROCbandwidthpost==bwii)))]);
                edges=[-0.5:0.05:0.5];
                pos2=[0.1 0.1 0.6 0.8];
                subplot('Position',pos2)
                hold on
                
                h2=histogram(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)),edges);
                h2.FaceColor='b';
                h1=histogram(auROCpost((ROCgroupNopost==grpost(grs))&(ROCbandwidthpost==bwii)),edges);
                h1.FaceColor='r';
                
                xlabel('auROC')
                ylabel('# of electrodes')
                legend('Pre','Laser')
                if grs==1
                    title(['auROC DBh Cre x halo for ' freq_names{bwii}])
                else
                    title(['auROC DBh Cre for ' freq_names{bwii}])
                end
                xlim([-0.3 0.6])
                ylim([0 45])
                ax=gca;
                ax.LineWidth=3;
                %                 if grs==1
                %                     ylim([0 30])
                %                 else
                %                     ylim([0 40])
                %                 end
                
                %Plot the single electrodes
                pos2=[0.8 0.1 0.1 0.8];
                subplot('Position',pos2)
                hold on
                for ii=1:length(auROCpre)
                    if (ROCgroupNopre(ii)==grpre(grs))&(ROCbandwidthpre(ii)==bwii)
                        plot([0 1],[auROCpre(ii) auROCpost(ii)],'-o', 'Color',[0.7 0.7 0.7])
                    end
                end
                
                
                plot([0 1],[mean(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))) mean(auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))],'-k','LineWidth', 3)
                CI = bootci(1000, @mean, auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                plot([0 0],CI,'-b','LineWidth',3)
                plot(0,mean(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))),'ob','MarkerSize', 10,'MarkerFace','b')
                CI = bootci(1000, @mean, auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                plot([1 1],CI,'-r','LineWidth',3)
                plot(1,mean(auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))),'or','MarkerSize', 10,'MarkerFace','r')
                ylabel('auROC')
                ylim([-0.2 0.5])
                ax=gca;
                ax.LineWidth=3;
                %Do the statistics for auROC differences
                %                 a={auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))' auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))'};
                %                 mode_statcond='perm';
                %                 [F df pval_auROCperm] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
                %
                pval_auROCperm=ranksum(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)), auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                
                if grs==1
                    fprintf(1, ['p value for premuted anovan for auROC DBH Cre x halo pre vs laser ' freq_names{bwii} '= %d\n'],  pval_auROCperm);
                else
                    fprintf(1, ['p value for premuted anovan for auROC DBH Cre pre vs laser ' freq_names{bwii} '= %d\n'],  pval_auROCperm);
                end
                pvals_auROCperm=[pvals_auROCperm pval_auROCperm];
                
                %Save the data for anovan interaction
                %Pre
                data_auROC=[data_auROC auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))];
                gr_auROC=[gr_auROC grs*ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
                pre_post_auROC=[pre_post_auROC ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
                
                %Post
                data_auROC=[data_auROC auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))];
                gr_auROC=[gr_auROC grs*ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
                pre_post_auROC=[pre_post_auROC 2*ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
            end
            figNo=figNo+2;
            x=x+3;
            
            %Calculate anovan for inteaction
            [p,tbl,stats]=anovan(data_auROC,{pre_post_auROC gr_auROC},'model','interaction','varnames',{'pre_vs_post','halo_vs_no_halo'},'display','off');
            fprintf(1, ['p value for anovan auROC interaction for ' freq_names{bwii} '= %d\n'],  p(3));
            p_aovan_int(bwii)=p(3);
            
        end
        
        pFDRauROC=drsFDRpval(pvals_auROCperm);
        fprintf(1, ['pFDR for auROC  = %d\n\n'],pFDRauROC);
        
        pFDRauROCint=drsFDRpval(p_aovan_int);
        fprintf(1, ['pFDR for auROC anovan interaction  = %d\n\n'],pFDRauROCint);
        
        
        
        %         p_perCorr=ranksum(perCorr_pre,perCorr_post);
        %         fprintf(1, '\np value for ranksum test for percent correct= %d\n\n',p_perCorr);
        %
        
        
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) '_out.mat'],'perCorr_pre','perCorr_post','group_pre', 'group_post');
        pfft=1;
        
    case 6
        %For the proficient mice in the first and last sessions
        %plot the ERP LFP spectrum for S+ vs S-, plot ERP LFP power for S+ vs S- for each electrode and plot ERP LFP auROCs
        %NOTE: This does the analysis in all the files and DOES not distinguish between groups!!!
        no_dBs=1;
        delta_dB_power=[];
        no_ROCs=0;
        ROCout=[];
        p_vals_ROC=[];
        delta_dB_powerEv1=[];
        no_Ev1=0;
        noWB=0;
        delta_dB_powerEv1WB=[];
        delta_dB_powerEv2WB=[];
        shift_ii=floor(length(handles_drgb.drgb.lfpevpair(1).out_times)/2)+1+shift_from_event;
        
        NoLicksEv1=[];
        NoLicksEv2=[];
        lickTimesEv1=[];
        lickTimesEv2=[];
        lickTrialsEv1=0;
        lickTrialsEv2=0;
        
        
        
        fprintf(1, ['Pairwise auROC analysis for Fig 1 of Daniel''s paper\n\n'])
        p_vals=[];
        for fileNo=1:length(files)
            %Analyze the licks
            elec=1;
            lfpodNo=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
            
            trials=length(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventERP(1,:));
            mask=logical([zeros(1,trials-trials_to_process) ones(1,trials_to_process)]);
            trials_in_eventEv1=(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventERP(event1,:)==1);
            trials_in_eventEv2=(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventERP(event2,:)==1);
            
            NoLicksEv1=[NoLicksEv1 handles_drgb.drgb.lfpevpair(lfpodNo).no_events_per_trial(trials_in_eventEv1&mask)];
            NoLicksEv2=[NoLicksEv2 handles_drgb.drgb.lfpevpair(lfpodNo).no_events_per_trial(trials_in_eventEv2&mask)];
            
            %Get times for Ev1
            for trialNo=1:length(trials_in_eventEv1)
                if trials_in_eventEv1(trialNo)==1
                    lickTimesEv1=[lickTimesEv1 handles_drgb.drgb.lfpevpair(lfpodNo).t_per_event_per_trial(trialNo,1:handles_drgb.drgb.lfpevpair(lfpodNo).no_events_per_trial(trialNo))];
                    lickTrialsEv1=lickTrialsEv1+1;
                end
                if trials_in_eventEv2(trialNo)==1
                    lickTimesEv2=[lickTimesEv2 handles_drgb.drgb.lfpevpair(lfpodNo).t_per_event_per_trial(trialNo,1:handles_drgb.drgb.lfpevpair(lfpodNo).no_events_per_trial(trialNo))];
                    lickTrialsEv2=lickTrialsEv2+1;
                end
            end
            
            pffft=1;
            
            for elec=1:16
                
                lfpodNo=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                
                if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo)))
                    
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP))
                        
                        if (length(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventERP(1,:))>=trials_to_process)
                            
                            trials=length(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventERP(1,:));
                            if front_mask==1
                                mask=logical([ones(1,trials_to_process) zeros(1,trials-trials_to_process)]);
                            else
                                mask=logical([zeros(1,trials-trials_to_process) ones(1,trials_to_process)]);
                            end
                            trials_in_eventEv1=(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventERP(event1,:)==1);
                            trials_in_eventEv2=(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventERP(event2,:)==1);
                            %Make sure there is at least one lick event
                            %within each trial
                            trials_with_event=(handles_drgb.drgb.lfpevpair(lfpodNo).no_events_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNo).no_events_per_trial<=max_events_per_sec)...
                                &(handles_drgb.drgb.lfpevpair(lfpodNo).no_ref_evs_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNo).no_ref_evs_per_trial<=max_events_per_sec);
                            
                            if (sum(trials_in_eventEv1&trials_with_event&mask)>=min_trials_per_event) & (sum(trials_in_eventEv2&trials_with_event&mask)>=min_trials_per_event)
                                
                                
                                
                                
                                this_dB_powerEv1=zeros(sum(trials_in_eventEv1&trials_with_event&mask),length(frequency));
                                this_dB_powerEv1(:,:)=handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP(trials_in_eventEv1&trials_with_event&mask,:,shift_ii);
                                
                                % Ev2
                                this_dB_powerEv2=zeros(sum(trials_in_eventEv2&trials_with_event&mask),length(frequency));
                                this_dB_powerEv2(:,:)=handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP(trials_in_eventEv2&trials_with_event&mask,:,shift_ii);
                                
                                
                                %Wide band spectrum
                                noWB=noWB+1;
                                
                                delta_dB_powerEv1WB(noWB,:)=mean(this_dB_powerEv1,1);
                                delta_dB_powerEv2WB(noWB,:)=mean(this_dB_powerEv2,1);
                                
                                
                                %Do per badwidth analysis
                                for bwii=1:no_bandwidths
                                    
                                    no_ROCs=no_ROCs+1;
                                    this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                    
                                    %Enter the  Ev1
                                    this_delta_dB_powerEv1=zeros(sum(trials_in_eventEv1&trials_with_event&mask),1);
                                    this_delta_dB_powerEv1=mean(this_dB_powerEv1(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:sum(trials_in_eventEv1&trials_with_event&mask),1)=this_delta_dB_powerEv1;
                                    roc_data(1:sum(trials_in_eventEv1&trials_with_event&mask),2)=zeros(sum(trials_in_eventEv1&trials_with_event&mask),1);
                                    
                                    %Enter  Ev2
                                    total_trials=sum(trials_in_eventEv1&trials_with_event&mask)+sum(trials_in_eventEv2&trials_with_event&mask);
                                    this_delta_dB_powerEv2=zeros(sum(trials_in_eventEv2&trials_with_event&mask),1);
                                    this_delta_dB_powerEv2=mean(this_dB_powerEv2(:,this_band),2);
                                    roc_data(sum(trials_in_eventEv1&trials_with_event&mask)+1:total_trials,1)=this_delta_dB_powerEv2;
                                    roc_data(sum(trials_in_eventEv1&trials_with_event&mask)+1:total_trials,2)=ones(sum(trials_in_eventEv2&trials_with_event&mask),1);
                                    
                                    
                                    %Find  ROC
                                    ROCout(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCout(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
                                    ROCgroupNo(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNo).fileNo);
                                    ROCout(no_ROCs).timeWindow=winNo;
                                    ROCbandwidth(no_ROCs)=bwii;
                                    auROC(no_ROCs)=ROCout(no_ROCs).roc.AUC-0.5;
                                    p_valROC(no_ROCs)=ROCout(no_ROCs).roc.p;
                                    
                                    p_vals_ROC=[p_vals_ROC ROCout(no_ROCs).roc.p];
                                    
                                    
                                    delta_dB_powerEv1(no_ROCs)=mean(this_delta_dB_powerEv1);
                                    delta_dB_powerEv2(no_ROCs)=mean(this_delta_dB_powerEv2);
                                    
                                    
                                    %Plot this point
                                    figure(bwii+1)
                                    pos2=[0.8 0.1 0.1 0.8];
                                    subplot('Position',pos2)
                                    hold on
                                    plot([1 0],[delta_dB_powerEv1(no_ROCs) delta_dB_powerEv2(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                    set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                                    
                                    
                                end
                                
                                
                                
                            else
                                
                                if sum(trials_in_eventEv1&trials_with_event&mask)<min_trials_per_event
                                    fprintf(1, ['%d trials with lick events in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_eventEv1&trials_with_event&mask),event1, min_trials_per_event,files(fileNo),elec);
                                end
                                
                                if sum(trials_in_eventEv2&trials_with_event&mask)<min_trials_per_event
                                    fprintf(1, ['%d trials with lick events in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_eventEv2&trials_with_event&mask),event2, min_trials_per_event,files(fileNo),elec);
                                end
                                
                            end
                            
                        else
                            
                            fprintf(1, ['%d trials fewer than %d trials to process for file No %d electrode %d\n'],length(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventERP(1,:)),trials_to_process,files(fileNo),elec);
                            
                        end
                    else
                        
                        fprintf(1, ['Empty allPower for file No %d electrode %d\n'],files(fileNo),elec);
                        
                    end
                    
                    
                else
                    fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],files(fileNo),elec);
                    
                    
                end
            end
            
        end
        fprintf(1, '\n\n')
        
        
        %Now plot the bounded line for
        
        %Calculate the mean and 95% CI for Ev1
        dB_Ev1_ci=zeros(length(frequency),2);
        for ifreq=1:length(frequency)
            %             pd=fitdist(delta_dB_powerEv1WB(:,ifreq),'Normal');
            %             ci=paramci(pd);
            %             dB_Ev1_ci(ifreq)=pd.mu-ci(1,1);
            dB_Ev1_mean(ifreq)=mean(delta_dB_powerEv1WB(:,ifreq));
            CI = bootci(1000, @mean, delta_dB_powerEv1WB(:,ifreq));
            dB_Ev1_ci(ifreq,1)=CI(2)-dB_Ev1_mean(ifreq);
            dB_Ev1_ci(ifreq,2)=-(CI(1)-dB_Ev1_mean(ifreq));
        end
        
        figure(1)
        [hl1, hp1] = boundedline(frequency,dB_Ev1_mean, dB_Ev1_ci, 'r');
        
        %Calculate the mean and 95% CI for Ev2
        dB_Ev2_ci=zeros(length(frequency),2);
        for ifreq=1:length(frequency)
            dB_Ev2_mean(ifreq)=mean(delta_dB_powerEv2WB(:,ifreq));
            CI = bootci(1000, @mean, delta_dB_powerEv2WB(:,ifreq));
            dB_Ev2_ci(ifreq,1)=CI(2)-dB_Ev2_mean(ifreq);
            dB_Ev2_ci(ifreq,2)=-(CI(1)-dB_Ev2_mean(ifreq));
        end
        
        hold on
        [hl2, hp2] = boundedline(frequency,dB_Ev2_mean, dB_Ev2_ci, 'b');
        xlabel('Frequency (Hz)')
        ylabel('delta Power (dB)')
        legend([hl1 hl2],'S+','S-')
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        
        %Now plot the histograms and the average
        for bwii=1:4
            %Plot the average
            figure(bwii+1)
            pos2=[0.8 0.1 0.1 0.8];
            subplot('Position',pos2)
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            plot([1 0],[mean(delta_dB_powerEv1(ROCbandwidth==bwii)) mean(delta_dB_powerEv2(ROCbandwidth==bwii))],'-k','LineWidth', 3)
            CI = bootci(1000, @mean, delta_dB_powerEv1(ROCbandwidth==bwii));
            plot([1 1],CI,'-r','LineWidth',3)
            plot(1,mean(delta_dB_powerEv1(ROCbandwidth==bwii)),'or','MarkerSize', 10,'MarkerFace','r')
            CI = bootci(1000, @mean, delta_dB_powerEv2(ROCbandwidth==bwii));
            plot([0 0],CI,'-b','LineWidth',3)
            plot(0,mean(delta_dB_powerEv2(ROCbandwidth==bwii)),'ob','MarkerSize', 10,'MarkerFace','b')
            ylabel('delta Power (dB)')
            ylim([-10 15])
            
            %Plot the histograms
            
            maxdB=max([max(delta_dB_powerEv1(ROCbandwidth==bwii)) max(delta_dB_powerEv2(ROCbandwidth==bwii))]);
            mindB=min([min(delta_dB_powerEv1(ROCbandwidth==bwii)) min(delta_dB_powerEv2(ROCbandwidth==bwii))]);
            edges=[-15:1:15];
            pos2=[0.1 0.1 0.6 0.8];
            subplot('Position',pos2)
            hold on
            
            h1=histogram(delta_dB_powerEv2(ROCbandwidth==bwii),edges);
            h1.FaceColor='b';
            h2=histogram(delta_dB_powerEv1(ROCbandwidth==bwii),edges);
            h2.FaceColor='r';
            xlabel('delta Power (dB)')
            ylabel('# of electrodes')
            legend('S-','S+')
            xlim([-12 12])
            ylim([0 70])
            title(freq_names{bwii})
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            
            
            a={ delta_dB_powerEv1(ROCbandwidth==bwii)' delta_dB_powerEv2(ROCbandwidth==bwii)'};
            mode_statcond='perm';
            [F df pvals_perm(bwii)] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for premuted anovan dB delta power S+ vs S- ' freq_names{bwii} '= %d\n'],  pvals_perm(bwii));
            
        end
        
        pFDRanovan=drsFDRpval(pvals_perm);
        fprintf(1, ['pFDR for premuted anovan p value  = %d\n\n'],pFDRanovan);
        
        
        
        fprintf(1, '\n\n')
        
        
        pFDRauROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for auROC  = %d\n\n'],pFDRauROC);
        %Plot cumulative histos for auROCs
        
        figNo=5;
        p_val_ROC=[];
        edges=-0.5:0.05:0.5;
        
        for bwii=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            n_cum=0;
            this_legend=[];
            
            histogram(auROC(( p_valROC>pFDRauROC)&(ROCbandwidth==bwii)),edges)
            histogram(auROC(( p_valROC<=pFDRauROC)&(ROCbandwidth==bwii)),edges)
            legend('auROC not singificant','auROC significant')
            title(['Histogram for ' freq_names{bwii} ' auROC for LFPs'])
            xlim([-0.2 0.6])
            ylim([0 30])
        end
        
        
        
        %Plot percent significant ROC
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        figure(figNo)
        
        hold on
        
        for bwii=1:4
            bar(bwii,100*sum(( p_valROC<=pFDRauROC)&(ROCbandwidth==bwii))/sum((ROCbandwidth==bwii)))
            auROC_sig.sig(bwii)=sum(( p_valROC<=pFDRauROC)&(ROCbandwidth==bwii));
            auROC_sig.not_sig(bwii)=sum((ROCbandwidth==bwii))-sum(( p_valROC<=pFDRauROC)&(ROCbandwidth==bwii));
        end
        title('Percent auROC significantly different from zero')
        ylim([0 100])
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        
        
        %Plot the lick time histogram
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        figure(figNo)
        
        hold on
        edges=[handles_drgb.drgb.lfpevpair(lfpodNo).timeStart:0.1:handles_drgb.drgb.lfpevpair(lfpodNo).timeEnd];
        h1=histogram(lickTimesEv1);
        h1.FaceColor='r';
        h2=histogram(lickTimesEv2);
        h2.FaceColor='b';
        xlabel('Time(sec)')
        ylabel('# of licks')
        legend(evTypeLabels{1},evTypeLabels{2})
        xlim([handles_drgb.drgb.lfpevpair(lfpodNo).timeStart handles_drgb.drgb.lfpevpair(lfpodNo).timeEnd])
        title('Lick time histogram')
        
        fprintf(1, ['Mean number of licks per trial for ' evTypeLabels{1} '= %d\n'],mean(NoLicksEv1))
        fprintf(1, ['Mean number of licks per trial for ' evTypeLabels{2} '= %d\n'],mean(NoLicksEv2))
        
        pffft=1;
        
    case 7
        %Compare auROC in the last few trials of the last session file with
        %first few trials of the first session
        %Generates Fig. 3 for Daniel's paper. first vs last.
        no_dBs=1;
        delta_dB_power_fp1=[];
        no_ROCs=0;
        ROCoutfp1=[];
        ROCoutfp2=[];
        p_vals_ROC=[];
        delta_dB_powerfp1Ev1=[];
        no_Ev1=0;
        pvals_auROCperm=[];
        pvals_dBperm=[];
        perCorr_fp1=[];
        perCorr_fp2=[];
        shift_ii=floor(length(handles_drgb.drgb.lfpevpair(1).out_times)/2)+1+shift_from_event;
        
        
        fprintf(1, ['Pairwise auROC log power LFP ERP analysis for ' evTypeLabels{1} ' and ' evTypeLabels{2} ' LFP power\n\n'])
        p_vals=[];
        
        no_file_pairs=length(file_pairs);
        
        for fps=1:no_file_pairs
            
            
            for elec=1:16
                
                lfpodNofp1=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                lfpodNofp2=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                
                
                if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNofp1)))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNofp2)))
                    
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNofp1).log_P_tERP))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNofp2).log_P_tERP))
                        
                        if (length(handles_drgb.drgb.lfpevpair(lfpodNofp1).which_eventERP(1,:))>=trials_to_process) &...
                                (length(handles_drgb.drgb.lfpevpair(lfpodNofp2).which_eventERP(1,:))>=trials_to_process)
                            
                            length_fp1=length(handles_drgb.drgb.lfpevpair(lfpodNofp1).which_eventERP(1,:));
                            if front_mask(1)==1
                                fp1_mask=logical([ones(1,trials_to_process) zeros(1,length_fp1-trials_to_process)]);
                            else
                                fp1_mask=logical([zeros(1,length_fp1-trials_to_process) ones(1,trials_to_process)]);
                            end
                            trials_in_event_fp1Ev1=(handles_drgb.drgb.lfpevpair(lfpodNofp1).which_eventERP(event1,:)==1);
                            trials_in_event_fp1Ev2=(handles_drgb.drgb.lfpevpair(lfpodNofp1).which_eventERP(event2,:)==1);
                            
                            trials_with_eventfp1=(handles_drgb.drgb.lfpevpair(lfpodNofp1).no_events_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNofp1).no_events_per_trial<=max_events_per_sec)...
                                &(handles_drgb.drgb.lfpevpair(lfpodNofp1).no_ref_evs_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNofp1).no_ref_evs_per_trial<=max_events_per_sec);
                            
                            
                            
                            length_fp2=length(handles_drgb.drgb.lfpevpair(lfpodNofp2).which_eventERP(1,:));
                            if front_mask(2)==1
                                fp2_mask=logical([ones(1,trials_to_process) zeros(1,length_fp2-trials_to_process)]);
                            else
                                fp2_mask=logical([zeros(1,length_fp2-trials_to_process) ones(1,trials_to_process)]);
                            end
                            trials_in_event_fp2Ev1=(handles_drgb.drgb.lfpevpair(lfpodNofp2).which_eventERP(event1,:)==1);
                            trials_in_event_fp2Ev2=(handles_drgb.drgb.lfpevpair(lfpodNofp2).which_eventERP(event2,:)==1);
                            
                            
                            trials_with_eventfp2=(handles_drgb.drgb.lfpevpair(lfpodNofp2).no_events_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNofp2).no_events_per_trial<=max_events_per_sec)...
                                &(handles_drgb.drgb.lfpevpair(lfpodNofp2).no_ref_evs_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNofp2).no_ref_evs_per_trial<=max_events_per_sec);
                            
                            if (sum(trials_in_event_fp1Ev1&fp1_mask&trials_with_eventfp1)>=min_trials_per_event) & (sum(trials_in_event_fp1Ev2&fp1_mask&trials_with_eventfp1)>=min_trials_per_event) & ...
                                    (sum(trials_in_event_fp2Ev1&fp2_mask&trials_with_eventfp2)>=min_trials_per_event) & (sum(trials_in_event_fp2Ev2&fp2_mask&trials_with_eventfp2)>=min_trials_per_event)
                            
                                
                                %fp1 Ev1
                                this_dB_powerfp1Ev1=zeros(sum(trials_in_event_fp1Ev1&fp1_mask&trials_with_eventfp1),length(frequency));
                                this_dB_powerfp1Ev1(:,:)=handles_drgb.drgb.lfpevpair(lfpodNofp1).log_P_tERP(trials_in_event_fp1Ev1&fp1_mask&trials_with_eventfp1,:,shift_ii);
                                
                                %fp1 Ev2
                                this_dB_powerfp1Ev2=zeros(sum(trials_in_event_fp1Ev2&fp1_mask&trials_with_eventfp1),length(frequency));
                                this_dB_powerfp1Ev2(:,:)=handles_drgb.drgb.lfpevpair(lfpodNofp1).log_P_tERP(trials_in_event_fp1Ev2&fp1_mask&trials_with_eventfp1,:,shift_ii);
                                
                                
                                %fp2 Ev1
                                this_dB_powerfp2Ev1=zeros(sum(trials_in_event_fp2Ev1&fp2_mask&trials_with_eventfp2),length(frequency));
                                this_dB_powerfp2Ev1(:,:)=handles_drgb.drgb.lfpevpair(lfpodNofp2).log_P_tERP(trials_in_event_fp2Ev1&fp2_mask&trials_with_eventfp2,:,shift_ii);
                                
                                %fp2 Ev2
                                this_dB_powerfp2Ev2=zeros(sum(trials_in_event_fp2Ev2&fp2_mask&trials_with_eventfp2),length(frequency));
                                this_dB_powerfp2Ev2(:,:)=handles_drgb.drgb.lfpevpair(lfpodNofp2).log_P_tERP(trials_in_event_fp2Ev2&fp2_mask&trials_with_eventfp2,:,shift_ii);
                                
                                for bwii=1:no_bandwidths
                                    
                                    no_ROCs=no_ROCs+1;
                                    this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
<<<<<<< HEAD
                                     
=======
                                    
>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855
                                    %Enter the fp1 Ev1
                                    this_delta_dB_powerfp1Ev1=zeros(sum(trials_in_event_fp1Ev1&fp1_mask&trials_with_eventfp1),1);
                                    this_delta_dB_powerfp1Ev1=mean(this_dB_powerfp1Ev1(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:sum(trials_in_event_fp1Ev1&fp1_mask&trials_with_eventfp1),1)=this_delta_dB_powerfp1Ev1;
                                    roc_data(1:sum(trials_in_event_fp1Ev1&fp1_mask&trials_with_eventfp1),2)=zeros(sum(trials_in_event_fp1Ev1&fp1_mask&trials_with_eventfp1),1);
                                    
                                    %Enter fp1 Ev2
                                    total_trials=sum(trials_in_event_fp1Ev2&fp1_mask&trials_with_eventfp1)+sum(trials_in_event_fp1Ev1&fp1_mask&trials_with_eventfp1);
                                    this_delta_dB_powerfp1Ev2=zeros(sum(trials_in_event_fp1Ev2&fp1_mask&trials_with_eventfp1),1);
                                    this_delta_dB_powerfp1Ev2=mean(this_dB_powerfp1Ev2(:,this_band),2);
                                    roc_data(sum(trials_in_event_fp1Ev1&fp1_mask&trials_with_eventfp1)+1:total_trials,1)=this_delta_dB_powerfp1Ev2;
                                    roc_data(sum(trials_in_event_fp1Ev1&fp1_mask&trials_with_eventfp1)+1:total_trials,2)=ones(sum(trials_in_event_fp1Ev2&fp1_mask&trials_with_eventfp1),1);
                                    
                                    
                                    %Find fp1 ROC
                                    ROCoutfp1(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCoutfp1(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNofp1).fileNo;
                                    ROCgroupNofp1(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNofp1).fileNo);
                                    ROCoutfp1(no_ROCs).timeWindow=winNo;
                                    ROCbandwidthfp1(no_ROCs)=bwii;
                                    auROCfp1(no_ROCs)=ROCoutfp1(no_ROCs).roc.AUC-0.5;
                                    p_valROCfp1(no_ROCs)=ROCoutfp1(no_ROCs).roc.p;
                                    
                                    p_vals_ROC=[p_vals_ROC ROCoutfp1(no_ROCs).roc.p];
                                    
                                    %Enter the fp2 Ev1
                                    this_delta_dB_powerfp2Ev1=zeros(sum(trials_in_event_fp2Ev1&fp2_mask&trials_with_eventfp2),1);
                                    this_delta_dB_powerfp2Ev1=mean(this_dB_powerfp2Ev1(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:sum(trials_in_event_fp2Ev1&fp2_mask&trials_with_eventfp2),1)=this_delta_dB_powerfp2Ev1;
                                    roc_data(1:sum(trials_in_event_fp2Ev1&fp2_mask&trials_with_eventfp2),2)=zeros(sum(trials_in_event_fp2Ev1&fp2_mask&trials_with_eventfp2),1);
                                    
                                    %Enter fp2 Ev2
                                    total_trials=sum(trials_in_event_fp2Ev1&fp2_mask&trials_with_eventfp2)+sum(trials_in_event_fp2Ev2&fp2_mask&trials_with_eventfp2);
                                    this_delta_dB_powerfp2Ev2=zeros(sum(trials_in_event_fp2Ev2&fp2_mask&trials_with_eventfp2),1);
                                    this_delta_dB_powerfp2Ev2=mean(this_dB_powerfp2Ev2(:,this_band),2);
                                    roc_data(sum(trials_in_event_fp2Ev1&fp2_mask&trials_with_eventfp2)+1:total_trials,1)=this_delta_dB_powerfp2Ev2;
                                    roc_data(sum(trials_in_event_fp2Ev1&fp2_mask&trials_with_eventfp2)+1:total_trials,2)=ones(sum(trials_in_event_fp2Ev2&fp2_mask&trials_with_eventfp2),1);
                                    
                                    
                                    %Find fp2 ROC
                                    ROCoutfp2(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCoutfp2(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNofp2).fileNo;
                                    ROCgroupNofp2(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNofp2).fileNo);
                                    ROCoutfp2(no_ROCs).timeWindow=winNo;
                                    ROCbandwidthfp2(no_ROCs)=bwii;
                                    auROCfp2(no_ROCs)=ROCoutfp2(no_ROCs).roc.AUC-0.5;
                                    p_valROCfp2(no_ROCs)=ROCoutfp2(no_ROCs).roc.p;
                                    
                                    p_vals_ROC=[p_vals_ROC ROCoutfp2(no_ROCs).roc.p];
                                    
                                    
                                    %Are the delta dB LFP's different?
                                    
                                    %Ev1
                                    
                                    p_val(no_dBs,bwii)=ranksum(this_delta_dB_powerfp1Ev1,this_delta_dB_powerfp2Ev1);
                                    p_vals=[p_vals p_val(no_dBs,bwii)];
                                    groupNofp1(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                    groupNofp2(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                    events(no_dBs)=1;
                                    
                                    
                                    %Ev2
                                    p_val(no_dBs+1,bwii)=ranksum(this_delta_dB_powerfp1Ev2,this_delta_dB_powerfp2Ev2);
                                    p_vals=[p_vals p_val(no_dBs+1,bwii)];
                                    groupNofp1(no_dBs+1)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                    groupNofp2(no_dBs+1)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                    events(no_dBs+1)=2;
                                    
                                    if p_val(no_dBs,bwii)<0.05
                                        dB_power_changeEv1(no_ROCs)=1;
                                    else
                                        dB_power_changeEv1(no_ROCs)=0;
                                    end
                                    
                                    if p_val(no_dBs+1,bwii)<0.05
                                        dB_power_changeEv2(no_ROCs)=1;
                                    else
                                        dB_power_changeEv2(no_ROCs)=0;
                                    end
                                    
                                    
                                    
                                    %Ev1, all points
                                    delta_dB_powerfp1Ev1(no_ROCs)=mean(this_delta_dB_powerfp1Ev1);
                                    delta_dB_powerfp2Ev1(no_ROCs)=mean(this_delta_dB_powerfp2Ev1);
                                    
                                    %Ev2, all points
                                    delta_dB_powerfp1Ev2(no_ROCs)=mean(this_delta_dB_powerfp1Ev2);
                                    delta_dB_powerfp2Ev2(no_ROCs)=mean(this_delta_dB_powerfp2Ev2);
                                    
                                    
                                end
                                
                                no_dBs=no_dBs+2;
                                
                            else
                                
                                if (sum(trials_in_event_fp1Ev1&fp1_mask&trials_with_eventfp1)<min_trials_per_event)
                                    fprintf(1, ['%d trials with lick events in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_fp1Ev1&fp1_mask&trials_with_eventfp1),event1, min_trials_per_event,file_pairs(fps,1),elec);
                                end
                                
                                if (sum(trials_in_event_fp1Ev2&fp1_mask&trials_with_eventfp1)<min_trials_per_event)
                                    fprintf(1, ['%d trials with lick events in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_fp1Ev2&fp1_mask&trials_with_eventfp1),event2, min_trials_per_event,file_pairs(fps,1),elec);
                                end
                                
                                if (sum(trials_in_event_fp2Ev1&fp2_mask&trials_with_eventfp2)<min_trials_per_event)
                                    fprintf(1, ['%d trials with lick events in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_fp2Ev1&fp2_mask&trials_with_eventfp2),event1, min_trials_per_event,file_pairs(fps,2),elec);
                                end
                                
                                if (sum(trials_in_event_fp2Ev2&fp2_mask&trials_with_eventfp2)<min_trials_per_event)
                                    fprintf(1, ['%d trials with lick events in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_event_fp2Ev2&fp2_mask&trials_with_eventfp2),event2, min_trials_per_event,file_pairs(fps,2),elec);
                                end
                                
                            end
                            
                        else
                            
                            if (length(handles_drgb.drgb.lfpevpair(lfpodNofp1).which_eventERP(1,:))<trials_to_process)
                                fprintf(1, ['%d trials fewer than %d trials to process for file No %d electrode %d\n'],length(handles_drgb.drgb.lfpevpair(lfpodNofp1).which_eventERP(1,:)),trials_to_process,file_pairs(fps,1),elec);
                            end
                            
                            if (length(handles_drgb.drgb.lfpevpair(lfpodNofp2).which_eventERP(1,:))<trials_to_process)
                                fprintf(1, ['%d trials fewer than %d trials to process for file No %d electrode %d\n'],length(handles_drgb.drgb.lfpevpair(lfpodNofp2).which_eventERP(1,:)),trials_to_process,file_pairs(fps,1),elec);
                            end
                            
                        end
                    else
                        
                        if isempty(handles_drgb.drgb.lfpevpair(lfpodNofp1).allPower)
                            fprintf(1, ['Empty allPower for file No %d electrode %d\n'],file_pairs(fps,1),elec);
                        end
                        
                        if isempty(handles_drgb.drgb.lfpevpair(lfpodNofp2).allPower)
                            fprintf(1, ['Empty allPower for file No %d electrode %d\n'],file_pairs(fps,2),elec);
                        end
                        
                    end
                    
                else
                    
                    if isempty(handles_drgb.drgb.lfpevpair(lfpodNofp1))
                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],file_pairs(fps,1),elec);
                    end
                    
                    if isempty(handles_drgb.drgb.lfpevpair(lfpodNofp2))
                        fprintf(1, ['Empty lfpevpairfor file No %d electrode %d\n'],file_pairs(fps,2),elec);
                    end
                    
                end
            end
            
        end
        fprintf(1, '\n\n')
        
        pFDRROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for significant difference of auROC p value from 0.5  = %d\n\n'],pFDRROC);
        
        
        
        fprintf(1, '\n\n')
        
        
        %Plot cumulative histos for auROCs
        dB_power_change=logical(dB_power_changeEv1+dB_power_changeEv2);
        figNo=0;
        p_val_ROC=[];
        
        try
            close(5)
        catch
        end
        figure(5)
        hold on
        x=0;
        
        for bwii=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            %Plot the histograms
            edges=[-0.5:0.05:0.5];
            pos2=[0.1 0.1 0.6 0.8];
            subplot('Position',pos2)
            hold on
            
            h2=histogram(auROCfp2(ROCbandwidthfp1==bwii),edges);
            h2.FaceColor='b';
            h1=histogram(auROCfp1(ROCbandwidthfp1==bwii),edges);
            h1.FaceColor='r';
            
            xlabel('auROC')
            ylabel('# of electrodes')
            legend(file_label{2},file_label{1})
            title(['auROC for ' freq_names{bwii}])
            xlim([-0.3 0.6])
            ylim([0 45])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Plot the single electrodes
            pos2=[0.8 0.1 0.1 0.8];
            subplot('Position',pos2)
            hold on
            for ii=1:length(auROCfp1)
                if ROCbandwidthfp1(ii)==bwii
                    plot([0 1],[auROCfp2(ii) auROCfp1(ii)],'-o', 'Color',[0.7 0.7 0.7])
                end
            end
            
            %PLot the mean and 95% CI
            plot([0 1],[mean(auROCfp2(ROCbandwidthfp1==bwii)) mean(auROCfp1(ROCbandwidthfp1==bwii))],'-k','LineWidth', 3)
            CI = bootci(1000, @mean, auROCfp2(ROCbandwidthfp1==bwii));
            plot([0 0],CI,'-b','LineWidth',3)
            plot(0,mean(auROCfp2(ROCbandwidthfp1==bwii)),'ob','MarkerSize', 10,'MarkerFace','b')
            CI = bootci(1000, @mean, auROCfp1(ROCbandwidthfp1==bwii));
            plot([1 1],CI,'-r','LineWidth',3)
            plot(1,mean(auROCfp1(ROCbandwidthfp1==bwii)),'or','MarkerSize', 10,'MarkerFace','r')
            ylabel('auROC')
            ylim([-0.2 0.5])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Do the statistics for auROC differences
            a={auROCfp2(ROCbandwidthfp1==bwii)' auROCfp1(ROCbandwidthfp1==bwii)'};
            mode_statcond='perm';
            [F df pval_auROCperm] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for permuted anovan for auROC S+ vs S- ' freq_names{bwii} '= %d\n\n'],  pval_auROCperm);
            pvals_auROCperm=[pvals_auROCperm pval_auROCperm];
            
            %Figure 5
            figure(5)
            
            percent_auROCfp2=100*sum(p_valROCfp2(ROCbandwidthfp1==bwii)<=pFDRROC)/sum(ROCbandwidthfp1==bwii);
            bar(x,percent_auROCfp2,'b')
            
            learn_sig(bwii)=sum(p_valROCfp2(ROCbandwidthfp1==bwii)<=pFDRROC);
            learn_not_sig(bwii)=sum(ROCbandwidthfp1==bwii)-sum(p_valROCfp2(ROCbandwidthfp1==bwii)<=pFDRROC);
            
            percent_auROCfp1=100*sum(p_valROCfp1(ROCbandwidthfp1==bwii)<=pFDRROC)/sum(ROCbandwidthfp1==bwii);
            bar(x+1,percent_auROCfp1,'r')
            
            prof_sig(bwii)=sum(p_valROCfp1(ROCbandwidthfp1==bwii)<=pFDRROC);
            prof_not_sig(bwii)=sum(ROCbandwidthfp1==bwii)-sum(p_valROCfp1(ROCbandwidthfp1==bwii)<=pFDRROC);
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            x=x+3;
            
        end
        
        figure(5)
        title('Percent singificant auROC')
        legend(file_label{2},file_label{1})
        ylim([0 100])
        
        pFDRanovanauROC=drsFDRpval(pval_auROCperm);
        fprintf(1, ['\npFDR for premuted anovan p value for difference between ' file_label{1} ' and ' file_label{2} ' for auROC = %d\n\n'],pFDRanovanauROC);
        
        
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix],'learn_sig','learn_not_sig','prof_sig','prof_not_sig');
        
        pffft=1;
        
        
    case 8
        
        %Compare auROC for ERP LFP in the last few trials of pre with first few trials of post
        %Used for New Fig. 7 of Daniel's paper
        
        fprintf(1, ['Pairwise auROC analysis for different shift times for ' evTypeLabels{1} ' and ' evTypeLabels{2} ' for ERP LFP power\n\n'],'perCorr_pre','perCorr_post')
        
        shift_ii=[1:4:41];
        
        delta_meanauROC=[];
        delta_meanauROCM=[];
        CIauROCpreLowM=[];
        CIauROCpreUppM=[];
        CIauROCpostLowM=[];
        CIauROCpostUppM=[];
        
        CIauROCdeltaLow=[];
        CIauROCdeltaUpp=[];
        
        for ii_shift=1:length(shift_ii)
            
            fprintf(1, '\nDelta t from event %d of %d\n',ii_shift, length(shift_ii));
            
            
            no_dBs=1;
            delta_dB_power_pre=[];
            no_ROCs=0;
            ROCoutpre=[];
            ROCoutpost=[];
            p_vals_ROC=[];
            delta_dB_powerpreHit=[];
            no_hits=0;
            perCorr_pre=[];
            perCorr_post=[];
            group_pre=[];
            group_post=[];
            sz_fps=size(file_pairs);
            no_file_pairs=sz_fps(1);
            
            p_vals=[];
            for fps=1:no_file_pairs
                for elec=1:16
                    
                    lfpodNopre=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                    lfpodNopost=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                    
                    if elec==1
                        %Find percent correct for pre
                        trials_Sp=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(event1,:)==1);
                        trials_Sm=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(event2,:)==1);
                        trials_Hit=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(perevent1,:)==1);
                        trials_CR=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(perevent2,:)==1);
                        perCorr_pre(fps)=100*(trials_Hit+trials_CR)/(trials_Sp+trials_Sm);
                        group_pre(fps)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopre).fileNo);
                        
                        %Find percent correct for post
                        trials_Sp=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(event1,:)==1);
                        trials_Sm=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(event2,:)==1);
                        trials_Hit=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(perevent1,:)==1);
                        trials_CR=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(perevent2,:)==1);
                        perCorr_post(fps)=100*(trials_Hit+trials_CR)/(trials_Sp+trials_Sm);
                        group_post(fps)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopost).fileNo);
                        
                    end
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre)))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost)))
                        
                        
                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre).log_P_tERP))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost).log_P_tERP))
                            
                            if (length(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(1,:))>=trials_to_process) &...
                                    (length(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(1,:))>=trials_to_process)
                                
                                length_pre=length(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(1,:));
                                pre_mask=logical([zeros(1,length_pre-trials_to_process) ones(1,trials_to_process)]);
                                trials_in_event_preHit=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(event1,:)==1);
                                trials_in_event_preCR=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(event2,:)==1);
                                
                                trials_with_event_pre=(handles_drgb.drgb.lfpevpair(lfpodNopre).no_events_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNopre).no_events_per_trial<=max_events_per_sec)...
                                    &(handles_drgb.drgb.lfpevpair(lfpodNopre).no_ref_evs_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNopre).no_ref_evs_per_trial<=max_events_per_sec);
                                
                                
                                length_post=length(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(1,:));
                                post_mask=logical([ones(1,trials_to_process) zeros(1,length_post-trials_to_process)]);
                                trials_in_event_postHit=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(event1,:)==1);
                                trials_in_event_postCR=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(event2,:)==1);
                                
                                trials_with_event_post=(handles_drgb.drgb.lfpevpair(lfpodNopost).no_events_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNopost).no_events_per_trial<=max_events_per_sec)...
                                    &(handles_drgb.drgb.lfpevpair(lfpodNopost).no_ref_evs_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNopost).no_ref_evs_per_trial<=max_events_per_sec);
                                
                                
                                if (sum(trials_in_event_preHit&trials_with_event_pre&pre_mask)>=min_trials_per_event) & (sum(trials_in_event_preCR&trials_with_event_pre&pre_mask)>=min_trials_per_event) & ...
                                        (sum(trials_in_event_postHit&trials_with_event_post&post_mask)>=min_trials_per_event) & (sum(trials_in_event_postCR&trials_with_event_post&post_mask)>=min_trials_per_event)
                                    
                                    
                                    %pre Hits
                                    this_dB_powerpreHit=zeros(sum(trials_in_event_preHit&pre_mask),length(frequency));
                                    this_dB_powerpreHit(:,:)=handles_drgb.drgb.lfpevpair(lfpodNopre).log_P_tERP(trials_in_event_preHit&pre_mask,:,shift_ii(ii_shift));
                                    
                                    
                                    %pre CRs
                                    this_dB_powerpreCR=zeros(sum(trials_in_event_preCR&pre_mask),length(frequency));
                                    this_dB_powerpreCR(:,:)=handles_drgb.drgb.lfpevpair(lfpodNopre).log_P_tERP(trials_in_event_preCR&pre_mask,:,shift_ii(ii_shift));
                                    
                                    
                                    %post Hits
                                    this_dB_powerpostHit=zeros(sum(trials_in_event_postHit&post_mask),length(frequency));
                                    this_dB_powerpostHit(:,:)=handles_drgb.drgb.lfpevpair(lfpodNopost).log_P_tERP(trials_in_event_postHit&post_mask,:,shift_ii(ii_shift));
                                    
                                    
                                    %post CRs
                                    this_dB_powerpostCR=zeros(sum(trials_in_event_postCR&post_mask),length(frequency));
                                    this_dB_powerpostCR(:,:)=handles_drgb.drgb.lfpevpair(lfpodNopost).log_P_tERP(trials_in_event_postCR&post_mask,:,shift_ii(ii_shift));
                                    
                                    for bwii=1:no_bandwidths
                                        
                                        no_ROCs=no_ROCs+1;
                                        this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                        
                                        %Enter the pre Hits
                                        this_delta_dB_powerpreHit=zeros(sum(trials_in_event_preHit&pre_mask),1);
                                        this_delta_dB_powerpreHit=mean(this_dB_powerpreHit(:,this_band),2);
                                        roc_data=[];
                                        roc_data(1:sum(trials_in_event_preHit&pre_mask),1)=this_delta_dB_powerpreHit;
                                        roc_data(1:sum(trials_in_event_preHit&pre_mask),2)=zeros(sum(trials_in_event_preHit&pre_mask),1);
                                        
                                        %Enter pre CR
                                        total_trials=sum(trials_in_event_preHit&pre_mask)+sum(trials_in_event_preCR&pre_mask);
                                        this_delta_dB_powerpreCR=zeros(sum(trials_in_event_preCR&pre_mask),1);
                                        this_delta_dB_powerpreCR=mean(this_dB_powerpreCR(:,this_band),2);
                                        roc_data(sum(trials_in_event_preHit&pre_mask)+1:total_trials,1)=this_delta_dB_powerpreCR;
                                        roc_data(sum(trials_in_event_preHit&pre_mask)+1:total_trials,2)=ones(sum(trials_in_event_preCR&pre_mask),1);
                                        
                                        
                                        %Find pre ROC
                                        ROCoutpre(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                        ROCoutpre(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopre).fileNo;
                                        ROCgroupNopre(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopre).fileNo);
                                        ROCoutpre(no_ROCs).timeWindow=winNo;
                                        ROCbandwidthpre(no_ROCs)=bwii;
                                        auROCpre(no_ROCs)=ROCoutpre(no_ROCs).roc.AUC-0.5;
                                        p_valROCpre(no_ROCs)=ROCoutpre(no_ROCs).roc.p;
                                        
                                        p_vals_ROC=[p_vals_ROC ROCoutpre(no_ROCs).roc.p];
                                        
                                        %Enter the post Hits
                                        this_delta_dB_powerpostHit=zeros(sum(trials_in_event_postHit&post_mask),1);
                                        this_delta_dB_powerpostHit=mean(this_dB_powerpostHit(:,this_band),2);
                                        roc_data=[];
                                        roc_data(1:sum(trials_in_event_postHit&post_mask),1)=this_delta_dB_powerpostHit;
                                        roc_data(1:sum(trials_in_event_postHit&post_mask),2)=zeros(sum(trials_in_event_postHit&post_mask),1);
                                        
                                        %Enter post CR
                                        total_trials=sum(trials_in_event_postHit&post_mask)+sum(trials_in_event_postCR&post_mask);
                                        this_delta_dB_powerpostCR=zeros(sum(trials_in_event_postCR&post_mask),1);
                                        this_delta_dB_powerpostCR=mean(this_dB_powerpostCR(:,this_band),2);
                                        roc_data(sum(trials_in_event_postHit&post_mask)+1:total_trials,1)=this_delta_dB_powerpostCR;
                                        roc_data(sum(trials_in_event_postHit&post_mask)+1:total_trials,2)=ones(sum(trials_in_event_postCR&post_mask),1);
                                        
                                        
                                        %Find post ROC
                                        ROCoutpost(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                        ROCoutpost(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopost).fileNo;
                                        ROCgroupNopost(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopost).fileNo);
                                        ROCoutpost(no_ROCs).timeWindow=winNo;
                                        ROCbandwidthpost(no_ROCs)=bwii;
                                        auROCpost(no_ROCs)=ROCoutpost(no_ROCs).roc.AUC-0.5;
                                        p_valROCpost(no_ROCs)=ROCoutpost(no_ROCs).roc.p;
                                        
                                        p_vals_ROC=[p_vals_ROC ROCoutpost(no_ROCs).roc.p];
                                        
                                        
                                        %Are the delta dB LFP's different?
                                        
                                        %Hit
                                        p_val(no_dBs,bwii)=ranksum(this_delta_dB_powerpreHit,this_delta_dB_powerpostHit);
                                        p_vals=[p_vals p_val(no_dBs,bwii)];
                                        groupNopre(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                        groupNopost(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                        events(no_dBs)=1;
                                        
                                        
                                        %CR
                                        p_val(no_dBs+1,bwii)=ranksum(this_delta_dB_powerpreCR,this_delta_dB_powerpostCR);
                                        p_vals=[p_vals p_val(no_dBs+1,bwii)];
                                        groupNopre(no_dBs+1)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                        groupNopost(no_dBs+1)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                        events(no_dBs+1)=2;
                                        
                                        if p_val(no_dBs,bwii)<0.05
                                            dB_power_changeHit(no_ROCs)=1;
                                        else
                                            dB_power_changeHit(no_ROCs)=0;
                                        end
                                        
                                        if p_val(no_dBs+1,bwii)<0.05
                                            dB_power_changeCR(no_ROCs)=1;
                                        else
                                            dB_power_changeCR(no_ROCs)=0;
                                        end
                                        
                                        %Plot the points and save the data
                                        if groupNopre(no_dBs)==1
                                            
                                            
                                            %Hit, all points
                                            delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                            delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                            
                                            
                                            %CR, all points
                                            delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                            delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                            
                                        else
                                            if groupNopre(no_dBs)==3
                                                
                                                %Hit, all points
                                                delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                                delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                                
                                                
                                                %CR, all points
                                                delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                                delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                                %                                             figure(bwii+4+12)
                                                %                                             hold on
                                                %                                             plot([3 4],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            end
                                        end
                                    end
                                    
                                    no_dBs=no_dBs+2;
                                    
                                    
                                end
                                
                                
                                
                            end
                            
                            
                            
                        end
                        
                        
                        
                    end
                end
                
            end
            
            
            
            pFDRROC=drsFDRpval(p_vals_ROC);
            
            %Now plot the bar graphs and do anovan for LFP power
            p_vals_anovan=[];
            
            
            %Plot cumulative histos for auROCs
            dB_power_change=logical(dB_power_changeHit+dB_power_changeCR);
            figNo=0;
            p_val_ROC=[];
            pvals_auROCperm=[];
            
            
            for bwii=1:4
                n_cum=0;
                this_legend=[];
                data_auROC=[];
                pre_post_auROC=[];
                gr_auROC=[];
                for grs=1:2
                    
                    this_deltaauROC=[];
                    this_deltaauROC=auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))-auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii));
                    delta_meanauROC(ii_shift,grs,bwii)=mean(this_deltaauROC);
                    CI= bootci(1000, @mean, this_deltaauROC);
                    CIauROCdeltaLow(ii_shift,grs,bwii)=mean(this_deltaauROC)-CI(1);
                    CIauROCdeltaUpp(ii_shift,grs,bwii)=CI(2)-mean(this_deltaauROC);
                    
                    
                    meanauROCpre(ii_shift,grs,bwii)=mean(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                    CI= bootci(1000, @mean, auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                    CIauROCpreLowM(ii_shift,grs,bwii)=mean(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))-CI(1);
                    CIauROCpreUppM(ii_shift,grs,bwii)=CI(2)-mean(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                    
                    
                    meanauROCpost(ii_shift,grs,bwii)=mean(auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                    CI = bootci(1000, @mean, auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                    CIauROCpostLowM(ii_shift,grs,bwii)=mean(auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))-CI(1);
                    CIauROCpostUppM(ii_shift,grs,bwii)=CI(2)-mean(auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                    
                    delta_meanauROCM(ii_shift,grs,bwii)=meanauROCpost(ii_shift,grs,bwii)-meanauROCpre(ii_shift,grs,bwii);
                    
                    
                    %Save the data for anovan interaction
                    %Pre
                    data_auROC=[data_auROC auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))];
                    gr_auROC=[gr_auROC grs*ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
                    pre_post_auROC=[pre_post_auROC ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
                    
                    %Post
                    data_auROC=[data_auROC auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))];
                    gr_auROC=[gr_auROC grs*ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
                    pre_post_auROC=[pre_post_auROC 2*ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
                end
                %             figNo=figNo+2;
                %             x=x+3;
                
                %Calculate anovan for inteaction
                [p,tbl,stats]=anovan(data_auROC,{pre_post_auROC gr_auROC},'model','interaction','varnames',{'pre_vs_post','halo_vs_no_halo'},'display','off');
                fprintf(1, ['p value for anovan auROC interaction for ii_shift= %d ' freq_names{bwii} '= %d\n'],  ii_shift,p(3));
                p_aovan_int(ii_shift,bwii)=p(3);
                
            end
            
        end
        
        pFDRauROCint=drsFDRpval(p_aovan_int(:));
        fprintf(1, ['pFDR for auROC anovan interaction  = %d\n\n'],pFDRauROCint);
        
        %Plot the effect of silencing NA fibers with halorhodopsin
        figure(1)
        hold on 
         plot([-0.6 0.6],[0 0],'-','LineWidth',3,'Color',[0.7 0.7 0.7])
        plot([0 0],[-0.25 0.25],'-','LineWidth',3,'Color',[0.7 0.7 0.7])
        delta_time=([0:10]*0.1-0.5);
        for bwii=1:4
            delta_auROC=zeros(1,length(shift_ii));
            delta_auROC=delta_meanauROC(:,1,bwii)-delta_meanauROC(:,2,bwii);
            delta_CIlow=sqrt(CIauROCdeltaLow(:,1,bwii).^2+CIauROCdeltaLow(:,2,bwii).^2);
            delta_CIupp=sqrt(CIauROCdeltaUpp(:,1,bwii).^2+CIauROCdeltaUpp(:,2,bwii).^2);
           
            for ii=1:length(delta_time)
                plot([delta_time(ii)+0.02*(bwii-2.5) delta_time(ii)+0.02*(bwii-2.5)],[delta_auROC(ii) delta_auROC(ii)+delta_CIupp(ii)],these_lines{bwii},'LineWidth',1)
                plot([delta_time(ii)+0.02*(bwii-2.5) delta_time(ii)+0.02*(bwii-2.5)],[delta_auROC(ii) delta_auROC(ii)-delta_CIlow(ii)],these_lines{bwii},'LineWidth',1)
            end
            plot(delta_time+0.02*(bwii-2.5),delta_auROC,'-','LineWidth',2,'Color',these_colors{bwii})
            plot(delta_time(p_aovan_int(:,bwii)<=pFDRauROCint)+0.02*(bwii-2.5),delta_auROC(p_aovan_int(:,bwii)<pFDRauROCint),'*','Color',these_colors{bwii},'MarkerSize',10)
            plot(delta_time(p_aovan_int(:,bwii)>pFDRauROCint)+0.02*(bwii-2.5),delta_auROC(p_aovan_int(:,bwii)>pFDRauROCint),'o','Color',these_colors{bwii},'MarkerFace',these_colors{bwii})
            
        end
       
        ylim([-0.25 0.25])
        xlim([-0.6 0.6])
        xlabel('dt to event (ms)')
        ylabel('delta auROC')
        title('Effect of silencing NA on auROC')
        legend('Theta','Beta','Low gamma','High gamma')
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2, 'Box', 'off')
        
case 9
    %Compare LFP power auROC in two percent windows for all of the files 
    no_dBs=1;
    delta_dB_power_fp1=[];
    no_ROCs=0;
    ROCoutfp1=[];
    ROCoutfp2=[];
    p_vals_ROC=[];
    delta_dB_powerfp1Ev1=[];
    no_Ev1=0;
    pvals_auROCperm=[];
    pvals_dBperm=[];
    perCorr_fp1=[];
    perCorr_fp2=[];
    
    fprintf(1, ['Pairwise auROC analysis for ' evTypeLabels{1} ' and ' evTypeLabels{2} ' LFP power\n\n'])
    p_vals=[];
    
    if exist('which_electrodes')==0
        which_electrodes=[1:16];
    end
    
    no_files=length(files);
    
    
    for fileNo=1:no_files
        
        
        for elec=1:16
            if sum(which_electrodes==elec)>0
                
                
                lfpodNo_ref=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                
                if ~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref))
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower))
                        
                        for per_ii=1:2
                            
                            percent_mask=[];
                            trials_in_event_Ev1=[];
                            trials_in_event_Ev2=[];
                            percent_mask=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower>=percent_windows(per_ii,1))&(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower<=percent_windows(per_ii,2));
                            trials_in_event_Ev1=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(event1,:)==1)&percent_mask;
                            trials_in_event_Ev2=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(event2,:)==1)&percent_mask;
                            
<<<<<<< HEAD
                         
=======
                            if (files(fileNo)==7)&(elec==5)
                                pffft=1;
                            end
>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855
                            if (sum(trials_in_event_Ev1)>=min_trials_per_event) & (sum( trials_in_event_Ev2)>=min_trials_per_event)
                                
                                fprintf(1, ['File no %d electrode %d for percent window No %d was processed succesfully\n'],files(fileNo),elec,per_ii)
                                lfpodNo=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                
                                %Ev1
                                this_dB_powerrefEv1=zeros(sum(trials_in_event_Ev1),length(frequency));
                                this_dB_powerrefEv1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower(trials_in_event_Ev1,:));
                                
                                this_dB_powerEv1=zeros(sum(trials_in_event_Ev1),length(frequency));
                                this_dB_powerEv1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_event_Ev1,:));
                                
                                %Ev2
                                this_dB_powerrefEv2=zeros(sum(trials_in_event_Ev2),length(frequency));
                                this_dB_powerrefEv2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower(trials_in_event_Ev2,:));
                                
                                this_dB_powerEv2=zeros(sum(trials_in_event_Ev2),length(frequency));
                                this_dB_powerEv2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_event_Ev2,:));
                                
                                for bwii=1:no_bandwidths
                                    
                                    no_ROCs=no_ROCs+1;
                                    this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                    
                                    %Enter Ev1
                                    this_delta_dB_powerEv1=zeros(sum(trials_in_event_Ev1),1);
                                    this_delta_dB_powerEv1=mean(this_dB_powerEv1(:,this_band)-this_dB_powerrefEv1(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:sum(trials_in_event_Ev1),1)=this_delta_dB_powerEv1;
                                    roc_data(1:sum(trials_in_event_Ev1),2)=zeros(sum(trials_in_event_Ev1),1);
                                    
                                    %Enter Ev2
                                    total_trials=sum(trials_in_event_Ev1)+sum(trials_in_event_Ev2);
                                    this_delta_dB_powerEv2=zeros(sum(trials_in_event_Ev2),1);
                                    this_delta_dB_powerEv2=mean(this_dB_powerEv2(:,this_band)-this_dB_powerrefEv2(:,this_band),2);
                                    roc_data(sum(trials_in_event_Ev1)+1:total_trials,1)=this_delta_dB_powerEv2;
                                    roc_data(sum(trials_in_event_Ev1)+1:total_trials,2)=ones(sum(trials_in_event_Ev2),1);
                                    
                                    
                                    %Find  ROC
                                    ROCout(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCout(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo_ref).fileNo;
                                    ROCgroupNo(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNo_ref).fileNo);
                                    ROCout(no_ROCs).timeWindow=winNo;
                                    ROCbandwidth(no_ROCs)=bwii;
                                    ROCper_ii(no_ROCs)=per_ii;
                                    auROC(no_ROCs)=ROCout(no_ROCs).roc.AUC-0.5;
                                    p_valROC(no_ROCs)=ROCout(no_ROCs).roc.p;
                                    
                                    p_vals_ROC=[p_vals_ROC ROCout(no_ROCs).roc.p];
                                    
                                    
                                    
                                end
                                
                                
                            else
                                
                                if (sum(trials_in_event_Ev1)<min_trials_per_event)
                                    fprintf(1, ['%d trials for ' evTypeLabels{1} ' fewer than minimum trials per event =%d for file No %d electrode %d\n'],sum(trials_in_event_Ev1), min_trials_per_event,files(fileNo),elec);
                                end
                                
                                if (sum(trials_in_event_Ev2)<min_trials_per_event)
                                    fprintf(1, ['%d trials for ' evTypeLabels{2} ' fewer than minimum trials per event =%d for file No %d electrode %d\n'],sum(trials_in_event_Ev2), min_trials_per_event,files(fileNo),elec);
                                end
                                
                                
                            end
                            
                        end
                        
                    else
                        fprintf(1, ['Empty allPower for file No %d electrode %d\n'],files(fileNo),elec);
                    end
                    
                else
                    fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],files(fileNo),elec);
                end
            end
        end
        
    end
    fprintf(1, '\n\n')
        
        pFDRROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for significant difference of auROC p value from 0.5  = %d\n\n'],pFDRROC);
        
        
        
        fprintf(1, '\n\n')
        
        
        %Plot cumulative histos for auROCs

        figNo=0;
        x=0;
        
        try
            close(5)
        catch
        end
        figure(5)
        hold on
        
        
        for bwii=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            %Plot the histograms
            edges=[-0.5:0.05:0.5];
            pos2=[0.1 0.1 0.6 0.8];
            subplot('Position',pos2)
            hold on
            
            h2=histogram(auROC((ROCbandwidth==bwii)&(ROCper_ii==1)),edges);
            h2.FaceColor='r';
            h1=histogram(auROC((ROCbandwidth==bwii)&(ROCper_ii==2)),edges);
            h1.FaceColor='b';
            
            xlabel('auROC')
            ylabel('# of electrodes')
            legend(file_label{2},file_label{1})
            title(['auROC for ' freq_names{bwii}])
            xlim([-0.3 0.6])
            ylim([0 80])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Plot the single electrodes
            pos2=[0.8 0.1 0.1 0.8];
            subplot('Position',pos2)
            hold on
            
            plot(ones(1,sum((ROCbandwidth==bwii)&(ROCper_ii==1))),auROC((ROCbandwidth==bwii)&(ROCper_ii==1)),'o', 'Color',[0.7 0.7 0.7])
            plot(zeros(1,sum((ROCbandwidth==bwii)&(ROCper_ii==2))),auROC((ROCbandwidth==bwii)&(ROCper_ii==2)),'o', 'Color',[0.7 0.7 0.7])
       
            
            %PLot the mean and 95% CI
            plot([0 1],[mean(auROC((ROCbandwidth==bwii)&(ROCper_ii==2))) mean(auROC((ROCbandwidth==bwii)&(ROCper_ii==1)))],'-k','LineWidth', 3)
            CI = bootci(1000, @mean, auROC((ROCbandwidth==bwii)&(ROCper_ii==2)));
            plot([0 0],CI,'-b','LineWidth',3)
            plot(0,mean(auROC((ROCbandwidth==bwii)&(ROCper_ii==2))),'ob','MarkerSize', 10,'MarkerFace','b')
            CI = bootci(1000, @mean, auROC((ROCbandwidth==bwii)&(ROCper_ii==1)));
            plot([1 1],CI,'-r','LineWidth',3)
            plot(1,mean(auROC((ROCbandwidth==bwii)&(ROCper_ii==1))),'or','MarkerSize', 10,'MarkerFace','r')
            ylabel('auROC')
            ylim([-0.2 0.5])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Do the statistics for auROC differences
            a={auROC((ROCbandwidth==bwii)&(ROCper_ii==1)) auROC((ROCbandwidth==bwii)&(ROCper_ii==2))};
            mode_statcond='perm';
            [F df pval_auROCperm] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for permuted anovan for auROC S+ vs S- ' freq_names{bwii} '= %d\n\n'],  pval_auROCperm);
            pvals_auROCperm=[pvals_auROCperm pval_auROCperm];
            
            %Figure 5
            figure(5)
            
            percent_auROCper2=100*sum(p_valROC((ROCbandwidth==bwii)&(ROCper_ii==2))<=pFDRROC)/sum((ROCbandwidth==bwii)&(ROCper_ii==2));
            bar(x,percent_auROCper2,'b')
            
            learn_sig(bwii)=sum(p_valROC((ROCbandwidth==bwii)&(ROCper_ii==2))<=pFDRROC);
            learn_not_sig(bwii)=sum((ROCbandwidth==bwii)&(ROCper_ii==2))-sum(p_valROC((ROCbandwidth==bwii)&(ROCper_ii==2))<=pFDRROC);
            
            percent_auROCper1=100*sum(p_valROC((ROCbandwidth==bwii)&(ROCper_ii==1))<=pFDRROC)/sum((ROCbandwidth==bwii)&(ROCper_ii==1));
            bar(x+1,percent_auROCper1,'r')
            
            prof_sig(bwii)=sum(p_valROC((ROCbandwidth==bwii)&(ROCper_ii==1))<=pFDRROC);
            prof_not_sig(bwii)=sum((ROCbandwidth==bwii)&(ROCper_ii==1))-sum(p_valROC((ROCbandwidth==bwii)&(ROCper_ii==1))<=pFDRROC);
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            x=x+3;
            
        end
        
        figure(5)
        title('Percent singificant auROC')
        legend(file_label{2},file_label{1})
        ylim([0 100])
        
        pFDRanovanauROC=drsFDRpval(pval_auROCperm);
        fprintf(1, ['\npFDR for premuted anovan p value for difference between ' file_label{1} ' and ' file_label{2} ' for auROC = %d\n\n'],pFDRanovanauROC);
        
        
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix],'learn_sig','learn_not_sig','prof_sig','prof_not_sig');
        
        
        pffft=1
        
        
    case 10
        %Compare LFP power auROC for two groups (i.e. NRG1 vs control) for a percent window
        no_dBs=1;
        delta_dB_power_fp1=[];
        no_ROCs=0;
        ROCoutfp1=[];
        ROCoutfp2=[];
        p_vals_ROC=[];
        delta_dB_powerfp1Ev1=[];
        no_Ev1=0;
        pvals_auROCperm=[];
        pvals_dBperm=[];
        perCorr_fp1=[];
        perCorr_fp2=[];
        
        fprintf(1, ['Pairwise auROC analysis for ' evTypeLabels{1} ' and ' evTypeLabels{2} ' LFP power\n\n'])
        p_vals=[];
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        no_files=length(files);
        
        
        for fileNo=1:no_files
            
            
            for elec=1:16
                if sum(which_electrodes==elec)>0
                    
                    
                    lfpodNo_ref=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                    
                    if ~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref))
                        
                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower))
                            
                            per_ii=1;
                            
                            percent_mask=[];
                            trials_in_event_Ev1=[];
                            trials_in_event_Ev2=[];
                            percent_mask=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower>=percent_windows(per_ii,1))...
                                &(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower<=percent_windows(per_ii,2));
                            trials_in_event_Ev1=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(event1,:)==1)&percent_mask;
                            trials_in_event_Ev2=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(event2,:)==1)&percent_mask;
                            if (files(fileNo)==7)&(elec==5)
                                pffft=1
                            end
                            if (sum(trials_in_event_Ev1)>=min_trials_per_event) & (sum( trials_in_event_Ev2)>=min_trials_per_event)
                                
                                fprintf(1, ['File no %d electrode %d for percent window No %d was processed succesfully\n'],files(fileNo),elec,per_ii)
                                lfpodNo=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                
                                %Ev1
                                this_dB_powerrefEv1=zeros(sum(trials_in_event_Ev1),length(frequency));
                                this_dB_powerrefEv1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower(trials_in_event_Ev1,:));
                                
                                this_dB_powerEv1=zeros(sum(trials_in_event_Ev1),length(frequency));
                                this_dB_powerEv1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_event_Ev1,:));
                                
                                %Ev2
                                this_dB_powerrefEv2=zeros(sum(trials_in_event_Ev2),length(frequency));
                                this_dB_powerrefEv2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower(trials_in_event_Ev2,:));
                                
                                this_dB_powerEv2=zeros(sum(trials_in_event_Ev2),length(frequency));
                                this_dB_powerEv2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_event_Ev2,:));
                                
                                for bwii=1:no_bandwidths
                                    
                                    no_ROCs=no_ROCs+1;
                                    this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                    
                                    %Enter Ev1
                                    this_delta_dB_powerEv1=zeros(sum(trials_in_event_Ev1),1);
                                    this_delta_dB_powerEv1=mean(this_dB_powerEv1(:,this_band)-this_dB_powerrefEv1(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:sum(trials_in_event_Ev1),1)=this_delta_dB_powerEv1;
                                    roc_data(1:sum(trials_in_event_Ev1),2)=zeros(sum(trials_in_event_Ev1),1);
                                    
                                    %Enter Ev2
                                    total_trials=sum(trials_in_event_Ev1)+sum(trials_in_event_Ev2);
                                    this_delta_dB_powerEv2=zeros(sum(trials_in_event_Ev2),1);
                                    this_delta_dB_powerEv2=mean(this_dB_powerEv2(:,this_band)-this_dB_powerrefEv2(:,this_band),2);
                                    roc_data(sum(trials_in_event_Ev1)+1:total_trials,1)=this_delta_dB_powerEv2;
                                    roc_data(sum(trials_in_event_Ev1)+1:total_trials,2)=ones(sum(trials_in_event_Ev2),1);
                                    
                                    
                                    %Find  ROC
                                    ROCout(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCout(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo_ref).fileNo;
                                    ROCgroupNo(no_ROCs)=groups(fileNo);
                                    ROCout(no_ROCs).timeWindow=winNo;
                                    ROCbandwidth(no_ROCs)=bwii;
                                    ROCper_ii(no_ROCs)=per_ii;
                                    auROC(no_ROCs)=ROCout(no_ROCs).roc.AUC-0.5;
                                    p_valROC(no_ROCs)=ROCout(no_ROCs).roc.p;
                                    
                                    p_vals_ROC=[p_vals_ROC ROCout(no_ROCs).roc.p];
                                    
                                    
                                    
                                end
                                
                                
                            else
                                
                                if (sum(trials_in_event_Ev1)<min_trials_per_event)
                                    fprintf(1, ['%d trials for ' evTypeLabels{1} ' fewer than minimum trials per event =%d for file No %d electrode %d\n'],sum(trials_in_event_Ev1), min_trials_per_event,fileNo,elec);
                                end
                                
                                if (sum(trials_in_event_Ev2)<min_trials_per_event)
                                    fprintf(1, ['%d trials for ' evTypeLabels{2} ' fewer than minimum trials per event =%d for file No %d electrode %d\n'],sum(trials_in_event_Ev2), min_trials_per_event,fileNo,elec);
                                end
                                
                                
                            end
                            
                            
                            
                        else
                            fprintf(1, ['Empty allPower for file No %d electrode %d\n'],fileNo,elec);
                        end
                        
                    else
                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],fileNo,elec);
                    end
                end
            end
            
        end
        fprintf(1, '\n\n')
        
        pFDRROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for significant difference of auROC p value from 0.5  = %d\n\n'],pFDRROC);
        
        
        
        fprintf(1, '\n\n')
        
        
        %Plot cumulative histos for auROCs
        
        figNo=0;
        x=0;
        
        try
            close(5)
        catch
        end
        figure(5)
        hold on
        
        
        for bwii=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            %Plot the histograms
            edges=[-0.5:0.05:0.5];
            pos2=[0.1 0.1 0.6 0.8];
            subplot('Position',pos2)
            hold on
            
            h2=histogram(auROC((ROCbandwidth==bwii)&(ROCgroupNo==1)),edges);
            h2.FaceColor='r';
            h1=histogram(auROC((ROCbandwidth==bwii)&(ROCgroupNo==2)),edges);
            h1.FaceColor='b';
            
            
            xlabel('auROC')
            ylabel('# of electrodes')
            legend(group_names{1},group_names{2})
            title(['auROC for ' freq_names{bwii}])
            xlim([-0.3 0.6])
            ylim([0 80])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Plot the single electrodes
            pos2=[0.8 0.1 0.1 0.8];
            subplot('Position',pos2)
            hold on
            
            plot(ones(1,sum((ROCbandwidth==bwii)&(ROCgroupNo==1))),auROC((ROCbandwidth==bwii)&(ROCgroupNo==1)),'o', 'Color',[0.7 0.7 0.7])
            plot(zeros(1,sum((ROCbandwidth==bwii)&(ROCgroupNo==2))),auROC((ROCbandwidth==bwii)&(ROCgroupNo==2)),'o', 'Color',[0.7 0.7 0.7])
            
            
            %PLot the mean and 95% CI
            plot([0 1],[mean(auROC((ROCbandwidth==bwii)&(ROCgroupNo==2))) mean(auROC((ROCbandwidth==bwii)&(ROCgroupNo==1)))],'-k','LineWidth', 3)
            CI = bootci(1000, @mean, auROC((ROCbandwidth==bwii)&(ROCgroupNo==2)));
            plot([0 0],CI,'-b','LineWidth',3)
            plot(0,mean(auROC((ROCbandwidth==bwii)&(ROCgroupNo==2))),'ob','MarkerSize', 10,'MarkerFace','b')
            CI = bootci(1000, @mean, auROC((ROCbandwidth==bwii)&(ROCgroupNo==1)));
            plot([1 1],CI,'-r','LineWidth',3)
            plot(1,mean(auROC((ROCbandwidth==bwii)&(ROCgroupNo==1))),'or','MarkerSize', 10,'MarkerFace','r')
            ylabel('auROC')
            ylim([-0.2 0.5])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Do the statistics for auROC differences
            a={auROC((ROCbandwidth==bwii)&(ROCgroupNo==1)) auROC((ROCbandwidth==bwii)&(ROCgroupNo==2))};
            mode_statcond='perm';
            [F df pval_auROCperm] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for permuted anovan for auROC S+ vs S- ' freq_names{bwii} '= %d\n\n'],  pval_auROCperm);
            pvals_auROCperm=[pvals_auROCperm pval_auROCperm];
            
            %Figure 5
            figure(5)
            
            percent_auROCper1=100*sum(p_valROC((ROCbandwidth==bwii)&(ROCgroupNo==1))<=pFDRROC)/sum((ROCbandwidth==bwii)&(ROCgroupNo==1));
            bar(x+1,percent_auROCper1,'r')
            
            percent_auROCper2=100*sum(p_valROC((ROCbandwidth==bwii)&(ROCgroupNo==2))<=pFDRROC)/sum((ROCbandwidth==bwii)&(ROCgroupNo==2));
            bar(x,percent_auROCper2,'b')
            
            learn_sig(bwii)=sum(p_valROC((ROCbandwidth==bwii)&(ROCgroupNo==2))<=pFDRROC);
            learn_not_sig(bwii)=sum((ROCbandwidth==bwii)&(ROCgroupNo==2))-sum(p_valROC((ROCbandwidth==bwii)&(ROCgroupNo==2))<=pFDRROC);
            
            prof_sig(bwii)=sum(p_valROC((ROCbandwidth==bwii)&(ROCgroupNo==1))<=pFDRROC);
            prof_not_sig(bwii)=sum((ROCbandwidth==bwii)&(ROCgroupNo==1))-sum(p_valROC((ROCbandwidth==bwii)&(ROCgroupNo==1))<=pFDRROC);
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            x=x+3;
            
        end
        
        figure(5)
        title('Percent singificant auROC')
        legend(group_names{1},group_names{2})
        ylim([0 100])
        
        pFDRanovanauROC=drsFDRpval(pval_auROCperm);
        fprintf(1, ['\npFDR for premuted anovan p value for difference between\n ' group_names{1} ' and ' group_names{2} ' for auROC = %d\n\n'],pFDRanovanauROC);
        
        
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix],'learn_sig','learn_not_sig','prof_sig','prof_not_sig');
        
        
        pffft=1
        
    case 11
        
<<<<<<< HEAD
        %Compare auROC in the last few trials of the last session file with
        %first few trials of the first session
        %Generates Fig. 3 for Daniel's paper. first vs last.
        no_dBs=1;
        delta_dB_power_fp1=[];
        no_ROCs=0;
        ROCoutfp1=[];
        ROCoutfp2=[];
        p_vals_ROC=[];
        delta_dB_powerfp1Ev1=[];
        no_Ev1=0;
        pvals_auROCperm=[];
        pvals_dBperm=[];
        perCorr_fp1=[];
        perCorr_fp2=[];
        shift_ii=floor(length(handles_drgb.drgb.lfpevpair(1).out_times)/2)+1+shift_from_event;
        
        
        fprintf(1, ['Pairwise auROC log power LFP ERP analysis for ' evTypeLabels{1} ' and ' evTypeLabels{2} ' LFP power\n\n'])
        p_vals=[];
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        no_files=length(files);
        
        for fileNo=1:no_files
            
            
            for elec=1:16
                if sum(which_electrodes==elec)>0
                    lfpodNo=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                    
                    
                    if ~isempty(handles_drgb.drgb.lfpevpair(lfpodNo))
                        
                        
                        if ~isempty(handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP)
                            
                            
                            for per_ii=1:2
                                
                                percent_mask=[];
                                trials_in_event_Ev1=[];
                                trials_in_event_Ev2=[];
                                percent_mask=(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrERP>=percent_windows(per_ii,1))&...
                                    (handles_drgb.drgb.lfpevpair(lfpodNo).perCorrERP<=percent_windows(per_ii,2));
                                trials_in_event_Ev1=(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventERP(event1,:)==1)&percent_mask;
                                trials_in_event_Ev2=(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventERP(event2,:)==1)&percent_mask;
                                
                                
                                if (sum(trials_in_event_Ev1)>=min_trials_per_event) & (sum( trials_in_event_Ev2)>=min_trials_per_event)
                                    
                                    fprintf(1, ['File no %d electrode %d for percent window No %d was processed succesfully\n'],files(fileNo),elec,per_ii)
                                    
                                    
                                    %Ev1
                                    this_dB_powerEv1=zeros(sum(trials_in_event_Ev1),length(frequency));
                                    this_dB_powerEv1(:,:)=handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP(trials_in_event_Ev1,:,shift_ii);
                                    
                                    %Ev2
                                    this_dB_powerEv2=zeros(sum(trials_in_event_Ev2),length(frequency));
                                    this_dB_powerEv2(:,:)=handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP(trials_in_event_Ev2,:,shift_ii);
                                    
=======
        %Compare auROC for ERP LFP between two groups for percent >80
        %Under construction (Diego)
        
        fprintf(1, ['Analysis of auROC differences between two groups for different shift times for ' evTypeLabels{1} ' and ' evTypeLabels{2} ' for ERP LFP power\n\n'])
        
        shift_ii=[1:4:41];
        
        delta_meanauROC=[];
        delta_meanauROCM=[];
        CIauROCpreLowM=[];
        CIauROCpreUppM=[];
        CIauROCpostLowM=[];
        CIauROCpostUppM=[];
        
        CIauROCdeltaLow=[];
        CIauROCdeltaUpp=[];
        
        for ii_shift=1:length(shift_ii)
            
            fprintf(1, '\nDelta t from event %d of %d\n',ii_shift, length(shift_ii));
            
            
            no_dBs=1;
            delta_dB_power_pre=[];
            no_ROCs=0;
            ROCoutpre=[];
            ROCoutpost=[];
            p_vals_ROC=[];
            delta_dB_powerpreHit=[];
            no_hits=0;
            perCorr_pre=[];
            perCorr_post=[];
            group_pre=[];
            group_post=[];
            sz_fps=size(file_pairs);
            no_file_pairs=sz_fps(1);
            
            p_vals=[];
            for fps=1:no_file_pairs
                for elec=1:16
                    
                    lfpodNopre=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                    lfpodNopost=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                    
                    if elec==1
                        %Find percent correct for pre
                        trials_Sp=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(event1,:)==1);
                        trials_Sm=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(event2,:)==1);
                        trials_Hit=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(perevent1,:)==1);
                        trials_CR=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(perevent2,:)==1);
                        perCorr_pre(fps)=100*(trials_Hit+trials_CR)/(trials_Sp+trials_Sm);
                        group_pre(fps)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopre).fileNo);
                        
                        %Find percent correct for post
                        trials_Sp=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(event1,:)==1);
                        trials_Sm=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(event2,:)==1);
                        trials_Hit=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(perevent1,:)==1);
                        trials_CR=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(perevent2,:)==1);
                        perCorr_post(fps)=100*(trials_Hit+trials_CR)/(trials_Sp+trials_Sm);
                        group_post(fps)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopost).fileNo);
                        
                    end
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre)))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost)))
                        
                        
                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre).log_P_tERP))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost).log_P_tERP))
                            
                            if (length(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(1,:))>=trials_to_process) &...
                                    (length(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(1,:))>=trials_to_process)
                                
                                length_pre=length(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(1,:));
                                pre_mask=logical([zeros(1,length_pre-trials_to_process) ones(1,trials_to_process)]);
                                trials_in_event_preHit=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(event1,:)==1);
                                trials_in_event_preCR=(handles_drgb.drgb.lfpevpair(lfpodNopre).which_eventERP(event2,:)==1);
                                
                                trials_with_event_pre=(handles_drgb.drgb.lfpevpair(lfpodNopre).no_events_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNopre).no_events_per_trial<=max_events_per_sec)...
                                    &(handles_drgb.drgb.lfpevpair(lfpodNopre).no_ref_evs_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNopre).no_ref_evs_per_trial<=max_events_per_sec);
                                
                                
                                length_post=length(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(1,:));
                                post_mask=logical([ones(1,trials_to_process) zeros(1,length_post-trials_to_process)]);
                                trials_in_event_postHit=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(event1,:)==1);
                                trials_in_event_postCR=(handles_drgb.drgb.lfpevpair(lfpodNopost).which_eventERP(event2,:)==1);
                                
                                trials_with_event_post=(handles_drgb.drgb.lfpevpair(lfpodNopost).no_events_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNopost).no_events_per_trial<=max_events_per_sec)...
                                    &(handles_drgb.drgb.lfpevpair(lfpodNopost).no_ref_evs_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNopost).no_ref_evs_per_trial<=max_events_per_sec);
                                
                                
                                if (sum(trials_in_event_preHit&trials_with_event_pre&pre_mask)>=min_trials_per_event) & (sum(trials_in_event_preCR&trials_with_event_pre&pre_mask)>=min_trials_per_event) & ...
                                        (sum(trials_in_event_postHit&trials_with_event_post&post_mask)>=min_trials_per_event) & (sum(trials_in_event_postCR&trials_with_event_post&post_mask)>=min_trials_per_event)
                                    
                                    
                                    %pre Hits
                                    this_dB_powerpreHit=zeros(sum(trials_in_event_preHit&pre_mask),length(frequency));
                                    this_dB_powerpreHit(:,:)=handles_drgb.drgb.lfpevpair(lfpodNopre).log_P_tERP(trials_in_event_preHit&pre_mask,:,shift_ii(ii_shift));
                                    
                                    
                                    %pre CRs
                                    this_dB_powerpreCR=zeros(sum(trials_in_event_preCR&pre_mask),length(frequency));
                                    this_dB_powerpreCR(:,:)=handles_drgb.drgb.lfpevpair(lfpodNopre).log_P_tERP(trials_in_event_preCR&pre_mask,:,shift_ii(ii_shift));
                                    
                                    
                                    %post Hits
                                    this_dB_powerpostHit=zeros(sum(trials_in_event_postHit&post_mask),length(frequency));
                                    this_dB_powerpostHit(:,:)=handles_drgb.drgb.lfpevpair(lfpodNopost).log_P_tERP(trials_in_event_postHit&post_mask,:,shift_ii(ii_shift));
                                    
                                    
                                    %post CRs
                                    this_dB_powerpostCR=zeros(sum(trials_in_event_postCR&post_mask),length(frequency));
                                    this_dB_powerpostCR(:,:)=handles_drgb.drgb.lfpevpair(lfpodNopost).log_P_tERP(trials_in_event_postCR&post_mask,:,shift_ii(ii_shift));
>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855
                                    
                                    for bwii=1:no_bandwidths
                                        
                                        no_ROCs=no_ROCs+1;
                                        this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                        
<<<<<<< HEAD
                                        %Enter Ev1
                                        this_delta_dB_powerEv1=zeros(sum(trials_in_event_Ev1),1);
                                        this_delta_dB_powerEv1=mean(this_dB_powerEv1(:,this_band),2);
                                        roc_data=[];
                                        roc_data(1:sum(trials_in_event_Ev1),1)=this_delta_dB_powerfp1Ev1;
                                        roc_data(1:sum(trials_in_event_Ev1),2)=zeros(sum(trials_in_event_Ev1),1);
                                        
                                        %Enter Ev2
                                        total_trials=sum(trials_in_event_Ev2)+sum(trials_in_event_Ev1);
                                        this_delta_dB_powerEv2=zeros(sum(trials_in_event_Ev2),1);
                                        this_delta_dB_powerEv2=mean(this_dB_powerEv2(:,this_band),2);
                                        roc_data(sum(trials_in_event_Ev1)+1:total_trials,1)=this_delta_dB_powerEv2;
                                        roc_data(sum(trials_in_event_Ev1)+1:total_trials,2)=ones(sum(trials_in_event_Ev2),1);
                                        
                                        
                                        %Find ROC
                                        ROCout(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                        ROCout(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
                                        ROCgroupNofp1(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNo).fileNo);
                                        ROCout(no_ROCs).timeWindow=winNo;
                                        ROCbandwidthfp1(no_ROCs)=bwii;
                                        ROCper_ii(no_ROCs)=per_ii;
                                        auROC(no_ROCs)=ROCout(no_ROCs).roc.AUC-0.5;
                                        p_valROC(no_ROCs)=ROCout(no_ROCs).roc.p;
                                        
                                        p_vals_ROC=[p_vals_ROC ROCout(no_ROCs).roc.p];
                                        
                                        
                                        
                                        
                                        %Ev1, all points
                                        delta_dB_powerEv1(no_ROCs)=mean(this_delta_dB_powerEv1);
                                        delta_dB_powerEv1(no_ROCs)=mean(this_delta_dB_powerEv1);
                                        
                                        
                                    end
                                    
                                    no_dBs=no_dBs+2;
                                    
                                    
                                    
                                else
                                    
                                    if (sum(trials_in_event_Ev1)<min_trials_per_event)
                                        fprintf(1, ['%d trials for ' evTypeLabels{1} ' fewer than minimum trials per event =%d for file No %d electrode %d\n'],sum(trials_in_event_Ev1), min_trials_per_event,files(fileNo),elec);
                                    end
                                    
                                    if (sum(trials_in_event_Ev2)<min_trials_per_event)
                                        fprintf(1, ['%d trials for ' evTypeLabels{2} ' fewer than minimum trials per event =%d for file No %d electrode %d\n'],sum(trials_in_event_Ev2), min_trials_per_event,files(fileNo),elec);
                                    end
                                    
                                end
                            end
                        else
                            
                            
                            fprintf(1, ['Empty allPower for file No %d electrode %d\n'],files(fileNo),elec);
                            
                            
                        end
                        
                    else
                        
                        
                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],files(fileNo),elec);
=======
                                        %Enter the pre Hits
                                        this_delta_dB_powerpreHit=zeros(sum(trials_in_event_preHit&pre_mask),1);
                                        this_delta_dB_powerpreHit=mean(this_dB_powerpreHit(:,this_band),2);
                                        roc_data=[];
                                        roc_data(1:sum(trials_in_event_preHit&pre_mask),1)=this_delta_dB_powerpreHit;
                                        roc_data(1:sum(trials_in_event_preHit&pre_mask),2)=zeros(sum(trials_in_event_preHit&pre_mask),1);
                                        
                                        %Enter pre CR
                                        total_trials=sum(trials_in_event_preHit&pre_mask)+sum(trials_in_event_preCR&pre_mask);
                                        this_delta_dB_powerpreCR=zeros(sum(trials_in_event_preCR&pre_mask),1);
                                        this_delta_dB_powerpreCR=mean(this_dB_powerpreCR(:,this_band),2);
                                        roc_data(sum(trials_in_event_preHit&pre_mask)+1:total_trials,1)=this_delta_dB_powerpreCR;
                                        roc_data(sum(trials_in_event_preHit&pre_mask)+1:total_trials,2)=ones(sum(trials_in_event_preCR&pre_mask),1);
                                        
                                        
                                        %Find pre ROC
                                        ROCoutpre(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                        ROCoutpre(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopre).fileNo;
                                        ROCgroupNopre(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopre).fileNo);
                                        ROCoutpre(no_ROCs).timeWindow=winNo;
                                        ROCbandwidthpre(no_ROCs)=bwii;
                                        auROCpre(no_ROCs)=ROCoutpre(no_ROCs).roc.AUC-0.5;
                                        p_valROCpre(no_ROCs)=ROCoutpre(no_ROCs).roc.p;
                                        
                                        p_vals_ROC=[p_vals_ROC ROCoutpre(no_ROCs).roc.p];
                                        
                                        %Enter the post Hits
                                        this_delta_dB_powerpostHit=zeros(sum(trials_in_event_postHit&post_mask),1);
                                        this_delta_dB_powerpostHit=mean(this_dB_powerpostHit(:,this_band),2);
                                        roc_data=[];
                                        roc_data(1:sum(trials_in_event_postHit&post_mask),1)=this_delta_dB_powerpostHit;
                                        roc_data(1:sum(trials_in_event_postHit&post_mask),2)=zeros(sum(trials_in_event_postHit&post_mask),1);
                                        
                                        %Enter post CR
                                        total_trials=sum(trials_in_event_postHit&post_mask)+sum(trials_in_event_postCR&post_mask);
                                        this_delta_dB_powerpostCR=zeros(sum(trials_in_event_postCR&post_mask),1);
                                        this_delta_dB_powerpostCR=mean(this_dB_powerpostCR(:,this_band),2);
                                        roc_data(sum(trials_in_event_postHit&post_mask)+1:total_trials,1)=this_delta_dB_powerpostCR;
                                        roc_data(sum(trials_in_event_postHit&post_mask)+1:total_trials,2)=ones(sum(trials_in_event_postCR&post_mask),1);
                                        
                                        
                                        %Find post ROC
                                        ROCoutpost(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                        ROCoutpost(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopost).fileNo;
                                        ROCgroupNopost(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopost).fileNo);
                                        ROCoutpost(no_ROCs).timeWindow=winNo;
                                        ROCbandwidthpost(no_ROCs)=bwii;
                                        auROCpost(no_ROCs)=ROCoutpost(no_ROCs).roc.AUC-0.5;
                                        p_valROCpost(no_ROCs)=ROCoutpost(no_ROCs).roc.p;
                                        
                                        p_vals_ROC=[p_vals_ROC ROCoutpost(no_ROCs).roc.p];
                                        
                                        
                                        %Are the delta dB LFP's different?
                                        
                                        %Hit
                                        p_val(no_dBs,bwii)=ranksum(this_delta_dB_powerpreHit,this_delta_dB_powerpostHit);
                                        p_vals=[p_vals p_val(no_dBs,bwii)];
                                        groupNopre(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                        groupNopost(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                        events(no_dBs)=1;
                                        
                                        
                                        %CR
                                        p_val(no_dBs+1,bwii)=ranksum(this_delta_dB_powerpreCR,this_delta_dB_powerpostCR);
                                        p_vals=[p_vals p_val(no_dBs+1,bwii)];
                                        groupNopre(no_dBs+1)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                        groupNopost(no_dBs+1)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                        events(no_dBs+1)=2;
                                        
                                        if p_val(no_dBs,bwii)<0.05
                                            dB_power_changeHit(no_ROCs)=1;
                                        else
                                            dB_power_changeHit(no_ROCs)=0;
                                        end
                                        
                                        if p_val(no_dBs+1,bwii)<0.05
                                            dB_power_changeCR(no_ROCs)=1;
                                        else
                                            dB_power_changeCR(no_ROCs)=0;
                                        end
                                        
                                        %Plot the points and save the data
                                        if groupNopre(no_dBs)==1
                                            
                                            
                                            %Hit, all points
                                            delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                            delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                            
                                            
                                            %CR, all points
                                            delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                            delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                            
                                        else
                                            if groupNopre(no_dBs)==3
                                                
                                                %Hit, all points
                                                delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                                delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                                
                                                
                                                %CR, all points
                                                delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                                delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                                %                                             figure(bwii+4+12)
                                                %                                             hold on
                                                %                                             plot([3 4],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            end
                                        end
                                    end
                                    
                                    no_dBs=no_dBs+2;
                                    
                                    
                                end
                                
                                
                                
                            end
                            
                            
                            
                        end
                        
>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855
                        
                        
                    end
                end
<<<<<<< HEAD
            end
        end
        fprintf(1, '\n\n')
        
        pFDRROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for significant difference of auROC p value from 0.5  = %d\n\n'],pFDRROC);
        
        
        
        fprintf(1, '\n\n')
        
        
        %Plot cumulative histos for auROCs
        dB_power_change=logical(dB_power_changeEv1+dB_power_changeEv2);
        figNo=0;
        p_val_ROC=[];
        
        try
            close(5)
        catch
        end
        figure(5)
        hold on
        x=0;
        
        for bwii=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            %Plot the histograms
            edges=[-0.5:0.05:0.5];
            pos2=[0.1 0.1 0.6 0.8];
            subplot('Position',pos2)
            hold on
            
            h2=histogram(auROCfp2(ROCbandwidthfp1==bwii),edges);
            h2.FaceColor='b';
            h1=histogram(auROCfp1(ROCbandwidthfp1==bwii),edges);
            h1.FaceColor='r';
            
            xlabel('auROC')
            ylabel('# of electrodes')
            legend(file_label{2},file_label{1})
            title(['auROC for ' freq_names{bwii}])
            xlim([-0.3 0.6])
            ylim([0 45])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Plot the single electrodes
            pos2=[0.8 0.1 0.1 0.8];
            subplot('Position',pos2)
            hold on
            for ii=1:length(auROCfp1)
                if ROCbandwidthfp1(ii)==bwii
                    plot([0 1],[auROCfp2(ii) auROCfp1(ii)],'-o', 'Color',[0.7 0.7 0.7])
                end
            end
            
            %PLot the mean and 95% CI
            plot([0 1],[mean(auROCfp2(ROCbandwidthfp1==bwii)) mean(auROCfp1(ROCbandwidthfp1==bwii))],'-k','LineWidth', 3)
            CI = bootci(1000, @mean, auROCfp2(ROCbandwidthfp1==bwii));
            plot([0 0],CI,'-b','LineWidth',3)
            plot(0,mean(auROCfp2(ROCbandwidthfp1==bwii)),'ob','MarkerSize', 10,'MarkerFace','b')
            CI = bootci(1000, @mean, auROCfp1(ROCbandwidthfp1==bwii));
            plot([1 1],CI,'-r','LineWidth',3)
            plot(1,mean(auROCfp1(ROCbandwidthfp1==bwii)),'or','MarkerSize', 10,'MarkerFace','r')
            ylabel('auROC')
            ylim([-0.2 0.5])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Do the statistics for auROC differences
            a={auROCfp2(ROCbandwidthfp1==bwii)' auROCfp1(ROCbandwidthfp1==bwii)'};
            mode_statcond='perm';
            [F df pval_auROCperm] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for permuted anovan for auROC S+ vs S- ' freq_names{bwii} '= %d\n\n'],  pval_auROCperm);
            pvals_auROCperm=[pvals_auROCperm pval_auROCperm];
            
            %Figure 5
            figure(5)
            
            percent_auROCfp2=100*sum(p_valROCfp2(ROCbandwidthfp1==bwii)<=pFDRROC)/sum(ROCbandwidthfp1==bwii);
            bar(x,percent_auROCfp2,'b')
            
            learn_sig(bwii)=sum(p_valROCfp2(ROCbandwidthfp1==bwii)<=pFDRROC);
            learn_not_sig(bwii)=sum(ROCbandwidthfp1==bwii)-sum(p_valROCfp2(ROCbandwidthfp1==bwii)<=pFDRROC);
            
            percent_auROCfp1=100*sum(p_valROCfp1(ROCbandwidthfp1==bwii)<=pFDRROC)/sum(ROCbandwidthfp1==bwii);
            bar(x+1,percent_auROCfp1,'r')
            
            prof_sig(bwii)=sum(p_valROCfp1(ROCbandwidthfp1==bwii)<=pFDRROC);
            prof_not_sig(bwii)=sum(ROCbandwidthfp1==bwii)-sum(p_valROCfp1(ROCbandwidthfp1==bwii)<=pFDRROC);
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            x=x+3;
            
        end
        
        figure(5)
        title('Percent singificant auROC')
        legend(file_label{2},file_label{1})
        ylim([0 100])
        
        pFDRanovanauROC=drsFDRpval(pval_auROCperm);
        fprintf(1, ['\npFDR for premuted anovan p value for difference between ' file_label{1} ' and ' file_label{2} ' for auROC = %d\n\n'],pFDRanovanauROC);
        
        
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix],'learn_sig','learn_not_sig','prof_sig','prof_not_sig');
        
        pffft=1;
        
=======
                
            end
            
            
            
            pFDRROC=drsFDRpval(p_vals_ROC);
            
            %Now plot the bar graphs and do anovan for LFP power
            p_vals_anovan=[];
            
            
            %Plot cumulative histos for auROCs
            dB_power_change=logical(dB_power_changeHit+dB_power_changeCR);
            figNo=0;
            p_val_ROC=[];
            pvals_auROCperm=[];
            
            
            for bwii=1:4
                n_cum=0;
                this_legend=[];
                data_auROC=[];
                pre_post_auROC=[];
                gr_auROC=[];
                for grs=1:2
                    
                    this_deltaauROC=[];
                    this_deltaauROC=auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))-auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii));
                    delta_meanauROC(ii_shift,grs,bwii)=mean(this_deltaauROC);
                    CI= bootci(1000, @mean, this_deltaauROC);
                    CIauROCdeltaLow(ii_shift,grs,bwii)=mean(this_deltaauROC)-CI(1);
                    CIauROCdeltaUpp(ii_shift,grs,bwii)=CI(2)-mean(this_deltaauROC);
                    
                    
                    meanauROCpre(ii_shift,grs,bwii)=mean(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                    CI= bootci(1000, @mean, auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                    CIauROCpreLowM(ii_shift,grs,bwii)=mean(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))-CI(1);
                    CIauROCpreUppM(ii_shift,grs,bwii)=CI(2)-mean(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                    
                    
                    meanauROCpost(ii_shift,grs,bwii)=mean(auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                    CI = bootci(1000, @mean, auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                    CIauROCpostLowM(ii_shift,grs,bwii)=mean(auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))-CI(1);
                    CIauROCpostUppM(ii_shift,grs,bwii)=CI(2)-mean(auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                    
                    delta_meanauROCM(ii_shift,grs,bwii)=meanauROCpost(ii_shift,grs,bwii)-meanauROCpre(ii_shift,grs,bwii);
                    
                    
                    %Save the data for anovan interaction
                    %Pre
                    data_auROC=[data_auROC auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))];
                    gr_auROC=[gr_auROC grs*ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
                    pre_post_auROC=[pre_post_auROC ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
                    
                    %Post
                    data_auROC=[data_auROC auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))];
                    gr_auROC=[gr_auROC grs*ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
                    pre_post_auROC=[pre_post_auROC 2*ones(1,sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)))];
                end
                %             figNo=figNo+2;
                %             x=x+3;
                
                %Calculate anovan for inteaction
                [p,tbl,stats]=anovan(data_auROC,{pre_post_auROC gr_auROC},'model','interaction','varnames',{'pre_vs_post','halo_vs_no_halo'},'display','off');
                fprintf(1, ['p value for anovan auROC interaction for ii_shift= %d ' freq_names{bwii} '= %d\n'],  ii_shift,p(3));
                p_aovan_int(ii_shift,bwii)=p(3);
                
            end
            
        end
        
        pFDRauROCint=drsFDRpval(p_aovan_int(:));
        fprintf(1, ['pFDR for auROC anovan interaction  = %d\n\n'],pFDRauROCint);
        
        %Plot the effect of silencing NA fibers with halorhodopsin
        figure(1)
        hold on
        plot([-0.6 0.6],[0 0],'-','LineWidth',3,'Color',[0.7 0.7 0.7])
        plot([0 0],[-0.25 0.25],'-','LineWidth',3,'Color',[0.7 0.7 0.7])
        delta_time=([0:10]*0.1-0.5);
        for bwii=1:4
            delta_auROC=zeros(1,length(shift_ii));
            delta_auROC=delta_meanauROC(:,1,bwii)-delta_meanauROC(:,2,bwii);
            delta_CIlow=sqrt(CIauROCdeltaLow(:,1,bwii).^2+CIauROCdeltaLow(:,2,bwii).^2);
            delta_CIupp=sqrt(CIauROCdeltaUpp(:,1,bwii).^2+CIauROCdeltaUpp(:,2,bwii).^2);
            
            for ii=1:length(delta_time)
                plot([delta_time(ii)+0.02*(bwii-2.5) delta_time(ii)+0.02*(bwii-2.5)],[delta_auROC(ii) delta_auROC(ii)+delta_CIupp(ii)],these_lines{bwii},'LineWidth',1)
                plot([delta_time(ii)+0.02*(bwii-2.5) delta_time(ii)+0.02*(bwii-2.5)],[delta_auROC(ii) delta_auROC(ii)-delta_CIlow(ii)],these_lines{bwii},'LineWidth',1)
            end
            plot(delta_time+0.02*(bwii-2.5),delta_auROC,'-','LineWidth',2,'Color',these_colors{bwii})
            plot(delta_time(p_aovan_int(:,bwii)<=pFDRauROCint)+0.02*(bwii-2.5),delta_auROC(p_aovan_int(:,bwii)<pFDRauROCint),'*','Color',these_colors{bwii},'MarkerSize',10)
            plot(delta_time(p_aovan_int(:,bwii)>pFDRauROCint)+0.02*(bwii-2.5),delta_auROC(p_aovan_int(:,bwii)>pFDRauROCint),'o','Color',these_colors{bwii},'MarkerFace',these_colors{bwii})
            
        end
        
        ylim([-0.25 0.25])
        xlim([-0.6 0.6])
        xlabel('dt to event (ms)')
        ylabel('delta auROC')
        title('Effect of silencing NA on auROC')
        legend('Theta','Beta','Low gamma','High gamma')
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2, 'Box', 'off')
>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855
        
    case 12
        %Justin
        
    case 13
        %Daniel
        %Generate Fig. 2  for Daniels' LFP power paper. For the proficient mice in the first and last sessions
        %plot the LFP spectrum for S+ vs S-, plot LFP power for S+ vs S- for each electrode and plot auROCs
        %NOTE: This does the analysis in all the files and DOES not distinguish between groups!!!
        no_dBs=1;
        delta_dB_power=[];
        no_ROCs=0;
        ROCout=[];
        p_vals_ROC=[];
        delta_dB_powerEv1=[];
        no_Ev1=0;
        noWB=0;
        delta_dB_powerEv1WB=[];
        delta_dB_powerEv2WB=[];
        
        fprintf(1, ['Pairwise auROC analysis for Fig 1 of Daniel''s paper\n\n'])
        p_vals=[];
        no_files=length(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        for fileNo=1:no_files
            for elec=1:16
                if sum(which_electrodes==elec)>0
                    lfpodNo_ref=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref)))
                        
                        
                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower))
                            
                            if (length(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(1,:))>=trials_to_process)
                                
                                trials=length(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(1,:));
                                mask=logical([zeros(1,trials-trials_to_process) ones(1,trials_to_process)]);
                                trials_in_eventEv1=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(event1,:)==1);
                                trials_in_eventEv2=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(event2,:)==1);
                                
                                if (sum(trials_in_eventEv1)>=min_trials_per_event) & (sum(trials_in_eventEv2)>=min_trials_per_event)
                                    
                                    lfpodNo=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                    
                                    % Ev1
                                    this_dB_powerrefEv1=zeros(sum(trials_in_eventEv1&mask),length(frequency));
                                    this_dB_powerrefEv1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower(trials_in_eventEv1&mask,:));
                                    
                                    
                                    this_dB_powerEv1=zeros(sum(trials_in_eventEv1&mask),length(frequency));
                                    this_dB_powerEv1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_eventEv1&mask,:));
                                    
                                    % Ev2
                                    this_dB_powerrefEv2=zeros(sum(trials_in_eventEv2&mask),length(frequency));
                                    this_dB_powerrefEv2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower(trials_in_eventEv2&mask,:));
                                    
                                    
                                    this_dB_powerEv2=zeros(sum(trials_in_eventEv2&mask),length(frequency));
                                    this_dB_powerEv2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_eventEv2&mask,:));
                                    
                                    
                                    %Wide band spectrum
                                    noWB=noWB+1;
                                    
                                    delta_dB_powerEv1WB(noWB,:)=mean(this_dB_powerEv1-this_dB_powerrefEv1,1);
                                    delta_dB_powerEv2WB(noWB,:)=mean(this_dB_powerEv2-this_dB_powerrefEv2,1);
                                    
                                    
                                    %Do per badwidth analysis
                                    for bwii=1:no_bandwidths
                                        
                                        no_ROCs=no_ROCs+1;
                                        this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                        
                                        %Enter the  Ev1
                                        this_delta_dB_powerEv1=zeros(sum(trials_in_eventEv1&mask),1);
                                        this_delta_dB_powerEv1=mean(this_dB_powerEv1(:,this_band)-this_dB_powerrefEv1(:,this_band),2);
                                        roc_data=[];
                                        roc_data(1:sum(trials_in_eventEv1&mask),1)=this_delta_dB_powerEv1;
                                        roc_data(1:sum(trials_in_eventEv1&mask),2)=zeros(sum(trials_in_eventEv1&mask),1);
                                        
                                        %Enter  Ev2
                                        total_trials=sum(trials_in_eventEv1&mask)+sum(trials_in_eventEv2&mask);
                                        this_delta_dB_powerEv2=zeros(sum(trials_in_eventEv2&mask),1);
                                        this_delta_dB_powerEv2=mean(this_dB_powerEv2(:,this_band)-this_dB_powerrefEv2(:,this_band),2);
                                        roc_data(sum(trials_in_eventEv1&mask)+1:total_trials,1)=this_delta_dB_powerEv2;
                                        roc_data(sum(trials_in_eventEv1&mask)+1:total_trials,2)=ones(sum(trials_in_eventEv2&mask),1);
                                        
                                        
                                        %Find  ROC
                                        ROCout(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                        ROCout(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo_ref).fileNo;
                                        ROCgroupNo(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNo_ref).fileNo);
                                        ROCout(no_ROCs).timeWindow=winNo;
                                        ROCbandwidth(no_ROCs)=bwii;
                                        auROC(no_ROCs)=ROCout(no_ROCs).roc.AUC-0.5;
                                        p_valROC(no_ROCs)=ROCout(no_ROCs).roc.p;
                                        
                                        p_vals_ROC=[p_vals_ROC ROCout(no_ROCs).roc.p];
                                        
                                        
                                        delta_dB_powerEv1(no_ROCs)=mean(this_delta_dB_powerEv1);
                                        delta_dB_powerEv2(no_ROCs)=mean(this_delta_dB_powerEv2);
                                        
                                        
                                        %Plot this point
                                        figure(bwii+1)
                                        pos2=[0.8 0.1 0.1 0.8];
                                        subplot('Position',pos2)
                                        hold on
                                        plot([1 0],[delta_dB_powerEv1(no_ROCs) delta_dB_powerEv2(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                                        
                                        
                                    end
                                    
                                    
                                    
                                else
                                    
                                    if (sum(trials_in_eventEv1)<min_trials_per_event)
                                        fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_eventEv1),event1, min_trials_per_event,files(fileNo),elec);
                                    end
                                    
                                    if (sum(trials_in_eventEv2)<min_trials_per_event)
                                        fprintf(1, ['%d trials in event No %d fewer than minimum trials per event %d for file No %d electrode %d\n'],sum(trials_in_eventEv2),event2, min_trials_per_event,files(fileNo),elec);
                                    end
                                    
                                end
                                
                            else
                                
                                fprintf(1, ['%d trials fewer than %d trials to process for file No %d electrode %d\n'],length(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(1,:)),trials_to_process,files(fileNo),elec);
                                
                            end
                        else
                            
                            fprintf(1, ['Empty allPower for file No %d electrode %d\n'],files(fileNo),elec);
                            
                        end
                        
                        
                    else
                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],files(fileNo),elec);
                        
                       
                    end
                end
            end
            
        end
        fprintf(1, '\n\n')
        
        
        %Now plot the bounded line for
        
        %Calculate the mean and 95% CI for Ev1
        dB_Ev1_ci=zeros(length(frequency),2);
        for ifreq=1:length(frequency)
            %             pd=fitdist(delta_dB_powerEv1WB(:,ifreq),'Normal');
            %             ci=paramci(pd);
            %             dB_Ev1_ci(ifreq)=pd.mu-ci(1,1);
            dB_Ev1_mean(ifreq)=mean(delta_dB_powerEv1WB(:,ifreq));
            CI = bootci(1000, @mean, delta_dB_powerEv1WB(:,ifreq));
            dB_Ev1_ci(ifreq,1)=CI(2)-dB_Ev1_mean(ifreq);
            dB_Ev1_ci(ifreq,2)=-(CI(1)-dB_Ev1_mean(ifreq));
        end
        
        figure(1)
        [hl1, hp1] = boundedline(frequency,dB_Ev1_mean, dB_Ev1_ci, 'r');
        
        %Calculate the mean and 95% CI for Ev2
        dB_Ev2_ci=zeros(length(frequency),2);
        for ifreq=1:length(frequency)
            dB_Ev2_mean(ifreq)=mean(delta_dB_powerEv2WB(:,ifreq));
            CI = bootci(1000, @mean, delta_dB_powerEv2WB(:,ifreq));
            dB_Ev2_ci(ifreq,1)=CI(2)-dB_Ev2_mean(ifreq);
            dB_Ev2_ci(ifreq,2)=-(CI(1)-dB_Ev2_mean(ifreq));
        end
        
        hold on
        [hl2, hp2] = boundedline(frequency,dB_Ev2_mean, dB_Ev2_ci, 'b');
        xlabel('Frequency (Hz)')
        ylabel('delta Power (dB)')
        legend([hl1 hl2],'S+','S-')
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        
        %Now plot the histograms and the average
        for bwii=1:4
            %Plot the average
            figure(bwii+1)
            pos2=[0.8 0.1 0.1 0.8];
            subplot('Position',pos2)
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            plot([1 0],[mean(delta_dB_powerEv1(ROCbandwidth==bwii)) mean(delta_dB_powerEv2(ROCbandwidth==bwii))],'-k','LineWidth', 3)
            CI = bootci(1000, @mean, delta_dB_powerEv1(ROCbandwidth==bwii));
            plot([1 1],CI,'-r','LineWidth',3)
            plot(1,mean(delta_dB_powerEv1(ROCbandwidth==bwii)),'or','MarkerSize', 10,'MarkerFace','r')
            CI = bootci(1000, @mean, delta_dB_powerEv2(ROCbandwidth==bwii));
            plot([0 0],CI,'-b','LineWidth',3)
            plot(0,mean(delta_dB_powerEv2(ROCbandwidth==bwii)),'ob','MarkerSize', 10,'MarkerFace','b')
            ylabel('delta Power (dB)')
            ylim([-10 15])
            
            %Plot the histograms
            
            maxdB=max([max(delta_dB_powerEv1(ROCbandwidth==bwii)) max(delta_dB_powerEv2(ROCbandwidth==bwii))]);
            mindB=min([min(delta_dB_powerEv1(ROCbandwidth==bwii)) min(delta_dB_powerEv2(ROCbandwidth==bwii))]);
            edges=[-15:1:15];
            pos2=[0.1 0.1 0.6 0.8];
            subplot('Position',pos2)
            hold on
            
            h1=histogram(delta_dB_powerEv2(ROCbandwidth==bwii),edges);
            h1.FaceColor='b';
            h2=histogram(delta_dB_powerEv1(ROCbandwidth==bwii),edges);
            h2.FaceColor='r';
            xlabel('delta Power (dB)')
            ylabel('# of electrodes')
            legend('S-','S+')
            xlim([-12 12])
            ylim([0 70])
            title(freq_names{bwii})
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            
            
            a={ delta_dB_powerEv1(ROCbandwidth==bwii)' delta_dB_powerEv2(ROCbandwidth==bwii)'};
            mode_statcond='perm';
            [F df pvals_perm(bwii)] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for premuted anovan dB delta power S+ vs S- ' freq_names{bwii} '= %d\n'],  pvals_perm(bwii));
            
        end
        
        pFDRanovan=drsFDRpval(pvals_perm);
        fprintf(1, ['pFDR for premuted anovan p value  = %d\n\n'],pFDRanovan);
        
        
        
        fprintf(1, '\n\n')
        
        
        pFDRauROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for auROC  = %d\n\n'],pFDRauROC);
        %Plot cumulative histos for auROCs
        
        figNo=5;
        p_val_ROC=[];
        edges=-0.5:0.05:0.5;
        
        for bwii=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            n_cum=0;
            this_legend=[];
            
            histogram(auROC(( p_valROC>pFDRauROC)&(ROCbandwidth==bwii)),edges)
            histogram(auROC(( p_valROC<=pFDRauROC)&(ROCbandwidth==bwii)),edges)
            legend('auROC not singificant','auROC significant')
            title(['Histogram for ' freq_names{bwii} ' auROC for LFPs'])
            xlim([-0.2 0.6])
            ylim([0 30])
        end
        
         
        
        %Plot percent significant ROC
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        figure(figNo)
        
        hold on
        
        for bwii=1:4
            bar(bwii,100*sum(( p_valROC<=pFDRauROC)&(ROCbandwidth==bwii))/sum((ROCbandwidth==bwii)))
            auROC_sig.sig(bwii)=sum(( p_valROC<=pFDRauROC)&(ROCbandwidth==bwii));
            auROC_sig.not_sig(bwii)=sum((ROCbandwidth==bwii))-sum(( p_valROC<=pFDRauROC)&(ROCbandwidth==bwii));
        end
        title('Percent auROC significantly different from zero')
        ylim([0 100])
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        pffft=1;
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix],'auROC_sig');
end


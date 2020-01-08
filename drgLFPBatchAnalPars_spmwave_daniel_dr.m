function handles_pars=drgLFPBatchAnalPars_spmwave_daniel_dr

%This file has the parameters specified for data analysis in
%drgAnalysisBatchLFP

%Initialize all variables
handles_pars.winNo=[];
handles_pars.trials_to_process=[];
handles_pars.which_display=[];
handles_pars.eventType=[];
handles_pars.evTypeLabels=[];
handles_pars.file_pairs=[];
handles_pars.trials_to_process=[];
handles_pars.min_trials_per_event=[];
handles_pars.shift_time=[];
handles_pars.shift_from_event=[];
handles_pars.grpre=[];
handles_pars.grpost=[];
handles_pars.output_suffix=[];
handles_pars.front_mask=[];
handles_pars.file_label=[];
handles_pars.refWin=[];
handles_pars.concs2=[];



handles_pars.winNo=2;
handles_pars.refWin=1;
% which_display=3;
handles_pars.which_display=24;

%
% eventType=[2 5];
% evTypeLabels={'Hit','CR'};
% comp_window=10;

% evTypeNos=[2 4 5 7];
% evTypeLabels={'Hit','Miss','CR','FA'};

% handles_pars.eventType=[8:13];
% handles_pars.evTypeLabels={'Hi1','Hi2','Hi3','Low4','Low5','Low6'};
% handles_pars.concs=[log10(10) log10(3.2) log10(1) log10(0.32) log10(0.1) log10(0.032)];
% handles_pars.concs2 = [0.032 0.1 0.32 1 3.2 10];

% handles_pars.eventType=[2 4 5 7];
% handles_pars.evTypeLabels={'Hit','Miss','CR','FA'};

handles_pars.eventType=[3 6];
handles_pars.evTypeLabels={'S+','S-',};
handles_pars.concs2 = {'S-', 'S+'};

% Enter the files to be processed
handles_pars.files=[1:22];     %EAPA Forward and reversal ethyl acetate 
handles_pars.output_suffix='_EAPA.mat';


% handles_pars.files=[53:80];     %APEB Forward aceto 
% handles_pars.output_suffix='_APEB.mat';

% 
% handles_pars.files=[95:112];     %Forward iso both control and exp
% handles_pars.output_suffix='_IAMO.mat';


% handles_pars.files=[105:112 118:121];     %IAMO iso experimental vs iso experimental laser

% handles_pars.files=[95:104 113:117];     %IAMO iso control vs iso control laser

% handles_pars.files=[65:80 87:94];     %aceto experimental vs aceto experimental laser

% handles_pars.files=[53:64 81:86];     %aceto control vs aceto control laser

% handles_pars.files=[15:22 49:52];     %ethyl experimental vs ethyl experimental laser

% handles_pars.files=[1:14 42:48];     %ethyl control vs ethyl control laser

handles_pars.trials_to_process=30;
handles_pars.min_trials_per_event=1;


handles_pars.percent_windows=[80 100;
    45 65];
handles_pars.per_lab = {'Proficient','Naive'};

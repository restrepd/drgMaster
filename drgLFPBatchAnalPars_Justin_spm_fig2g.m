function handles_pars=drgLFPBatchAnalPars_Justin_spm_dr

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
handles_pars.which_display=19;

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

% handles_pars.files=[1:5 6:7 8:10 12:15 17:18 42:44 46:49];     %IAAP 

% handles_pars.files=[1:8 12 42:43 46];     %IAAP fwd only
% handles_pars.output_suffix='_IAAP.mat';

handles_pars.files=[9:10 13:15 44];     %IAAP reverse only
handles_pars.output_suffix='_IAAPr.mat';
% 
% handles_pars.files=[8 12 42:43];     %IAAP fwd matched to reverse only
% handles_pars.output_suffix='_IAAPfr.mat';

% handles_pars.files=[11 16 19:23 27:29 33:36 40:41 45 50 51];     %EAPA2


% handles_pars.files=[11 16 19 21:23 27 33 40 45 50];     %EAPA fwd only
% handles_pars.output_suffix='_EAPA.mat';

% handles_pars.files=[17 18 24:26 30:32 37:39];     %IAMO fwd only
% handles_pars.output_suffix='_IAMO.mat';

handles_pars.trials_to_process=30;
handles_pars.min_trials_per_event=1;


handles_pars.percent_windows=[80 100;
    45 65];
handles_pars.per_lab = {'Proficient','Naive'};

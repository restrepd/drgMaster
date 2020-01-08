function handles=drgbChoicesDiscriminantJustin_spm_wavephase_05172019_IAAP
%Output file name 
handles.drgb.outFileName='spm_discriminant_LFP_wavephase_05172019_IAAP.mat';
handles.drgb.outPathName='C:\DATA\Justin\spmout\';

%Which event do you want as time start reference? This should the event that
%encompasses all events you want to obtain results on. All analysis will be
%timed with respect to this event. Thus, if you are interested in S+, S-,
%Hit, CR, etc for spm you should make this event OdorOn=2
%If you are interested in tstart, this should be tstart=1
handles.drgbchoices.referenceEvent=2; 

%Which events do you want evaluated? 
%Remember these choices are different depending on how you saved the
%jt_times file in drta. e.g. 5 and 11 are S+ and S- if you saved the file
%under dropcspm
% For spm these are the events you can choose:
%Enter the event type
%   Events 1 through 6
%     'TStart'    'OdorOn'    'Hit'    'HitE'    'S+'    'S+E'
%   Events 7 through 13
%     'Miss'    'MissE'    'CR'    'CRE'    'S-'    'S-E'    'FA'
%   Events 14 through 19
%     'FAE'    'Reinf'    'L+'    'L-' 'S+TStart' 'S-TStart'
%   'S+TStart' = 18

%The program will keep track of all of these events
handles.drgbchoices.evTypeNos=[3 5 7 9 11 13]; 
%OdorOn, Hit, S+, Miss, CR, S-, FA

%The program will use these events to do the discrimination analysis
handles.drgbchoices.events_to_discriminate=[5 11]; %These are S+ and S-

%Which electrodes should be used
handles.drgbchoices.which_electrodes=[1:16];

%Which time window do you want to be evaluated?
handles.drgbchoices.timeStart=-2;
handles.drgbchoices.timeEnd=5;

%Reference LFP
handles.drgbchoices.subtractRef=1;


%which discriminant analysis
%1 Perceptron for power LFP
%2 Linear discriminant analysis for power LFP
handles.drgbchoices.which_discriminant=[10 11];

%Percent windows
handles.drgbchoices.percent_windows=[80 100;
    45 65];
handles.drgbchoices.per_lab = {'proficient','naive'};


%First file to process
handles.drgb.first_file = 1;

%Which jt_times_ files do you want to process?
% handles.drgbchoices.no_files=8;

handles.drgbchoices.PathName='C:\DATA\Justin\SPM\';

%spm
handles.drgbchoices.MouseName{1}='M4';
handles.drgbchoices.mouse_no(1:5)=1;
handles.drgbchoices.group_no(1:5)=1; % fwd
handles.drgbchoices.session_no(1:5) = (1:5);
handles.drgbchoices.FileName{1}='jt_times_M4_spm_iso_ace_170809_080253.mat';
handles.drgbchoices.FileName{2}='jt_times_M4_spm_iso_ace_170810_073545.mat';
handles.drgbchoices.FileName{3}='jt_times_M4_spm_ace_iso_170830_085055.mat';
handles.drgbchoices.FileName{4}='jt_times_M4_spm_ace_iso_170831_090444.mat';
handles.drgbchoices.FileName{5}='jt_times_M4_spm_iso_ace_170907_080617.mat';

handles.drgbchoices.MouseName{2}='M5';
handles.drgbchoices.mouse_no(6:7)=2;
handles.drgbchoices.group_no(6:7)=1; % fwd
handles.drgbchoices.session_no(6:7) = (1:2);
handles.drgbchoices.FileName{6}='jt_times_M5_spm_iso_ace_170809_091037.mat';
handles.drgbchoices.FileName{7}='jt_times_M5_spm_iso_ace_170810_085523.mat';

handles.drgbchoices.MouseName{3}='T0';
handles.drgbchoices.mouse_no(8:10)=3;
handles.drgbchoices.group_no(8)=1; % fwd
handles.drgbchoices.group_no([9 10])=2; % rev
handles.drgbchoices.session_no(8:10) = (1:3);
handles.drgbchoices.FileName{8}='jt_times_T0_spm_iso_ace_180130_111647.mat';
handles.drgbchoices.FileName{9}='jt_times_T0_spmr_ace_iso_180131_103513.mat';
handles.drgbchoices.FileName{10}='jt_times_T0_spmr_ace_iso_180201_104706.mat';

handles.drgbchoices.MouseName{4}='T2';
handles.drgbchoices.mouse_no(11:14)=4;
handles.drgbchoices.group_no(11)=1; % fwd
handles.drgbchoices.group_no(12:14)=2; % rev
handles.drgbchoices.session_no(11:14) = (1:4);
handles.drgbchoices.FileName{11}='jt_times_T2_spm_iso_ace_180130_122900.mat';
handles.drgbchoices.FileName{12}='jt_times_T2_spmr_ace_iso_180131_115342.mat';
handles.drgbchoices.FileName{13}='jt_times_T2_spmr_ace_iso_180201_113431.mat';
handles.drgbchoices.FileName{14}='jt_times_T2_spmr_ace_iso_180202_114037.mat';

handles.drgbchoices.MouseName{5}='R8';
handles.drgbchoices.mouse_no(15:17)=5;
handles.drgbchoices.group_no(15:16)=1; % fwd
handles.drgbchoices.group_no(17)=2; % rev
handles.drgbchoices.session_no(15:17) = (1:3);
handles.drgbchoices.FileName{15}='jt_times_R8_spm_iso_ace_190222_162007.mat';
handles.drgbchoices.FileName{16}='jt_times_R8_spm_iso_ace_190224_104344.mat';
handles.drgbchoices.FileName{17}='jt_times_R8_spm_ace_iso_190225_083746.mat';

handles.drgbchoices.MouseName{6}='R9';
handles.drgbchoices.mouse_no(18:21)=6;
handles.drgbchoices.group_no(18)=1; % fwd
handles.drgbchoices.group_no(19:21)=2; % rev
handles.drgbchoices.session_no(18:21) = (1:4);
handles.drgbchoices.FileName{18}='jt_times_R9_spm_iso_ace_190222_165359.mat';
handles.drgbchoices.FileName{19}='jt_times_R9_spm_ace_iso_190224_112151.mat';
handles.drgbchoices.FileName{20}='jt_times_R9_spm_ace_iso_190225_091732.mat';
handles.drgbchoices.FileName{21}='jt_times_R9_spm_ace_iso_190226_083131.mat';

handles.drgbchoices.no_files=max(size(handles.drgbchoices.FileName(1,:)));

% %Enter the mouse number for each file
% handles.drgbchoices.mouse_no=[1 1 1 1 1 2 2 3 3 3 3 4 4 4 4 4 5 5 5 5 3 3 3 7 7 7 7 7 7 8 8 8 8 8 8 8];
% 
% %Enter the group number for each file
% handles.drgbchoices.group_no=[1 1 1 1 1 1 1 1 2 2 1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 2 2 2];
% 
% %Enter the session number for each file
% handles.drgbchoices.session_no=[1 2 3 4 5 1 2 1 2 3 4 1 2 3 4 5 1 2 3 4 1 2 3 1 2 3 4 5 6 1 2 3 4 5 6 7];

%Enter the name of each group No
handles.drgbchoices.group_no_names{1}='fwd';
handles.drgbchoices.group_no_names{2}='rev';

%Theta
handles.drgbchoices.lowF(1)=6;
handles.drgbchoices.highF(1)=14;

%Beta
handles.drgbchoices.lowF(2)=15;
handles.drgbchoices.highF(2)=30;

%Low gamma
handles.drgbchoices.lowF(3)=35;
handles.drgbchoices.highF(3)=55;

%High gamma
handles.drgbchoices.lowF(4)=65;
handles.drgbchoices.highF(4)=95;

% %Wideband
% handles.drgbchoices.lowF(5)=4;
% handles.drgbchoices.highF(5)=95;


handles.drgbchoices.bwlabels{1}='Theta';
handles.drgbchoices.bwlabels{2}='Beta';
handles.drgbchoices.bwlabels{3}='Low gamma';
handles.drgbchoices.bwlabels{4}='High gamma';
% handles.drgbchoices.bwlabels{5}='Wideband';


%PAC parameters
handles.drgbchoices.no_PACpeaks=3;
handles.drgbchoices.PACpeakLowF=6;
handles.drgbchoices.PACpeakHighF=14;
handles.drgbchoices.PACburstLowF=[15 35 65];
handles.drgbchoices.PACburstHighF=[30 55 95];
handles.drgbchoices.PACnames{1}='Beta';
handles.drgbchoices.PACnames{2}='Low Gamma';
handles.drgbchoices.PACnames{3}='High Gamma';

handles.drgbchoices.n_phase_bins=50;



%LFP power parameters
handles.LFPPowerSpectrumLowF=4;
handles.LFPPowerSpectrumHighF=95;



%All other choices
handles.analysisNoOsc=1;
handles.analysisNoSpikes=1;
handles.peakLFPNo=1; %3 for Anan, 2 for Justin
handles.burstLFPNo=1;
handles.evTypeNo=5;
handles.unitNo=1;
handles.time_start=-0.2;
handles.time_end=2.2; %2.2 for Shane
handles.sessionNo=1;
handles.trialNo=1;
handles.lastTrialNo=60;
handles.data_vs_simulate=0;
handles.window=1; %This is the FFT window in sec
handles.noverlap=handles.window*0.9;
handles.deltaLowF_PAC=1;  %2
handles.deltaHighF_PAC=5; %30 
handles.bandwidth_lowF=3;
handles.bandwidth_highF=25;
handles.which_method=1;
handles.amplitudeHz=80;
handles.phaseHz=10;
handles.unitNo=1;
handles.analysisNoBeh=1;
handles.notch60=1;
handles.displayData=0;
handles.save_drgb=1;
handles.save_events=1;
handles.perTetrode=0;
handles.max_dt_between_events=3.5;
handles.read_entire_file=0;
handles.time_pad=0.2;
handles.startRef=-2-handles.time_pad; %Remember these have the pad
handles.endRef=0+handles.time_pad;
handles.dt_lick=0.1;
handles.smallest_inter_lick_interval=0.02;  %Note: this is used to reject "lick" events due to noise





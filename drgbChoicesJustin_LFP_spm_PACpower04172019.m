function handles=drgbChoicesJustin_LFP_spm_PACpower04172019
%Output file name 
% handles.drgb.outFileName='spm_LFP_20180809.mat';
% handles.drgb.outPathName='D:\data\spm\';

handles.drgb.outFileName='spm_LFP_PACpower04172019.mat';
handles.drgb.outPathName='C:\Data\Justin\spmout\';


%test the batch? 1=yes, 0=no
handles.drgbchoices.test_batch=0;


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

handles.drgbchoices.evTypeNos=[2 3 5 7 9 11 13]; 
%OdorOn, Hit, S+, Miss, CR, S-, FA

%Number of LFP electrodes
%Anan 1, others 16
handles.drgbchoices.no_LFP_elect=16;

%Which time windows do you want to be evaluated?
%For example if you want -2 to 0 and 0 to 2 sec:
% handles.drgbchoices.timeStart=[-2 0];
% handles.drgbchoices.timeEnd=[0 2];

handles.drgbchoices.noWindows=2;
handles.drgbchoices.timeStart=[-1.5 0.5]; %3.15 to 3.6 is reinforcement, 3.8 to 4.8 something interesing happens
handles.drgbchoices.timeEnd=[0 2.5];

%Reference LFP
handles.drgbchoices.subtractRef = 1;
handles.subtractRef=1;
handles.time_pad=0.2;
handles.startRef=-2-handles.time_pad; %Remember these have the pad
handles.endRef=0+handles.time_pad;



%What analysis do you want to run?
%
% 1=drgThetaAmpPhase
% 2=drgGetLFPPowerForThisEvTypeNo

handles.drgbchoices.analyses=[1 2];



%First file to process
handles.drgb.first_file = 1;

%Which jt_times_ files do you want to process?
% handles.drgbchoices.no_files=8;

handles.drgbchoices.PathName='C:\DATA\Justin\SPM\';

%spm
handles.drgbchoices.MouseName{1}='M4';
handles.drgbchoices.mouse_no(1:10)=1;
handles.drgbchoices.group_no(1:5)=1; % fwd
% handles.drgbchoices.group_no(6:10)=2; % rev
handles.drgbchoices.session_no(1:10) = (1:10);
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
handles.drgbchoices.mouse_no(8:11)=3;
handles.drgbchoices.group_no([8 11])=1; % fwd
handles.drgbchoices.group_no([9 10])=2; % rev
handles.drgbchoices.session_no(8:11) = (1:4);
handles.drgbchoices.FileName{8}='jt_times_T0_spm_iso_ace_180130_111647.mat';
handles.drgbchoices.FileName{9}='jt_times_T0_spmr_ace_iso_180131_103513.mat';
handles.drgbchoices.FileName{10}='jt_times_T0_spmr_ace_iso_180201_104706.mat';
handles.drgbchoices.FileName{11}='jt_times_T0_spm_EAPA_180208_115640.mat';

handles.drgbchoices.MouseName{4}='T2';
handles.drgbchoices.mouse_no(12:16)=4;
handles.drgbchoices.group_no([12 16])=1; % fwd
handles.drgbchoices.group_no(13:15)=2; % rev
handles.drgbchoices.session_no(12:16) = (1:5);
handles.drgbchoices.FileName{12}='jt_times_T2_spm_iso_ace_180130_122900.mat';
handles.drgbchoices.FileName{13}='jt_times_T2_spmr_ace_iso_180131_115342.mat';
handles.drgbchoices.FileName{14}='jt_times_T2_spmr_ace_iso_180201_113431.mat';
handles.drgbchoices.FileName{15}='jt_times_T2_spmr_ace_iso_180202_114037.mat';
handles.drgbchoices.FileName{16}='jt_times_T2_spm_EAPA_180208_122330.mat';

handles.drgbchoices.MouseName{5}='R2';
handles.drgbchoices.mouse_no(17:20)=5;
handles.drgbchoices.group_no(17:20)=1; % fwd
handles.drgbchoices.session_no(17:20) = (1:4);
handles.drgbchoices.FileName{17}='jt_times_R2_spm_iso_mo_180507_100258.mat';
handles.drgbchoices.FileName{18}='jt_times_R2_spm_iso_mo_180716_102103.mat';
handles.drgbchoices.FileName{19}='jt_times_R2_spm_EAPA_180717_101444.mat';
handles.drgbchoices.FileName{20}='jt_times_R2_spm_PAEA_180718_101418.mat';

handles.drgbchoices.MouseName{6}='R3';
handles.drgbchoices.mouse_no(21:23)=6;
handles.drgbchoices.group_no(21:23)=1; % fwd
handles.drgbchoices.session_no(21:23) = (1:3);
handles.drgbchoices.FileName{21}='jt_times_R3_spm_EAPA_180726_104350.mat';
handles.drgbchoices.FileName{22}='jt_times_R3_spm_EAPA_180727_113406.mat';
handles.drgbchoices.FileName{23}='jt_times_R3_spm_EAPA_180814_131722.mat';

handles.drgbchoices.MouseName{7}='R4';
handles.drgbchoices.mouse_no(24:29)=7;
handles.drgbchoices.group_no(24:27)=1; % fwd
handles.drgbchoices.group_no(28:29)=2; % rev
handles.drgbchoices.session_no(24:29) = (1:6);
handles.drgbchoices.FileName{24}='jt_times_R4_spm_isomo_180801_103554.mat';
handles.drgbchoices.FileName{25}='jt_times_R4_spm_isomo_180802_100655.mat';
handles.drgbchoices.FileName{26}='jt_times_R4_spm_isomo_180806_100637.mat';
handles.drgbchoices.FileName{27}='jt_times_R4_spm_PAEA_180807_101303.mat';
handles.drgbchoices.FileName{28}='jt_times_R4_spmr_EAPA_180808_100929.mat';
handles.drgbchoices.FileName{29}='jt_times_R4_spmr_EAPA_180809_101536.mat';

handles.drgbchoices.MouseName{8}='R5';
handles.drgbchoices.mouse_no(30:36)=8;
handles.drgbchoices.group_no(30:33)=1; % fwd
handles.drgbchoices.group_no(34:36)=2; % rev
handles.drgbchoices.session_no(30:36) = (1:7);
handles.drgbchoices.FileName{30}='jt_times_R5_spm_isomo_180801_112402.mat';
handles.drgbchoices.FileName{31}='jt_times_R5_spm_isomo_180802_111225.mat';
handles.drgbchoices.FileName{32}='jt_times_R5_spm_isomo_180806_104729.mat';
%jt_times_R5_spm_PAEA_180807_111935.mat will be used for simulations
handles.drgbchoices.FileName{33}='jt_times_R5_spm_PAEA_180807_111935.mat';
handles.drgbchoices.FileName{34}='jt_times_R5_spmr_EAPA_180808_111500.mat';
handles.drgbchoices.FileName{35}='jt_times_R5_spmr_EAPA_180809_105230.mat';
handles.drgbchoices.FileName{36}='jt_times_R5_spmr_EAPA_180810_111557.mat';
handles.drgbchoices.FileName{37}='jt_times_R5_spm_iso_mo_181119_101055.mat';

handles.drgbchoices.MouseName{9}='R6';
handles.drgbchoices.mouse_no(38:41)=9;
handles.drgbchoices.group_no(38:40)=1; % fwd
handles.drgbchoices.group_no(41)=2; % rev
handles.drgbchoices.session_no(38:41) = (1:4);
handles.drgbchoices.FileName{38}='jt_times_R6_spm_1iso_MO_181101_105120.mat';
handles.drgbchoices.FileName{39}='jt_times_R6_spm_iso_mo_181119_094140.mat';
handles.drgbchoices.FileName{40}='jt_times_R6_spm_PAEA_190210_143426.mat';
handles.drgbchoices.FileName{41}='jt_times_R6_spm_EAPA_190211_080729.mat';

handles.drgbchoices.MouseName{10}='R8';
handles.drgbchoices.mouse_no(42:45)=10;
handles.drgbchoices.group_no([42 43 45])=1; % fwd
handles.drgbchoices.group_no(44)=2; % rev
handles.drgbchoices.session_no(42:45) = (1:4);
handles.drgbchoices.FileName{42}='jt_times_R8_spm_iso_ace_190222_162007.mat';
handles.drgbchoices.FileName{43}='jt_times_R8_spm_iso_ace_190224_104344.mat';
handles.drgbchoices.FileName{44}='jt_times_R8_spm_ace_iso_190225_083746.mat';
handles.drgbchoices.FileName{45}='jt_times_R8_spm_PAEA_190226_085848.mat';

handles.drgbchoices.MouseName{11}='R9';
handles.drgbchoices.mouse_no(46:51)=11;
handles.drgbchoices.group_no([46 50])=1; % fwd
handles.drgbchoices.group_no([47:49 51])=2; % rev
handles.drgbchoices.session_no(46:51) = (1:6);
handles.drgbchoices.FileName{46}='jt_times_R9_spm_iso_ace_190222_165359.mat';
handles.drgbchoices.FileName{47}='jt_times_R9_spm_ace_iso_190224_112151.mat';
handles.drgbchoices.FileName{48}='jt_times_R9_spm_ace_iso_190225_091732.mat';
handles.drgbchoices.FileName{49}='jt_times_R9_spm_ace_iso_190226_083131.mat';
handles.drgbchoices.FileName{50}='jt_times_R9_spm_PAEA_190227_090029.mat';
handles.drgbchoices.FileName{51}='jt_times_R9_spm_EAPA_190227_083042.mat';

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

%PAC parameters
handles.no_PACpeaks=3;
handles.PACpeakLowF=6;
handles.PACpeakHighF=14;
handles.PACburstLowF=[15 35 65];
handles.PACburstHighF=[30 55 95];
handles.PACnames{1}='Beta PAC';
handles.PACnames{2}='Low Gamma PAC';
handles.PACnames{3}='High Gamma PAC';

handles.n_phase_bins=50;

%LFP power parameters
handles.LFPPowerSpectrumLowF=4;
handles.LFPPowerSpectrumHighF=95;

% handles.no_lfp_bands=4;
% 
% %Theta
% handles.lowF(1)=2;
% handles.highF(1)=12;
% 
% %Beta
% handles.lowF(2)=15;
% handles.highF(2)=30;
% 
% %Low gamma
% handles.lowF(3)=35;
% handles.highF(3)=55;
% 
% %High gamma
% handles.lowF(4)=65;
% handles.highF(4)=95;


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
handles.max_events_per_sec=20;
handles.dt_lick=0.1;
handles.smallest_inter_lick_interval=0.02;  %Note: this is used to reject "lick" events due to noise



function handles=drgbChoicesLFPexample
%Output file name 
handles.drgb.outFileName='acetoethylben_electrodetest.mat';
handles.drgb.outPathName='C:\Users\Daniel\Desktop\acetotheylben\';

%Number of LFP electrodes
% Anan: please enter handles.drgbchoices.no_LFP_elect=1;
% Restrepo lab enter handles.drgbchoices.no_LFP_elect=16;
handles.drgbchoices.no_LFP_elect=16;

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

%Which time windows do you want to be evaluated?
%For example if you want -2 to 0 and 0 to 2 sec:
% handles.drgbchoices.timeStart=[-2 0];
% handles.drgbchoices.timeEnd=[0 2];

handles.drgbchoices.noWindows=3;
handles.drgbchoices.timeStart=[-1.5 0.5 3.15]; %3.15 to 3.6 is reinforcement, 3.8 to 4.8 something interesing happens
handles.drgbchoices.timeEnd=[0 2.5 3.6];

%Reference LFP
handles.subtractRef=1;
handles.time_pad=0.2;
handles.startRef=-2-handles.time_pad; %Remember these have the pad
handles.endRef=0+handles.time_pad;



%What analysis do you want to run?
%
% 1=drgThetaAmpPhase
% 2=drgGetLFPPowerForThisEvTypeNo

handles.drgbchoices.analyses=[1 2];

%Which is the first file to process?
handles.drgb.first_file=1;

%Which jt_times_ files do you want to process?
handles.drgbchoices.no_files=16;

%no laser
% handles.drgbchoices.PathName{1}='C:\Users\Daniel\Desktop\acetotheylben\nolaser\159863\';
% handles.drgbchoices.FileName{1}='jt_times_103016female159863acetoethylben.mat';

%For troubleshooting in Diego's computer
handles.drgbchoices.PathName{1}='/Users/restrepd/Documents/Projects/Daniel/matlabprocessdatabatchprogram/acetotheylben/nolaser/159863/';
handles.drgbchoices.FileName{1}='jt_times_103016female159863acetoethylben_nl.mat';


handles.drgbchoices.PathName{2}='C:\Users\Daniel\Desktop\acetotheylben\nolaser\159865\';
handles.drgbchoices.FileName{2}='jt_times_103116male159865acetoethylben.mat';
handles.drgbchoices.PathName{3}='C:\Users\Daniel\Desktop\acetotheylben\nolaser\159866\';
handles.drgbchoices.FileName{3}='jt_times_102816male159866acetoethylben.mat';
handles.drgbchoices.PathName{4}='C:\Users\Daniel\Desktop\acetotheylben\nolaser\159867\';
handles.drgbchoices.FileName{4}='jt_times_103016male159867acetoethylben.mat';
handles.drgbchoices.PathName{5}='C:\Users\Daniel\Desktop\acetotheylben\nolaser\159869\';
handles.drgbchoices.FileName{5}='jt_times_103116male159869acetoethylben.mat';

%no laser control
handles.drgbchoices.PathName{6}='C:\Users\Daniel\Desktop\acetotheylben\nolaser\control\';
handles.drgbchoices.FileName{6}='jt_times_103116male159868acetoethylben.mat';

%laser
handles.drgbchoices.PathName{7}='C:\Users\Daniel\Desktop\acetotheylben\laser\159863\';
handles.drgbchoices.FileName{7}='jt_times_103016female159863acetoethlben.mat';
handles.drgbchoices.PathName{8}='C:\Users\Daniel\Desktop\acetotheylben\laser\159865\';
handles.drgbchoices.FileName{8}='jt_times_103116male159865acetoethlben.mat';
handles.drgbchoices.PathName{9}='C:\Users\Daniel\Desktop\acetotheylben\laser\159866\';
handles.drgbchoices.FileName{9}='jt_times_10286male1598661acetoethylben.mat';
handles.drgbchoices.PathName{10}='C:\Users\Daniel\Desktop\acetotheylben\laser\159867\';
handles.drgbchoices.FileName{10}='jt_times_11216male159867acetoethlben.mat';
handles.drgbchoices.PathName{11}='C:\Users\Daniel\Desktop\acetotheylben\laser\159869\';
handles.drgbchoices.FileName{11}='jt_times_11116female159869acetoethylben.mat';

%laser control
handles.drgbchoices.PathName{12}='C:\Users\Daniel\Desktop\acetotheylben\laser\control\';
handles.drgbchoices.FileName{12}='jt_times_103116male159868acetoethlben.mat';

%Let's add

%159863
handles.drgbchoices.PathName{13}='C:\Users\Daniel\Desktop\acetotheylben\nolaser\159863\';
handles.drgbchoices.FileName{13}='jt_times_102916female159863acetoethylben_nl.mat';
handles.drgbchoices.PathName{14}='C:\Users\Daniel\Desktop\acetotheylben\laser\159863\';
handles.drgbchoices.FileName{14}='jt_times_102916female159863acetoethlben_l.mat';

%146679
handles.drgbchoices.PathName{15}='C:\Users\Daniel\Desktop\acetotheylben\nolaser\146679\';
handles.drgbchoices.FileName{15}='jt_times_102916female146679acetoethylben_nl.mat';
handles.drgbchoices.PathName{16}='C:\Users\Daniel\Desktop\acetotheylben\laser\146679\';
handles.drgbchoices.FileName{16}='jt_times_102916female1466791peracetoethlben_l.mat';


%Enter the group number for each file
handles.drgbchoices.group_no=[1 1 1 1 1 3 2 2 2 2 2 4 1 2 1 2];

%Enter the name of each group No
handles.drgbchoices.group_no_names{1}='NL';
handles.drgbchoices.group_no_names{2}='L';
handles.drgbchoices.group_no_names{3}='NLc';
handles.drgbchoices.group_no_names{4}='Lc';

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

handles.no_lfp_bands=4;

%Theta
handles.lowF(1)=6;
handles.highF(1)=14;

%Beta
handles.lowF(2)=15;
handles.highF(2)=30;

%Low gamma
handles.lowF(3)=35;
handles.highF(3)=55;

%High gamma
handles.lowF(4)=65;
handles.highF(4)=95;




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




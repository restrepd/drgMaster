function handles=drgbChoicesDiscriminantDanielolfa_spmcall_PCA_LFP_Aceto42919
%Output file name 
handles.drgb.outFileName='spmc_discriminantolfac_all_PCA_LFP_aceto42919.mat';
handles.drgb.outPathName='C:\DATA\Daniel\olfactory paper\discriminant\fwd\';

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
% handles.drgbchoices.which_electrodes=[5:12]; %Hippocampus
% handles.drgbchoices.which_electrodes=[1:4,13:16]; %Prefrontal

%Which time window do you want to be evaluated?
handles.drgbchoices.timeStart=-2;
handles.drgbchoices.timeEnd=5;

%Reference LFP
handles.drgbchoices.subtractRef=1;


%which discriminant analysis
%1 Perceptron for power LFP
%2 Linear discriminant analysis for power LFP
%10 LDA for PAC wavelet power
%11 PCA for PAC wavelet power
handles.drgbchoices.which_discriminant=[10 11];

%Percent windows
handles.drgbchoices.percent_windows=[80 100;
    45 65];
handles.drgbchoices.per_lab = {'proficient','naive'};


%First file to process
handles.drgb.first_file = 1;

%Which jt_times_ files do you want to process?
handles.drgbchoices.no_files=28;

for ii=1:handles.drgbchoices.no_files
    handles.drgbchoices.epoch(ii)=1;
end

handles.drgbchoices.epoch_names{1}=' ';


%aceto control
handles.drgbchoices.MouseName{2}='261366a';
handles.drgbchoices.group_no(1:2)=1;
handles.drgbchoices.mouse_no(1:2)=2;
handles.drgbchoices.session_no(1:2)=[1:2];
handles.drgbchoices.PathName{1}='C:\DATA\Daniel\olfactory paper\261366\aceto\';
handles.drgbchoices.FileName{1}='jt_times_72517261366acetoethylbennl.mat';
handles.drgbchoices.PathName{2}='C:\DATA\Daniel\olfactory paper\261366\aceto\';
handles.drgbchoices.FileName{2}='jt_times_71217261366acetoethylbennl.mat';

handles.drgbchoices.MouseName{3}='261367a';
handles.drgbchoices.group_no(3:4)=1;
handles.drgbchoices.mouse_no(3:4)=3;
handles.drgbchoices.session_no(3:4)=[1:2];
handles.drgbchoices.PathName{3}='C:\DATA\Daniel\olfactory paper\261367\aceto\';
handles.drgbchoices.FileName{3}='jt_times_72517261367acetoethylbennl.mat';
handles.drgbchoices.PathName{4}='C:\DATA\Daniel\olfactory paper\261367\aceto\';
handles.drgbchoices.FileName{4}='jt_times_72417261367acetoethylben_170724_172207.mat';

handles.drgbchoices.MouseName{1}='261365a';
handles.drgbchoices.group_no(5:6)=1;
handles.drgbchoices.mouse_no(5:6)=1;
handles.drgbchoices.session_no(5:6)=[1:2];
handles.drgbchoices.PathName{5}='C:\DATA\Daniel\olfactory paper\261365\aceto\';
handles.drgbchoices.FileName{5}='jt_times_73017261365acetoethylbennl.mat';
handles.drgbchoices.PathName{6}='C:\DATA\Daniel\olfactory paper\261365\aceto\';
handles.drgbchoices.FileName{6}='jt_times_71117261365acetoethylbennl.mat';

handles.drgbchoices.MouseName{4}='261368a';
handles.drgbchoices.group_no(7:8)=1;
handles.drgbchoices.mouse_no(7:8)=4;
handles.drgbchoices.session_no(7:8)=[1:2];
handles.drgbchoices.PathName{7}='C:\DATA\Daniel\olfactory paper\261368\aceto\';
handles.drgbchoices.FileName{7}='jt_times_72817261368acetoethylbennl.mat';
handles.drgbchoices.PathName{8}='C:\DATA\Daniel\olfactory paper\261368\aceto\';
handles.drgbchoices.FileName{8}='jt_times_71217261368acetoethylbennl.mat';

handles.drgbchoices.MouseName{13}='159868a';
handles.drgbchoices.group_no(9:10)=1;
handles.drgbchoices.mouse_no(9:10)=13;
handles.drgbchoices.session_no(9:10)=[1:2];
handles.drgbchoices.PathName{9}='C:\DATA\Daniel\olfactory paper\159868\';
handles.drgbchoices.FileName{9}='jt_times_103116male159868acetoethylbennl.mat';
handles.drgbchoices.PathName{10}='C:\DATA\Daniel\olfactory paper\159868\';
handles.drgbchoices.FileName{10}='jt_times_101216male159868acetoethylben.mat';


handles.drgbchoices.MouseName{5}='262641a';
handles.drgbchoices.group_no(11:12)=1;
handles.drgbchoices.mouse_no(11:12)=5;
handles.drgbchoices.session_no(11:12)=[1:2];
handles.drgbchoices.PathName{11}='C:\DATA\Daniel\olfactory paper\262641\aceto\';
handles.drgbchoices.FileName{11}='jt_times_8117262641acetoethylbennl.mat';
handles.drgbchoices.PathName{12}='C:\DATA\Daniel\olfactory paper\262641\aceto\';
handles.drgbchoices.FileName{12}='jt_times_72016262641acetoethylben.mat';

%aceto experimental
handles.drgbchoices.MouseName{14}='159863a';
handles.drgbchoices.group_no(13:14)=1;
handles.drgbchoices.mouse_no(13:14)=14;
handles.drgbchoices.session_no(13:14)=[1:2];
handles.drgbchoices.PathName{13}='C:\DATA\Daniel\olfactory paper\159863\';
handles.drgbchoices.FileName{13}='jt_times_103016female159863acetoethylbennl.mat';
handles.drgbchoices.PathName{14}='C:\DATA\Daniel\olfactory paper\159863\';
handles.drgbchoices.FileName{14}='jt_times_102416female1598631peracetoethylben.mat';

handles.drgbchoices.MouseName{15}='159865a';
handles.drgbchoices.group_no(15:16)=1;
handles.drgbchoices.mouse_no(15:16)=15;
handles.drgbchoices.session_no(15:16)=[1:2];
handles.drgbchoices.PathName{15}='C:\DATA\Daniel\olfactory paper\159865\';
handles.drgbchoices.FileName{15}='jt_times_103116male159865acetoethylbennl.mat';
handles.drgbchoices.PathName{16}='C:\DATA\Daniel\olfactory paper\159865\';
handles.drgbchoices.FileName{16}='jt_times_101216male1598651peracetoethylben.mat';

handles.drgbchoices.MouseName{16}='159866a';
handles.drgbchoices.group_no(17:18)=1;
handles.drgbchoices.mouse_no(17:18)=16;
handles.drgbchoices.session_no(17:18)=[1:2];
handles.drgbchoices.PathName{17}='C:\DATA\Daniel\olfactory paper\159866\';
handles.drgbchoices.FileName{17}='jt_times_102816male159866acetoethylben.mat';
handles.drgbchoices.PathName{18}='C:\DATA\Daniel\olfactory paper\159866\';
handles.drgbchoices.FileName{18}='jt_times_101316male159866acetoethylben.mat';

handles.drgbchoices.MouseName{17}='159867a';
handles.drgbchoices.group_no(19:20)=1;
handles.drgbchoices.mouse_no(19:20)=17;
handles.drgbchoices.session_no(19:20)=[1:2];
handles.drgbchoices.PathName{19}='C:\DATA\Daniel\olfactory paper\159867\';
handles.drgbchoices.FileName{19}='jt_times_11216male159867acetoethylbennl.mat';
handles.drgbchoices.PathName{20}='C:\DATA\Daniel\olfactory paper\159867\';
handles.drgbchoices.FileName{20}='jt_times_101216male1598671peracetoethylben.mat';

handles.drgbchoices.MouseName{18}='159869a';
handles.drgbchoices.group_no(21:22)=1;
handles.drgbchoices.mouse_no(21:22)=18;
handles.drgbchoices.session_no(21:22)=[1:2];
handles.drgbchoices.PathName{21}='C:\DATA\Daniel\olfactory paper\159869\';
handles.drgbchoices.FileName{21}='jt_times_11116female159869acetoethylbennl.mat';
handles.drgbchoices.PathName{22}='C:\DATA\Daniel\olfactory paper\159869\';
handles.drgbchoices.FileName{22}='jt_times_102916female159869acetoethylben_stch.mat';

handles.drgbchoices.MouseName{6}='254071a';
handles.drgbchoices.group_no(23:24)=1;
handles.drgbchoices.mouse_no(23:24)=6;
handles.drgbchoices.session_no(23:24)=[1:2];
handles.drgbchoices.PathName{23}='C:\DATA\Daniel\olfactory paper\254071\aceto\';
handles.drgbchoices.FileName{23}='jt_times_72517254071acetoethylbennl.mat';
handles.drgbchoices.PathName{24}='C:\DATA\Daniel\olfactory paper\254071\aceto\';
handles.drgbchoices.FileName{24}='jt_times_71117254071acetoethylbennl_170711_210112.mat';


handles.drgbchoices.MouseName{7}='261032a';
handles.drgbchoices.group_no(25:26)=1;
handles.drgbchoices.mouse_no(25:26)=7;
handles.drgbchoices.session_no(25:26)=[1:2];
handles.drgbchoices.PathName{25}='C:\DATA\Daniel\olfactory paper\261032\aceto\';
handles.drgbchoices.FileName{25}='jt_times_71117261032acetoethylbennl.mat';
handles.drgbchoices.PathName{26}='C:\DATA\Daniel\olfactory paper\261032\aceto\';
handles.drgbchoices.FileName{26}='jt_times_7917261032acetoethylbennl_170709_101250.mat';


handles.drgbchoices.MouseName{19}='254072a';
handles.drgbchoices.group_no(27:28)=1;
handles.drgbchoices.mouse_no(27:28)=19;
handles.drgbchoices.session_no(27:28)=[1:2];
handles.drgbchoices.PathName{27}='C:\DATA\Daniel\olfactory paper\254072\';
handles.drgbchoices.FileName{27}='jt_times_72017254072acetoethylbennl.mat';
handles.drgbchoices.PathName{28}='C:\DATA\Daniel\olfactory paper\254072\';
handles.drgbchoices.FileName{28}='jt_times_71717254072acetoethylbennl.mat';
% 
% 
% 
% %laser control aceto
% handles.drgbchoices.MouseName{13}='159868La';
% handles.drgbchoices.group_no(81)=3;
% handles.drgbchoices.mouse_no(81)=13;
% handles.drgbchoices.session_no(81)=[1];
% handles.drgbchoices.PathName{81}='C:\DATA\Daniel\olfactory paper\159868\';
% handles.drgbchoices.FileName{81}='jt_times_103116male159868acetoethlbenl.mat';
% 
% handles.drgbchoices.MouseName{2}='261366La';
% handles.drgbchoices.group_no(82)=3;
% handles.drgbchoices.mouse_no(82)=2;
% handles.drgbchoices.session_no(82)=[1];
% handles.drgbchoices.PathName{82}='C:\DATA\Daniel\olfactory paper\261366\laser\';
% handles.drgbchoices.FileName{82}='jt_times_72517261366acetoethylbenl_stch.mat';
% 
% handles.drgbchoices.MouseName{3}='261367La';
% handles.drgbchoices.group_no(83)=3;
% handles.drgbchoices.mouse_no(83)=3;
% handles.drgbchoices.session_no(83)=[1];
% handles.drgbchoices.PathName{83}='C:\DATA\Daniel\olfactory paper\261367\laser\';
% handles.drgbchoices.FileName{83}='jt_times_72517261367acetoethylbenl.mat';
% 
% handles.drgbchoices.MouseName{5}='262641La';
% handles.drgbchoices.group_no(84)=3;
% handles.drgbchoices.mouse_no(84)=5;
% handles.drgbchoices.session_no(84)=[1];
% handles.drgbchoices.PathName{84}='C:\DATA\Daniel\olfactory paper\262641\laser\';
% handles.drgbchoices.FileName{84}='jt_times_8117262641acetoethylbenl.mat';
% 
% handles.drgbchoices.MouseName{1}='261365La';
% handles.drgbchoices.group_no(85)=3;
% handles.drgbchoices.mouse_no(85)=1;
% handles.drgbchoices.session_no(85)=[1];
% handles.drgbchoices.PathName{85}='C:\DATA\Daniel\olfactory paper\261365\laser\';
% handles.drgbchoices.FileName{85}='jt_times_73017261365acetoethylbenl.mat';
% 
% handles.drgbchoices.MouseName{4}='261368La';
% handles.drgbchoices.group_no(86)=3;
% handles.drgbchoices.mouse_no(86)=4;
% handles.drgbchoices.session_no(86)=[1];
% handles.drgbchoices.PathName{86}='C:\DATA\Daniel\olfactory paper\261368\laser\';
% handles.drgbchoices.FileName{86}='jt_times_72817261368acetoethylbenl_stch.mat';
% 
% 
% %laser experimental aceto
% handles.drgbchoices.MouseName{14}='159863La';
% handles.drgbchoices.group_no(87)=3;
% handles.drgbchoices.mouse_no(87)=14;
% handles.drgbchoices.session_no(87)=[1];
% handles.drgbchoices.PathName{87}='C:\DATA\Daniel\olfactory paper\159863\';
% handles.drgbchoices.FileName{87}='jt_times_103016female159863acetoethlbenl.mat';
% 
% handles.drgbchoices.MouseName{15}='159865La';
% handles.drgbchoices.group_no(88)=3;
% handles.drgbchoices.mouse_no(88)=15;
% handles.drgbchoices.session_no(88)=[1];
% handles.drgbchoices.PathName{88}='C:\DATA\Daniel\olfactory paper\159865\';
% handles.drgbchoices.FileName{88}='jt_times_103116male159865acetoethlbenl.mat';
% 
% handles.drgbchoices.MouseName{16}='159866La';
% handles.drgbchoices.group_no(89)=3;
% handles.drgbchoices.mouse_no(89)=16;
% handles.drgbchoices.session_no(89)=[1];
% handles.drgbchoices.PathName{89}='C:\DATA\Daniel\olfactory paper\159866\';
% handles.drgbchoices.FileName{89}='jt_times_10286male1598661acetoethylbenl.mat';
% 
% handles.drgbchoices.MouseName{17}='159867La';
% handles.drgbchoices.group_no(90)=3;
% handles.drgbchoices.mouse_no(90)=17;
% handles.drgbchoices.session_no(90)=[1];
% handles.drgbchoices.PathName{90}='C:\DATA\Daniel\olfactory paper\159867\';
% handles.drgbchoices.FileName{90}='jt_times_11216male159867acetoethlbenl.mat';
% 
% handles.drgbchoices.MouseName{18}='159869La';
% handles.drgbchoices.group_no(91)=3;
% handles.drgbchoices.mouse_no(91)=18;
% handles.drgbchoices.session_no(91)=[1];
% handles.drgbchoices.PathName{91}='C:\DATA\Daniel\olfactory paper\159869\';
% handles.drgbchoices.FileName{91}='jt_times_11116female159869acetoethylbenl.mat';
% 
% handles.drgbchoices.MouseName{6}='254071La';
% handles.drgbchoices.group_no(92)=3;
% handles.drgbchoices.mouse_no(92)=6;
% handles.drgbchoices.session_no(92)=[1];
% handles.drgbchoices.PathName{92}='C:\DATA\Daniel\olfactory paper\254071\laser\';
% handles.drgbchoices.FileName{92}='jt_times_72517254071acetoethylbenl.mat';
% 
% handles.drgbchoices.MouseName{7}='261032La';
% handles.drgbchoices.group_no(93)=3;
% handles.drgbchoices.mouse_no(93)=7;
% handles.drgbchoices.session_no(93)=[1];
% handles.drgbchoices.PathName{93}='C:\DATA\Daniel\olfactory paper\261032\laser\';
% handles.drgbchoices.FileName{93}='jt_times_71117261032acetoethylbenl.mat';
% 
% 
% handles.drgbchoices.MouseName{19}='254072La';
% handles.drgbchoices.group_no(94)=3;
% handles.drgbchoices.mouse_no(94)=19;
% handles.drgbchoices.session_no(94)=[1];
% handles.drgbchoices.PathName{94}='C:\DATA\Daniel\olfactory paper\254072\';
% handles.drgbchoices.FileName{94}='jt_times_72017254072acetoethylbenl.mat';
% 
% %iso control
% handles.drgbchoices.MouseName{1}='261365I';
% handles.drgbchoices.group_no(95:96)=1;
% handles.drgbchoices.mouse_no(95:96)=1;
% handles.drgbchoices.session_no(95:96)=[1:2];
% handles.drgbchoices.PathName{95}='C:\DATA\Daniel\olfactory paper\261365\iso\';
% handles.drgbchoices.FileName{95}='jt_times_71017261365isominnl_stch.mat';
% handles.drgbchoices.PathName{96}='C:\DATA\Daniel\olfactory paper\261365\iso\';
% handles.drgbchoices.FileName{96}='jt_times_7717261365isominnl_stch.mat';
% 
% handles.drgbchoices.MouseName{2}='261366I';
% handles.drgbchoices.group_no(97:98)=1;
% handles.drgbchoices.mouse_no(97:98)=2;
% handles.drgbchoices.session_no(97:98)=[1:2];
% handles.drgbchoices.PathName{97}='C:\DATA\Daniel\olfactory paper\261366\iso\';
% handles.drgbchoices.FileName{97}='jt_times_71117261366isominnl.mat';
% handles.drgbchoices.PathName{98}='C:\DATA\Daniel\olfactory paper\261366\iso\';
% handles.drgbchoices.FileName{98}='jt_times_7717261366isominnl.mat';
% 
% 
% 
% handles.drgbchoices.MouseName{3}='261367I';
% handles.drgbchoices.group_no(99:100)=1;
% handles.drgbchoices.mouse_no(99:100)=3;
% handles.drgbchoices.session_no(99:100)=[1:2];
% handles.drgbchoices.PathName{99}='C:\DATA\Daniel\olfactory paper\261367\iso\';
% handles.drgbchoices.FileName{99}='jt_times_71117261367isominnl.mat';
% handles.drgbchoices.PathName{100}='C:\DATA\Daniel\olfactory paper\261367\iso\';
% handles.drgbchoices.FileName{100}='jt_times_7617261367isominnl.mat';
% 
% handles.drgbchoices.MouseName{4}='261368I';
% handles.drgbchoices.group_no(101:102)=1;
% handles.drgbchoices.mouse_no(101:102)=4;
% handles.drgbchoices.session_no(101:102)=[1:2];
% handles.drgbchoices.PathName{101}='C:\DATA\Daniel\olfactory paper\261368\iso\';
% handles.drgbchoices.FileName{101}='jt_times_71117261368isominnl.mat';
% handles.drgbchoices.PathName{102}='C:\DATA\Daniel\olfactory paper\261368\iso\';
% handles.drgbchoices.FileName{102}='jt_times_7517261368isominnl.mat';
% 
% 
% handles.drgbchoices.MouseName{5}='262641I';
% handles.drgbchoices.group_no(103:104)=1;
% handles.drgbchoices.mouse_no(103:104)=5;
% handles.drgbchoices.session_no(103:104)=[1:2];
% handles.drgbchoices.PathName{103}='C:\DATA\Daniel\olfactory paper\262641\iso\';
% handles.drgbchoices.FileName{103}='jt_times_71717262641isominnl.mat';
% handles.drgbchoices.PathName{104}='C:\DATA\Daniel\olfactory paper\262641\iso\';
% handles.drgbchoices.FileName{104}='jt_times_71517262641isominnl_170715_203504.mat';
% 
% %iso experimental
% handles.drgbchoices.MouseName{6}='254071I';
% handles.drgbchoices.group_no(105:106)=1;
% handles.drgbchoices.mouse_no(105:106)=6;
% handles.drgbchoices.session_no(105:106)=[1:2];
% handles.drgbchoices.PathName{105}='C:\DATA\Daniel\olfactory paper\254071\iso\';
% handles.drgbchoices.FileName{105}='jt_times_71017254071isominnl.mat';
% handles.drgbchoices.PathName{106}='C:\DATA\Daniel\olfactory paper\254071\iso\';
% handles.drgbchoices.FileName{106}='jt_times_7617254071isominnl_stch.mat';
% 
% handles.drgbchoices.MouseName{7}='261032I';
% handles.drgbchoices.group_no(107:108)=1;
% handles.drgbchoices.mouse_no(107:108)=7;
% handles.drgbchoices.session_no(107:108)=[1:2];
% handles.drgbchoices.PathName{107}='C:\DATA\Daniel\olfactory paper\261032\iso\';
% handles.drgbchoices.FileName{107}='jt_times_7817261032isominnl.mat';
% handles.drgbchoices.PathName{108}='C:\DATA\Daniel\olfactory paper\261032\iso\';
% handles.drgbchoices.FileName{108}='jt_times_7617261032isominnl.mat';
% 
% handles.drgbchoices.MouseName{19}='254072I';
% handles.drgbchoices.group_no(109:110)=1;
% handles.drgbchoices.mouse_no(109:110)=19;
% handles.drgbchoices.session_no(109:110)=[1:2];
% handles.drgbchoices.PathName{109}='C:\DATA\Daniel\olfactory paper\254072\';
% handles.drgbchoices.FileName{109}='jt_times_71617254072isominnl.mat';
% handles.drgbchoices.PathName{110}='C:\DATA\Daniel\olfactory paper\254072\';
% handles.drgbchoices.FileName{110}='jt_times_71517254072isominnl.mat';
% 
% handles.drgbchoices.MouseName{9}='261036I';
% handles.drgbchoices.group_no(111:112)=1;
% handles.drgbchoices.mouse_no(111:112)=9;
% handles.drgbchoices.session_no(111:112)=[1:2];
% handles.drgbchoices.PathName{111}='C:\DATA\Daniel\olfactory paper\261036\';
% handles.drgbchoices.FileName{111}='jt_times_71117261036isominnl.mat';
% handles.drgbchoices.PathName{112}='C:\DATA\Daniel\olfactory paper\261036\';
% handles.drgbchoices.FileName{112}='jt_times_7617261036isominnl_stch.mat';
% 
% 
% %laser control iso
% handles.drgbchoices.MouseName{1}='261365LI';
% handles.drgbchoices.group_no(113)=3;
% handles.drgbchoices.mouse_no(113)=1;
% handles.drgbchoices.session_no(113)=[1];
% handles.drgbchoices.PathName{113}='C:\DATA\Daniel\olfactory paper\261365\laser\';
% handles.drgbchoices.FileName{113}='jt_times_71017261365isominl_stch.mat';
% 
% handles.drgbchoices.MouseName{2}='261366LI';
% handles.drgbchoices.group_no(114)=3;
% handles.drgbchoices.mouse_no(114)=2;
% handles.drgbchoices.session_no(114)=[1];
% handles.drgbchoices.PathName{114}='C:\DATA\Daniel\olfactory paper\261366\laser\';
% handles.drgbchoices.FileName{114}='jt_times_71117261366isominl.mat';
% 
% handles.drgbchoices.MouseName{3}='261367LI';
% handles.drgbchoices.group_no(115)=3;
% handles.drgbchoices.mouse_no(115)=3;
% handles.drgbchoices.session_no(115)=[1];
% handles.drgbchoices.PathName{115}='C:\DATA\Daniel\olfactory paper\261367\laser\';
% handles.drgbchoices.FileName{115}='jt_times_71117261367isominl.mat';
% 
% handles.drgbchoices.MouseName{4}='261368LI';
% handles.drgbchoices.group_no(116)=3;
% handles.drgbchoices.mouse_no(116)=4;
% handles.drgbchoices.session_no(116)=[1];
% handles.drgbchoices.PathName{116}='C:\DATA\Daniel\olfactory paper\261368\laser\';
% handles.drgbchoices.FileName{116}='jt_times_71117261368isominl.mat';
% 
% handles.drgbchoices.MouseName{5}='262641LI';
% handles.drgbchoices.group_no(117)=3;
% handles.drgbchoices.mouse_no(117)=5;
% handles.drgbchoices.session_no(117)=[1];
% handles.drgbchoices.PathName{117}='C:\DATA\Daniel\olfactory paper\262641\laser\';
% handles.drgbchoices.FileName{117}='jt_times_71717262641isominl.mat';
% 
% %laser expe iso
% handles.drgbchoices.MouseName{6}='254071LI';
% handles.drgbchoices.group_no(118)=3;
% handles.drgbchoices.mouse_no(118)=6;
% handles.drgbchoices.session_no(118)=[1];
% handles.drgbchoices.PathName{118}='C:\DATA\Daniel\olfactory paper\254071\laser\';
% handles.drgbchoices.FileName{118}='jt_times_71017254071isominl.mat';
% 
% handles.drgbchoices.MouseName{7}='261032LI';
% handles.drgbchoices.group_no(119)=3;
% handles.drgbchoices.mouse_no(119)=7;
% handles.drgbchoices.session_no(119)=[1];
% handles.drgbchoices.PathName{119}='C:\DATA\Daniel\olfactory paper\261032\laser\';
% handles.drgbchoices.FileName{119}='jt_times_7817261032isominl.mat';
% 
% handles.drgbchoices.MouseName{19}='254072LI';
% handles.drgbchoices.group_no(120)=3;
% handles.drgbchoices.mouse_no(120)=19;
% handles.drgbchoices.session_no(120)=[1];
% handles.drgbchoices.PathName{120}='C:\DATA\Daniel\olfactory paper\254072\';
% handles.drgbchoices.FileName{120}='jt_times_71617254072isominl_stch.mat';
% 
% handles.drgbchoices.MouseName{9}='261036LI';
% handles.drgbchoices.group_no(121)=3;
% handles.drgbchoices.mouse_no(121)=9;
% handles.drgbchoices.session_no(121)=[1];
% handles.drgbchoices.PathName{121}='C:\DATA\Daniel\olfactory paper\261036\';
% handles.drgbchoices.FileName{121}='jt_times_71117261036isominl.mat';


%Enter the name of each group No
handles.drgbchoices.group_no_names{1}='Fwd';
handles.drgbchoices.group_no_names{2}='Rev';
% handles.drgbchoices.group_no_names{3}='laser';



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

%Wideband
handles.drgbchoices.lowF(5)=4;
handles.drgbchoices.highF(5)=95;



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





function handles=drgChR2_choices_test7252023


%run with drgCalmAn_batch_dropc_fsdz
% handles.PathName_out='/home/simoesdf/Grin6and7_spm_raw/';
% handles.FileName_out='entire_session_deepcad6_02052022.mat';

%The choices below are for drgCaImAn_pval_rev_batch
handles.suffix_out='_LFPp.mat';
handles.PathName_out='/Volumes/Diego HD/Joe/';

%The choices below are for both drgCaImAn_warp_ROIs_batch and drgCaImAn_pval_rev_batch
handles.first_file=1;

handles.genotype{1}='OMPChR2';
handles.genotype{2}='C57BL/6';

handles.electrode_labels{1}='Right CA1';
handles.electrode_labels{2}='Left CA1';
handles.electrode_labels{3}='Right OB';
handles.electrode_labels{4}='Left OB';

handles.burstLowF=1;
handles.burstHighF=100;

handles.mouse{1}='Barholomew';
handles.mouse_no(1)=1;
handles.frequency(1)=40;
handles.group(1)=1;
handles.PathName_in{1}='/Volumes/Diego HD/Joe/Bartholomew_OB_Hippo_Recordings/20230821_Bar_40Hz_230821_103911/';
handles.FileName_in{1}='jt_times_20230821_Bar_40Hz_230821_103911.mat';
handles.file(1).peakLFPNo=[2 7 10 15];

handles.mouse{2}='Barholomew';
handles.mouse_no(2)=1;
handles.frequency(2)=20;
handles.group(2)=1;
handles.PathName_in{2}='/Volumes/Diego HD/Joe/Bartholomew_OB_Hippo_Recordings/20230821_Bar_20Hz_230821_103114/';
handles.FileName_in{2}='jt_times_20230821_Bar_20Hz_230821_103114.mat';
handles.file(2).peakLFPNo=[2 7 10 15];


handles.mouse{3}='Oscar';
handles.mouse_no(3)=2;
handles.frequency(3)=40;
handles.group(3)=1;
handles.PathName_in{3}='/Volumes/Diego HD/Joe/Bartholomew_OB_Hippo_Recordings/20230821_Bar_40Hz_230821_103911/';
handles.FileName_in{3}='jt_times_20230821_Bar_40Hz_230821_103911.mat';
handles.file(3).peakLFPNo=[2 7 10 15];

handles.mouse{4}='Oscar';
handles.mouse_no(4)=1;
handles.frequency(4)=20;
handles.group(4)=1;
handles.PathName_in{4}='/Volumes/Diego HD/Joe/Bartholomew_OB_Hippo_Recordings/20230821_Bar_20Hz_230821_103114/';
handles.FileName_in{4}='jt_times_20230821_Bar_20Hz_230821_103114.mat';
handles.file(4).peakLFPNo=[2 7 10 15];

handles.no_files=length(handles.FileName_in);

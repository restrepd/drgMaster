function drgAnalysisBatchLFPCaMKIIv2(handles)

%drgAnalysisBatchLFPCaMKIIv2 displays the LFP power and PAC data processed by drgRunBatchLFPpar

%Which analysis is performed is determined by the value enterd in the
%variable which_display:
%
% 1 Oscillatory wavelet power calculated at the peak and trough of the low frequency PAC phase
% calculated from a full time window drgRunBatchLFPpar (e.g. -1.5 to 5
% sec). We calculate phase referenced power (PRP) for peak and trough and
% lick-referenced power (LRP).
% For PRP this yields the same graphs of case 24 of drgAnalysisBatchLFPCaMKII
% (below) that was used for Losacco, Ramirez-Gordillo et al., eLife 2020 Figure 3
%
% 2 Oscillatory PRP and LRP timecourses are converted to a z score to allow comparing 
% S+ vs S- using all data in all electrodes. This yields p value
% timecourses and discrimination times
%
% 3 Yields histograms for autocorrelation and cross correlation of lick and
% peak/trough times and performs z score analysis for discrimination for
% PRP for peaks chosen in relation to licks. This case also yields lick
% rate plots
%
% 19 PAC MI analysis for events (concentrations or S+/S-) for naive and proficient
% Analyzed per mouse for groups defined by the user
% Used for Figure 2 of Losacco, Ramirez-Gordillo et al., eLife 2020
%
% 20  Multiclass ROC analysis of LFP power differences for naive and proficient
% mice for different epochs (concentrations or S+ vs. S-) and different groups. Analyzed per mouse
%
% 21  Multiclass ROC analysis of coherence for naive and proficient
%  mice for different epochs (concentrations or S+ vs. S-) and different groups. Analyzed per mouse
%
% 22 ERWA analysis for LFP power
%
% 23 Oscillatory power calculated from the amplitude of the Hilbert envelope of the high
% frequency at the peak and trough of the low frequency phase. This was NOT
% used for Figure 3 of the Losacco, Ramirez-Gordillo et al., eLife 2020 paper. The code has not been vetted
%
% 24 Oscillatory wavelet power calculated at the peak and trough of the low frequency PAC phase
% This was used for Losacco, Ramirez-Gordillo et al., eLife 2020 Figure 3
%

%% Read the BatchParameters
[parsFileName,parsPathName] = uigetfile({'drgLFPBatchAnalPars*.m'},'Select the .m file with all the parameters for LFP batch analysis');
fprintf(1, ['\ndrgAnalysisBatchLFP run for ' parsFileName '\n\n']);

addpath(parsPathName)
eval(['handles_pars=' parsFileName(1:end-2) ';'])
handles.parsFileName=parsFileName;
handles.parsPathName=parsPathName;

if isfield(handles_pars,'winNo')
    winNo=handles_pars.winNo;
end
if isfield(handles_pars,'which_display')
    which_display=handles_pars.which_display;
end
if isfield(handles_pars,'eventType')
    eventType=handles_pars.eventType;
end
if isfield(handles_pars,'evTypeLabels')
    evTypeLabels=handles_pars.evTypeLabels;
end
if isfield(handles_pars,'file_pairs')
    file_pairs=handles_pars.file_pairs;
end
if isfield(handles_pars,'trials_to_process')
    trials_to_process=handles_pars.trials_to_process;
end
if isfield(handles_pars,'min_trials_per_event')
    min_trials_per_event=handles_pars.min_trials_per_event;
end
if isfield(handles_pars,'shift_time')
    shift_time=handles_pars.shift_time;
end
if isfield(handles_pars,'shift_from_event')
    shift_from_event=handles_pars.shift_from_event;
end
if isfield(handles_pars,'grpost')
    grpost=handles_pars.grpost;
end
if isfield(handles_pars,'file_label')
    file_label=handles_pars.file_label;
end
if isfield(handles_pars,'front_mask')
    front_mask=handles_pars.front_mask;
end
if isfield(handles_pars,'output_suffix')
    output_suffix=handles_pars.output_suffix;
end
if isfield(handles_pars,'percent_windows')
    percent_windows=handles_pars.percent_windows;
end
if isfield(handles_pars,'delta_t_ii')
    delta_t_ii=handles_pars.delta_t_ii;
end
if isfield(handles_pars,'which_electrodes')
    which_electrodes=handles_pars.which_electrodes;
end
if isfield(handles_pars,'files')
    files=handles_pars.files;
end
if isfield(handles_pars,'concs2')
    concs2=handles_pars.concs2;
end
if isfield(handles_pars,'concs')
    concs=handles_pars.concs;
end
if isfield(handles_pars,'per_lab')
    per_lab=handles_pars.per_lab;
end

if ~isfield(handles_pars,'no_bandwidths')
    no_bandwidths=4;
    low_freq=[6 15 35 65];
    high_freq=[14 30 55 95];
else
    no_bandwidths=handles_pars.no_bandwidths;
    low_freq=handles_pars.low_freq;
    high_freq=handles_pars.high_freq;
end

if ~isfield(handles_pars,'freq_names')
    freq_names={'Theta','Beta','Low gamma','High gamma'};
else
    freq_names=handles_pars.freq_names;
end
if ~isfield(handles_pars,'no_pacii')
    no_pacii=3;
else
    no_pacii=handles_pars.no_pacii;
end

refWin=handles_pars.refWin;


%% The code processing pairwise batch LFP starts here

close all
warning('off')


%Bandwidths

if exist('no_bandwidths')==0
    no_bandwidths=4;
    low_freq=[6 15 35 65];
    high_freq=[12 30 55 95];
    freq_names={'Theta','Beta','Low gamma','High gamma'};
end


event1=eventType(1);
event2=eventType(2);


%Ask user for the drgb output .mat file and load those data
[handles.drgb.outFileName,handles.PathName] = uigetfile('*.mat','Select the drgb output file');
% [handles.drgb.outFileName,handles.PathName] = uigetfile('*');
% handles.PathName=handles_pars.PathName;
% handles.drgb.outFileName=handles_pars.outFileName;
load([handles.PathName handles.drgb.outFileName])

if (~isfield(handles_pars,'files'))||(isempty(handles_pars.files))
    files=[1:handles_drgb.drgbchoices.no_files];
end

if isfield(handles_pars,'groupNo')
    handles_drgb.drgbchoices.group_no=handles_pars.groupNo;
end

fprintf(1, ['\ndrgAnalysisBatchLFP run for ' handles.drgb.outFileName '\nwhich_display= %d\n\n'],which_display);

switch which_display
    case 21
        frequency=handles_drgb.drgb.lfpevpair(1).f_coh;
    case 22
        frequency=handles_drgb.drgb.lfpevpair(1).wave_fERWA;
    case {1,2,3,24}
        %Do nothing
    otherwise
        if isfield(handles_drgb.drgb,'freq_for_LFPpower')
            frequency=handles_drgb.drgb.freq_for_LFPpower;
        else
            frequency=handles_drgb.drgb.file(1).freq_for_LFPpower;
        end
end

if isfield(handles_drgb.drgb.file,'eventlabels')
    if ~isempty(handles_drgb.drgb.file(1).eventlabels)
        %Overwrite the event labels
        for ii=1:length(eventType)
            evTypeLabels{ii}=handles_drgb.drgb.file(1).eventlabels{handles_drgb.drgbchoices.evTypeNos(eventType(ii))};
        end
    end
end

%These are the colors for the different lines

these_colors{1}='b';
these_colors{2}='r';
these_colors{3}='m';
these_colors{8}='g';
these_colors{5}='y';
these_colors{6}='k';
these_colors{7}='c';
these_colors{4}='k';

these_light_colors{1}=[0.7 0.7 1];
these_light_colors{2}=[0.7 0.7 1];

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
    window_per_lfp(lfpodNo)=handles_drgb.drgb.lfpevpair(lfpodNo).timeWindow;
    switch which_display
        case {21}
            elec_pair_No(lfpodNo)=handles_drgb.drgb.lfpevpair(lfpodNo).elec_pair_No;
            elec1(lfpodNo)=handles_drgb.drgb.lfpevpair(lfpodNo).elec1;
            elec2(lfpodNo)=handles_drgb.drgb.lfpevpair(lfpodNo).elec2;
        otherwise
            elec_per_lfp(lfpodNo)=handles_drgb.drgb.lfpevpair(lfpodNo).elecNo;
    end
    
end

switch which_display
    
    case 1
        %1 Oscillatory wavelet power timecourse at the peak and trough of the PAC
        
        mean_PACpower_No_per_mouse=0;
        mean_deltaLickFreq_No_per_mouse=0;
        
        mean_PACpower_No=0;
        mean_peakPACpower=[];
        mean_troughPACpower=[];
        mean_lickPACpower=[];
        mean_PACpower_perii=[];
        mean_PACpower_evNo=[];
        mean_PACpower_pacii=[];
        mean_PACpower_fileNo=[];
        per_session_group_no=[];
        mean_VL=[];
        mean_VA=[];
        mean_PA=[];
        
        handles_out=[];
        handles_out.PRP_ii=0;
        handles_out.LickF_ii=0;
        handles_out.AUC_ii=0;
        handles_out.AUClick_ii=0;
        
        prof_naive_leg{1}='Proficient';
        prof_naive_leg{2}='Naive';
        
        %Initialize ROC
        no_ROCs=0;
        ROCfileNo=[];
        ROCelec=[];
        ROCgroups=[];
        ROCmouse=[];
        ROCpacii=[];
        ROCper_ii=[];
        ROCEvNo1=[];
        ROCEvNo2=[];
        ROC_between=[];
        auROCpeak=[];
        auROCtrough=[];
        p_valROCpeak=[];
        p_vals_ROCpeak=[];
        p_valROCtrough=[];
        p_vals_ROCtrough=[];
        
        no_lick_ROCs=0;
        ROClick_fileNo=[];
        ROClick_elec=[];
        ROClick_groups=[];
        ROClick_mouse=[];
        ROClick_pacii=[];
        ROClick_per_ii=[];
        ROClick_EvNo1=[];
        ROClick_EvNo2=[];
        ROClick_between=[];
        auROClick=[];
        p_valROClick=[];
        p_vals_ROClick=[];
        
        fprintf(1, ['PAC power analysis using wavelet power for Justin''s paper\n\n'])
        p_vals=[];
        no_files=max(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        
        szpc=size(percent_windows);
        for per_ii=1:szpc(1)
            
            
            
            for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
                mouse_has_files=0;
                for elec=1:16
                    if sum(which_electrodes==elec)>0
                        
                        theseEvNos_thisMouse_thisElec=[];
                        for evNo=1:length(eventType)
                            for pacii=1:no_pacii
                                for group_no=1:max(handles_drgb.drgbchoices.group_no)
                                    theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv=0;
                                    theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).nolickEv=0;
                                end
                            end
                        end
                        
                        
                        for fileNo=1:no_files
                            %If this file is in the list of files the user wants to process in drgAnalysisBatchLFP continue
                            if sum(files==fileNo)>0
                                if handles_drgb.drgbchoices.mouse_no(fileNo)==mouseNo
                                    lfpodNo=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                    %                                     lfpodRefNo=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                                    
                                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo)))
                                        
                                        
                                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo).PAC))
                                            
                                            percent_mask=[];
                                            percent_mask=logical((handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).perCorrPAC>=percent_windows(per_ii,1))...
                                                &(handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).perCorrPAC<=percent_windows(per_ii,2)));
                                            
                                            if ~isempty(percent_mask)
                                                
                                                
                                                for evNo=1:length(eventType)
                                                    
                                                    trials_in_event_Ev=[];
                                                    trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).which_eventPAC(eventType(evNo),:)==1)&percent_mask;
                                                    
                                                    if (sum(trials_in_event_Ev)>=1)
                                                        
                                                        %Do per bandwidth analysis
                                                        for pacii=1:no_pacii
                                                            
                                                            group_no=handles_drgb.drgbchoices.group_no(fileNo);
                                                            
                                                            %Enter the PACpower
                                                            t_pac=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.t_pac;
                                                            
                                                            ref_mean_peakPACpower=zeros(1, length(trials_in_event_Ev));
                                                            analysis_mean_peakPACpower=zeros(1, length(trials_in_event_Ev));
                                                            for ww=1:length(trials_in_event_Ev)
                                                                this_p=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(ww).peakPower;
                                                                ref_mean_peakPACpower(1,ww)=mean(this_p((t_pac>=handles_pars.refWin_times(1))&(t_pac<=handles_pars.refWin_times(2))));
                                                                analysis_mean_peakPACpower(1,ww)=mean(this_p((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2))));
                                                            end
                                                            
                                                            ref_mean_troughPACpower=zeros(1, length(trials_in_event_Ev));
                                                            analysis_mean_troughPACpower=zeros(1, length(trials_in_event_Ev));
                                                            for ww=1:length(trials_in_event_Ev)
                                                                this_p=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(ww).troughPower;
                                                                ref_mean_troughPACpower(1,ww)=mean(this_p((t_pac>=handles_pars.refWin_times(1))&(t_pac<=handles_pars.refWin_times(2))));
                                                                analysis_mean_troughPACpower(1,ww)=mean(this_p((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2))));
                                                            end
                                                            
                                                            trials_in_lick_event=zeros(1,length(handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(1).PACwave.lickPower_trials));
                                                            trials_in_lick_event(1,:)=trials_in_event_Ev(handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(1).PACwave.lickPower_trials);
                                                            ref_mean_lickPACpower=zeros(1, length(trials_in_lick_event));
                                                            analysis_mean_lickPACpower=zeros(1, length(trials_in_lick_event));
                                                            for ww=1:length(trials_in_lick_event)
                                                                this_p=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(ww).lickPower;
                                                                ref_mean_lickPACpower(1,ww)=mean(this_p((t_pac>=handles_pars.refWin_times(1))&(t_pac<=handles_pars.refWin_times(2))));
                                                                analysis_mean_lickPACpower(1,ww)=mean(this_p((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2))));
                                                            end
                                                            
                                                            this_peakPACpower_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_peakPACpower_Ev=analysis_mean_peakPACpower(trials_in_event_Ev)-ref_mean_peakPACpower(trials_in_event_Ev);
                                                            this_troughPACpower_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_troughPACpower_Ev=analysis_mean_troughPACpower(trials_in_event_Ev)-ref_mean_troughPACpower(trials_in_event_Ev);
                                                            this_lickPACpower_Ev=zeros(sum(trials_in_lick_event),1);
                                                            this_lickPACpower_Ev=analysis_mean_lickPACpower(logical(trials_in_lick_event))-ref_mean_lickPACpower(logical(trials_in_lick_event));
                                                            %                                                             if (elec==which_electrodes(1))&(pacii==1)
                                                            %                                                                 this_deltaLickFreq_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            %                                                                 this_deltaLickFreq_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.mean_lick_freq(trials_in_event_Ev)-handles_drgb.drgb.lfpevpair(lfpodRefNo).PACwave(pacii).PACwave.mean_lick_freq(trials_in_event_Ev);
                                                            %                                                             end
                                                            
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_peakPACpower_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1 ...
                                                                :theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=this_peakPACpower_Ev;
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_troughPACpower_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1 ...
                                                                :theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=this_troughPACpower_Ev;
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_lickPACpower_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).nolickEv+1 ...
                                                                :theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).nolickEv+sum(trials_in_lick_event))=this_lickPACpower_Ev;
                                                            
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).whichMouse(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=mouseNo*ones(1,length(this_peakPACpower_Ev));
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev);
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).nolickEv=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).nolickEv+sum(trials_in_lick_event);
                                                            
                                                            
                                                            %                                                             if (elec==which_electrodes(1))&(pacii==1)
                                                            %                                                                  theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_deltaLickFreq_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1 ...
                                                            %                                                                 :theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=this_deltaLickFreq_Ev;
                                                            %                                                             end
                                                            
                                                            %Save per session value for peak power
                                                            mean_PACpower_No=mean_PACpower_No+1;
                                                            mean_peakPACpower(mean_PACpower_No)=mean(this_peakPACpower_Ev);
                                                            mean_troughPACpower(mean_PACpower_No)=mean(this_troughPACpower_Ev);
                                                            mean_lickPACpower(mean_PACpower_No)=mean(this_lickPACpower_Ev);
                                                            mean_PACpower_perii(mean_PACpower_No)=per_ii;
                                                            mean_PACpower_evNo(mean_PACpower_No)=evNo;
                                                            mean_PACpower_pacii(mean_PACpower_No)=pacii;
                                                            mean_PACpower_fileNo(mean_PACpower_No)=fileNo;
                                                            mean_PACpower_mouse(mean_PACpower_No)=mouseNo;
                                                            per_session_group_no(mean_PACpower_No)=handles_drgb.drgbchoices.group_no(fileNo);
                                                            
                                                            
                                                            mouse_has_files=1;
                                                            
                                                        end
                                                        
                                                        
                                                        fprintf(1, ['%d trials in event No %d succesfully processed for file No %d electrode %d\n'],sum(trials_in_event_Ev), evNo,fileNo,elec);
                                                        
                                                    else
                                                        
                                                        
                                                        fprintf(1, ['%d trials in event No %d fewer than minimum trials per event ' evTypeLabels{evNo} ' for file No %d electrode %d\n'],sum(trials_in_event_Ev), min_trials_per_event,fileNo,elec);
                                                        
                                                        
                                                    end
                                                    
                                                    
                                                end
                                                
                                            else
                                                
                                                fprintf(1, ['Empty percent_mask for file No %d electrode %d\n'],fileNo,elec);
                                                
                                            end
                                            
                                        else
                                            
                                            fprintf(1, ['Empty PAC for file No %d electrode %d\n'],fileNo,elec);
                                            
                                        end
                                        
                                        
                                    else
                                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],fileNo,elec);
                                        
                                        
                                    end
                                end %if mouseNo
                            end
                        end %fileNo
                        
                        if theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv>0
                            
                            %Calculate per mouse PAC power
                            for evNo=1:length(eventType)
                                for pacii=1:no_pacii
                                    for group_no=1:max(handles_drgb.drgbchoices.group_no)
                                        %Calculate per mouse PAC power
                                        mean_PACpower_No_per_mouse=mean_PACpower_No_per_mouse+1;
                                        this_mouse_peakPACpower=[];
                                        this_mouse_peakPACpower=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_peakPACpower_Ev;
                                        this_mouse_troughPACpower=[];
                                        this_mouse_troughPACpower=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_troughPACpower_Ev;
                                        this_mouse_lickPACpower=[];
                                        this_mouse_lickPACpower=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_lickPACpower_Ev;
                                        
                                        if ~isempty(this_mouse_peakPACpower)
                                            mean_peakPACpower_per_mouse(mean_PACpower_No_per_mouse)=mean(this_mouse_peakPACpower);
                                            mean_troughPACpower_per_mouse(mean_PACpower_No_per_mouse)=mean(this_mouse_troughPACpower);
                                        else
                                            mean_peakPACpower_per_mouse(mean_PACpower_No_per_mouse)=NaN;
                                            mean_troughPACpower_per_mouse(mean_PACpower_No_per_mouse)=NaN;
                                        end
                                        
                                        if ~isempty(this_mouse_lickPACpower)
                                            mean_lickPACpower_per_mouse(mean_PACpower_No_per_mouse)=mean(this_mouse_lickPACpower);
                                        else
                                            mean_lickPACpower_per_mouse(mean_PACpower_No_per_mouse)=NaN;
                                        end
                                        
                                        mean_PACpower_perii_per_mouse(mean_PACpower_No_per_mouse)=per_ii;
                                        mean_PACpower_evNo_per_mouse(mean_PACpower_No_per_mouse)=evNo;
                                        mean_PACpower_pacii_per_mouse(mean_PACpower_No_per_mouse)=pacii;
                                        mean_PACpower_mouseNo_per_mouse(mean_PACpower_No_per_mouse)=mouseNo;
                                        mean_PACpower_electNo_per_mouse(mean_PACpower_No_per_mouse)=elec;
                                        mean_PACpower_group_no_per_mouse(mean_PACpower_No_per_mouse)=group_no;
                                        
                                        
                                        
                                    end
                                end
                            end
                            
                            
                            
                            %Calculate per electrode ROC
                            can_calculate_auroc=1;
                            if can_calculate_auroc==1
                                for group_no=1:max(handles_drgb.drgbchoices.group_no)
                                    for evNo1=1:length(eventType)
                                        %                                         if theseEvNos(evNo1).noEv>0
                                        for evNo2=evNo1+1:length(eventType)
                                            %                                                 if theseEvNos(evNo2).noEv>0
                                            if (theseEvNos_thisMouse_thisElec(evNo1,pacii,group_no).noEv>0)&...
                                                    (theseEvNos_thisMouse_thisElec(evNo2,pacii,group_no).noEv>0)
                                                for pacii=1:no_pacii
                                                    
                                                    if (theseEvNos_thisMouse_thisElec(evNo1,pacii,group_no).noEv>=5)&...
                                                            (theseEvNos_thisMouse_thisElec(evNo2,pacii,group_no).noEv>=5)
                                                        
                                                        %Enter Ev1
                                                        this_mouse_peakPACpowerEv1=[];
                                                        this_mouse_peakPACpowerEv1=theseEvNos_thisMouse_thisElec(evNo1,pacii,group_no).this_peakPACpower_Ev;
                                                        this_mouse_troughPACpowerEv1=[];
                                                        this_mouse_troughPACpowerEv1=theseEvNos_thisMouse_thisElec(evNo1,pacii,group_no).this_troughPACpower_Ev;
                                                        
                                                        trials_in_event_Ev1=theseEvNos_thisMouse_thisElec(evNo1,pacii,group_no).noEv;
                                                        
                                                        roc_data_peak=[];
                                                        roc_data_peak(1:trials_in_event_Ev1,1)=this_mouse_peakPACpowerEv1;
                                                        roc_data_peak(1:trials_in_event_Ev1,2)=zeros(trials_in_event_Ev1,1);
                                                        roc_data_trough=[];
                                                        roc_data_trough(1:trials_in_event_Ev1,1)=this_mouse_troughPACpowerEv1;
                                                        roc_data_trough(1:trials_in_event_Ev1,2)=zeros(trials_in_event_Ev1,1);
                                                        
                                                        
                                                        %Enter Ev2
                                                        this_mouse_peakPACpowerEv2=[];
                                                        this_mouse_peakPACpowerEv2=theseEvNos_thisMouse_thisElec(evNo2,pacii,group_no).this_peakPACpower_Ev;
                                                        this_mouse_troughPACpowerEv2=[];
                                                        this_mouse_troughPACpowerEv2=theseEvNos_thisMouse_thisElec(evNo2,pacii,group_no).this_troughPACpower_Ev;
                                                        
                                                        trials_in_event_Ev2=theseEvNos_thisMouse_thisElec(evNo2,pacii,group_no).noEv;
                                                        
                                                        roc_data_peak(trials_in_event_Ev1+1:trials_in_event_Ev2+trials_in_event_Ev1,1)=this_mouse_peakPACpowerEv2;
                                                        roc_data_peak(trials_in_event_Ev1+1:trials_in_event_Ev2+trials_in_event_Ev1,2)=ones(trials_in_event_Ev2,1);
                                                        roc_data_trough(trials_in_event_Ev1+1:trials_in_event_Ev2+trials_in_event_Ev1,1)=this_mouse_troughPACpowerEv2;
                                                        roc_data_trough(trials_in_event_Ev1+1:trials_in_event_Ev2+trials_in_event_Ev1,2)=ones(trials_in_event_Ev2,1);
                                                        
                                                        %Find  per electrode ROC
                                                        
                                                        no_ROCs=no_ROCs+1;
                                                        
                                                        
                                                        %                                                         ROCfileNo(no_ROCs)=handles_drgb.drgb.lfpevpair(lfpodNo_ref).fileNo;
                                                        ROCelec(no_ROCs)=elec;
                                                        ROCgroups(no_ROCs)=group_no;
                                                        ROCmouse(no_ROCs)=mouseNo;
                                                        ROCpacii(no_ROCs)=pacii;
                                                        ROCper_ii(no_ROCs)=per_ii;
                                                        ROCEvNo1(no_ROCs)=evNo1;
                                                        ROCEvNo2(no_ROCs)=evNo2;
                                                        
                                                        if ((abs(evNo1-2)<=1)&(abs(evNo2-5)<=1))||((abs(evNo1-5)<=1)&(abs(evNo2-2)<=1))
                                                            ROC_between(no_ROCs)=1;
                                                        else
                                                            ROC_between(no_ROCs)=0;
                                                        end
                                                        
                                                        roc=[];
                                                        roc=roc_calc(roc_data_peak,0,0.05,0);
                                                        auROCpeak(no_ROCs)=roc.AUC-0.5;
                                                        p_valROCpeak(no_ROCs)=roc.p;
                                                        p_vals_ROCpeak=[p_vals_ROCpeak roc.p];
                                                        
                                                        roc=[];
                                                        roc=roc_calc(roc_data_trough,0,0.05,0);
                                                        auROCtrough(no_ROCs)=roc.AUC-0.5;
                                                        p_valROCtrough(no_ROCs)=roc.p;
                                                        p_vals_ROCtrough=[p_vals_ROCtrough roc.p];
                                                        
                                                        
                                                        
                                                        %I have this code here to plot the ROC
                                                        
                                                        show_roc=0;
                                                        if (show_roc==1)&(mouseNo==4)
                                                            %I have this code here to plot the ROC
                                                            roc=roc_calc(roc_data_peak,0,0.05,1);
                                                            
                                                            %Do the histograms
                                                            try
                                                                close(2)
                                                            catch
                                                            end
                                                            figure(2)
                                                            
                                                            hold on
                                                            
                                                            max_dB=max([max(this_mouse_peakPACpowerEv1) max(this_mouse_peakPACpowerEv2)]);
                                                            min_dB=min([min(this_mouse_peakPACpowerEv1) min(this_mouse_peakPACpowerEv2)]);
                                                            
                                                            edges=[min_dB-0.1*(max_dB-min_dB):(max_dB-min_dB)/20:max_dB+0.1*(max_dB-min_dB)];
                                                            histogram(this_mouse_peakPACpowerEv1,edges,'FaceColor','b','EdgeColor','b')
                                                            histogram(this_mouse_peakPACpowerEv2,edges,'FaceColor','r','EdgeColor','r')
                                                            xlabel('delta power dB')
                                                            title(['Histogram for concentrations ' num2str(concs2(evNo1)) ' and ' num2str(concs2(evNo2))])
                                                            pffft=1;
                                                        end
                                                        
                                                    end
                                                end
                                                %
                                            end
                                        end
                                    end
                                    %                                         end
                                end
                                
                                
                                %Calculate ROC for licks if appropriate
                                for group_no=1:max(handles_drgb.drgbchoices.group_no)
                                    for evNo1=1:length(eventType)
                                        %                                         if theseEvNos(evNo1).noEv>0
                                        for evNo2=evNo1+1:length(eventType)
                                            %                                                 if theseEvNos(evNo2).noEv>0
                                            if (theseEvNos_thisMouse_thisElec(evNo1,pacii,group_no).nolickEv>0)&...
                                                    (theseEvNos_thisMouse_thisElec(evNo2,pacii,group_no).nolickEv>0)
                                                for pacii=1:no_pacii
                                                    
                                                    if (theseEvNos_thisMouse_thisElec(evNo1,pacii,group_no).nolickEv>=5)&...
                                                            (theseEvNos_thisMouse_thisElec(evNo2,pacii,group_no).nolickEv>=5)
                                                        
                                                        %Enter Ev1
                                                        this_mouse_lickPACpowerEv1=[];
                                                        this_mouse_lickPACpowerEv1=theseEvNos_thisMouse_thisElec(evNo1,pacii,group_no).this_lickPACpower_Ev;
                                                        trials_in_event_Ev1=theseEvNos_thisMouse_thisElec(evNo1,pacii,group_no).nolickEv;
                                                        
                                                        roc_data_lick=[];
                                                        roc_data_lick(1:trials_in_event_Ev1,1)=this_mouse_lickPACpowerEv1;
                                                        roc_data_lick(1:trials_in_event_Ev1,2)=zeros(trials_in_event_Ev1,1);
                                                        
                                                        %Enter Ev2
                                                        this_mouse_lickPACpowerEv2=[];
                                                        this_mouse_lickPACpowerEv2=theseEvNos_thisMouse_thisElec(evNo2,pacii,group_no).this_lickPACpower_Ev;
                                                        trials_in_event_Ev2=theseEvNos_thisMouse_thisElec(evNo2,pacii,group_no).nolickEv;
                                                        
                                                        
                                                        roc_data_lick(trials_in_event_Ev1+1:trials_in_event_Ev2+trials_in_event_Ev1,1)=this_mouse_lickPACpowerEv2;
                                                        roc_data_lick(trials_in_event_Ev1+1:trials_in_event_Ev2+trials_in_event_Ev1,2)=ones(trials_in_event_Ev2,1);
                                                        
                                                        %Find  per electrode ROC
                                                        no_lick_ROCs=no_lick_ROCs+1;
                                                        
                                                        
                                                        %                                                         ROCfileNo(no_ROCs)=handles_drgb.drgb.lfpevpair(lfpodNo_ref).fileNo;
                                                        ROClick_elec(no_lick_ROCs)=elec;
                                                        ROClick_groups(no_lick_ROCs)=group_no;
                                                        ROClick_mouse(no_lick_ROCs)=mouseNo;
                                                        ROClick_pacii(no_lick_ROCs)=pacii;
                                                        ROClick_per_ii(no_lick_ROCs)=per_ii;
                                                        ROClick_EvNo1(no_lick_ROCs)=evNo1;
                                                        ROClick_EvNo2(no_lick_ROCs)=evNo2;
                                                        
                                                        if ((abs(evNo1-2)<=1)&(abs(evNo2-5)<=1))||((abs(evNo1-5)<=1)&(abs(evNo2-2)<=1))
                                                            ROC_between(no_lick_ROCs)=1;
                                                        else
                                                            ROC_between(no_lick_ROCs)=0;
                                                        end
                                                        
                                                        roc=[];
                                                        roc=roc_calc(roc_data_lick,0,0.05,0);
                                                        auROClick(no_lick_ROCs)=roc.AUC-0.5;
                                                        p_valROClick(no_lick_ROCs)=roc.p;
                                                        p_vals_ROClick=[p_vals_ROClick roc.p];
                                                        
                                                        %I have this code here to plot the ROC
                                                        
                                                        show_roc=0;
                                                        if (show_roc==1)&(mouseNo==4)
                                                            %I have this code here to plot the ROC
                                                            roc=roc_calc(roc_data_peak,0,0.05,1);
                                                            
                                                            %Do the histograms
                                                            try
                                                                close(2)
                                                            catch
                                                            end
                                                            figure(2)
                                                            
                                                            hold on
                                                            
                                                            max_dB=max([max(this_mouse_peakPACpowerEv1) max(this_mouse_peakPACpowerEv2)]);
                                                            min_dB=min([min(this_mouse_peakPACpowerEv1) min(this_mouse_peakPACpowerEv2)]);
                                                            
                                                            edges=[min_dB-0.1*(max_dB-min_dB):(max_dB-min_dB)/20:max_dB+0.1*(max_dB-min_dB)];
                                                            histogram(this_mouse_peakPACpowerEv1,edges,'FaceColor','b','EdgeColor','b')
                                                            histogram(this_mouse_peakPACpowerEv2,edges,'FaceColor','r','EdgeColor','r')
                                                            xlabel('delta power dB')
                                                            title(['Histogram for concentrations ' num2str(concs2(evNo1)) ' and ' num2str(concs2(evNo2))])
                                                            pffft=1;
                                                        end
                                                        
                                                    end
                                                end
                                                %
                                            end
                                        end
                                    end
                                    %                                         end
                                end
                                
                            end %Can calculate auROC
                        end
                        
                    end
                    
                end
                
                if mouse_has_files==1
                    mouse_included(mouseNo)=1;
                else
                    mouse_included(mouseNo)=0;
                end
            end
            
            
        end
        fprintf(1, '\n\n')
        
        fprintf(1, ['The number of mice included in the PAC analysis for this odor pair is %d\n\n\n'], sum(mouse_included))
        
        figureNo = 0;
        %Now plot the average peakPACpower for each electrode calculated per mouse
        %(including all sessions for each mouse)
        edges=[-25:0.5:15];
        
        rand_offset=0.8;
        
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            data_PACpower=[];
            prof_naive=[];
            events=[];
            groups=[];
            mice=[];
            
            %Plot the average
            figureNo = figureNo + 1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            set(hFig, 'units','normalized','position',[.1 .2 .7 .7])
            %             subplot(2,1,1)
            ax=gca;ax.LineWidth=3;
            
            hold on
            
            bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            bar_offset = 0;
            
            %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
            
            include_group=0;
            
            for evNo=1:length(eventType)
                
                for per_ii=2:-1:1
                    
                    for grNo=1:max(handles_drgb.drgbchoices.group_no)
                        bar_offset = bar_offset +1;
                        %                         if sum(eventType==3)>0
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
                        %                         else
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        %                         end
                        
                        %                         these_offsets(per_ii)=bar_offset;
                        bar_offset = bar_offset + 1;
                        
                        if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
                            
                            include_group=1;
                            
                            switch grNo
                                case 1
                                    bar(bar_offset,mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'g','LineWidth', 3,'EdgeColor','none')
                                case 2
                                    bar(bar_offset,mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                                case 3
                                    bar(bar_offset,mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'y','LineWidth', 3,'EdgeColor','none')
                            end
                            
                            %Save data
                            handles_out.PRP_ii=handles_out.PRP_ii+1;
                            handles_out.PRP_values(handles_out.PRP_ii).pacii=pacii;
                            handles_out.PRP_values(handles_out.PRP_ii).evNo=evNo;
                            handles_out.PRP_values(handles_out.PRP_ii).per_ii=per_ii;
                            handles_out.PRP_values(handles_out.PRP_ii).groupNo=grNo;
                            handles_out.PRP_values(handles_out.PRP_ii).peak=1;
                            handles_out.PRP_values(handles_out.PRP_ii).PRP=mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)));
                            
                            %Save the mean per mouse
                            these_mice=mean_PACpower_mouseNo_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo));
                            handles_out.PRP_values(handles_out.PRP_ii).noMice=0;
                            for iiMice=min(these_mice):max(these_mice)
                                if sum(these_mice==iiMice)>0
                                    handles_out.PRP_values(handles_out.PRP_ii).noMice=handles_out.PRP_values(handles_out.PRP_ii).noMice+1;
                                    handles_out.PRP_values(handles_out.PRP_ii).mouseNo(handles_out.PRP_values(handles_out.PRP_ii).noMice)=iiMice;
                                    handles_out.PRP_values(handles_out.PRP_ii).PRP_per_mouse(handles_out.PRP_values(handles_out.PRP_ii).noMice)=mean(mean_peakPACpower_per_mouse((mean_PACpower_mouseNo_per_mouse==iiMice)&(~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)));
                                end
                            end
                            
                            %Violin plot
                            [mean_out, CIout]=drgViolinPoint(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))...
                                ,edges,bar_offset,rand_offset,'k','k',3);
                            
                            %                             plot(bar_offset,mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
                            %                             plot((bar_offset)*ones(1,sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),...
                            %                                 mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)),'o',...
                            %                                 'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            %
                            %                             if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>=2
                            %                                 CI = bootci(1000, {@mean, mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))},'type','cper');
                            %                                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            %                             end
                            
                            % %Save data for anovan
                            % data_PACpower=[data_PACpower mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
                            % prof_naive=[prof_naive per_ii*ones(1,sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            % events=[events evNo*ones(1,sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            % groups=[groups grNo*ones(1,sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            % mice=[mice mean_PACpower_mouseNo_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
                            
                        end
                    end
                    bar_offset = bar_offset + 2;
                    %                     if include_group==1
                    %                         bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                    %                         no_ev_labels=no_ev_labels+1;
                    %                         if sum(eventType==3)>0
                    %                             bar_labels{no_ev_labels}=evTypeLabels{evNo};
                    %                         else
                    %                             bar_labels{no_ev_labels}=num2str(concs(evNo));
                    %                         end
                    %                     end
                end
                bar_offset = bar_offset + 3;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%aqui
                %                 if include_group==1
                %                     ii_gr_included=ii_gr_included+1;
                %                     groups_included(ii_gr_included)=grNo;
                %                 end
                
            end
            
            title(['Peak PAC power theta/' freq_names{pacii+1}])
            
            
            %Annotations identifying groups
            x_interval=0.8/ii_gr_included;
            for ii=1:ii_gr_included
                annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %             %x labels
            %             to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
            %             sorted_A=sortrows(to_sort);
            %             sorted_bar_lab_loc=sorted_A(:,1);
            %             for ii=1:length(bar_lab_loc)
            %                 sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
            %             end
            xticks([2 4 6 10 12 14 21 23 25 29 31 33])
            xticklabels({'nwtS+', 'nHETS+', 'nKOS+', 'pwtS+', 'pHETS+', 'pKOS+', 'nwtS-', 'nHETS-', 'nKOS-', 'pwtS-', 'pHETS-', 'pKOS-'})
            
            if sum(eventType==3)==0
                xlabel('Concentration (%)')
            end
            
            ylabel('dB')
            
            %ylim1(pacii,:)=ylim;
            ylim([-7 5])
            
            
        end
        
        
        
        
        
        
        
        %Now plot the average through PACpower for each electrode calculated per mouse
        %(including all sessions for each mouse)
        edges=[-25:0.5:15];
        
        rand_offset=0.8;
        
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            data_PACpower=[];
            prof_naive=[];
            events=[];
            groups=[];
            mice=[];
            
            %Plot the average
            figureNo = figureNo + 1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            set(hFig, 'units','normalized','position',[.1 .2 .7 .7])
            %             subplot(2,1,1)
            
            hold on
            
            bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            bar_offset = 0;
            
            %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
            
            include_group=0;
            
            for evNo=1:length(eventType)
                
                for per_ii=2:-1:1
                    
                    for grNo=1:max(handles_drgb.drgbchoices.group_no)
                        bar_offset = bar_offset +1;
                        %                         if sum(eventType==3)>0
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
                        %                         else
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        %                         end
                        
                        %                         these_offsets(per_ii)=bar_offset;
                        bar_offset = bar_offset + 1;
                        
                        if sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
                            
                            include_group=1;
                            
                            switch grNo
                                case 1
                                    bar(bar_offset,mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'g','LineWidth', 3,'EdgeColor','none')
                                case 2
                                    bar(bar_offset,mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                                case 3
                                    bar(bar_offset,mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'y','LineWidth', 3,'EdgeColor','none')
                            end
                            
                            %Save data
                            
                            %Save data
                            handles_out.PRP_ii=handles_out.PRP_ii+1;
                            handles_out.PRP_values(handles_out.PRP_ii).pacii=pacii;
                            handles_out.PRP_values(handles_out.PRP_ii).evNo=evNo;
                            handles_out.PRP_values(handles_out.PRP_ii).per_ii=per_ii;
                            handles_out.PRP_values(handles_out.PRP_ii).groupNo=grNo;
                            handles_out.PRP_values(handles_out.PRP_ii).peak=0;
                            handles_out.PRP_values(handles_out.PRP_ii).PRP=mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)));
                            
                            %Save the mean per mouse
                            these_mice=mean_PACpower_mouseNo_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo));
                            handles_out.PRP_values(handles_out.PRP_ii).noMice=0;
                            for iiMice=min(these_mice):max(these_mice)
                                if sum(these_mice==iiMice)>0
                                    handles_out.PRP_values(handles_out.PRP_ii).noMice=handles_out.PRP_values(handles_out.PRP_ii).noMice+1;
                                    handles_out.PRP_values(handles_out.PRP_ii).mouseNo(handles_out.PRP_values(handles_out.PRP_ii).noMice)=iiMice;
                                    handles_out.PRP_values(handles_out.PRP_ii).PRP_per_mouse(handles_out.PRP_values(handles_out.PRP_ii).noMice)=mean(mean_troughPACpower_per_mouse((mean_PACpower_mouseNo_per_mouse==iiMice)&(~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)));
                                end
                            end
                            
                            %Violin plot
                            [mean_out, CIout]=drgViolinPoint(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))...
                                ,edges,bar_offset,rand_offset,'k','k',3);
                            
                            
                            %                             plot(bar_offset,mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
                            %                             plot((bar_offset)*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),...
                            %                                 mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)),'o',...
                            %                                 'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            %
                            %                             if sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>=2
                            %                                 CI = bootci(1000, {@mean, mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))},'type','cper');
                            %                                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            %                             end
                            
                            %                             %Save data for anovan
                            %                             data_PACpower=[data_PACpower mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
                            %                             prof_naive=[prof_naive per_ii*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            %                             events=[events evNo*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            %                             groups=[groups grNo*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            %                             mice=[mice mean_PACpower_mouseNo_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
                            %
                        end
                    end
                    bar_offset = bar_offset + 2;
                    %                     if include_group==1
                    %                         bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                    %                         no_ev_labels=no_ev_labels+1;
                    %                         if sum(eventType==3)>0
                    %                             bar_labels{no_ev_labels}=evTypeLabels{evNo};
                    %                         else
                    %                             bar_labels{no_ev_labels}=num2str(concs(evNo));
                    %                         end
                    %                     end
                end
                bar_offset = bar_offset + 3;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%aqui
                %                 if include_group==1
                %                     ii_gr_included=ii_gr_included+1;
                %                     groups_included(ii_gr_included)=grNo;
                %                 end
                
            end
            
            title(['Trough PAC power theta/' freq_names{pacii+1}])
            
            
            %Annotations identifying groups
            x_interval=0.8/ii_gr_included;
            for ii=1:ii_gr_included
                annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %x labels
            to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
            sorted_A=sortrows(to_sort);
            sorted_bar_lab_loc=sorted_A(:,1);
            for ii=1:length(bar_lab_loc)
                sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
            end
            xticks(sorted_bar_lab_loc)
            %             xticklabels(sorted_bar_labels)
            
            if sum(eventType==3)==0
                xlabel('Concentration (%)')
            end
            
            ylabel('dB')
            %             ylim1(pacii,:)=ylim;
            ylim([-7 5])
            
        end
        
        
        
        %Now plot the average lick PACpower for each electrode calculated per mouse
        %(including all sessions for each mouse)
        edges=[-25:0.5:15];
        
        rand_offset=0.8;
        
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            data_PACpower=[];
            prof_naive=[];
            events=[];
            groups=[];
            mice=[];
            
            %Plot the average
            figureNo = figureNo + 1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            set(hFig, 'units','normalized','position',[.1 .2 .7 .7])
            %             subplot(2,1,1)
            
            hold on
            
            bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            bar_offset = 0;
            
            %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
            
            include_group=0;
            
            for evNo=1:length(eventType)
                
                for per_ii=2:-1:1
                    
                    for grNo=1:max(handles_drgb.drgbchoices.group_no)
                        bar_offset = bar_offset +1;
                        %                         if sum(eventType==3)>0
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
                        %                         else
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        %                         end
                        
                        %                         these_offsets(per_ii)=bar_offset;
                        bar_offset = bar_offset + 1;
                        
                        if sum((~isnan(mean_lickPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
                            
                            include_group=1;
                            
                            switch grNo
                                case 1
                                    bar(bar_offset,mean(mean_lickPACpower_per_mouse((~isnan(mean_lickPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'g','LineWidth', 3,'EdgeColor','none')
                                case 2
                                    bar(bar_offset,mean(mean_lickPACpower_per_mouse((~isnan(mean_lickPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                                case 3
                                    bar(bar_offset,mean(mean_lickPACpower_per_mouse((~isnan(mean_lickPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'y','LineWidth', 3,'EdgeColor','none')
                            end
                            
                            %Save data
                            
                            %Save data
                            handles_out.PRP_ii=handles_out.PRP_ii+1;
                            handles_out.PRP_values(handles_out.PRP_ii).pacii=pacii;
                            handles_out.PRP_values(handles_out.PRP_ii).evNo=evNo;
                            handles_out.PRP_values(handles_out.PRP_ii).per_ii=per_ii;
                            handles_out.PRP_values(handles_out.PRP_ii).groupNo=grNo;
                            handles_out.PRP_values(handles_out.PRP_ii).peak=0;
                            handles_out.PRP_values(handles_out.PRP_ii).PRP=mean(mean_lickPACpower_per_mouse((~isnan(mean_lickPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)));
                            
                            %Save the mean per mouse
                            these_mice=mean_PACpower_mouseNo_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo));
                            handles_out.PRP_values(handles_out.PRP_ii).noMice=0;
                            for iiMice=min(these_mice):max(these_mice)
                                if sum(these_mice==iiMice)>0
                                    handles_out.PRP_values(handles_out.PRP_ii).noMice=handles_out.PRP_values(handles_out.PRP_ii).noMice+1;
                                    handles_out.PRP_values(handles_out.PRP_ii).mouseNo(handles_out.PRP_values(handles_out.PRP_ii).noMice)=iiMice;
                                    handles_out.PRP_values(handles_out.PRP_ii).PRP_per_mouse(handles_out.PRP_values(handles_out.PRP_ii).noMice)=mean(mean_lickPACpower_per_mouse((mean_PACpower_mouseNo_per_mouse==iiMice)&(~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)));
                                end
                            end
                            
                            %Violin plot
                            [mean_out, CIout]=drgViolinPoint(mean_lickPACpower_per_mouse((~isnan(mean_lickPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))...
                                ,edges,bar_offset,rand_offset,'k','k',3);
                            
                            
                            %                             plot(bar_offset,mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
                            %                             plot((bar_offset)*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),...
                            %                                 mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)),'o',...
                            %                                 'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            %
                            %                             if sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>=2
                            %                                 CI = bootci(1000, {@mean, mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))},'type','cper');
                            %                                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            %                             end
                            
                            %                             %Save data for anovan
                            %                             data_PACpower=[data_PACpower mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
                            %                             prof_naive=[prof_naive per_ii*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            %                             events=[events evNo*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            %                             groups=[groups grNo*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            %                             mice=[mice mean_PACpower_mouseNo_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
                            %
                        end
                    end
                    bar_offset = bar_offset + 2;
                    %                     if include_group==1
                    %                         bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                    %                         no_ev_labels=no_ev_labels+1;
                    %                         if sum(eventType==3)>0
                    %                             bar_labels{no_ev_labels}=evTypeLabels{evNo};
                    %                         else
                    %                             bar_labels{no_ev_labels}=num2str(concs(evNo));
                    %                         end
                    %                     end
                end
                bar_offset = bar_offset + 3;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%aqui
                %                 if include_group==1
                %                     ii_gr_included=ii_gr_included+1;
                %                     groups_included(ii_gr_included)=grNo;
                %                 end
                
            end
            
            title(['Lick PAC power theta/' freq_names{pacii+1}])
            
            
            %Annotations identifying groups
            x_interval=0.8/ii_gr_included;
            for ii=1:ii_gr_included
                annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %x labels
            to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
            sorted_A=sortrows(to_sort);
            sorted_bar_lab_loc=sorted_A(:,1);
            for ii=1:length(bar_lab_loc)
                sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
            end
            xticks(sorted_bar_lab_loc)
            %             xticklabels(sorted_bar_labels)
            
            if sum(eventType==3)==0
                xlabel('Concentration (%)')
            end
            
            ylabel('dB')
            %             ylim1(pacii,:)=ylim;
            ylim([-7 5])
            
        end
        
        
        %Display auROC
        edges=[-0.3:0.05:0.5];
        rand_offset=0.8;
        for pacii=1:no_pacii    %for different PACs
            
            ii_roc=0;
            roc_data=[];
            glm_roc=[];
            glm_roc_ii=0;
            
            ii_rocpk=0;
            rocpk_data=[];
            glm_rocpk=[];
            glm_rocpk_ii=0;
            
            ii_roctr=0;
            roctr_data=[];
            glm_roctr=[];
            glm_roctr_ii=0;
            
            ii_roclk=0;
            roclk_data=[];
            glm_roclk=[];
            glm_roclk_ii=0;
            
            %Display the peak auROC
            figureNo=figureNo+1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            
            
            bar_offset=0;
            
            for per_ii=2:-1:1
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    
                    if sum((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))>0
                        
                        bar_offset=bar_offset+1;
                        
                        switch grNo
                            case 1
                                bar(bar_offset,mean(auROCpeak((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))),'g','LineWidth', 3,'EdgeColor','none')
                            case 2
                                bar(bar_offset,mean(auROCpeak((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                            case 3
                                bar(bar_offset,mean(auROCpeak((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))),'y','LineWidth', 3,'EdgeColor','none')
                        end
                        
                        
                        %Save data
                        handles_out.AUC_ii=handles_out.AUC_ii+1;
                        handles_out.AUC_values(handles_out.AUC_ii).pacii=pacii;
                        handles_out.AUC_values(handles_out.AUC_ii).evNo=evNo;
                        handles_out.AUC_values(handles_out.AUC_ii).per_ii=per_ii;
                        handles_out.AUC_values(handles_out.AUC_ii).groupNo=grNo;
                        handles_out.AUC_values(handles_out.AUC_ii).peak=1;
                        handles_out.AUC_values(handles_out.AUC_ii).AUC=mean(auROCpeak((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo)));
                        
                        %Violin plot
                        [mean_out, CIout]=drgViolinPoint(auROCpeak((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))...
                            ,edges,bar_offset,rand_offset,'k','k',3);
                        
                        
                        %Enter data for glm for peak only
                        these_data=auROCpeak((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))';
                        glm_rocpk.data(glm_rocpk_ii+1:glm_rocpk_ii+length(these_data))=these_data;
                        glm_rocpk.group(glm_rocpk_ii+1:glm_rocpk_ii+length(these_data))=grNo;
                        glm_rocpk.perCorr(glm_rocpk_ii+1:glm_rocpk_ii+length(these_data))=per_ii;
                        glm_rocpk_ii=glm_rocpk_ii+length(these_data);
                        
                        %Enter the data for t-test/ranksum
                        ii_rocpk=ii_rocpk+1;
                        rocpk_data(ii_rocpk).data=these_data;
                        rocpk_data(ii_rocpk).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' prof_naive_leg{per_ii}];
                        
                        
                        %Enter data for glm for peak and trough only
                        glm_roc.data(glm_roc_ii+1:glm_roc_ii+length(these_data))=these_data;
                        glm_roc.group(glm_roc_ii+1:glm_roc_ii+length(these_data))=grNo;
                        glm_roc.perCorr(glm_roc_ii+1:glm_roc_ii+length(these_data))=per_ii;
                        glm_roc.peak_trough(glm_roc_ii+1:glm_roc_ii+length(these_data))=1;
                        glm_roc_ii=glm_roc_ii+length(these_data);
                        
                        %Enter the data for t-test/ranksum
                        ii_roc=ii_roc+1;
                        roc_data(ii_roc).data=these_data;
                        roc_data(ii_roc).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' prof_naive_leg{per_ii} ' peak'];
                        
                        
                        
                        
                    end
                    
                    
                end
                
                
                bar_offset=bar_offset+1;
            end
            
            title(['auROC per mouse, per electrode peak-referenced power theta/' freq_names{pacii+1}])
            ylim([-0.2 0.5])
            
            
            ylabel('auROC')
            
            pffft=1;
            
            
            
            %Display the trough auROC for all trials per mouse (per electrode) for within vs
            %betweeen
            
            %Plot the average
            figureNo=figureNo+1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            
            
            bar_offset=0;
            
            for per_ii=2:-1:1
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    
                    if sum((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))>0
                        
                        bar_offset=bar_offset+1;
                        
                        switch grNo
                            case 1
                                bar(bar_offset,mean(auROCtrough((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))),'g','LineWidth', 3,'EdgeColor','none')
                            case 2
                                bar(bar_offset,mean(auROCtrough((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                            case 3
                                bar(bar_offset,mean(auROCtrough((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))),'y','LineWidth', 3,'EdgeColor','none')
                        end
                        
                        
                        %Save data
                        handles_out.AUC_ii=handles_out.AUC_ii+1;
                        handles_out.AUC_values(handles_out.AUC_ii).pacii=pacii;
                        handles_out.AUC_values(handles_out.AUC_ii).evNo=evNo;
                        handles_out.AUC_values(handles_out.AUC_ii).per_ii=per_ii;
                        handles_out.AUC_values(handles_out.AUC_ii).groupNo=grNo;
                        handles_out.AUC_values(handles_out.AUC_ii).peak=0;
                        handles_out.AUC_values(handles_out.AUC_ii).AUC=mean(auROCtrough((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo)));
                        
                        
                        %Violin plot
                        
                        [mean_out, CIout]=drgViolinPoint(auROCtrough((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))...
                            ,edges,bar_offset,rand_offset,'k','k',3);
                        
                        
                        %Enter data for glm for peak only
                        these_data=auROCtrough((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))';
                        glm_roctr.data(glm_roctr_ii+1:glm_roctr_ii+length(these_data))=these_data;
                        glm_roctr.group(glm_roctr_ii+1:glm_roctr_ii+length(these_data))=grNo;
                        glm_roctr.perCorr(glm_roctr_ii+1:glm_roctr_ii+length(these_data))=per_ii;
                        glm_roctr_ii=glm_roctr_ii+length(these_data);
                        
                        %Enter the data for t-test/ranksum
                        ii_roctr=ii_roctr+1;
                        roctr_data(ii_roctr).data=these_data;
                        roctr_data(ii_roctr).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' prof_naive_leg{per_ii}];
                        
                        
                        %Enter data for glm for peak and trough only
                        glm_roc.data(glm_roc_ii+1:glm_roc_ii+length(these_data))=these_data;
                        glm_roc.group(glm_roc_ii+1:glm_roc_ii+length(these_data))=grNo;
                        glm_roc.perCorr(glm_roc_ii+1:glm_roc_ii+length(these_data))=per_ii;
                        glm_roc.peak_trough(glm_roc_ii+1:glm_roc_ii+length(these_data))=2;
                        glm_roc_ii=glm_roc_ii+length(these_data);
                        
                        %Enter the data for t-test/ranksum
                        ii_roc=ii_roc+1;
                        roc_data(ii_roc).data=these_data;
                        roc_data(ii_roc).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' prof_naive_leg{per_ii} ' peak'];
                        
                        
                        
                    end
                    
                    
                end
                
                
                bar_offset=bar_offset+1;
            end
            
            title(['auROC per mouse, per electrode trough-referenced power theta/' freq_names{pacii+1}])
            ylim([-0.2 0.5])
            
            
            ylabel('auROC')
            
            %Display the lick-refrenced power auROC for all trials per mouse (per electrode) for within vs
            %betweeen
            
            %Plot the average
            figureNo=figureNo+1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            
            
            bar_offset=0;
            
            for per_ii=2:-1:1
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    
                    if sum((ROClick_per_ii==per_ii)&(ROClick_pacii==pacii)&(ROClick_groups==grNo))>0
                        
                        bar_offset=bar_offset+1;
                        
                        switch grNo
                            case 1
                                bar(bar_offset,mean(auROClick((ROClick_per_ii==per_ii)&(ROClick_pacii==pacii)&(ROClick_groups==grNo))),'g','LineWidth', 3,'EdgeColor','none')
                            case 2
                                bar(bar_offset,mean(auROClick((ROClick_per_ii==per_ii)&(ROClick_pacii==pacii)&(ROClick_groups==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                            case 3
                                bar(bar_offset,mean(auROClick((ROClick_per_ii==per_ii)&(ROClick_pacii==pacii)&(ROClick_groups==grNo))),'y','LineWidth', 3,'EdgeColor','none')
                        end
                        
                        
                        %Save data
                        handles_out.AUClick_ii=handles_out.AUClick_ii+1;
                        handles_out.AUClick_values(handles_out.AUClick_ii).pacii=pacii;
                        handles_out.AUClick_values(handles_out.AUClick_ii).evNo=evNo;
                        handles_out.AUClick_values(handles_out.AUClick_ii).per_ii=per_ii;
                        handles_out.AUClick_values(handles_out.AUClick_ii).groupNo=grNo;
                        handles_out.AUClick_values(handles_out.AUClick_ii).peak=0;
                        handles_out.AUClick_values(handles_out.AUClick_ii).AUC=mean(auROClick((ROClick_per_ii==per_ii)&(ROClick_pacii==pacii)&(ROClick_groups==grNo)));
                        
                        
                        %Violin plot
                        
                        [mean_out, CIout]=drgViolinPoint(auROClick((ROClick_per_ii==per_ii)&(ROClick_pacii==pacii)&(ROClick_groups==grNo))...
                            ,edges,bar_offset,rand_offset,'k','k',3);
                        
                        
                        %Enter data for glm for peak only
                        these_data=auROClick((ROClick_per_ii==per_ii)&(ROClick_pacii==pacii)&(ROClick_groups==grNo))';
                        glm_roclk.data(glm_roclk_ii+1:glm_roclk_ii+length(these_data))=these_data;
                        glm_roclk.group(glm_roclk_ii+1:glm_roclk_ii+length(these_data))=grNo;
                        glm_roclk.perCorr(glm_roclk_ii+1:glm_roclk_ii+length(these_data))=per_ii;
                        glm_roclk_ii=glm_roclk_ii+length(these_data);
                        
                        %Enter the data for t-test/ranksum
                        ii_roclk=ii_roclk+1;
                        roclk_data(ii_roclk).data=these_data;
                        roclk_data(ii_roclk).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' prof_naive_leg{per_ii}];
                        
                        
                        %Enter data for glm for peak and trough only
                        glm_roc.data(glm_roc_ii+1:glm_roc_ii+length(these_data))=these_data;
                        glm_roc.group(glm_roc_ii+1:glm_roc_ii+length(these_data))=grNo;
                        glm_roc.perCorr(glm_roc_ii+1:glm_roc_ii+length(these_data))=per_ii;
                        glm_roc.peak_trough(glm_roc_ii+1:glm_roc_ii+length(these_data))=3;
                        glm_roc_ii=glm_roc_ii+length(these_data);
                        
                        %Enter the data for t-test/ranksum
                        ii_roc=ii_roc+1;
                        roc_data(ii_roc).data=these_data;
                        roc_data(ii_roc).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' prof_naive_leg{per_ii} ' lick'];
                        
                        
                        
                    end
                    
                    
                end
                
                
                bar_offset=bar_offset+1;
            end
            
            title(['auROC per mouse, per electrode lick-referenced power theta/' freq_names{pacii+1}])
            ylim([-0.2 0.5])
            
            
            ylabel('auROC')
            
            %Do glm for peak/trough
            fprintf(1, ['\n\nglm for auROC peak/trough for Theta/' freq_names{pacii+1} '\n'])
            tbl = table(glm_roc.data',glm_roc.group',glm_roc.perCorr',glm_roc.peak_trough',...
                'VariableNames',{'Peak_wave','group','perCorr','peak_trough'});
            mdl = fitglm(tbl,'Peak_wave~group+perCorr+peak_trough+group*perCorr*peak_trough'...
                ,'CategoricalVars',[2,3,4])
            
            
            fprintf(1, ['\n\nRanksum or t-test for auROC peak for theta ' freq_names{pacii+1} '\n'])
            %Now do the ranksums
            output_data = drgMutiRanksumorTtest(rocpk_data);
            
            fprintf(1, ['\n\nRanksum or t-test for auROC trough for theta ' freq_names{pacii+1} '\n'])
            %Now do the ranksums
            output_data = drgMutiRanksumorTtest(roctr_data);
            
            fprintf(1, ['\n\nRanksum or t-test for auROC trough for lick ' freq_names{pacii+1} '\n'])
            %Now do the ranksums
            output_data = drgMutiRanksumorTtest(roclk_data);
            
            pffft=1;
            
        end
        
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) handles_pars.output_suffix],'handles_out')
        pfft=1;
        
    case 2
        %1 Oscillatory wavelet power timecourse divergence after odorant addition
        
        t_odor_on=0.1045;  %From Supplementary Fig 1 of Losacco et al 2020 DOI: 10.7554/eLife.52583
        use_rs=0; %If this is 1 then ranksums are used for the decision time p value
                    %ranksum does not work well because it often goes below 0.05 before the
                    %odor
    
        handles_out=[];
        handles_out.PRP_ii=0;
        handles_out.LickF_ii=0;
        handles_out.AUC_ii=0;
        handles_out.AUClick_ii=0;
       
        
        group_legend{1}='WT';
        group_legend{2}='Het';
        group_legend{3}='KO';
        
        prof_naive_leg{1}='Proficient';
        prof_naive_leg{2}='Naive';
        
      
        
        PRPtimecourse=[];
        
        for ii_mouse=1:max(handles_drgb.drgbchoices.mouse_no)
            PRPtimecourse.mouse(ii_mouse).no_sessions=0;
            PRPtimecourse.mouse(ii_mouse).files=[];
        end
        
        fprintf(1, ['PAC power analysis using wavelet power for Justin''s paper\n\n'])
        p_vals=[];
        no_files=max(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        
        szpc=size(percent_windows);
        for per_ii=1:szpc(1)
            
            
            
            for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
                mouse_has_files=0;
                for elec=1:16
                    if sum(which_electrodes==elec)>0
                  
                        for fileNo=1:no_files
                            %If this file is in the list of files the user wants to process in drgAnalysisBatchLFP continue
                            if sum(files==fileNo)>0
                                if handles_drgb.drgbchoices.mouse_no(fileNo)==mouseNo
                                    
                                    if PRPtimecourse.mouse(mouseNo).no_sessions==0
                                        PRPtimecourse.mouse(mouseNo).no_sessions=1;
                                        PRPtimecourse.mouse(mouseNo).files(1)=fileNo;
                                        this_session=1;
                                    else
                                        this_session=find(PRPtimecourse.mouse(mouseNo).files==fileNo);
                                        if isempty(this_session)
                                            PRPtimecourse.mouse(mouseNo).no_sessions=PRPtimecourse.mouse(mouseNo).no_sessions+1;
                                            PRPtimecourse.mouse(mouseNo).files(PRPtimecourse.mouse(mouseNo).no_sessions)=fileNo;
                                            this_session=PRPtimecourse.mouse(mouseNo).no_sessions;
                                        end
                                    end
                                    
                                    PRPtimecourse.mouse(mouseNo).group_no=handles_drgb.drgbchoices.group_no(fileNo);
                                    
                                    lfpodNo=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                    %                                     lfpodRefNo=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                                    if (mouseNo==2)&(this_session==4)
                                        pffft=1;
                                    end
                                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo)))
                                        
                                        
                                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo).PAC))
                                            
                                            percent_mask=[];
                                            percent_mask=logical((handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).perCorrPAC>=percent_windows(per_ii,1))...
                                                &(handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).perCorrPAC<=percent_windows(per_ii,2)));
                                            
                                            if ~isempty(percent_mask)
                                                
                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed=1;
                                                
                                                for evNo=1:length(eventType)
                                                    
                                                    trials_in_event_Ev=[];
                                                    trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).which_eventPAC(eventType(evNo),:)==1)&percent_mask;
                                                    for pacii=1:no_pacii
                                                        PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).no_events=sum(trials_in_event_Ev);
                                                    end
                                                    
                                                    if (sum(trials_in_event_Ev)>=1)
                                                        
                                                           trialNos_in_event_Ev=zeros(1,sum(trials_in_event_Ev));
                                                        ii_no=0;
                                                        for ii=1:length(trials_in_event_Ev)
                                                            if trials_in_event_Ev(ii)==1
                                                                ii_no=ii_no+1;
                                                                trialNos_in_event_Ev(ii_no)=ii;
                                                            end
                                                        end

                                                        %Do per bandwidth analysis
                                                        for pacii=1:no_pacii
                                                            
                                                            group_no=handles_drgb.drgbchoices.group_no(fileNo);
                                                              
                                                            %Enter the PACpower
                                                            t_pac=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.t_pac;
                                                            PRPtimecourse.t_pac=t_pac((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
                                                            
                                                            %Get peak power
                                                            these_ref_power=[];
                                                            for ww=1:length(trialNos_in_event_Ev)
                                                                this_p=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_event_Ev(ww)).peakPower;
                                                                these_ref_power=[these_ref_power this_p((t_pac>=handles_pars.refWin_times(1))&(t_pac<=handles_pars.refWin_times(2)))];
                                                            end
                                                            ref_power_SD=std(these_ref_power);
                                                            ref_power_mean=mean(these_ref_power);
                                                            
                                                           
                                                            for ww=1:length(trialNos_in_event_Ev)
                                                                this_p=(handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_event_Ev(ww)).peakPower-ref_power_mean)/ref_power_SD;
                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peak=this_p((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
                                                            end
                                                            PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).no_z_power=sum(trials_in_event_Ev);
                                                            
                                                            %Save the times for peaks
                                                            for ww=1:length(trialNos_in_event_Ev)
                                                                these_peak_times=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_event_Ev(ww)).peakPower_times;
                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peak_times=these_peak_times;
                                                            end
                                                            
                                                            %Get trough power
                                                            these_ref_power=[];
                                                            for ww=1:length(trialNos_in_event_Ev)
                                                                this_p=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_event_Ev(ww)).troughPower;
                                                                these_ref_power=[these_ref_power this_p((t_pac>=handles_pars.refWin_times(1))&(t_pac<=handles_pars.refWin_times(2)))];
                                                            end
                                                            ref_power_SD=std(these_ref_power);
                                                            ref_power_mean=mean(these_ref_power);
                                                            
                                                            for ww=1:length(trialNos_in_event_Ev)
                                                                this_p=(handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_event_Ev(ww)).troughPower-ref_power_mean)/ref_power_SD;
                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).trough=this_p((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
                                                            end
                                                            
                                                            %Save the times for troughs
                                                            for ww=1:length(trialNos_in_event_Ev)
                                                                these_trough_times=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_event_Ev(ww)).troughPower_times;
                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).trough_times=these_trough_times;
                                                            end
                                                            
                                                            
                                                            %Now do the lick power and other lick computations
                                                            if isfield(handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(1).PACwave,'lickPower_trials')
                                                                if length(handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(1).PACwave.lickPower_trials)>0
                                                                    
                                                                    lick_trials=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.lickPower_trials;
                                                                    
                                                                    %Do lick referenced power
                                                                    lickPower_trials=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(1).PACwave.lickPower_trials;
                                                                    trials_in_lick_event=zeros(1,length(lickPower_trials(lickPower_trials<=length(trials_in_event_Ev))));
                                                                    trials_in_lick_event(1,:)=trials_in_event_Ev(lickPower_trials(lickPower_trials<=length(trials_in_event_Ev)));
                                                                    
                                                                    trialNos_in_lick_event_Ev=zeros(1,sum(trials_in_lick_event));
                                                                    ii_no=0;
                                                                    for ii=1:length(trials_in_lick_event)
                                                                        if trials_in_lick_event(ii)==1
                                                                            ii_no=ii_no+1;
                                                                            trialNos_in_lick_event_Ev(ii_no)=ii;
                                                                        end
                                                                    end
                                                                    
                                                                    these_ref_power=[];
                                                                    for ww=1:length(trialNos_in_lick_event_Ev)
                                                                        this_p=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_lick_event_Ev(ww)).lickPower;
                                                                        these_ref_power=[these_ref_power this_p((t_pac>=handles_pars.refWin_times(1))&(t_pac<=handles_pars.refWin_times(2)))];
                                                                    end
                                                                    ref_power_SD=std(these_ref_power);
                                                                    ref_power_mean=mean(these_ref_power);
                                                                    
                                                                    
                                                                    for ww=1:length(trialNos_in_lick_event_Ev)
                                                                        this_p=(handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_lick_event_Ev(ww)).lickPower-ref_power_mean)/ref_power_SD;
                                                                        PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).lickp=this_p((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
                                                                    end
                                                                    PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).no_z_lick_power=sum(trials_in_lick_event);
                                                                    PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).lick_trials=lick_trials;
                                                                    
                                                                    %Do lick accounting for lick p value
                                                                    ii_t_first=find((t_pac>=handles_pars.analysisWin_times(1)),1,'first');
                                                                    ii_t_last=find(t_pac<=handles_pars.analysisWin_times(2),1,'last');
                                                                    for ww=1:length(trialNos_in_lick_event_Ev)
                                                                        licks_this_trial=zeros(1,ii_t_last-ii_t_first+1);
                                                                        these_lick_times=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_lick_event_Ev(ww)).lickPower_times;
                                                                        for ii_lick=1:length(these_lick_times)
                                                                            if (these_lick_times(ii_lick)>=t_pac(ii_t_first))&(these_lick_times(ii_lick)<t_pac(ii_t_last))
                                                                                this_ii=find((t_pac<=these_lick_times(ii_lick))&(t_pac+(t_pac(2)-t_pac(1))>these_lick_times(ii_lick)),1,'first')-ii_t_first+1;
                                                                                licks_this_trial(this_ii)=1;
                                                                            end
                                                                        end
                                                                        PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).binary_lick=licks_this_trial;
                                                                    end
                                                                    
                                                                    %Save the times for licks
                                                                    for ww=1:length(trialNos_in_lick_event_Ev)
                                                                        these_lick_times=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_lick_event_Ev(ww)).lickPower_times;
                                                                        PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).lick_times=these_lick_times;
                                                                    end
                                                                else
                                                                    PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).no_z_lick_power=0;
                                                                end
                                                            else
                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).no_z_lick_power=0;
                                                            end
                                                            
                                                            
                                                            %Save lick analysis
                                                            if (pacii==1)&(which_electrodes(1)==elec)
                                                                %Do lick referenced power
                                                                lick_trials_included=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).lickStuff.lick_trials_included;
                                                                these_lick_trials=zeros(1,length(lick_trials_included));
                                                                for trNo=1:length(trials_in_event_Ev)
                                                                    if trials_in_event_Ev(trNo)==1
                                                                        this_lick_tr=find(lick_trials_included==trNo);
                                                                        if ~isempty(this_lick_tr)
                                                                            these_lick_trials(1,this_lick_tr)=1;
                                                                        end
                                                                    end
                                                                end
                                                                
                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).no_lick_events=sum(these_lick_trials);
                                                                last_ii=0;
                                                                for ww=1:sum(these_lick_trials)
                                                                    this_delta_ii=find(these_lick_trials(last_ii+1:end)==1,1,'first');
                                                                    PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).lick_tr(ww).binary_lick_per_t=...
                                                                        handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).lickStuff.binary_lick_per_t(last_ii+this_delta_ii,(t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
                                                                    last_ii=last_ii+this_delta_ii;
                                                                end
                                                            end
                                                        end
                                                        
                                                        
                                                        fprintf(1, ['%d trials in event No %d succesfully processed for file No %d electrode %d\n'],sum(trials_in_event_Ev), evNo,fileNo,elec);
                                                        
                                                    else
                                                        
                                                        
                                                        fprintf(1, ['%d trials in event No %d fewer than minimum trials per event ' evTypeLabels{evNo} ' for file No %d electrode %d\n'],sum(trials_in_event_Ev), min_trials_per_event,fileNo,elec);
                                                        
                                                        
                                                    end
                                                    
                                                    
                                                end
                                                
                                            else
                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed=0;
                                                fprintf(1, ['Empty percent_mask for file No %d electrode %d\n'],fileNo,elec);
                                                
                                            end
                                            
                                        else
                                            
                                            fprintf(1, ['Empty PAC for file No %d electrode %d\n'],fileNo,elec);
                                            
                                        end
                                        
                                        
                                    else
                                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],fileNo,elec);
                                        
                                        
                                    end
                                end %if mouseNo
                            end
                        end %fileNo

                    end
                    
                end
                
                if mouse_has_files==1
                    mouse_included(mouseNo)=1;
                else
                    mouse_included(mouseNo)=0;
                end
            end
            
            
        end
        fprintf(1, '\n\n')
        
        elec=[];
        
        figureNo = 0;
        anal_t_pac=t_pac((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
        %Let's try a per mouse z score t test (using the z scores for all
        %electrodes?)
        handles_out.RPtimecourse=[];
        handles_out.RP_ii=0;
        handles_out.anal_t_pac=anal_t_pac;
        
        for pacii=1:no_pacii
            for group_no=1:3
                handles_out.disc_t.pacii(pacii).groupNo(group_no).no_peak=0;
                handles_out.disc_t.pacii(pacii).groupNo(group_no).no_licks=0;
                handles_out.disc_t.pacii(pacii).groupNo(group_no).no_trough=0;
            end
        end
        
        figureNo = figureNo + 1;
        for pacii=1:no_pacii
            for per_ii=2:-1:1
                
                for mouseNo=1:length(PRPtimecourse.mouse)
                    handles_out.RP_ii=handles_out.RP_ii+1;
                    handles_out.RPtimecourse(handles_out.RP_ii).pacii=pacii;
                    handles_out.RPtimecourse(handles_out.RP_ii).per_ii=per_ii;
                    handles_out.RPtimecourse(handles_out.RP_ii).group_no=PRPtimecourse.mouse(mouseNo).group_no;
                    handles_out.RPtimecourse(handles_out.RP_ii).mouseNo=mouseNo;
                    
                    data_calculated=1;
                    
                    %Peak p values
                    SpPRPtimecourse_peak=[];
                    ii_Sp=0;
                    SmPRPtimecourse_peak=[];
                    ii_Sm=0;
                    for elecNo=which_electrodes
                        for this_session=1:PRPtimecourse.mouse(mouseNo).no_sessions
                            if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed==1
                                
                                %Splus
                                evNo=1;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).peak)
                                            this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                            this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).peak;
                                            ii_Sp=ii_Sp+1;
                                            SpPRPtimecourse_peak(ii_Sp,:)=this_PRPtimecourse;
                                        end
                                    end
                                end
                                
                                %Sminus
                                evNo=2;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).peak)
                                            this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                            this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).peak;
                                            ii_Sm=ii_Sm+1;
                                            SmPRPtimecourse_peak(ii_Sm,:)=this_PRPtimecourse;
                                        end
                                    end
                                end
                            end
                        end
                    end
                    handles_out.RPtimecourse(handles_out.RP_ii).p_vals_peak=[];
                    handles_out.RPtimecourse(handles_out.RP_ii).t_sig_peak=[];
                    
                    if (size(SpPRPtimecourse_peak,1)>10)&(size(SmPRPtimecourse_peak,1)>10)
                        p_vals_peak=zeros(1,length(anal_t_pac));
                        for ii_t=1:length(anal_t_pac)
                            if use_rs==1
                                p_vals_peak(ii_t)=ranksum(SpPRPtimecourse_peak(:,ii_t),SmPRPtimecourse_peak(:,ii_t));
                            else
                                [h,p_vals_peak(ii_t)]=ttest2(SpPRPtimecourse_peak(:,ii_t),SmPRPtimecourse_peak(:,ii_t));
                            end
                        end
                        
                        handles_out.RPtimecourse(handles_out.RP_ii).p_vals_peak=p_vals_peak;
                        ii_zero=find(anal_t_pac>=t_odor_on,1,'first');
                        delta_ii_sig=find(p_vals_peak(ii_zero:end)<=0.05,1,'first');
                        handles_out.RPtimecourse(handles_out.RP_ii).t_sig_peak=anal_t_pac(ii_zero+delta_ii_sig-1)-t_odor_on;
                    else
                        data_calculated=0;
                    end
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSmPRPtimecourse_peak=mean(SmPRPtimecourse_peak,1)';
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSpPRPtimecourse_peak=mean(SpPRPtimecourse_peak,1)';
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).SmPRPtimecourse_peak=SmPRPtimecourse_peak;
                    handles_out.RPtimecourse(handles_out.RP_ii).SpPRPtimecourse_peak=SpPRPtimecourse_peak;
                    
                    
                    %Trough p values
                    SpPRPtimecourse_trough=[];
                    ii_Sp=0;
                    SmPRPtimecourse_trough=[];
                    ii_Sm=0;
                    for elecNo=which_electrodes
                        for this_session=1:PRPtimecourse.mouse(mouseNo).no_sessions
                            if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed==1
                                %Splus
                                evNo=1;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0 
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).trough)
                                            this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                            this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).trough;
                                            ii_Sp=ii_Sp+1;
                                            SpPRPtimecourse_trough(ii_Sp,:)=this_PRPtimecourse;
                                        end
                                    end
                                end
                                %Sminus
                                evNo=2;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).trough)
                                            this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                            this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).trough;
                                            ii_Sm=ii_Sm+1;
                                            SmPRPtimecourse_trough(ii_Sm,:)=this_PRPtimecourse;
                                        end
                                    end
                                end
                           
                            end
                        end
                    end
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).p_vals_trough=[];
                    handles_out.RPtimecourse(handles_out.RP_ii).t_sig_trough=[];
                    
                    if (size(SpPRPtimecourse_trough,1)>10)&(size(SmPRPtimecourse_trough,1)>10)
                        p_vals_trough=zeros(1,length(anal_t_pac));
                        for ii_t=1:length(anal_t_pac)
                            if use_rs==1
                                p_vals_trough(ii_t)=ranksum(SpPRPtimecourse_trough(:,ii_t),SmPRPtimecourse_trough(:,ii_t));
                            else
                                [h,p_vals_trough(ii_t)]=ttest2(SpPRPtimecourse_trough(:,ii_t),SmPRPtimecourse_trough(:,ii_t));
                            end
                        end
                        
                        handles_out.mouse(mouseNo).per_ii(per_ii).pacii(pacii).p_vals_trough=p_vals_trough;
                        
                        handles_out.RPtimecourse(handles_out.RP_ii).p_vals_trough=p_vals_trough;
                        ii_zero=find(anal_t_pac>=t_odor_on,1,'first');
                        delta_ii_sig=find(p_vals_trough(ii_zero:end)<=0.05,1,'first');
                        handles_out.RPtimecourse(handles_out.RP_ii).t_sig_trough=anal_t_pac(ii_zero+delta_ii_sig-1)-t_odor_on;
                    else
                        data_calculated=0;
                    end
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSmPRPtimecourse_trough=mean(SmPRPtimecourse_trough,1)';
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSpPRPtimecourse_trough=mean(SpPRPtimecourse_trough,1)';
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).SmPRPtimecourse_trough=SmPRPtimecourse_trough;
                    handles_out.RPtimecourse(handles_out.RP_ii).SpPRPtimecourse_trough=SpPRPtimecourse_trough;
                    
                    
                    %Lick power p values
                    SpPRPtimecourse_lickp=[];
                    ii_Sp=0;
                    SmPRPtimecourse_lickp=[];
                    ii_Sm=0;
                    for elecNo=which_electrodes
                        for this_session=1:PRPtimecourse.mouse(mouseNo).no_sessions
                            if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed==1
                                %Splus
                                evNo=1;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'lickp')
                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lickp)
                                                this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                                this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lickp;
                                                if sum(isnan(this_PRPtimecourse))==0
                                                    ii_Sp=ii_Sp+1;
                                                    SpPRPtimecourse_lickp(ii_Sp,:)=this_PRPtimecourse;
                                                end
                                            end
                                        end
                                    end
                                end
                                %Sminus
                                evNo=2;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'lickp')
                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lickp)
                                                this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                                this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lickp;
                                                if sum(isnan(this_PRPtimecourse))==0
                                                    ii_Sm=ii_Sm+1;
                                                    SmPRPtimecourse_lickp(ii_Sm,:)=this_PRPtimecourse;
                                                end
                                            end
                                        end
                                    end
                                end
                                
                            end
                        end
                    end
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).p_vals_lickp=[];
                    handles_out.RPtimecourse(handles_out.RP_ii).t_sig_lickp=[];
                    
                    handles_out.mouse(mouseNo).per_ii(per_ii).pacii(pacii).p_vals_lickp=[];
                    if (size(SpPRPtimecourse_lickp,1)>10)&(size(SmPRPtimecourse_lickp,1)>10)
                        p_vals_lickp=zeros(1,length(anal_t_pac));
                        for ii_t=1:length(anal_t_pac)
                            if use_rs==1
                                p_vals_lickp(ii_t)=ranksum(SpPRPtimecourse_lickp(:,ii_t),SmPRPtimecourse_lickp(:,ii_t));
                            else
                                [h,p_vals_lickp(ii_t)]=ttest2(SpPRPtimecourse_lickp(:,ii_t),SmPRPtimecourse_lickp(:,ii_t));
                            end
                        end
                        
                        handles_out.RPtimecourse(handles_out.RP_ii).p_vals_lickp=p_vals_lickp;
                        ii_zero=find(anal_t_pac>=t_odor_on,1,'first');
                        delta_ii_sig=find(p_vals_lickp(ii_zero:end)<=0.05,1,'first');
                        handles_out.RPtimecourse(handles_out.RP_ii).t_sig_lickp=anal_t_pac(ii_zero+delta_ii_sig-1)-t_odor_on;
                    else
                        data_calculated=0;
                    end
                    
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSmPRPtimecourse_lickp=mean(SmPRPtimecourse_lickp,1)';
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSpPRPtimecourse_lickp=mean(SpPRPtimecourse_lickp,1)';
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).SmPRPtimecourse_lickp=SmPRPtimecourse_lickp;
                    handles_out.RPtimecourse(handles_out.RP_ii).SpPRPtimecourse_lickp=SpPRPtimecourse_lickp;
                    
                    pffft=1;
                    
                    %Lick p values
                    SpPRPtimecourse=[];
                    ii_Sp=0;
                    SmPRPtimecourse=[];
                    ii_Sm=0;
                    %                     for elecNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode)
                    for this_session=1:PRPtimecourse.mouse(mouseNo).no_sessions
                        if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed==1
                            if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii),'event')
                                if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event)
                                    if length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event)>1
                                    %Splus
                                    evNo=1;
                                    if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).no_lick_events>0
                                        for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).lick_tr)
                                            this_lick_bin_timecourse=zeros(1,length(anal_t_pac));
                                            if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).lick_tr(trNo),'binary_lick_per_t')
                                                this_lick_bin_timecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).lick_tr(trNo).binary_lick_per_t;
                                                ii_Sp=ii_Sp+1;
                                                SpPRPtimecourse(ii_Sp,:)=this_lick_bin_timecourse;
                                            end
                                        end
                                    end
                                    %Sminus
                                    evNo=2;
                                    if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).no_lick_events>0
                                        for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).lick_tr)
                                            this_lick_bin_timecourse=zeros(1,length(anal_t_pac));
                                            if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).lick_tr(trNo),'binary_lick_per_t')
                                                this_lick_bin_timecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).lick_tr(trNo).binary_lick_per_t;
                                                ii_Sm=ii_Sm+1;
                                                SmPRPtimecourse(ii_Sm,:)=this_lick_bin_timecourse;
                                            end
                                        end
                                    end
                                end
                                end
                            end
                        end
                    end
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).p_vals_licks=[];
                    handles_out.RPtimecourse(handles_out.RP_ii).t_sig_licks=[];
                    
                    if (size(SpPRPtimecourse,1)>10)&(size(SmPRPtimecourse,1)>10)
                        p_vals_licks=zeros(1,length(anal_t_pac));
                        for ii_t=1:length(anal_t_pac)
                            p_vals_licks(ii_t)=ranksum(SpPRPtimecourse(:,ii_t),SmPRPtimecourse(:,ii_t));
                        end
                        
                        handles_out.RPtimecourse(handles_out.RP_ii).p_vals_licks=p_vals_licks;
                        ii_zero=find(anal_t_pac>=t_odor_on,1,'first');
                        delta_ii_sig=find(p_vals_licks(ii_zero:end)<=0.05,1,'first');
                        handles_out.RPtimecourse(handles_out.RP_ii).t_sig_licks=anal_t_pac(ii_zero+delta_ii_sig-1)-t_odor_on;
                    else
                        data_calculated=0;
                    end
                    
                    
                    %Plot the per mouse figures and calculate the decision making times
                    plot_figure=0;
                    if (data_calculated==1)&(per_ii==1)
                        
                        if plot_figure==1
                            %Plot the p values
                            
                            try
                                close(figureNo)
                            catch
                            end
                            hFig=figure(figureNo);
                            
                            set(hFig, 'units','normalized','position',[.1 .1 .7 .7])
                            %             subplot(2,1,1)
                            ax=gca;ax.LineWidth=3;
                            
                            hold on
                            plot(anal_t_pac,log10(p_vals_licks),'-k','LineWidth',3)
                            %                         plot(anal_t_pac,log10(p_vals_lickp),'-m','LineWidth',3)
                            plot(anal_t_pac,log10(p_vals_trough),'-b','LineWidth',3)
                            plot(anal_t_pac,log10(p_vals_peak),'-r','LineWidth',3)
                            
                            plot([anal_t_pac(1) anal_t_pac(end)],[log10(0.05) log10(0.05)],'-r','LineWidth', 2)
%                             pFDR=drsFDRpval(p_vals_licks);
%                             plot([anal_t_pac(1) anal_t_pac(end)],[log10(pFDR) log10(pFDR)],'-m','LineWidth', 2)
                        end
                        
                        %Find discrimination times and p value timecourses
                        
                        %Licks
                        ii_t0=find(anal_t_pac>=t_odor_on,1,'first');
                        ii=ii_t0+1;
                        not_found=1;
                        ii_cross=ii;
                        delta_ii_dt=ceil(0.3/0.0333);
                        while (ii<=ii_t0+50)&(not_found==1)
                            if (log10(p_vals_licks(ii-1))>=log10(0.05))&(log10(p_vals_licks(ii))<=log10(0.05))
                               if sum(log10(p_vals_licks(ii:ii+delta_ii_dt))>log10(0.05))==0
                                   ii_cross=ii;
                                   not_found=0;
                               end
                            end
                            ii=ii+1;
                        end
                        
                        if not_found==0
                            group_no=PRPtimecourse.mouse(mouseNo).group_no;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).no_licks=handles_out.disc_t.pacii(pacii).groupNo(PRPtimecourse.mouse(mouseNo).group_no).no_licks+1;
                            no_licks=handles_out.disc_t.pacii(pacii).groupNo(group_no).no_licks;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).licks(no_licks).disc_t=anal_t_pac(ii_cross)-t_odor_on;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).licks(no_licks).p_vals=log10(p_vals_licks);
                            
                            if plot_figure==1
                                plot(anal_t_pac(ii_cross),log10(0.05),'ok','MarkerFaceColor','k','MarkerSize',8)
                            end
                        end
                         
                        %Peak PRP
                        ii_t0=find(anal_t_pac>=t_odor_on,1,'first');
                        ii=ii_t0+1;
                        not_found=1;
                        ii_cross=ii;
                        delta_ii_dt=ceil(0.3/0.0333);
                        while (ii<=ii_t0+50)&(not_found==1)
                            if (log10(p_vals_peak(ii-1))>=log10(0.05))&(log10(p_vals_peak(ii))<=log10(0.05))
                               if sum(log10(p_vals_peak(ii:ii+delta_ii_dt))>log10(0.05))==0
                                   ii_cross=ii;
                                   not_found=0;
                               end
                            end
                            ii=ii+1;
                        end
                        
                        if not_found==0
                            group_no=PRPtimecourse.mouse(mouseNo).group_no;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).no_peak=handles_out.disc_t.pacii(pacii).groupNo(PRPtimecourse.mouse(mouseNo).group_no).no_peak+1;
                            no_peak=handles_out.disc_t.pacii(pacii).groupNo(group_no).no_peak;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).peak(no_peak).disc_t=anal_t_pac(ii_cross)-t_odor_on;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).peak(no_peak).p_vals=log10(p_vals_peak);
                            
                            if plot_figure==1
                                plot(anal_t_pac(ii_cross),log10(0.05),'or','MarkerFaceColor','r','MarkerSize',10)
                            end
                            
                        end
                        
                        %Trough PRP
                        ii_t0=find(anal_t_pac>=t_odor_on,1,'first');
                        ii=ii_t0+1;
                        not_found=1;
                        ii_cross=ii;
                        delta_ii_dt=ceil(0.3/0.0333);
                        while (ii<=ii_t0+50)&(not_found==1)
                            if (log10(p_vals_trough(ii-1))>=log10(0.05))&(log10(p_vals_trough(ii))<=log10(0.05))
                               if sum(log10(p_vals_trough(ii:ii+delta_ii_dt))>log10(0.05))==0
                                   ii_cross=ii;
                                   not_found=0;
                               end
                            end
                            ii=ii+1;
                        end  
                        
                        if not_found==0
                            group_no=PRPtimecourse.mouse(mouseNo).group_no;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).no_trough=handles_out.disc_t.pacii(pacii).groupNo(PRPtimecourse.mouse(mouseNo).group_no).no_trough+1;
                            no_trough=handles_out.disc_t.pacii(pacii).groupNo(group_no).no_trough;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).trough(no_trough).disc_t=anal_t_pac(ii_cross)-t_odor_on;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).trough(no_trough).p_vals=log10(p_vals_trough);
                            
                            if plot_figure==1
                                plot(anal_t_pac(ii_cross),log10(0.05),'ob','MarkerFaceColor','b','MarkerSize',7)
                            end
                        end
                        
                        if plot_figure==1
                            xlabel('Time (sec)')
                            ylabel('log10(p)')
                             xlim([-0.4 0.7])
                            ylim([-20 0])
                            title(['z score t test log10(p) for mouse no ' num2str(mouseNo) ' ' freq_names{pacii+1} ' ' prof_naive_leg{per_ii} ' group No ' num2str(PRPtimecourse.mouse(mouseNo).group_no)])
                            
                            %Plot the bounded line averages
                             
                            try
                                close(figureNo+1)
                            catch
                            end
                            hFig=figure(figureNo+1);
                            
                            set(hFig, 'units','normalized','position',[.1 .2 .8 .3])
                            
                            %Show the average peak PRP dynamics
                            %S+, S-
                            subplot(1,3,1)
                            hold on
                            ax=gca;ax.LineWidth=3;
                            
                            CIsm = bootci(1000, @mean, SmPRPtimecourse_peak);
                            meansm=mean(SmPRPtimecourse_peak,1);
                            CIsm(1,:)=meansm-CIsm(1,:);
                            CIsm(2,:)=CIsm(2,:)-meansm;
                            
                            CIsp = bootci(1000, @mean, SpPRPtimecourse_peak);
                            meansp=mean(SpPRPtimecourse_peak,1);
                            CIsp(1,:)=meansp-CIsp(1,:);
                            CIsp(2,:)=CIsp(2,:)-meansp;
                            
                            
                            [hlsm, hpsm] = boundedline(anal_t_pac',mean(SmPRPtimecourse_peak,1)', CIsm', 'b');
                            [hlsp, hpsp] = boundedline(anal_t_pac',mean(SpPRPtimecourse_peak,1)', CIsp', 'r');
                            
                            ylim([-1.6 0.8])
                            xlabel('Time (sec)')
                            ylabel('PRP z score')
                            
                            
                            title(['PRP peak z score for mouse no ' num2str(mouseNo) ' ' freq_names{pacii+1} ' ' prof_naive_leg{per_ii} ' group No ' num2str(PRPtimecourse.mouse(mouseNo).group_no)])
                            
                            %Show the average trough PRP dynamics
                            %S+, S-
                            subplot(1,3,2)
                            hold on
                            ax=gca;ax.LineWidth=3;
                            
                            CIsm = bootci(1000, @mean, SmPRPtimecourse_trough);
                            meansm=mean(SmPRPtimecourse_trough,1);
                            CIsm(1,:)=meansm-CIsm(1,:);
                            CIsm(2,:)=CIsm(2,:)-meansm;
                            
                            CIsp = bootci(1000, @mean, SpPRPtimecourse_trough);
                            meansp=mean(SpPRPtimecourse_trough,1);
                            CIsp(1,:)=meansp-CIsp(1,:);
                            CIsp(2,:)=CIsp(2,:)-meansp;
                            
                            
                            [hlsm, hpsm] = boundedline(anal_t_pac',mean(SmPRPtimecourse_trough,1)', CIsm', 'b');
                            [hlsp, hpsp] = boundedline(anal_t_pac',mean(SpPRPtimecourse_trough,1)', CIsp', 'r');
                            
                            ylim([-1.6 0.8])
                            xlabel('Time (sec)')
                            ylabel('PRP z score')
                            
                            
                            title(['PRP trough z score for mouse no ' num2str(mouseNo) ' ' freq_names{pacii+1} ' ' prof_naive_leg{per_ii} ' group No ' num2str(PRPtimecourse.mouse(mouseNo).group_no)])
                            
                            %Show the average LRP dynamics
                            %S+, S-
                            subplot(1,3,3)
                            hold on
                            ax=gca;ax.LineWidth=3;
                            
                            CIsm = bootci(1000, @mean, SmPRPtimecourse_lickp);
                            meansm=mean(SmPRPtimecourse_lickp,1);
                            CIsm(1,:)=meansm-CIsm(1,:);
                            CIsm(2,:)=CIsm(2,:)-meansm;
                            
                            CIsp = bootci(1000, @mean, SpPRPtimecourse_lickp);
                            meansp=mean(SpPRPtimecourse_lickp,1);
                            CIsp(1,:)=meansp-CIsp(1,:);
                            CIsp(2,:)=CIsp(2,:)-meansp;
                            
                            
                            [hlsm, hpsm] = boundedline(anal_t_pac',mean(SmPRPtimecourse_lickp,1)', CIsm', 'b');
                            [hlsp, hpsp] = boundedline(anal_t_pac',mean(SpPRPtimecourse_lickp,1)', CIsp', 'r');
                            
                            ylim([-1.6 0.8])
                            xlabel('Time (sec)')
                            ylabel('LRP z score')
                            
                            
                            title(['LRP z score for mouse no ' num2str(mouseNo) ' ' freq_names{pacii+1} ' ' prof_naive_leg{per_ii} ' group No ' num2str(PRPtimecourse.mouse(mouseNo).group_no)])
                        end
                         pffft=1;
                    end
                   
                end
            end
        end
        
        pffft=1;
        
        %Peak PRP
        for pacii=1:no_pacii
            figureNo = figureNo + 1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            for per_ii=2:-1:1
                
                for groupNo=1:3
                    subplot(2,3,groupNo + 3*(2-per_ii))
                    hold on
                    
                    all_meanSmPRPtimecourse_peak=[];
                    all_meanSpPRPtimecourse_peak=[];
                    ii_tcs=0;
                    for ii=1:handles_out.RP_ii
                        if (handles_out.RPtimecourse(ii).pacii==pacii)&(handles_out.RPtimecourse(ii).per_ii==per_ii)&(handles_out.RPtimecourse(ii).group_no==groupNo)
                            if (length(handles_out.RPtimecourse(ii).meanSmPRPtimecourse_peak)>1)&(length(handles_out.RPtimecourse(ii).meanSpPRPtimecourse_peak)>1)
                                ii_tcs=ii_tcs+1;
                                all_meanSmPRPtimecourse_peak(ii_tcs,:)=handles_out.RPtimecourse(ii).meanSmPRPtimecourse_peak;
                                all_meanSpPRPtimecourse_peak(ii_tcs,:)=handles_out.RPtimecourse(ii).meanSpPRPtimecourse_peak;
                            end
                        end
                    end
                    fprintf(1, ['The number of mice included for ' prof_naive_leg{per_ii} ' ' group_legend{groupNo} ' is %d\n'], ii_tcs)
                    
                    CIsp = bootci(1000, @mean, all_meanSpPRPtimecourse_peak);
                    meansp=mean(all_meanSpPRPtimecourse_peak,1);
                    CIsp(1,:)=meansp-CIsp(1,:);
                    CIsp(2,:)=CIsp(2,:)-meansp;
                    
                    
                    [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_meanSpPRPtimecourse_peak,1)', CIsp', 'r');
                    
                    CIsm = bootci(1000, @mean, all_meanSmPRPtimecourse_peak);
                    meansm=mean(all_meanSmPRPtimecourse_peak,1);
                    CIsm(1,:)=meansm-CIsm(1,:);
                    CIsm(2,:)=CIsm(2,:)-meansm;
                    
                    [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_meanSmPRPtimecourse_peak,1)', CIsm', 'b');
                    
                    title([group_legend{groupNo} ' ' prof_naive_leg{per_ii}])
                    xlabel('Time(sec)')
                    ylabel('z')
                    ylim([-1.5 1])
                    
                end
            end
            
            sgtitle(['Peak PRP ' freq_names{pacii+1}])
        end
        
        %Trough PRP
        for pacii=1:no_pacii
            figureNo = figureNo + 1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            for per_ii=2:-1:1
                
                for groupNo=1:3
                    subplot(2,3,groupNo + 3*(2-per_ii))
                    hold on
                    
                    all_meanSmPRPtimecourse_trough=[];
                    all_meanSpPRPtimecourse_trough=[];
                    ii_tcs=0;
                    for ii=1:handles_out.RP_ii
                        if (handles_out.RPtimecourse(ii).pacii==pacii)&(handles_out.RPtimecourse(ii).per_ii==per_ii)&(handles_out.RPtimecourse(ii).group_no==groupNo)
                            if (length(handles_out.RPtimecourse(ii).meanSmPRPtimecourse_trough)>1)&(length(handles_out.RPtimecourse(ii).meanSpPRPtimecourse_trough)>1)
                                ii_tcs=ii_tcs+1;
                                all_meanSmPRPtimecourse_trough(ii_tcs,:)=handles_out.RPtimecourse(ii).meanSmPRPtimecourse_trough;
                                all_meanSpPRPtimecourse_trough(ii_tcs,:)=handles_out.RPtimecourse(ii).meanSpPRPtimecourse_trough;
                            end
                        end
                    end
                    fprintf(1, ['The number of mice included for ' prof_naive_leg{per_ii} ' ' group_legend{groupNo} ' is %d\n'], ii_tcs)
                    
                    CIsp = bootci(1000, @mean, all_meanSpPRPtimecourse_trough);
                    meansp=mean(all_meanSpPRPtimecourse_trough,1);
                    CIsp(1,:)=meansp-CIsp(1,:);
                    CIsp(2,:)=CIsp(2,:)-meansp;
                    
                    
                    [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_meanSpPRPtimecourse_trough,1)', CIsp', 'r');
                    
                    CIsm = bootci(1000, @mean, all_meanSmPRPtimecourse_trough);
                    meansm=mean(all_meanSmPRPtimecourse_trough,1);
                    CIsm(1,:)=meansm-CIsm(1,:);
                    CIsm(2,:)=CIsm(2,:)-meansm;
                    
                    [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_meanSmPRPtimecourse_trough,1)', CIsm', 'b');
                    
                    title([group_legend{groupNo} ' ' prof_naive_leg{per_ii}])
                    xlabel('Time(sec)')
                    ylabel('z')
                    ylim([-1.5 1])
                    
                end
            end
            
            sgtitle(['Trough PRP ' freq_names{pacii+1}])
        end
        
        %LRP
        for pacii=1:no_pacii
            figureNo = figureNo + 1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            for per_ii=2:-1:1
                
                for groupNo=1:3
                    subplot(2,3,groupNo + 3*(2-per_ii))
                    hold on
                    
                    all_meanSmPRPtimecourse_lickp=[];
                    all_meanSpPRPtimecourse_lickp=[];
                    ii_tcs=0;
                    for ii=1:handles_out.RP_ii
                        if (handles_out.RPtimecourse(ii).pacii==pacii)&(handles_out.RPtimecourse(ii).per_ii==per_ii)&(handles_out.RPtimecourse(ii).group_no==groupNo)
                            if (length(handles_out.RPtimecourse(ii).meanSmPRPtimecourse_lickp)>1)&(length(handles_out.RPtimecourse(ii).meanSpPRPtimecourse_lickp)>1)
                                ii_tcs=ii_tcs+1;
                                all_meanSmPRPtimecourse_lickp(ii_tcs,:)=handles_out.RPtimecourse(ii).meanSmPRPtimecourse_lickp;
                                all_meanSpPRPtimecourse_lickp(ii_tcs,:)=handles_out.RPtimecourse(ii).meanSpPRPtimecourse_lickp;
                            end
                        end
                    end
                    fprintf(1, ['The number of mice included for ' prof_naive_leg{per_ii} ' ' group_legend{groupNo} ' is %d\n'], ii_tcs)
                    
                    CIsp = bootci(1000, @mean, all_meanSpPRPtimecourse_lickp);
                    meansp=mean(all_meanSpPRPtimecourse_lickp,1);
                    CIsp(1,:)=meansp-CIsp(1,:);
                    CIsp(2,:)=CIsp(2,:)-meansp;
                    
                    
                    [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_meanSpPRPtimecourse_lickp,1)', CIsp', 'r');
                    
                    CIsm = bootci(1000, @mean, all_meanSmPRPtimecourse_lickp);
                    meansm=mean(all_meanSmPRPtimecourse_lickp,1);
                    CIsm(1,:)=meansm-CIsm(1,:);
                    CIsm(2,:)=CIsm(2,:)-meansm;
                    
                    [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_meanSmPRPtimecourse_lickp,1)', CIsm', 'b');
                    
                    title([group_legend{groupNo} ' ' prof_naive_leg{per_ii}])
                    xlabel('Time(sec)')
                    ylabel('z')
                    ylim([-1.5 1])
                    
                end
            end
            
            sgtitle(['Lick RP ' freq_names{pacii+1}])
        end
        
        %Now let's plot the p value timecourses
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        for pacii=1:no_pacii
            
            
            per_ii=1;
                
                for groupNo=1:3
                    subplot(3,3,groupNo + 3*(3-pacii))
                    hold on
                    
                    
                    all_mean_p_vals_peak=[];
                    no_peak=handles_out.disc_t.pacii(pacii).groupNo(groupNo).no_peak
                    for ii=1:no_peak
                        this_pval=handles_out.disc_t.pacii(pacii).groupNo(groupNo).peak(ii).p_vals;
                        log10_this_pval=zeros(1,length(this_pval));
                        log10_this_pval(1,:)=this_pval;
                        log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                        all_mean_p_vals_peak(ii,:)=log10_this_pval;
                    end
                    fprintf(1, ['The number of mice included for peak pval ' prof_naive_leg{per_ii} ' ' group_legend{groupNo} ' ' freq_names{pacii+1} ' is %d\n'], ii)
                    
                    all_mean_p_vals_trough=[];
                    no_trough=handles_out.disc_t.pacii(pacii).groupNo(groupNo).no_trough
                    for ii=1:no_trough
                        this_pval=handles_out.disc_t.pacii(pacii).groupNo(groupNo).trough(ii).p_vals;
                        log10_this_pval=zeros(1,length(this_pval));
                        log10_this_pval(1,:)=this_pval;
                        log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                        all_mean_p_vals_trough(ii,:)=log10_this_pval;
                    end
                    fprintf(1, ['The number of mice included for trough pval ' prof_naive_leg{per_ii} ' ' group_legend{groupNo} ' ' freq_names{pacii+1} ' is %d\n'], ii)
                    
                    all_mean_p_vals_licks=[];
                    no_licks=handles_out.disc_t.pacii(pacii).groupNo(groupNo).no_licks;
                    for ii=1:no_licks
                        this_pval=handles_out.disc_t.pacii(pacii).groupNo(groupNo).licks(ii).p_vals;
                        if sum(isnan(this_pval))~=0
                            ii_nan=find(isnan(this_pval));
                            for jj=ii_nan
                                this_pval(jj)=this_pval(jj-1);
                            end
                        end
                        log10_this_pval=zeros(1,length(this_pval));
                        log10_this_pval(1,:)=this_pval;
                        log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                        all_mean_p_vals_licks(ii,:)=log10_this_pval;
                    end
                    fprintf(1, ['The number of mice included for lick pval ' prof_naive_leg{per_ii} ' ' group_legend{groupNo} ' ' freq_names{pacii+1} ' is %d\n'], ii)
                    
                    fprintf(1,'\n')
                     
                  
                    CIpv = bootci(1000, @mean, all_mean_p_vals_licks);
                    meanpv=mean(all_mean_p_vals_licks,1);
                    CIpv(1,:)=meanpv-CIpv(1,:);
                    CIpv(2,:)=CIpv(2,:)-meanpv;
                    
                    
                    [hlpvl, hppvl] = boundedline(anal_t_pac',mean(all_mean_p_vals_licks,1)', CIpv', 'k');
                    
                   
                    CIpv = bootci(1000, @mean, all_mean_p_vals_trough);
                    meanpv=mean(all_mean_p_vals_trough,1);
                    CIpv(1,:)=meanpv-CIpv(1,:);
                    CIpv(2,:)=CIpv(2,:)-meanpv;
                    
                    
                    [hlpvt, hppvt] = boundedline(anal_t_pac',mean(all_mean_p_vals_trough,1)', CIpv', 'b');
                    
                     CIpv = bootci(1000, @mean, all_mean_p_vals_peak);
                    meanpv=mean(all_mean_p_vals_peak,1);
                    CIpv(1,:)=meanpv-CIpv(1,:);
                    CIpv(2,:)=CIpv(2,:)-meanpv;
                    
                    
                    [hlpvp, hppvp] = boundedline(anal_t_pac',mean(all_mean_p_vals_peak,1)', CIpv', 'r');
                  
                    
                    plot(anal_t_pac',mean(all_mean_p_vals_licks,1)',  'k');
                    plot(anal_t_pac',mean(all_mean_p_vals_trough,1)',  'b');
                    plot(anal_t_pac',mean(all_mean_p_vals_peak,1)',  'r');
                    
                    plot([anal_t_pac(1) anal_t_pac(end)],[log10(0.05) log10(0.05)],'-r','LineWidth', 2)
                    
                    ylim([-400 0])
                    title([group_legend{groupNo} ' ' freq_names{pacii+1}])
                    xlabel('Time(sec)')
                    ylabel('log10(p)')

                end
                fprintf(1, '\n')
                
                
        end

        sgtitle('p values for proficient mice')
        
        if pacii==3
            %Plot the p value mean for licks
            figureNo = figureNo + 1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            ax=gca;ax.LineWidth=3;
            set(hFig, 'units','normalized','position',[.3 .3 .4 .25])

            hold on
            
            CIpv = bootci(1000, @mean, all_mean_p_vals_licks);
            meanpv=mean(all_mean_p_vals_licks,1);
            CIpv(1,:)=meanpv-CIpv(1,:);
            CIpv(2,:)=CIpv(2,:)-meanpv;
            
            
            [hlpvl, hppvl] = boundedline(anal_t_pac',mean(all_mean_p_vals_licks,1)', CIpv', 'k');
            
            plot(anal_t_pac',mean(all_mean_p_vals_licks,1)','-k','LineWidth',3);
            
            plot([anal_t_pac(1) anal_t_pac(end)],[log10(0.05) log10(0.05)],'-r','LineWidth', 2)
            
            plot([0 0],[-10 0],'-k','LineWidth',2)
            
            ylim([-10 0])
            title('Lick p value')
            xlabel('Time(sec)')
            ylabel('log10(p)')
            
            handles_out.lick_rate.all_mean_p_vals_licks_p=all_mean_p_vals_licks;
            handles_out.lick_rate.mean_all_mean_p_vals_licks_p=mean(all_mean_p_vals_licks,1)';
            handles_out.lick_rate.anal_t_pac=anal_t_pac';
           
            pffft=1;
        end
        
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        for pacii=1:no_pacii
            
            
            per_ii=1;
                
                for groupNo=1:3
                    subplot(3,3,groupNo + 3*(3-pacii))
                    hold on
                    
                    
                    all_mean_p_vals_peak=[];
                    no_peak=handles_out.disc_t.pacii(pacii).groupNo(groupNo).no_peak
                    for ii=1:no_peak
                        this_pval=handles_out.disc_t.pacii(pacii).groupNo(groupNo).peak(ii).p_vals;
                        log10_this_pval=zeros(1,length(this_pval));
                        log10_this_pval(1,:)=this_pval;
                        log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                        all_mean_p_vals_peak(ii,:)=log10_this_pval;
                    end
                    fprintf(1, ['The number of mice included for peak pval ' prof_naive_leg{per_ii} ' ' group_legend{groupNo} ' ' freq_names{pacii+1} ' is %d\n'], ii)
                    
                    all_mean_p_vals_trough=[];
                    no_trough=handles_out.disc_t.pacii(pacii).groupNo(groupNo).no_trough
                    for ii=1:no_trough
                        this_pval=handles_out.disc_t.pacii(pacii).groupNo(groupNo).trough(ii).p_vals;
                        log10_this_pval=zeros(1,length(this_pval));
                        log10_this_pval(1,:)=this_pval;
                        log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                        all_mean_p_vals_trough(ii,:)=log10_this_pval;
                    end
                    fprintf(1, ['The number of mice included for trough pval ' prof_naive_leg{per_ii} ' ' group_legend{groupNo} ' ' freq_names{pacii+1} ' is %d\n'], ii)
                    
                    all_mean_p_vals_licks=[];
                    no_licks=handles_out.disc_t.pacii(pacii).groupNo(groupNo).no_licks;
                    for ii=1:no_licks
                        this_pval=handles_out.disc_t.pacii(pacii).groupNo(groupNo).licks(ii).p_vals;
                        log10_this_pval=zeros(1,length(this_pval));
                        log10_this_pval(1,:)=this_pval;
                        log10_this_pval(isinf(log10_this_pval))=min(log10_this_pval(~isinf(log10_this_pval)));
                        all_mean_p_vals_licks(ii,:)=log10_this_pval;
                    end
                    fprintf(1, ['The number of mice included for lick pval ' prof_naive_leg{per_ii} ' ' group_legend{groupNo} ' ' freq_names{pacii+1} ' is %d\n'], ii)
                    
                    fprintf(1,'\n')
                    
                    %                     CIpv = bootci(1000, @mean, all_mean_p_vals_licks);
                    %                     meanpv=mean(all_mean_p_vals_licks,1);
                    %                     CIpv(1,:)=meanpv-CIpv(1,:);
                    %                     CIpv(2,:)=CIpv(2,:)-meanpv;
                    %
                    %
                    %                     [hlpvl, hppvl] = boundedline(anal_t_pac',mean(all_mean_p_vals_licks,1)', CIpv', 'k');
                    %
                    %
                    %                     CIpv = bootci(1000, @mean, all_mean_p_vals_trough);
                    %                     meanpv=mean(all_mean_p_vals_trough,1);
                    %                     CIpv(1,:)=meanpv-CIpv(1,:);
                    %                     CIpv(2,:)=CIpv(2,:)-meanpv;
                    %
                    %
                    %                     [hlpvt, hppvt] = boundedline(anal_t_pac',mean(all_mean_p_vals_trough,1)', CIpv', 'b');
                    %
                    %                      CIpv = bootci(1000, @mean, all_mean_p_vals_peak);
                    %                     meanpv=mean(all_mean_p_vals_peak,1);
                    %                     CIpv(1,:)=meanpv-CIpv(1,:);
                    %                     CIpv(2,:)=CIpv(2,:)-meanpv;
                    %
                    %
                    %                     [hlpvp, hppvp] = boundedline(anal_t_pac',mean(all_mean_p_vals_peak,1)', CIpv', 'r');
                    
                    
                    plot(anal_t_pac',mean(all_mean_p_vals_licks,1)',  'k');
                    plot(anal_t_pac',mean(all_mean_p_vals_trough,1)',  'b');
                    plot(anal_t_pac',mean(all_mean_p_vals_peak,1)',  'r');
                    
                    plot([anal_t_pac(1) anal_t_pac(end)],[log10(0.05) log10(0.05)],'-r','LineWidth', 2)
                    
                    ylim([-5 0])
                    xlim([-0.05 0.5966])
                    title([group_legend{groupNo} ' ' freq_names{pacii+1}])
                    xlabel('Time(sec)')
                    ylabel('log10(p)')
                    
                end
                fprintf(1, '\n')
                
                
        end
        
        sgtitle('p values for proficient mice closeup')
        
        %Now let's plot decision times
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        for pacii=1:no_pacii
            
            
            per_ii=1;
            
            for groupNo=1:3
                subplot(3,3,groupNo + 3*(pacii-1))
                hold on
                
                all_mean_decision_dt_peak=[];
                all_mean_decision_dt_trough=[];
                all_mean_decision_dt_licks=[];
                all_mean_decision_dt_lickp=[];
                
                
                all_decision_dt_licks=[];
                for ii=1:handles_out.disc_t.pacii(pacii).groupNo(group_no).no_licks
                    all_decision_dt_licks(ii)=handles_out.disc_t.pacii(pacii).groupNo(group_no).licks(ii).disc_t;
                end
                
                all_decision_dt_peak=[];
                for ii=1:handles_out.disc_t.pacii(pacii).groupNo(group_no).no_peak
                    all_decision_dt_peak(ii)=handles_out.disc_t.pacii(pacii).groupNo(group_no).peak(ii).disc_t;
                end
                
                all_decision_dt_trough=[];
                for ii=1:handles_out.disc_t.pacii(pacii).groupNo(group_no).no_trough
                    all_decision_dt_trough(ii)=handles_out.disc_t.pacii(pacii).groupNo(group_no).trough(ii).disc_t;
                end
                
                
                edges=[0:0.033:0.5];
                rand_offset=0.8;
                
                %Peak PRP
                bar_offset=1;
                bar(bar_offset,mean(all_decision_dt_peak),'r','LineWidth', 3,'EdgeColor','none')
                
                %Violin plot
                [mean_out, CIout]=drgViolinPoint(all_decision_dt_peak...
                    ,edges,bar_offset,rand_offset,'k','k',3);
                
                %Trough PRP
                bar_offset=bar_offset+1;
                bar(bar_offset,mean(all_decision_dt_trough),'b','LineWidth', 3,'EdgeColor','none')
                
                %Violin plot
                [mean_out, CIout]=drgViolinPoint(all_decision_dt_trough...
                    ,edges,bar_offset,rand_offset,'k','k',3);
                
                
                %Licks
                bar_offset=bar_offset+1;
                bar(bar_offset,mean(all_decision_dt_licks),'FaceColor',[0.7 0.7 0.7],'LineWidth', 3,'EdgeColor','k')
                
                %Violin plot
                [mean_out, CIout]=drgViolinPoint(all_decision_dt_licks...
                    ,edges,bar_offset,rand_offset,'k','k',3);
                
                
                title([group_legend{groupNo} ' ' freq_names{pacii+1}])
                
                ylabel('dt')
                
                
            end
            fprintf(1, '\n')
        end
        
        sgtitle(['Decision time (sec)'])
          
        if use_rs==1
            save([handles.PathName handles.drgb.outFileName(1:end-4) '_rs_case' num2str(which_display) '_' handles_pars.output_suffix],'handles_out')
        else
            save([handles.PathName handles.drgb.outFileName(1:end-4) '_case' num2str(which_display) '_' handles_pars.output_suffix],'handles_out')
        end
        
        
    case 3
        %1 Oscillatory wavelet power timecourse divergence after odorant addition
        
        t_odor_on=0.1045;  %From Supplementary Fig 1 of Losacco et al 2020 DOI: 10.7554/eLife.52583
%         handles_pars.lick_analysis_times=[0.5 2.5]; %Used by case 3
        lick_analysis_window=0.3; %2;
        lick_analysis_dt=0.005; %0.02;
        edges=[-lick_analysis_window:lick_analysis_dt:lick_analysis_window];
        
        win_lick_start=[-1.5 0.1 0.6];
        win_lick_end=[0 0.6 2.6];
        
        handles_out=[];
        handles_out.PRP_ii=0;
        handles_out.LickF_ii=0;
        handles_out.AUC_ii=0;
        handles_out.AUClick_ii=0;
        
        
        group_legend{1}='WT';
        group_legend{2}='Het';
        group_legend{3}='KO';
        
        prof_naive_leg{1}='Proficient';
        prof_naive_leg{2}='Naive';
        
        
        
        PRPtimecourse=[];
        
        for ii_mouse=1:max(handles_drgb.drgbchoices.mouse_no)
            PRPtimecourse.mouse(ii_mouse).no_sessions=0;
            PRPtimecourse.mouse(ii_mouse).files=[];
        end
        
        fprintf(1, ['PAC power analysis using wavelet power for Justin''s paper\n\n'])
        p_vals=[];
        no_files=max(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        
        szpc=size(percent_windows);
        for per_ii=1:szpc(1)
            
            
            
            for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
                mouse_has_files=0;
                for elec=1:16
                    if sum(which_electrodes==elec)>0
                        
                        for fileNo=1:no_files
                            %If this file is in the list of files the user wants to process in drgAnalysisBatchLFP continue
                            if sum(files==fileNo)>0
                                if handles_drgb.drgbchoices.mouse_no(fileNo)==mouseNo
                                    
                                    if PRPtimecourse.mouse(mouseNo).no_sessions==0
                                        PRPtimecourse.mouse(mouseNo).no_sessions=1;
                                        PRPtimecourse.mouse(mouseNo).files(1)=fileNo;
                                        this_session=1;
                                    else
                                        this_session=find(PRPtimecourse.mouse(mouseNo).files==fileNo);
                                        if isempty(this_session)
                                            PRPtimecourse.mouse(mouseNo).no_sessions=PRPtimecourse.mouse(mouseNo).no_sessions+1;
                                            PRPtimecourse.mouse(mouseNo).files(PRPtimecourse.mouse(mouseNo).no_sessions)=fileNo;
                                            this_session=PRPtimecourse.mouse(mouseNo).no_sessions;
                                        end
                                    end
                                    
                                    PRPtimecourse.mouse(mouseNo).group_no=handles_drgb.drgbchoices.group_no(fileNo);
                                    
                                    lfpodNo=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                    %                                     lfpodRefNo=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                                    if (mouseNo==2)&(this_session==4)
                                        pffft=1;
                                    end
                                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo)))
                                        
                                        
                                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo).PAC))
                                            
                                            percent_mask=[];
                                            percent_mask=logical((handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).perCorrPAC>=percent_windows(per_ii,1))...
                                                &(handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).perCorrPAC<=percent_windows(per_ii,2)));
                                            
                                            if ~isempty(percent_mask)
                                                
                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed=1;
                                                
                                                for evNo=1:length(eventType)
                                                    
                                                    trials_in_event_Ev=[];
                                                    trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).which_eventPAC(eventType(evNo),:)==1)&percent_mask;
                                                    for pacii=1:no_pacii
                                                        PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).no_events=sum(trials_in_event_Ev);
                                                    end
                                                    
                                                    if (sum(trials_in_event_Ev)>=1)
                                                        
                                                        trialNos_in_event_Ev=zeros(1,sum(trials_in_event_Ev));
                                                        ii_no=0;
                                                        for ii=1:length(trials_in_event_Ev)
                                                            if trials_in_event_Ev(ii)==1
                                                                ii_no=ii_no+1;
                                                                trialNos_in_event_Ev(ii_no)=ii;
                                                            end
                                                        end
                                                        
                                                        %Do per bandwidth analysis
                                                        for pacii=1:no_pacii
                                                            
                                                            group_no=handles_drgb.drgbchoices.group_no(fileNo);
                                                            
                                                            %Enter the PACpower
                                                            t_pac=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.t_pac;
                                                            PRPtimecourse.t_pac=t_pac((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
                                                             
                                                            %Get peak power
                                                            these_ref_power=[];
                                                            for ww=1:length(trialNos_in_event_Ev)
                                                                this_p=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_event_Ev(ww)).peakPower;
                                                                these_ref_power=[these_ref_power this_p((t_pac>=handles_pars.refWin_times(1))&(t_pac<=handles_pars.refWin_times(2)))];
                                                            end
                                                            ref_power_SD=std(these_ref_power);
                                                            ref_power_mean=mean(these_ref_power);
                                                            

                                                            for ww=1:length(trialNos_in_event_Ev)
                                                                this_p=(handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_event_Ev(ww)).peakPower-ref_power_mean)/ref_power_SD;
                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peak=this_p((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
                                                            end
                                                            PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).no_z_power=sum(trials_in_event_Ev);
                                                             
                                                            %Save the times for peaks
                                                            for ww=1:length(trialNos_in_event_Ev)
                                                                these_peak_times=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_event_Ev(ww)).peakPower_times;
                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peak_times=these_peak_times;
                                                            end
                                                            
                                                            %Get trough power
                                                            these_ref_power=[];
                                                            for ww=1:length(trialNos_in_event_Ev)
                                                                this_p=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_event_Ev(ww)).troughPower;
                                                                these_ref_power=[these_ref_power this_p((t_pac>=handles_pars.refWin_times(1))&(t_pac<=handles_pars.refWin_times(2)))];
                                                            end
                                                            ref_power_SD=std(these_ref_power);
                                                            ref_power_mean=mean(these_ref_power);
                                                            
                                                            for ww=1:length(trialNos_in_event_Ev)
                                                                this_p=(handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_event_Ev(ww)).troughPower-ref_power_mean)/ref_power_SD;
                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).trough=this_p((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
                                                            end
                                                            
                                                            %Save the times for troughs
                                                            for ww=1:length(trialNos_in_event_Ev)
                                                                these_trough_times=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_event_Ev(ww)).troughPower_times;
                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).trough_times=these_trough_times;
                                                            end
                                                            
                                                            
                                                            %Now do the lick power and other lick computations
                                                            if isfield(handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(1).PACwave,'lickPower_trials')
                                                                if length(handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(1).PACwave.lickPower_trials)>0
                                                                    
                                                                    lick_trials=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.lickPower_trials;
                                                                    
                                                                    %Do lick referenced power
                                                                    
                                                                    these_lick_trialNos=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(1).PACwave.lickPower_trials;
                                                                    these_lick_trialNos=these_lick_trialNos(these_lick_trialNos<=length(trials_in_event_Ev));
                                                                    trials_in_lick_event=zeros(1,length(these_lick_trialNos));
                                                                    trials_in_lick_event(1,:)=trials_in_event_Ev(these_lick_trialNos);
                                                                    
                                                                    trialNos_in_lick_event_Ev=zeros(1,sum(trials_in_lick_event));
                                                                    ii_no=0;
                                                                    for ii=1:length(trials_in_lick_event)
                                                                        if trials_in_lick_event(ii)==1
                                                                            ii_no=ii_no+1;
                                                                            trialNos_in_lick_event_Ev(ii_no)=ii;
                                                                        end
                                                                    end
                                                                    
                                                                    these_ref_power=[];
                                                                    for ww=1:length(trialNos_in_lick_event_Ev)
                                                                        this_p=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_lick_event_Ev(ww)).lickPower;
                                                                        these_ref_power=[these_ref_power this_p((t_pac>=handles_pars.refWin_times(1))&(t_pac<=handles_pars.refWin_times(2)))];
                                                                    end
                                                                    ref_power_SD=std(these_ref_power);
                                                                    ref_power_mean=mean(these_ref_power);
                                                                    
                                                                    
                                                                    for ww=1:length(trialNos_in_lick_event_Ev)
                                                                        this_p=(handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_lick_event_Ev(ww)).lickPower-ref_power_mean)/ref_power_SD;
                                                                        PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).lickp=this_p((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
                                                                    end
                                                                    PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).no_z_lick_power=sum(trials_in_lick_event);
                                                                    PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).lick_trials=lick_trials;
                                                                    
                                                                    %Do lick accounting for lick p value
                                                                    ii_t_first=find((t_pac>=handles_pars.analysisWin_times(1)),1,'first');
                                                                    ii_t_last=find(t_pac<=handles_pars.analysisWin_times(2),1,'last');
                                                                    for ww=1:length(trialNos_in_lick_event_Ev)
                                                                        licks_this_trial=zeros(1,ii_t_last-ii_t_first+1);
                                                                        these_lick_times=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_lick_event_Ev(ww)).lickPower_times;
                                                                        for ii_lick=1:length(these_lick_times)
                                                                            if (these_lick_times(ii_lick)>=t_pac(ii_t_first))&(these_lick_times(ii_lick)<t_pac(ii_t_last))
                                                                                this_ii=find((t_pac<=these_lick_times(ii_lick))&(t_pac+(t_pac(2)-t_pac(1))>these_lick_times(ii_lick)),1,'first')-ii_t_first+1;
                                                                                licks_this_trial(this_ii)=1;
                                                                            end
                                                                        end
                                                                        PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).binary_lick=licks_this_trial;
                                                                    end
                                                                    
                                                                    %Save the times for licks
                                                                    for ww=1:sum(trials_in_lick_event)
                                                                        these_lick_times=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_lick_event_Ev(ww)).lickPower_times;
                                                                        these_peak_times=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_lick_event_Ev(ww)).peakPower_times;
                                                                        these_peak_powers=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_lick_event_Ev(ww)).Power_per_peak;
                                                                        these_trough_times=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.PACtimecourse(trialNos_in_lick_event_Ev(ww)).troughPower_times;
                                                                        
                                                                        %Autocorrelogram for licks
                                                                        these_deltas=[];
                                                                        these_d_lick_times=[];
                                                                        xx=0;
                                                                        if length(these_lick_times)>1
                                                                            for zz1=1:length(these_lick_times)
                                                                                if zz1>1
                                                                                    xx=xx+1;
                                                                                    these_deltas(xx)=these_lick_times(zz1-1)-these_lick_times(zz1);
                                                                                    these_d_lick_times(xx)=these_lick_times(zz1);
                                                                                end
                                                                                if zz1<length(these_lick_times)
                                                                                    xx=xx+1;
                                                                                    these_deltas(xx)=these_lick_times(zz1+1)-these_lick_times(zz1);
                                                                                    these_d_lick_times(xx)=these_lick_times(zz1);
                                                                                end
                                                                                
                                                                            end
                                                                        end
                                                                        
                                                                        for ii=1:3
                                                                            PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).lick_auto(ii).dt=...
                                                                                these_deltas((these_d_lick_times>=win_lick_start(ii))&(these_d_lick_times<=win_lick_end(ii)));
                                                                        end
                                                                        
                                                                        %Lick-peak crosscorrelogram and lick-peak delta-gated peak power timecourse
                                                                        if ~isempty(these_lick_times)
                                                                            these_deltas=[];
                                                                            xx=0;
                                                                            for zz1=1:length(these_peak_times)
                                                                                [mint,ii_near]=min(abs(these_peak_times(zz1)-these_lick_times));
                                                                                xx=xx+1;
                                                                                these_deltas(xx)=these_peak_times(zz1)-these_lick_times(ii_near);
                                                                            end
                                                                            
                                                                            for ii=1:3
                                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).lick_peak_cross(ii).dt=...
                                                                                    these_deltas((these_peak_times>=win_lick_start(ii))&(these_peak_times<=win_lick_end(ii)));
                                                                            end
                                                                            
                                                                            min_times=20;
                                                                            
                                                                            %Peak power timecourse below
                                                                            this_peakPower_times=these_peak_times(these_deltas<0);
                                                                            if length(this_peakPower_times)>=min_times
                                                                                this_peakPower=these_peak_powers(these_deltas<0);
                                                                                this_peakPower_timecourse=zeros(1,length(t_pac));
                                                                                for ii_t=1:length(t_pac)
                                                                                    if t_pac(ii_t)<=this_peakPower_times(1)
                                                                                        this_peakPower_timecourse(ii_t)=this_peakPower(1);
                                                                                        %this_peakPower_t_pac(ii_t)=this_peakPower(1);
                                                                                    else
                                                                                        if t_pac(ii_t)>=this_peakPower_times(end)
                                                                                            this_peakPower_timecourse(ii_t)=this_peakPower(end);
                                                                                            %this_peakPower_t_pac(ii_t)=this_peakPower(end);
                                                                                        else
                                                                                            ii_pt=find(this_peakPower_times>=t_pac(ii_t),1,'first');
                                                                                            this_peakPower_timecourse(ii_t)=this_peakPower(ii_pt-1)+(t_pac(ii_t)-this_peakPower_times(ii_pt-1))*((this_peakPower(ii_pt)-this_peakPower(ii_pt-1))/(this_peakPower_times(ii_pt)-this_peakPower_times(ii_pt-1)));
                                                                                        end
                                                                                    end
                                                                                end
                                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_below=this_peakPower_timecourse;
                                                                            else
                                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_below=[];
                                                                            end
                                                                            
                                                                            %Peak power timecourse above
                                                                            this_peakPower_times=these_peak_times(these_deltas>0);
                                                                            if length(this_peakPower_times)>=min_times
                                                                                this_peakPower=these_peak_powers(these_deltas>0);
                                                                                this_peakPower_timecourse=zeros(1,length(t_pac));
                                                                                for ii_t=1:length(t_pac)
                                                                                    if t_pac(ii_t)<=this_peakPower_times(1)
                                                                                        this_peakPower_timecourse(ii_t)=this_peakPower(1);
                                                                                        %this_peakPower_t_pac(ii_t)=this_peakPower(1);
                                                                                    else
                                                                                        if t_pac(ii_t)>=this_peakPower_times(end)
                                                                                            this_peakPower_timecourse(ii_t)=this_peakPower(end);
                                                                                            %this_peakPower_t_pac(ii_t)=this_peakPower(end);
                                                                                        else
                                                                                            ii_pt=find(this_peakPower_times>=t_pac(ii_t),1,'first');
                                                                                            this_peakPower_timecourse(ii_t)=this_peakPower(ii_pt-1)+(t_pac(ii_t)-this_peakPower_times(ii_pt-1))*((this_peakPower(ii_pt)-this_peakPower(ii_pt-1))/(this_peakPower_times(ii_pt)-this_peakPower_times(ii_pt-1)));
                                                                                        end
                                                                                    end
                                                                                end
                                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_above=this_peakPower_timecourse;
                                                                            else
                                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_above=[];
                                                                            end
                                                                            
                                                                            %Peak power timecourse all
                                                                            this_peakPower_times=these_peak_times;
                                                                            if length(this_peakPower_times)>=min_times
                                                                                this_peakPower=these_peak_powers;
                                                                                this_peakPower_timecourse=zeros(1,length(t_pac));
                                                                                for ii_t=1:length(t_pac)
                                                                                    if t_pac(ii_t)<=this_peakPower_times(1)
                                                                                        this_peakPower_timecourse(ii_t)=this_peakPower(1);
                                                                                        %this_peakPower_t_pac(ii_t)=this_peakPower(1);
                                                                                    else
                                                                                        if t_pac(ii_t)>=this_peakPower_times(end)
                                                                                            this_peakPower_timecourse(ii_t)=this_peakPower(end);
                                                                                            %this_peakPower_t_pac(ii_t)=this_peakPower(end);
                                                                                        else
                                                                                            ii_pt=find(this_peakPower_times>=t_pac(ii_t),1,'first');
                                                                                            this_peakPower_timecourse(ii_t)=this_peakPower(ii_pt-1)+(t_pac(ii_t)-this_peakPower_times(ii_pt-1))*((this_peakPower(ii_pt)-this_peakPower(ii_pt-1))/(this_peakPower_times(ii_pt)-this_peakPower_times(ii_pt-1)));
                                                                                        end
                                                                                    end
                                                                                end
                                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_all=this_peakPower_timecourse;
                                                                            else
                                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_all=[];
                                                                            end
                                                                        end
                                                                        
                                                                        %Lick-trough crosscorrelogram
                                                                        if ~isempty(these_lick_times)
                                                                            these_deltas=[];
                                                                            xx=0;
                                                                            for zz1=1:length(these_trough_times)
                                                                                [mint,ii_near]=min(abs(these_trough_times(zz1)-these_lick_times));
                                                                                xx=xx+1;
                                                                                these_deltas(xx)=these_trough_times(zz1)-these_lick_times(ii_near);
                                                                            end
                                                                            
                                                                            for ii=1:3
                                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).lick_trough_cross(ii).dt=...
                                                                                    these_deltas((these_trough_times>=win_lick_start(ii))&(these_trough_times<=win_lick_end(ii)));
                                                                            end
                                                                            
                                                                        end
                                                                    end
                                                                    
                                                                    %Get z scored peak power below
                                                                    if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power,'peakPower_timecourse_below')
                                                                        
                                                                        these_ref_power=[];
                                                                        for ww=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power)
                                                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_below)
                                                                                this_p=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_below;
                                                                                these_ref_power=[these_ref_power this_p((t_pac>=handles_pars.refWin_times(1))&(t_pac<=handles_pars.refWin_times(2)))];
                                                                            end
                                                                        end
                                                                        ref_power_SD=std(these_ref_power);
                                                                        ref_power_mean=mean(these_ref_power);
                                                                        
                                                                        
                                                                        for ww=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power)
                                                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_below)
                                                                                this_p=(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_below-ref_power_mean)/ref_power_SD;
                                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).z_peakPower_timecourse_below=this_p((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
                                                                            end
                                                                        end
                                                                    end
                                                                    
                                                                    %Get z scored peak power above
                                                                    if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power,'peakPower_timecourse_above')
                                                                        
                                                                        these_ref_power=[];
                                                                        for ww=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power)
                                                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_above)
                                                                                this_p=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_above;
                                                                                these_ref_power=[these_ref_power this_p((t_pac>=handles_pars.refWin_times(1))&(t_pac<=handles_pars.refWin_times(2)))];
                                                                            end
                                                                        end
                                                                        ref_power_SD=std(these_ref_power);
                                                                        ref_power_mean=mean(these_ref_power);
                                                                        
                                                                        
                                                                        for ww=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power)
                                                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_above)
                                                                                this_p=(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_above-ref_power_mean)/ref_power_SD;
                                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).z_peakPower_timecourse_above=this_p((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
                                                                            end
                                                                        end
                                                                    end
                                                                    
                                                                    %Get z scored peak power all
                                                                    if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power,'peakPower_timecourse_all')
                                                                        
                                                                        these_ref_power=[];
                                                                        for ww=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power)
                                                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_all)
                                                                                this_p=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_all;
                                                                                these_ref_power=[these_ref_power this_p((t_pac>=handles_pars.refWin_times(1))&(t_pac<=handles_pars.refWin_times(2)))];
                                                                            end
                                                                        end
                                                                        ref_power_SD=std(these_ref_power);
                                                                        ref_power_mean=mean(these_ref_power);
                                                                        
                                                                        
                                                                        for ww=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power)
                                                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_all)
                                                                                this_p=(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).peakPower_timecourse_all-ref_power_mean)/ref_power_SD;
                                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).z_power(ww).z_peakPower_timecourse_all=this_p((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
                                                                            end
                                                                        end
                                                                    end
                                                                else
                                                                    PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).no_z_lick_power=0;
                                                                end
                                                            else
                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elec).event(evNo).no_z_lick_power=0;
                                                            end
                                                            
                                                            %Save lick analysis
                                                            if (pacii==1)&(which_electrodes(1)==elec)
                                                                %Do lick referenced power
                                                                lick_trials_included=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).lickStuff.lick_trials_included;
                                                                these_lick_trials=zeros(1,length(lick_trials_included));
                                                                for trNo=1:length(trials_in_event_Ev)
                                                                    if trials_in_event_Ev(trNo)==1
                                                                        this_lick_tr=find(lick_trials_included==trNo);
                                                                        if ~isempty(this_lick_tr)
                                                                            these_lick_trials(1,this_lick_tr)=1;
                                                                        end
                                                                    end
                                                                end
                                                                
                                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).no_lick_events=sum(these_lick_trials);
                                                                last_ii=0;
                                                                for ww=1:sum(these_lick_trials)
                                                                    this_delta_ii=find(these_lick_trials(last_ii+1:end)==1,1,'first');
                                                                    PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).lick_tr(ww).binary_lick_per_t=...
                                                                        handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).lickStuff.binary_lick_per_t(last_ii+this_delta_ii,(t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
                                                                    last_ii=last_ii+this_delta_ii;
                                                                end
                                                            end
                                                        end
                                                        
                                                        
                                                        fprintf(1, ['%d trials in event No %d succesfully processed for file No %d electrode %d\n'],sum(trials_in_event_Ev), evNo,fileNo,elec);
                                                        
                                                    else
                                                        
                                                        
                                                        fprintf(1, ['%d trials in event No %d fewer than minimum trials per event ' evTypeLabels{evNo} ' for file No %d electrode %d\n'],sum(trials_in_event_Ev), min_trials_per_event,fileNo,elec);
                                                        
                                                        
                                                    end
                                                    
                                                    
                                                end
                                                
                                            else
                                                PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed=0;
                                                fprintf(1, ['Empty percent_mask for file No %d electrode %d\n'],fileNo,elec);
                                                
                                            end
                                            
                                        else
                                            
                                            fprintf(1, ['Empty PAC for file No %d electrode %d\n'],fileNo,elec);
                                            
                                        end
                                        
                                        
                                    else
                                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],fileNo,elec);
                                        
                                        
                                    end
                                end %if mouseNo
                            end
                        end %fileNo
                        
                    end
                    
                end
                
                if mouse_has_files==1
                    mouse_included(mouseNo)=1;
                else
                    mouse_included(mouseNo)=0;
                end
            end
            
            
        end
        fprintf(1, '\n\n')
        
        elec=[];
        
        figureNo = 0;
        anal_t_pac=t_pac((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
        
        %Let's try a per mouse z score t test (using the z scores for all
        %electrodes?)
        handles_out.correl=[];
        handles_out.RP_ii=0;
        handles_out.anal_t_pac=anal_t_pac;
        
        
        
        for pacii=1:no_pacii
            
                for per_ii=2:-1:1
                    
                    %For the lick autocorrelogram we only use the first electrode and first
                    %pacii
                    if pacii==1
                        elecNo=which_electrodes(1);
                        for mouseNo=1:length(PRPtimecourse.mouse)
                            
                            for ii_dt=1:length(win_lick_start)
                                handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_auto(ii_dt).dt=zeros(1,length(edges)-1);
                                handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_auto(ii_dt).dt=zeros(1,length(edges)-1);
                            end
                            handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_auto_calc=0;
                            handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_auto_calc=0;
                            
                            for this_session=1:PRPtimecourse.mouse(mouseNo).no_sessions
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed==1
                                    
                                    
                                    
                                    %Splus
                                    evNo=1;
                                    this_calc=0;
                                    if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                        for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                            if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'lick_auto')
                                                if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_auto)
                                                    for ii_dt=1:length(win_lick_start)
                                                        if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_auto(ii_dt).dt)
                                                            this_lick_auto=[];
                                                            this_lick_auto=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_auto(ii_dt).dt;
                                                            counts=histcounts(this_lick_auto,edges);
                                                            handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_auto(ii_dt).dt(:,:)=...
                                                                handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_auto(ii_dt).dt+counts;
                                                            this_calc=1;
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_auto_calc=this_calc;
                                    
                                    %Sminus
                                    evNo=2;
                                    this_calc=0;
                                    if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                        for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                            if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'lick_auto')
                                                if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_auto)
                                                    for ii_dt=1:length(win_lick_start)
                                                        if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_auto(ii_dt).dt)
                                                            this_lick_auto=[];
                                                            this_lick_auto=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_auto(ii_dt).dt;
                                                            counts=histcounts(this_lick_auto,edges);
                                                            handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_auto(ii_dt).dt(:,:)=...
                                                                handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_auto(ii_dt).dt+counts;
                                                            this_calc=1;
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_auto_calc=this_calc;
                                    
                                end
                            end
                        end
                        
                    end
                    
                    %Now do the peaks
                    
                    for mouseNo=1:length(PRPtimecourse.mouse)
                        for ii_dt=1:length(win_lick_start)
                            handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_peak_cross(ii_dt).dt=zeros(1,length(edges)-1);
                            handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_peak_cross(ii_dt).dt=zeros(1,length(edges)-1);
                        end
                        handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_peak_cross_calc=0;
                        handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_peak_cross_calc=0;
                        
                        for elecNo=which_electrodes
                            for this_session=1:PRPtimecourse.mouse(mouseNo).no_sessions
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed==1
                                    
                                    
                                    
                                    %Splus
                                    evNo=1;
                                    this_calc=0;
                                    if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                        for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                            if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'lick_peak_cross')
                                                if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_peak_cross)
                                                    for ii_dt=1:length(win_lick_start)
                                                        if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_peak_cross(ii_dt).dt)
                                                            this_lick_peak_cross=[];
                                                            this_lick_peak_cross=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_peak_cross(ii_dt).dt;
                                                            counts=histcounts(this_lick_peak_cross,edges);
                                                            handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_peak_cross(ii_dt).dt(:,:)=...
                                                                handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_peak_cross(ii_dt).dt+counts;
                                                            this_calc=1;
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_peak_cross_calc=this_calc;
                                    
                                    %Sminus
                                    evNo=2;
                                    this_calc=0;
                                    if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                        for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                            if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'lick_peak_cross')
                                                if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_peak_cross)
                                                    for ii_dt=1:length(win_lick_start)
                                                        if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_peak_cross(ii_dt).dt)
                                                            this_lick_peak_cross=[];
                                                            this_lick_peak_cross=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_peak_cross(ii_dt).dt;
                                                            counts=histcounts(this_lick_peak_cross,edges);
                                                            handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_peak_cross(ii_dt).dt(:,:)=...
                                                                handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_peak_cross(ii_dt).dt+counts;
                                                            this_calc=1;
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_peak_cross_calc=this_calc;
                                    
                                end
                            end
                        end
                    end
                    
                    %Now do the troughs
                    
                    for mouseNo=1:length(PRPtimecourse.mouse)
                         for ii_dt=1:length(win_lick_start)
                            handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_trough_cross(ii_dt).dt=zeros(1,length(edges)-1);
                            handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_trough_cross(ii_dt).dt=zeros(1,length(edges)-1);
                         end
                         
                        handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_trough_cross_calc=0;
                        handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_trough_cross_calc=0;
                        for elecNo=which_electrodes
                            for this_session=1:PRPtimecourse.mouse(mouseNo).no_sessions
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed==1
                                    
                                    
                                    
                                    %Splus
                                    evNo=1;
                                    this_calc=0;
                                    if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                        for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                            if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'lick_trough_cross')
                                                if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_trough_cross)
                                                    for ii_dt=1:length(win_lick_start)
                                                        if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_trough_cross(ii_dt).dt)
                                                            this_lick_trough_cross=[];
                                                            this_lick_trough_cross=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_trough_cross(ii_dt).dt;
                                                            counts=histcounts(this_lick_trough_cross,edges);
                                                            handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_trough_cross(ii_dt).dt(:,:)=...
                                                                handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_trough_cross(ii_dt).dt+counts;
                                                            this_calc=1;
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_trough_cross_calc=this_calc;
                                    
                                    %Sminus
                                    evNo=2;
                                    this_calc=0;
                                    if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                        for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                            if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'lick_trough_cross')
                                                if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_trough_cross)
                                                    for ii_dt=1:length(win_lick_start)
                                                        if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_trough_cross(ii_dt).dt)
                                                            this_lick_trough_cross=[];
                                                            this_lick_trough_cross=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lick_trough_cross(ii_dt).dt;
                                                            counts=histcounts(this_lick_trough_cross,edges);
                                                            handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_trough_cross(ii_dt).dt(:,:)=...
                                                                handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_trough_cross(ii_dt).dt+counts;
                                                            this_calc=1;
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_trough_cross_calc=this_calc;
                                    
                                end
                            end
                        end
                    end
     
                    pffft=1;
                end
            
        end
        
        %Now do lick frequency timecourse
        handles_out.lickf_timecourse=[];
        handles_out.lf_ii=0;
        
        pacii=1;
        
        for per_ii=2:-1:1
            
            %For the lick autocorrelogram we only use the first electrode and first
            %pacii
            
            elecNo=which_electrodes(1);
            for mouseNo=1:length(PRPtimecourse.mouse)
                
                handles_out.lf_ii=handles_out.lf_ii+1;
                handles_out.lickf_timecourse(handles_out.lf_ii).pacii=pacii;
                handles_out.lickf_timecourse(handles_out.lf_ii).per_ii=per_ii;
                handles_out.lickf_timecourse(handles_out.lf_ii).group_no=PRPtimecourse.mouse(mouseNo).group_no;
                handles_out.lickf_timecourse(handles_out.lf_ii).mouseNo=mouseNo;
                
                this_Sp_binary_lick=[];
                ii_sp=0;
                this_Sm_binary_lick=[];
                ii_sm=0;
                
                for this_session=1:PRPtimecourse.mouse(mouseNo).no_sessions
                    if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed==1
                        
                        
                        %Splus
                        evNo=1;
                        if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                            for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'binary_lick')
                                    if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).binary_lick)
                                        this_binary_lick=[];
                                        this_binary_lick=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).binary_lick;
                                        ii_sp=ii_sp+1;
                                        this_Sp_binary_lick(ii_sp,:)=this_binary_lick;
                                    end
                                end
                            end
                        end
                        
                        
                        %Sminus
                        evNo=2;
                        this_calc=0;
                        if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                            for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'binary_lick')
                                    if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).binary_lick)
                                        this_binary_lick=[];
                                        this_binary_lick=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).binary_lick;
                                        ii_sm=ii_sm+1;
                                        this_Sm_binary_lick(ii_sm,:)=this_binary_lick;
                                    end
                                end
                            end
                        end
                        
                    end
                end
                
                handles_out.lickf_timecourse(handles_out.lf_ii).Sp_binary_lick=this_Sp_binary_lick;
                handles_out.lickf_timecourse(handles_out.lf_ii).Sp_lick_freq=sum(this_Sp_binary_lick,1)/(size(this_Sp_binary_lick,1)*(anal_t_pac(2)-anal_t_pac(1)));
                handles_out.lickf_timecourse(handles_out.lf_ii).Sm_binary_lick=this_Sm_binary_lick;
                handles_out.lickf_timecourse(handles_out.lf_ii).Sm_lick_freq=sum(this_Sm_binary_lick,1)/(size(this_Sm_binary_lick,1)*(anal_t_pac(2)-anal_t_pac(1)));
            end
            
        end
        
        
        
        dt_legend{1}='pre';
        dt_legend{2}='decision';
        dt_legend{3}='post';
        
        figureNo=0;
        
        between_edges=(edges(2:end)+edges(1:end-1))/2;
        
        %Plot the histograms for peak and trough relative to licks
        for groupNo=1:3
            for pacii=1:no_pacii
                for per_ii=2:-1:1
                    figureNo=figureNo+1;
                    try
                        close(figureNo)
                    catch
                    end
                    hFig=figure(figureNo);
                    
                    set(hFig, 'units','normalized','position',[.1 .1 .7 .7])
                    %             subplot(2,1,1)
                    ax=gca;ax.LineWidth=3;
                    
                    hold on
                    
                   %First do lick auto
                    for ii_dt=1:3
                        this_Sp_lick_auto=[];
                        this_Sm_lick_auto=[];
                        this_Sp_lick_auto_asym=[];
                        this_Sm_lick_auto_asym=[];
                        ii_mouse=0;
                        for mouseNo=1:length(PRPtimecourse.mouse)
                             
                            group_no=PRPtimecourse.mouse(mouseNo).group_no;
                            
                            if group_no==groupNo
                               
                                if (isfield(handles_out.correl.mouseNo(mouseNo).pacii(1).per_ii(per_ii),'Sp_lick_auto'))&(isfield(handles_out.correl.mouseNo(mouseNo).pacii(1).per_ii(per_ii),'Sm_lick_auto'))
                                    this_Sp=handles_out.correl.mouseNo(mouseNo).pacii(1).per_ii(per_ii).Sp_lick_auto(ii_dt).dt;
                                    this_Sm=handles_out.correl.mouseNo(mouseNo).pacii(1).per_ii(per_ii).Sm_lick_auto(ii_dt).dt;
                                    if (sum(this_Sp)~=0)&(sum(this_Sm)~=0)
                                        ii_mouse=ii_mouse+1;
                                        this_Sp=this_Sp/sum(this_Sp);
                                        this_Sm=this_Sm/sum(this_Sm);
                                        this_Sp_lick_auto(ii_mouse,:)=this_Sp;
                                        this_Sm_lick_auto(ii_mouse,:)=this_Sm;
                                        this_Sp_lick_auto_asym(ii_mouse)=((sum(this_Sp(between_edges>0))/sum(this_Sp))-0.5)/0.5;
                                        this_Sm_lick_auto_asym(ii_mouse)=((sum(this_Sm(between_edges>0))/sum(this_Sm))-0.5)/0.5;
                                    end
                                end
                            end
                        end
                        subplot(3,3,ii_dt)
                        hold on
                        
                        handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).pSp_lick_auto=this_Sp_lick_auto;
                        handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).pSm_lick_auto=this_Sm_lick_auto;
                        handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).ii_mouse_lick=ii_mouse;
                        handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).Sp_lick_auto_asym=this_Sp_lick_auto_asym;
                        handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).Sm_lick_auto_asym=this_Sm_lick_auto_asym;
                        
                        try
                            CIsm = bootci(1000, @mean, this_Sm_lick_auto);
                            meansm=mean(this_Sm_lick_auto,1);
                            CIsm(1,:)=meansm-CIsm(1,:);
                            CIsm(2,:)=CIsm(2,:)-meansm;
                            
                            [hlsm, hpsm] = boundedline(between_edges',mean(this_Sm_lick_auto,1)', CIsm', 'b');
                            
                            CIsp = bootci(1000, @mean, this_Sp_lick_auto);
                            meansp=mean(this_Sp_lick_auto,1);
                            CIsp(1,:)=meansp-CIsp(1,:);
                            CIsp(2,:)=CIsp(2,:)-meansp;
                            
                            
                            [hlsp, hpsp] = boundedline(between_edges',mean(this_Sp_lick_auto,1)', CIsp', 'r');
                        catch
                        end
                        
                        plot([0 0],[0 0.35],'-k')
                        
                        title(['licks ' dt_legend{ii_dt}])
                        xlabel('Time(sec)')
                        ylabel('p')
                        xlim([-0.3 0.3])
                        ylim([0 0.05])
                        
                        
                    end
                    
                    
                    %Now do peak cross
                    for ii_dt=1:3
                        this_Sp_lick_peak_cross=[];
                        this_Sm_lick_peak_cross=[];
                        this_Sp_lick_peak_cross_asym=[];
                        this_Sm_lick_peak_cross_asym=[];
                        ii_mouse=0;
                        for mouseNo=1:length(PRPtimecourse.mouse)
                            
                            group_no=PRPtimecourse.mouse(mouseNo).group_no;
                            
                            if group_no==groupNo
                                
                                if (isfield(handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii),'Sp_lick_peak_cross'))...
                                        &(isfield(handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii),'Sm_lick_peak_cross')) 
                                    this_Sp=handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_peak_cross(ii_dt).dt;
                                    this_Sm=handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_peak_cross(ii_dt).dt;
                                    if (sum(this_Sp)~=0)&(sum(this_Sm)~=0)
                                        ii_mouse=ii_mouse+1;
                                        this_Sp=this_Sp/sum(this_Sp);
                                        this_Sm=this_Sm/sum(this_Sm);
                                        this_Sp_lick_peak_cross(ii_mouse,:)=this_Sp;
                                        this_Sm_lick_peak_cross(ii_mouse,:)=this_Sm;
                                        this_Sp_lick_peak_cross_asym(ii_mouse)=((sum(this_Sp(between_edges>0))/sum(this_Sp))-0.5)/0.5;
                                        this_Sm_lick_peak_cross_asym(ii_mouse)=((sum(this_Sm(between_edges>0))/sum(this_Sm))-0.5)/0.5;
                                    end
                                end
                            end
                        end
                        
                        handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).pSp_lick_peak_cross=this_Sp_lick_peak_cross;
                        handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).pSm_lick_peak_cross=this_Sm_lick_peak_cross;
                        handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).ii_mouse_peak=ii_mouse;
                        handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).Sp_lick_peak_cross_asym=this_Sp_lick_peak_cross_asym;
                        handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).Sm_lick_peak_cross_asym=this_Sm_lick_peak_cross_asym;
                        
                        
                        subplot(3,3,3+ii_dt)
                        hold on
                        
                        try
                            CIsm = bootci(1000, @mean, this_Sm_lick_peak_cross);
                            meansm=mean(this_Sm_lick_peak_cross,1);
                            CIsm(1,:)=meansm-CIsm(1,:);
                            CIsm(2,:)=CIsm(2,:)-meansm;
                            
                            [hlsm, hpsm] = boundedline(between_edges',mean(this_Sm_lick_peak_cross,1)', CIsm', 'b');
                            
                            CIsp = bootci(1000, @mean, this_Sp_lick_peak_cross);
                            meansp=mean(this_Sp_lick_peak_cross,1);
                            CIsp(1,:)=meansp-CIsp(1,:);
                            CIsp(2,:)=CIsp(2,:)-meansp;
                            
                            
                            [hlsp, hpsp] = boundedline(between_edges',mean(this_Sp_lick_peak_cross,1)', CIsp', 'r');
                        catch
                        end
                        
                        plot([0 0],[0 0.05],'-k')
                        
                        title(['peak ' dt_legend{ii_dt}])
                        xlabel('Time(sec)')
                        ylabel('p')
                        xlim([-0.3 0.3])
                        ylim([0 0.05])
                        
                    end
                           
                    %Now do trough cross
                    for ii_dt=1:3
                        this_Sp_lick_trough_cross=[];
                        this_Sm_lick_trough_cross=[];
                        this_Sp_lick_trough_cross_asym=[];
                        this_Sm_lick_trough_cross_asym=[];
                        ii_mouse=0;
                        for mouseNo=1:length(PRPtimecourse.mouse)
                            
                            group_no=PRPtimecourse.mouse(mouseNo).group_no;
                            
                            if group_no==groupNo
                                
                                if (isfield(handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii),'Sp_lick_peak_cross'))&(isfield(handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii),'Sm_lick_peak_cross'))
                                    this_Sp=handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sp_lick_trough_cross(ii_dt).dt;
                                    this_Sm=handles_out.correl.mouseNo(mouseNo).pacii(pacii).per_ii(per_ii).Sm_lick_trough_cross(ii_dt).dt;
                                    if (sum(this_Sp)~=0)&(sum(this_Sm)~=0)
                                        ii_mouse=ii_mouse+1;
                                        this_Sp=this_Sp/sum(this_Sp);
                                        this_Sm=this_Sm/sum(this_Sm);
                                        this_Sp_lick_trough_cross(ii_mouse,:)=this_Sp;
                                        this_Sm_lick_trough_cross(ii_mouse,:)=this_Sm;
                                        this_Sp_lick_trough_cross_asym(ii_mouse)=((sum(this_Sp(between_edges>0))/sum(this_Sp))-0.5)/0.5;
                                        this_Sm_lick_trough_cross_asym(ii_mouse)=((sum(this_Sm(between_edges>0))/sum(this_Sm))-0.5)/0.5;
                                    end
                                end
                            end
                        end
                        
                         handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).pSp_lick_peak_cross=this_Sp_lick_trough_cross;
                        handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).pSm_lick_peak_cross=this_Sm_lick_trough_cross;
                        handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).ii_mouse_trough=ii_mouse;
                        handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).Sp_lick_trough_cross_asym=this_Sp_lick_trough_cross_asym;
                        handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).Sm_lick_trough_cross_asym=this_Sm_lick_trough_cross_asym;
                        
                        subplot(3,3,6+ii_dt)
                        hold on 
                        try
                            CIsm = bootci(1000, @mean, this_Sm_lick_trough_cross);
                            meansm=mean(this_Sm_lick_trough_cross,1);
                            CIsm(1,:)=meansm-CIsm(1,:);
                            CIsm(2,:)=CIsm(2,:)-meansm;
                            
                            
                            [hlsm, hpsm] = boundedline(between_edges',mean(this_Sm_lick_trough_cross,1)', CIsm', 'b');
                        catch
                        end
                        
                        try
                            CIsp = bootci(1000, @mean, this_Sp_lick_trough_cross);
                            meansp=mean(this_Sp_lick_trough_cross,1);
                            CIsp(1,:)=meansp-CIsp(1,:);
                            CIsp(2,:)=CIsp(2,:)-meansp;
                            
                            
                            [hlsp, hpsp] = boundedline(between_edges',mean(this_Sp_lick_trough_cross,1)', CIsp', 'r');
                        catch
                        end
                        
                        plot([0 0],[0 0.35],'-k')
                        
                        title(['trough ' dt_legend{ii_dt}])
                        xlabel('Time(sec)')
                        ylabel('p')
                        xlim([-0.3 0.3])
                        ylim([0 0.35])
                        
                    end
                    
                    sgtitle(['Cross ' freq_names{pacii+1} ' ' prof_naive_leg{per_ii} ' ' group_legend{groupNo} ])
                    
                    pfft=1;
                end
            end
        end
        pffft=1;
        
        these_edges=-1:0.1:1;
        rand_offset=0.8;
        %Bar graph for asymmetry
        for pacii=1:3
            figureNo=figureNo+1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            set(hFig, 'units','normalized','position',[.1 .1 .7 .7])
            
            ax=gca;ax.LineWidth=3;
            
            
            %Proficient Splus
           
            
            bar_offset=0;
            
  
            for per_ii=2:-1:1
                subplot(2,1,per_ii)
                hold on
                for grNo=1:3
                    
                    for ii_dt=1:3
                        these_asym=handles_out.p_correl.groupNo(group_no).pacii(pacii).per_ii(per_ii).dt(ii_dt).Sp_lick_peak_cross_asym;
                           switch grNo
                                case 1
                                    bar(bar_offset,mean(these_asym),'g','LineWidth', 3,'EdgeColor','none')
                                case 2
                                    bar(bar_offset,mean(these_asym),'b','LineWidth', 3,'EdgeColor','none')
                                case 3
                                    bar(bar_offset,mean(these_asym),'y','LineWidth', 3,'EdgeColor','none')
                           end
                             %Violin plot
                            [mean_out, CIout]=drgViolinPoint(these_asym,these_edges,bar_offset,rand_offset,'k','k',3);
                            bar_offset=bar_offset+1;
                    end
                    bar_offset=bar_offset+1;
                end
                title(prof_naive_leg{per_ii})
            end
            sgtitle(['Asymmetry ' freq_names{pacii+1} ' S+'])
                    
        end
        
        %Now plot the lick frequency timecourse
        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        
        set(hFig, 'units','normalized','position',[.1 .1 .7 .7])
        
        ax=gca;ax.LineWidth=3;
        
        for per_ii=2:-1:1
            for groupNo=1:3
                
                
                
                Sp_lick_freq=[];
                Sm_lick_freq=[];
                ii_mouse=0;
                for lf_ii=1:handles_out.lf_ii
                    if (handles_out.lickf_timecourse(lf_ii).group_no==groupNo)&...
                            (handles_out.lickf_timecourse(lf_ii).per_ii==per_ii)
                        
                        if (isfield(handles_out.lickf_timecourse(lf_ii),'Sp_lick_freq'))&(isfield(handles_out.lickf_timecourse(lf_ii),'Sm_lick_freq'))
                            if (~isempty(handles_out.lickf_timecourse(lf_ii).Sp_lick_freq))&(~isempty(handles_out.lickf_timecourse(lf_ii).Sm_lick_freq))
                                this_Sp_lick_freq=handles_out.lickf_timecourse(lf_ii).Sp_lick_freq;
                                this_Sm_lick_freq=handles_out.lickf_timecourse(lf_ii).Sm_lick_freq;
                                
                                ii_mouse=ii_mouse+1;
                                Sp_lick_freq(ii_mouse,:)=this_Sp_lick_freq;
                                Sm_lick_freq(ii_mouse,:)=this_Sm_lick_freq;
                            end
                        end
                    end
                end
                
                subplot(2,3,3*(per_ii-1)+groupNo)
                hold on
                
                CIsm = bootci(1000, @mean, Sm_lick_freq);
                meansm=mean(Sm_lick_freq,1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;
                
                
                [hlsm, hpsm] = boundedline(anal_t_pac',mean(Sm_lick_freq,1)', CIsm', 'b');
                
                CIsp = bootci(1000, @mean, Sm_lick_freq);
                meansp=mean(Sm_lick_freq,1);
                CIsp(1,:)=meansp-CIsp(1,:);
                CIsp(2,:)=CIsp(2,:)-meansp;
                
                
                [hlsp, hpsp] = boundedline(anal_t_pac',mean(Sp_lick_freq,1)', CIsp', 'r');
                
                
                plot([0 0],[0 0.35],'-k')
                
                title([prof_naive_leg{per_ii} ' ' group_legend{groupNo} ])
                xlabel('Time(sec)')
                ylabel('licks/sec')
               
            end
        end
        
        sgtitle('Lick frequency')
        
        %For WT plot the lick frequency timecourse for naive
        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        
        set(hFig, 'units','normalized','position',[.3 .3 .4 .25])
        
        ax=gca;ax.LineWidth=3;
        hold on
        
        per_ii=2;
        groupNo=1;
        
        
        
        Sp_lick_freq=[];
        Sm_lick_freq=[];
        ii_mouse=0;
        for lf_ii=1:handles_out.lf_ii
            if (handles_out.lickf_timecourse(lf_ii).group_no==groupNo)&...
                    (handles_out.lickf_timecourse(lf_ii).per_ii==per_ii)
                
                if (isfield(handles_out.lickf_timecourse(lf_ii),'Sp_lick_freq'))&(isfield(handles_out.lickf_timecourse(lf_ii),'Sm_lick_freq'))
                    if (~isempty(handles_out.lickf_timecourse(lf_ii).Sp_lick_freq))&(~isempty(handles_out.lickf_timecourse(lf_ii).Sm_lick_freq))
                        this_Sp_lick_freq=handles_out.lickf_timecourse(lf_ii).Sp_lick_freq;
                        this_Sm_lick_freq=handles_out.lickf_timecourse(lf_ii).Sm_lick_freq;
                        
                        ii_mouse=ii_mouse+1;
                        Sp_lick_freq(ii_mouse,:)=this_Sp_lick_freq;
                        Sm_lick_freq(ii_mouse,:)=this_Sm_lick_freq;
                    end
                end
            end
        end
        
 
        
        CIsm = bootci(1000, @mean, Sm_lick_freq);
        meansm=mean(Sm_lick_freq,1);
        CIsm(1,:)=meansm-CIsm(1,:);
        CIsm(2,:)=CIsm(2,:)-meansm;
        
        
        [hlsm, hpsm] = boundedline(anal_t_pac',mean(Sm_lick_freq,1)', CIsm', 'cmap',[80/255 194/255 255/255]);
        
        CIsp = bootci(1000, @mean, Sm_lick_freq);
        meansp=mean(Sm_lick_freq,1);
        CIsp(1,:)=meansp-CIsp(1,:);
        CIsp(2,:)=CIsp(2,:)-meansp;
        
        
        [hlsp, hpsp] = boundedline(anal_t_pac',mean(Sp_lick_freq,1)', CIsp', 'cmap',[238/255 111/255 179/255]);
        
        plot(anal_t_pac',mean(Sm_lick_freq,1)','-','Color',[80/255 194/255 255/255]);
        plot(anal_t_pac',mean(Sp_lick_freq,1)','-','Color',[238/255 111/255 179/255]);
        
        ylim([0 10])
        xlim([-1.45 4.5])
        
        title('Lick rate for wild type naive')
        xlabel('Time(sec)')
        ylabel('licks/sec')
        
        handles_out.lick_rate.Sp_lick_freq_n=Sp_lick_freq;
        handles_out.lick_rate.mean_Sp_lick_freq_n=mean(Sp_lick_freq,1)';
        handles_out.lick_rate.Sm_lick_freq_n=Sm_lick_freq;
        handles_out.lick_rate.mean_Sm_lick_freq_n=mean(Sm_lick_freq,1)';
        handles_out.lick_rate.anal_t_pac=anal_t_pac';
        %For WT plot the lick frequency timecourse for proficient
        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        
        set(hFig, 'units','normalized','position',[.3 .3 .4 .25])
        
        ax=gca;ax.LineWidth=3;
        hold on
        
        per_ii=1;
        groupNo=1;
        
        
        
        Sp_lick_freq=[];
        Sm_lick_freq=[];
        ii_mouse=0;
        for lf_ii=1:handles_out.lf_ii
            if (handles_out.lickf_timecourse(lf_ii).group_no==groupNo)&...
                    (handles_out.lickf_timecourse(lf_ii).per_ii==per_ii)
                
                if (isfield(handles_out.lickf_timecourse(lf_ii),'Sp_lick_freq'))&(isfield(handles_out.lickf_timecourse(lf_ii),'Sm_lick_freq'))
                    if (~isempty(handles_out.lickf_timecourse(lf_ii).Sp_lick_freq))&(~isempty(handles_out.lickf_timecourse(lf_ii).Sm_lick_freq))
                        this_Sp_lick_freq=handles_out.lickf_timecourse(lf_ii).Sp_lick_freq;
                        this_Sm_lick_freq=handles_out.lickf_timecourse(lf_ii).Sm_lick_freq;
                        
                        ii_mouse=ii_mouse+1;
                        Sp_lick_freq(ii_mouse,:)=this_Sp_lick_freq;
                        Sm_lick_freq(ii_mouse,:)=this_Sm_lick_freq;
                    end
                end
            end
        end
        
 
        
        CIsm = bootci(1000, @mean, Sm_lick_freq);
        meansm=mean(Sm_lick_freq,1);
        CIsm(1,:)=meansm-CIsm(1,:);
        CIsm(2,:)=CIsm(2,:)-meansm;
        
        
        [hlsm, hpsm] = boundedline(anal_t_pac',mean(Sm_lick_freq,1)', CIsm', 'cmap',[0 114/255 178/255]);
        
        CIsp = bootci(1000, @mean, Sm_lick_freq);
        meansp=mean(Sm_lick_freq,1);
        CIsp(1,:)=meansp-CIsp(1,:);
        CIsp(2,:)=CIsp(2,:)-meansp;
        
        
        [hlsp, hpsp] = boundedline(anal_t_pac',mean(Sp_lick_freq,1)', CIsp', 'cmap',[158/255 31/255 99/255]);
        
        plot(anal_t_pac',mean(Sm_lick_freq,1)','-','Color',[0 114/255 178/255]);
        plot(anal_t_pac',mean(Sp_lick_freq,1)','-','Color',[158/255 31/255 99/255]);
        
        ylim([0 10])
        xlim([-1.45 4.5])
        
        title('Lick rate for wild type proficient')
        xlabel('Time(sec)')
        ylabel('licks/sec')
        
        handles_out.lick_rate.Sp_lick_freq_p=Sp_lick_freq;
        handles_out.lick_rate.mean_Sp_lick_freq_p=mean(Sp_lick_freq,1)';
        handles_out.lick_rate.Sm_lick_freq_p=Sm_lick_freq;
        handles_out.lick_rate.mean_Sm_lick_freq_p=mean(Sm_lick_freq,1)';
        
        %
        %          anal_t_pac=t_pac((t_pac>=handles_pars.analysisWin_times(1))&(t_pac<=handles_pars.analysisWin_times(2)));
        %Let's try a per mouse z score t test (using the z scores for all
        %electrodes?)
        handles_out.RPtimecourse=[];
        handles_out.RP_ii=0;
        handles_out.anal_t_pac=anal_t_pac;
        
        for pacii=1:no_pacii
            for group_no=1:3
                handles_out.disc_t.pacii(pacii).groupNo(group_no).no_peak=0;
                handles_out.disc_t.pacii(pacii).groupNo(group_no).no_licks=0;
                handles_out.disc_t.pacii(pacii).groupNo(group_no).no_trough=0;
            end
        end
        
        figureNo = figureNo + 1;
        for pacii=1:no_pacii
            for per_ii=2:-1:1
                
                for mouseNo=1:length(PRPtimecourse.mouse)
                    handles_out.RP_ii=handles_out.RP_ii+1;
                    handles_out.RPtimecourse(handles_out.RP_ii).pacii=pacii;
                    handles_out.RPtimecourse(handles_out.RP_ii).per_ii=per_ii;
                    handles_out.RPtimecourse(handles_out.RP_ii).group_no=PRPtimecourse.mouse(mouseNo).group_no;
                    handles_out.RPtimecourse(handles_out.RP_ii).mouseNo=mouseNo;
                    
                    data_calculated=1;
                    
                    %Peak p values
                    SpPRPtimecourse_peak=[];
                    ii_Sp=0;
                    SmPRPtimecourse_peak=[];
                    ii_Sm=0;
                    for elecNo=which_electrodes
                        for this_session=1:PRPtimecourse.mouse(mouseNo).no_sessions
                            if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed==1
                                
                                %Splus
                                evNo=1;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).peak)
                                            this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                            this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).peak;
                                            ii_Sp=ii_Sp+1;
                                            SpPRPtimecourse_peak(ii_Sp,:)=this_PRPtimecourse;
                                        end
                                    end
                                end
                                
                                %Sminus
                                evNo=2;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).peak)
                                            this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                            this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).peak;
                                            ii_Sm=ii_Sm+1;
                                            SmPRPtimecourse_peak(ii_Sm,:)=this_PRPtimecourse;
                                        end
                                    end
                                end
                            end
                        end
                    end
                    handles_out.RPtimecourse(handles_out.RP_ii).p_vals_peak=[];
                    handles_out.RPtimecourse(handles_out.RP_ii).t_sig_peak=[];
                    
                    if (size(SpPRPtimecourse_peak,1)>10)&(size(SmPRPtimecourse_peak,1)>10)
                        p_vals_peak=zeros(1,length(anal_t_pac));
                        for ii_t=1:length(anal_t_pac)
                            [h,p_vals_peak(ii_t)]=ttest2(SpPRPtimecourse_peak(:,ii_t),SmPRPtimecourse_peak(:,ii_t));
                        end
                        
                        handles_out.RPtimecourse(handles_out.RP_ii).p_vals_peak=p_vals_peak;
                        ii_zero=find(anal_t_pac>=t_odor_on,1,'first');
                        delta_ii_sig=find(p_vals_peak(ii_zero:end)<=0.05,1,'first');
                        handles_out.RPtimecourse(handles_out.RP_ii).t_sig_peak=anal_t_pac(ii_zero+delta_ii_sig-1)-t_odor_on;
                    else
                        data_calculated=0;
                    end
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSmPRPtimecourse_peak=mean(SmPRPtimecourse_peak,1)';
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSpPRPtimecourse_peak=mean(SpPRPtimecourse_peak,1)';
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).SmPRPtimecourse_peak=SmPRPtimecourse_peak;
                    handles_out.RPtimecourse(handles_out.RP_ii).SpPRPtimecourse_peak=SpPRPtimecourse_peak;
                    
                    %Peak below p values
                    SpPRPtimecourse_peak_below=[];
                    ii_Spb=0;
                    SmPRPtimecourse_peak_below=[];
                    ii_Smb=0;
                    for elecNo=which_electrodes
                        for this_session=1:PRPtimecourse.mouse(mouseNo).no_sessions
                            if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed==1
                                
                                %Splus
                                evNo=1;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'z_peakPower_timecourse_below')
                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_below)
                                                if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_below)
                                                    this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                                    this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_below;
                                                    ii_Spb=ii_Spb+1;
                                                    SpPRPtimecourse_peak_below(ii_Spb,:)=this_PRPtimecourse;
                                                end
                                            end
                                        end
                                    end
                                end
                                
                                %Sminus
                                evNo=2;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'z_peakPower_timecourse_below')
                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_below)
                                                if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_below)
                                                    this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                                    this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_below;
                                                    ii_Smb=ii_Smb+1;
                                                    SmPRPtimecourse_peak_below(ii_Smb,:)=this_PRPtimecourse;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    handles_out.RPtimecourse(handles_out.RP_ii).p_vals_peak_below=[];
                    handles_out.RPtimecourse(handles_out.RP_ii).t_sig_peak_below=[];
                    
                    if (size(SpPRPtimecourse_peak_below,1)>10)&(size(SmPRPtimecourse_peak_below,1)>10)
                        p_vals_peak_below=zeros(1,length(anal_t_pac));
                        for ii_t=1:length(anal_t_pac)
                            [h,p_vals_peak_below(ii_t)]=ttest2(SpPRPtimecourse_peak_below(:,ii_t),SmPRPtimecourse_peak_below(:,ii_t));
                        end
                        
                        handles_out.RPtimecourse(handles_out.RP_ii).p_vals_peak_below=p_vals_peak_below;
                        ii_zero=find(anal_t_pac>=t_odor_on,1,'first');
                        delta_ii_sig=find(p_vals_peak_below(ii_zero:end)<=0.05,1,'first');
                        handles_out.RPtimecourse(handles_out.RP_ii).t_sig_peak_below=anal_t_pac(ii_zero+delta_ii_sig-1)-t_odor_on;
                    else
                        data_calculated=0;
                    end
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSmPRPtimecourse_peak_below=mean(SmPRPtimecourse_peak_below,1)';
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSpPRPtimecourse_peak_below=mean(SpPRPtimecourse_peak_below,1)';
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).SmPRPtimecourse_peak_below=SmPRPtimecourse_peak_below;
                    handles_out.RPtimecourse(handles_out.RP_ii).SpPRPtimecourse_peak_below=SpPRPtimecourse_peak_below;
                    
                    %Peak above p values
                    SpPRPtimecourse_peak_above=[];
                    ii_Spb=0;
                    SmPRPtimecourse_peak_above=[];
                    ii_Smb=0;
                    for elecNo=which_electrodes
                        for this_session=1:PRPtimecourse.mouse(mouseNo).no_sessions
                            if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed==1
                                
                                %Splus
                                evNo=1;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'z_peakPower_timecourse_above')
                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_above)
                                                if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_above)
                                                    this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                                    this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_above;
                                                    if sum(isnan(this_PRPtimecourse))==0
                                                        ii_Spb=ii_Spb+1;
                                                        SpPRPtimecourse_peak_above(ii_Spb,:)=this_PRPtimecourse;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                                
                                %Sminus
                                evNo=2;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'z_peakPower_timecourse_above')
                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_above)
                                                if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_above)
                                                    this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                                    this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_above;
                                                    if sum(isnan(this_PRPtimecourse))==0
                                                        ii_Smb=ii_Smb+1;
                                                        SmPRPtimecourse_peak_above(ii_Smb,:)=this_PRPtimecourse;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    handles_out.RPtimecourse(handles_out.RP_ii).p_vals_peak_above=[];
                    handles_out.RPtimecourse(handles_out.RP_ii).t_sig_peak_above=[];
                     
                    if (size(SpPRPtimecourse_peak_above,1)>10)&(size(SmPRPtimecourse_peak_above,1)>10)
                        p_vals_peak_above=zeros(1,length(anal_t_pac));
                        for ii_t=1:length(anal_t_pac)
                            [h,p_vals_peak_above(ii_t)]=ttest2(SpPRPtimecourse_peak_above(:,ii_t),SmPRPtimecourse_peak_above(:,ii_t));
                        end
                        
                        handles_out.RPtimecourse(handles_out.RP_ii).p_vals_peak_above=p_vals_peak_above;
                        ii_zero=find(anal_t_pac>=t_odor_on,1,'first');
                        delta_ii_sig=find(p_vals_peak_above(ii_zero:end)<=0.05,1,'first');
                        handles_out.RPtimecourse(handles_out.RP_ii).t_sig_peak_above=anal_t_pac(ii_zero+delta_ii_sig-1)-t_odor_on;
                    else
                        data_calculated=0;
                    end
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSmPRPtimecourse_peak_above=mean(SmPRPtimecourse_peak_above,1)';
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSpPRPtimecourse_peak_above=mean(SpPRPtimecourse_peak_above,1)';

                    
                    handles_out.RPtimecourse(handles_out.RP_ii).SmPRPtimecourse_peak_above=SmPRPtimecourse_peak_above;
                    handles_out.RPtimecourse(handles_out.RP_ii).SpPRPtimecourse_peak_above=SpPRPtimecourse_peak_above;
                    
                    %Peak all p values
                    SpPRPtimecourse_peak_all=[];
                    ii_Spb=0;
                    SmPRPtimecourse_peak_all=[];
                    ii_Smb=0;
                    for elecNo=which_electrodes
                        for this_session=1:PRPtimecourse.mouse(mouseNo).no_sessions
                            if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed==1
                                
                                %Splus
                                evNo=1;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'z_peakPower_timecourse_all')
                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_all)
                                                if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_all)
                                                    this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                                    this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_all;
                                                    if sum(isnan(this_PRPtimecourse))==0
                                                        ii_Spb=ii_Spb+1;
                                                        SpPRPtimecourse_peak_all(ii_Spb,:)=this_PRPtimecourse;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                                
                                %Sminus
                                evNo=2;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'z_peakPower_timecourse_all')
                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_all)
                                                if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_all)
                                                    this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                                    this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).z_peakPower_timecourse_all;
                                                    if sum(isnan(this_PRPtimecourse))==0
                                                        ii_Smb=ii_Smb+1;
                                                        SmPRPtimecourse_peak_all(ii_Smb,:)=this_PRPtimecourse;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    handles_out.RPtimecourse(handles_out.RP_ii).p_vals_peak_all=[];
                    handles_out.RPtimecourse(handles_out.RP_ii).t_sig_peak_all=[];
                     
                    if (size(SpPRPtimecourse_peak_all,1)>10)&(size(SmPRPtimecourse_peak_all,1)>10)
                        p_vals_peak_all=zeros(1,length(anal_t_pac));
                        for ii_t=1:length(anal_t_pac)
                            [h,p_vals_peak_all(ii_t)]=ttest2(SpPRPtimecourse_peak_all(:,ii_t),SmPRPtimecourse_peak_all(:,ii_t));
                        end
                        
                        handles_out.RPtimecourse(handles_out.RP_ii).p_vals_peak_all=p_vals_peak_all;
                        ii_zero=find(anal_t_pac>=t_odor_on,1,'first');
                        delta_ii_sig=find(p_vals_peak_all(ii_zero:end)<=0.05,1,'first');
                        handles_out.RPtimecourse(handles_out.RP_ii).t_sig_peak_all=anal_t_pac(ii_zero+delta_ii_sig-1)-t_odor_on;
                    else
                        data_calculated=0;
                    end
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSmPRPtimecourse_peak_all=mean(SmPRPtimecourse_peak_all,1)';
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSpPRPtimecourse_peak_all=mean(SpPRPtimecourse_peak_all,1)';

                    
                    handles_out.RPtimecourse(handles_out.RP_ii).SmPRPtimecourse_peak_all=SmPRPtimecourse_peak_all;
                    handles_out.RPtimecourse(handles_out.RP_ii).SpPRPtimecourse_peak_all=SpPRPtimecourse_peak_all;
                    
                    
                    %Trough p values
                    SpPRPtimecourse_trough=[];
                    ii_Sp=0;
                    SmPRPtimecourse_trough=[];
                    ii_Sm=0;
                    for elecNo=which_electrodes
                        for this_session=1:PRPtimecourse.mouse(mouseNo).no_sessions
                            if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed==1
                                %Splus
                                evNo=1;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0 
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).trough)
                                            this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                            this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).trough;
                                            ii_Sp=ii_Sp+1;
                                            SpPRPtimecourse_trough(ii_Sp,:)=this_PRPtimecourse;
                                        end
                                    end
                                end
                                %Sminus
                                evNo=2;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).trough)
                                            this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                            this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).trough;
                                            ii_Sm=ii_Sm+1;
                                            SmPRPtimecourse_trough(ii_Sm,:)=this_PRPtimecourse;
                                        end
                                    end
                                end
                           
                            end
                        end
                    end
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).p_vals_trough=[];
                    handles_out.RPtimecourse(handles_out.RP_ii).t_sig_trough=[];
                    
                    if (size(SpPRPtimecourse_trough,1)>10)&(size(SmPRPtimecourse_trough,1)>10)
                        p_vals_trough=zeros(1,length(anal_t_pac));
                        for ii_t=1:length(anal_t_pac)
                            [h,p_vals_trough(ii_t)]=ttest2(SpPRPtimecourse_trough(:,ii_t),SmPRPtimecourse_trough(:,ii_t));
                        end
                        
                        handles_out.mouse(mouseNo).per_ii(per_ii).pacii(pacii).p_vals_trough=p_vals_trough;
                        
                        handles_out.RPtimecourse(handles_out.RP_ii).p_vals_trough=p_vals_trough;
                        ii_zero=find(anal_t_pac>=t_odor_on,1,'first');
                        delta_ii_sig=find(p_vals_trough(ii_zero:end)<=0.05,1,'first');
                        handles_out.RPtimecourse(handles_out.RP_ii).t_sig_trough=anal_t_pac(ii_zero+delta_ii_sig-1)-t_odor_on;
                    else
                        data_calculated=0;
                    end
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSmPRPtimecourse_trough=mean(SmPRPtimecourse_trough,1)';
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSpPRPtimecourse_trough=mean(SpPRPtimecourse_trough,1)';
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).SmPRPtimecourse_trough=SmPRPtimecourse_trough;
                    handles_out.RPtimecourse(handles_out.RP_ii).SpPRPtimecourse_trough=SpPRPtimecourse_trough;
                    
                    
                    %Lick power p values
                    SpPRPtimecourse_lickp=[];
                    ii_Sp=0;
                    SmPRPtimecourse_lickp=[];
                    ii_Sm=0;
                    for elecNo=which_electrodes
                        for this_session=1:PRPtimecourse.mouse(mouseNo).no_sessions
                            if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed==1
                                %Splus
                                evNo=1;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'lickp')
                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lickp)
                                                this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                                this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lickp;
                                                if sum(isnan(this_PRPtimecourse))==0
                                                    ii_Sp=ii_Sp+1;
                                                    SpPRPtimecourse_lickp(ii_Sp,:)=this_PRPtimecourse;
                                                end
                                            end
                                        end
                                    end
                                end
                                %Sminus
                                evNo=2;
                                if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).no_events>0
                                    
                                    for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power)
                                        if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo),'lickp')
                                            if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lickp)
                                                this_PRPtimecourse=zeros(1,length(anal_t_pac));
                                                this_PRPtimecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode(elecNo).event(evNo).z_power(trNo).lickp;
                                                if sum(isnan(this_PRPtimecourse))==0
                                                    ii_Sm=ii_Sm+1;
                                                    SmPRPtimecourse_lickp(ii_Sm,:)=this_PRPtimecourse;
                                                end
                                            end
                                        end
                                    end
                                end
                                
                            end
                        end
                    end
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).p_vals_lickp=[];
                    handles_out.RPtimecourse(handles_out.RP_ii).t_sig_lickp=[];
                    
                    handles_out.mouse(mouseNo).per_ii(per_ii).pacii(pacii).p_vals_lickp=[];
                    if (size(SpPRPtimecourse_lickp,1)>10)&(size(SmPRPtimecourse_lickp,1)>10)
                        p_vals_lickp=zeros(1,length(anal_t_pac));
                        for ii_t=1:length(anal_t_pac)
                            [h,p_vals_lickp(ii_t)]=ttest2(SpPRPtimecourse_lickp(:,ii_t),SmPRPtimecourse_lickp(:,ii_t));
                        end
                        
                        handles_out.RPtimecourse(handles_out.RP_ii).p_vals_lickp=p_vals_lickp;
                        ii_zero=find(anal_t_pac>=t_odor_on,1,'first');
                        delta_ii_sig=find(p_vals_lickp(ii_zero:end)<=0.05,1,'first');
                        handles_out.RPtimecourse(handles_out.RP_ii).t_sig_lickp=anal_t_pac(ii_zero+delta_ii_sig-1)-t_odor_on;
                    else
                        data_calculated=0;
                    end
                    
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSmPRPtimecourse_lickp=mean(SmPRPtimecourse_lickp,1)';
                    handles_out.RPtimecourse(handles_out.RP_ii).meanSpPRPtimecourse_lickp=mean(SpPRPtimecourse_lickp,1)';
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).SmPRPtimecourse_lickp=SmPRPtimecourse_lickp;
                    handles_out.RPtimecourse(handles_out.RP_ii).SpPRPtimecourse_lickp=SpPRPtimecourse_lickp;
                    
                    pffft=1;
                    
                    %Lick p values
                    SpPRPtimecourse=[];
                    ii_Sp=0;
                    SmPRPtimecourse=[];
                    ii_Sm=0;
                    %                     for elecNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).pacii(pacii).electrode)
                    for this_session=1:PRPtimecourse.mouse(mouseNo).no_sessions
                        if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).per_ii_processed==1
                            if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii),'event')
                                if ~isempty(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event)
                                    if length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event)>1
                                    %Splus
                                    evNo=1;
                                    if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).no_lick_events>0
                                        for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).lick_tr)
                                            this_lick_bin_timecourse=zeros(1,length(anal_t_pac));
                                            if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).lick_tr(trNo),'binary_lick_per_t')
                                                this_lick_bin_timecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).lick_tr(trNo).binary_lick_per_t;
                                                ii_Sp=ii_Sp+1;
                                                SpPRPtimecourse(ii_Sp,:)=this_lick_bin_timecourse;
                                            end
                                        end
                                    end
                                    %Sminus
                                    evNo=2;
                                    if PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).no_lick_events>0
                                        for trNo=1:length(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).lick_tr)
                                            this_lick_bin_timecourse=zeros(1,length(anal_t_pac));
                                            if isfield(PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).lick_tr(trNo),'binary_lick_per_t')
                                                this_lick_bin_timecourse(1,:)=PRPtimecourse.mouse(mouseNo).session(this_session).per_ii(per_ii).event(evNo).lick_tr(trNo).binary_lick_per_t;
                                                ii_Sm=ii_Sm+1;
                                                SmPRPtimecourse(ii_Sm,:)=this_lick_bin_timecourse;
                                            end
                                        end
                                    end
                                end
                                end
                            end
                        end
                    end
                    
                    handles_out.RPtimecourse(handles_out.RP_ii).p_vals_licks=[];
                    handles_out.RPtimecourse(handles_out.RP_ii).t_sig_licks=[];
                    
                    if (size(SpPRPtimecourse,1)>10)&(size(SmPRPtimecourse,1)>10)
                        p_vals_licks=zeros(1,length(anal_t_pac));
                        for ii_t=1:length(anal_t_pac)
                            p_vals_licks(ii_t)=ranksum(SpPRPtimecourse(:,ii_t),SmPRPtimecourse(:,ii_t));
                        end
                        
                        handles_out.RPtimecourse(handles_out.RP_ii).p_vals_licks=p_vals_licks;
                        ii_zero=find(anal_t_pac>=t_odor_on,1,'first');
                        delta_ii_sig=find(p_vals_licks(ii_zero:end)<=0.05,1,'first');
                        handles_out.RPtimecourse(handles_out.RP_ii).t_sig_licks=anal_t_pac(ii_zero+delta_ii_sig-1)-t_odor_on;
                    else
                        data_calculated=0;
                    end
                    
                    
                    %Plot the per mouse figures and calculate the decision making times
                    plot_figure=0;
                    if (data_calculated==1)&(per_ii==1)
                        
                        if plot_figure==1
                            %Plot the p values
                            
                            try
                                close(figureNo)
                            catch
                            end
                            hFig=figure(figureNo);
                            
                            set(hFig, 'units','normalized','position',[.1 .1 .7 .7])
                            %             subplot(2,1,1)
                            ax=gca;ax.LineWidth=3;
                            
                            hold on
                            plot(anal_t_pac,log10(p_vals_licks),'-k','LineWidth',3)
                            %                         plot(anal_t_pac,log10(p_vals_lickp),'-m','LineWidth',3)
                            plot(anal_t_pac,log10(p_vals_trough),'-b','LineWidth',3)
                            plot(anal_t_pac,log10(p_vals_peak),'-r','LineWidth',3)
                            
                            plot([anal_t_pac(1) anal_t_pac(end)],[log10(0.05) log10(0.05)],'-r','LineWidth', 2)
                        end
                        
                        %Find discrimination times and p value timecourses
                        
                        %Licks
                        ii_t0=find(anal_t_pac>=t_odor_on,1,'first');
                        ii=ii_t0+1;
                        not_found=1;
                        ii_cross=ii;
                        delta_ii_dt=ceil(0.3/0.0333);
                        while (ii<=ii_t0+50)&(not_found==1)
                            if (log10(p_vals_licks(ii-1))>=log10(0.05))&(log10(p_vals_licks(ii))<=log10(0.05))
                               if sum(log10(p_vals_licks(ii:ii+delta_ii_dt))>log10(0.05))==0
                                   ii_cross=ii;
                                   not_found=0;
                               end
                            end
                            ii=ii+1;
                        end
                        
                        if not_found==0
                            group_no=PRPtimecourse.mouse(mouseNo).group_no;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).no_licks=handles_out.disc_t.pacii(pacii).groupNo(PRPtimecourse.mouse(mouseNo).group_no).no_licks+1;
                            no_licks=handles_out.disc_t.pacii(pacii).groupNo(group_no).no_licks;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).licks(no_licks).disc_t=anal_t_pac(ii_cross)-t_odor_on;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).licks(no_licks).p_vals=log10(p_vals_licks);
                            
                            if plot_figure==1
                                plot(anal_t_pac(ii_cross),log10(0.05),'ok','MarkerFaceColor','k','MarkerSize',8)
                            end
                        end
                        
                        %Peak PRP
                        ii_t0=find(anal_t_pac>=t_odor_on,1,'first');
                        ii=ii_t0+1;
                        not_found=1;
                        ii_cross=ii;
                        delta_ii_dt=ceil(0.3/0.0333);
                        while (ii<=ii_t0+50)&(not_found==1)
                            if (log10(p_vals_peak(ii-1))>=log10(0.05))&(log10(p_vals_peak(ii))<=log10(0.05))
                               if sum(log10(p_vals_peak(ii:ii+delta_ii_dt))>log10(0.05))==0
                                   ii_cross=ii;
                                   not_found=0;
                               end
                            end
                            ii=ii+1;
                        end
                        
                        if not_found==0
                            group_no=PRPtimecourse.mouse(mouseNo).group_no;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).no_peak=handles_out.disc_t.pacii(pacii).groupNo(PRPtimecourse.mouse(mouseNo).group_no).no_peak+1;
                            no_peak=handles_out.disc_t.pacii(pacii).groupNo(group_no).no_peak;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).peak(no_peak).disc_t=anal_t_pac(ii_cross)-t_odor_on;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).peak(no_peak).p_vals=log10(p_vals_peak);
                            
                            if plot_figure==1
                                plot(anal_t_pac(ii_cross),log10(0.05),'or','MarkerFaceColor','r','MarkerSize',10)
                            end
                        end
                        
                        %Trough PRP
                        ii_t0=find(anal_t_pac>=t_odor_on,1,'first');
                        ii=ii_t0+1;
                        not_found=1;
                        ii_cross=ii;
                        delta_ii_dt=ceil(0.3/0.0333);
                        while (ii<=ii_t0+50)&(not_found==1)
                            if (log10(p_vals_trough(ii-1))>=log10(0.05))&(log10(p_vals_trough(ii))<=log10(0.05))
                               if sum(log10(p_vals_trough(ii:ii+delta_ii_dt))>log10(0.05))==0
                                   ii_cross=ii;
                                   not_found=0;
                               end
                            end
                            ii=ii+1;
                        end  
                        
                        if not_found==0
                            group_no=PRPtimecourse.mouse(mouseNo).group_no;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).no_trough=handles_out.disc_t.pacii(pacii).groupNo(PRPtimecourse.mouse(mouseNo).group_no).no_trough+1;
                            no_trough=handles_out.disc_t.pacii(pacii).groupNo(group_no).no_trough;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).trough(no_trough).disc_t=anal_t_pac(ii_cross)-t_odor_on;
                            handles_out.disc_t.pacii(pacii).groupNo(group_no).trough(no_trough).p_vals=log10(p_vals_trough);
                            
                            if plot_figure==1
                                plot(anal_t_pac(ii_cross),log10(0.05),'ob','MarkerFaceColor','b','MarkerSize',7)
                            end
                        end
                        
%                         if plot_figure==1
%                             xlabel('Time (sec)')
%                             ylabel('log10(p)')
%                             ylim([-10 0.1])
%                             title(['z score t test log10(p) for mouse no ' num2str(mouseNo) ' ' freq_names{pacii+1} ' ' prof_naive_leg{per_ii} ' group No ' num2str(PRPtimecourse.mouse(mouseNo).group_no)])
%                             
%                             %Plot the bounded line averages
%                             
%                             try
%                                 close(figureNo+1)
%                             catch
%                             end
%                             hFig=figure(figureNo+1);
%                             
%                             set(hFig, 'units','normalized','position',[.1 .2 .8 .3])
%                             
%                             %Show the average peak PRP dynamics
%                             %S+, S-
%                             subplot(1,3,1)
%                             hold on
%                             ax=gca;ax.LineWidth=3;
%                             
%                             CIsm = bootci(1000, @mean, SmPRPtimecourse_peak);
%                             meansm=mean(SmPRPtimecourse_peak,1);
%                             CIsm(1,:)=meansm-CIsm(1,:);
%                             CIsm(2,:)=CIsm(2,:)-meansm;
%                             
%                             CIsp = bootci(1000, @mean, SpPRPtimecourse_peak);
%                             meansp=mean(SpPRPtimecourse_peak,1);
%                             CIsp(1,:)=meansp-CIsp(1,:);
%                             CIsp(2,:)=CIsp(2,:)-meansp;
%                             
%                             
%                             [hlsm, hpsm] = boundedline(anal_t_pac',mean(SmPRPtimecourse_peak,1)', CIsm', 'b');
%                             [hlsp, hpsp] = boundedline(anal_t_pac',mean(SpPRPtimecourse_peak,1)', CIsp', 'r');
%                             
%                             ylim([-1.6 0.8])
%                             xlabel('Time (sec)')
%                             ylabel('PRP z score')
%                             
%                             
%                             title(['PRP peak z score for mouse no ' num2str(mouseNo) ' ' freq_names{pacii+1} ' ' prof_naive_leg{per_ii} ' group No ' num2str(PRPtimecourse.mouse(mouseNo).group_no)])
%                             
%                             %Show the average trough PRP dynamics
%                             %S+, S-
%                             subplot(1,3,2)
%                             hold on
%                             ax=gca;ax.LineWidth=3;
%                             
%                             CIsm = bootci(1000, @mean, SmPRPtimecourse_trough);
%                             meansm=mean(SmPRPtimecourse_trough,1);
%                             CIsm(1,:)=meansm-CIsm(1,:);
%                             CIsm(2,:)=CIsm(2,:)-meansm;
%                             
%                             CIsp = bootci(1000, @mean, SpPRPtimecourse_trough);
%                             meansp=mean(SpPRPtimecourse_trough,1);
%                             CIsp(1,:)=meansp-CIsp(1,:);
%                             CIsp(2,:)=CIsp(2,:)-meansp;
%                             
%                             
%                             [hlsm, hpsm] = boundedline(anal_t_pac',mean(SmPRPtimecourse_trough,1)', CIsm', 'b');
%                             [hlsp, hpsp] = boundedline(anal_t_pac',mean(SpPRPtimecourse_trough,1)', CIsp', 'r');
%                             
%                             ylim([-1.6 0.8])
%                             xlabel('Time (sec)')
%                             ylabel('PRP z score')
%                             
%                             
%                             title(['PRP trough z score for mouse no ' num2str(mouseNo) ' ' freq_names{pacii+1} ' ' prof_naive_leg{per_ii} ' group No ' num2str(PRPtimecourse.mouse(mouseNo).group_no)])
%                             
%                             %Show the average LRP dynamics
%                             %S+, S-
%                             subplot(1,3,3)
%                             hold on
%                             ax=gca;ax.LineWidth=3;
%                             
%                             CIsm = bootci(1000, @mean, SmPRPtimecourse_lickp);
%                             meansm=mean(SmPRPtimecourse_lickp,1);
%                             CIsm(1,:)=meansm-CIsm(1,:);
%                             CIsm(2,:)=CIsm(2,:)-meansm;
%                             
%                             CIsp = bootci(1000, @mean, SpPRPtimecourse_lickp);
%                             meansp=mean(SpPRPtimecourse_lickp,1);
%                             CIsp(1,:)=meansp-CIsp(1,:);
%                             CIsp(2,:)=CIsp(2,:)-meansp;
%                             
%                             
%                             [hlsm, hpsm] = boundedline(anal_t_pac',mean(SmPRPtimecourse_lickp,1)', CIsm', 'b');
%                             [hlsp, hpsp] = boundedline(anal_t_pac',mean(SpPRPtimecourse_lickp,1)', CIsp', 'r');
%                             
%                             ylim([-1.6 0.8])
%                             xlabel('Time (sec)')
%                             ylabel('LRP z score')
%                             
%                             
%                             title(['LRP z score for mouse no ' num2str(mouseNo) ' ' freq_names{pacii+1} ' ' prof_naive_leg{per_ii} ' group No ' num2str(PRPtimecourse.mouse(mouseNo).group_no)])
%                         end
                         pffft=1;
                    end
                   
                end
            end
        end
        
        %Peak PRP
        for pacii=1:no_pacii
            figureNo = figureNo + 1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            for per_ii=2:-1:1
                
                for groupNo=1:3
                    subplot(2,3,groupNo + 3*(2-per_ii))
                    hold on
                    
                    all_meanSmPRPtimecourse_peak=[];
                    all_meanSpPRPtimecourse_peak=[];
                    ii_tcs=0;
                    for ii=1:handles_out.RP_ii
                        if (handles_out.RPtimecourse(ii).pacii==pacii)&(handles_out.RPtimecourse(ii).per_ii==per_ii)&(handles_out.RPtimecourse(ii).group_no==groupNo)
                            if (length(handles_out.RPtimecourse(ii).meanSmPRPtimecourse_peak)>1)&(length(handles_out.RPtimecourse(ii).meanSpPRPtimecourse_peak)>1)
                                ii_tcs=ii_tcs+1;
                                all_meanSmPRPtimecourse_peak(ii_tcs,:)=handles_out.RPtimecourse(ii).meanSmPRPtimecourse_peak;
                                all_meanSpPRPtimecourse_peak(ii_tcs,:)=handles_out.RPtimecourse(ii).meanSpPRPtimecourse_peak;
                            end
                        end
                    end
                    fprintf(1, ['The number of mice included for ' prof_naive_leg{per_ii} ' ' group_legend{groupNo} ' is %d\n'], ii_tcs)
                    
                    CIsp = bootci(1000, @mean, all_meanSpPRPtimecourse_peak);
                    meansp=mean(all_meanSpPRPtimecourse_peak,1);
                    CIsp(1,:)=meansp-CIsp(1,:);
                    CIsp(2,:)=CIsp(2,:)-meansp;
                    
                    
                    [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_meanSpPRPtimecourse_peak,1)', CIsp', 'r');
                    
                    CIsm = bootci(1000, @mean, all_meanSmPRPtimecourse_peak);
                    meansm=mean(all_meanSmPRPtimecourse_peak,1);
                    CIsm(1,:)=meansm-CIsm(1,:);
                    CIsm(2,:)=CIsm(2,:)-meansm;
                    
                    [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_meanSmPRPtimecourse_peak,1)', CIsm', 'b');
                    
                    title([group_legend{groupNo} ' ' prof_naive_leg{per_ii}])
                    xlabel('Time(sec)')
                    ylabel('z')
                    ylim([-2 2])
                    
                end
            end
            
            sgtitle(['Peak PRP ' freq_names{pacii+1}])
        end
        
        %Peak PRP previous to lick
        for pacii=1:no_pacii
            figureNo = figureNo + 1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            for per_ii=2:-1:1
                
                for groupNo=1:3
                    subplot(2,3,groupNo + 3*(2-per_ii))
                    hold on
                    
                    all_meanSmPRPtimecourse_peak_below=[];
                    all_meanSpPRPtimecourse_peak_below=[];
                    ii_tcs=0;
                    for ii=1:handles_out.RP_ii
                        if (handles_out.RPtimecourse(ii).pacii==pacii)&(handles_out.RPtimecourse(ii).per_ii==per_ii)&(handles_out.RPtimecourse(ii).group_no==groupNo)
                            if (length(handles_out.RPtimecourse(ii).meanSmPRPtimecourse_peak_below)>1)&(length(handles_out.RPtimecourse(ii).meanSpPRPtimecourse_peak_below)>1)
                                ii_tcs=ii_tcs+1;
                                all_meanSmPRPtimecourse_peak_below(ii_tcs,:)=handles_out.RPtimecourse(ii).meanSmPRPtimecourse_peak_below;
                                all_meanSpPRPtimecourse_peak_below(ii_tcs,:)=handles_out.RPtimecourse(ii).meanSpPRPtimecourse_peak_below;
                            end
                        end
                    end
                    fprintf(1, ['The number of mice included for ' prof_naive_leg{per_ii} ' ' group_legend{groupNo} ' is %d\n'], ii_tcs)
                    
                    CIsp = bootci(1000, @mean, all_meanSpPRPtimecourse_peak_below);
                    meansp=mean(all_meanSpPRPtimecourse_peak_below,1);
                    CIsp(1,:)=meansp-CIsp(1,:);
                    CIsp(2,:)=CIsp(2,:)-meansp;
                    
                    
                    [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_meanSpPRPtimecourse_peak_below,1)', CIsp', 'r');
                    
                    CIsm = bootci(1000, @mean, all_meanSmPRPtimecourse_peak_below);
                    meansm=mean(all_meanSmPRPtimecourse_peak_below,1);
                    CIsm(1,:)=meansm-CIsm(1,:);
                    CIsm(2,:)=CIsm(2,:)-meansm;
                    
                    [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_meanSmPRPtimecourse_peak_below,1)', CIsm', 'b');
                    
                    title([group_legend{groupNo} ' ' prof_naive_leg{per_ii}])
                    xlabel('Time(sec)')
                    ylabel('z')
                    ylim([-2 2])
                    
                end
            end
            
            sgtitle(['Peak PRP for peaks before licks ' freq_names{pacii+1}])
        end
        
        %Peak PRP after lick
        for pacii=1:no_pacii
            figureNo = figureNo + 1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            for per_ii=2:-1:1
                
                for groupNo=1:3
                    subplot(2,3,groupNo + 3*(2-per_ii))
                    hold on
                    
                    all_meanSmPRPtimecourse_peak_above=[];
                    all_meanSpPRPtimecourse_peak_above=[];
                    ii_tcs=0;
                    for ii=1:handles_out.RP_ii
                        if (handles_out.RPtimecourse(ii).pacii==pacii)&(handles_out.RPtimecourse(ii).per_ii==per_ii)&(handles_out.RPtimecourse(ii).group_no==groupNo)
                            if (length(handles_out.RPtimecourse(ii).meanSmPRPtimecourse_peak_above)>1)&(length(handles_out.RPtimecourse(ii).meanSpPRPtimecourse_peak_above)>1)
                                ii_tcs=ii_tcs+1;
                                all_meanSmPRPtimecourse_peak_above(ii_tcs,:)=handles_out.RPtimecourse(ii).meanSmPRPtimecourse_peak_above;
                                all_meanSpPRPtimecourse_peak_above(ii_tcs,:)=handles_out.RPtimecourse(ii).meanSpPRPtimecourse_peak_above;
                            end
                        end
                    end
                    fprintf(1, ['The number of mice included for ' prof_naive_leg{per_ii} ' ' group_legend{groupNo} ' is %d\n'], ii_tcs)
                     
                    try
                        if size(all_meanSpPRPtimecourse_peak_above,1)>2
                            CIsp = bootci(1000, @mean, all_meanSpPRPtimecourse_peak_above);
                            meansp=mean(all_meanSpPRPtimecourse_peak_above,1);
                            CIsp(1,:)=meansp-CIsp(1,:);
                            CIsp(2,:)=CIsp(2,:)-meansp;
                            
                            
                            [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_meanSpPRPtimecourse_peak_above,1)', CIsp', 'r');
                            
                        else
                            plot(anal_t_pac',mean(all_meanSpPRPtimecourse_peak_above,1)', 'r');
                        end
                        
                        if size(all_meanSmPRPtimecourse_peak_above,1)>2
                            CIsm = bootci(1000, @mean, all_meanSmPRPtimecourse_peak_above);
                            meansm=mean(all_meanSmPRPtimecourse_peak_above,1);
                            CIsm(1,:)=meansm-CIsm(1,:);
                            CIsm(2,:)=CIsm(2,:)-meansm;
                            
                            [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_meanSmPRPtimecourse_peak_above,1)', CIsm', 'b');
                        else
                            plot(anal_t_pac',mean(all_meanSmPRPtimecourse_peak_above,1)', 'b');
                        end
                    catch
                    end
                    
                  
                    title([group_legend{groupNo} ' ' prof_naive_leg{per_ii}])
                    xlabel('Time(sec)')
                    ylabel('z')
                    ylim([-2 2])
                    
                end
            end
            
            sgtitle(['Peak PRP for peaks after licks ' freq_names{pacii+1}])
        end
        
        %Peak PRP all lick
        for pacii=1:no_pacii
            figureNo = figureNo + 1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            for per_ii=2:-1:1
                
                for groupNo=1:3
                    subplot(2,3,groupNo + 3*(2-per_ii))
                    hold on
                    
                    all_meanSmPRPtimecourse_peak_all=[];
                    all_meanSpPRPtimecourse_peak_all=[];
                    ii_tcs=0;
                    for ii=1:handles_out.RP_ii
                        if (handles_out.RPtimecourse(ii).pacii==pacii)&(handles_out.RPtimecourse(ii).per_ii==per_ii)&(handles_out.RPtimecourse(ii).group_no==groupNo)
                            if (length(handles_out.RPtimecourse(ii).meanSmPRPtimecourse_peak_all)>1)&(length(handles_out.RPtimecourse(ii).meanSpPRPtimecourse_peak_all)>1)
                                ii_tcs=ii_tcs+1;
                                all_meanSmPRPtimecourse_peak_all(ii_tcs,:)=handles_out.RPtimecourse(ii).meanSmPRPtimecourse_peak_all;
                                all_meanSpPRPtimecourse_peak_all(ii_tcs,:)=handles_out.RPtimecourse(ii).meanSpPRPtimecourse_peak_all;
                            end
                        end
                    end
                    fprintf(1, ['The number of mice included for ' prof_naive_leg{per_ii} ' ' group_legend{groupNo} ' is %d\n'], ii_tcs)
                     
                    
                    if size(all_meanSpPRPtimecourse_peak_all,1)>2
                        CIsp = bootci(1000, @mean, all_meanSpPRPtimecourse_peak_all);
                        meansp=mean(all_meanSpPRPtimecourse_peak_all,1);
                        CIsp(1,:)=meansp-CIsp(1,:);
                        CIsp(2,:)=CIsp(2,:)-meansp;
                        
                        
                        [hlsp, hpsp] = boundedline(anal_t_pac',mean(all_meanSpPRPtimecourse_peak_all,1)', CIsp', 'r');
                        
                    else
                        plot(anal_t_pac',mean(all_meanSpPRPtimecourse_peak_all,1)', 'r');
                    end
                    
                    if size(all_meanSmPRPtimecourse_peak_all,1)>2
                        CIsm = bootci(1000, @mean, all_meanSmPRPtimecourse_peak_all);
                        meansm=mean(all_meanSmPRPtimecourse_peak_all,1);
                        CIsm(1,:)=meansm-CIsm(1,:);
                        CIsm(2,:)=CIsm(2,:)-meansm;
                        
                        [hlsm, hpsm] = boundedline(anal_t_pac',mean(all_meanSmPRPtimecourse_peak_all,1)', CIsm', 'b');
                    else
                        plot(anal_t_pac',mean(all_meanSmPRPtimecourse_peak_all,1)', 'b');
                    end
                    
                   
                  
                    title([group_legend{groupNo} ' ' prof_naive_leg{per_ii}])
                    xlabel('Time(sec)')
                    ylabel('z')
                    ylim([-2 2])
                    
                end
            end
            
            sgtitle(['Peak PRP for all peaks referenced to licks ' freq_names{pacii+1}])
        end
        save([handles.PathName handles.drgb.outFileName(1:end-4) '_case' num2str(which_display) '_' handles_pars.output_suffix],'handles_out')
        pffft=1; 
        
    case 19
        % 19 PAC MI analysis for events (concentrations or S+/S-) for naive and proficient
        % Analyzed per mouse for groups defined by the user
        
        mean_MI_No_per_mouse=0;
        
        mean_MI_No=0;
        mean_MI=[];
        mean_MI_perii=[];
        mean_MI_evNo=[];
        mean_MI_pacii=[];
        mean_MI_fileNo=[];
        per_session_group_no=[];
        mean_VL=[];
        mean_VA=[];
        mean_PA=[];
        
        handles_out=[];
        handles_out.mi_ii=0;
        handles_out.PA_ii=0;
        
        prof_naive_leg{1}='Proficient';
        prof_naive_leg{2}='Naive';
        
        group_legend{1}='WT';
        group_legend{2}='Het';
        group_legend{3}='KO';
        
        fprintf(1, ['PAC analysis for Justin''s paper\n\n'])
        p_vals=[];
        no_files=max(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        
        szpc=size(percent_windows);
        for per_ii=1:szpc(1)
            
            
            
            for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
                mouse_has_files=0;
                for elec=1:16
                    if sum(which_electrodes==elec)>0
                        
                        theseEvNos_thisMouse_thisElec=[];
                        for evNo=1:length(eventType)
                            for pacii=1:no_pacii
                                for group_no=1:max(handles_drgb.drgbchoices.group_no)
                                    theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv=0;
                                end
                            end
                        end
                        
                        
                        for fileNo=1:no_files
                            %If this file is in the list of files the user wants to process in drgAnalysisBatchLFP continue
                            if sum(files==fileNo)>0
                                if handles_drgb.drgbchoices.mouse_no(fileNo)==mouseNo
                                    lfpodNo=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                    
                                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo)))
                                        
                                        
                                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo).PAC))
                                            
                                            percent_mask=[];
                                            percent_mask=logical((handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).perCorrPAC>=percent_windows(per_ii,1))...
                                                &(handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).perCorrPAC<=percent_windows(per_ii,2)));
                                            
                                            if ~isempty(percent_mask)
                                                
                                                
                                                for evNo=1:length(eventType)
                                                    
                                                    noWB_for_evNo(evNo)=-1;
                                                    
                                                    trials_in_event_Ev=[];
                                                    trials_in_event_Ev=logical((handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).which_eventPAC(eventType(evNo),:)==1)&percent_mask);
                                                    
                                                    if (sum(trials_in_event_Ev)>=1)
                                                        
                                                        %Do per bandwidth analysis
                                                        for pacii=1:no_pacii
                                                            
                                                            group_no=handles_drgb.drgbchoices.group_no(fileNo);
                                                            
                                                            %Enter the modulation index
                                                            this_MI_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_MI_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).mod_indx(trials_in_event_Ev);
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_MI_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=this_MI_Ev;
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).whichMouse(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=mouseNo*ones(1,length(this_MI_Ev));
                                                            
                                                            %Save per session
                                                            %value for MI
                                                            mean_MI_No=mean_MI_No+1;
                                                            mean_MI(mean_MI_No)=mean(this_MI_Ev);
                                                            mean_MI_perii(mean_MI_No)=per_ii;
                                                            mean_MI_evNo(mean_MI_No)=evNo;
                                                            mean_MI_pacii(mean_MI_No)=pacii;
                                                            mean_MI_fileNo(mean_MI_No)=fileNo;
                                                            mean_MI_mouse(mean_MI_No)=mouseNo;
                                                            per_session_group_no(mean_MI_No)=handles_drgb.drgbchoices.group_no(fileNo);
                                                            
                                                            
                                                            %                                                             if mean_MI(mean_MI_No)>=0.035
                                                            %                                                                 fprintf(1, ['MI larger than 0.035 for mouse no %d, file no %d, electrode, %d, pac no %d, perii %d, conc, %d\n'],mouseNo, fileNo, elec, pacii, per_ii, evNo);
                                                            %                                                             end
                                                            
                                                            %Enter the meanVectorLength
                                                            this_meanVectorLength_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_meanVectorLength_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).meanVectorLength(trials_in_event_Ev);
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_meanVectorLength_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=this_meanVectorLength_Ev;
                                                            
                                                            mean_VL(mean_MI_No)=mean(this_meanVectorLength_Ev);
                                                            
                                                            %Enter the meanVectorAngle
                                                            this_meanVectorAngle_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_meanVectorAngle_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).meanVectorAngle(trials_in_event_Ev);
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_meanVectorAngle_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=this_meanVectorAngle_Ev;
                                                            
                                                            mean_VA(mean_MI_No)=mean(this_meanVectorAngle_Ev);
                                                            
                                                            %Enter the peakAngle
                                                            this_peakAngle_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_peakAngle_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).peakAngle(trials_in_event_Ev);
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_peakAngle_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=this_peakAngle_Ev;
                                                            
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev);
                                                            
                                                            mean_PA(mean_MI_No)=mean(this_peakAngle_Ev);
                                                            
                                                            mouse_has_files=1;
                                                            
                                                        end
                                                        
                                                        
                                                        fprintf(1, ['%d trials in event No %d succesfully processed for file No %d electrode %d\n'],sum(trials_in_event_Ev), evNo,fileNo,elec);
                                                        
                                                    else
                                                        
                                                        
                                                        fprintf(1, ['%d trials in event No %d fewer than minimum trials per event ' evTypeLabels{evNo} ' for file No %d electrode %d\n'],sum(trials_in_event_Ev), min_trials_per_event,fileNo,elec);
                                                        
                                                        
                                                    end
                                                    
                                                    
                                                end
                                                
                                            else
                                                
                                                fprintf(1, ['Empty percent_mask for file No %d electrode %d\n'],fileNo,elec);
                                                
                                            end
                                            
                                        else
                                            
                                            fprintf(1, ['Empty allPower for file No %d electrode %d\n'],fileNo,elec);
                                            
                                        end
                                        
                                        
                                    else
                                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],fileNo,elec);
                                        
                                        
                                    end
                                end %if mouseNo
                            end
                        end %fileNo
                        
                        if theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv>0
                            
                            %Calculate per mouse PAC measures
                            for evNo=1:length(eventType)
                                for pacii=1:no_pacii
                                    for group_no=1:max(handles_drgb.drgbchoices.group_no)
                                        %Calculate per mouse MI
                                        mean_MI_No_per_mouse=mean_MI_No_per_mouse+1;
                                        this_mouse_MI=[];
                                        this_mouse_MI=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_MI_Ev;
                                        if ~isempty(this_mouse_MI)
                                            mean_MI_per_mouse(mean_MI_No_per_mouse)=mean(this_mouse_MI);
                                            if mean(this_mouse_MI)>0.01
                                                fprintf(1, ['MI larger than 0.01 for mouse no %d, electrode, %d, pac no %d, perii %d, conc, %d\n'],mouseNo, elec, pacii, per_ii, evNo);
                                            end
                                        else
                                            mean_MI_per_mouse(mean_MI_No_per_mouse)=NaN;
                                        end
                                        
                                        mean_MI_perii_per_mouse(mean_MI_No_per_mouse)=per_ii;
                                        mean_MI_evNo_per_mouse(mean_MI_No_per_mouse)=evNo;
                                        mean_MI_pacii_per_mouse(mean_MI_No_per_mouse)=pacii;
                                        mean_MI_mouseNo_per_mouse(mean_MI_No_per_mouse)=mouseNo;
                                        mean_MI_electNo_per_mouse(mean_MI_No_per_mouse)=elec;
                                        mean_MI_group_no_per_mouse(mean_MI_No_per_mouse)=group_no;
                                        
                                        %Calculate per mouse meanVectorLength
                                        this_mouse_meanVectorLength=[];
                                        this_mouse_meanVectorLength=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_meanVectorLength_Ev;
                                        if ~isempty(this_mouse_meanVectorLength)
                                            mean_meanVectorLength_per_mouse(mean_MI_No_per_mouse)=mean(this_mouse_meanVectorLength);
                                        else
                                            mean_meanVectorLength_per_mouse(mean_MI_No_per_mouse)=NaN;
                                        end
                                        
                                        %Calculate per mouse meanVectorAngle
                                        this_mouse_meanVectorAngle=[];
                                        this_mouse_meanVectorAngle=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_meanVectorAngle_Ev;
                                        if ~isempty(this_mouse_meanVectorAngle)
                                            mean_meanVectorAngle_per_mouse(mean_MI_No_per_mouse)=(180/pi)*circ_axial(circ_mean(this_mouse_meanVectorAngle'*pi/180)');
                                        else
                                            mean_meanVectorAngle_per_mouse(mean_MI_No_per_mouse)=NaN;
                                        end
                                        
                                        %Calculate per mouse peakAngle
                                        this_mouse_peakAngle=[];
                                        this_mouse_peakAngle=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_peakAngle_Ev;
                                        if ~isempty(this_mouse_peakAngle)
                                            mean_peakAngle_per_mouse(mean_MI_No_per_mouse)=(180/pi)*circ_axial(circ_mean(this_mouse_peakAngle'*pi/180)');
                                        else
                                            mean_peakAngle_per_mouse(mean_MI_No_per_mouse)=NaN;
                                        end
                                        
                                        %Calculate per mouse peakAngleVar
                                        this_mouse_peakAngle=[];
                                        this_mouse_peakAngle=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_peakAngle_Ev;
                                        if ~isempty(this_mouse_peakAngle)
                                            mean_peakAngleVar_per_mouse(mean_MI_No_per_mouse)=((180/pi)^2)*circ_var(this_mouse_peakAngle'*pi/180)';
                                        else
                                            mean_peakAngleVar_per_mouse(mean_MI_No_per_mouse)=NaN;
                                        end
                                    end
                                end
                            end
                        end
                        
                    end
                    
                end
                
                if mouse_has_files==1
                    mouse_included(mouseNo)=1;
                else
                    mouse_included(mouseNo)=0;
                end
            end
            
            
        end
        fprintf(1, '\n\n')
        
        fprintf(1, ['The number of mice included in the PAC analysis for this odor pair is %d\n\n\n'], sum(mouse_included))
        
        figNo=0;
        
        %Now plot the average MI for each electrode calculated per mouse
        %(including all sessions for each mouse)
        edges=[-15:0.5:100];
        rand_offset=0.8;
        
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            glm_mi=[];
            glm_ii=0;
            
            id_ii=0;
            input_data=[];
            
            %Plot the average
            figNo = figNo +1;
            
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            %             try
            %                 close(figNo+pacii)
            %             catch
            %             end
            %             hFig=figure(figNo+pacii);
            
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            hold on
            
            bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            bar_offset = 0;
            
            %             for grNo=1:max(handles_drgb.drgbchoices.group_no) lo puse
            %             abajo
            
            include_group=0;
            
            for evNo=1:length(eventType)
                
                for per_ii=2:-1:1
                    
                    for grNo=1:max(handles_drgb.drgbchoices.group_no)
                        bar_offset = bar_offset +1;
                        
                        %                         if sum(eventType==3)>0
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
                        %                         else
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        %                         end
                        %
                        %                         these_offsets(per_ii)=bar_offset;
                        bar_offset = bar_offset + 1;
                        
                        if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
                            
                            include_group=1;
                            
                            switch grNo
                                case 1
                                    bar(bar_offset,mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),'g','LineWidth', 3,'EdgeColor','none')
                                case 2
                                    bar(bar_offset,mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                                case 3
                                    bar(bar_offset,mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),'y','LineWidth', 3,'EdgeColor','none')
                            end
                            
                            handles_out.mi_ii=handles_out.mi_ii+1;
                            handles_out.mi_values(handles_out.mi_ii).pacii=pacii;
                            handles_out.mi_values(handles_out.mi_ii).evNo=evNo;
                            handles_out.mi_values(handles_out.mi_ii).per_ii=per_ii;
                            handles_out.mi_values(handles_out.mi_ii).groupNo=grNo;
                            handles_out.mi_values(handles_out.mi_ii).MI=mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)));
                            
                            %Save the mean per mouse
                            these_mice=mean_MI_mouseNo_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
                            handles_out.mi_values(handles_out.mi_ii).noMice=0;
                            for iiMice=min(these_mice):max(these_mice)
                                if sum(these_mice==iiMice)>0
                                    handles_out.mi_values(handles_out.mi_ii).noMice=handles_out.mi_values(handles_out.mi_ii).noMice+1;
                                    handles_out.mi_values(handles_out.mi_ii).mouseNo(handles_out.mi_values(handles_out.mi_ii).noMice)=iiMice;
                                    handles_out.mi_values(handles_out.mi_ii).MI_per_mouse(handles_out.mi_values(handles_out.mi_ii).noMice)=mean(mean_MI_per_mouse((mean_MI_mouseNo_per_mouse==iiMice)&(~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)));
                                end
                            end
                            
                            %Violin plot
                            [mean_out, CIout]=drgViolinPoint(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))...
                                ,edges,bar_offset,rand_offset,'k','k',1);
                            
                            %                                 %Save data for glm and ranksum
                            x_mi=mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
                            glm_mi.data(glm_ii+1:glm_ii+length(x_mi))=x_mi;
                            glm_mi.group(glm_ii+1:glm_ii+length(x_mi))=grNo*ones(1,length(x_mi));
                            glm_mi.perCorr(glm_ii+1:glm_ii+length(x_mi))=per_ii*ones(1,length(x_mi));
                            glm_mi.event(glm_ii+1:glm_ii+length(x_mi))=evNo*ones(1,length(x_mi));
                            glm_ii=glm_ii+length(x_mi);
                            
                            id_ii=id_ii+1;
                            input_data(id_ii).data=x_mi;
                            input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                            
                            %                             if per_ii==1
                            %                                 bar(bar_offset,mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),'r','LineWidth', 3,'EdgeColor','none')
                            %                             else
                            %                                 bar(bar_offset,mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                            %                             end
                            
                            
                            %                             plot(bar_offset,mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
                            %                             plot((bar_offset)*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),...
                            %                                 mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)),'o',...
                            %                                 'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            %
                            %                             if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>=2
                            %                                 CI = bootci(1000, {@mean, mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))},'type','cper');
                            %                                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            %                             end
                            
                            %Save data for anovan
                            %                             data_MI=[data_MI mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))];
                            %                             prof_naive=[prof_naive per_ii*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))];
                            %                             events=[events evNo*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))];
                            %                             groups=[groups grNo*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))];
                            %                             mice=[mice mean_MI_mouseNo_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))];
                            %
                        end
                    end
                    bar_offset = bar_offset + 2;
                    %                     if include_group==1
                    %                         bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                    %                         no_ev_labels=no_ev_labels+1;
                    %                         if sum(eventType==3)>0
                    %                             bar_labels{no_ev_labels}=evTypeLabels{evNo};
                    %                         else
                    %                             bar_labels{no_ev_labels}=num2str(concs(evNo));
                    %                         end
                    %                     end
                end
                bar_offset = bar_offset + 3;
                %                 if include_group==1
                %                     ii_gr_included=ii_gr_included+1;
                %                     groups_included(ii_gr_included)=grNo;
                %                 end
            end
            
            title(['Average MI for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1}])
            
            %              ylim([0 1.2*max(all_bar)])
            %Annotations identifying groups
            x_interval=0.8/ii_gr_included;
            for ii=1:ii_gr_included
                annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %             %x labels
            %             to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
            %             sorted_A=sortrows(to_sort);
            %             sorted_bar_lab_loc=sorted_A(:,1);
            %             for ii=1:length(bar_lab_loc)
            %                 sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
            %             end
            xticks([2 4 6 10 12 14 21 23 25 29 31 33])
            xticklabels({'nwtS+', 'nHS+', 'nKOS+', 'pwtS+', 'pHS+', 'pKOS+', 'nwtS-', 'nHS-', 'nKOS-', 'pwtS-', 'pHS-', 'pKOS-'})
            
            
            if sum(eventType==3)==0
                xlabel('Concentration (%)')
            end
            
            ylabel('Modulation Index')
            
            %Perform the glm
            fprintf(1, ['glm for average MI for each electrode calculated per mouse for PAC theta' freq_names{pacii+1} '\n'])
            
            fprintf(1, ['\n\nglm for MI for Theta/' freq_names{pacii+1} '\n'])
            tbl = table(glm_mi.data',glm_mi.group',glm_mi.perCorr',glm_mi.event',...
                'VariableNames',{'MI','group','perCorr','event'});
            mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
                ,'CategoricalVars',[2,3,4])
            
            
            %Do the ranksum/t-test
            fprintf(1, ['\n\nRanksum or t-test p values for average MI for each electrode calculated per mouse for PAC theta' freq_names{pacii+1} '\n'])
            [output_data] = drgMutiRanksumorTtest(input_data);
            
            
        end
        
        pfft=1;
        
        %         %Now do the cumulative histograms and ranksums for MI for each electrode calculated with all sessons per mouse
        %
        %         pvals=[];
        %         legends=[];
        %         out_mi_rank=[];
        %
        %         for pacii=1:no_pacii
        %
        %             glm_ii=0;
        %             glm_mi=[];
        %             glm_perm_ii=0;
        %             glm_mi_perm=[];
        %             ii_rank=0;
        %             maxmi=-200;
        %             minmi=200;
        %             input_data=[];
        %             ii=0;
        %             prof_naive_leg{1}='Proficient';
        %             prof_naive_leg{2}='Naive';
        %
        %             %Find out which groups the user is including
        %             group_included=zeros(1,max(handles_drgb.drgbchoices.group_no));
        %             for evNo=1:length(eventType)
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %                             group_included(grNo)=1;
        %                         end
        %                     end
        %                 end
        %             end
        %
        %             %If this is only one group plot for Justin's paper
        %             if sum(group_included)==1
        %                 figNo = figNo + 1;
        %                 try
        %                     close(figNo)
        %                 catch
        %                 end
        %                 hFig=figure(figNo);
        %
        %                 if length(eventType)>2
        %                     set(hFig, 'units','normalized','position',[.1 .1 .7 .8])
        %                 else
        %                     set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
        %                 end
        %                 hold on
        %
        %                 for evNo=1:length(eventType)
        %
        %
        %
        %                     for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                         for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                             if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                                 [f_mi,x_mi] = drg_ecdf(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)));
        %                                 cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_mi=f_mi;
        %                                 cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_mi=x_mi;
        %
        %                              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I
        %                              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%added
        %                              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%this
        %                                  switch grNo
        %                                 case 1
        %                                     if per_ii==1
        %                                         legends.pacii(pacii).evNo(evNo).per_ii(per_ii).p1=plot(x_mi,f_mi,'Color',[0 0 1],'LineWidth',3);
        %                                     else
        %                                         legends.pacii(pacii).evNo(evNo).per_ii(per_ii).p2=plot(x_mi,f_mi,'Color',[0.7 0.7 1],'LineWidth',3);
        %                                     end
        %                                 case 2
        %                                     if per_ii==1
        %                                         legends.pacii(pacii).evNo(evNo).per_ii(per_ii).p3=plot(x_mi,f_mi,'Color',[1 0 0],'LineWidth',3);
        %                                     else
        %                                         legends.pacii(pacii).evNo(evNo).per_ii(per_ii).p4=plot(x_mi,f_mi,'Color',[1 0.7 0.7],'LineWidth',3);
        %                                     end
        %                                 case 3
        %                                     if per_ii==1
        %                                         legends.pacii(pacii).evNo(evNo).per_ii(per_ii).p5=plot(x_mi,f_mi,'Color',[0 1 0],'LineWidth',3);
        %                                     else
        %                                         legends.pacii(pacii).evNo(evNo).per_ii(per_ii).p6=plot(x_mi,f_mi,'Color',[0.7 1 0.7],'LineWidth',3);
        %                                     end
        % %                                 if evNo==1
        % %                                     if per_ii==1
        % %                                         legends.pacii(pacii).evNo(evNo).per_ii(per_ii).p=plot(x_mi,f_mi,'Color',[1 0 0],'LineWidth',3);
        % %                                     else
        % %                                         legends.pacii(pacii).evNo(evNo).per_ii(per_ii).p=plot(x_mi,f_mi,'Color',[0 0 1],'LineWidth',3);
        % %                                     end
        % %                                 else
        % %                                     if per_ii==1
        % %                                         legends.pacii(pacii).evNo(evNo).per_ii(per_ii).p=plot(x_mi,f_mi,'Color',[1 0.7 0.7],'LineWidth',3);
        % %                                     else
        % %                                         legends.pacii(pacii).evNo(evNo).per_ii(per_ii).p=plot(x_mi,f_mi,'Color',[0.7 0.7 1],'LineWidth',3);
        % %                                     end
        %                                 end
        %
        %
        %                                 %Save data for glm and ranksum
        %                                 glm_mi.data(glm_ii+1:glm_ii+length(x_mi))=x_mi;
        %                                 glm_mi.group(glm_ii+1:glm_ii+length(x_mi))=grNo*ones(1,length(x_mi));
        %                                 glm_mi.perCorr(glm_ii+1:glm_ii+length(x_mi))=per_ii*ones(1,length(x_mi));
        %                                 glm_mi.event(glm_ii+1:glm_ii+length(x_mi))=evNo*ones(1,length(x_mi));
        %                                 glm_ii=glm_ii+length(x_mi);
        %
        %                                 ii=ii+1;
        %                                 input_data(ii).data=x_mi;
        %                                 input_data(ii).description=[evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
        %
        %                                 ii_rank=ii_rank+1;
        %                                 mi_rank(ii_rank).mi=mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %                                 mi_rank(ii_rank).per_ii=per_ii;
        %                                 mi_rank(ii_rank).grNo=grNo;
        %                                 mi_rank(ii_rank).evNo=evNo;
        %
        %                                 maxmi=max([maxmi max(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))]);
        %                                 minmi=min([minmi min(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))]);
        %                             end
        %                         end
        %
        %                     end
        %
        %                     xlabel('MI')
        %                     ylabel('Probability')
        %                 end
        %
        %
        %
        %
        %                     xlim([minmi-0.1*(maxmi-minmi) maxmi+0.1*(maxmi-minmi)])
        %
        %
        %                 %suptitle(['Average MI for each electrode calculated per  for PAC theta/' freq_names{pacii+1}])
        %
        %                 %Now do the ranksums
        %
        %
        %
        %
        %                 %Perform the glm
        %                 fprintf(1, ['glm for average MI for each electrode calculated per mouse for PAC theta' freq_names{pacii+1} '\n'])
        %
        %                 if sum(glm_mi.group==1)==length(glm_mi.group)
        %                     %There is only one group here (e.g. for Justin's paper we only include
        %                     %forward)
        %                     fprintf(1, ['\n\nglm for MI for Theta/' freq_names{pacii+1} '\n'])
        %                     tbl = table(glm_mi.data',glm_mi.perCorr',glm_mi.event',...
        %                         'VariableNames',{'MI','perCorr','event'});
        %                     mdl = fitglm(tbl,'MI~perCorr+event+perCorr*event'...
        %                         ,'CategoricalVars',[2,3])
        %                 else
        %
        %                     fprintf(1, ['\n\nglm for MI for Theta/' freq_names{pacii+1} '\n'])
        %                     tbl = table(glm_mi.data',glm_mi.group',glm_mi.perCorr',glm_mi.event',...
        %                         'VariableNames',{'MI','group','perCorr','event'});
        %                     mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
        %                         ,'CategoricalVars',[2,3,4])
        %                 end
        %
        %                 %Do the ranksum/t-test
        %                 fprintf(1, ['\n\nRanksum or t-test p values for average MI for each electrode calculated per mouse for PAC theta' freq_names{pacii+1} '\n'])
        %                 [output_data] = drgMutiRanksumorTtest(input_data);
        %
        %                 out_mi_rank(pacii).mi_rank=mi_rank;
        %             else
        %                 figNo = figNo + 1;
        %                 try
        %                     close(figNo)
        %                 catch
        %                 end
        %                 hFig=figure(figNo);
        %
        %                 if length(eventType)>2
        %                     set(hFig, 'units','normalized','position',[.1 .1 .7 .8])
        %                 else
        %                     set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
        %                 end
        %                 for evNo=1:length(eventType)
        %
        %                     subplot(ceil(length(eventType)/2),2,evNo)
        %                     hold on
        %
        %                     for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                         for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                             if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                                 [f_mi,x_mi] = drg_ecdf(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)));
        %                                 cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_mi=f_mi;
        %                                 cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_mi=x_mi;
        %                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I
        %                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%added
        %                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%this
        %                                  switch grNo
        %                                 case 1
        %                                     if per_ii==1
        %                                         legends.pacii(pacii).evNo(evNo).p1=plot(x_mi,f_mi,'Color',[0 0 1],'LineWidth',3);
        %                                     else
        %                                         legends.pacii(pacii).evNo(evNo).p2=plot(x_mi,f_mi,'Color',[0.7 0.7 1],'LineWidth',3);
        %                                     end
        %                                 case 2
        %                                     if per_ii==1
        %                                         legends.pacii(pacii).evNo(evNo).p3=plot(x_mi,f_mi,'Color',[1 0 0],'LineWidth',3);
        %                                     else
        %                                         legends.pacii(pacii).evNo(evNo).p4=plot(x_mi,f_mi,'Color',[1 0.7 0.7],'LineWidth',3);
        %                                     end
        %                                 case 3
        %                                     if per_ii==1
        %                                         legends.pacii(pacii).evNo(evNo).p5=plot(x_mi,f_mi,'Color',[0 1 0],'LineWidth',3);
        %                                     else
        %                                         legends.pacii(pacii).evNo(evNo).p6=plot(x_mi,f_mi,'Color',[0.7 1 0.7],'LineWidth',3);
        %                                     end
        % %                                 if grNo==1
        % %                                     if per_ii==1
        % %                                         legends.pacii(pacii).evNo(evNo).p1=plot(x_mi,f_mi,'Color',[1 0 0],'LineWidth',3);
        % %                                     else
        % %                                         legends.pacii(pacii).evNo(evNo).p2=plot(x_mi,f_mi,'Color',[0 0 1],'LineWidth',3);
        % %                                     end
        % %                                 else
        % %                                     if per_ii==1
        % %                                         legends.pacii(pacii).evNo(evNo).p3=plot(x_mi,f_mi,'Color',[1 0.7 0.7],'LineWidth',3);
        % %                                     else
        % %                                         legends.pacii(pacii).evNo(evNo).p4=plot(x_mi,f_mi,'Color',[0.7 0.7 1],'LineWidth',3);
        % %                                     end
        %                                 end
        %
        %
        %                                 %Save data for ranksum
        %                                 ii_rank=ii_rank+1;
        %                                 mi_rank(ii_rank).mi=mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %                                 mi_rank(ii_rank).per_ii=per_ii;
        %                                 mi_rank(ii_rank).grNo=grNo;
        %                                 mi_rank(ii_rank).evNo=evNo;
        %                                 maxmi=max([maxmi max(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))]);
        %                                 minmi=min([minmi min(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))]);
        %                             end
        %                         end
        %
        %                     end
        %
        %                     title(evTypeLabels{evNo})
        %                     xlabel('MI')
        %                     ylabel('Probability')
        %                 end
        %
        %
        %
        %             for evNo=1:length(eventType)
        %                 subplot(ceil(length(eventType)/2),2,evNo)
        %                 xlim([minmi-0.1*(maxmi-minmi) maxmi+0.1*(maxmi-minmi)])
        %             end
        %
        %             %suptitle(['Average MI for each electrode calculated per  for PAC theta/' freq_names{pacii+1}])
        %
        %             %Now do the ranksums
        %             prof_naive_leg{1}='Proficient';
        %             prof_naive_leg{2}='Naive';
        %
        %             input_data=[];
        %             for ii=1:ii_rank
        %                 input_data(ii).data=mi_rank(ii).mi;
        %                 input_data(ii).description=[handles_drgb.drgbchoices.group_no_names{mi_rank(ii).grNo} ' ' evTypeLabels{mi_rank(ii).evNo} ' ' prof_naive_leg{mi_rank(ii).per_ii}];
        %
        %                 glm_mi.data(glm_ii+1:glm_ii+length(mi_rank(ii).mi))=mi_rank(ii).mi;
        %                 glm_mi.group(glm_ii+1:glm_ii+length(mi_rank(ii).mi))=mi_rank(ii).grNo;
        %                 glm_mi.perCorr(glm_ii+1:glm_ii+length(mi_rank(ii).mi))=mi_rank(ii).per_ii;
        %                 glm_mi.event(glm_ii+1:glm_ii+length(mi_rank(ii).mi))=mi_rank(ii).evNo;
        %                 glm_ii=glm_ii+length(mi_rank(ii).mi);
        %             end
        %
        %             %Perform the glm
        %             fprintf(1, ['glm for average MI for each electrode calculated per mouse for PAC theta' freq_names{pacii+1} '\n'])
        %
        %             if sum(glm_mi.group==1)==length(glm_mi.group)
        %                 %There is only one group here (e.g. for Justin's paper we only include
        %                 %forward)
        %                 fprintf(1, ['\n\nglm for MI for Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_mi.data',glm_mi.perCorr',glm_mi.event',...
        %                     'VariableNames',{'MI','perCorr','event'});
        %                 mdl = fitglm(tbl,'MI~perCorr+event+perCorr*event'...
        %                     ,'CategoricalVars',[2,3])
        %             else
        %
        %                 fprintf(1, ['\n\nglm for MI for Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_mi.data',glm_mi.group',glm_mi.perCorr',glm_mi.event',...
        %                     'VariableNames',{'MI','group','perCorr','event'});
        %                 mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
        %                     ,'CategoricalVars',[2,3,4])
        %             end
        %
        %             %Do the ranksum/t-test
        %             fprintf(1, ['\n\nRanksum or t-test p values for average MI for each electrode calculated per mouse for PAC theta' freq_names{pacii+1} '\n'])
        %             [output_data] = drgMutiRanksumorTtest(input_data);
        %
        %             out_mi_rank(pacii).mi_rank=mi_rank;
        %             end
        % %             for ii=1:ii_rank
        % %                 for jj=ii+1:ii_rank
        % %                     [p, r_or_t]=drg_ranksum_or_ttest(mi_rank(ii).mi,mi_rank(jj).mi);
        % %                     if r_or_t==0
        % %                         fprintf(1, ['p value ranksum for ' handles_drgb.drgbchoices.group_no_names{mi_rank(ii).grNo} ' ' evTypeLabels{mi_rank(ii).evNo} ' ' prof_naive_leg{mi_rank(ii).per_ii} ' vs ' ...
        % %                             handles_drgb.drgbchoices.group_no_names{mi_rank(jj).grNo} ' ' evTypeLabels{mi_rank(jj).evNo} ' ' prof_naive_leg{mi_rank(jj).per_ii} ' =  %d\n'],p)
        % %                     else
        % %                         fprintf(1, ['p value t-test for ' handles_drgb.drgbchoices.group_no_names{mi_rank(ii).grNo} ' ' evTypeLabels{mi_rank(ii).evNo} ' ' prof_naive_leg{mi_rank(ii).per_ii} ' vs ' ...
        % %                             handles_drgb.drgbchoices.group_no_names{mi_rank(jj).grNo} ' ' evTypeLabels{mi_rank(jj).evNo} ' ' prof_naive_leg{mi_rank(jj).per_ii} ' =  %d\n'],p)
        % %                     end
        % %                     pvals=[pvals p];
        % %                 end
        % %             end
        % %             fprintf(1, ['\n\n'])
        %         end
        % %
        % %         pFDR_mi_rank=drsFDRpval( pvals);
        % %         fprintf(1, ['pFDR for per mi per mouse, per electrode  = %d\n\n'],pFDR_mi_rank);
        %
        %         %Now plot the average MI per mouse averaged over electrodes
        %
        %         for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
        %
        %             data_MI_per_mouse=[];
        %             prof_naive_per_mouse=[];
        %             events_per_mouse=[];
        %             groups_per_mouse=[];
        %
        %
        %             glm_perm_ii=0;
        %             glm_mi_perm=[];
        %
        %             %Plot the average
        %             figNo=figNo+3;
        %             try
        %                 close(figNo)
        %             catch
        %             end
        %             hFig=figure(figNo);
        %
        %             set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
        %             hold on
        %
        %             %             bar_lab_loc=[];
        %             no_ev_labels=0;
        %             ii_gr_included=0;
        %
        %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %                 include_group=0;
        %
        %                 for evNo=1:length(eventType)
        %
        %                     for per_ii=1:2
        %
        %                         if sum(eventType==3)>0
        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
        %                         else
        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
        %                         end
        %
        %                         these_offsets(per_ii)=bar_offset;
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             include_group=1;
        %
        %                             %Compute per mouse avearge
        %                             each_mouse_average_MI=[];
        %                             no_mice_included=0;
        %                             for mouseNo=1:max(mean_MI_mouseNo_per_mouse)
        %                                 if mouse_included(mouseNo)==1
        %                                     if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
        %                                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo))>0
        %                                             no_mice_included=no_mice_included+1;
        %                                             each_mouse_average_MI(no_mice_included)=mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo)));
        %                                         end
        %                                     end
        %                                 end
        %                             end
        %                             if no_mice_included>0
        %
        %                                 %Save the per mouse averages
        %                                 glm_mi_perm.data(glm_perm_ii+1:glm_perm_ii+no_mice_included)=each_mouse_average_MI;
        %                                 glm_mi_perm.group(glm_perm_ii+1:glm_perm_ii+no_mice_included)=grNo;
        %                                 glm_mi_perm.perCorr(glm_perm_ii+1:glm_perm_ii+no_mice_included)=per_ii;
        %                                 glm_mi_perm.event(glm_perm_ii+1:glm_perm_ii+no_mice_included)=evNo;
        %                                 glm_perm_ii=glm_perm_ii+no_mice_included;
        %
        %                                 include_group=1;
        %
        %                                 if per_ii==1
        %                                     bar(bar_offset,mean(each_mouse_average_MI),'r','LineWidth', 3,'EdgeColor','none')
        %                                 else
        %                                     bar(bar_offset,mean(each_mouse_average_MI),'b','LineWidth', 3,'EdgeColor','none')
        %                                 end
        %
        %
        %                                 plot(bar_offset,mean(each_mouse_average_MI),'ok','LineWidth', 3)
        %                                 plot((bar_offset)*ones(1,length(each_mouse_average_MI)),each_mouse_average_MI,'o',...
        %                                     'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
        %
        %                                 if length(each_mouse_average_MI)>2
        %                                     CI = bootci(1000, {@mean, each_mouse_average_MI},'type','cper');
        %                                     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        %                                 end
        %
        %                                 %Show the mean in the cumulative histos
        %                                 if sum(group_included)==1
        %
        %                                     figure(figNo-3)
        %                                     hold on
        %
        %
        %                                     for jj=1:length(each_mouse_average_MI)
        %                                         this_f_mi=[];
        %                                         this_f_mi=cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_mi;
        %
        %                                         this_x_mi=[];
        %                                         this_x_mi=cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_mi;
        %
        %                                         xii_below=find(this_x_mi<each_mouse_average_MI(jj),1,'last');
        %                                         xii_above=find(this_x_mi>each_mouse_average_MI(jj),1,'first');
        %
        %                                         slope=(this_f_mi(xii_above)-this_f_mi(xii_below))/(this_x_mi(xii_above)-this_x_mi(xii_below));
        %                                         intercept=this_f_mi(xii_above)-slope*this_x_mi(xii_above);
        %
        %                                         this_f=slope*each_mouse_average_MI(jj)+intercept;
        %                                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I added this
        %                                              switch grNo
        %                                         case 1
        %                                             if per_ii==1
        %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        %                                             else
        %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        %                                             end
        %                                         case 2
        %                                             if per_ii==1
        %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        %                                             else
        %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[1 0.7 .7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
        %                                             end
        %                                         case 3
        %
        %                                             if per_ii==1
        %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[0 1 0],'MarkerEdge',[0 1 0],'MarkerSize',10)
        %                                             else
        %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[0.7 1 0.7],'MarkerEdge',[0.7 1 0.7 ],'MarkerSize',10)
        %                                             end
        %
        % %                                         if evNo==1
        % %                                             if per_ii==1
        % %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        % %                                             else
        % %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        % %                                             end
        % %                                         else
        % %                                             if per_ii==1
        % %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
        % %                                             else
        % %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        % %                                             end
        %                                         end
        %
        %                                     end
        %
        %                                 else
        %                                     figure(figNo-3)
        %                                     subplot(ceil(length(eventType)/2),2,evNo)
        %                                     hold on
        %
        %
        %                                     for jj=1:length(each_mouse_average_MI)
        %                                         this_f_mi=[];
        %                                         this_f_mi=cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_mi;
        %
        %                                         this_x_mi=[];
        %                                         this_x_mi=cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_mi;
        %
        %                                         xii_below=find(this_x_mi<each_mouse_average_MI(jj),1,'last');
        %                                         xii_above=find(this_x_mi>each_mouse_average_MI(jj),1,'first');
        %
        %                                         slope=(this_f_mi(xii_above)-this_f_mi(xii_below))/(this_x_mi(xii_above)-this_x_mi(xii_below));
        %                                         intercept=this_f_mi(xii_above)-slope*this_x_mi(xii_above);
        %
        %                                         this_f=slope*each_mouse_average_MI(jj)+intercept;
        %
        %                                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I added this
        %                                            switch grNo
        %                                 case 1
        %                                     if per_ii==1
        % %                                         legends.pacii(pacii).evNo(evNo).p1=plot(x_mi,f_mi,'Color',[0 0 1],'LineWidth',3);
        % %                                     else
        % %                                         legends.pacii(pacii).evNo(evNo).p2=plot(x_mi,f_mi,'Color',[0.7 0.7 1],'LineWidth',3);
        % %                                     end
        %                                         plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        %                                             else
        %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        %                                             end
        %
        %                                 case 2
        %                                     if per_ii==1
        % %                                         legends.pacii(pacii).evNo(evNo).p3=plot(x_mi,f_mi,'Color',[1 0 0],'LineWidth',3);
        % %                                     else
        % %                                         legends.pacii(pacii).evNo(evNo).p4=plot(x_mi,f_mi,'Color',[1 0.7 0.7],'LineWidth',3);
        % %                                     end
        %                                         plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        %                                             else
        %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
        %                                             end
        %                                 case 3
        %                                     if per_ii==1
        % %                                         legends.pacii(pacii).evNo(evNo).p5=plot(x_mi,f_mi,'Color',[0 1 0],'LineWidth',3);
        % %                                     else
        % %                                         legends.pacii(pacii).evNo(evNo).p6=plot(x_mi,f_mi,'Color',[0.7 1 0.7],'LineWidth',3);
        % %                                     end
        %                                         plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[0 1 0],'MarkerEdge',[0 1 0],'MarkerSize',10)
        %                                             else
        %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[0.7 1 0.7],'MarkerEdge',[0.7 1 0.7],'MarkerSize',10)
        %                                             end
        %
        %
        % %                                         if grNo==1
        % %                                             if per_ii==1
        % %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        % %                                             else
        % %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        % %                                             end
        % %                                         else
        % %                                             if per_ii==1
        % %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
        % %                                             else
        % %                                                 plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        % %                                             end
        %                                         end
        %
        %                                     end
        %                                 end
        %
        %                                 figure(figNo)
        %                                 hold on
        %
        %                                 %                                 %Save data for anovan
        %                                 %                                 data_MI_per_mouse=[data_MI_per_mouse mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))];
        %                                 %                                 prof_naive_per_mouse=[prof_naive_per_mouse per_ii*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))];
        %                                 %                                 events_per_mouse=[events_per_mouse evNo*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))];
        %                                 %                                 groups_per_mouse=[groups_per_mouse grNo*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))];
        %                             end
        %                         end
        %                     end
        %                 end
        %
        %
        %                 if include_group==1
        %                     ii_gr_included=ii_gr_included+1;
        %                     groups_included(ii_gr_included)=grNo;
        %                 end
        %             end
        %
        %             title(['MI per mouse averaged over all electrodes for PAC theta/' freq_names{pacii+1}])
        %
        %
        %             %Annotations identifying groups
        %             x_interval=0.8/ii_gr_included;
        %             for ii=1:ii_gr_included
        %                 annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
        %             end
        %
        %             %Proficient/Naive annotations
        %             annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
        %             annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
        %
        %             %             %x labels
        %             %             to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
        %             %             sorted_A=sortrows(to_sort);
        %             %             sorted_bar_lab_loc=sorted_A(:,1);
        %             %             for ii=1:length(bar_lab_loc)
        %             %                 sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
        %             %             end
        % %             xticks(sorted_bar_lab_loc)
        % %             xticklabels(sorted_bar_labels)
        %
        %             if sum(eventType==3)==0
        %                 xlabel('Concentration (%)')
        %             end
        %
        %             ylabel('Modulation Index')
        %
        %             %             %Calculate anovan for inteaction
        %             %             [p,tbl,stats]=anovan(data_MI_per_mouse,{prof_naive_per_mouse events_per_mouse, groups_per_mouse},'model','interaction','varnames',{'proficient_vs_naive','events','groups'},'display','off');
        %             %             fprintf(1, ['anovan MI per mouse averaged over electrodes for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(1));
        %             %
        %             %             fprintf(1, ['anovan p value for naive vs proficient for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(1));
        %             %             fprintf(1, ['anovan p value for events for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(2));
        %             %             fprintf(1, ['anovan p value for groups for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(3));
        %             %             fprintf(1, ['anovan p value for naive vs proficient * events for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(1));
        %             %             fprintf(1, ['anovan p value for naive vs proficient * groups for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(2));
        %             %             fprintf(1, ['anovan p value for events * groups for PAC theta/' freq_names{pacii+1} '= %d \n\n'],  p(3));
        %             %
        %             %Perform the glm
        %             fprintf(1, ['glm for average MI calculated with per mouse for PAC theta' freq_names{pacii+1} '\n'])
        %
        %             if sum(glm_mi_perm.group==1)==length(glm_mi_perm.group)
        %                 %There is only one group here (e.g. for Justin's paper we only include
        %                 %forward)
        %                 fprintf(1, ['\n\nglm for MI calculated per mouse for Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_mi_perm.data',glm_mi_perm.perCorr',glm_mi_perm.event',...
        %                     'VariableNames',{'MI','perCorr','event'});
        %                 mdl = fitglm(tbl,'MI~perCorr+event+perCorr*event'...
        %                     ,'CategoricalVars',[2,3])
        %             else
        %
        %                 fprintf(1, ['\n\nglm for MI calculated per mouse for Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_mi_perm.data',glm_mi_perm.group',glm_mi_perm.perCorr',glm_mi_perm.event',...
        %                     'VariableNames',{'MI','group','perCorr','event'});
        %                 mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
        %                     ,'CategoricalVars',[2,3,4])
        %             end
        %         end
        %
        %
        %      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I
        %      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Stop
        %      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%here
        %
        %         if sum(group_included)~=1
        %
        %             for pacii=1:no_pacii
        %                 for evNo=1:length(eventType)
        %                     figure(figNo-3)
        %
        %                     subplot(ceil(length(eventType)/2),2,evNo)
        %
        %                     hold on
        %                     try
        %                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I
        %                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%added
        %                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%this
        %                            legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4 legends.pacii(pacii).evNo(evNo).p5 legends.pacii(pacii).evNo(evNo).p6],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
        %                         [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'],...
        %                         [handles_drgb.drgbchoices.group_no_names{3} ' proficient'],[handles_drgb.drgbchoices.group_no_names{3} ' naive'])
        % %                         legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
        % %                             [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
        %                     catch
        %                     end
        %
        %                 end
        %                 suptitle(['Average MI for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1}])
        %             end
        %         else
        % %             for pacii=1:no_pacii
        % %                 figure(figNo-3)
        % %                 try
        % %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I added this
        % %                             legend([legends.pacii(pacii).evNo(1).per_ii(1).p1 legends.pacii(pacii).evNo(2).per_ii(2).p2 legends.pacii(pacii).evNo(3).per_ii(3).p3 legends.pacii(pacii).evNo(4).per_ii(4).p4 legends.pacii(pacii).evNo(5).per_ii(5).p5 legends.pacii(pacii).evNo(6).per_ii(6).p6],[evTypeLabels{1} ' proficient'],[evTypeLabels{1} ' naive'],...
        % %                         [evTypeLabels{2} ' proficient'],[evTypeLabels{2} ' naive'],...
        % %                         [evTypeLabels{3} ' proficient'],[evTypeLabels{3} ' naive'])
        % %
        % % %                     legend([legends.pacii(pacii).evNo(1).per_ii(1).p legends.pacii(pacii).evNo(1).per_ii(2).p legends.pacii(pacii).evNo(2).per_ii(1).p legends.pacii(pacii).evNo(2).per_ii(2).p],[evTypeLabels{1} ' proficient'],[evTypeLabels{1} ' naive'],...
        % % %                         [evTypeLabels{2} ' proficient'],[evTypeLabels{2} ' naive'])
        % %                 catch
        % %                 end
        % %                 title(['Average MI for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1}])
        % %             end
        %         end
        %
        %         %Now do the cumulative histograms and ranksums for meanVectorLength per electrode per mouse
        % %         mean_meanVectorLength_per_mouse
        %
        %         %Compute the meanVectorLength average per mouse
        %         for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
        %
        %
        %             ii_gr_included=0;
        %
        %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %                 include_group=0;
        %
        %                 for evNo=1:length(eventType)
        %
        %                     for per_ii=1:2
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             include_group=1;
        %
        %                             %Compute per mouse avearge
        %                             each_mouse_average_VL=[];
        %                             no_mice_included=0;
        %                             for mouseNo=1:max(mean_MI_mouseNo_per_mouse)
        %                                 if mouse_included(mouseNo)==1
        %                                     if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
        %                                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo))>0
        %                                             no_mice_included=no_mice_included+1;
        %                                             these_VL=[];
        %                                             these_VL=mean_meanVectorLength_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo));
        %                                             each_mouse_average_VL(no_mice_included)=mean(these_VL);
        %                                         end
        %                                     end
        %                                 end
        %                             end
        %                             if (pacii==1)&(evNo==2)
        %                                 pffft=1;
        %                             end
        %                             if no_mice_included>0
        %
        %                                 include_group=1;
        %
        %                                 mouse_avg_VL.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).each_mouse_average_VL=each_mouse_average_VL;
        %
        %
        %                             end
        %                         end
        %                     end
        %                 end
        %                 if include_group==1
        %                     ii_gr_included=ii_gr_included+1;
        %                     groups_included(ii_gr_included)=grNo;
        %                 end
        %             end
        %         end
        %
        %
        %
        %         pvals=[];
        %         legends=[];
        %         cum_histoVL=[];
        %         mean_all_meansVL=[];
        %         max_all_shifted_meansVL=[];
        %         all_means_shifted_meanVLs=[];
        %         all_CIs_shifted_meanVLs=[];
        %
        %         for pacii=1:no_pacii
        %
        %             ii_rank=0;
        %             VL_rank=[];
        %             glm_VL=[];
        %             glm_ii=0;
        %             maxVL=-2000;
        %             minVL=2000;
        %             figNo = figNo + 1;
        %             try
        %                 close(figNo)
        %             catch
        %             end
        %             hFig=figure(figNo);
        %
        %             if length(eventType)>2
        %                 set(hFig, 'units','normalized','position',[.1 .1 .7 .8])
        %             else
        %                 set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
        %             end
        %
        %             %Figure out minVL and maxVL
        %             for evNo=1:length(eventType)
        %                 subplot(ceil(length(eventType)/2),2,evNo)
        %                 hold on
        %
        %                 %Calculate the mean of the mean of each distribution
        %                 all_meansVL=[];
        %                 no_means=0;
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             if (evNo==2)&(pacii==2)
        %                                 pffft=1;
        %                             end
        %                             these_meanVLs=[];
        %                             these_meanVLs=mean_meanVectorLength_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %
        %                             sztmVL=size(these_meanVLs);
        %                             shifted_meanVLs=zeros(sztmVL(1),sztmVL(2));
        %
        %                             this_meanVL=[];
        %                             this_meanVL=mean(these_meanVLs);
        %
        %                             no_means=no_means+1;
        %                             all_meansVL(no_means)=this_meanVL;
        %                         end
        %                     end
        %                 end
        %
        %
        %                 mean_all_meansVL(pacii,evNo)=mean(all_meansVL);
        %
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             these_meanVLs=[];
        %                             these_meanVLs=mean_meanVectorLength_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %
        %                             sztmVL=size(these_meanVLs);
        %                             shifted_meanVLs=zeros(sztmVL(1),sztmVL(2));
        %
        %                             this_meanVL=[];
        %                             this_meanVL=mean(these_meanVLs);
        %
        %                             max_all_meansVL(pacii,evNo,grNo,per_ii)=max(these_meanVLs);
        %                             if length(eventType)>2
        %                                 all_means_meanVLs(pacii,evNo,grNo,per_ii)=mean(these_meanVLs);
        %                                 all_CIs_meanVLs(pacii,evNo,grNo,per_ii,1:2) = bootci(1000, {@mean, these_meanVLs},'type','cper');
        %                             end
        %
        %
        %                             maxVL=max([maxVL max(these_meanVLs)]);
        %                             minVL=min([minVL min(these_meanVLs)]);
        %                         end
        %                     end
        %
        %                 end
        %
        %
        %             end
        %
        %             for evNo=1:length(eventType)
        %                 subplot(ceil(length(eventType)/2),2,evNo)
        %                 hold on
        %
        %                 %Calculate the mean of the mean of each distribution
        %                 all_meansVL=[];
        %                 no_means=0;
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             these_meanVLs=[];
        %                             these_meanVLs=mean_meanVectorLength_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %
        %                             this_meanVL=[];
        %                             this_meanVL=mean(these_meanVLs);
        %
        %                             no_means=no_means+1;
        %                             all_meansVL(no_means)=this_meanVL;
        %                         end
        %                     end
        %                 end
        %
        %
        %                 mean_all_meansVL(pacii,evNo)=mean(all_meansVL);
        %
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             these_meanVLs=[];
        %                             these_meanVLs=mean_meanVectorLength_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %
        %                             this_meanVL=[];
        %                             this_meanVL=mean(these_meanVLs);
        %
        %                             max_all_meansVL(pacii,evNo,grNo,per_ii)=max(these_meanVLs');
        %                             if length(eventType)>2
        %                                 all_means_meanVLs(pacii,evNo,grNo,per_ii)=mean(these_meanVLs');
        %                                 all_CIs_meanVLs(pacii,evNo,grNo,per_ii,1:2) = bootci(1000, {@mean, these_meanVLs'},'type','cper');
        %                             end
        %
        %                             [f_VL,x_VL] = drg_ecdf(these_meanVLs);
        %
        %                             cum_histoVL.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_VL=f_VL;
        %                             cum_histoVL.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_VL=x_VL;
        %                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I added this
        %                               switch grNo
        %                                 case 1
        %                                 if per_ii==1
        %                                     legends.pacii(pacii).evNo(evNo).p1=plot(x_VL,f_VL,'Color',[0 0 1],'LineWidth',3);
        %                                 else
        %                                     legends.pacii(pacii).evNo(evNo).p2=plot(x_VL,f_VL,'Color',[0.7 0.7 1],'LineWidth',3);
        %                                 end
        %                                 case 2
        %                                 if per_ii==1
        %                                     legends.pacii(pacii).evNo(evNo).p3=plot(x_VL,f_VL,'Color',[1 0 0],'LineWidth',3);
        %                                 else
        %                                     legends.pacii(pacii).evNo(evNo).p4=plot(x_VL,f_VL,'Color',[1 0.7 0.7],'LineWidth',3);
        %                                 end
        %                                   case 3
        %                                 if per_ii==1
        %                                     legends.pacii(pacii).evNo(evNo).p5=plot(x_VL,f_VL,'Color',[0 1 0],'LineWidth',3);
        %                                 else
        %                                     legends.pacii(pacii).evNo(evNo).p6=plot(x_VL,f_VL,'Color',[0.7 1 0.7],'LineWidth',3);
        %                                 end
        %                             end
        % %                             if grNo==1
        % %                                 if per_ii==1
        % %                                     legends.pacii(pacii).evNo(evNo).p1=plot(x_VL,f_VL,'Color',[1 0 0],'LineWidth',3);
        % %                                 else
        % %                                     legends.pacii(pacii).evNo(evNo).p2=plot(x_VL,f_VL,'Color',[0 0 1],'LineWidth',3);
        % %                                 end
        % %                             else
        % %                                 if per_ii==1
        % %                                     legends.pacii(pacii).evNo(evNo).p3=plot(x_VL,f_VL,'Color',[1 0.7 0.7],'LineWidth',3);
        % %                                 else
        % %                                     legends.pacii(pacii).evNo(evNo).p4=plot(x_VL,f_VL,'Color',[0.7 0.7 1],'LineWidth',3);
        % %                                 end
        % %                             end
        % %
        %                             %Now plot the average per mouse
        %                             each_mouse_average_VL=[];
        %                             each_mouse_average_VL=mouse_avg_VL.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).each_mouse_average_VL;
        %
        %                             %Show the mean in the cumulative histos
        %                             for jj=1:length(each_mouse_average_VL)
        %                                 this_f_VL=[];
        %                                 this_f_VL=cum_histoVL.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_VL;
        %
        %                                 this_x_VL=[];
        %                                 this_x_VL=cum_histoVL.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_VL;
        %
        %                                 xii_below=find(this_x_VL<each_mouse_average_VL(jj),1,'last');
        %                                 xii_above=find(this_x_VL>each_mouse_average_VL(jj),1,'first');
        %
        %                                 slope=(this_f_VL(xii_above)-this_f_VL(xii_below))/(this_x_VL(xii_above)-this_x_VL(xii_below));
        %                                 intercept=this_f_VL(xii_above)-slope*this_x_VL(xii_above);
        %
        %                                 this_f=slope*each_mouse_average_VL(jj)+intercept;
        %
        %                                 if each_mouse_average_VL(jj)>max(this_x_VL)
        %                                     this_f=1;
        %                                 end
        %
        %                                 if each_mouse_average_VL(jj)<min(this_x_VL)
        %                                     this_f=0;
        %                                 end
        %                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I added this
        %                                     switch grNo
        %                                    case 1
        %                                     if per_ii==1
        %                                         plot(each_mouse_average_VL(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        %                                     else
        %                                         plot(each_mouse_average_VL(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        %                                     end
        %                                    case 2
        %                                     if per_ii==1
        %                                         plot(each_mouse_average_VL(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        %                                     else
        %                                         plot(each_mouse_average_VL(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
        %                                     end
        %                                    case 3
        %                                     if per_ii==1
        %                                         plot(each_mouse_average_VL(jj),this_f,'o','MarkerFace',[0 1 0],'MarkerEdge',[0 1 0],'MarkerSize',10)
        %                                     else
        %                                         plot(each_mouse_average_VL(jj),this_f,'o','MarkerFace',[0.7 1 0.7],'MarkerEdge',[0.7 1 0.7],'MarkerSize',10)
        %                                     end
        %                                 end
        %
        %                             end
        %
        %
        % %                                 if grNo==1
        % %                                     if per_ii==1
        % %                                         plot(each_mouse_average_VL(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        % %                                     else
        % %                                         plot(each_mouse_average_VL(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        % %                                     end
        % %                                 else
        % %                                     if per_ii==1
        % %                                         plot(each_mouse_average_VL(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
        % %                                     else
        % %                                         plot(each_mouse_average_VL(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        % %                                     end
        % %                                 end
        % %
        % %                             end
        %
        %                             %Save data for ranksum
        %                             ii_rank=ii_rank+1;
        %                             VL_rank(ii_rank).meanVL=these_meanVLs;
        %                             VL_rank(ii_rank).per_ii=per_ii;
        %                             VL_rank(ii_rank).grNo=grNo;
        %                             VL_rank(ii_rank).evNo=evNo;
        %                             maxVL=max([maxVL max(these_meanVLs)]);
        %                             minVL=min([minVL min(these_meanVLs)]);
        %                         end
        %                     end
        %
        %                 end
        %
        %                 title(evTypeLabels{evNo})
        %                 xlabel('mean vector length')
        %                 ylabel('Probability')
        %                 try
        %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I added
        %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%this
        %                       legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4 legends.pacii(pacii).evNo(evNo).p5 legends.pacii(pacii).evNo(evNo).p6],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
        %                         [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'],...
        %                         [handles_drgb.drgbchoices.group_no_names{3} ' proficient'],[handles_drgb.drgbchoices.group_no_names{3} ' naive'])
        % %                     legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
        % %                         [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
        %                 catch
        %                 end
        %                  xlim([minVL-0.1*(maxVL-minVL) maxVL+0.1*(maxVL-minVL)])
        %             end
        %
        %             %Now do the ranksum or t-test
        %             prof_naive_leg{1}='Proficient';
        %             prof_naive_leg{2}='Naive';
        %
        %             input_data=[];
        %             for ii=1:ii_rank
        %                 input_data(ii).data=VL_rank(ii).meanVL;
        %                 input_data(ii).description=[handles_drgb.drgbchoices.group_no_names{VL_rank(ii).grNo} ' ' evTypeLabels{VL_rank(ii).evNo} ' ' prof_naive_leg{VL_rank(ii).per_ii}];
        %
        %                 glm_VL.data(glm_ii+1:glm_ii+length(VL_rank(ii).meanVL))=VL_rank(ii).meanVL;
        %                 glm_VL.group(glm_ii+1:glm_ii+length(VL_rank(ii).meanVL))=VL_rank(ii).grNo;
        %                 glm_VL.perCorr(glm_ii+1:glm_ii+length(VL_rank(ii).meanVL))=VL_rank(ii).per_ii;
        %                 glm_VL.event(glm_ii+1:glm_ii+length(VL_rank(ii).meanVL))=VL_rank(ii).evNo;
        %                 glm_ii=glm_ii+length(VL_rank(ii).meanVL);
        %             end
        %
        %             %Perform the glm
        %             fprintf(1, ['glm for mean vector length for each electrode calculated per mouse for PAC theta ' freq_names{pacii+1} '\n'])
        %
        %             if sum(glm_VL.group==1)==length(glm_VL.group)
        %                 %There is only one group here (e.g. for Justin's paper we only include
        %                 %forward)
        %                 fprintf(1, ['\n\nglm for vector length for Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_VL.data',glm_VL.perCorr',glm_VL.event',...
        %                     'VariableNames',{'Vector_length','perCorr','event'});
        %                 mdl = fitglm(tbl,'Vector_length~perCorr+event+perCorr*event'...
        %                     ,'CategoricalVars',[2,3])
        %             else
        %
        %                 fprintf(1, ['\n\nglm for vector length for Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_VL.data',glm_VL.group',glm_VL.perCorr',glm_VL.event',...
        %                     'VariableNames',{'Vector_length','group','perCorr','event'});
        %                 mdl = fitglm(tbl,'Vector_length~group+perCorr+event+perCorr*group*event'...
        %                     ,'CategoricalVars',[2,3,4])
        %             end
        %
        %             fprintf(1, ['\n\nRanksum or t-test for mean vector length for each electrode calculated per mouse for PAC theta ' freq_names{pacii+1} '\n'])
        %             input_data=[];
        %             for ii=1:ii_rank
        %                 input_data(ii).data=VL_rank(ii).meanVL;
        %                 input_data(ii).description=[handles_drgb.drgbchoices.group_no_names{VL_rank(ii).grNo} ' ' evTypeLabels{VL_rank(ii).evNo} ' ' prof_naive_leg{VL_rank(ii).per_ii}];
        %             end
        %             [output_data] = drgMutiRanksumorTtest(input_data);
        %
        %
        %
        %             fprintf(1, ['\n\n'])
        %
        %             suptitle(['Mean vector length per mouse for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1} ])
        %
        %         end
        %
        % %   Now do the cumulative histograms and ranksums for circuar variance of the peak angle per electrode per mouse
        %
        % % (180/pi)*circ_axial(circ_mean(these_meanPAs'*pi/180))'
        
        %         %Compute the peakAngle variance average per mouse
        %     for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
        %
        %
        %
        %             ii_gr_included=0;
        %
        %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %                 include_group=0;
        %
        %                 for evNo=1:length(eventType)
        %
        %                     for per_ii=1:2
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             include_group=1;
        %
        %                             %Compute per mouse avearge
        %                             each_mouse_average_PAvar=[];
        %                             no_mice_included=0;
        %                             for mouseNo=1:max(mean_MI_mouseNo_per_mouse)
        %                                 if mouse_included(mouseNo)==1
        %                                     if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
        %                                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo))>0
        %                                             no_mice_included=no_mice_included+1;
        %                                             these_PAvar=[];
        %                                             these_PAvar=mean_peakAngleVar_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo));
        %                                             each_mouse_average_PAvar(no_mice_included)=mean(these_PAvar);
        %                                         end
        %                                     end
        %                                 end
        %                             end
        %                             if (pacii==1)&(evNo==2)
        %                                 pffft=1;
        %                             end
        %                             if no_mice_included>0
        %
        %                                 include_group=1;
        %
        %                                 mouse_avg_PAvar.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).each_mouse_average_PAvar=each_mouse_average_PAvar;
        %
        %
        %                             end
        %                         end
        %                     end
        %                 end
        %                 if include_group==1
        %                     ii_gr_included=ii_gr_included+1;
        %                     groups_included(ii_gr_included)=grNo;
        %                 end
        %             end
        %         end
        %
        %
        %         figNo=figNo+3;
        %         pvals=[];
        %         legends=[];
        %         cum_histoPAvar=[];
        %         mean_all_meansPAvar=[];
        %         max_all_shifted_meansPAvar=[];
        %         all_means_shifted_meanPAvars=[];
        %         all_CIs_shifted_meanPAvars=[];
        %         out_PAvar_rank=[];
        %
        %         for pacii=1:no_pacii
        %
        %             ii_rank=0;
        %             PAvar_rank=[];
        %             glm_PAvar=[];
        %             glm_ii=0;
        %             glm_perm_ii=0;
        %             glm_PAvar_perm=[];
        %             maxPAvar=-2000;
        %             minPAvar=2000;
        %             figNo = figNo + 1;
        %             try
        %                 close(figNo)
        %             catch
        %             end
        %             hFig=figure(figNo);
        %
        %             if length(eventType)>2
        %                 set(hFig, 'units','normalized','position',[.1 .1 .7 .8])
        %             else
        %                 set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
        %             end
        %
        %             %Figure out minPAvar and maxPAvar
        %             for evNo=1:length(eventType)
        %                 subplot(ceil(length(eventType)/2),2,evNo)
        %                 hold on
        %
        %                 %Calculate the mean of the mean of each distribution
        %                 all_meansPAvar=[];
        %                 no_means=0;
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             if (evNo==2)&(pacii==2)
        %                                 pffft=1;
        %                             end
        %                             these_meanPAvars=[];
        %                             these_meanPAvars=mean_peakAngleVar_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %
        %                             sztmPAvar=size(these_meanPAvars);
        %                             shifted_meanPAvars=zeros(sztmPAvar(1),sztmPAvar(2));
        %
        %                             this_meanPAvar=[];
        %                             this_meanPAvar=mean(these_meanPAvars);
        %
        %                             no_means=no_means+1;
        %                             all_meansPAvar(no_means)=this_meanPAvar;
        %                         end
        %                     end
        %                 end
        %
        %
        %                 mean_all_meansPAvar(pacii,evNo)=mean(all_meansPAvar);
        %
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             these_meanPAvars=[];
        %                             these_meanPAvars=mean_peakAngleVar_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %
        %                             sztmPAvar=size(these_meanPAvars);
        %                             shifted_meanPAvars=zeros(sztmPAvar(1),sztmPAvar(2));
        %
        %                             this_meanPAvar=[];
        %                             this_meanPAvar=mean(these_meanPAvars);
        %
        %                             max_all_meansPAvar(pacii,evNo,grNo,per_ii)=max(these_meanPAvars);
        %                             if length(eventType)>2
        %                                 all_means_meanPAvars(pacii,evNo,grNo,per_ii)=mean(these_meanPAvars);
        %                                 all_CIs_meanPAvars(pacii,evNo,grNo,per_ii,1:2) = bootci(1000, {@mean, these_meanPAvars},'type','cper');
        %                             end
        %
        %
        %                             maxPAvar=max([maxPAvar max(these_meanPAvars)]);
        %                             minPAvar=min([minPAvar min(these_meanPAvars)]);
        %                         end
        %                     end
        %
        %                 end
        %
        %
        %             end
        %
        %             for evNo=1:length(eventType)
        %                 subplot(ceil(length(eventType)/2),2,evNo)
        %                 hold on
        %
        %                 %Calculate the mean of the mean of each distribution
        %                 all_meansPAvar=[];
        %                 no_means=0;
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             these_meanPAvars=[];
        %                             these_meanPAvars=mean_peakAngleVar_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %
        %                             this_meanPAvar=[];
        %                             this_meanPAvar=mean(these_meanPAvars);
        %
        %                             no_means=no_means+1;
        %                             all_meansPAvar(no_means)=this_meanPAvar;
        %                         end
        %                     end
        %                 end
        %
        %
        %                 mean_all_meansPAvar(pacii,evNo)=mean(all_meansPAvar);
        %
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             these_meanPAvars=[];
        %                             these_meanPAvars=mean_peakAngleVar_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %
        %                             this_meanPAvar=[];
        %                             this_meanPAvar=mean(these_meanPAvars);
        %
        %                             max_all_meansPAvar(pacii,evNo,grNo,per_ii)=max(these_meanPAvars');
        %                             if length(eventType)>2
        %                                 all_means_meanPAvars(pacii,evNo,grNo,per_ii)=mean(these_meanPAvars');
        %                                 all_CIs_meanPAvars(pacii,evNo,grNo,per_ii,1:2) = bootci(1000, {@mean, these_meanPAvars'},'type','cper');
        %                             end
        %
        %                             [f_PAvar,x_PAvar] = drg_ecdf(these_meanPAvars);
        %
        %                             cum_histoPAvar.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_PAvar=f_PAvar;
        %                             cum_histoPAvar.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_PAvar=x_PAvar;
        %                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I added this
        %                                    switch grNo
        %                                 case 1
        %                                 if per_ii==1
        %                                     legends.pacii(pacii).evNo(evNo).p1=plot(x_PAvar,f_PAvar,'Color',[0 0 1],'LineWidth',3);
        %                                 else
        %                                     legends.pacii(pacii).evNo(evNo).p2=plot(x_PAvar,f_PAvar,'Color',[0.7 0.7 1],'LineWidth',3);
        %                                 end
        %                                 case 2
        %                                 if per_ii==1
        %                                     legends.pacii(pacii).evNo(evNo).p3=plot(x_PAvar,f_PAvar,'Color',[1 0 0],'LineWidth',3);
        %                                 else
        %                                     legends.pacii(pacii).evNo(evNo).p4=plot(x_PAvar,f_PAvar,'Color',[1 0.7 0.7],'LineWidth',3);
        %                                 end
        %                                 case 3
        %                                 if per_ii==1
        %                                     legends.pacii(pacii).evNo(evNo).p5=plot(x_PAvar,f_PAvar,'Color',[0 1 0],'LineWidth',3);
        %                                 else
        %                                     legends.pacii(pacii).evNo(evNo).p6=plot(x_PAvar,f_PAvar,'Color',[0.7 1 0.7],'LineWidth',3);
        %                                 end
        %                             end
        % %                             if grNo==1
        % %                                 if per_ii==1
        % %                                     legends.pacii(pacii).evNo(evNo).p1=plot(x_PAvar,f_PAvar,'Color',[1 0 0],'LineWidth',3);
        % %                                 else
        % %                                     legends.pacii(pacii).evNo(evNo).p2=plot(x_PAvar,f_PAvar,'Color',[0 0 1],'LineWidth',3);
        % %                                 end
        % %                             else
        % %                                 if per_ii==1
        % %                                     legends.pacii(pacii).evNo(evNo).p3=plot(x_PAvar,f_PAvar,'Color',[1 0.7 0.7],'LineWidth',3);
        % %                                 else
        % %                                     legends.pacii(pacii).evNo(evNo).p4=plot(x_PAvar,f_PAvar,'Color',[0.7 0.7 1],'LineWidth',3);
        % %                                 end
        % %                             end
        %
        %                             %Now plot the average per mouse
        %                             each_mouse_average_PAvar=[];
        %                             each_mouse_average_PAvar=mouse_avg_PAvar.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).each_mouse_average_PAvar;
        %
        %                             %Save the per mouse averages
        %                             glm_PAvar_perm.data(glm_perm_ii+1:glm_perm_ii+length(each_mouse_average_PAvar))=each_mouse_average_PAvar;
        %                             glm_PAvar_perm.group(glm_perm_ii+1:glm_perm_ii+length(each_mouse_average_PAvar))=grNo;
        %                             glm_PAvar_perm.perCorr(glm_perm_ii+1:glm_perm_ii+length(each_mouse_average_PAvar))=per_ii;
        %                             glm_PAvar_perm.event(glm_perm_ii+1:glm_perm_ii+length(each_mouse_average_PAvar))=evNo;
        %                             glm_perm_ii=glm_perm_ii+length(each_mouse_average_PAvar);
        %
        %                             %Show the mean in the cumulative histos
        %                             for jj=1:length(each_mouse_average_PAvar)
        %                                 this_f_PAvar=[];
        %                                 this_f_PAvar=cum_histoPAvar.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_PAvar;
        %
        %                                 this_x_PAvar=[];
        %                                 this_x_PAvar=cum_histoPAvar.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_PAvar;
        %
        %                                 xii_below=find(this_x_PAvar<each_mouse_average_PAvar(jj),1,'last');
        %                                 xii_above=find(this_x_PAvar>each_mouse_average_PAvar(jj),1,'first');
        %
        %                                 slope=(this_f_PAvar(xii_above)-this_f_PAvar(xii_below))/(this_x_PAvar(xii_above)-this_x_PAvar(xii_below));
        %                                 intercept=this_f_PAvar(xii_above)-slope*this_x_PAvar(xii_above);
        %
        %                                 this_f=slope*each_mouse_average_PAvar(jj)+intercept;
        %
        %                                 if each_mouse_average_PAvar(jj)>max(this_x_PAvar)
        %                                     this_f=1;
        %                                 end
        %
        %                                 if each_mouse_average_PAvar(jj)<min(this_x_PAvar)
        %                                     this_f=0;
        %                                 end
        %                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I added this
        %                                        switch grNo
        %                                     case 1
        %                                     if per_ii==1
        %                                         plot(each_mouse_average_PAvar(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        %                                     else
        %                                         plot(each_mouse_average_PAvar(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        %                                     end
        %                                     case 2
        %                                     if per_ii==1
        %                                         plot(each_mouse_average_PAvar(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        %                                     else
        %                                         plot(each_mouse_average_PAvar(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
        %                                     end
        %                                     case 3
        %                                     if per_ii==1
        %                                         plot(each_mouse_average_PAvar(jj),this_f,'o','MarkerFace',[0 1 0],'MarkerEdge',[0 1 0],'MarkerSize',10)
        %                                     else
        %                                         plot(each_mouse_average_PAvar(jj),this_f,'o','MarkerFace',[0.7 1 0.7],'MarkerEdge',[0.7 1 0.7],'MarkerSize',10)
        %                                     end
        %                                 end
        %
        %                             end
        %
        % %                                 if grNo==1
        % %                                     if per_ii==1
        % %                                         plot(each_mouse_average_PAvar(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        % %                                     else
        % %                                         plot(each_mouse_average_PAvar(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        % %                                     end
        % %                                 else
        % %                                     if per_ii==1
        % %                                         plot(each_mouse_average_PAvar(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
        % %                                     else
        % %                                         plot(each_mouse_average_PAvar(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        % %                                     end
        % %                                 end
        % %
        % %                             end
        % %
        %                             %Save data for ranksum
        %                             ii_rank=ii_rank+1;
        %                             PAvar_rank(ii_rank).meanPAvar=these_meanPAvars;
        %                             PAvar_rank(ii_rank).per_ii=per_ii;
        %                             PAvar_rank(ii_rank).grNo=grNo;
        %                             PAvar_rank(ii_rank).evNo=evNo;
        %                             maxPAvar=max([maxPAvar max(these_meanPAvars)]);
        %                             minPAvar=min([minPAvar min(these_meanPAvars)]);
        %                         end
        %                     end
        %
        %                 end
        %
        %                 title(evTypeLabels{evNo})
        %                 xlabel('Peak angle variance (deg^2)')
        %                 ylabel('Probability')
        %                 try
        %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I added this
        %                      legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4 legends.pacii(pacii).evNo(evNo).p5 legends.pacii(pacii).evNo(evNo).p6],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
        %                         [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'],...
        %                     [handles_drgb.drgbchoices.group_no_names{3} ' proficient'],[handles_drgb.drgbchoices.group_no_names{3} ' naive'])
        %
        % %                     legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
        % %                         [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
        %                 catch
        %                 end
        %                  xlim([minPAvar-0.1*(maxPAvar-minPAvar) maxPAvar+0.1*(maxPAvar-minPAvar)])
        %             end
        %
        %             input_data=[];
        %             for ii=1:ii_rank
        %                 input_data(ii).data=PAvar_rank(ii).meanPAvar;
        %                 input_data(ii).description=[handles_drgb.drgbchoices.group_no_names{PAvar_rank(ii).grNo} ' ' evTypeLabels{PAvar_rank(ii).evNo} ' ' prof_naive_leg{PAvar_rank(ii).per_ii}];
        %
        %                 glm_PAvar.data(glm_ii+1:glm_ii+length(PAvar_rank(ii).meanPAvar))=PAvar_rank(ii).meanPAvar;
        %                 glm_PAvar.group(glm_ii+1:glm_ii+length(PAvar_rank(ii).meanPAvar))=PAvar_rank(ii).grNo;
        %                 glm_PAvar.perCorr(glm_ii+1:glm_ii+length(PAvar_rank(ii).meanPAvar))=PAvar_rank(ii).per_ii;
        %                 glm_PAvar.event(glm_ii+1:glm_ii+length(PAvar_rank(ii).meanPAvar))=PAvar_rank(ii).evNo;
        %                 glm_ii=glm_ii+length(PAvar_rank(ii).meanPAvar);
        %             end
        %
        %             %Perform the glm
        %             fprintf(1, ['glm for  peak angle variance for each electrode calculated per mouse for PAC theta ' freq_names{pacii+1} '\n'])
        %
        %             if sum(glm_PAvar.group==1)==length(glm_PAvar.group)
        %                 %There is only one group here (e.g. for Justin's paper we only include
        %                 %forward)
        %                 fprintf(1, ['\n\nglm for peak angle variance for Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_PAvar.data',glm_PAvar.perCorr',glm_PAvar.event',...
        %                     'VariableNames',{'Peak_angle_var','perCorr','event'});
        %                 mdl = fitglm(tbl,'Peak_angle_var~perCorr+event+perCorr*event'...
        %                     ,'CategoricalVars',[2,3])
        %             else
        %
        %                 fprintf(1, ['\n\nglm for peak angle variance for Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_PAvar.data',glm_PAvar.group',glm_PAvar.perCorr',glm_PAvar.event',...
        %                     'VariableNames',{'Peak_angle_var','group','perCorr','event'});
        %                 mdl = fitglm(tbl,'Peak_angle_var~group+perCorr+event+perCorr*group*event'...
        %                     ,'CategoricalVars',[2,3,4])
        %             end
        %
        %             %Now do the ranksum or t-test
        %             fprintf(1, ['\n\nRanksum or t-test for peak angle variance for each electrode calculated per mouse for PAC theta ' freq_names{pacii+1} '\n'])
        %             prof_naive_leg{1}='Proficient';
        %             prof_naive_leg{2}='Naive';
        %
        %
        %             input_data=[];
        %             for ii=1:ii_rank
        %                 input_data(ii).data=PAvar_rank(ii).meanPAvar;
        %                 input_data(ii).description=[handles_drgb.drgbchoices.group_no_names{PAvar_rank(ii).grNo} ' ' evTypeLabels{PAvar_rank(ii).evNo} ' ' prof_naive_leg{PAvar_rank(ii).per_ii}];
        %             end
        %             [output_data] = drgMutiRanksumorTtest(input_data);
        %
        %             out_PAvar_rank(pacii).PAvar_rank=PAvar_rank;
        %
        %             fprintf(1, ['\n\n'])
        %
        %             suptitle(['Peak angle variance for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1} ])
        %
        %             fprintf(1, ['glm for average PA variance calculated with per mouse for PAC theta' freq_names{pacii+1} '\n'])
        %
        %             if sum(glm_PAvar_perm.group==1)==length(glm_PAvar_perm.group)
        %                 %There is only one group here (e.g. for Justin's paper we only include
        %                 %forward)
        %                 fprintf(1, ['\n\nglm for average PA variance calculated with per mouse for Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_PAvar_perm.data',glm_PAvar_perm.perCorr',glm_PAvar_perm.event',...
        %                     'VariableNames',{'MI','perCorr','event'});
        %                 mdl = fitglm(tbl,'MI~perCorr+event+perCorr*event'...
        %                     ,'CategoricalVars',[2,3])
        %             else
        %
        %                 fprintf(1, ['\n\nglm for average PA variance calculated with per mouse for Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_PAvar_perm.data',glm_PAvar_perm.group',glm_PAvar_perm.perCorr',glm_PAvar_perm.event',...
        %                     'VariableNames',{'MI','group','perCorr','event'});
        %                 mdl = fitglm(tbl,'MI~group+perCorr+event+perCorr*group*event'...
        %                     ,'CategoricalVars',[2,3,4])
        %             end
        %         end
        %
        %Now plot the peak angle variance for each electrode calculated per mouse
        %(including all sessions for each mouse)
        edges=[0:100:3500];
        rand_offset=0.8;
        
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            %             data_MI=[];
            %             prof_naive=[];
            %             events=[];
            %             groups=[];
            %             mice=[];
            
            glm_PAvar=[];
            glm_ii=0;
            id_ii=0;
            input_data=[];
            
            %
            %Plot the average
            figNo = figNo +1;
            
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            %             try
            %                 close(figNo+pacii)
            %             catch
            %             end
            %             hFig=figure(figNo+pacii);
            
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            hold on
            
            bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            bar_offset = 0;
            
            %             for grNo=1:max(handles_drgb.drgbchoices.group_no) lo puse
            %             abajo
            
            include_group=0;
            
            for evNo=1:length(eventType)
                
                for per_ii=2:-1:1
                    
                    for grNo=1:max(handles_drgb.drgbchoices.group_no)
                        bar_offset = bar_offset +1;
                        
                        %                         if sum(eventType==3)>0
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
                        %                         else
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        %                         end
                        %
                        %                         these_offsets(per_ii)=bar_offset;
                        bar_offset = bar_offset + 1;
                        
                        if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
                            
                            include_group=1;
                            
                            switch grNo
                                case 1
                                    bar(bar_offset,mean(mean_peakAngleVar_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),'g','LineWidth', 3,'EdgeColor','none')
                                case 2
                                    bar(bar_offset,mean(mean_peakAngleVar_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                                case 3
                                    bar(bar_offset,mean(mean_peakAngleVar_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),'y','LineWidth', 3,'EdgeColor','none')
                            end
                            
                            handles_out.PA_ii=handles_out.PA_ii+1;
                            handles_out.PA_values(handles_out.PA_ii).pacii=pacii;
                            handles_out.PA_values(handles_out.PA_ii).evNo=evNo;
                            handles_out.PA_values(handles_out.PA_ii).per_ii=per_ii;
                            handles_out.PA_values(handles_out.PA_ii).groupNo=grNo;
                            handles_out.PA_values(handles_out.PA_ii).PA_var=mean(mean_peakAngleVar_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)));
                            
                            %Save the mean per mouse
                            these_mice=mean_MI_mouseNo_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
                            handles_out.PA_values(handles_out.PA_ii).noMice=0;
                            for iiMice=min(these_mice):max(these_mice)
                                if sum(these_mice==iiMice)>0
                                    handles_out.PA_values(handles_out.PA_ii).noMice=handles_out.PA_values(handles_out.PA_ii).noMice+1;
                                    handles_out.PA_values(handles_out.PA_ii).mouseNo(handles_out.PA_values(handles_out.PA_ii).noMice)=iiMice;
                                    handles_out.PA_values(handles_out.PA_ii).PA_var_per_mouse(handles_out.PA_values(handles_out.PA_ii).noMice)=mean(mean_peakAngleVar_per_mouse((mean_MI_mouseNo_per_mouse==iiMice)&(~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)));
                                end
                            end
                            
                            x_PAvar=mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
                            glm_PAvar.data(glm_ii+1:glm_ii+length(x_PAvar))=x_PAvar;
                            glm_PAvar.group(glm_ii+1:glm_ii+length(x_PAvar))=grNo*ones(1,length(x_PAvar));
                            glm_PAvar.perCorr(glm_ii+1:glm_ii+length(x_PAvar))=per_ii*ones(1,length(x_PAvar));
                            glm_PAvar.event(glm_ii+1:glm_ii+length(x_PAvar))=evNo*ones(1,length(x_PAvar));
                            glm_ii=glm_ii+length(x_PAvar);
                            
                            id_ii=id_ii+1;
                            input_data(id_ii).data=x_PAvar;
                            input_data(id_ii).description=[group_legend{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                            
                            
                            %Violin plot
                            [mean_out, CIout]=drgViolinPoint(mean_peakAngleVar_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))...
                                ,edges,bar_offset,rand_offset,'k','k',1);
                            
                            %                             if per_ii==1
                            %                                 bar(bar_offset,mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),'r','LineWidth', 3,'EdgeColor','none')
                            %                             else
                            %                                 bar(bar_offset,mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                            %                             end
                            
                            
                            %                             plot(bar_offset,mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
                            %                             plot((bar_offset)*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),...
                            %                                 mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)),'o',...
                            %                                 'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            %
                            %                             if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>=2
                            %                                 CI = bootci(1000, {@mean, mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))},'type','cper');
                            %                                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            %                             end
                            
                            %Save data for anovan
                            %                             data_MI=[data_MI mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))];
                            %                             prof_naive=[prof_naive per_ii*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))];
                            %                             events=[events evNo*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))];
                            %                             groups=[groups grNo*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))];
                            %                             mice=[mice mean_MI_mouseNo_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))];
                            %
                        end
                    end
                    bar_offset = bar_offset + 2;
                    %                     if include_group==1
                    %                         bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                    %                         no_ev_labels=no_ev_labels+1;
                    %                         if sum(eventType==3)>0
                    %                             bar_labels{no_ev_labels}=evTypeLabels{evNo};
                    %                         else
                    %                             bar_labels{no_ev_labels}=num2str(concs(evNo));
                    %                         end
                    %                     end
                end
                bar_offset = bar_offset + 3;
                %                 if include_group==1
                %                     ii_gr_included=ii_gr_included+1;
                %                     groups_included(ii_gr_included)=grNo;
                %                 end
            end
            
            title(['Average peak angle variance for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1}])
            
            %              ylim([0 1.2*max(all_bar)])
            %Annotations identifying groups
            x_interval=0.8/ii_gr_included;
            for ii=1:ii_gr_included
                annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %             %x labels
            %             to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
            %             sorted_A=sortrows(to_sort);
            %             sorted_bar_lab_loc=sorted_A(:,1);
            %             for ii=1:length(bar_lab_loc)
            %                 sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
            %             end
            xticks([2 4 6 10 12 14 21 23 25 29 31 33])
            xticklabels({'nwtS+', 'nHS+', 'nKOS+', 'pwtS+', 'pHS+', 'pKOS+', 'nwtS-', 'nHS-', 'nKOS-', 'pwtS-', 'pHS-', 'pKOS-'})
            
            
            if sum(eventType==3)==0
                xlabel('Concentration (%)')
            end
            
            ylabel('Peak angle variance')
            
            
            %Perform the glm
            fprintf(1, ['glm for peak angle variance for each electrode calculated per mouse for PAC theta' freq_names{pacii+1} '\n'])
            
            fprintf(1, ['\n\nglm for peak angle variance for Theta/' freq_names{pacii+1} '\n'])
            tbl = table(glm_PAvar.data',glm_PAvar.group',glm_PAvar.perCorr',glm_PAvar.event',...
                'VariableNames',{'PAvar','group','perCorr','event'});
            mdl = fitglm(tbl,'PAvar~group+perCorr+event+perCorr*group*event'...
                ,'CategoricalVars',[2,3,4])
            
            
            %Do the ranksum/t-test
            fprintf(1, ['\n\nRanksum or t-test p values for peak angle variance for each electrode calculated per mouse for PAC theta' freq_names{pacii+1} '\n'])
            [output_data] = drgMutiRanksumorTtest(input_data);
            
            
            
        end
        
        pffft=1;
        %Now do the cumulative histograms and ranksums for meanVectorAngle per electrode per mouse
        
        %         %Compute the meanVectorAngle average per mouse
        %         for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
        %
        %
        %             ii_gr_included=0;
        %
        %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %                 include_group=0;
        %
        %                 for evNo=1:length(eventType)
        %
        %                     for per_ii=1:2
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             include_group=1;
        %
        %                             %Compute per mouse avearge
        %                             each_mouse_average_VA=[];
        %                             no_mice_included=0;
        %                             for mouseNo=1:max(mean_MI_mouseNo_per_mouse)
        %                                 if mouse_included(mouseNo)==1
        %                                     if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
        %                                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo))>0
        %                                             no_mice_included=no_mice_included+1;
        %                                             these_VA=[];
        %                                             these_VA=mean_meanVectorAngle_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo));
        %                                             each_mouse_average_VA(no_mice_included)=(180/pi)*circ_axial(circ_mean(these_VA'*pi/180)');
        %                                         end
        %                                     end
        %                                 end
        %                             end
        %                             if (pacii==1)&(evNo==2)
        %                                 pffft=1;
        %                             end
        %                             if no_mice_included>0
        %
        %                                 include_group=1;
        %
        %                                 mouse_avg_VA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).each_mouse_average_VA=each_mouse_average_VA;
        %
        %
        %                             end
        %                         end
        %                     end
        %                 end
        %                 if include_group==1
        %                     ii_gr_included=ii_gr_included+1;
        %                     groups_included(ii_gr_included)=grNo;
        %                 end
        %             end
        %         end
        %
        %
        %         figNo=figNo+3;
        %         pvals=[];
        %         legends=[];
        %         cum_histoVA=[];
        %         mean_all_meansVA=[];
        %         max_all_shifted_meansVA=[];
        %         all_means_shifted_meanVAs=[];
        %         all_CIs_shifted_meanVAs=[];
        %
        %         for pacii=1:no_pacii
        %
        %             ii_rank=0;
        %             VA_rank=[];
        %             maxVA=-2000;
        %             minVA=2000;
        %             figNo = figNo + 1;
        %             try
        %                 close(figNo)
        %             catch
        %             end
        %             hFig=figure(figNo);
        %
        %             if length(eventType)>2
        %                 set(hFig, 'units','normalized','position',[.1 .1 .7 .8])
        %             else
        %                 set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
        %             end
        %
        %             %Figure out minVA and maxVA
        %             for evNo=1:length(eventType)
        %                 subplot(ceil(length(eventType)/2),2,evNo)
        %                 hold on
        %
        %                 %Calculate the mean of the mean of each distribution
        %                 all_meansVA=[];
        %                 no_means=0;
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             %Note: we have to shift the mean to accomodate
        %                             %all cumulative histograms in one angle axis
        %                             if (evNo==2)&(pacii==2)
        %                                 pffft=1;
        %                             end
        %                             these_meanVAs=[];
        %                             these_meanVAs=mean_meanVectorAngle_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %
        %                             sztmVA=size(these_meanVAs);
        %                             shifted_meanVAs=zeros(sztmVA(1),sztmVA(2));
        %
        %                             this_meanVA=[];
        %                             this_meanVA=(180/pi)*circ_axial(circ_mean(these_meanVAs'*pi/180))';
        %
        %                             no_means=no_means+1;
        %                             all_meansVA(no_means)=this_meanVA;
        %                         end
        %                     end
        %                 end
        %
        %
        %                 mean_all_meansVA(pacii,evNo)=(180/pi)*circ_axial(circ_mean(all_meansVA'*pi/180))';
        %
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             %Note: we have to shift the mean to accomodate
        %                             %all cumulative histograms in one angle axis
        %                             if (evNo==2)&(pacii==1)
        %                                 pffft=1;
        %                             end
        %                             these_meanVAs=[];
        %                             these_meanVAs=mean_meanVectorAngle_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %
        %                             sztmVA=size(these_meanVAs);
        %                             shifted_meanVAs=zeros(sztmVA(1),sztmVA(2));
        %
        %                             this_meanVA=[];
        %                             this_meanVA=(180/pi)*circ_axial(circ_mean(these_meanVAs'*pi/180))';
        %
        %                             shifted_meanVAs(these_meanVAs>180+this_meanVA)=-(360-these_meanVAs(these_meanVAs>180+this_meanVA));
        %                             shifted_meanVAs(these_meanVAs<this_meanVA-180)=360+these_meanVAs(these_meanVAs<this_meanVA-180);
        %                             shifted_meanVAs((these_meanVAs<=180+this_meanVA)&(these_meanVAs>=this_meanVA-180))=these_meanVAs((these_meanVAs<=180+this_meanVA)&(these_meanVAs>=this_meanVA-180));
        %
        %                             %Make sure they are all grouped on the same
        %                             %side of the cumulative histogram
        %                             if mean_all_meansVA(pacii,evNo)<180
        %                                 if abs(this_meanVA-mean_all_meansVA(pacii,evNo))>abs(this_meanVA-360+mean_all_meansVA(pacii,evNo))
        %                                     shifted_meanVAs=shifted_meanVAs-360;
        %                                 end
        %                             else
        %                                 if abs(this_meanVA-mean_all_meansVA(pacii,evNo))<abs(this_meanVA-360+mean_all_meansVA(pacii,evNo))
        %                                     shifted_meanVAs=shifted_meanVAs-360;
        %                                 end
        %                             end
        %
        %                             max_all_shifted_meansVA(pacii,evNo,grNo,per_ii)=max(shifted_meanVAs');
        %                             if length(eventType)>2
        %                                 all_means_shifted_meanVAs(pacii,evNo,grNo,per_ii)=mean(shifted_meanVAs');
        %                                 all_CIs_shifted_meanVAs(pacii,evNo,grNo,per_ii,1:2) = bootci(1000, {@mean, shifted_meanVAs'},'type','cper');
        %                             end
        %
        %
        %                             maxVA=max([maxVA max(shifted_meanVAs)]);
        %                             minVA=min([minVA min(shifted_meanVAs)]);
        %                         end
        %                     end
        %
        %                 end
        %
        %
        %             end
        %
        %             for evNo=1:length(eventType)
        %                 subplot(ceil(length(eventType)/2),2,evNo)
        %                 hold on
        %
        %                 %Calculate the mean of the mean of each distribution
        %                 all_meansVA=[];
        %                 no_means=0;
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             %Note: we have to shift the mean to accomodate
        %                             %all cumulative histograms in one angle axis
        %                             if (evNo==2)&(pacii==2)
        %                                 pffft=1;
        %                             end
        %                             these_meanVAs=[];
        %                             these_meanVAs=mean_meanVectorAngle_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %
        %                             sztmVA=size(these_meanVAs);
        %                             shifted_meanVAs=zeros(sztmVA(1),sztmVA(2));
        %
        %                             this_meanVA=[];
        %                             this_meanVA=(180/pi)*circ_axial(circ_mean(these_meanVAs'*pi/180))';
        %
        %                             no_means=no_means+1;
        %                             all_meansVA(no_means)=this_meanVA;
        %                         end
        %                     end
        %                 end
        %
        %
        %                 mean_all_meansVA(pacii,evNo)=(180/pi)*circ_axial(circ_mean(all_meansVA'*pi/180))';
        %
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             %Note: we have to shift the mean to accomodate
        %                             %all cumulative histograms in one angle axis
        %                             if (evNo==2)&(pacii==1)
        %                                 pffft=1;
        %                             end
        %                             these_meanVAs=[];
        %                             these_meanVAs=mean_meanVectorAngle_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %
        %                             sztmVA=size(these_meanVAs);
        %                             shifted_meanVAs=zeros(sztmVA(1),sztmVA(2));
        %
        %                             this_meanVA=[];
        %                             this_meanVA=(180/pi)*circ_axial(circ_mean(these_meanVAs'*pi/180))';
        %
        %                             shifted_meanVAs(these_meanVAs>180+this_meanVA)=-(360-these_meanVAs(these_meanVAs>180+this_meanVA));
        %                             shifted_meanVAs(these_meanVAs<this_meanVA-180)=360+these_meanVAs(these_meanVAs<this_meanVA-180);
        %                             shifted_meanVAs((these_meanVAs<=180+this_meanVA)&(these_meanVAs>=this_meanVA-180))=these_meanVAs((these_meanVAs<=180+this_meanVA)&(these_meanVAs>=this_meanVA-180));
        %
        %                             %Make sure they are all grouped on the same
        %                             %side of the cumulative histogram
        %                             if mean_all_meansVA(pacii,evNo)<180
        %                                 if abs(this_meanVA-mean_all_meansVA(pacii,evNo))>abs(this_meanVA-360+mean_all_meansVA(pacii,evNo))
        %                                     shifted_meanVAs=shifted_meanVAs-360;
        %                                 end
        %                             else
        %                                 if abs(this_meanVA-mean_all_meansVA(pacii,evNo))<abs(this_meanVA-360+mean_all_meansVA(pacii,evNo))
        %                                     shifted_meanVAs=shifted_meanVAs-360;
        %                                 end
        %                             end
        %
        %                             max_all_shifted_meansVA(pacii,evNo,grNo,per_ii)=max(shifted_meanVAs');
        %                             if length(eventType)>2
        %                                 all_means_shifted_meanVAs(pacii,evNo,grNo,per_ii)=mean(shifted_meanVAs');
        %                                 all_CIs_shifted_meanVAs(pacii,evNo,grNo,per_ii,1:2) = bootci(1000, {@mean, shifted_meanVAs'},'type','cper');
        %                             end
        %
        %                             [f_VA,x_VA] = drg_ecdf(shifted_meanVAs);
        %
        %
        %
        %                             cum_histoVA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_VA=f_VA;
        %                             cum_histoVA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_VA=x_VA;
        %
        %                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I added this
        %                             switch grNo
        %                                 case 1
        %                                 if per_ii==1
        %                                     legends.pacii(pacii).evNo(evNo).p1=plot(x_VA,f_VA,'Color',[0 0 1],'LineWidth',3);
        %                                 else
        %                                     legends.pacii(pacii).evNo(evNo).p2=plot(x_VA,f_VA,'Color',[0.7 0.7 1],'LineWidth',3);
        %                                 end
        %                                 case 2
        %                                 if per_ii==1
        %                                     legends.pacii(pacii).evNo(evNo).p3=plot(x_VA,f_VA,'Color',[1 0 0],'LineWidth',3);
        %                                 else
        %                                     legends.pacii(pacii).evNo(evNo).p4=plot(x_VA,f_VA,'Color',[1 0.7 0.7],'LineWidth',3);
        %                                 end
        %                                 case 3
        %                                 if per_ii==1
        %                                     legends.pacii(pacii).evNo(evNo).p3=plot(x_VA,f_VA,'Color',[0 1 0],'LineWidth',3);
        %                                 else
        %                                     legends.pacii(pacii).evNo(evNo).p4=plot(x_VA,f_VA,'Color',[0.7 1 0.7],'LineWidth',3);
        %                                 end
        %                             end
        % %                             if grNo==1
        % %                                 if per_ii==1
        % %                                     legends.pacii(pacii).evNo(evNo).p1=plot(x_VA,f_VA,'Color',[1 0 0],'LineWidth',3);
        % %                                 else
        % %                                     legends.pacii(pacii).evNo(evNo).p2=plot(x_VA,f_VA,'Color',[0 0 1],'LineWidth',3);
        % %                                 end
        % %                             else
        % %                                 if per_ii==1
        % %                                     legends.pacii(pacii).evNo(evNo).p3=plot(x_VA,f_VA,'Color',[1 0.7 0.7],'LineWidth',3);
        % %                                 else
        % %                                     legends.pacii(pacii).evNo(evNo).p4=plot(x_VA,f_VA,'Color',[0.7 0.7 1],'LineWidth',3);
        % %                                 end
        % %                             end
        %
        %                             %Now plot the average per mouse
        %                             each_mouse_average_VA=[];
        %                             each_mouse_average_VA=mouse_avg_VA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).each_mouse_average_VA;
        %
        %                             for no_mice=1:length(each_mouse_average_VA)
        %                                 if each_mouse_average_VA(no_mice)>max_all_shifted_meansVA(pacii,evNo,grNo,per_ii)
        %                                     each_mouse_average_VA(no_mice)=each_mouse_average_VA(no_mice)-360;
        %                                 end
        %                             end
        %
        %
        %                             %Show the mean in the cumulative histos
        %                             for jj=1:length(each_mouse_average_VA)
        %                                 this_f_VA=[];
        %                                 this_f_VA=cum_histoVA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_VA;
        %
        %                                 this_x_VA=[];
        %                                 this_x_VA=cum_histoVA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_VA;
        %
        %                                 xii_below=find(this_x_VA<each_mouse_average_VA(jj),1,'last');
        %                                 xii_above=find(this_x_VA>each_mouse_average_VA(jj),1,'first');
        %
        %                                 slope=(this_f_VA(xii_above)-this_f_VA(xii_below))/(this_x_VA(xii_above)-this_x_VA(xii_below));
        %                                 intercept=this_f_VA(xii_above)-slope*this_x_VA(xii_above);
        %
        %                                 this_f=slope*each_mouse_average_VA(jj)+intercept;
        %
        %                                 if each_mouse_average_VA(jj)>max(this_x_VA)
        %                                     this_f=1;
        %                                 end
        %
        %                                 if each_mouse_average_VA(jj)<min(this_x_VA)
        %                                     this_f=0;
        %                                 end
        %
        %                                        switch grNo
        %                                     case 1
        %                                     if per_ii==1
        %                                         plot(each_mouse_average_VA(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        %                                     else
        %                                         plot(each_mouse_average_VA(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        %                                     end
        %                                     case 2
        %                                     if per_ii==1
        %                                         plot(each_mouse_average_VA(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        %                                     else
        %                                         plot(each_mouse_average_VA(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
        %                                     end
        %                                     case 3
        %                                     if per_ii==1
        %                                         plot(each_mouse_average_VA(jj),this_f,'o','MarkerFace',[0 1 0],'MarkerEdge',[0 1 0],'MarkerSize',10)
        %                                     else
        %                                         plot(each_mouse_average_VA(jj),this_f,'o','MarkerFace',[0.7 1  0.7],'MarkerEdge',[0.7 1 0.7],'MarkerSize',10)
        %                                     end
        %                                 end
        %
        %                             end
        %
        % %                                 if grNo==1
        % %                                     if per_ii==1
        % %                                         plot(each_mouse_average_VA(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        % %                                     else
        % %                                         plot(each_mouse_average_VA(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        % %                                     end
        % %                                 else
        % %                                     if per_ii==1
        % %                                         plot(each_mouse_average_VA(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
        % %                                     else
        % %                                         plot(each_mouse_average_VA(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        % %                                     end
        % %                                 end
        % %
        % %                             end
        %
        %                             %Save data for ranksum
        %                             ii_rank=ii_rank+1;
        %                             VA_rank(ii_rank).meanVA=shifted_meanVAs;
        %                             VA_rank(ii_rank).per_ii=per_ii;
        %                             VA_rank(ii_rank).grNo=grNo;
        %                             VA_rank(ii_rank).evNo=evNo;
        %                             maxVA=max([maxVA max(shifted_meanVAs)]);
        %                             minVA=min([minVA min(shifted_meanVAs)]);
        %                         end
        %                     end
        %
        %                 end
        %
        %                 title(evTypeLabels{evNo})
        %                 xlabel('mean vector angle')
        %                 ylabel('Probability')
        %                 try
        %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I
        %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%added
        %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%this
        %                         legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4 legends.pacii(pacii).evNo(evNo).p5 legends.pacii(pacii).evNo(evNo).p6],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
        %                         [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'],...
        %                     [handles_drgb.drgbchoices.group_no_names{3} ' proficient'],[handles_drgb.drgbchoices.group_no_names{3} ' naive'])
        % %                     legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
        % %                         [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
        %                 catch
        %                 end
        %                  xlim([minVA-0.1*(maxVA-minVA) maxVA+0.1*(maxVA-minVA)])
        %             end
        %
        %             %Now do the Watson-Williams test
        %             fprintf(1, ['\n\nWatson-Williams test p values for mean vector angle for each electrode calculated per mouse for PAC theta ' freq_names{pacii+1} '\n'])
        %             prof_naive_leg{1}='Proficient';
        %             prof_naive_leg{2}='Naive';
        %
        %             ww_test_out.ii_out=0;
        %             pvals=[];
        %             for ii=1:ii_rank
        %                 for jj=ii+1:ii_rank
        %                     %p=ranksum(VA_rank(ii).meanVA,VA_rank(jj).meanVA);
        %                     p=circ_wwtest(pi*VA_rank(ii).meanVA/180,pi*VA_rank(jj).meanVA/180);
        % %                     fprintf(1, ['p value for ' handles_drgb.drgbchoices.group_no_names{VA_rank(ii).grNo} ' ' evTypeLabels{VA_rank(ii).evNo} ' ' prof_naive_leg{VA_rank(ii).per_ii} ' vs ' ...
        % %                         handles_drgb.drgbchoices.group_no_names{VA_rank(jj).grNo} ' ' evTypeLabels{VA_rank(jj).evNo} ' ' prof_naive_leg{VA_rank(jj).per_ii} ' =  %d\n'],p)
        %                     pvals=[pvals p];
        %                     ww_test_out.ii_out=ww_test_out.ii_out+1;
        %                     ww_test_out.p(ww_test_out.ii_out)=p;
        %                     ww_test_out.label{ww_test_out.ii_out}=[handles_drgb.drgbchoices.group_no_names{VA_rank(ii).grNo} ' ' evTypeLabels{VA_rank(ii).evNo} ' ' prof_naive_leg{VA_rank(ii).per_ii} ' vs ' ...
        %                         handles_drgb.drgbchoices.group_no_names{VA_rank(jj).grNo} ' ' evTypeLabels{VA_rank(jj).evNo} ' ' prof_naive_leg{VA_rank(jj).per_ii}];
        %                 end
        %             end
        %
        %             pFDR_VA_rank=drsFDRpval(pvals);
        %
        %             %Now sort the data
        %             these_ii_pairs=[1:ww_test_out.ii_out];
        %             to_sort=[pvals' these_ii_pairs'];
        %             sorted_rows=sortrows(to_sort);
        %             ww_test_out.sorted_ii_pairs=sorted_rows(:,2);
        %
        %             first_above=1;
        %             for ii_out=1:ww_test_out.ii_out
        %                 if first_above==1
        %                     if ww_test_out.p(ww_test_out.sorted_ii_pairs(ii_out))>pFDR_VA_rank
        %                        first_above=0;
        %                        fprintf(1, ['\n\npFDR for per mean vector angle for each electrode calculated per mouse  = %d\n\n'],pFDR_VA_rank);
        %                     end
        %                 end
        %                 this_label=ww_test_out.label(ww_test_out.sorted_ii_pairs(ii_out));
        %                 fprintf(1, ['p value for ' this_label{1} ' =  %d\n'],ww_test_out.p(ww_test_out.sorted_ii_pairs(ii_out)));
        %             end
        %
        %             fprintf(1, ['\n\n'])
        %
        %             suptitle(['Mean vector angle per mouse for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1} ])
        %
        %         end
        %
        %
        %
        %         %Now do the cumulative histograms and ranksums for peak angle per electrode per mouse
        %         pvals=[];
        %         legends=[];
        %         cum_histoPA=[];
        %         mean_all_meansPA=[];
        %         max_all_shifted_meansPA=[];
        %         all_means_shifted_meanPAs=[];
        %         all_CIs_shifted_meanPAs=[];
        %         figNo=figNo+3;
        %
        %         for pacii=1:no_pacii
        %
        %             ii_rank=0;
        %             PA_rank=[];
        %             glm_PA=[];
        %             glm_ii=0;
        %             maxPA=-2000;
        %             minPA=2000;
        %
        %             try
        %                 close(figNo)
        %             catch
        %             end
        %             hFig=figure(figNo);
        %
        %             if length(eventType)>2
        %                 set(hFig, 'units','normalized','position',[.1 .1 .7 .8])
        %             else
        %                 set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
        %             end
        %
        %
        %             for evNo=1:length(eventType)
        %                 subplot(ceil(length(eventType)/2),2,evNo)
        %                 hold on
        %
        %                 %Calculate the mean of the mean of each distribution
        %                 all_meansPA=[];
        %                 no_means=0;
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             %Note: we have to shift the mean to accomodate
        %                             %all cumulative histograms in one angle axis
        %                             if (evNo==2)&(pacii==2)
        %                                 pffft=1;
        %                             end
        %                             these_meanPAs=[];
        %                             these_meanPAs=mean_peakAngle_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %
        %                             sztmPA=size(these_meanPAs);
        %                             shifted_meanPAs=zeros(sztmPA(1),sztmPA(2));
        %
        %                             this_meanPA=[];
        %                             this_meanPA=(180/pi)*circ_axial(circ_mean(these_meanPAs'*pi/180))';
        %
        %                             no_means=no_means+1;
        %                             all_meansPA(no_means)=this_meanPA;
        %                         end
        %                     end
        %                 end
        %
        %
        %                 mean_all_meansPA(pacii,evNo)=(180/pi)*circ_axial(circ_mean(all_meansPA'*pi/180))';
        %
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             %Note: we have to shift the mean to accomodate
        %                             %all cumulative histograms in one angle axis
        %                             if (evNo==2)&(pacii==1)
        %                                 pffft=1;
        %                             end
        %                             these_meanPAs=[];
        %                             these_meanPAs=mean_peakAngle_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
        %
        %                             sztmPA=size(these_meanPAs);
        %                             shifted_meanPAs=zeros(sztmPA(1),sztmPA(2));
        %
        %                             this_meanPA=[];
        %                             this_meanPA=(180/pi)*circ_axial(circ_mean(these_meanPAs'*pi/180))';
        %
        %                             shifted_meanPAs(these_meanPAs>180+this_meanPA)=-(360-these_meanPAs(these_meanPAs>180+this_meanPA));
        %                             shifted_meanPAs(these_meanPAs<this_meanPA-180)=360+these_meanPAs(these_meanPAs<this_meanPA-180);
        %                             shifted_meanPAs((these_meanPAs<=180+this_meanPA)&(these_meanPAs>=this_meanPA-180))=these_meanPAs((these_meanPAs<=180+this_meanPA)&(these_meanPAs>=this_meanPA-180));
        %
        %                             %Make sure they are all grouped on the same
        %                             %side of the cumulative histogram
        %                             if mean_all_meansPA(pacii,evNo)<180
        %                                 if abs(this_meanPA-mean_all_meansPA(pacii,evNo))>abs(this_meanPA-360+mean_all_meansPA(pacii,evNo))
        %                                     shifted_meanPAs=shifted_meanPAs-360;
        %                                 end
        %                             else
        %                                 if abs(this_meanPA-mean_all_meansPA(pacii,evNo))<abs(this_meanPA-360+mean_all_meansPA(pacii,evNo))
        %                                     shifted_meanPAs=shifted_meanPAs-360;
        %                                 end
        %                             end
        %
        %                             max_all_shifted_meansPA(pacii,evNo,grNo,per_ii)=max(shifted_meanPAs');
        %                             if length(eventType)>2
        %                                 all_means_shifted_meanPAs(pacii,evNo,grNo,per_ii)=mean(shifted_meanPAs');
        %                                 all_CIs_shifted_meanPAs(pacii,evNo,grNo,per_ii,1:2) = bootci(1000, {@mean, shifted_meanPAs'},'type','cper');
        %                             end
        %
        %                             [f_PA,x_PA] = drg_ecdf(shifted_meanPAs);
        %
        %
        %
        %                             cum_histoPA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_PA=f_PA;
        %                             cum_histoPA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_PA=x_PA;
        %
        %                                     switch grNo
        %                                 case 1
        %                                 if per_ii==1
        %                                     legends.pacii(pacii).evNo(evNo).p1=plot(x_PA,f_PA,'Color',[0 0 1],'LineWidth',3);
        %                                 else
        %                                     legends.pacii(pacii).evNo(evNo).p2=plot(x_PA,f_PA,'Color',[0.7 0.7 1],'LineWidth',3);
        %                                 end
        %                                 case 2
        %                                 if per_ii==1
        %                                     legends.pacii(pacii).evNo(evNo).p3=plot(x_PA,f_PA,'Color',[1 0 0],'LineWidth',3);
        %                                 else
        %                                     legends.pacii(pacii).evNo(evNo).p4=plot(x_PA,f_PA,'Color',[1 0.7 0.7],'LineWidth',3);
        %                                 end
        %                                  case 3
        %                                 if per_ii==1
        %                                     legends.pacii(pacii).evNo(evNo).p5=plot(x_PA,f_PA,'Color',[0 1 0],'LineWidth',3);
        %                                 else
        %                                     legends.pacii(pacii).evNo(evNo).p6=plot(x_PA,f_PA,'Color',[0.7 1 0.7],'LineWidth',3);
        %                                 end
        %                             end
        % %                             if grNo==1
        % %                                 if per_ii==1
        % %                                     legends.pacii(pacii).evNo(evNo).p1=plot(x_PA,f_PA,'Color',[1 0 0],'LineWidth',3);
        % %                                 else
        % %                                     legends.pacii(pacii).evNo(evNo).p2=plot(x_PA,f_PA,'Color',[0 0 1],'LineWidth',3);
        % %                                 end
        % %                             else
        % %                                 if per_ii==1
        % %                                     legends.pacii(pacii).evNo(evNo).p3=plot(x_PA,f_PA,'Color',[1 0.7 0.7],'LineWidth',3);
        % %                                 else
        % %                                     legends.pacii(pacii).evNo(evNo).p4=plot(x_PA,f_PA,'Color',[0.7 0.7 1],'LineWidth',3);
        % %                                 end
        % %                             end
        %
        %
        %                             %Save data for ranksum
        %                             ii_rank=ii_rank+1;
        %                             PA_rank(ii_rank).meanPA=shifted_meanPAs;
        %                             PA_rank(ii_rank).per_ii=per_ii;
        %                             PA_rank(ii_rank).grNo=grNo;
        %                             PA_rank(ii_rank).evNo=evNo;
        %                             maxPA=max([maxPA max(shifted_meanPAs)]);
        %                             minPA=min([minPA min(shifted_meanPAs)]);
        %                         end
        %                     end
        %
        %                 end
        %
        %                 title(evTypeLabels{evNo})
        %                 xlabel('Peak angle')
        %                 ylabel('Probability')
        %             end
        %
        %
        %             for evNo=1:length(eventType)
        %                 subplot(ceil(length(eventType)/2),2,evNo)
        %                 xlim([minPA-0.1*(maxPA-minPA) maxPA+0.1*(maxPA-minPA)])
        %             end
        %
        %             %suptitle(['Mean vector angle per mouse for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1} ])
        %
        %             %Do a glm
        %
        %             for ii=1:ii_rank
        %                 glm_PA.data(glm_ii+1:glm_ii+length(PA_rank(ii).meanPA))=PA_rank(ii).meanPA;
        %                 glm_PA.group(glm_ii+1:glm_ii+length(PA_rank(ii).meanPA))=PA_rank(ii).grNo;
        %                 glm_PA.perCorr(glm_ii+1:glm_ii+length(PA_rank(ii).meanPA))=PA_rank(ii).per_ii;
        %                 glm_PA.event(glm_ii+1:glm_ii+length(PA_rank(ii).meanPA))=PA_rank(ii).evNo;
        %                 glm_ii=glm_ii+length(PA_rank(ii).meanPA);
        %             end
        %
        % %             %Perform the glm
        % %             fprintf(1, ['glm forpeak angle for each electrode calculated per mouse for PAC theta ' freq_names{pacii+1} '\n'])
        % %
        % %             if sum(glm_PA.group==1)==length(glm_PA.group)
        % %                 %There is only one group here (e.g. for Justin's paper we only include
        % %                 %forward)
        % %                 fprintf(1, ['\n\nglm for peak angle for Theta/' freq_names{pacii+1} '\n'])
        % %                 tbl = table(glm_PA.data',glm_PA.perCorr',glm_PA.event',...
        % %                     'VariableNames',{'Peak_angle_var','perCorr','event'});
        % %                 mdl = fitglm(tbl,'Peak_angle_var~perCorr+event+perCorr*event'...
        % %                     ,'CategoricalVars',[2,3])
        % %             else
        % %
        % %                 fprintf(1, ['\n\nglm for peak angle for Theta/' freq_names{pacii+1} '\n'])
        % %                 tbl = table(glm_PA.data',glm_PA.group',glm_PA.perCorr',glm_PA.event',...
        % %                     'VariableNames',{'Peak_angle_var','group','perCorr','event'});
        % %                 mdl = fitglm(tbl,'Peak_angle_var~group+perCorr+event+perCorr*group*event'...
        % %                     ,'CategoricalVars',[2,3,4])
        % %             end
        %
        %             %Now do the Watson-Williams test
        %             fprintf(1, ['Watson-Williams test p values for peak angle for each electrode calculated per mouse for PAC theta ' freq_names{pacii+1} '\n'])
        %             prof_naive_leg{1}='Proficient';
        %             prof_naive_leg{2}='Naive';
        %
        %
        %
        %             ww_test_out.ii_out=0;
        %             pvals=[];
        %             for ii=1:ii_rank
        %                 for jj=ii+1:ii_rank
        %                     %p=ranksum(PA_rank(ii).meanPA,PA_rank(jj).meanPA);
        %                     p=circ_wwtest(pi*PA_rank(ii).meanPA/180,pi*PA_rank(jj).meanPA/180);
        % %                     fprintf(1, ['p value for ' handles_drgb.drgbchoices.group_no_names{PA_rank(ii).grNo} ' ' evTypeLabels{PA_rank(ii).evNo} ' ' prof_naive_leg{PA_rank(ii).per_ii} ' vs ' ...
        % %                         handles_drgb.drgbchoices.group_no_names{PA_rank(jj).grNo} ' ' evTypeLabels{PA_rank(jj).evNo} ' ' prof_naive_leg{PA_rank(jj).per_ii} ' =  %d\n'],p)
        %                     pvals=[pvals p];
        %                     ww_test_out.ii_out=ww_test_out.ii_out+1;
        %                     ww_test_out.p(ww_test_out.ii_out)=p;
        %                     ww_test_out.label{ww_test_out.ii_out}=[handles_drgb.drgbchoices.group_no_names{PA_rank(ii).grNo} ' ' evTypeLabels{PA_rank(ii).evNo} ' ' prof_naive_leg{PA_rank(ii).per_ii} ' vs ' ...
        %                         handles_drgb.drgbchoices.group_no_names{PA_rank(jj).grNo} ' ' evTypeLabels{PA_rank(jj).evNo} ' ' prof_naive_leg{PA_rank(jj).per_ii}];
        %                 end
        %             end
        %
        %             pFDR_PA_rank=drsFDRpval(pvals);
        %
        %             %Now sort the data
        %             these_ii_pairs=[1:ww_test_out.ii_out];
        %             to_sort=[pvals' these_ii_pairs'];
        %             sorted_rows=sortrows(to_sort);
        %             ww_test_out.sorted_ii_pairs=sorted_rows(:,2);
        %
        %             first_above=1;
        %             for ii_out=1:ww_test_out.ii_out
        %                 if first_above==1
        %                     if ww_test_out.p(ww_test_out.sorted_ii_pairs(ii_out))>pFDR_PA_rank
        %                        first_above=0;
        %                        fprintf(1, ['\n\npFDR for peak angle for each electrode calculated per mouse  = %d\n\n'],pFDR_PA_rank);
        %                     end
        %                 end
        %                 this_label=ww_test_out.label(ww_test_out.sorted_ii_pairs(ii_out));
        %                 fprintf(1, ['p value for ' this_label{1} ' =  %d\n'],ww_test_out.p(ww_test_out.sorted_ii_pairs(ii_out)));
        %             end
        %
        %             fprintf(1, ['\n\n'])
        %
        %
        %
        %         end
        %
        %         %Now plot the average PA per mouse averaged over electrodes
        %         for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
        %
        %
        %             ii_gr_included=0;
        %
        %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %                 include_group=0;
        %
        %                 for evNo=1:length(eventType)
        %
        %                     for per_ii=1:2
        %
        %
        %                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
        %
        %                             include_group=1;
        %
        %                             %Compute per mouse avearge
        %                             each_mouse_average_PA=[];
        %                             no_mice_included=0;
        %                             for mouseNo=1:max(mean_MI_mouseNo_per_mouse)
        %                                 if mouse_included(mouseNo)==1
        %                                     if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
        %                                         if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo))>0
        %                                             no_mice_included=no_mice_included+1;
        %                                             these_PA=[];
        %                                             these_PA=mean_peakAngle_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo));
        %                                             each_mouse_average_PA(no_mice_included)=(180/pi)*circ_axial(circ_mean(these_PA'*pi/180)');
        %                                         end
        %                                     end
        %                                 end
        %                             end
        %                             if (pacii==1)&(evNo==2)
        %                                 pffft=1;
        %                             end
        %                             if no_mice_included>0
        %
        %                                 include_group=1;
        %
        %                                 for no_mice=1:no_mice_included
        %                                     if each_mouse_average_PA(no_mice)>max_all_shifted_meansPA(pacii,evNo,grNo,per_ii)
        %                                         each_mouse_average_PA(no_mice)=each_mouse_average_PA(no_mice)-360;
        %                                     end
        %                                 end
        %
        %                                 %Show the mean in the cumulative histos
        %                                 figure(figNo)
        %                                 subplot(ceil(length(eventType)/2),2,evNo)
        %                                 hold on
        %
        %
        %                                 for jj=1:length(each_mouse_average_PA)
        %                                     this_f_PA=[];
        %                                     this_f_PA=cum_histoPA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_PA;
        %
        %                                     this_x_PA=[];
        %                                     this_x_PA=cum_histoPA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_PA;
        %
        %                                     xii_below=find(this_x_PA<each_mouse_average_PA(jj),1,'last');
        %                                     xii_above=find(this_x_PA>each_mouse_average_PA(jj),1,'first');
        %
        %                                     slope=(this_f_PA(xii_above)-this_f_PA(xii_below))/(this_x_PA(xii_above)-this_x_PA(xii_below));
        %                                     intercept=this_f_PA(xii_above)-slope*this_x_PA(xii_above);
        %
        %                                     this_f=slope*each_mouse_average_PA(jj)+intercept;
        %
        %                                     if each_mouse_average_PA(jj)>max(this_x_PA)
        %                                         this_f=1;
        %                                     end
        %
        %                                     if each_mouse_average_PA(jj)<min(this_x_PA)
        %                                         this_f=0;
        %                                     end
        %                                     %%%%%%%%%%%%%%%%%%%%%%%%%I added this
        %                                                switch grNo
        %                                         case 1
        %                                         if per_ii==1
        %                                             plot(each_mouse_average_PA(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        %                                         else
        %                                             plot(each_mouse_average_PA(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        %                                         end
        %                                         case 2
        %                                         if per_ii==1
        %                                             plot(each_mouse_average_PA(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        %                                         else
        %                                             plot(each_mouse_average_PA(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
        %                                         end
        %                                          case 3
        %                                         if per_ii==1
        %                                             plot(each_mouse_average_PA(jj),this_f,'o','MarkerFace',[0 1 0],'MarkerEdge',[0 1 0],'MarkerSize',10)
        %                                         else
        %                                             plot(each_mouse_average_PA(jj),this_f,'o','MarkerFace',[0.7 1 0.7],'MarkerEdge',[0.7 1 0.7],'MarkerSize',10)
        %                                         end
        %                                     end
        %
        %                                 end
        %
        %
        %                             end
        %                         end
        %                     end
        %                 end
        %
        % %                                     if grNo==1
        % %                                         if per_ii==1
        % %                                             plot(each_mouse_average_PA(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        % %                                         else
        % %                                             plot(each_mouse_average_PA(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        % %                                         end
        % %                                     else
        % %                                         if per_ii==1
        % %                                             plot(each_mouse_average_PA(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
        % %                                         else
        % %                                             plot(each_mouse_average_PA(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        % %                                         end
        % %                                     end
        % %
        % %                                 end
        % %
        % %
        % %                             end
        % %                         end
        % %                     end
        % %                 end
        %                 if include_group==1
        %                     ii_gr_included=ii_gr_included+1;
        %                     groups_included(ii_gr_included)=grNo;
        %                 end
        %             end
        %         end
        %
        %
        %         for pacii=1:no_pacii
        %             %             figure(pacii+3)
        %             %             for evNo=1:length(eventType)
        %             %                 subplot(1,2,evNo)
        %             %                 hold on
        %             %                 try
        %             %                 legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
        %             %                     [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
        %             %                 catch
        %             %                 end
        %             %
        %             %             end
        %             %             suptitle(['Average per mouse MI for each electrode calculated per for PAC theta/' freq_names{pacii+1}])
        %
        %             figure(figNo)
        %             for evNo=1:length(eventType)
        %                 subplot(ceil(length(eventType)/2),2,evNo)
        %                 hold on
        %                 try
        %
        %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I dded this
        %                      legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4 legends.pacii(pacii).evNo(evNo).p5 legends.pacii(pacii).evNo(evNo).p6],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
        %                         [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'],...
        %                     [handles_drgb.drgbchoices.group_no_names{3} ' proficient'],[handles_drgb.drgbchoices.group_no_names{3} ' naive'])
        %
        % % %                     legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
        % % %                         [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
        %                 catch
        %                 end
        %
        %             end
        %             suptitle(['Peak angle per mouse for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1} ])
        %         end
        %
        %
        %         %If this is for concentration plot as a function of concentration
        %         if length(eventType)>2
        %
        %
        %             for pacii=1:no_pacii
        %                 try
        %                     close(figNo)
        %                 catch
        %                 end
        %                 hFig=figure(figNo);
        %
        %                 set(hFig, 'units','normalized','position',[.1 .1 .4 .4])
        %                 hold on
        %
        %
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %                     for per_ii=1:2
        %
        %                         these_means=zeros(length(eventType),1);
        %                         these_means(:,1)=all_means_shifted_meanVAs(pacii,:,grNo,per_ii);
        %
        %                         these_CIs=zeros(length(eventType),2);
        %                         these_CIs(:,:)=all_CIs_shifted_meanVAs(pacii,:,grNo,per_ii,1:2);
        %
        %                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I added this
        %                             switch grNo
        %                             case 1
        %                             if per_ii==1
        %                                 p1=plot(concs,these_means,'-o','Color',[0 0 1], 'MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        %                             else
        %                                 p2=plot(concs,these_means,'-o','Color',[0.7 0.7 1],'MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        %                             end
        %                             case 2
        %                             if per_ii==1
        %                                 p3=plot(concs,these_means,'-o','Color',[1 0 0],'MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        %                             else
        %                                 p4=plot(concs,these_means,'-o','Color',[1 0 0],'MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        %                             end
        %                              case 3
        %                             if per_ii==1
        %                                 p5=plot(concs,these_means,'-o','Color',[0 1 0],'MarkerFace',[0 1 0],'MarkerEdge',[0 1 0],'MarkerSize',10)
        %                             else
        %                                 p6=plot(concs,these_means,'-o','Color',[0.7  1 0.7],'MarkerFace',[0.7 1 0.7],'MarkerEdge',[0.7 1 0.7],'MarkerSize',10)
        %                             end
        %                         end
        %
        % %                         if grNo==1
        % %                             if per_ii==1
        % %                                 p1=plot(concs,these_means,'-o','Color',[1 0 0], 'MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        % %                             else
        % %                                 p2=plot(concs,these_means,'-o','Color',[0 0 1],'MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        % %                             end
        % %                         else
        % %                             if per_ii==1
        % %                                 p3=plot(concs,these_means,'-o','Color',[1 0.7 0.7],'MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
        % %                             else
        % %                                 p4=plot(concs,these_means,'-o','Color',[0.7 0.7 1],'MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        % %                             end
        % %                         end
        %
        %                         for evTy=1:length(eventType)
        %                             plot([concs(evTy) concs(evTy)],these_CIs(evTy,:),'-k')
        %                         end
        %
        %
        %                     end
        %                 end
        %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I added this
        %                  legend([p1 p2 p3 p4 p5 p6],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
        %                     [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'],...
        %                 [handles_drgb.drgbchoices.group_no_names{3} ' proficient'],[handles_drgb.drgbchoices.group_no_names{3} ' naive'])
        % %                 legend([p1 p2 p3 p4],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
        % %                     [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
        %                 xlabel('Percent odor in mineral oil')
        %                 ylabel('Phase (deg)')
        %                 title(['Mean vector angle per mouse for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1} ])
        %
        %             end
        %         end
        pfft=1;
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) handles_pars.output_suffix],'handles_out')
        
        pffft=1;
        
        
    case 20
        % Multiclass ROC analysis of LFP power differences for naive and proficient
        % mice for different epochs (concentrations or S+ vs. S-) and different groups. Analyzed per mouse
        handles_out=[];
        no_dBs=1;
        delta_dB_power=[];
        no_ROCs=0;
        ROCout=[];
        p_vals_ROC=[];
        per_mouse_no_ROCs=0;
        per_mouse_ROCout=[];
        per_mouse_p_vals_ROC=[];
        delta_dB_powerEv1=[];
        no_Ev1=0;
        for evNo=1:length(eventType)
            evNo_out(evNo).noWB=0;
        end
        delta_dB_powerEv1WB=[];
        delta_dB_powerEv2WB=[];
        delta_dB_No_per_mouse=0;
        delta_dB_per_mouse=[];
        delta_dB_perii_per_mouse=[];
        delta_dB_evNo_per_mouse=[];
        delta_dB_bwii_per_mouse=[];
        delta_dB_mouseNo_per_mouse=[];
        delta_dB_electrode_per_mouse=[];
        mouse_included=[];
        
        
        delta_dB_WB_No_per_mouse=0;
        delta_dB_WB_per_mouse=[];
        delta_dB_WB_perii_per_mouse=[];
        delta_dB_WB_evNo_per_mouse=[];
        delta_dB_WB_mouseNo_per_mouse=[];
        delta_dB_WB_electrode_per_mouse=[];
        delta_dB_WB_group_no_per_mouse=[];
        
        
        
        fprintf(1, ['Pairwise auROC analysis for Fig 1 of Justin''s paper\n\n'])
        p_vals=[];
        no_files=max(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        no_ROCs=0;
        szpc=size(percent_windows);
        for per_ii=1:szpc(1)
            
            for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
                mouse_has_files=0;
                theseEvNosPerEl=[];
                for evNo=1:length(eventType)
                    for bwii=1:no_bandwidths
                        for elec=1:16
                            theseEvNosPerEl(evNo,bwii,elec).noEv=0;
                        end
                    end
                end
                
                for elec=1:16
                    if sum(which_electrodes==elec)>0
                        
                        theseEvNos=[];
                        for evNo=1:length(eventType)
                            for bwii=1:no_bandwidths
                                theseEvNos(evNo,bwii).noEv=0;
                            end
                        end
                        
                        theseEvNosWB=[];
                        for evNo=1:length(eventType)
                            theseEvNosWB(evNo).noWB=0;
                        end
                        
                        for fileNo=1:no_files
                            if sum(files==fileNo)>0
                                if handles_drgb.drgbchoices.mouse_no(fileNo)==mouseNo
                                    
                                    
                                    
                                    lfpodNo_ref=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                                    
                                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref)))
                                        
                                        
                                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower))
                                            
                                            percent_mask=[];
                                            percent_mask=logical((handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower>=percent_windows(per_ii,1))...
                                                &(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower<=percent_windows(per_ii,2)));
                                            
                                            
                                            for evNo=1:length(eventType)
                                                
                                                noWB_for_evNo(evNo)=-1;
                                                
                                                trials_in_event_Ev=[];
                                                trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(eventType(evNo),:)==1)&percent_mask;
                                                
                                                if (sum(trials_in_event_Ev)>=1)
                                                    
                                                    lfpodNo=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                                    
                                                    % Ev1
                                                    this_dB_powerref=zeros(sum(trials_in_event_Ev),length(frequency));
                                                    this_dB_powerref(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower(trials_in_event_Ev,:));
                                                    
                                                    
                                                    this_dB_power=zeros(sum(trials_in_event_Ev),length(frequency));
                                                    this_dB_power(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_event_Ev,:));
                                                    
                                                    %Wide band spectrum
                                                    this_delta_dB_powerEvWB=zeros(sum(trials_in_event_Ev),length(frequency));
                                                    this_delta_dB_powerEvWB=this_dB_power-this_dB_powerref;
                                                    theseEvNosWB(evNo).delta_dB_powerEvWB(theseEvNosWB(evNo).noWB+1:theseEvNosWB(evNo).noWB+sum(trials_in_event_Ev),1:length(frequency))=this_delta_dB_powerEvWB;
                                                    dB_powerEvWB=zeros(sum(trials_in_event_Ev),length(frequency));
                                                    dB_powerEvWB=this_dB_power;
                                                    theseEvNosWB(evNo).dB_powerEvWB(theseEvNosWB(evNo).noWB+1:theseEvNosWB(evNo).noWB+sum(trials_in_event_Ev),1:length(frequency))=dB_powerEvWB;
                                                    theseEvNosWB(evNo).groupNo=handles_drgb.drgbchoices.group_no(fileNo)*ones(1,sum(trials_in_event_Ev));
                                                    theseEvNosWB(evNo).noWB=theseEvNosWB(evNo).noWB+sum(trials_in_event_Ev);
                                                    %                                                     noWB_for_evNo(evNo)=evNo_out(evNo).noWB;
                                                    
                                                    
                                                    %Do per bandwidth analysis
                                                    for bwii=1:no_bandwidths
                                                        
                                                        theseEvNos(evNo,bwii).groupNo(theseEvNos(evNo,bwii).noEv+1:theseEvNos(evNo,bwii).noEv+sum(trials_in_event_Ev))=handles_drgb.drgbchoices.group_no(fileNo)*ones(1,sum(trials_in_event_Ev));
                                                        
                                                        this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                                        
                                                        %Enter the  Ev1
                                                        this_delta_dB_powerEv=zeros(sum(trials_in_event_Ev),1);
                                                        this_delta_dB_powerEv=mean(this_dB_power(:,this_band)-this_dB_powerref(:,this_band),2);
                                                        
                                                        theseEvNos(evNo,bwii).this_delta_dB_powerEv(1,theseEvNos(evNo,bwii).noEv+1:theseEvNos(evNo,bwii).noEv+sum(trials_in_event_Ev))=this_delta_dB_powerEv';
                                                        theseEvNos(evNo,bwii).noEv=theseEvNos(evNo,bwii).noEv+sum(trials_in_event_Ev);
                                                        
                                                        if theseEvNos(evNo,bwii).noEv~=length(theseEvNos(evNo,bwii).this_delta_dB_powerEv)
                                                            pffft=1;
                                                        end
                                                        
                                                        theseEvNosPerEl(evNo,bwii,elec).this_delta_dB_powerEv(1,theseEvNosPerEl(evNo,bwii,elec).noEv+1:theseEvNosPerEl(evNo,bwii,elec).noEv+sum(trials_in_event_Ev))=this_delta_dB_powerEv';
                                                        theseEvNosPerEl(evNo,bwii,elec).groupNo(1,theseEvNosPerEl(evNo,bwii,elec).noEv+1:theseEvNosPerEl(evNo,bwii,elec).noEv+sum(trials_in_event_Ev))=handles_drgb.drgbchoices.group_no(fileNo)*ones(1,sum(trials_in_event_Ev));
                                                        
                                                        theseEvNosPerEl(evNo,bwii,elec).noEv=theseEvNosPerEl(evNo,bwii,elec).noEv+sum(trials_in_event_Ev);
                                                        
                                                        %                                                         evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).noWB,bwii)=mean(this_delta_dB_powerEv);
                                                        
                                                        mouse_has_files=1;
                                                    end
                                                    
                                                    
                                                    fprintf(1, ['%d trials in event No %d succesfully processed for file No %d electrode %d\n'],sum(trials_in_event_Ev), evNo,fileNo,elec);
                                                    
                                                    try
                                                        close(1)
                                                    catch
                                                    end
                                                    figure(1)
                                                    hold on
                                                    plot(this_delta_dB_powerEvWB')
                                                    plot(mean(this_delta_dB_powerEvWB)','-k','LineWidth',3)
                                                    plot([4 95],[mean(this_delta_dB_powerEv) mean(this_delta_dB_powerEv)],'-b','LineWidth',2)
                                                    
                                                    pffft=1;
                                                else
                                                    
                                                    
                                                    fprintf(1, ['%d trials in event No %d fewer than minimum trials per event ' evTypeLabels{evNo} ' for file No %d electrode %d\n'],sum(trials_in_event_Ev), min_trials_per_event,fileNo,elec);
                                                    
                                                    
                                                end
                                                
                                                
                                            end
                                            
                                            
                                        else
                                            
                                            fprintf(1, ['Empty allPower for file No %d electrode %d\n'],fileNo,elec);
                                            
                                        end
                                        
                                        
                                    else
                                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],fileNo,elec);
                                        
                                        
                                    end
                                end %if mouseNo
                            end
                        end %fileNo
                        
                        if mouse_has_files==1
                            mouse_included(mouseNo)=1;
                            %Calculate per mouse per electrode delta_dB
                            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                                for evNo=1:length(eventType)
                                    if theseEvNos(evNo).noEv>0
                                        if sum(theseEvNosWB(evNo).groupNo==grNo)>0
                                            %Wide band spectrum
                                            delta_dB_WB_No_per_mouse=delta_dB_WB_No_per_mouse+1;
                                            delta_dB_WB_per_mouse(delta_dB_WB_No_per_mouse,:)=mean(theseEvNosWB(evNo).delta_dB_powerEvWB(theseEvNosWB(evNo).groupNo==grNo,:),1);
                                            dB_WB_per_mouse(delta_dB_WB_No_per_mouse,:)=mean(theseEvNosWB(evNo).dB_powerEvWB(theseEvNosWB(evNo).groupNo==grNo,:),1);
                                            delta_dB_WB_perii_per_mouse(delta_dB_WB_No_per_mouse)=per_ii;
                                            delta_dB_WB_evNo_per_mouse(delta_dB_WB_No_per_mouse)=evNo;
                                            delta_dB_WB_mouseNo_per_mouse(delta_dB_WB_No_per_mouse)=mouseNo;
                                            delta_dB_WB_electrode_per_mouse(delta_dB_WB_No_per_mouse)=elec;
                                            delta_dB_WB_group_no_per_mouse(delta_dB_WB_No_per_mouse)=grNo;
                                        end
                                        for bwii=1:no_bandwidths
                                            if sum(theseEvNos(evNo,1).groupNo==grNo)>0
                                                delta_dB_No_per_mouse=delta_dB_No_per_mouse+1;
                                                delta_dB_per_mouse(delta_dB_No_per_mouse)=mean(theseEvNos(evNo,bwii).this_delta_dB_powerEv(1,logical(theseEvNos(evNo,bwii).groupNo==grNo)));
                                                delta_dB_perii_per_mouse(delta_dB_No_per_mouse)=per_ii;
                                                delta_dB_evNo_per_mouse(delta_dB_No_per_mouse)=evNo;
                                                delta_dB_bwii_per_mouse(delta_dB_No_per_mouse)=bwii;
                                                delta_dB_mouseNo_per_mouse(delta_dB_No_per_mouse)=mouseNo;
                                                delta_dB_electrode_per_mouse(delta_dB_No_per_mouse)=elec;
                                                delta_dB_group_no_per_mouse(delta_dB_No_per_mouse)=grNo;
                                            end
                                        end
                                    end
                                end
                            end
                            
                            %Calculate per electrode ROC
                            can_calculate_auroc=1;
                            if can_calculate_auroc==1
                                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                                    for evNo1=1:length(eventType)
                                        if theseEvNos(evNo1).noEv>0
                                            for evNo2=evNo1+1:length(eventType)
                                                if theseEvNos(evNo2).noEv>0
                                                    for bwii=1:no_bandwidths
                                                        
                                                        
                                                        %Enter Ev1
                                                        trials_in_event_Ev1=length(theseEvNos(evNo1,bwii).this_delta_dB_powerEv(theseEvNos(evNo1).groupNo==grNo));
                                                        this_delta_dB_powerEv1=zeros(trials_in_event_Ev1,1);
                                                        this_delta_dB_powerEv1=theseEvNos(evNo1,bwii).this_delta_dB_powerEv(theseEvNos(evNo1).groupNo==grNo);
                                                        roc_data=[];
                                                        roc_data(1:trials_in_event_Ev1,1)=this_delta_dB_powerEv1;
                                                        roc_data(1:trials_in_event_Ev1,2)=zeros(trials_in_event_Ev1,1);
                                                        
                                                        %Enter Ev2
                                                        trials_in_event_Ev2=length(theseEvNos(evNo2,bwii).this_delta_dB_powerEv(theseEvNos(evNo2).groupNo==grNo));
                                                        total_trials=trials_in_event_Ev1+trials_in_event_Ev2;
                                                        this_delta_dB_powerEv2=zeros(trials_in_event_Ev2,1);
                                                        this_delta_dB_powerEv2=theseEvNos(evNo2,bwii).this_delta_dB_powerEv(theseEvNos(evNo2).groupNo==grNo);
                                                        roc_data(trials_in_event_Ev1+1:total_trials,1)=this_delta_dB_powerEv2;
                                                        roc_data(trials_in_event_Ev1+1:total_trials,2)=ones(trials_in_event_Ev2,1);
                                                        
                                                        
                                                        %Find  per electrode ROC
                                                        if (trials_in_event_Ev1>=5)&(trials_in_event_Ev2>=5)
                                                            no_ROCs=no_ROCs+1;
                                                            roc=roc_calc(roc_data,0,0.05,0);
                                                            ROCout(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo_ref).fileNo;
                                                            ROCelec(no_ROCs)=elec;
                                                            ROCgroups(no_ROCs)=grNo;
                                                            ROCmouse(no_ROCs)=mouseNo;
                                                            ROCbandwidth(no_ROCs)=bwii;
                                                            ROCper_ii(no_ROCs)=per_ii;
                                                            ROCEvNo1(no_ROCs)=evNo1;
                                                            ROCEvNo2(no_ROCs)=evNo2;
                                                            if sum(eventType==3)==0
                                                                if ((abs(evNo1-2)<=1)&(abs(evNo2-5)<=1))||((abs(evNo1-5)<=1)&(abs(evNo2-2)<=1))
                                                                    ROC_between(no_ROCs)=1;
                                                                else
                                                                    ROC_between(no_ROCs)=0;
                                                                end
                                                                ROC_neighbor(no_ROCs)=abs(evNo1-evNo2);
                                                            else
                                                                %This is S+/S-,
                                                                %these values are
                                                                %assigned
                                                                %arbitrarily so
                                                                %that plotting
                                                                %auROC works
                                                                ROC_between(no_ROCs)=1;
                                                                ROC_neighbor(no_ROCs)=2;
                                                            end
                                                            
                                                            auROC(no_ROCs)=roc.AUC-0.5;
                                                            p_valROC(no_ROCs)=roc.p;
                                                            p_vals_ROC=[p_vals_ROC roc.p];
                                                            
                                                            %I have this code here to plot the ROC
                                                            if (per_ii==1)&(bwii==4)&(roc.AUC-0.5>0.3)
                                                                show_roc=0;
                                                                if show_roc==1
                                                                    %I have this code here to plot the ROC
                                                                    roc=roc_calc(roc_data,0,0.05,1);
                                                                    
                                                                    %Do the histograms
                                                                    try
                                                                        close(2)
                                                                    catch
                                                                    end
                                                                    figure(2)
                                                                    
                                                                    hold on
                                                                    
                                                                    max_dB=max([max(this_delta_dB_powerEv1) max(this_delta_dB_powerEv2)]);
                                                                    min_dB=min([min(this_delta_dB_powerEv1) min(this_delta_dB_powerEv2)]);
                                                                    
                                                                    edges=[min_dB-0.1*(max_dB-min_dB):(max_dB-min_dB)/20:max_dB+0.1*(max_dB-min_dB)];
                                                                    histogram(this_delta_dB_powerEv1,edges,'FaceColor','b','EdgeColor','b')
                                                                    histogram(this_delta_dB_powerEv2,edges,'FaceColor','r','EdgeColor','r')
                                                                    xlabel('delta power dB')
                                                                    title(['Histogram for conentrations ' num2str(concs2(evNo1)) ' and ' num2str(concs2(evNo2))])
                                                                    pffft=1;
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            
                        else
                            mouse_included(mouseNo)=0;
                        end
                        
                        
                    end
                    
                end
                
                %Calculate per mouse ROC
                if mouse_has_files==1
                    for grNo=1:max(handles_drgb.drgbchoices.group_no)
                        for evNo1=1:length(eventType)
                            for evNo2=evNo1+1:length(eventType)
                                
                                
                                for bwii=1:no_bandwidths
                                    
                                    %Enter Ev1
                                    trials_in_event_Ev1=length(theseEvNosPerEl(evNo1,bwii,which_electrodes(1)).this_delta_dB_powerEv(theseEvNosPerEl(evNo1,bwii,which_electrodes(1)).groupNo==grNo));
                                    
                                    this_delta_dB_powerEv1=zeros(trials_in_event_Ev1,1);
                                    
                                    %                                     for elec=which_electrodes
                                    %                                         this_delta_dB_powerEv1=this_delta_dB_powerEv1+(theseEvNosPerEl(evNo1,bwii,elec).this_delta_dB_powerEv((theseEvNosPerEl(evNo1,bwii,elec).groupNo==grNo))')/length(which_electrodes);
                                    %                                     end
                                    %
                                    no_elects=0;
                                    for elec=which_electrodes
                                        if length(this_delta_dB_powerEv1)==sum(theseEvNosPerEl(evNo1,bwii,elec).groupNo==grNo)
                                            this_delta_dB_powerEv1=this_delta_dB_powerEv1+(theseEvNosPerEl(evNo1,bwii,elec).this_delta_dB_powerEv((theseEvNosPerEl(evNo1,bwii,elec).groupNo==grNo))');
                                            no_elects=no_elects+1;
                                        end
                                    end
                                    this_delta_dB_powerEv1=this_delta_dB_powerEv1/no_elects;
                                    
                                    roc_data=[];
                                    roc_data(1:trials_in_event_Ev1,1)=this_delta_dB_powerEv1;
                                    roc_data(1:trials_in_event_Ev1,2)=zeros(trials_in_event_Ev1,1);
                                    
                                    %Enter Ev2
                                    trials_in_event_Ev2=length(theseEvNosPerEl(evNo2,bwii,which_electrodes(1)).this_delta_dB_powerEv(theseEvNosPerEl(evNo2,bwii,which_electrodes(1)).groupNo==grNo));
                                    total_trials=trials_in_event_Ev1+trials_in_event_Ev2;
                                    this_delta_dB_powerEv2=zeros(trials_in_event_Ev2,1);
                                    
                                    no_elects=0;
                                    for elec=which_electrodes
                                        if length(this_delta_dB_powerEv2)==sum((theseEvNosPerEl(evNo2,bwii,elec).groupNo==grNo))
                                            this_delta_dB_powerEv2=this_delta_dB_powerEv2+(theseEvNosPerEl(evNo2,bwii,elec).this_delta_dB_powerEv((theseEvNosPerEl(evNo2,bwii,elec).groupNo==grNo))');
                                            no_elects=no_elects+1;
                                        end
                                    end
                                    this_delta_dB_powerEv2=this_delta_dB_powerEv2/no_elects;
                                    
                                    roc_data(trials_in_event_Ev1+1:total_trials,1)=this_delta_dB_powerEv2;
                                    roc_data(trials_in_event_Ev1+1:total_trials,2)=ones(trials_in_event_Ev2,1);
                                    
                                    
                                    %Find  per electrode ROC
                                    if (trials_in_event_Ev1>=5)&(trials_in_event_Ev2>=5)
                                        per_mouse_no_ROCs=per_mouse_no_ROCs+1;
                                        try
                                            roc=roc_calc(roc_data,0,0.05,0);
                                        catch
                                            pffft=1
                                        end
                                        per_mouse_ROCout(per_mouse_no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo_ref).fileNo;
                                        per_mouse_ROCbandwidth(per_mouse_no_ROCs)=bwii;
                                        per_mouse_ROCper_ii(per_mouse_no_ROCs)=per_ii;
                                        per_mouse_ROCgroup(per_mouse_no_ROCs)=grNo;
                                        per_mouse_ROCEvNo1(per_mouse_no_ROCs)=evNo1;
                                        per_mouse_ROCEvNo2(per_mouse_no_ROCs)=evNo2;
                                        if ((abs(evNo1-2)<=1)&(abs(evNo2-5)<=1))||((abs(evNo1-5)<=1)&(abs(evNo2-2)<=1))
                                            per_mouse_ROC_between(per_mouse_no_ROCs)=1;
                                        else
                                            per_mouse_ROC_between(per_mouse_no_ROCs)=0;
                                        end
                                        per_mouse_ROC_neighbor(per_mouse_no_ROCs)=abs(evNo1-evNo2);
                                        per_mouse_auROC(per_mouse_no_ROCs)=roc.AUC-0.5;
                                        per_mouse_p_valROC(per_mouse_no_ROCs)=roc.p;
                                        per_mouse_p_vals_ROC=[p_vals_ROC roc.p];
                                        
                                        %I have this code here to plot the ROC
                                        if roc.AUC-0.5>0.3
                                            show_roc=0;
                                            if show_roc==1
                                                %I have this code here to plot the ROC
                                                roc=roc_calc(roc_data,0,0.05,1);
                                                
                                                %Do the histograms
                                                try
                                                    close(2)
                                                catch
                                                end
                                                figure(2)
                                                
                                                hold on
                                                
                                                max_dB=max([max(this_delta_dB_powerEv1) max(this_delta_dB_powerEv2)]);
                                                min_dB=min([min(this_delta_dB_powerEv1) min(this_delta_dB_powerEv2)]);
                                                
                                                edges=[min_dB-0.1*(max_dB-min_dB):(max_dB-min_dB)/20:max_dB+0.1*(max_dB-min_dB)];
                                                histogram(this_delta_dB_powerEv1,edges,'FaceColor','b','EdgeColor','b')
                                                histogram(this_delta_dB_powerEv2,edges,'FaceColor','r','EdgeColor','r')
                                                xlabel('delta power dB')
                                                title(['Histogram for conentrations ' num2str(concs2(evNo1)) ' and ' num2str(concs2(evNo2))])
                                                pffft=1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                if mouse_has_files==1
                    mouse_included(mouseNo)=1;
                else
                    mouse_included(mouseNo)=0;
                end
            end
            
            
            
        end
        fprintf(1, '\n\n')
        
        
        
        %Now plot the per mouse delta LFP power computed for each electrode
        pvals=[];
        for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            
            %Plot the average
            try
                close(bwii)
            catch
            end
            hFig=figure(bwii);
            
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            data_delta_dB=[];
            prof_naive=[];
            events=[];
            mice=[];
            electrodes=[];
            groups=[];
            
            bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            
            
            %             fprintf(1, ['ANOVAN for delta dB power per mouse per electrode, mouse as random factor\n\n'])
            %
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                
                for evNo=1:length(eventType)
                    
                    per_included=0;
                    these_dB_per_e=[];
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        %bar_offset=21-evNo*3+(2-per_ii);
                        if sum(eventType==3)>0
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(evNo-1);
                        else
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        end
                        
                        these_offsets(per_ii)=bar_offset;
                        
                        if sum((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo))>0
                            
                            include_group=1;
                            
                            per_included=per_included+1;
                            
                            if per_ii==1
                                bar(bar_offset,mean(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo))),'r','LineWidth', 3)
                            else
                                bar(bar_offset,mean(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo))),'b','LineWidth', 3)
                            end
                            
                            %Individual points; in the future add lines linking the
                            %points?
                            plot((bar_offset)*ones(1,sum((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo))),...
                                delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo)),'o',...
                                'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            data_for_lines(per_ii).these_dB_per_e(1:length(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo))))=...
                                delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo));
                            data_for_lines(per_ii).these_mice(1:length(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo))))=...
                                delta_dB_mouseNo_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo));
                            data_for_lines(per_ii).these_bar_offsets=bar_offset;
                            
                            
                            %Average and CI
                            plot(bar_offset,mean(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
                            if sum((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo))>2
                                CI = bootci(1000, {@mean, delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo))},'type','cper');
                                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            end
                            
                            %Save data for anovan
                            data_delta_dB=[data_delta_dB delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo))];
                            prof_naive=[prof_naive per_ii*ones(1,sum((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo)))];
                            events=[events evNo*ones(1,sum((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo)))];
                            mice=[mice delta_dB_mouseNo_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo))];
                            electrodes=[electrodes delta_dB_electrode_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo))];
                            groups=[groups delta_dB_group_no_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo))];
                        end
                    end
                    if per_included==2
                        for mouseNo=1:length(mouse_included)
                            if mouse_included(mouseNo)==1
                                if (sum(data_for_lines(1).these_mice==mouseNo)>0)&(sum(data_for_lines(2).these_mice==mouseNo)>0)
                                    try
                                        plot([data_for_lines(1).these_bar_offsets*ones(1,sum(data_for_lines(1).these_mice==mouseNo)); data_for_lines(2).these_bar_offsets*ones(1,sum(data_for_lines(1).these_mice==mouseNo))],...
                                            [data_for_lines(1).these_dB_per_e(data_for_lines(1).these_mice==mouseNo); data_for_lines(2).these_dB_per_e(data_for_lines(2).these_mice==mouseNo)],'-','Color',[0.7 0.7 0.7])
                                    catch
                                    end
                                end
                            end
                        end
                    end
                    if include_group==1
                        bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                        no_ev_labels=no_ev_labels+1;
                        if sum(eventType==3)>0
                            bar_labels{no_ev_labels}=evTypeLabels{evNo};
                        else
                            bar_labels{no_ev_labels}=num2str(concs(evNo));
                        end
                    end
                end
                if include_group==1
                    ii_gr_included=ii_gr_included+1;
                    groups_included(ii_gr_included)=grNo;
                end
            end
            
            title([freq_names{bwii} ' average delta dB power per mouse, per electrode'])
            
            %Annotations identifying groups
            x_interval=0.8/ii_gr_included;
            for ii=1:ii_gr_included
                annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.8 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %x labels
            to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
            sorted_A=sortrows(to_sort);
            sorted_bar_lab_loc=sorted_A(:,1);
            for ii=1:length(bar_lab_loc)
                sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
            end
            xticks(sorted_bar_lab_loc)
            xticklabels(sorted_bar_labels)
            
            if sum(eventType==3)==0
                xlabel('Concentration (%)')
            end
            
            ylabel('Delta power (dB)')
            
            
            
            %             %Calculate anovan for inteaction
            %             [p,tbl,stats]=anovan(data_delta_dB,{prof_naive events mice electrodes},'varnames',{'proficient_vs_naive','events','groups','mice'},'display','off','random',4);
            %             fprintf(1, ['p value for anovan delta dB power per mouse per electrode for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
            %             fprintf(1, ['p value for anovan delta dB power  per mouse per electrode for events ' freq_names{bwii} '= %d \n'],  p(2));
            %             fprintf(1, ['p value for anovan delta dB power  per mouse per electrode for groups for ' freq_names{bwii} '= %d \n\n'],  p(3));
            %
            
            %Plot the cumulative histos and do ranksum
            %Plot the average
            try
                close(bwii+4)
            catch
            end
            hFig=figure(bwii+4);
            
            set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            prof_naive_leg{1}='Proficient';
            prof_naive_leg{2}='Naive';
            ii_rank=0;
            dBpower_rank=[];
            glm_dBpower_rank=[];
            glm_ii=0;
            maxdB=-200;
            mindB=200;
            for evNo=1:length(eventType)
                subplot(length(eventType),1,evNo)
                hold on
                
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        
                        
                        if sum((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo))>0
                            
                            
                            if per_ii==1
                                if grNo==1
                                    [f_aic,x_aic] = drg_ecdf(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo)));
                                    plot(x_aic,f_aic,'Color',[1 0 0],'LineWidth',3)
                                else
                                    [f_aic,x_aic] = drg_ecdf(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo)));
                                    plot(x_aic,f_aic,'Color',[0 0 1],'LineWidth',3)
                                end
                            else
                                if grNo==1
                                    [f_aic,x_aic] = drg_ecdf(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo)));
                                    plot(x_aic,f_aic,'Color',[1 0 0])
                                else
                                    [f_aic,x_aic] = drg_ecdf(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo)));
                                    plot(x_aic,f_aic,'Color',[0 0 1])
                                end
                            end
                            
                            
                            
                            %Save data for ranksum
                            ii_rank=ii_rank+1;
                            dBpower_rank(ii_rank).data=delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo));
                            dBpower_rank(ii_rank).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                            dBpower_rank(ii_rank).per_ii=per_ii;
                            dBpower_rank(ii_rank).grNo=grNo;
                            dBpower_rank(ii_rank).evNo=evNo;
                            
                            no_points=length(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo)));
                            glm_dBpower_rank.data(glm_ii+1:glm_ii+no_points)=delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo));
                            glm_dBpower_rank.per_ii(glm_ii+1:glm_ii+no_points)=per_ii;
                            glm_dBpower_rank.grNo(glm_ii+1:glm_ii+no_points)=grNo;
                            glm_dBpower_rank.evNo(glm_ii+1:glm_ii+no_points)=evNo;
                            glm_ii=glm_ii+no_points;
                            
                            
                            maxdB=max([maxdB max(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo)))]);
                            mindB=min([mindB min(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo)))]);
                            
                            
                        end
                    end
                    
                end
                legend([handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
                    [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
                title([freq_names{bwii} ' delta dB power per mouse, per electrode for ' evTypeLabels{evNo}])
                xlabel('delta power (dB)')
                ylabel('Probability')
            end
            
            for evNo=1:length(eventType)
                subplot(length(eventType),1,evNo)
                xlim([mindB-0.1*(maxdB-mindB) maxdB+0.1*(maxdB-mindB)])
            end
            
            %Perform the glm
            fprintf(1, ['\n\nglm for delta dB power per electrode for ' freq_names{bwii} '\n'])
            tbl = table(glm_dBpower_rank.data',glm_dBpower_rank.grNo',glm_dBpower_rank.per_ii',glm_dBpower_rank.evNo',...
                'VariableNames',{'dBpower','group','perCorr','event'});
            mdl = fitglm(tbl,'dBpower~group+perCorr+event+group*event*perCorr'...
                ,'CategoricalVars',[2,3,4])
            
            %Do ranksum/t test
            fprintf(1, ['\n\nRanksum or t-test p values for delta dB power per electrode for ' freq_names{bwii} '\n'])
            try
                [output_data] = drgMutiRanksumorTtest(dBpower_rank);
                fprintf(1, '\n\n')
            catch
            end
            
            
        end
        fprintf(1, ['\n\n'])
        
        %Now plot the histograms and the average per mouse LFP power
        %computed per mouse
        
        for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            %Plot the average
            
            try
                close(bwii+8)
            catch
            end
            hFig=figure(bwii+8);
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            
            
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            per_mouse_data_delta_dB=[];
            per_mouse_prof_naive=[];
            per_mouse_events=[];
            per_mouse_groups=[];
            
            bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                
                for evNo=1:length(eventType)
                    
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        %                         bar_offset=21-evNo*3+(2-per_ii);
                        
                        if sum(eventType==3)>0
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(evNo-1);
                        else
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        end
                        
                        these_offsets(per_ii)=bar_offset;
                        
                        %Compute per mouse avearge for this group
                        no_mice_for_this_group=0;
                        each_mouse_average_delta_dB=[];
                        for mouseNo=1:max(delta_dB_mouseNo_per_mouse)
                            if sum((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_mouseNo_per_mouse==mouseNo)&(delta_dB_group_no_per_mouse==grNo))>0
                                no_mice_for_this_group=no_mice_for_this_group+1;
                                each_mouse_average_delta_dB(no_mice_for_this_group)=mean(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_mouseNo_per_mouse==mouseNo)&(delta_dB_group_no_per_mouse==grNo)));
                            end
                        end
                        
                        if no_mice_for_this_group>0
                            
                            include_group=1;
                            
                            if per_ii==1
                                bar(bar_offset,mean(each_mouse_average_delta_dB),'r','LineWidth', 3)
                            else
                                bar(bar_offset,mean(each_mouse_average_delta_dB),'b','LineWidth', 3)
                            end
                            
                            
                            %In the future add lines linking the points
                            plot((bar_offset)*ones(1,no_mice_for_this_group),each_mouse_average_delta_dB,'o',...
                                'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            
                            %Average and CI
                            plot(bar_offset,mean(each_mouse_average_delta_dB),'ok','LineWidth', 3)
                            if no_mice_for_this_group>2
                                CI = bootci(1000, {@mean, each_mouse_average_delta_dB},'type','cper');
                                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            end
                            
                            
                            %Save data for anovan
                            per_mouse_data_delta_dB=[per_mouse_data_delta_dB each_mouse_average_delta_dB];
                            per_mouse_prof_naive=[per_mouse_prof_naive per_ii*ones(1,no_mice_for_this_group)];
                            per_mouse_events=[per_mouse_events evNo*ones(1,no_mice_for_this_group)];
                            per_mouse_groups=[per_mouse_groups grNo*ones(1,no_mice_for_this_group)];
                        end
                    end
                end
                if include_group==1
                    ii_gr_included=ii_gr_included+1;
                    groups_included(ii_gr_included)=grNo;
                end
            end
            if sum(eventType==3)>0
                title([freq_names{bwii} ' auROC per mouse, electrode average'])
            else
                title([freq_names{bwii} ' auROC per mouse, electrode avearage concentrations two steps appart'])
            end
            
            
            
            %Annotations identifying groups
            x_interval=0.8/ii_gr_included;
            for ii=1:ii_gr_included
                annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.8 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %x labels
            to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
            sorted_A=sortrows(to_sort);
            sorted_bar_lab_loc=sorted_A(:,1);
            for ii=1:length(bar_lab_loc)
                sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
            end
            xticks(sorted_bar_lab_loc)
            xticklabels(sorted_bar_labels)
            
            if sum(eventType==3)==0
                xlabel('Concentration (%)')
            end
            
            
            ylabel('Delta power (dB)')
            
            
            
        end
        fprintf(1, ['\n\n'])
        
        %Display the auROC for all trials per mouse (per electrode) for for concentrations separated by two log steps
        for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            %Plot the average
            
            try
                close(bwii+12)
            catch
            end
            hFig=figure(bwii+12);
            
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            data_auROC=[];
            prof_naive=[];
            within_between=[];
            mice=[];
            groups=[];
            
            bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            
            
            
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                
                for between=0:1
                    
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        %bar_offset=21-evNo*3+(2-per_ii);
                        if sum(eventType==3)>0
                            bar_offset=(grNo-1)*(3.5*2)+(2-(per_ii-1))+3*(between-1);
                        else
                            bar_offset=(grNo-1)*(3.5*2)+(2-(per_ii-1))+3*(2-between);
                        end
                        
                        these_offsets(per_ii)=bar_offset;
                        
                        if sum((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))>0
                            
                            include_group=1;
                            
                            
                            if per_ii==1
                                bar(bar_offset,mean(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))),'r','LineWidth', 3)
                            else
                                bar(bar_offset,mean(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))),'b','LineWidth', 3)
                            end
                            
                            %Individual points; in the future add lines linking the
                            %points?
                            plot((bar_offset)*ones(1,sum((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))),...
                                auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2)),'o',...
                                'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            
                            %Average and CI
                            plot(bar_offset,mean(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))),'ok','LineWidth', 3)
                            CI = bootci(1000, {@mean, auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))},'type','cper');
                            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            
                            %Save data for anovan
                            data_auROC=[data_auROC auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))];
                            prof_naive=[prof_naive per_ii*ones(1,sum((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2)))];
                            within_between=[within_between between*ones(1,sum((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2)))];
                            mice=[mice ROCmouse((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))];
                            groups=[groups grNo*ones(1,sum((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2)))];
                        end
                    end
                    if include_group==1
                        bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                        no_ev_labels=no_ev_labels+1;
                        if sum(eventType==3)==0
                            if between==0
                                bar_labels{no_ev_labels}='within';
                            else
                                bar_labels{no_ev_labels}='between';
                            end
                            
                        end
                    end
                end
                if include_group==1
                    ii_gr_included=ii_gr_included+1;
                    groups_included(ii_gr_included)=grNo;
                end
            end
            
            title([freq_names{bwii} ' auROC per mouse, per electrode'])
            
            %Annotations identifying groups
            
            if sum(eventType==3)==0
                x_interval=0.8/ii_gr_included;
                for ii=1:ii_gr_included
                    annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.8 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
                end
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %x labels
            if sum(eventType==3)==0
                to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
                sorted_A=sortrows(to_sort);
                sorted_bar_lab_loc=sorted_A(:,1);
                for ii=1:length(bar_lab_loc)
                    sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
                end
                xticks(sorted_bar_lab_loc)
                xticklabels(sorted_bar_labels)
            else
                for ii=1:ii_gr_included
                    bar_labels{ii}=handles_drgb.drgbchoices.group_no_names{ groups_included(ii)};
                end
                xticks(bar_lab_loc)
                xticklabels(bar_labels)
            end
            
            
            
            ylabel('auROC')
            
            
        end
        
        
        
        
        %Display cumulative histograms for the auROC for all trials per mouse (per electrode) for for concentrations separated by two log steps
        %This only works for Daniel's NRG, we have to modify in the future
        pvals=[];
        if sum(eventType==3)>0
            for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
                %Plot the average
                
                try
                    close(bwii+16)
                catch
                end
                hFig=figure(bwii+16);
                
                set(hFig, 'units','normalized','position',[.2 .2 .6 .6])
                
                set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                hold on
                
                ii_rank=0;
                glm_auROC=[];
                glm_ii=0;
                
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    include_group=0;
                    
                    
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        if sum((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo))>0
                            
                            include_group=1;
                            
                            
                            if per_ii==1
                                if grNo==1
                                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo)));
                                    plot(x_aic,f_aic,'Color',[1 0 0],'LineWidth',3)
                                else
                                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo)));
                                    plot(x_aic,f_aic,'Color',[0 0 1],'LineWidth',3)
                                end
                            else
                                if grNo==1
                                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo)));
                                    plot(x_aic,f_aic,'Color',[1 0 0])
                                else
                                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo)));
                                    plot(x_aic,f_aic,'Color',[0 0 1])
                                end
                            end
                            
                            
                            %Save data for ranksum
                            ii_rank=ii_rank+1;
                            ranksum_auROC(ii_rank).data=auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo));
                            ranksum_auROC(ii_rank).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' prof_naive_leg{per_ii}];
                            ranksum_auROC(ii_rank).per_ii=per_ii;
                            ranksum_auROC(ii_rank).grNo=grNo;
                            
                            glm_auROC.data(glm_ii+1:glm_ii+length(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo))))=auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo));
                            glm_auROC.perCorr(glm_ii+1:glm_ii+length(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo))))=per_ii;
                            glm_auROC.grNo(glm_ii+1:glm_ii+length(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo))))=grNo;
                            glm_ii=glm_ii+length(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo)));
                            
                        end
                    end
                    
                    if include_group==1
                        ii_gr_included=ii_gr_included+1;
                        groups_included(ii_gr_included)=grNo;
                    end
                end
                
                title([freq_names{bwii} ' auROC per mouse, per electrode'])
                
                legend([handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive']...
                    ,[handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
                
                xlabel('auROC')
                ylabel('Cumulative probability')
                
                
                %Perform the glm
                fprintf(1, ['\n\nglm for auROC for ' freq_names{bwii} '\n'])
                tbl = table(glm_auROC.data',glm_auROC.grNo',glm_auROC.perCorr',...
                    'VariableNames',{'auROC','group','perCorr'});
                mdl = fitglm(tbl,'auROC~group+perCorr+perCorr*group','CategoricalVars',[2,3])
                
                %Do ranksum/t test
                fprintf(1, ['\n\nRanksum or t-test p values for auROC for ' freq_names{bwii} '\n'])
                try
                    [output_data] = drgMutiRanksumorTtest(ranksum_auROC);
                    fprintf(1, '\n\n')
                catch
                end
                
                
            end
        end
        fprintf(1, ['\n\n'])
        
        %Plot the WB delta dB
        maxdB=-100000;
        mindB=1000000;
        
        
        glm_dbWB=[];
        glm_dbWB_ii=0;
        
        no_out=0;
        
        for grNo=1:max(handles_drgb.drgbchoices.group_no)
            
            try
                close(grNo+20)
            catch
            end
            hFig=figure(grNo+20);
            
            set(hFig, 'units','normalized','position',[.1 .5 .35 .25])
            
            %set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            for per_ii=2:-1:1      %performance bins. blue = naive, red = proficient
                
                for evNo=1:2
                    
                    this_mean_dbWB=zeros(1,length(frequency));
                    this_mean_dbWB(1,:)=mean(delta_dB_WB_per_mouse((delta_dB_WB_perii_per_mouse==per_ii)&(delta_dB_WB_evNo_per_mouse==evNo)&(delta_dB_WB_group_no_per_mouse==grNo),:),1);
                    
                    no_out=no_out+1;
                    handles_out.delta_dB_WB(no_out).mean_dbWB=this_mean_dbWB;
                    handles_out.delta_dB_WB(no_out).evNo=evNo;
                    handles_out.delta_dB_WB(no_out).grNo=grNo;
                    handles_out.delta_dB_WB(no_out).per_ii=per_ii;
                    
                    no_dbWB=sum((delta_dB_WB_perii_per_mouse==per_ii)&(delta_dB_WB_evNo_per_mouse==evNo)&(delta_dB_WB_group_no_per_mouse==grNo));
                    these_dbWB=zeros(no_dbWB,length(frequency));
                    these_dbWB(:,:)=delta_dB_WB_per_mouse((delta_dB_WB_perii_per_mouse==per_ii)&(delta_dB_WB_evNo_per_mouse==evNo)&(delta_dB_WB_group_no_per_mouse==grNo),:);
                    
                    CI=[];
                    CI = bootci(1000, {@mean, these_dbWB})';
                    CI(:,1)= this_mean_dbWB'-CI(:,1);
                    CI(:,2)=CI(:,2)- this_mean_dbWB';
                    
                    
                    
                    if evNo==1
                        if per_ii==1
                            %S+ Proficient
                            [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[158/255 31/255 99/255]);
                        else
                            %S+ Naive
                            [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[238/255 111/255 179/255]);
                        end
                    else
                        if per_ii==1
                            %S- Proficient
                            [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[0 114/255 178/255]);
                        else
                            %S- naive
                            [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[80/255 194/255 255/255]);
                        end
                    end
                    
                    
                    for iif=1:length(frequency)
                        glm_dbWB.data(glm_dbWB_ii+1:glm_dbWB_ii+no_dbWB)=these_dbWB(:,iif);
                        glm_dbWB.freq(glm_dbWB_ii+1:glm_dbWB_ii+no_dbWB)=frequency(iif);
                        glm_dbWB.per_ii(glm_dbWB_ii+1:glm_dbWB_ii+no_dbWB)=per_ii;
                        glm_dbWB.grNo(glm_dbWB_ii+1:glm_dbWB_ii+no_dbWB)=grNo;
                        glm_dbWB.evNo(glm_dbWB_ii+1:glm_dbWB_ii+no_dbWB)=evNo;
                        glm_dbWB_ii=glm_dbWB_ii+no_dbWB;
                    end
                    
                    maxdB=max([maxdB max(this_mean_dbWB)]);
                    mindB=min([mindB min(this_mean_dbWB)]);
                end
            end
            
            xlim([4 95])
            ylim([mindB-0.05*(maxdB-mindB) maxdB+0.05*(maxdB-mindB)])
            
            title(['Delta power vs. frequency ' handles_drgb.drgbchoices.group_no_names{grNo}])
            xlabel('Frequency (Hz)')
            ylabel('dB')
            
        end
        
        %Plot the WB dB
        maxdB=-100000;
        mindB=1000000;
        
        
        glm_dbWB=[];
        glm_dbWB_ii=0;
        
        no_out=0;
        
        for grNo=1:max(handles_drgb.drgbchoices.group_no)
            
            try
                close(grNo+23)
            catch
            end
            hFig=figure(grNo+23);
            
            set(hFig, 'units','normalized','position',[.1 .5 .35 .25])
            
            
            %set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            for per_ii=2:-1:1      %performance bins. blue = naive, red = proficient
                
                for evNo=1:2
                    
                    this_mean_dbWB=zeros(1,length(frequency));
                    this_mean_dbWB(1,:)=mean(dB_WB_per_mouse((delta_dB_WB_perii_per_mouse==per_ii)&(delta_dB_WB_evNo_per_mouse==evNo)&(delta_dB_WB_group_no_per_mouse==grNo),:),1);
                    
                    no_out=no_out+1;
                    handles_out.dB_WB(no_out).mean_dbWB=this_mean_dbWB;
                    handles_out.dB_WB(no_out).evNo=evNo;
                    handles_out.dB_WB(no_out).grNo=grNo;
                    handles_out.dB_WB(no_out).per_ii=per_ii;
                    
                    no_dbWB=sum((delta_dB_WB_perii_per_mouse==per_ii)&(delta_dB_WB_evNo_per_mouse==evNo)&(delta_dB_WB_group_no_per_mouse==grNo));
                    these_dbWB=zeros(no_dbWB,length(frequency));
                    these_dbWB(:,:)=dB_WB_per_mouse((delta_dB_WB_perii_per_mouse==per_ii)&(delta_dB_WB_evNo_per_mouse==evNo)&(delta_dB_WB_group_no_per_mouse==grNo),:);
                    
                    CI=[];
                    CI = bootci(1000, {@mean, these_dbWB})';
                    CI(:,1)= this_mean_dbWB'-CI(:,1);
                    CI(:,2)=CI(:,2)- this_mean_dbWB';
                    
                    
                    
                    if evNo==1
                        if per_ii==1
                            %S+ Proficient
                            [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[158/255 31/255 99/255]);
                        else
                            %S+ Naive
                            [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[238/255 111/255 179/255]);
                        end
                    else
                        if per_ii==1
                            %S- Proficient
                            [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[0 114/255 178/255]);
                        else
                            %S- naive
                            [hlCR, hpCR] = boundedline(frequency',this_mean_dbWB', CI, 'cmap',[80/255 194/255 255/255]);
                        end
                    end
                    
                    
                    for iif=1:length(frequency)
                        glm_dbWB.data(glm_dbWB_ii+1:glm_dbWB_ii+no_dbWB)=these_dbWB(:,iif);
                        glm_dbWB.freq(glm_dbWB_ii+1:glm_dbWB_ii+no_dbWB)=frequency(iif);
                        glm_dbWB.per_ii(glm_dbWB_ii+1:glm_dbWB_ii+no_dbWB)=per_ii;
                        glm_dbWB.grNo(glm_dbWB_ii+1:glm_dbWB_ii+no_dbWB)=grNo;
                        glm_dbWB.evNo(glm_dbWB_ii+1:glm_dbWB_ii+no_dbWB)=evNo;
                        glm_dbWB_ii=glm_dbWB_ii+no_dbWB;
                    end
                    
                    maxdB=max([maxdB max(this_mean_dbWB)]);
                    mindB=min([mindB min(this_mean_dbWB)]);
                end
            end
            
            xlim([4 95])
            ylim([mindB-0.05*(maxdB-mindB) maxdB+0.05*(maxdB-mindB)])
            
            title(['dB power vs. frequency ' handles_drgb.drgbchoices.group_no_names{grNo}])
            xlabel('Frequency (Hz)')
            ylabel('dB')
            
        end
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) handles_pars.output_suffix],'handles_out')
        
        pffft=1;
        
        
    case 21
        % Multiclass ROC analysis of coherence for naive and proficient
        % mice for different epochs (concentrations or S+ vs. S-) and different groups. Analyzed per mouse
        
        deltaCxy=[];        %Change in coherence after stimulation with the odor
        Cxy=[];             %This is the coherence before odor on
        p_vals_ROC=[];
        per_mouse_no_ROCs=0;
        per_mouse_ROCout=[];
        per_mouse_p_vals_ROC=[];
        %         deltaCxy_Ev1=[];
        %         no_Ev1=0;
        for evNo=1:length(eventType)
            evNo_out(evNo).noWB=0;
        end
        %         deltaCxy_Ev1WB=[];
        %         deltaCxy_Ev2WB=[];
        deltaCxy_No_per_mouse=0;
        deltaCxy_per_mouse=[];
        deltaCxy_perii_per_mouse=[];
        deltaCxy_evNo_per_mouse=[];
        deltaCxy_bwii_per_mouse=[];
        deltaCxy_mouseNo_per_mouse=[];
        deltaCxy_elec_pair_per_mouse=[];
        mouse_included=[];
        
        deltaCxy_af_No_per_mouse=0;
        deltaCxy_af_per_mouse=[];
        deltaCxy_af_perii_per_mouse=[];
        deltaCxy_af_evNo_per_mouse=[];
        deltaCxy_af_mouseNo_per_mouse=[];
        deltaCxy_af_elec_pair_per_mouse=[];
        deltaCxy_af_group_no_per_mouse=[];
        
        prof_naive_leg{1}='Proficient';
        prof_naive_leg{2}='Naive';
        
        handles_out=[];
        handles_out.dcoh_ii=0;
        handles_out.auc_ii=0;
        handles_out.dcohaf_ii=0;
        
        
        fprintf(1, ['Pairwise auROC analysis for coherence analysis of LFP\n\n'])
        %         p_vals=[];
        no_files=max(files);
        
        %Initialize ROC
        no_ROCs=0;
        ROCelec_pair=[];
        ROCgroups=[];
        ROCmouse=[];
        ROCbwii=[];
        ROCper_ii=[];
        ROCEvNo1=[];
        ROCEvNo2=[];
        auROC=[];
        p_valROC=[];
        p_vals_ROC=[];
        
        
        
        szpc=size(percent_windows);
        for per_ii=1:szpc(1)
            
            for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
                mouse_has_files=0;
                theseEvNosPerEl=[];
                for evNo=1:length(eventType)
                    for bwii=1:no_bandwidths
                        for elec=1:16
                            theseEvNosPerEl(evNo,bwii,elec).noEv=0;
                        end
                    end
                end
                
                for elec_pair=1:max(elec_pair_No)
                    
                    theseEvNos=[];
                    theseEvNos_af=[];
                    for evNo=1:length(eventType)
                        for bwii=1:no_bandwidths
                            theseEvNos(evNo,bwii).noEv=0;
                        end
                        theseEvNos_af(evNo).noEv=0;
                    end
                    
                    this_file_is_nan=0;
                    for fileNo=1:no_files
                        if sum(files==fileNo)>0
                            if handles_drgb.drgbchoices.mouse_no(fileNo)==mouseNo
                                group_no_per_mouse(mouseNo)=handles_drgb.drgbchoices.group_no(fileNo);
                                lfppairNo=find((files_per_lfp==fileNo)&(elec_pair_No==elec_pair)&(window_per_lfp==winNo));
                                lfppairNo_ref=find((files_per_lfp==fileNo)&(elec_pair_No==elec_pair)&(window_per_lfp==refWin));
                                
                                if (~isempty(handles_drgb.drgb.lfpevpair(lfppairNo).all_Cxy_timecourse))&(~isempty(handles_drgb.drgb.lfpevpair(lfppairNo_ref).all_Cxy_timecourse))
                                    
                                    percent_mask=[];
                                    percent_mask=logical((handles_drgb.drgb.lfpevpair(lfppairNo).perCorrCoh>=percent_windows(per_ii,1))...
                                        &(handles_drgb.drgb.lfpevpair(lfppairNo).perCorrCoh<=percent_windows(per_ii,2)));
                                    
                                    
                                    for evNo=1:length(eventType)
                                        
                                        noWB_for_evNo(evNo)=-1;
                                        
                                        trials_in_event_refEv=[];
                                        trials_in_event_refEv=(handles_drgb.drgb.lfpevpair(lfppairNo_ref).which_event_coh(eventType(evNo),:)==1)&percent_mask;
                                        
                                        
                                        trials_in_event_Ev=[];
                                        trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfppairNo).which_event_coh(eventType(evNo),:)==1)&percent_mask;
                                        
                                        if (sum(trials_in_event_Ev)>=1)
                                            
                                            
                                            
                                            % Ev1
                                            this_deltaCxy=zeros(sum(trials_in_event_Ev),length(frequency));
                                            
                                            %Use this if you want coherence in the presence of the odorant
                                            %                                             this_deltaCxy(:,:)=mean(handles_drgb.drgb.lfpevpair(lfppairNo).all_Cxy_timecourse(trials_in_event_Ev,:,:),3);
                                            
                                            %Use this if you want delta coherence
                                            this_deltaCxy(:,:)=mean(handles_drgb.drgb.lfpevpair(lfppairNo).all_Cxy_timecourse(trials_in_event_Ev,:,:)...
                                                -handles_drgb.drgb.lfpevpair(lfppairNo_ref).all_Cxy_timecourse(trials_in_event_Ev,:,:),3);
                                            
                                            if sum(isnan(this_deltaCxy))>0
                                                this_file_is_nan=1;
                                            end
                                            
                                            this_Cxy=zeros(sum(trials_in_event_Ev),length(frequency));
                                            this_Cxy(:,:)=mean(handles_drgb.drgb.lfpevpair(lfppairNo_ref).all_Cxy_timecourse(trials_in_event_Ev,:,:),3);
                                            
                                            
                                            %Wide band spectrum
                                            evNo_out(evNo).noWB=evNo_out(evNo).noWB+1;
                                            evNo_out(evNo).deltaCxy_EvWB(evNo_out(evNo).noWB,1:length(frequency))=mean(this_deltaCxy,1);
                                            evNo_out(evNo).per_ii(evNo_out(evNo).noWB)=per_ii;
                                            evNo_out(evNo).groupNo(evNo_out(evNo).noWB)=handles_drgb.drgbchoices.group_no(fileNo);
                                            
                                            noWB_for_evNo(evNo)=evNo_out(evNo).noWB;
                                            
                                            
                                            %Do per bandwidth analysis
                                            for bwii=1:no_bandwidths
                                                
                                                this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                                
                                                %Enter the  Ev1
                                                this_deltaCxy_Ev=zeros(sum(trials_in_event_Ev),1);
                                                this_deltaCxy_Ev=mean(this_deltaCxy(:,this_band),2);
                                                
                                                theseEvNos(evNo,bwii).this_deltaCxy_Ev(1,theseEvNos(evNo,bwii).noEv+1:theseEvNos(evNo,bwii).noEv+sum(trials_in_event_Ev))=this_deltaCxy_Ev';
                                                %                                                         theseEvNos(evNo,bwii).noEv=theseEvNos(evNo,bwii).noEv+sum(trials_in_event_Ev);
                                                
                                                %Enter the  Ev1
                                                this_Cxy_Ev=zeros(sum(trials_in_event_Ev),1);
                                                this_Cxy_Ev=mean(this_Cxy(:,this_band),2);
                                                
                                                theseEvNos(evNo,bwii).this_Cxy_Ev(1,theseEvNos(evNo,bwii).noEv+1:theseEvNos(evNo,bwii).noEv+sum(trials_in_event_Ev))=this_Cxy_Ev';
                                                theseEvNos(evNo,bwii).groupNo(1,theseEvNos(evNo,bwii).noEv+1:theseEvNos(evNo,bwii).noEv+sum(trials_in_event_Ev))=handles_drgb.drgbchoices.group_no(fileNo)*ones(1,sum(trials_in_event_Ev));
                                                theseEvNos(evNo,bwii).noEv=theseEvNos(evNo,bwii).noEv+sum(trials_in_event_Ev);
                                                
                                                mouse_has_files=1;
                                            end
                                            
                                            %Now save all frequencies
                                            %Enter the  Ev1
                                            theseEvNos_af(evNo).this_deltaCxy_Ev(theseEvNos_af(evNo).noEv+1:theseEvNos_af(evNo).noEv+sum(trials_in_event_Ev),:)=this_deltaCxy;
                                            
                                            %Enter the  Ev1
                                            theseEvNos_af(evNo).this_Cxy_Ev(1,theseEvNos_af(evNo).noEv+1:theseEvNos_af(evNo).noEv+sum(trials_in_event_Ev))=this_Cxy_Ev';
                                            theseEvNos_af(evNo).groupNo(1,theseEvNos_af(evNo).noEv+1:theseEvNos_af(evNo).noEv+sum(trials_in_event_Ev))=handles_drgb.drgbchoices.group_no(fileNo)*ones(1,sum(trials_in_event_Ev));
                                            theseEvNos_af(evNo).noEv=theseEvNos_af(evNo).noEv+sum(trials_in_event_Ev);
                                            
                                            fprintf(1, ['%d trials in event No %d succesfully processed for file No %d electrode pair %d\n'],sum(trials_in_event_Ev), evNo,fileNo,elec_pair);
                                            
                                        else
                                            
                                            
                                            fprintf(1, ['%d trials in event No %d fewer than minimum trials per event ' evTypeLabels{evNo} ' for file No %d electrode pair %d\n'],sum(trials_in_event_Ev), min_trials_per_event,fileNo,elec_pair);
                                            
                                            
                                        end
                                        
                                        
                                    end
                                    
                                    
                                else
                                    
                                    fprintf(1, ['Empty allPower for file No %d electrode %d\n'],fileNo,elec_pair);
                                    
                                end
                                
                                
                            end %if mouseNo
                        end
                    end %fileNo
                    
                    if this_file_is_nan==1
                        fprintf(1, ['WARNING file No %d has NaN delta Cxy\n\n'])
                    end
                    
                    if theseEvNos(evNo,1).noEv>0
                        
                        %Calculate coherence
                        if mouse_has_files==1
                            mouse_included(mouseNo)=1;
                            %Calculate per mouse per electrode delta_dB
                            for evNo=1:length(eventType)
                                for bwii=1:no_bandwidths
                                    if theseEvNos(evNo,bwii).noEv>0
                                        deltaCxy_No_per_mouse=deltaCxy_No_per_mouse+1;
                                        deltaCxy_per_mouse(deltaCxy_No_per_mouse)=mean(theseEvNos(evNo,bwii).this_deltaCxy_Ev);
                                        Cxy_per_mouse(deltaCxy_No_per_mouse)=mean(theseEvNos(evNo,bwii).this_Cxy_Ev);
                                        deltaCxy_perii_per_mouse(deltaCxy_No_per_mouse)=per_ii;
                                        deltaCxy_evNo_per_mouse(deltaCxy_No_per_mouse)=evNo;
                                        deltaCxy_bwii_per_mouse(deltaCxy_No_per_mouse)=bwii;
                                        deltaCxy_mouseNo_per_mouse(deltaCxy_No_per_mouse)=mouseNo;
                                        deltaCxy_elec_pair_per_mouse(deltaCxy_No_per_mouse)=elec_pair;
                                        deltaCxy_group_no_per_mouse(deltaCxy_No_per_mouse)=group_no_per_mouse(mouseNo);
                                    end
                                end
                                if theseEvNos_af(evNo).noEv>0
                                    deltaCxy_af_No_per_mouse=deltaCxy_af_No_per_mouse+1;
                                    deltaCxy_af_per_mouse(deltaCxy_af_No_per_mouse,:)=mean(theseEvNos_af(evNo).this_deltaCxy_Ev,1);
                                    deltaCxy_af_perii_per_mouse(deltaCxy_af_No_per_mouse)=per_ii;
                                    deltaCxy_af_evNo_per_mouse(deltaCxy_af_No_per_mouse)=evNo;
                                    deltaCxy_af_mouseNo_per_mouse(deltaCxy_af_No_per_mouse)=mouseNo;
                                    deltaCxy_af_elec_pair_per_mouse(deltaCxy_af_No_per_mouse)=elec_pair;
                                    deltaCxy_af_group_no_per_mouse(deltaCxy_af_No_per_mouse)=group_no_per_mouse(mouseNo);
                                end
                                
                            end
                            
                            %Calculate per electrode ROC
                            can_calculate_auroc=1;
                            if can_calculate_auroc==1
                                
                                for evNo1=1:length(eventType)
                                    for evNo2=evNo1+1:length(eventType)
                                        for bwii=1:no_bandwidths
                                            if (theseEvNos(evNo1,bwii).noEv>5)&...
                                                    (theseEvNos(evNo2,bwii).noEv>5)
                                                
                                                %Enter Ev1
                                                this_mouse_deltaCxyEv1=[];
                                                this_mouse_deltaCxyEv1=theseEvNos(evNo1,bwii).this_deltaCxy_Ev;
                                                trials_in_event_Ev1=theseEvNos(evNo1,bwii).noEv;
                                                
                                                roc_data=[];
                                                roc_data(1:trials_in_event_Ev1,1)=this_mouse_deltaCxyEv1;
                                                roc_data(1:trials_in_event_Ev1,2)=zeros(trials_in_event_Ev1,1);
                                                
                                                
                                                %Enter Ev2
                                                this_mouse_deltaCxyEv2=[];
                                                this_mouse_deltaCxyEv2=theseEvNos(evNo2,bwii).this_deltaCxy_Ev;
                                                trials_in_event_Ev2=theseEvNos(evNo2,bwii).noEv;
                                                
                                                
                                                roc_data(trials_in_event_Ev1+1:trials_in_event_Ev2+trials_in_event_Ev1,1)=this_mouse_deltaCxyEv2;
                                                roc_data(trials_in_event_Ev1+1:trials_in_event_Ev2+trials_in_event_Ev1,2)=ones(trials_in_event_Ev2,1);
                                                
                                                
                                                %Find  per electrode ROC
                                                no_ROCs=no_ROCs+1;
                                                
                                                
                                                ROCelec_pair(no_ROCs)=elec_pair;
                                                ROCgroups(no_ROCs)=group_no_per_mouse(mouseNo);
                                                ROCmouse(no_ROCs)=mouseNo;
                                                ROCbwii(no_ROCs)=bwii;
                                                ROCper_ii(no_ROCs)=per_ii;
                                                ROCEvNo1(no_ROCs)=evNo1;
                                                ROCEvNo2(no_ROCs)=evNo2;
                                                
                                                
                                                roc=[];
                                                roc=roc_calc(roc_data,0,0.05,0);
                                                auROC(no_ROCs)=roc.AUC-0.5;
                                                p_valROC(no_ROCs)=roc.p;
                                                p_vals_ROC=[p_vals_ROC roc.p];
                                                
                                                
                                                %I have this code here to plot the ROC
                                                
                                                show_roc=0;
                                                if (show_roc==1)&(mouseNo==4)
                                                    %I have this code here to plot the ROC
                                                    roc=roc_calc(roc_data,0,0.05,1);
                                                    
                                                    %Do the histograms
                                                    try
                                                        close(2)
                                                    catch
                                                    end
                                                    figure(2)
                                                    
                                                    hold on
                                                    
                                                    max_deltaCxy=max([max(this_mouse_deltaCxyEv1) max(this_mouse_deltaCxyEv2)]);
                                                    min_deltaCxy=min([min(this_mouse_deltaCxyEv1) min(this_mouse_deltaCxyEv2)]);
                                                    
                                                    edges=[min_deltaCxy-0.1*(max_deltaCxy-min_deltaCxy):(max_deltaCxy-min_deltaCxy)/20:max_deltaCxy+0.1*(max_deltaCxy-min_deltaCxy)];
                                                    histogram(this_mouse_deltaCxyEv1,edges,'FaceColor','b','EdgeColor','b')
                                                    histogram(this_mouse_deltaCxyEv2,edges,'FaceColor','r','EdgeColor','r')
                                                    xlabel('delta coherence')
                                                    title(['Histogram for S+ vs S-'])
                                                    pffft=1;
                                                end
                                                
                                            end
                                            
                                            %
                                        end
                                    end
                                end
                                %                                         end
                                
                                
                            end
                            
                        else
                            mouse_included(mouseNo)=0;
                        end
                    end
                end
                
                
                
                if mouse_has_files==1
                    mouse_included(mouseNo)=1;
                else
                    mouse_included(mouseNo)=0;
                end
            end
            
            
            
        end
        %         fprintf(1, '\n\n')
        %
        %         pFDRauROC=drsFDRpval(p_vals_ROC);
        %         fprintf(1, ['pFDR for per electrode pair auROC  = %d\n\n'],pFDRauROC);
        %
        %         per_mouse_pFDRauROC=drsFDRpval(per_mouse_p_vals_ROC);
        %         fprintf(1, ['pFDR for per mouse auROC  = %d\n\n'],per_mouse_pFDRauROC);
        
        
        if length(eventType) ==4
            %do Hits, etc for proficient
            
            figureNo = 0;
            %Now plot the average per electrode, all sessions for each mouse
            %         edges=[-25:0.5:15];
            edges=[-0.5:0.05:0.5];
            rand_offset=0.8;
            
            for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
                %Plot the average
                figureNo = figureNo + 1;
                try
                    close(figureNo)
                catch
                end
                hFig=figure(figureNo);
                
                
                %             try
                %                 close(bwii)
                %             catch
                %             end
                %             hFig=figure(bwii);
                
                set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
                
                set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                hold on
                
                data_delta_dB=[];
                prof_naive=[];
                events=[];
                mice=[];
                electrodes=[];
                groups=[];
                
                bar_lab_loc=[];
                no_ev_labels=0;
                ii_gr_included=0;
                %defining and initinalizing all bar
                all_bar=[];
                bar_offset = 0;
                
                
                input_data=[];
                ii_rank=0;
                glm_coh=[];
                glm_ii=0;
                
                %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                
                for evNo=1:length(eventType)
                    
                    per_included=0;
                    these_dB_per_e=[];
                    there_are_NaNs=0;
                    
                    
                    per_ii=1;      %proficient
                    for grNo=1:max(handles_drgb.drgbchoices.group_no)
                        bar_offset = bar_offset +1;
                        
                        %                         if sum(eventType==3)>0
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(evNo-1);
                        %                         else
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        %                         end
                        %
                        %                         these_offsets(per_ii)=bar_offset;
                        bar_offset = bar_offset + 1;
                        if sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))>0
                            
                            include_group=1;
                            
                            switch grNo
                                case 1
                                    bar(bar_offset,mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))),'g','LineWidth', 3,'EdgeColor','none')
                                case 2
                                    bar(bar_offset,mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))),'LineWidth', 3,'EdgeColor','none','FaceColor',[0.3010 0.7450 0.9330])
                                case 3
                                    bar(bar_offset,mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))),'y','LineWidth', 3,'EdgeColor','none')
                            end
                            
                            handles_out.dcoh_ii=handles_out.dcoh_ii+1;
                            handles_out.dcoh_values(handles_out.dcoh_ii).bwii=bwii;
                            handles_out.dcoh_values(handles_out.dcoh_ii).evNo=evNo;
                            handles_out.dcoh_values(handles_out.dcoh_ii).per_ii=per_ii;
                            handles_out.dcoh_values(handles_out.dcoh_ii).groupNo=grNo;
                            handles_out.dcoh_values(handles_out.dcoh_ii).dcoh=mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
                            
                            %Save the mean per mouse
                            these_mice=deltaCxy_mouseNo_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo));
                            handles_out.dcoh_values(handles_out.dcoh_ii).noMice=0;
                            for iiMice=min(these_mice):max(these_mice)
                                if sum(these_mice==iiMice)>0
                                    handles_out.dcoh_values(handles_out.dcoh_ii).noMice=handles_out.dcoh_values(handles_out.dcoh_ii).noMice+1;
                                    handles_out.dcoh_values(handles_out.dcoh_ii).mouseNo(handles_out.dcoh_values(handles_out.dcoh_ii).noMice)=iiMice;
                                    handles_out.dcoh_values(handles_out.dcoh_ii).dcoh_per_mouse(handles_out.dcoh_values(handles_out.dcoh_ii).noMice)=mean( deltaCxy_per_mouse((deltaCxy_mouseNo_per_mouse==iiMice)&(deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)) );
                                end
                            end
                            
                            %Violin plot
                            
                            [mean_out, CIout]=drgViolinPoint(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))...
                                ,edges,bar_offset,rand_offset,'k','k',1);
                            
                            %Save data for glm and t test/ranksum
                            these_data=deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo));
                            glm_coh.data(glm_ii+1:glm_ii+length(these_data))=these_data;
                            glm_coh.group(glm_ii+1:glm_ii+length(these_data))=grNo;
                            glm_coh.spm(glm_ii+1:glm_ii+length(these_data))=evNo;
                            glm_ii=glm_ii+length(these_data);
                            
                            %Enter the data for t-test/ranksum
                            ii_rank=ii_rank+1;
                            input_data(ii_rank).data=these_data;
                            input_data(ii_rank).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                            
                            
                            %                             %I added this to make the bars smaller
                            %                             all_bar = [all_bar mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)))];
                            %
                            %Individual points; in the future add lines linking the
                            %points?
                            
                            %                             ylim([0 1.2*max(all_bar)])
                            %                             plot((bar_offset)*ones(1,sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))),...
                            %                                 deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)),'o',...
                            %                                 'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            %                             data_for_lines(per_ii).these_dB_per_e(1:length(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))))=...
                            %                                 deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo));
                            %                             data_for_lines(per_ii).these_mice(1:length(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))))=...
                            %                                 deltaCxy_mouseNo_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo));
                            %                             data_for_lines(per_ii).these_bar_offsets=bar_offset;
                            %
                            %
                            %                             %Average and CI
                            %                             plot(bar_offset,mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
                            %                             if sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))>2
                            %                                 CI = bootci(1000, {@mean, deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))},'type','cper');
                            %                                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            %                             end
                            
                            %                             %Save data for anovan
                            %                             data_delta_dB=[data_delta_dB deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))];
                            %                             prof_naive=[prof_naive per_ii*ones(1,sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)))];
                            %                             events=[events evNo*ones(1,sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)))];
                            %                             mice=[mice deltaCxy_mouseNo_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))];
                            %                             electrodes=[electrodes deltaCxy_elec_pair_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))];
                            %                             groups=[groups deltaCxy_group_no_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))];
                            %
                        end
                    end
                    bar_offset = bar_offset + 2;
                    try
                        if per_included==2
                            for mouseNo=1:length(mouse_included)
                                if mouse_included(mouseNo)==1
                                    if (sum(data_for_lines(1).these_mice==mouseNo)>0)&(sum(data_for_lines(2).these_mice==mouseNo)>0)
                                        plot([data_for_lines(1).these_bar_offsets*ones(1,sum(data_for_lines(1).these_mice==mouseNo)); data_for_lines(2).these_bar_offsets*ones(1,sum(data_for_lines(1).these_mice==mouseNo))],...
                                            [data_for_lines(1).these_dB_per_e(data_for_lines(1).these_mice==mouseNo); data_for_lines(2).these_dB_per_e(data_for_lines(2).these_mice==mouseNo)],'-','Color',[0.7 0.7 0.7])
                                    end
                                end
                            end
                        end
                    catch
                        pffft=1
                    end
                    %                     if include_group==1
                    %                         bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                    %                         no_ev_labels=no_ev_labels+1;
                    %                         if sum(eventType==3)>0
                    %                             bar_labels{no_ev_labels}=evTypeLabels{evNo};
                    %                         else
                    %                             bar_labels{no_ev_labels}=num2str(concs(evNo));
                    %                         end
                    %                     end
                    
                    bar_offset = bar_offset + 3;
                    %                 if include_group==1
                    %                     ii_gr_included=ii_gr_included+1;
                    %                     groups_included(ii_gr_included)=grNo;
                    %                 end
                end
                ylim([-0.5 0.5])
                %             title([freq_names{bwii} ' average delta coherence per mouse, per electrode'])
                title([freq_names{bwii} ' average coherence per mouse, per electrode during odor'])
                
                %Annotations identifying groups
                x_interval=0.8/ii_gr_included;
                for ii=1:ii_gr_included
                    annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.8 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
                end
                
                %Proficient/Naive annotations
                annotation('textbox',[0.15 0.70 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
                annotation('textbox',[0.15 0.65 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
                
                %x labels
                %             to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
                %             sorted_A=sortrows(to_sort);
                %             sorted_bar_lab_loc=sorted_A(:,1);
                %             for ii=1:length(bar_lab_loc)
                %                 sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
                %             end
                
                
                xticks([2 4 6 13 15 17 24 26 28 35 37 39])
                xticklabels({'PwH', 'PHH', 'PKOH', 'PwM', 'PHM', 'PKOM', 'PwCR', 'PHCR', 'PKOCR', 'PwFA', 'PHFA', 'PKOFA'})
                
                
                if sum(eventType==3)==0
                    xlabel('Concentration (%)')
                end
                
                %             ylabel('Delta coherence')
                ylabel('Delta coherence')
                
                
                %             %Calculate anovan for inteaction
                %             [p,tbl,stats]=anovan(data_delta_dB,{prof_naive events mice electrodes},'varnames',{'proficient_vs_naive','events','groups','mice'},'display','off','random',4);
                %             fprintf(1, ['p value for anovan delta dB power per mouse per electrode for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
                %             fprintf(1, ['p value for anovan delta dB power  per mouse per electrode for events ' freq_names{bwii} '= %d \n'],  p(2));
                %             fprintf(1, ['p value for anovan delta dB power  per mouse per electrode for groups for ' freq_names{bwii} '= %d \n\n'],  p(3));
                %
                
                %Now do the GLM
                fprintf(1, ['\n\nglm for odor-elicited change in coherence for ' freq_names{bwii} '\n'])
                tbl = table(glm_coh.data',glm_coh.group',glm_coh.spm',...
                    'VariableNames',{'delta_coherence','group','event'});
                mdl = fitglm(tbl,'delta_coherence~group+event+group*event'...
                    ,'CategoricalVars',[2,3])
                
                
                fprintf(1, ['\n\nRanksum or t-test for delta coherence ' freq_names{bwii} '\n'])
                %Now do the ranksums
                output_data = drgMutiRanksumorTtest(input_data);
                
                fprintf(1, ['\n\n'])
                
                pffft=1;
            end
            fprintf(1, ['\n\n'])
            
        else
            figureNo = 0;
            %Now plot the average per electrode, all sessions for each mouse
            %         edges=[-25:0.5:15];
            edges=[-0.5:0.05:0.5];
            rand_offset=0.8;
            
            for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
                %Plot the average
                figureNo = figureNo + 1;
                try
                    close(figureNo)
                catch
                end
                hFig=figure(figureNo);
                
                
                %             try
                %                 close(bwii)
                %             catch
                %             end
                %             hFig=figure(bwii);
                
                set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
                
                set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                hold on
                
                data_delta_dB=[];
                prof_naive=[];
                events=[];
                mice=[];
                electrodes=[];
                groups=[];
                
                bar_lab_loc=[];
                no_ev_labels=0;
                ii_gr_included=0;
                %defining and initinalizing all bar
                all_bar=[];
                bar_offset = 0;
                
                
                input_data=[];
                ii_rank=0;
                glm_coh=[];
                glm_ii=0;
                
                %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                
                for evNo=1:length(eventType)
                    
                    per_included=0;
                    these_dB_per_e=[];
                    there_are_NaNs=0;
                    
                    
                    for per_ii=2:-1:1      %performance bins. blue = naive, red = proficient
                        for grNo=1:max(handles_drgb.drgbchoices.group_no)
                            bar_offset = bar_offset +1;
                            
                            %                         if sum(eventType==3)>0
                            %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(evNo-1);
                            %                         else
                            %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                            %                         end
                            %
                            %                         these_offsets(per_ii)=bar_offset;
                            bar_offset = bar_offset + 1;
                            if sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))>0
                                
                                include_group=1;
                                
                                switch grNo
                                    case 1
                                        bar(bar_offset,mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))),'g','LineWidth', 3,'EdgeColor','none')
                                    case 2
                                        bar(bar_offset,mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))),'LineWidth', 3,'EdgeColor','none','FaceColor',[0.3010 0.7450 0.9330])
                                    case 3
                                        bar(bar_offset,mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))),'y','LineWidth', 3,'EdgeColor','none')
                                end
                                
                                handles_out.dcoh_ii=handles_out.dcoh_ii+1;
                                handles_out.dcoh_values(handles_out.dcoh_ii).pacii=bwii;
                                handles_out.dcoh_values(handles_out.dcoh_ii).evNo=evNo;
                                handles_out.dcoh_values(handles_out.dcoh_ii).per_ii=per_ii;
                                handles_out.dcoh_values(handles_out.dcoh_ii).groupNo=grNo;
                                handles_out.dcoh_values(handles_out.dcoh_ii).dcoh=mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
                                handles_out.dcoh_values(handles_out.dcoh_ii).all_dcoh=deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo));
                                
                                %Save the mean per mouse
                                these_mice=deltaCxy_mouseNo_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo));
                                handles_out.dcoh_values(handles_out.dcoh_ii).noMice=0;
                                for iiMice=min(these_mice):max(these_mice)
                                    if sum(these_mice==iiMice)>0
                                        handles_out.dcoh_values(handles_out.dcoh_ii).noMice=handles_out.dcoh_values(handles_out.dcoh_ii).noMice+1;
                                        handles_out.dcoh_values(handles_out.dcoh_ii).mouseNo(handles_out.dcoh_values(handles_out.dcoh_ii).noMice)=iiMice;
                                        handles_out.dcoh_values(handles_out.dcoh_ii).dcoh_per_mouse(handles_out.dcoh_values(handles_out.dcoh_ii).noMice)=mean(deltaCxy_per_mouse((deltaCxy_mouseNo_per_mouse==iiMice)&(deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
                                    end
                                end
                                
                                %Violin plot
                                
                                [mean_out, CIout]=drgViolinPoint(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))...
                                    ,edges,bar_offset,rand_offset,'k','k',1);
                                
                                %Save data for glm and t test/ranksum
                                these_data=deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo));
                                glm_coh.data(glm_ii+1:glm_ii+length(these_data))=these_data;
                                glm_coh.group(glm_ii+1:glm_ii+length(these_data))=grNo;
                                glm_coh.perCorr(glm_ii+1:glm_ii+length(these_data))=per_ii;
                                glm_coh.spm(glm_ii+1:glm_ii+length(these_data))=evNo;
                                glm_ii=glm_ii+length(these_data);
                                
                                %Enter the data for t-test/ranksum
                                ii_rank=ii_rank+1;
                                input_data(ii_rank).data=these_data;
                                input_data(ii_rank).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                                
                            end
                        end
                        bar_offset = bar_offset + 2;
                        try
                            if per_included==2
                                for mouseNo=1:length(mouse_included)
                                    if mouse_included(mouseNo)==1
                                        if (sum(data_for_lines(1).these_mice==mouseNo)>0)&(sum(data_for_lines(2).these_mice==mouseNo)>0)
                                            plot([data_for_lines(1).these_bar_offsets*ones(1,sum(data_for_lines(1).these_mice==mouseNo)); data_for_lines(2).these_bar_offsets*ones(1,sum(data_for_lines(1).these_mice==mouseNo))],...
                                                [data_for_lines(1).these_dB_per_e(data_for_lines(1).these_mice==mouseNo); data_for_lines(2).these_dB_per_e(data_for_lines(2).these_mice==mouseNo)],'-','Color',[0.7 0.7 0.7])
                                        end
                                    end
                                end
                            end
                        catch
                            pffft=1
                        end
                        %                     if include_group==1
                        %                         bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                        %                         no_ev_labels=no_ev_labels+1;
                        %                         if sum(eventType==3)>0
                        %                             bar_labels{no_ev_labels}=evTypeLabels{evNo};
                        %                         else
                        %                             bar_labels{no_ev_labels}=num2str(concs(evNo));
                        %                         end
                        %                     end
                    end
                    bar_offset = bar_offset + 3;
                    %                 if include_group==1
                    %                     ii_gr_included=ii_gr_included+1;
                    %                     groups_included(ii_gr_included)=grNo;
                    %                 end
                end
                ylim([-0.5 0.5])
                %             title([freq_names{bwii} ' average delta coherence per mouse, per electrode'])
                title([freq_names{bwii} ' average coherence per mouse, per electrode during odor'])
                
                %Annotations identifying groups
                x_interval=0.8/ii_gr_included;
                for ii=1:ii_gr_included
                    annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.8 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
                end
                
                %Proficient/Naive annotations
                annotation('textbox',[0.15 0.70 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
                annotation('textbox',[0.15 0.65 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
                
                %x labels
                %             to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
                %             sorted_A=sortrows(to_sort);
                %             sorted_bar_lab_loc=sorted_A(:,1);
                %             for ii=1:length(bar_lab_loc)
                %                 sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
                %             end
                
                
                xticks([2 4 6 10 12 14 21 23 25 29 31 33])
                xticklabels({'NwS+', 'NHS+', 'NKOS+', 'PwS+', 'PHS+', 'PKOS+', 'NwS-', 'NHS-', 'NKOS-', 'PwS-', 'PHS-', 'PKOS-'})
                
                
                if sum(eventType==3)==0
                    xlabel('Concentration (%)')
                end
                
                %             ylabel('Delta coherence')
                ylabel('Delta coherence')
                
                
                %             %Calculate anovan for inteaction
                %             [p,tbl,stats]=anovan(data_delta_dB,{prof_naive events mice electrodes},'varnames',{'proficient_vs_naive','events','groups','mice'},'display','off','random',4);
                %             fprintf(1, ['p value for anovan delta dB power per mouse per electrode for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
                %             fprintf(1, ['p value for anovan delta dB power  per mouse per electrode for events ' freq_names{bwii} '= %d \n'],  p(2));
                %             fprintf(1, ['p value for anovan delta dB power  per mouse per electrode for groups for ' freq_names{bwii} '= %d \n\n'],  p(3));
                %
                
                %Now do the GLM
                fprintf(1, ['\n\nglm for odor-elicited change in coherence for ' freq_names{bwii} '\n'])
                tbl = table(glm_coh.data',glm_coh.group',glm_coh.perCorr',glm_coh.spm',...
                    'VariableNames',{'delta_coherence','group','perCorr','spm'});
                mdl = fitglm(tbl,'delta_coherence~group+perCorr+spm+group*perCorr*spm'...
                    ,'CategoricalVars',[2,3,4])
                
                
                fprintf(1, ['\n\nRanksum or t-test for delta coherence ' freq_names{bwii} '\n'])
                %Now do the ranksums
                output_data = drgMutiRanksumorTtest(input_data);
                
                fprintf(1, ['\n\n'])
                
            end
            fprintf(1, ['\n\n'])
            
            
            %Now plot the frequency spectra
            
            %Plot the bounded lines
            
            
            
            maxlP=-200000;
            minlP=200000;
            
            for evNo=1:length(eventType)
                
                per_included=0;
                these_dB_per_e=[];
                there_are_NaNs=0;
                
                
                for per_ii=2:-1:1      %performance bins. blue = naive, red = proficient
                    
                    figureNo = figureNo + 1;
                    try
                        close(figureNo)
                    catch
                    end
                    hFig=figure(figureNo);
                    
                    set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
                    
                    set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                    hold on
                    
                    
                    for grNo=max(handles_drgb.drgbchoices.group_no):-1:0
                        
                        if sum((deltaCxy_af_perii_per_mouse==per_ii)&(deltaCxy_af_evNo_per_mouse==evNo)&(deltaCxy_af_group_no_per_mouse==grNo))>0
                            
                            mean_deltaCxy=[];
                            mean_deltaCxy=mean(deltaCxy_af_per_mouse((deltaCxy_af_perii_per_mouse==per_ii)&(deltaCxy_af_evNo_per_mouse==evNo)&(deltaCxy_af_group_no_per_mouse==grNo),:),1);
                            
                            CI=[];
                            CI = bootci(1000, {@mean, deltaCxy_af_per_mouse((deltaCxy_af_perii_per_mouse==per_ii)&(deltaCxy_af_evNo_per_mouse==evNo)&(deltaCxy_af_group_no_per_mouse==grNo),:)})';
                            maxlP=max([maxlP max(CI(:))]);
                            minlP=min([minlP min(CI(:))]);
                            CI(:,1)= mean_deltaCxy'-CI(:,1);
                            CI(:,2)=CI(:,2)- mean_deltaCxy';
                            
                            
                            switch grNo
                                case 1
                                    [hlCR, hpCR] = boundedline(frequency',mean_deltaCxy', CI, 'g');
                                case 2
                                    [hlCR, hpCR] = boundedline(frequency',mean_deltaCxy', CI, 'cmap',[0.3010 0.7450 0.9330]);
                                case 3
                                    [hlCR, hpCR] = boundedline(frequency',mean_deltaCxy', CI, 'y');
                            end
                            
                            %Save the coherence spectrum
                            handles_out.dcohaf_ii=handles_out.dcohaf_ii+1;
                            handles_out.dcohaf_values(handles_out.dcohaf_ii).evNo=evNo;
                            handles_out.dcohaf_values(handles_out.dcohaf_ii).per_ii=per_ii;
                            handles_out.dcohaf_values(handles_out.dcohaf_ii).groupNo=grNo;
                            handles_out.dcohaf_values(handles_out.dcohaf_ii).dcohaf=mean_deltaCxy;
                            handles_out.frequency=frequency;
                            
                            %Save the mean per mouse
                            these_mice=deltaCxy_af_mouseNo_per_mouse((deltaCxy_af_perii_per_mouse==per_ii)&(deltaCxy_af_evNo_per_mouse==evNo)&(deltaCxy_af_group_no_per_mouse==grNo));
                            these_deltaCxy_af_per_mouse=deltaCxy_af_per_mouse((deltaCxy_af_perii_per_mouse==per_ii)&(deltaCxy_af_evNo_per_mouse==evNo)&(deltaCxy_af_group_no_per_mouse==grNo),:);
                            handles_out.dcohaf_values(handles_out.dcohaf_ii).noMice=0;
                            handles_out.dcohaf_values(handles_out.dcohaf_ii).dcoh_per_mouse=[];
                            for iiMice=min(these_mice):max(these_mice)
                                if sum(these_mice==iiMice)>0
                                    handles_out.dcohaf_values(handles_out.dcohaf_ii).noMice=handles_out.dcohaf_values(handles_out.dcohaf_ii).noMice+1;
                                    handles_out.dcohaf_values(handles_out.dcohaf_ii).mouseNo(handles_out.dcohaf_values(handles_out.dcohaf_ii).noMice)=iiMice;
                                    handles_out.dcohaf_values(handles_out.dcohaf_ii).dcoh_per_mouse(handles_out.dcohaf_values(handles_out.dcohaf_ii).noMice,:)=mean(these_deltaCxy_af_per_mouse((these_mice==iiMice),:));
                                end
                            end
                        end
                        
                        
                    end
                    
                    title(['delta coherence per mouse, per electrode during odor ' prof_naive_leg{per_ii} ' ' evTypeLabels{evNo}])
                    xlabel('Frequency (Hz')
                    ylabel('delta coherence')
                end
                
            end
            
            %Set the same ylim for all figures
            maxyl=maxlP+0.1*(maxlP-minlP);
            minyl=minlP-0.1*(maxlP-minlP);
            
            fNo=figureNo-4;
            for evNo=1:length(eventType)
                for per_ii=2:-1:1
                    fNo=fNo+1;
                    hFig=figure(fNo);
                    ylim([minyl maxyl])
                end
            end
            
            
            %         %Now plot the average per mouse LFP power
            %         for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            %             %Plot the average
            %
            %             figureNo=figureNo+1;
            %             try
            %                 close(figureNo)
            %             catch
            %             end
            %             hFig=figure(figureNo);
            %             set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            %
            %
            %
            %             set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            %             hold on
            %
            %             per_mouse_data_delta_dB=[];
            %             per_mouse_prof_naive=[];
            %             per_mouse_events=[];
            %             per_mouse_groups=[];
            %
            %             bar_lab_loc=[];
            %             no_ev_labels=0;
            %             ii_gr_included=0;
            %             all_bar=[];
            %
            %
            %
            % %             fprintf(1, ['\n\n'])
            % %             fprintf(1, ['ANOVAN for delta dB power per mouse, electrode average\n\n'])
            %
            %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
            %
            %                 include_group=0;
            %
            %                 for evNo=1:length(eventType)
            %
            %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
            %
            %                         %                         bar_offset=21-evNo*3+(2-per_ii);
            %
            %                         if sum(eventType==3)>0
            %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(evNo-1);
            %                         else
            %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
            %                         end
            %
            %                         these_offsets(per_ii)=bar_offset;
            %                         %I added this
            % %                         all_bar = [all_bar mean(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo)];
            %
            %                         %Compute per mouse avearge for this group
            %                         no_mice_for_this_group=0;
            %                         each_mouse_average_delta_dB=[];
            %                         for mouseNo=1:max(deltaCxy_mouseNo_per_mouse)
            %                             if sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_mouseNo_per_mouse==mouseNo)&(deltaCxy_group_no_per_mouse==grNo))>0
            %                                 no_mice_for_this_group=no_mice_for_this_group+1;
            %                                 each_mouse_average_delta_dB(no_mice_for_this_group)=mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_mouseNo_per_mouse==mouseNo)&(deltaCxy_group_no_per_mouse==grNo)));
            %                             end
            %                         end
            %
            %                         if no_mice_for_this_group>0
            %
            %                             include_group=1;
            %
            %                             if per_ii==1
            %                                 bar(bar_offset,mean(each_mouse_average_delta_dB),'r','LineWidth', 3)
            %                             else
            %                                 bar(bar_offset,mean(each_mouse_average_delta_dB),'b','LineWidth', 3)
            %                             end
            %                             all_bar = [all_bar mean(each_mouse_average_delta_dB)];
            %
            %                             %In the future add lines linking the points
            %                             plot((bar_offset)*ones(1,no_mice_for_this_group),each_mouse_average_delta_dB,'o',...
            %                                 'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
            %
            %                             %Average and CI
            %                             plot(bar_offset,mean(each_mouse_average_delta_dB),'ok','LineWidth', 3)
            %                             if no_mice_for_this_group>2
            %                                 CI = bootci(1000, {@mean, each_mouse_average_delta_dB},'type','cper');
            %                                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            %                             end
            %
            %
            %                             %Save data for anovan
            %                             per_mouse_data_delta_dB=[per_mouse_data_delta_dB each_mouse_average_delta_dB];
            %                             per_mouse_prof_naive=[per_mouse_prof_naive per_ii*ones(1,no_mice_for_this_group)];
            %                             per_mouse_events=[per_mouse_events evNo*ones(1,no_mice_for_this_group)];
            %                             per_mouse_groups=[per_mouse_groups grNo*ones(1,no_mice_for_this_group)];
            %                         end
            %                     end
            %                 end
            %                 if include_group==1
            %                     ii_gr_included=ii_gr_included+1;
            %                     groups_included(ii_gr_included)=grNo;
            %                 end
            %             end
            % %             if sum(eventType==3)>0
            % %                 title([freq_names{bwii} ' delta coherence per mouse, electrode average'])
            % %             else
            % %                 title([freq_names{bwii} ' delta coherence per mouse, electrode avearage concentrations two steps appart'])
            % %             end
            %
            %              if sum(eventType==3)>0
            %                 title([freq_names{bwii} ' coherence during odor per mouse, electrode average'])
            %             else
            %                 title([freq_names{bwii} ' coherence during odor per mouse, electrode avearage concentrations two steps appart'])
            %              end
            %             % I added
            %             ylim([-0.5 0.5])
            %             %Annotations identifying groups
            %             x_interval=0.8/ii_gr_included;
            %             for ii=1:ii_gr_included
            %                 annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.8 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            %             end
            %
            %             %Proficient/Naive annotations
            %             annotation('textbox',[0.15 0.70 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            %             annotation('textbox',[0.15 0.65 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            %
            %             %x labels
            %             to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
            %             sorted_A=sortrows(to_sort);
            %             sorted_bar_lab_loc=sorted_A(:,1);
            %             for ii=1:length(bar_lab_loc)
            %                 sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
            %             end
            % %             xticks(sorted_bar_lab_loc)
            % %             xticklabels(sorted_bar_labels)
            %
            %             if sum(eventType==3)==0
            %                 xlabel('Concentration (%)')
            %             end
            %
            %
            % %             ylabel('Delta coherence')
            %             ylabel('Coherence')
            %
            %
            %
            % %             %Calculate anovan for inteaction
            % %
            % %
            % %             [p,tbl,stats]=anovan(per_mouse_data_delta_dB,{per_mouse_prof_naive per_mouse_events per_mouse_groups},'varnames',{'proficient_vs_naive','events','groups'},'display','off');
            % %             fprintf(1, ['p value for anovan delta dB histogram per mouse, electrode avearage for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
            % %             fprintf(1, ['p value for anovan delta dB histogram  per mouse, electrode avearage for events ' freq_names{bwii} '= %d \n'],  p(2));
            % %             fprintf(1, ['p value for anovan delta dB histogram  per mouse, electrode avearage for groups ' freq_names{bwii} '= %d \n\n'],  p(2));
            % %
            % %
            %
            %         end
            
            
            
            %         %Display cumulative histograms for delta coherence for average per electrode per mouse and do ranksum
            %         pvals=[];
            %         ranksum_deltaCxy=[];
            %         if sum(eventType==3)>0
            %             for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            %
            %                 %fprintf for ranksums
            % %                 fprintf(1, ['Ranksum or t-test p values for delta coherence for ' freq_names{bwii} '\n'])
            %
            %                 ii_rank=0;
            %                 glm_coh=[];
            %                 glm_ii=0;
            %                 for evNo=1:length(eventType)
            %                     %Plot the average
            %
            %                     figureNo=figureNo+1;
            %                     try
            %                         close(figureNo)
            %                     catch
            %                     end
            %                     hFig=figure(figureNo);
            %
            %                     set(hFig, 'units','normalized','position',[.2 .2 .6 .6])
            %
            %                     set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            %                     hold on
            %                     for grNo=1:max(handles_drgb.drgbchoices.group_no)
            %
            %                         include_group=0;
            %
            %
            %                         for per_ii=1:2      %performance bins. blue = naive, red = proficient
            %
            %                             if sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))>0
            %
            %                                 include_group=1;
            %
            %                                 [f_aic,x_aic] = drg_ecdf(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
            %                                 if per_ii==1
            %                                     switch grNo
            %                                         case 1
            %                                             p1=plot(x_aic,f_aic,'Color',[0 0 1],'LineWidth',3);
            %                                         case 2
            %                                             p2=plot(x_aic,f_aic,'Color',[1 0 0],'LineWidth',3);
            %                                         case 3
            %                                             p3=plot(x_aic,f_aic,'Color',[0 1 0],'LineWidth',3);
            %                                     end
            %                                 else
            %                                     switch grNo
            %                                         case 1
            %                                             p4=plot(x_aic,f_aic,'Color',[0.7 0.7 1],'LineWidth',3);
            %                                         case 2
            %                                             p5=plot(x_aic,f_aic,'Color',[1 0.7 0.7],'LineWidth',3);
            %                                         case 3
            %                                             p6=plot(x_aic,f_aic,'Color',[0.7 1 0.7],'LineWidth',3);
            %                                     end
            %                                 end
            %
            %                                 %Compute and plot per mouse avearge for this group
            %                                 no_mice_for_this_group=0;
            %                                 each_mouse_average_delta_Cxy=[];
            %                                 for mouseNo=1:max(deltaCxy_mouseNo_per_mouse)
            %                                     if sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_mouseNo_per_mouse==mouseNo)&(deltaCxy_group_no_per_mouse==grNo))>0
            %                                         no_mice_for_this_group=no_mice_for_this_group+1;
            %                                         each_mouse_average_delta_Cxy(no_mice_for_this_group)=mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_mouseNo_per_mouse==mouseNo)&(deltaCxy_group_no_per_mouse==grNo)));
            %                                     end
            %                                 end
            %
            %                                 for jj=1:length(each_mouse_average_delta_Cxy)
            %
            %                                     xii_below=find(x_aic<each_mouse_average_delta_Cxy(jj),1,'last');
            %                                     xii_above=find(x_aic>each_mouse_average_delta_Cxy(jj),1,'first');
            %
            %                                     slope=(f_aic(xii_above)-f_aic(xii_below))/(x_aic(xii_above)-x_aic(xii_below));
            %                                     intercept=f_aic(xii_above)-slope*x_aic(xii_above);
            %
            %                                     this_f=slope*each_mouse_average_delta_Cxy(jj)+intercept;
            %
            %                                     switch grNo
            %                                         case 1
            %                                             if per_ii==1
            %                                                 plot(each_mouse_average_delta_Cxy(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
            %                                             else
            %                                                 plot(each_mouse_average_delta_Cxy(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
            %                                             end
            %                                         case 2
            %                                             if per_ii==1
            %                                                 plot(each_mouse_average_delta_Cxy(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
            %                                             else
            %                                                 plot(each_mouse_average_delta_Cxy(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
            %                                             end
            %                                         case 3
            %                                             if per_ii==1
            %                                                 plot(each_mouse_average_delta_Cxy(jj),this_f,'o','MarkerFace',[0 1 0],'MarkerEdge',[0 1 0],'MarkerSize',10)
            %                                             else
            %                                                 plot(each_mouse_average_delta_Cxy(jj),this_f,'o','MarkerFace',[0.7 1 0.7],'MarkerEdge',[0.7 1 0.7],'MarkerSize',10)
            %                                             end
            %                                     end
            %
            %                                 end
            %
            %                                 %Save data for ranksum
            %                                 ii_rank=ii_rank+1;
            %                                 ranksum_deltaCxy(ii_rank).deltaCxy=deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo));
            %                                 ranksum_deltaCxy(ii_rank).per_ii=per_ii;
            %                                 ranksum_deltaCxy(ii_rank).grNo=grNo;
            %                                 ranksum_deltaCxy(ii_rank).evNo=evNo;
            %                             end
            %                         end
            %
            %                         if include_group==1
            %                             ii_gr_included=ii_gr_included+1;
            %                             groups_included(ii_gr_included)=grNo;
            %                         end
            %                     end
            %                     title([freq_names{bwii} ' delta coherence per electrode (per mouse) for ' evTypeLabels{evNo}])
            % %                     title([freq_names{bwii} ' coherence per electrode (per mouse) for ' evTypeLabels{evNo}])
            %
            %                     legend([p1 p2 p3 p4 p5 p6],{[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' proficient'] ,[handles_drgb.drgbchoices.group_no_names{3} ' proficient'],...
            %                         [handles_drgb.drgbchoices.group_no_names{1} ' naive'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'],[handles_drgb.drgbchoices.group_no_names{3} ' naive']})
            %
            %                     xlabel('Delta coherence')
            % %                      xlabel('Coherence')
            %                     ylabel('Cumulative probability')
            %
            %
            %                     prof_naive_leg{1}='Proficient';
            %                     prof_naive_leg{2}='Naive';
            %
            %                 end
            %
            %                 %Now do the ranksums
            %                 prof_naive_leg{1}='Proficient';
            %                 prof_naive_leg{2}='Naive';
            %
            %                 input_data=[];
            %                 for ii=1:ii_rank
            %                     input_data(ii).data=ranksum_deltaCxy(ii).deltaCxy;
            %                     input_data(ii).description=[handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(ii).grNo} ' ' evTypeLabels{ranksum_deltaCxy(ii).evNo} ' ' prof_naive_leg{ranksum_deltaCxy(ii).per_ii}];
            %                     input_data(ii).per_ii=ranksum_deltaCxy(ii).per_ii;
            %                     input_data(ii).grNo=ranksum_deltaCxy(ii).grNo;
            %                     input_data(ii).evNo=ranksum_deltaCxy(ii).evNo;
            %
            %                     glm_coh.data(glm_ii+1:glm_ii+length(ranksum_deltaCxy(ii).deltaCxy))=ranksum_deltaCxy(ii).deltaCxy;
            %                     glm_coh.group(glm_ii+1:glm_ii+length(ranksum_deltaCxy(ii).deltaCxy))=ranksum_deltaCxy(ii).grNo;
            %                     glm_coh.perCorr(glm_ii+1:glm_ii+length(ranksum_deltaCxy(ii).deltaCxy))=ranksum_deltaCxy(ii).per_ii;
            %                     glm_coh.event(glm_ii+1:glm_ii+length(ranksum_deltaCxy(ii).deltaCxy))=ranksum_deltaCxy(ii).evNo;
            %                     glm_ii=glm_ii+length(ranksum_deltaCxy(ii).deltaCxy);
            %                 end
            %
            %                 %Perform the glm
            % %                 fprintf(1, ['glm for delta coherence for each electrode pair calculated per mouse for ' freq_names{bwii} '\n'])
            %                 fprintf(1, ['glm for ' freq_names{bwii} ' coherence during odor for each electrode pair calculated per mouse\n'])
            %
            %
            %                 if sum(glm_coh.group==1)==length(glm_coh.group)
            %                     %There is only one group here (e.g. for Justin's paper we only include
            %                     %forward)
            %                     fprintf(1, ['\n\nglm for odor coherence for ' freq_names{bwii} '\n'])
            %                     tbl = table(glm_coh.data',glm_coh.perCorr',glm_coh.event',...
            %                         'VariableNames',{'delta_coherence','perCorr','event'});
            %                     mdl = fitglm(tbl,'delta_coherence~perCorr+event+perCorr*event'...
            %                         ,'CategoricalVars',[2,3])
            %                 else
            %
            %                     fprintf(1, ['\n\nglm for odor coherence for ' freq_names{bwii} '\n'])
            %                     tbl = table(glm_coh.data',glm_coh.group',glm_coh.perCorr',glm_coh.event',...
            %                         'VariableNames',{'delta_coherence','group','perCorr','event'});
            %                     mdl = fitglm(tbl,'delta_coherence~group+perCorr+event+perCorr*group*event'...
            %                         ,'CategoricalVars',[2,3,4])
            %                 end
            %
            %                 %Do the ranksum/t-test
            %                 fprintf(1, ['\n\nRanksum or t-test p values for odor coherence for each electrode pair calculated per mouse for ' freq_names{bwii} '\n'])
            %                 [output_data] = drgMutiRanksumorTtest(input_data);
            %
            %
            % %                 for ii=1:ii_rank
            % %                     for jj=ii+1:ii_rank
            % %                         [p, r_or_t]=drg_ranksum_or_ttest(ranksum_deltaCxy(ii).deltaCxy,ranksum_deltaCxy(jj).deltaCxy);
            % %                         if r_or_t==0
            % %                             fprintf(1, ['p value ranksum for ' handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(ii).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(ii).per_ii} ' ' evTypeLabels{ranksum_deltaCxy(ii).evNo} ' vs ' ...
            % %                                 handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(jj).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(jj).per_ii} ' ' evTypeLabels{ranksum_deltaCxy(jj).evNo} ' =  %d\n'],p)
            % %                         else
            % %                             fprintf(1, ['p value t-test for ' handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(ii).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(ii).per_ii} ' ' evTypeLabels{ranksum_deltaCxy(ii).evNo} ' vs ' ...
            % %                                 handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(jj).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(jj).per_ii} ' ' evTypeLabels{ranksum_deltaCxy(jj).evNo} ' =  %d\n'],p)
            % %                         end
            % %
            % %                         pvals=[pvals p];
            % %                     end
            % %                 end
            %
            %                 fprintf(1, ['\n\n'])
            %
            %             end
            %         end
            %         fprintf(1, ['\n\n'])
            %         pFDR = drsFDRpval(pvals);
            %         fprintf(1, ['pFDR = %d \n\n'],pFDR)
            %         fprintf(1, ['\n\n'])
            
            
            %Display cumulative histograms for  coherence before odor on for average per electrode per mouse and do ranksum
            %         pvals=[];
            %         ranksum_Cxy=[];
            %         if sum(eventType==3)>0
            %             for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            %
            %                 %fprintf for ranksums
            % %                 fprintf(1, ['Ranksum or t-test p values for coherence before odor for ' freq_names{bwii} '\n'])
            %
            %                 ii_rank=0;
            %                 glm_coh_pre=[];
            %                 glm_ii=0;
            %
            %                 figureNo=figureNo+1;
            %                 try
            %                     close(figureNo)
            %                 catch
            %                 end
            %                 hFig=figure(figureNo);
            %
            %                 set(hFig, 'units','normalized','position',[.2 .2 .6 .6])
            %
            %                 set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            %                 hold on
            %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
            %
            %                     include_group=0;
            %
            %
            %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
            %
            %                         if sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))>0
            %
            %                             include_group=1;
            %
            % %                             [f_aic,x_aic] = drg_ecdf(Cxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
            % %
            % %                             switch grNo
            % %                             case 1
            % %                             if per_ii==1
            % % %                                 if grNo==1
            % %
            % %                                     plot(x_aic,f_aic,'Color',[0 0 1],'LineWidth',3)
            % %                                 else
            % %
            % %                                     plot(x_aic,f_aic,'Color',[0.7 0.7 1],'LineWidth',3)
            % %                                 end
            % %                                 case 2
            % %                                 if per_ii==1
            % %
            % %                                     plot(x_aic,f_aic,'Color',[1 0 0])
            % %                                 else
            % %
            % %                                     plot(x_aic,f_aic,'Color',[1 0.7 0.7])
            % %                                 end
            % %                                  case 3
            % %                                 if per_ii==1
            % %
            % %                                     plot(x_aic,f_aic,'Color',[0 1 0])
            % %                                 else
            % %
            % %                                     plot(x_aic,f_aic,'Color',[0.7 1 0.7])
            % %                                 end
            % %                             end
            %
            %                             [f_aic,x_aic] = drg_ecdf(Cxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
            %
            %                             if per_ii==1
            %                                 switch  grNo
            %                                     case 1
            %
            %                                     p1=plot(x_aic,f_aic,'Color',[0 0 1],'LineWidth',3);
            %
            %                                     case 2
            %
            %                                     p2=plot(x_aic,f_aic,'Color',[1 0 0],'LineWidth',3)
            %
            %                                     case 3
            %
            %                                     p3=plot(x_aic,f_aic,'Color',[0 1 0],'LineWidth',3)
            %                                 end
            %                             else
            %                                 switch  grNo
            %                                     case 1
            %
            %                                     p4=plot(x_aic,f_aic,'Color',[0.7 0.7 1],'LineWidth',3)
            %
            %                                     case 2
            %
            %                                     p5=plot(x_aic,f_aic,'Color',[1 0.7 0.7],'LineWidth',3)
            %
            %                                     case 3
            %
            %                                     p6=plot(x_aic,f_aic,'Color',[0.7 1 0.7],'LineWidth',3)
            %                                 end
            % %                                   if grNo==1
            % %                                     [f_aic,x_aic] = drg_ecdf(Cxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
            % %                                     plot(x_aic,f_aic,'Color',[0 1 0])
            % %                                 else
            % %                                     [f_aic,x_aic] = drg_ecdf(Cxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
            % %                                     plot(x_aic,f_aic,'Color',[0 1 0])
            % %                                 end
            %                             end
            %
            %                               %Compute and plot per mouse avearge for this group
            %                                 no_mice_for_this_group=0;
            %                                 each_mouse_average_Cxy_per_mouse=[];
            %                                 for mouseNo=1:max(deltaCxy_mouseNo_per_mouse)
            %                                     if sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_mouseNo_per_mouse==mouseNo)&(deltaCxy_group_no_per_mouse==grNo))>0
            %                                         no_mice_for_this_group=no_mice_for_this_group+1;
            %                                         each_mouse_average_Cxy_per_mouse(no_mice_for_this_group)=mean(Cxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_mouseNo_per_mouse==mouseNo)&(deltaCxy_group_no_per_mouse==grNo)));
            %                                     end
            %                                 end
            %
            %                                 for jj=1:length(each_mouse_average_Cxy_per_mouse)
            %
            %                                     xii_below=find(x_aic<each_mouse_average_Cxy_per_mouse(jj),1,'last');
            %                                     xii_above=find(x_aic>each_mouse_average_Cxy_per_mouse(jj),1,'first');
            %
            %                                     slope=(f_aic(xii_above)-f_aic(xii_below))/(x_aic(xii_above)-x_aic(xii_below));
            %                                     intercept=f_aic(xii_above)-slope*x_aic(xii_above);
            %
            %                                     this_f=slope*each_mouse_average_Cxy_per_mouse(jj)+intercept;
            %
            %                                     switch grNo
            %                                         case 1
            %                                             if per_ii==1
            %                                                 plot(each_mouse_average_Cxy_per_mouse(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
            %                                             else
            %                                                 plot(each_mouse_average_Cxy_per_mouse(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
            %                                             end
            %                                         case 2
            %                                             if per_ii==1
            %                                                 plot(each_mouse_average_Cxy_per_mouse(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
            %                                             else
            %                                                 plot(each_mouse_average_Cxy_per_mouse(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
            %                                             end
            %                                         case 3
            %                                             if per_ii==1
            %                                                 plot(each_mouse_average_Cxy_per_mouse(jj),this_f,'o','MarkerFace',[0 1 0],'MarkerEdge',[0 1 0],'MarkerSize',10)
            %                                             else
            %                                                 plot(each_mouse_average_Cxy_per_mouse(jj),this_f,'o','MarkerFace',[0.7 1 0.7],'MarkerEdge',[0.7 1 0.7],'MarkerSize',10)
            %                                             end
            %                                     end
            %
            %                                 end
            %
            %                             %Save data for ranksum
            %                             ii_rank=ii_rank+1;
            %                             ranksum_Cxy(ii_rank).Cxy=Cxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo));
            %                             ranksum_Cxy(ii_rank).per_ii=per_ii;
            %                             ranksum_Cxy(ii_rank).grNo=grNo;
            %
            %                         end
            %                     end
            %
            %                     if include_group==1
            %                         ii_gr_included=ii_gr_included+1;
            %                         groups_included(ii_gr_included)=grNo;
            %                     end
            %                 end
            %                 title([freq_names{bwii} ' coherence before odor per electrode (per mouse)'])
            %
            %                legend([p1 p2 p3 p4 p5 p6],{[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' proficient'] ,[handles_drgb.drgbchoices.group_no_names{3} ' proficient'],...
            %                         [handles_drgb.drgbchoices.group_no_names{1} ' naive'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'],[handles_drgb.drgbchoices.group_no_names{3} ' naive']})
            %
            % %                 legend([handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive']...
            % %                     ,[handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'],...
            % %                  [handles_drgb.drgbchoices.group_no_names{3} ' proficient'],[handles_drgb.drgbchoices.group_no_names{3} ' naive'])
            %
            %                 xlabel('coherence')
            %                 ylabel('Cumulative probability')
            %
            %
            %                 prof_naive_leg{1}='Proficient';
            %                 prof_naive_leg{2}='Naive';
            %
            %
            %
            %                 input_data=[];
            %                 for ii=1:ii_rank
            %                     input_data(ii).data=ranksum_Cxy(ii).Cxy;
            %                     input_data(ii).description=[handles_drgb.drgbchoices.group_no_names{ranksum_Cxy(ii).grNo} ' '  prof_naive_leg{ranksum_Cxy(ii).per_ii}];
            %                     input_data(ii).per_ii=ranksum_Cxy(ii).per_ii;
            %                     input_data(ii).grNo=ranksum_Cxy(ii).grNo;
            %
            %                     glm_coh_pre.data(glm_ii+1:glm_ii+length(ranksum_Cxy(ii).Cxy))=ranksum_Cxy(ii).Cxy;
            %                     glm_coh_pre.group(glm_ii+1:glm_ii+length(ranksum_Cxy(ii).Cxy))=ranksum_Cxy(ii).grNo;
            %                     glm_coh_pre.perCorr(glm_ii+1:glm_ii+length(ranksum_Cxy(ii).Cxy))=ranksum_Cxy(ii).per_ii;
            %                     glm_ii=glm_ii+length(ranksum_Cxy(ii).Cxy);
            %                 end
            %
            %                 %Perform the glm
            %                 fprintf(1, ['glm for pre-odor coherence for each electrode pair calculated per mouse for ' freq_names{bwii} '\n'])
            %
            %                 if sum(glm_coh.group==1)==length(glm_coh.group)
            %                     %There is only one group here (e.g. for Justin's paper we only include
            %                     %forward)
            %                     fprintf(1, ['\n\nglm for pre-odor coherence for ' freq_names{bwii} '\n'])
            %                     tbl = table(glm_coh.data',glm_coh.perCorr',...
            %                         'VariableNames',{'pre_odor_coherence','perCorr'});
            %                     mdl = fitglm(tbl,'pre_odor_coherence~perCorr'...
            %                         ,'CategoricalVars',[2])
            %                 else
            %
            %                     fprintf(1, ['\n\nglm for pre-odor coherence for ' freq_names{bwii} '\n'])
            %                     tbl = table(glm_coh.data',glm_coh.group',glm_coh.perCorr',...
            %                         'VariableNames',{'pre_odor_coherence','group','perCorr'});
            %                     mdl = fitglm(tbl,'pre_odor_coherence~group+perCorr+perCorr*group'...
            %                         ,'CategoricalVars',[2,3])
            %                 end
            %
            %                 %Do the ranksum/t-test
            %                 fprintf(1, ['\n\nRanksum or t-test p values for pre-odor coherence for each electrode pair calculated per mouse for ' freq_names{bwii} '\n'])
            %                 [output_data] = drgMutiRanksumorTtest(input_data);
            %
            %
            % %                 for ii=1:ii_rank
            % %                     for jj=ii+1:ii_rank
            % %                         [p, r_or_t]=drg_ranksum_or_ttest(ranksum_Cxy(ii).Cxy,ranksum_Cxy(jj).Cxy);
            % %                         if r_or_t==0
            % %                             fprintf(1, ['p value ranksum for ' handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(ii).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(ii).per_ii}  ' vs ' ...
            % %                                 handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(jj).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(jj).per_ii}  ' =  %d\n'],p)
            % %                         else
            % %                             fprintf(1, ['p value t-test for ' handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(ii).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(ii).per_ii}  ' vs ' ...
            % %                                 handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(jj).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(jj).per_ii}  ' =  %d\n'],p)
            % %                         end
            % %                         pvals=[pvals p];
            % %                     end
            % %                 end
            % %
            %                 fprintf(1, ['\n\n'])
            %
            %             end
            %         end
            % %         fprintf(1, ['\n\n'])
            % %         pFDR = drsFDRpval(pvals);
            % %         fprintf(1, ['pFDR = %d \n\n'],pFDR)
            %         fprintf(1, ['\n\n'])
            %
            %Display auROC
            edges=[-0.3:0.05:0.5];
            rand_offset=0.8;
            coh_auROC_per_mouse=[];
            coh_group_no_per_mouse=[];
            
            for bwii=1:no_bandwidths    %for different bandwidths
                
                ii_roc=0;
                roc_data=[];
                glm_roc=[];
                glm_roc_ii=0;
                
                
                
                %Display the  auROC
                figureNo=figureNo+1;
                try
                    close(figureNo)
                catch
                end
                hFig=figure(figureNo);
                
                set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
                
                set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                hold on
                
                
                
                bar_offset=0;
                
                for per_ii=2:-1:1
                    for grNo=1:max(handles_drgb.drgbchoices.group_no)
                        
                        
                        if sum((ROCper_ii==per_ii)&(ROCbwii==bwii)&(ROCgroups==grNo))>0
                            
                            bar_offset=bar_offset+1;
                            
                            switch grNo
                                case 1
                                    bar(bar_offset,mean(auROC((ROCper_ii==per_ii)&(ROCbwii==bwii)&(ROCgroups==grNo))),'g','LineWidth', 3,'EdgeColor','none')
                                case 2
                                    bar(bar_offset,mean(auROC((ROCper_ii==per_ii)&(ROCbwii==bwii)&(ROCgroups==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                                case 3
                                    bar(bar_offset,mean(auROC((ROCper_ii==per_ii)&(ROCbwii==bwii)&(ROCgroups==grNo))),'y','LineWidth', 3,'EdgeColor','none')
                            end
                            
                            handles_out.auc_ii=handles_out.auc_ii+1;
                            handles_out.auc_values(handles_out.auc_ii).pacii=bwii;
                            handles_out.auc_values(handles_out.auc_ii).per_ii=per_ii;
                            handles_out.auc_values(handles_out.auc_ii).groupNo=grNo;
                            handles_out.auc_values(handles_out.auc_ii).auc_coh=mean(auROC((ROCper_ii==per_ii)&(ROCbwii==bwii)&(ROCgroups==grNo)));
                            
                            %Save the mean per mouse
                            these_mice=ROCmouse((ROCper_ii==per_ii)&(ROCbwii==bwii)&(ROCgroups==grNo));
                            handles_out.auc_values(handles_out.auc_ii).noMice=0;
                            for iiMice=min(these_mice):max(these_mice)
                                if sum(these_mice==iiMice)>0
                                    handles_out.auc_values(handles_out.auc_ii).noMice=handles_out.auc_values(handles_out.auc_ii).noMice+1;
                                    handles_out.auc_values(handles_out.auc_ii).mouseNo(handles_out.auc_values(handles_out.auc_ii).noMice)=iiMice;
                                    handles_out.auc_values(handles_out.auc_ii).auROC_per_mouse(handles_out.auc_values(handles_out.auc_ii).noMice)=mean(auROC((ROCmouse==iiMice)&(ROCper_ii==per_ii)&(ROCbwii==bwii)&(ROCgroups==grNo)));
                                    coh_auROC_per_mouse(iiMice,bwii,per_ii)=mean(auROC((ROCmouse==iiMice)&(ROCper_ii==per_ii)&(ROCbwii==bwii)&(ROCgroups==grNo)));
                                    coh_group_no_per_mouse(iiMice)=grNo;
                                end
                            end
                            
                            %Violin plot
                            [mean_out, CIout]=drgViolinPoint(auROC((ROCper_ii==per_ii)&(ROCbwii==bwii)&(ROCgroups==grNo))...
                                ,edges,bar_offset,rand_offset,'k','k',3);
                            
                            
                            %Enter data for glm for peak and trough only
                            these_data=auROC((ROCper_ii==per_ii)&(ROCbwii==bwii)&(ROCgroups==grNo));
                            glm_roc.data(glm_roc_ii+1:glm_roc_ii+length(these_data))=these_data;
                            glm_roc.group(glm_roc_ii+1:glm_roc_ii+length(these_data))=grNo;
                            glm_roc.perCorr(glm_roc_ii+1:glm_roc_ii+length(these_data))=per_ii;
                            glm_roc_ii=glm_roc_ii+length(these_data);
                            
                            %Enter the data for t-test/ranksum
                            ii_roc=ii_roc+1;
                            roc_data(ii_roc).data=these_data;
                            roc_data(ii_roc).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' prof_naive_leg{per_ii}];
                            
                        end
                        
                        
                    end
                    
                    
                    bar_offset=bar_offset+1;
                end
                
                title(['auROC per mouse, per electrode for ' freq_names{bwii} ' coherence'])
                
                
                
                ylabel('auROC')
                
                pffft=1;
                
                
                
                %Do glm for auROC coherence
                fprintf(1, ['\n\nglm for auROC coherence for ' freq_names{bwii} '\n'])
                tbl = table(glm_roc.data',glm_roc.group',glm_roc.perCorr',...
                    'VariableNames',{'Peak_wave','group','perCorr'});
                mdl = fitglm(tbl,'Peak_wave~group+perCorr+group*perCorr'...
                    ,'CategoricalVars',[2,3])
                
                
                fprintf(1, ['\n\nRanksum or t-test for auROC coherence for ' freq_names{bwii} '\n'])
                %Now do the ranksums
                output_data = drgMutiRanksumorTtest(roc_data);
                
                
                
            end
            
        end
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) handles_pars.output_suffix],'handles_out','coh_auROC_per_mouse','coh_group_no_per_mouse')
        
    case 22
        % Multiclass ROC analysis of LFP power differences for naive and proficient
        % mice for different epochs (concentrations or S+ vs. S-) and different groups. Analyzed per mouse
        %         no_dBs=1;
        %         wave_logP_power=[];
        no_ROCs=0;
        ROCout=[];
        p_vals_ROC=[];
        per_mouse_no_ROCs=0;
        per_mouse_ROCout=[];
        per_mouse_p_vals_ROC=[];
        %         wave_logP_powerEv1=[];
        %         no_Ev1=0;
        %         for evNo=1:length(eventType)
        %             evNo_out(evNo).noWB=0;
        %         end
        %         wave_logP_powerEv1WB=[];
        %         wave_logP_powerEv2WB=[];
        wave_logP_No_per_mouse=0;
        wave_logPpeak_per_mouse=[];
        wave_logP_perii_per_mouse=[];
        wave_logP_evNo_per_mouse=[];
        wave_logP_bwii_per_mouse=[];
        wave_logP_mouseNo_per_mouse=[];
        wave_logP_electrode_per_mouse=[];
        wave_logP_per_mouse=[];
        mouse_included=[];
        
        out_times=handles_drgb.drgb.lfpevpair(1).wave_out_times;
        
        fprintf(1, ['Pairwise ERWA analysis for Justin''s paper\n\n'])
        %         p_vals=[];
        no_files=max(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        
        szpc=size(percent_windows);
        for per_ii=1:szpc(1)
            
            for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
                mouse_has_files=0;
                theseEvNosPerEl=[];
                for evNo=1:length(eventType)
                    for bwii=1:no_bandwidths
                        for elec=1:16
                            theseEvNosPerEl(evNo,bwii,elec).noEv=0;
                        end
                    end
                end
                
                for elec=1:16
                    if sum(which_electrodes==elec)>0
                        
                        theseEvNos=[];
                        for evNo=1:length(eventType)
                            for bwii=1:no_bandwidths
                                theseEvNos(evNo,bwii).noEv=0;
                            end
                        end
                        
                        for fileNo=1:no_files
                            if sum(files==fileNo)>0
                                if handles_drgb.drgbchoices.mouse_no(fileNo)==mouseNo
                                    
                                    
                                    
                                    
                                    %Note: FOr ERWA the reference is
                                    %already subtracted
                                    lfpodNo=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                    
                                    
                                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo).wave_log_P_t_lick_referenced))
                                        
                                        percent_mask=[];
                                        percent_mask=logical((handles_drgb.drgb.lfpevpair(lfpodNo).wave_perCorrERWA>=percent_windows(per_ii,1))...
                                            &(handles_drgb.drgb.lfpevpair(lfpodNo).wave_perCorrERWA<=percent_windows(per_ii,2)));
                                        
                                        
                                        for evNo=1:length(eventType)
                                            
                                            noWB_for_evNo(evNo)=-1;
                                            
                                            trials_in_event_Ev=[];
                                            trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfpodNo).wave_which_eventERWA(eventType(evNo),:)==1)&percent_mask;
                                            
                                            if (sum(trials_in_event_Ev)>=1)
                                                
                                                
                                                this_wave_logP=zeros(sum(trials_in_event_Ev),length(frequency),length(out_times));
                                                this_wave_logP(:,:,:)=handles_drgb.drgb.lfpevpair(lfpodNo).wave_log_P_t_lick_referenced(trials_in_event_Ev,:,:);
                                                
                                                
                                                %Do per bandwidth analysis
                                                for bwii=1:no_bandwidths
                                                    theseEvNos(evNo,bwii).groupNo(1,theseEvNos(evNo,bwii).noEv+1:theseEvNos(evNo,bwii).noEv+sum(trials_in_event_Ev))=handles_drgb.drgbchoices.group_no(fileNo)*ones(1,sum(trials_in_event_Ev));
                                                    theseEvNos(evNo,bwii).lick_freq(1,theseEvNos(evNo,bwii).noEv+1:theseEvNos(evNo,bwii).noEv+sum(trials_in_event_Ev))=handles_drgb.drgb.lfpevpair(lfpodNo).wave_no_licks_per_trial(1,trials_in_event_Ev)/...
                                                        (handles_drgb.drgbchoices.timeEnd(winNo)-handles_drgb.drgbchoices.timeStart(winNo));
                                                    
                                                    this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                                    
                                                    %Enter the  Ev1
                                                    this_wave_logPEv=zeros(sum(trials_in_event_Ev),length(out_times));
                                                    this_wave_logPEv(:,:)=mean(this_wave_logP(:,this_band,:),2);
                                                    
                                                    theseEvNos(evNo,bwii).this_wave_logPEv(theseEvNos(evNo,bwii).noEv+1:theseEvNos(evNo,bwii).noEv+sum(trials_in_event_Ev),1:length(out_times))=this_wave_logPEv;
                                                    theseEvNos(evNo,bwii).this_wave_logPEvpeak(1,theseEvNos(evNo,bwii).noEv+1:theseEvNos(evNo,bwii).noEv+sum(trials_in_event_Ev))=prctile(this_wave_logPEv',95);
                                                    theseEvNos(evNo,bwii).noEv=theseEvNos(evNo,bwii).noEv+sum(trials_in_event_Ev);
                                                    
                                                    theseEvNosPerEl(evNo,bwii,elec).this_wave_logPEv(theseEvNosPerEl(evNo,bwii,elec).noEv+1:theseEvNosPerEl(evNo,bwii,elec).noEv+sum(trials_in_event_Ev),1:length(out_times))=this_wave_logPEv;
                                                    theseEvNosPerEl(evNo,bwii,elec).this_wave_logPEvpeak(1,theseEvNosPerEl(evNo,bwii,elec).noEv+1:theseEvNosPerEl(evNo,bwii,elec).noEv+sum(trials_in_event_Ev))=prctile(this_wave_logPEv',95);
                                                    theseEvNosPerEl(evNo,bwii,elec).groupNo(1,theseEvNosPerEl(evNo,bwii,elec).noEv+1:theseEvNosPerEl(evNo,bwii,elec).noEv+sum(trials_in_event_Ev))=handles_drgb.drgbchoices.group_no(fileNo)*ones(1,sum(trials_in_event_Ev));
                                                    
                                                    theseEvNosPerEl(evNo,bwii,elec).noEv=theseEvNosPerEl(evNo,bwii,elec).noEv+sum(trials_in_event_Ev);
                                                    
                                                    mouse_has_files=1;
                                                end
                                                
                                                
                                                fprintf(1, ['%d trials in event No %d succesfully processed for file No %d electrode %d\n'],sum(trials_in_event_Ev), evNo,fileNo,elec);
                                                
                                            else
                                                
                                                
                                                fprintf(1, ['%d trials in event No %d fewer than minimum trials per event ' evTypeLabels{evNo} ' for file No %d electrode %d\n'],sum(trials_in_event_Ev), min_trials_per_event,fileNo,elec);
                                                
                                                
                                            end
                                            
                                            
                                        end
                                        
                                        
                                    else
                                        
                                        fprintf(1, ['Empty allPower for file No %d electrode %d\n'],fileNo,elec);
                                        
                                    end
                                    
                                    
                                    
                                end %if mouseNo
                            end
                        end %fileNo
                        
                        if mouse_has_files==1
                            mouse_included(mouseNo)=1;
                            %Calculate per mouse per electrode wave_logP
                            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                                for evNo=1:length(eventType)
                                    for bwii=1:no_bandwidths
                                        if theseEvNos(evNo,bwii).noEv>0
                                            if sum(theseEvNos(evNo,bwii).groupNo==grNo)>0
                                                wave_logP_No_per_mouse=wave_logP_No_per_mouse+1;
                                                wave_logPpeak_per_mouse(wave_logP_No_per_mouse)=mean(theseEvNos(evNo,bwii).this_wave_logPEvpeak(theseEvNos(evNo,bwii).groupNo==grNo));
                                                wave_logP_per_mouse(wave_logP_No_per_mouse,1:length(out_times))=mean(theseEvNos(evNo,bwii).this_wave_logPEv(theseEvNos(evNo,bwii).groupNo==grNo,:),1);
                                                wave_lick_freq_per_mouse(wave_logP_No_per_mouse)=mean(theseEvNos(evNo,bwii).lick_freq(theseEvNos(evNo,bwii).groupNo==grNo));
                                                wave_logP_perii_per_mouse(wave_logP_No_per_mouse)=per_ii;
                                                wave_logP_evNo_per_mouse(wave_logP_No_per_mouse)=evNo;
                                                wave_logP_bwii_per_mouse(wave_logP_No_per_mouse)=bwii;
                                                wave_logP_mouseNo_per_mouse(wave_logP_No_per_mouse)=mouseNo;
                                                wave_logP_electrode_per_mouse(wave_logP_No_per_mouse)=elec;
                                                wave_logP_group_no_per_mouse(wave_logP_No_per_mouse)=grNo;
                                            end
                                        end
                                    end
                                end
                            end
                            
                            %Calculate per electrode ROC
                            can_calculate_auroc=1;
                            if can_calculate_auroc==1
                                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                                    for evNo1=1:length(eventType)
                                        if theseEvNos(evNo1).noEv>0
                                            for evNo2=evNo1+1:length(eventType)
                                                if theseEvNos(evNo2).noEv>0
                                                    for bwii=1:no_bandwidths
                                                        
                                                        
                                                        %Enter Ev1
                                                        trials_in_event_Ev1=length(theseEvNos(evNo1,bwii).this_wave_logPEv(theseEvNos(evNo1).groupNo==grNo));
                                                        this_wave_logPEv1=zeros(trials_in_event_Ev1,1);
                                                        this_wave_logPEv1=theseEvNos(evNo1,bwii).this_wave_logPEv(theseEvNos(evNo1).groupNo==grNo);
                                                        roc_data=[];
                                                        roc_data(1:trials_in_event_Ev1,1)=this_wave_logPEv1;
                                                        roc_data(1:trials_in_event_Ev1,2)=zeros(trials_in_event_Ev1,1);
                                                        
                                                        %Enter Ev2
                                                        trials_in_event_Ev2=length(theseEvNos(evNo2,bwii).this_wave_logPEv(theseEvNos(evNo2).groupNo==grNo));
                                                        total_trials=trials_in_event_Ev1+trials_in_event_Ev2;
                                                        this_wave_logPEv2=zeros(trials_in_event_Ev2,1);
                                                        this_wave_logPEv2=theseEvNos(evNo2,bwii).this_wave_logPEv(theseEvNos(evNo2).groupNo==grNo);
                                                        roc_data(trials_in_event_Ev1+1:total_trials,1)=this_wave_logPEv2;
                                                        roc_data(trials_in_event_Ev1+1:total_trials,2)=ones(trials_in_event_Ev2,1);
                                                        
                                                        
                                                        %Find  per electrode ROC
                                                        if (trials_in_event_Ev1>=5)&(trials_in_event_Ev2>=5)
                                                            no_ROCs=no_ROCs+1;
                                                            roc=roc_calc(roc_data,0,0.05,0);
                                                            ROCout(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
                                                            ROCelec(no_ROCs)=elec;
                                                            ROCgroups(no_ROCs)=grNo;
                                                            ROCmouse(no_ROCs)=mouseNo;
                                                            ROCbandwidth(no_ROCs)=bwii;
                                                            ROCper_ii(no_ROCs)=per_ii;
                                                            ROCEvNo1(no_ROCs)=evNo1;
                                                            ROCEvNo2(no_ROCs)=evNo2;
                                                            if sum(eventType==3)==0
                                                                if ((abs(evNo1-2)<=1)&(abs(evNo2-5)<=1))||((abs(evNo1-5)<=1)&(abs(evNo2-2)<=1))
                                                                    ROC_between(no_ROCs)=1;
                                                                else
                                                                    ROC_between(no_ROCs)=0;
                                                                end
                                                                ROC_neighbor(no_ROCs)=abs(evNo1-evNo2);
                                                            else
                                                                %This is S+/S-,
                                                                %these values are
                                                                %assigned
                                                                %arbitrarily so
                                                                %that plotting
                                                                %auROC works
                                                                ROC_between(no_ROCs)=1;
                                                                ROC_neighbor(no_ROCs)=2;
                                                            end
                                                            
                                                            auROC(no_ROCs)=roc.AUC-0.5;
                                                            p_valROC(no_ROCs)=roc.p;
                                                            p_vals_ROC=[p_vals_ROC roc.p];
                                                            
                                                            %I have this code here to plot the ROC
                                                            if (per_ii==1)&(bwii==4)&(roc.AUC-0.5>0.3)
                                                                show_roc=0;
                                                                if show_roc==1
                                                                    %I have this code here to plot the ROC
                                                                    roc=roc_calc(roc_data,0,0.05,1);
                                                                    
                                                                    %Do the histograms
                                                                    try
                                                                        close(2)
                                                                    catch
                                                                    end
                                                                    figure(2)
                                                                    
                                                                    hold on
                                                                    
                                                                    max_dB=max([max(this_wave_logPEv1) max(this_wave_logPEv2)]);
                                                                    min_dB=min([min(this_wave_logPEv1) min(this_wave_logPEv2)]);
                                                                    
                                                                    edges=[min_dB-0.1*(max_dB-min_dB):(max_dB-min_dB)/20:max_dB+0.1*(max_dB-min_dB)];
                                                                    histogram(this_wave_logPEv1,edges,'FaceColor','b','EdgeColor','b')
                                                                    histogram(this_wave_logPEv2,edges,'FaceColor','r','EdgeColor','r')
                                                                    xlabel('delta power dB')
                                                                    title(['Histogram for conentrations ' num2str(concs2(evNo1)) ' and ' num2str(concs2(evNo2))])
                                                                    pffft=1;
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            
                        else
                            mouse_included(mouseNo)=0;
                        end
                        
                        
                    end
                    
                end
                
                %Calculate per mouse ROC
                if mouse_has_files==1
                    for grNo=1:max(handles_drgb.drgbchoices.group_no)
                        for evNo1=1:length(eventType)
                            for evNo2=evNo1+1:length(eventType)
                                
                                
                                for bwii=1:no_bandwidths
                                    
                                    %Enter Ev1
                                    trials_in_event_Ev1=length(theseEvNosPerEl(evNo1,bwii,which_electrodes(1)).this_wave_logPEv(theseEvNosPerEl(evNo1,bwii,which_electrodes(1)).groupNo==grNo));
                                    
                                    this_wave_logPEv1=zeros(trials_in_event_Ev1,1);
                                    
                                    %                                     for elec=which_electrodes
                                    %                                         this_wave_logPEv1=this_wave_logPEv1+(theseEvNosPerEl(evNo1,bwii,elec).this_wave_logPEv((theseEvNosPerEl(evNo1,bwii,elec).groupNo==grNo))')/length(which_electrodes);
                                    %                                     end
                                    %
                                    no_elects=0;
                                    for elec=which_electrodes
                                        if length(this_wave_logPEv1)==length(theseEvNosPerEl(evNo1,bwii,elec).groupNo)
                                            this_wave_logPEv1=this_wave_logPEv1+(theseEvNosPerEl(evNo1,bwii,elec).this_wave_logPEv((theseEvNosPerEl(evNo1,bwii,elec).groupNo==grNo))');
                                            no_elects=no_elects+1;
                                        end
                                    end
                                    this_wave_logPEv1=this_wave_logPEv1/no_elects;
                                    
                                    roc_data=[];
                                    roc_data(1:trials_in_event_Ev1,1)=this_wave_logPEv1;
                                    roc_data(1:trials_in_event_Ev1,2)=zeros(trials_in_event_Ev1,1);
                                    
                                    %Enter Ev2
                                    trials_in_event_Ev2=length(theseEvNosPerEl(evNo2,bwii,which_electrodes(1)).this_wave_logPEv(theseEvNosPerEl(evNo2,bwii,which_electrodes(1)).groupNo==grNo));
                                    total_trials=trials_in_event_Ev1+trials_in_event_Ev2;
                                    this_wave_logPEv2=zeros(trials_in_event_Ev2,1);
                                    
                                    no_elects=0;
                                    for elec=which_electrodes
                                        if length(this_wave_logPEv2)==length(theseEvNosPerEl(evNo2,bwii,elec).groupNo)
                                            this_wave_logPEv2=this_wave_logPEv2+(theseEvNosPerEl(evNo2,bwii,elec).this_wave_logPEv((theseEvNosPerEl(evNo2,bwii,elec).groupNo==grNo))');
                                            no_elects=no_elects+1;
                                        end
                                    end
                                    this_wave_logPEv2=this_wave_logPEv2/no_elects;
                                    
                                    roc_data(trials_in_event_Ev1+1:total_trials,1)=this_wave_logPEv2;
                                    roc_data(trials_in_event_Ev1+1:total_trials,2)=ones(trials_in_event_Ev2,1);
                                    
                                    
                                    %Find  per electrode ROC
                                    if (trials_in_event_Ev1>=5)&(trials_in_event_Ev2>=5)
                                        per_mouse_no_ROCs=per_mouse_no_ROCs+1;
                                        try
                                            roc=roc_calc(roc_data,0,0.05,0);
                                        catch
                                            pffft=1
                                        end
                                        per_mouse_ROCout(per_mouse_no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
                                        per_mouse_ROCbandwidth(per_mouse_no_ROCs)=bwii;
                                        per_mouse_ROCper_ii(per_mouse_no_ROCs)=per_ii;
                                        per_mouse_ROCgroup(per_mouse_no_ROCs)=grNo;
                                        per_mouse_ROCEvNo1(per_mouse_no_ROCs)=evNo1;
                                        per_mouse_ROCEvNo2(per_mouse_no_ROCs)=evNo2;
                                        if ((abs(evNo1-2)<=1)&(abs(evNo2-5)<=1))||((abs(evNo1-5)<=1)&(abs(evNo2-2)<=1))
                                            per_mouse_ROC_between(per_mouse_no_ROCs)=1;
                                        else
                                            per_mouse_ROC_between(per_mouse_no_ROCs)=0;
                                        end
                                        per_mouse_ROC_neighbor(per_mouse_no_ROCs)=abs(evNo1-evNo2);
                                        per_mouse_auROC(per_mouse_no_ROCs)=roc.AUC-0.5;
                                        per_mouse_p_valROC(per_mouse_no_ROCs)=roc.p;
                                        per_mouse_p_vals_ROC=[p_vals_ROC roc.p];
                                        
                                        %I have this code here to plot the ROC
                                        if roc.AUC-0.5>0.3
                                            show_roc=0;
                                            if show_roc==1
                                                %I have this code here to plot the ROC
                                                roc=roc_calc(roc_data,0,0.05,1);
                                                
                                                %Do the histograms
                                                try
                                                    close(2)
                                                catch
                                                end
                                                figure(2)
                                                
                                                hold on
                                                
                                                max_dB=max([max(this_wave_logPEv1) max(this_wave_logPEv2)]);
                                                min_dB=min([min(this_wave_logPEv1) min(this_wave_logPEv2)]);
                                                
                                                edges=[min_dB-0.1*(max_dB-min_dB):(max_dB-min_dB)/20:max_dB+0.1*(max_dB-min_dB)];
                                                histogram(this_wave_logPEv1,edges,'FaceColor','b','EdgeColor','b')
                                                histogram(this_wave_logPEv2,edges,'FaceColor','r','EdgeColor','r')
                                                xlabel('delta power dB')
                                                title(['Histogram for conentrations ' num2str(concs2(evNo1)) ' and ' num2str(concs2(evNo2))])
                                                pffft=1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                if mouse_has_files==1
                    mouse_included(mouseNo)=1;
                else
                    mouse_included(mouseNo)=0;
                end
            end
            
            
            
        end
        fprintf(1, '\n\n')
        
        pFDRauROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for per electrode auROC  = %d\n\n'],pFDRauROC);
        
        per_mouse_pFDRauROC=drsFDRpval(per_mouse_p_vals_ROC);
        fprintf(1, ['pFDR for per mouse auROC  = %d\n\n'],per_mouse_pFDRauROC);
        
        %Plot the lick frequency
        
        %Plot the average
        try
            close(1)
        catch
        end
        hFig=figure(1);
        
        set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
        
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        hold on
        
        data_wave_lickf=[];
        bwii=1;
        
        bar_lab_loc=[];
        no_ev_labels=0;
        ii_gr_included=0;
        
        ii_for_test=0;
        
        fprintf(1, ['\n\n'])
        
        
        for grNo=1:max(handles_drgb.drgbchoices.group_no)
            
            include_group=0;
            
            for evNo=1:length(eventType)
                
                per_included=0;
                these_dB_per_e=[];
                for per_ii=1:2      %performance bins. blue = naive, red = proficient
                    
                    %bar_offset=21-evNo*3+(2-per_ii);
                    if sum(eventType==3)>0
                        bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(evNo-1);
                    else
                        bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                    end
                    
                    these_offsets(per_ii)=bar_offset;
                    
                    if sum((wave_logP_electrode_per_mouse==1)&(wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))>0
                        
                        include_group=1;
                        
                        per_included=per_included+1;
                        
                        if per_ii==1
                            bar(bar_offset,mean(wave_lick_freq_per_mouse((wave_logP_electrode_per_mouse==1)&(wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))),'r','LineWidth', 3)
                        else
                            bar(bar_offset,mean(wave_lick_freq_per_mouse((wave_logP_electrode_per_mouse==1)&(wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))),'b','LineWidth', 3)
                        end
                        
                        %Individual points; in the future add lines linking the
                        %points?
                        plot((bar_offset)*ones(1,sum((wave_logP_electrode_per_mouse==1)&(wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))),...
                            wave_lick_freq_per_mouse((wave_logP_electrode_per_mouse==1)&(wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo)),'o',...
                            'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                        data_for_lines(per_ii).these_dB_per_e(1:length(wave_lick_freq_per_mouse((wave_logP_electrode_per_mouse==1)&(wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))))=...
                            wave_lick_freq_per_mouse((wave_logP_electrode_per_mouse==1)&(wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo));
                        data_for_lines(per_ii).these_mice(1:length(wave_lick_freq_per_mouse((wave_logP_electrode_per_mouse==1)&(wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))))=...
                            wave_logP_mouseNo_per_mouse((wave_logP_electrode_per_mouse==1)&(wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo));
                        data_for_lines(per_ii).these_bar_offsets=bar_offset;
                        
                        
                        %Average and CI
                        plot(bar_offset,mean(wave_lick_freq_per_mouse((wave_logP_electrode_per_mouse==1)&(wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
                        if sum((wave_logP_electrode_per_mouse==1)&(wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))>2
                            CI = bootci(1000, {@mean, wave_lick_freq_per_mouse((wave_logP_electrode_per_mouse==1)&(wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))},'type','cper');
                            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                        end
                        
                        %Save data for ranksum/ttest
                        ii_for_test=ii_for_test+1;
                        data_wave_lickf(ii_for_test).data=wave_lick_freq_per_mouse((wave_logP_electrode_per_mouse==1)&(wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo));
                        data_wave_lickf(ii_for_test).prof_naive=per_ii;
                        data_wave_lickf(ii_for_test).evNo=evNo;
                        data_wave_lickf(ii_for_test).groupNo=grNo;
                    end
                end
                if per_included==2
                    for mouseNo=1:length(mouse_included)
                        if mouse_included(mouseNo)==1
                            if (sum(data_for_lines(1).these_mice==mouseNo)>0)&(sum(data_for_lines(2).these_mice==mouseNo)>0)
                                try
                                    plot([data_for_lines(1).these_bar_offsets*ones(1,sum(data_for_lines(1).these_mice==mouseNo)); data_for_lines(2).these_bar_offsets*ones(1,sum(data_for_lines(1).these_mice==mouseNo))],...
                                        [data_for_lines(1).these_dB_per_e(data_for_lines(1).these_mice==mouseNo); data_for_lines(2).these_dB_per_e(data_for_lines(2).these_mice==mouseNo)],'-','Color',[0.7 0.7 0.7])
                                catch
                                end
                            end
                        end
                    end
                end
                if include_group==1
                    bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                    no_ev_labels=no_ev_labels+1;
                    if sum(eventType==3)>0
                        bar_labels{no_ev_labels}=evTypeLabels{evNo};
                    else
                        bar_labels{no_ev_labels}=num2str(concs(evNo));
                    end
                end
            end
            if include_group==1
                ii_gr_included=ii_gr_included+1;
                groups_included(ii_gr_included)=grNo;
            end
        end
        
        title(['Average lick freqency per mouse'])
        
        %Annotations identifying groups
        x_interval=0.8/ii_gr_included;
        for ii=1:ii_gr_included
            annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.8 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
        end
        
        %Proficient/Naive annotations
        annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
        annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
        
        %x labels
        to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
        sorted_A=sortrows(to_sort);
        sorted_bar_lab_loc=sorted_A(:,1);
        for ii=1:length(bar_lab_loc)
            sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
        end
        xticks(sorted_bar_lab_loc)
        xticklabels(sorted_bar_labels)
        
        if sum(eventType==3)==0
            xlabel('Concentration (%)')
        end
        
        ylabel('Lick frequency (Hz)')
        
        %Now do the ranksums/t tests
        prof_naive_leg{1}='Proficient';
        prof_naive_leg{2}='Naive';
        for ii=1:ii_for_test
            data_wave_lickf(ii).description=[handles_drgb.drgbchoices.group_no_names{data_wave_lickf(ii).groupNo} ' ' evTypeLabels{data_wave_lickf(ii).evNo} ' ' prof_naive_leg{data_wave_lickf(ii).prof_naive} ];
        end
        fprintf(1, ['Ranksum or t-test p values for lick rate\n'])
        output_data = drgMutiRanksumorTtest(data_wave_lickf);
        
        
        
        %Now plot the per mouse delta peak wavelet LFP power computed for each electrode
        pvals=[];
        figOffset=1;
        for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            
            %Plot the average
            try
                close(bwii+figOffset)
            catch
            end
            hFig=figure(bwii+figOffset);
            
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            data_wave_logP=[];
            prof_naive=[];
            events=[];
            mice=[];
            electrodes=[];
            groups=[];
            
            bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            
            fprintf(1, ['\n\n'])
            fprintf(1, ['ANOVAN for peak ERWA power per mouse per electrode, mouse as random factor\n\n'])
            
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                
                for evNo=1:length(eventType)
                    
                    per_included=0;
                    these_dB_per_e=[];
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        %bar_offset=21-evNo*3+(2-per_ii);
                        if sum(eventType==3)>0
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(evNo-1);
                        else
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        end
                        
                        these_offsets(per_ii)=bar_offset;
                        
                        if sum((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))>0
                            
                            include_group=1;
                            
                            per_included=per_included+1;
                            
                            if per_ii==1
                                bar(bar_offset,mean(wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))),'r','LineWidth', 3)
                            else
                                bar(bar_offset,mean(wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))),'b','LineWidth', 3)
                            end
                            
                            %Individual points; in the future add lines linking the
                            %points?
                            plot((bar_offset)*ones(1,sum((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))),...
                                wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo)),'o',...
                                'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            data_for_lines(per_ii).these_dB_per_e(1:length(wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))))=...
                                wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo));
                            data_for_lines(per_ii).these_mice(1:length(wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))))=...
                                wave_logP_mouseNo_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo));
                            data_for_lines(per_ii).these_bar_offsets=bar_offset;
                            
                            
                            %Average and CI
                            plot(bar_offset,mean(wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
                            if sum((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))>2
                                CI = bootci(1000, {@mean, wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))},'type','cper');
                                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            end
                            
                            %Save data for anovan
                            data_wave_logP=[data_wave_logP wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))];
                            prof_naive=[prof_naive per_ii*ones(1,sum((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo)))];
                            events=[events evNo*ones(1,sum((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo)))];
                            mice=[mice wave_logP_mouseNo_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))];
                            electrodes=[electrodes wave_logP_electrode_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))];
                            groups=[groups wave_logP_group_no_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))];
                        end
                    end
                    if per_included==2
                        for mouseNo=1:length(mouse_included)
                            if mouse_included(mouseNo)==1
                                if (sum(data_for_lines(1).these_mice==mouseNo)>0)&(sum(data_for_lines(2).these_mice==mouseNo)>0)
                                    try
                                        plot([data_for_lines(1).these_bar_offsets*ones(1,sum(data_for_lines(1).these_mice==mouseNo)); data_for_lines(2).these_bar_offsets*ones(1,sum(data_for_lines(1).these_mice==mouseNo))],...
                                            [data_for_lines(1).these_dB_per_e(data_for_lines(1).these_mice==mouseNo); data_for_lines(2).these_dB_per_e(data_for_lines(2).these_mice==mouseNo)],'-','Color',[0.7 0.7 0.7])
                                    catch
                                    end
                                end
                            end
                        end
                    end
                    if include_group==1
                        bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                        no_ev_labels=no_ev_labels+1;
                        if sum(eventType==3)>0
                            bar_labels{no_ev_labels}=evTypeLabels{evNo};
                        else
                            bar_labels{no_ev_labels}=num2str(concs(evNo));
                        end
                    end
                end
                if include_group==1
                    ii_gr_included=ii_gr_included+1;
                    groups_included(ii_gr_included)=grNo;
                end
            end
            
            title([freq_names{bwii} ' average peak ERWA delta dB power per mouse, per electrode'])
            
            %Annotations identifying groups
            x_interval=0.8/ii_gr_included;
            for ii=1:ii_gr_included
                annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.8 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %x labels
            to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
            sorted_A=sortrows(to_sort);
            sorted_bar_lab_loc=sorted_A(:,1);
            for ii=1:length(bar_lab_loc)
                sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
            end
            xticks(sorted_bar_lab_loc)
            xticklabels(sorted_bar_labels)
            
            if sum(eventType==3)==0
                xlabel('Concentration (%)')
            end
            
            ylabel('Peak ERWA power (dB)')
            
            
            
            %Calculate anovan for inteaction
            [p,tbl,stats]=anovan(data_wave_logP,{prof_naive events mice electrodes},'varnames',{'proficient_vs_naive','events','groups','mice'},'display','off','random',4);
            fprintf(1, ['p value for anovan peak ERWA dB power per mouse per electrode for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan peak ERWA dB power  per mouse per electrode for events ' freq_names{bwii} '= %d \n'],  p(2));
            fprintf(1, ['p value for anovan peak ERWA dB power  per mouse per electrode for groups for ' freq_names{bwii} '= %d \n\n'],  p(3));
            
            
            %Plot the cumulative histos and do ranksum
            %Plot the average
            
            try
                close(bwii+figOffset+4)
            catch
            end
            hFig=figure(bwii+figOffset+4);
            
            set(hFig, 'units','normalized','position',[.1 .5 .4 .4])
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            ii_rank=0;
            dBpower_rank=[];
            maxdB=-200;
            mindB=200;
            for evNo=1:length(eventType)
                subplot(length(eventType),1,evNo)
                hold on
                
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        
                        
                        if sum((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))>0
                            
                            
                            if per_ii==1
                                if grNo==1
                                    [f_aic,x_aic] = drg_ecdf(wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo)));
                                    plot(x_aic,f_aic,'Color',[1 0 0],'LineWidth',3)
                                else
                                    [f_aic,x_aic] = drg_ecdf(wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo)));
                                    plot(x_aic,f_aic,'Color',[0 0 1],'LineWidth',3)
                                end
                            else
                                if grNo==1
                                    [f_aic,x_aic] = drg_ecdf(wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo)));
                                    plot(x_aic,f_aic,'Color',[1 0 0])
                                else
                                    [f_aic,x_aic] = drg_ecdf(wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo)));
                                    plot(x_aic,f_aic,'Color',[0 0 1])
                                end
                            end
                            
                            
                            
                            %Save data for ranksum
                            ii_rank=ii_rank+1;
                            dBpower_rank(ii_rank).wave_logPpower=wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo));
                            dBpower_rank(ii_rank).per_ii=per_ii;
                            dBpower_rank(ii_rank).grNo=grNo;
                            dBpower_rank(ii_rank).evNo=evNo;
                            maxdB=max([maxdB max(wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo)))]);
                            mindB=min([mindB min(wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo)))]);
                        end
                    end
                    
                end
                legend([handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
                    [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
                title([freq_names{bwii} ' peak ERWA dB power per mouse, per electrode for ' evTypeLabels{evNo}])
                xlabel('peak ERWA power (dB)')
                ylabel('Probability')
            end
            
            for evNo=1:length(eventType)
                subplot(length(eventType),1,evNo)
                xlim([mindB-0.1*(maxdB-mindB) maxdB+0.1*(maxdB-mindB)])
            end
            
            %Now do the ranksums
            fprintf(1, ['Ranksum or t-test p values for peak ERWA dB power per electrode for ' freq_names{bwii} '\n'])
            prof_naive_leg{1}='Proficient';
            prof_naive_leg{2}='Naive';
            for ii=1:ii_rank
                for jj=ii+1:ii_rank
                    [p, r_or_t]=drg_ranksum_or_ttest(dBpower_rank(ii).wave_logPpower,dBpower_rank(jj).wave_logPpower);
                    if r_or_t==0
                        fprintf(1, ['p value ranksum for ' handles_drgb.drgbchoices.group_no_names{dBpower_rank(ii).grNo} ' ' evTypeLabels{dBpower_rank(ii).evNo} ' ' prof_naive_leg{dBpower_rank(ii).per_ii} ' vs ' ...
                            handles_drgb.drgbchoices.group_no_names{dBpower_rank(jj).grNo} ' ' evTypeLabels{dBpower_rank(jj).evNo} ' ' prof_naive_leg{dBpower_rank(jj).per_ii} ' =  %d\n'],p)
                    else
                        fprintf(1, ['p value t-test for ' handles_drgb.drgbchoices.group_no_names{dBpower_rank(ii).grNo} ' ' evTypeLabels{dBpower_rank(ii).evNo} ' ' prof_naive_leg{dBpower_rank(ii).per_ii} ' vs ' ...
                            handles_drgb.drgbchoices.group_no_names{dBpower_rank(jj).grNo} ' ' evTypeLabels{dBpower_rank(jj).evNo} ' ' prof_naive_leg{dBpower_rank(jj).per_ii} ' =  %d\n'],p)
                    end
                    pvals=[pvals p];
                end
            end
            fprintf(1, ['\n\n'])
            
        end
        pFDR = drsFDRpval(pvals);
        fprintf(1, ['pFDR = %d \n\n'],pFDR)
        fprintf(1, ['\n\n'])
        
        
        pvals=[];
        
        figOffset=figOffset+8;
        maxlP=-200;
        minlP=200;
        
        for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                for evNo=1:length(eventType)
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        if sum((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))>0
                            CI=[];
                            CI = bootci(1000, {@mean, wave_logP_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo),:)})';
                            maxlP=max([maxlP max(CI(:))]);
                            minlP=min([minlP min(CI(:))]);
                        end
                    end
                end
            end
        end
        
        for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            
            %Plot the average
            try
                close(bwii+figOffset)
            catch
            end
            hFig=figure(bwii+figOffset);
            
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            
            ii_gr_included=0;
            
            fprintf(1, ['\n\n'])
            fprintf(1, ['ANOVAN for delta dB power per mouse per electrode, mouse as random factor\n\n'])
            
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                subplot(1,2,grNo)
                hold on
                for evNo=1:length(eventType)
                    
                    per_included=0;
                    
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        if sum((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo))>0
                            
                            include_group=1;
                            
                            per_included=per_included+1;
                            
                            mean_wave_logP=[];
                            mean_wave_logP=mean(wave_logP_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo),:),1);
                            
                            CI=[];
                            CI = bootci(1000, {@mean, wave_logP_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo),:)})';
                            maxlP=max([maxlP max(CI(:))]);
                            minlP=min([minlP min(CI(:))]);
                            CI(:,1)= mean_wave_logP'-CI(:,1);
                            CI(:,2)=CI(:,2)- mean_wave_logP';
                            
                            if evNo==1
                                if per_ii==1
                                    [hlCR, hpCR] = boundedline(out_times',mean_wave_logP', CI, 'r');
                                else
                                    [hlCR, hpCR] = boundedline(out_times',mean_wave_logP', CI, 'm');
                                end
                            else
                                if per_ii==1
                                    [hlCR, hpCR] = boundedline(out_times',mean_wave_logP', CI, 'b');
                                else
                                    [hlCR, hpCR] = boundedline(out_times',mean_wave_logP', CI, 'c');
                                end
                            end
                            this_wave_logP=[];
                            this_wave_logP=wave_logP_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_group_no_per_mouse==grNo),:);
                            
                            
                        end
                    end
                    
                    
                end
                
                xlim([-0.2 0.2])
                xlabel('lag (sec)')
                ylabel('ERWA logP (dB)')
                ylim([minlP-0.1*(maxlP-minlP) maxlP+0.1*(maxlP-minlP)])
                title([handles_drgb.drgbchoices.group_no_names{grNo}])
                if grNo==1
                    annotation('textbox',[0.15 0.75 0.3 0.1],'String','Proficient S+','FitBoxToText','on','Color','r','LineStyle','none');
                    annotation('textbox',[0.15 0.72 0.3 0.1],'String','Naive S+','FitBoxToText','on','Color','m','LineStyle','none');
                    annotation('textbox',[0.15 0.69 0.3 0.1],'String','Proficient S-','FitBoxToText','on','Color','b','LineStyle','none');
                    annotation('textbox',[0.15 0.66 0.3 0.1],'String','Naive S-','FitBoxToText','on','Color','c','LineStyle','none');
                end
            end
            
            suptitle([freq_names{bwii} ' average peak ERWA delta dB power per mouse, per electrode'])
            
            
        end
        
        
        
        %Now plot the histograms and the average per mouse LFP power
        %computed per mouse
        figOffset=figOffset+4;
        for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            %Plot the average
            
            try
                close(bwii+figOffset)
            catch
            end
            hFig=figure(bwii+figOffset);
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            
            
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            per_mouse_data_wave_logP=[];
            per_mouse_prof_naive=[];
            per_mouse_events=[];
            per_mouse_groups=[];
            
            bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            
            fprintf(1, ['\n\n'])
            fprintf(1, ['ANOVAN for peak ERWA dB power per mouse, electrode average\n\n'])
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                
                for evNo=1:length(eventType)
                    
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        %                         bar_offset=21-evNo*3+(2-per_ii);
                        
                        if sum(eventType==3)>0
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(evNo-1);
                        else
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        end
                        
                        these_offsets(per_ii)=bar_offset;
                        
                        %Compute per mouse avearge for this group
                        no_mice_for_this_group=0;
                        each_mouse_average_wave_logP=[];
                        for mouseNo=1:max(wave_logP_mouseNo_per_mouse)
                            if sum((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_mouseNo_per_mouse==mouseNo)&(wave_logP_group_no_per_mouse==grNo))>0
                                no_mice_for_this_group=no_mice_for_this_group+1;
                                each_mouse_average_wave_logP(no_mice_for_this_group)=mean(wave_logPpeak_per_mouse((wave_logP_perii_per_mouse==per_ii)&(wave_logP_evNo_per_mouse==evNo)&(wave_logP_bwii_per_mouse==bwii)&(wave_logP_mouseNo_per_mouse==mouseNo)&(wave_logP_group_no_per_mouse==grNo)));
                            end
                        end
                        
                        if no_mice_for_this_group>0
                            
                            include_group=1;
                            
                            if per_ii==1
                                bar(bar_offset,mean(each_mouse_average_wave_logP),'r','LineWidth', 3)
                            else
                                bar(bar_offset,mean(each_mouse_average_wave_logP),'b','LineWidth', 3)
                            end
                            
                            
                            %In the future add lines linking the points
                            plot((bar_offset)*ones(1,no_mice_for_this_group),each_mouse_average_wave_logP,'o',...
                                'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            
                            %Average and CI
                            plot(bar_offset,mean(each_mouse_average_wave_logP),'ok','LineWidth', 3)
                            if no_mice_for_this_group>2
                                CI = bootci(1000, {@mean, each_mouse_average_wave_logP},'type','cper');
                                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            end
                            
                            
                            %Save data for anovan
                            per_mouse_data_wave_logP=[per_mouse_data_wave_logP each_mouse_average_wave_logP];
                            per_mouse_prof_naive=[per_mouse_prof_naive per_ii*ones(1,no_mice_for_this_group)];
                            per_mouse_events=[per_mouse_events evNo*ones(1,no_mice_for_this_group)];
                            per_mouse_groups=[per_mouse_groups grNo*ones(1,no_mice_for_this_group)];
                        end
                    end
                end
                if include_group==1
                    ii_gr_included=ii_gr_included+1;
                    groups_included(ii_gr_included)=grNo;
                end
            end
            if sum(eventType==3)>0
                title([freq_names{bwii} ' auROC peak ERWA dB power per mouse, electrode average'])
            else
                title([freq_names{bwii} ' auROC peak ERWA dB power per mouse, electrode avearage concentrations two steps appart'])
            end
            
            
            
            %Annotations identifying groups
            x_interval=0.8/ii_gr_included;
            for ii=1:ii_gr_included
                annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.8 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %x labels
            to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
            sorted_A=sortrows(to_sort);
            sorted_bar_lab_loc=sorted_A(:,1);
            for ii=1:length(bar_lab_loc)
                sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
            end
            xticks(sorted_bar_lab_loc)
            xticklabels(sorted_bar_labels)
            
            if sum(eventType==3)==0
                xlabel('Concentration (%)')
            end
            
            
            ylabel('peak ERWA power (dB)')
            
            
            
            %Calculate anovan for inteaction
            
            
            [p,tbl,stats]=anovan(per_mouse_data_wave_logP,{per_mouse_prof_naive per_mouse_events per_mouse_groups},'varnames',{'proficient_vs_naive','events','groups'},'display','off');
            fprintf(1, ['p value for anovan peak ERWA dB power histogram per mouse, electrode avearage for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan peak ERWA dB power histogram  per mouse, electrode avearage for events ' freq_names{bwii} '= %d \n'],  p(2));
            fprintf(1, ['p value for anovan peak ERWA dB power histogram  per mouse, electrode avearage for groups ' freq_names{bwii} '= %d \n\n'],  p(2));
            
            
            
        end
        fprintf(1, ['\n\n'])
        
        
        
        
        %Display the auROC for all trials per mouse (per electrode) for for concentrations separated by two log steps
        figOffset=figOffset+4;
        for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            %Plot the average
            
            try
                close(bwii+figOffset)
            catch
            end
            hFig=figure(bwii+figOffset);
            
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            data_auROC=[];
            prof_naive=[];
            within_between=[];
            mice=[];
            groups=[];
            
            bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            
            fprintf(1, ['\n\n'])
            fprintf(1, ['ANOVAN for auROC peak ERWA dB power per mouse per electrode\n\n'])
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                
                for between=0:1
                    
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        %bar_offset=21-evNo*3+(2-per_ii);
                        if sum(eventType==3)>0
                            bar_offset=(grNo-1)*(3.5*2)+(2-(per_ii-1))+3*(between-1);
                        else
                            bar_offset=(grNo-1)*(3.5*2)+(2-(per_ii-1))+3*(2-between);
                        end
                        
                        these_offsets(per_ii)=bar_offset;
                        
                        if sum((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))>0
                            
                            include_group=1;
                            
                            
                            if per_ii==1
                                bar(bar_offset,mean(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))),'r','LineWidth', 3)
                            else
                                bar(bar_offset,mean(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))),'b','LineWidth', 3)
                            end
                            
                            %Individual points; in the future add lines linking the
                            %points?
                            plot((bar_offset)*ones(1,sum((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))),...
                                auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2)),'o',...
                                'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            
                            %Average and CI
                            plot(bar_offset,mean(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))),'ok','LineWidth', 3)
                            CI = bootci(1000, {@mean, auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))},'type','cper');
                            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            
                            %Save data for anovan
                            data_auROC=[data_auROC auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))];
                            prof_naive=[prof_naive per_ii*ones(1,sum((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2)))];
                            within_between=[within_between between*ones(1,sum((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2)))];
                            mice=[mice ROCmouse((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2))];
                            groups=[groups grNo*ones(1,sum((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROC_between==between)&(ROCgroups==grNo)&(ROC_neighbor==2)))];
                        end
                    end
                    if include_group==1
                        bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                        no_ev_labels=no_ev_labels+1;
                        if sum(eventType==3)==0
                            if between==0
                                bar_labels{no_ev_labels}='within';
                            else
                                bar_labels{no_ev_labels}='between';
                            end
                            
                        end
                    end
                end
                if include_group==1
                    ii_gr_included=ii_gr_included+1;
                    groups_included(ii_gr_included)=grNo;
                end
            end
            
            title([freq_names{bwii} ' auROC peak ERWA dB power per mouse, per electrode'])
            
            %Annotations identifying groups
            
            if sum(eventType==3)==0
                x_interval=0.8/ii_gr_included;
                for ii=1:ii_gr_included
                    annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.8 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
                end
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %x labels
            if sum(eventType==3)==0
                to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
                sorted_A=sortrows(to_sort);
                sorted_bar_lab_loc=sorted_A(:,1);
                for ii=1:length(bar_lab_loc)
                    sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
                end
                xticks(sorted_bar_lab_loc)
                xticklabels(sorted_bar_labels)
            else
                for ii=1:ii_gr_included
                    bar_labels{ii}=handles_drgb.drgbchoices.group_no_names{ groups_included(ii)};
                end
                xticks(bar_lab_loc)
                xticklabels(bar_labels)
            end
            
            
            
            ylabel('auROC')
            
            if sum(eventType==3)>0
                %S+ S-
                %Calculate anovan for inteaction
                fprintf(1, ['ANOVAN for auROC peak ERWA dB power per mouse per electrode for concentrations separated by two steps, mouse random factor\n'])
                [p,tbl,stats]=anovan(data_auROC,{prof_naive groups mice},'varnames',{'proficient_vs_naive','groups','mice'},'display','off','random',3);
                fprintf(1, ['p value for anovan auROC peak ERWA dB power per mouse per electrode for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
                fprintf(1, ['p value for anovan auROC peak ERWA dB power per mouse per electrode for groups for ' freq_names{bwii} '= %d \n\n'],  p(2));
                
                %Calculate anovan for inteaction
                fprintf(1, ['ANOVAN for auROC peak ERWA dB power per mouse per electrode for concentrations separated by two steps\n'])
                [p,tbl,stats]=anovan(data_auROC,{prof_naive groups},'varnames',{'proficient_vs_naive','groups'},'display','off');
                fprintf(1, ['p value for anovan auROC peak ERWA dB power per mouse per electrode for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
                fprintf(1, ['p value for anovan auROC peak ERWA dB power per mouse per electrode for groups for ' freq_names{bwii} '= %d \n\n'],  p(2));
            else
                %Calculate anovan for inteaction
                fprintf(1, ['ANOVAN for auROC per mouse per electrode for concentrations separated by two steps, mouse random factor\n'])
                [p,tbl,stats]=anovan(data_auROC,{prof_naive within_between groups mice},'varnames',{'proficient_vs_naive','within_between','groups','mice'},'display','off','random',4);
                fprintf(1, ['p value for anovan auROC peak ERWA dB power per mouse per electrode for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
                fprintf(1, ['p value for anovan auROC peak ERWA dB power per mouse per electrode for within vs between ' freq_names{bwii} '= %d \n'],  p(2));
                fprintf(1, ['p value for anovan auROC peak ERWA dB power per mouse per electrode for groups for ' freq_names{bwii} '= %d \n\n'],  p(3));
                
                
                %Calculate anovan for inteaction
                fprintf(1, ['ANOVAN for auROC per mouse per electrode for concentrations separated by two steps\n'])
                [p,tbl,stats]=anovan(data_auROC,{prof_naive within_between groups},'varnames',{'proficient_vs_naive','within_between','groups'},'display','off');
                fprintf(1, ['p value for anovan auROC peak ERWA dB power per mouse per electrode for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
                fprintf(1, ['p value for anovan auROC peak ERWA dB power per mouse per electrode for within vs between ' freq_names{bwii} '= %d \n'],  p(2));
                fprintf(1, ['p value for anovan auROC peak ERWA dB power per mouse per electrode for groups for ' freq_names{bwii} '= %d \n\n'],  p(3));
                
            end
        end
        
        fprintf(1, ['\n\n'])
        
        
        
        %Display cumulative histograms for the auROC for all trials per mouse (per electrode) for for concentrations separated by two log steps
        %This only works for Daniel's NRG, we have to modify in the future
        pvals=[];
        if sum(eventType==3)>0
            figOffset=figOffset+4;
            for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
                %Plot the average
                
                try
                    close(bwii+figOffset)
                catch
                end
                hFig=figure(bwii+figOffset);
                
                set(hFig, 'units','normalized','position',[.2 .2 .6 .6])
                
                set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                hold on
                
                ii_rank=0;
                
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    include_group=0;
                    
                    
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        if sum((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo))>0
                            
                            include_group=1;
                            
                            
                            if per_ii==1
                                if grNo==1
                                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo)));
                                    plot(x_aic,f_aic,'Color',[1 0 0],'LineWidth',3)
                                else
                                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo)));
                                    plot(x_aic,f_aic,'Color',[0 0 1],'LineWidth',3)
                                end
                            else
                                if grNo==1
                                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo)));
                                    plot(x_aic,f_aic,'Color',[1 0 0])
                                else
                                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo)));
                                    plot(x_aic,f_aic,'Color',[0 0 1])
                                end
                            end
                            
                            
                            %Save data for ranksum
                            ii_rank=ii_rank+1;
                            ranksum_auROC(ii_rank).auROC=auROC((ROCper_ii==per_ii)&(ROCbandwidth==bwii)&(ROCgroups==grNo));
                            ranksum_auROC(ii_rank).per_ii=per_ii;
                            ranksum_auROC(ii_rank).grNo=grNo;
                        end
                    end
                    
                    if include_group==1
                        ii_gr_included=ii_gr_included+1;
                        groups_included(ii_gr_included)=grNo;
                    end
                end
                
                title([freq_names{bwii} ' auROC peak ERWA dB power per mouse, per electrode'])
                
                legend([handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive']...
                    ,[handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
                
                xlabel('auROC')
                ylabel('Cumulative probability')
                
                %Now do the ranksums
                fprintf(1, ['Ranksum or t-test p values for ' freq_names{bwii} '\n'])
                prof_naive_leg{1}='Proficient';
                prof_naive_leg{2}='Naive';
                for ii=1:ii_rank
                    for jj=ii+1:ii_rank
                        [p, r_or_t]=drg_ranksum_or_ttest(ranksum_auROC(ii).auROC,ranksum_auROC(jj).auROC);
                        if r_or_t==0
                            fprintf(1, ['p value ranksum for ' handles_drgb.drgbchoices.group_no_names{ranksum_auROC(ii).grNo} ' ' prof_naive_leg{ranksum_auROC(ii).per_ii} ' vs ' ...
                                handles_drgb.drgbchoices.group_no_names{ranksum_auROC(jj).grNo} ' ' prof_naive_leg{ranksum_auROC(jj).per_ii} ' =  %d\n'],p)
                        else
                            fprintf(1, ['p value t-test for ' handles_drgb.drgbchoices.group_no_names{ranksum_auROC(ii).grNo} ' ' prof_naive_leg{ranksum_auROC(ii).per_ii} ' vs ' ...
                                handles_drgb.drgbchoices.group_no_names{ranksum_auROC(jj).grNo} ' ' prof_naive_leg{ranksum_auROC(jj).per_ii} ' =  %d\n'],p)
                        end
                        pvals=[pvals p];
                    end
                end
                fprintf(1, ['\n\n'])
                
            end
        end
        fprintf(1, ['\n\n'])
        pFDR = drsFDRpval(pvals);
        fprintf(1, ['pFDR = %d \n\n'],pFDR)
        
        
        pffft=1;
        
    case 23
        %23 Oscillatory power at the peak and trough of the PAC
        
        mean_PACpower_No_per_mouse=0;
        
        mean_PACpower_No=0;
        mean_peakPACpower=[];
        mean_troughPACpower=[];
        mean_PACpower_perii=[];
        mean_PACpower_evNo=[];
        mean_PACpower_pacii=[];
        mean_PACpower_fileNo=[];
        per_session_group_no=[];
        mean_VL=[];
        mean_VA=[];
        mean_PA=[];
        
        
        fprintf(1, ['PAC power analysis using Hilbert transform was not used for Justin''s paper\n\n'])
        p_vals=[];
        no_files=max(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        
        szpc=size(percent_windows);
        for per_ii=1:szpc(1)
            
            
            
            for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
                mouse_has_files=0;
                for elec=1:16
                    if sum(which_electrodes==elec)>0
                        
                        theseEvNos_thisMouse_thisElec=[];
                        for evNo=1:length(eventType)
                            for pacii=1:no_pacii
                                for group_no=1:max(handles_drgb.drgbchoices.group_no)
                                    theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv=0;
                                end
                            end
                        end
                        
                        
                        for fileNo=1:no_files
                            %If this file is in the list of files the user wants to process in drgAnalysisBatchLFP continue
                            if sum(files==fileNo)>0
                                if handles_drgb.drgbchoices.mouse_no(fileNo)==mouseNo
                                    lfpodNo=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                    lfpodRefNo=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                                    
                                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo)))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodRefNo)))
                                        
                                        
                                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo).PAC))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodRefNo).PAC))
                                            
                                            percent_mask=[];
                                            percent_mask=logical((handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).perCorrPAC>=percent_windows(per_ii,1))...
                                                &(handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).perCorrPAC<=percent_windows(per_ii,2)));
                                            
                                            if ~isempty(percent_mask)
                                                
                                                
                                                for evNo=1:length(eventType)
                                                    
                                                    trials_in_event_Ev=[];
                                                    trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).which_eventPAC(eventType(evNo),:)==1)&percent_mask;
                                                    
                                                    if (sum(trials_in_event_Ev)>=1)
                                                        
                                                        %Do per bandwidth analysis
                                                        for pacii=1:no_pacii
                                                            
                                                            group_no=handles_drgb.drgbchoices.group_no(fileNo);
                                                            
                                                            %Enter the PACpower
                                                            this_peakPACpower_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_peakPACpower_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).meanPeakPower(trials_in_event_Ev)-handles_drgb.drgb.lfpevpair(lfpodRefNo).PAC(pacii).meanPeakPower(trials_in_event_Ev);
                                                            this_troughPACpower_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_troughPACpower_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).meanTroughPower(trials_in_event_Ev)-handles_drgb.drgb.lfpevpair(lfpodRefNo).PAC(pacii).meanTroughPower(trials_in_event_Ev);
                                                            
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_peakPACpower_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1 ...
                                                                :theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=this_peakPACpower_Ev;
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_troughPACpower_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1 ...
                                                                :theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=this_troughPACpower_Ev;
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).whichMouse(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=mouseNo*ones(1,length(this_peakPACpower_Ev));
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev);
                                                            
                                                            %Save per session value for peak power
                                                            mean_PACpower_No=mean_PACpower_No+1;
                                                            mean_peakPACpower(mean_PACpower_No)=mean(this_peakPACpower_Ev);
                                                            mean_troughPACpower(mean_PACpower_No)=mean(this_troughPACpower_Ev);
                                                            mean_PACpower_perii(mean_PACpower_No)=per_ii;
                                                            mean_PACpower_evNo(mean_PACpower_No)=evNo;
                                                            mean_PACpower_pacii(mean_PACpower_No)=pacii;
                                                            mean_PACpower_fileNo(mean_PACpower_No)=fileNo;
                                                            mean_PACpower_mouse(mean_PACpower_No)=mouseNo;
                                                            per_session_group_no(mean_PACpower_No)=handles_drgb.drgbchoices.group_no(fileNo);
                                                            
                                                            
                                                            mouse_has_files=1;
                                                            
                                                        end
                                                        
                                                        
                                                        fprintf(1, ['%d trials in event No %d succesfully processed for file No %d electrode %d\n'],sum(trials_in_event_Ev), evNo,fileNo,elec);
                                                        
                                                    else
                                                        
                                                        
                                                        fprintf(1, ['%d trials in event No %d fewer than minimum trials per event ' evTypeLabels{evNo} ' for file No %d electrode %d\n'],sum(trials_in_event_Ev), min_trials_per_event,fileNo,elec);
                                                        
                                                        
                                                    end
                                                    
                                                    
                                                end
                                                
                                            else
                                                
                                                fprintf(1, ['Empty percent_mask for file No %d electrode %d\n'],fileNo,elec);
                                                
                                            end
                                            
                                        else
                                            
                                            fprintf(1, ['Empty PAC for file No %d electrode %d\n'],fileNo,elec);
                                            
                                        end
                                        
                                        
                                    else
                                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],fileNo,elec);
                                        
                                        
                                    end
                                end %if mouseNo
                            end
                        end %fileNo
                        
                        if theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv>0
                            
                            %Calculate per mouse PAC power
                            for evNo=1:length(eventType)
                                for pacii=1:no_pacii
                                    for group_no=1:max(handles_drgb.drgbchoices.group_no)
                                        %Calculate per mouse PAC power
                                        mean_PACpower_No_per_mouse=mean_PACpower_No_per_mouse+1;
                                        this_mouse_peakPACpower=[];
                                        this_mouse_peakPACpower=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_peakPACpower_Ev;
                                        this_mouse_troughPACpower=[];
                                        this_mouse_troughPACpower=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_troughPACpower_Ev;
                                        if ~isempty(this_mouse_peakPACpower)
                                            mean_peakPACpower_per_mouse(mean_PACpower_No_per_mouse)=mean(this_mouse_peakPACpower);
                                            mean_troughPACpower_per_mouse(mean_PACpower_No_per_mouse)=mean(this_mouse_troughPACpower);
                                        else
                                            mean_peakPACpower_per_mouse(mean_PACpower_No_per_mouse)=NaN;
                                            mean_troughPACpower_per_mouse(mean_PACpower_No_per_mouse)=NaN;
                                        end
                                        
                                        mean_PACpower_perii_per_mouse(mean_PACpower_No_per_mouse)=per_ii;
                                        mean_PACpower_evNo_per_mouse(mean_PACpower_No_per_mouse)=evNo;
                                        mean_PACpower_pacii_per_mouse(mean_PACpower_No_per_mouse)=pacii;
                                        mean_PACpower_mouseNo_per_mouse(mean_PACpower_No_per_mouse)=mouseNo;
                                        mean_PACpower_electNo_per_mouse(mean_PACpower_No_per_mouse)=elec;
                                        mean_PACpower_group_no_per_mouse(mean_PACpower_No_per_mouse)=group_no;
                                        
                                    end
                                end
                            end
                        end
                        
                    end
                    
                end
                
                if mouse_has_files==1
                    mouse_included(mouseNo)=1;
                else
                    mouse_included(mouseNo)=0;
                end
            end
            
            
        end
        fprintf(1, '\n\n')
        
        
        %Now plot the average peakPACpower for each electrode calculated per mouse
        %(including all sessions for each mouse)
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            data_PACpower=[];
            prof_naive=[];
            events=[];
            groups=[];
            mice=[];
            
            %Plot the average
            try
                close(pacii)
            catch
            end
            hFig=figure(pacii);
            
            set(hFig, 'units','normalized','position',[.1 .2 .7 .7])
            subplot(2,1,1)
            hold on
            
            bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                
                for evNo=1:length(eventType)
                    
                    for per_ii=1:2
                        
                        if sum(eventType==3)>0
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
                        else
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        end
                        
                        these_offsets(per_ii)=bar_offset;
                        
                        if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
                            
                            include_group=1;
                            
                            if per_ii==1
                                bar(bar_offset,mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'r','LineWidth', 3,'EdgeColor','none')
                            else
                                bar(bar_offset,mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                            end
                            
                            
                            plot(bar_offset,mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
                            plot((bar_offset)*ones(1,sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),...
                                mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)),'o',...
                                'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            
                            if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>=2
                                CI = bootci(1000, {@mean, mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))},'type','cper');
                                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            end
                            
                            %Save data for anovan
                            data_PACpower=[data_PACpower mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
                            prof_naive=[prof_naive per_ii*ones(1,sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            events=[events evNo*ones(1,sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            groups=[groups grNo*ones(1,sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            mice=[mice mean_PACpower_mouseNo_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
                            
                        end
                    end
                    if include_group==1
                        bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                        no_ev_labels=no_ev_labels+1;
                        if sum(eventType==3)>0
                            bar_labels{no_ev_labels}=evTypeLabels{evNo};
                        else
                            bar_labels{no_ev_labels}=num2str(concs(evNo));
                        end
                    end
                end
                if include_group==1
                    ii_gr_included=ii_gr_included+1;
                    groups_included(ii_gr_included)=grNo;
                end
                
            end
            
            title('Peak PAC power')
            
            
            %Annotations identifying groups
            x_interval=0.8/ii_gr_included;
            for ii=1:ii_gr_included
                annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %x labels
            to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
            sorted_A=sortrows(to_sort);
            sorted_bar_lab_loc=sorted_A(:,1);
            for ii=1:length(bar_lab_loc)
                sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
            end
            xticks(sorted_bar_lab_loc)
            xticklabels(sorted_bar_labels)
            
            if sum(eventType==3)==0
                xlabel('Concentration (%)')
            end
            
            ylabel('dB')
            ylim1(pacii,:)=ylim;
            
            %Calculate anovan for inteaction
            fprintf(1, ['ANOVAN for average peak PAC power for each electrode calculated per mouse with mouse as random factor PAC theta/' freq_names{pacii+1} '\n'])
            [p,tbl,stats]=anovan(data_PACpower,{prof_naive, events, groups, mice},'varnames',{'proficient_vs_naive','events','groups','mice'},'display','off','random',4);
            fprintf(1, ['p value for anovan peak PAC power per mouse per electrode for naive vs proficient for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan peak PAC power per mouse per electrode for events for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(2));
            fprintf(1, ['p value for anovan peak PAC power per mouse per electrode for groups for PAC theta/' freq_names{pacii+1} '= %d \n\n'],  p(3));
            
            %Calculate anovan for inteaction
            fprintf(1, ['ANOVAN for average peak PAC power for each electrode calculated per mouse without mouse as a factor for PAC theta/' freq_names{pacii+1} '\n'])
            [p,tbl,stats]=anovan(data_PACpower,{prof_naive, events, groups},'varnames',{'proficient_vs_naive','events','groups'},'display','off');
            fprintf(1, ['p value for anovan peak PAC power per mouse per electrode for naive vs proficient for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan peak PAC power per mouse per electrode for events for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(2));
            fprintf(1, ['p value for anovan peak PAC power per mouse per electrode for groups for PAC theta/' freq_names{pacii+1} '= %d \n\n'],  p(3));
            
        end
        
        
        %Now plot the average troughPACpower for each electrode calculated per mouse
        %(including all sessions for each mouse)
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            data_PACpower=[];
            prof_naive=[];
            events=[];
            groups=[];
            mice=[];
            
            %Plot the average
            
            hFig=figure(pacii);
            subplot(2,1,2)
            hold on
            
            bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                
                for evNo=1:length(eventType)
                    
                    for per_ii=1:2
                        
                        if sum(eventType==3)>0
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
                        else
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        end
                        
                        these_offsets(per_ii)=bar_offset;
                        
                        if sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
                            
                            include_group=1;
                            
                            if per_ii==1
                                bar(bar_offset,mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'r','LineWidth', 3,'EdgeColor','none')
                            else
                                bar(bar_offset,mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                            end
                            
                            
                            plot(bar_offset,mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
                            plot((bar_offset)*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),...
                                mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)),'o',...
                                'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            
                            if sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>=2
                                CI = bootci(1000, {@mean, mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))},'type','cper');
                                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            end
                            
                            %Save data for anovan
                            data_PACpower=[data_PACpower mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
                            prof_naive=[prof_naive per_ii*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            events=[events evNo*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            groups=[groups grNo*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            mice=[mice mean_PACpower_mouseNo_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
                            
                        end
                    end
                    if include_group==1
                        bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                        no_ev_labels=no_ev_labels+1;
                        if sum(eventType==3)>0
                            bar_labels{no_ev_labels}=evTypeLabels{evNo};
                        else
                            bar_labels{no_ev_labels}=num2str(concs(evNo));
                        end
                    end
                end
                if include_group==1
                    ii_gr_included=ii_gr_included+1;
                    groups_included(ii_gr_included)=grNo;
                end
            end
            
            title('Trough PAC power')
            
            suptitle(['Average PAC power for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1}])
            
            
            %Annotations identifying groups
            x_interval=0.8/ii_gr_included;
            for ii=1:ii_gr_included
                annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %x labels
            to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
            sorted_A=sortrows(to_sort);
            sorted_bar_lab_loc=sorted_A(:,1);
            for ii=1:length(bar_lab_loc)
                sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
            end
            xticks(sorted_bar_lab_loc)
            xticklabels(sorted_bar_labels)
            
            if sum(eventType==3)==0
                xlabel('Concentration (%)')
            end
            
            ylabel('dB')
            ylim2=ylim;
            ylim([min([ylim1(pacii,1) ylim2(1)]) max([ylim1(pacii,2) ylim2(2)])])
            
        end
        
        %Now do the cumulative histograms and ranksums for PAC power for each electrode calculated with all sessons per mouse
        pvals=[];
        legends=[];
        prof_naive_leg{1}='Proficient';
        prof_naive_leg{2}='Naive';
        
        
        for pacii=1:no_pacii
            
            %Find max and min to use for xlim
            maxpower=-200;
            minpower=200;
            for evNo=1:length(eventType)
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
                            maxpower=max([maxpower max(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))]);
                            minpower=min([minpower min(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))]);
                        end
                    end
                end
            end
            
            try
                close(pacii+3)
            catch
            end
            hFig=figure(pacii+3);
            
            
            set(hFig, 'units','normalized','position',[.1 .1 .7 .8])
            
            
            %Plot the peak PAC power
            ii_rank=0;
            input_data=[];
            for evNo=1:length(eventType)
                
                subplot(2,2,evNo)
                hold on
                
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    
                    
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        
                        if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
                            
                            [f_mi,x_mi] = drg_ecdf(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)));
                            cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_mi=f_mi;
                            cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_mi=x_mi;
                            if grNo==1
                                if per_ii==1
                                    legends.pacii(pacii).evNo(evNo).p1=plot(x_mi,f_mi,'Color',[1 0 0],'LineWidth',3);
                                else
                                    legends.pacii(pacii).evNo(evNo).p2=plot(x_mi,f_mi,'Color',[0 0 1],'LineWidth',3);
                                end
                            else
                                if per_ii==1
                                    legends.pacii(pacii).evNo(evNo).p3=plot(x_mi,f_mi,'Color',[1 0.7 0.7],'LineWidth',3);
                                else
                                    legends.pacii(pacii).evNo(evNo).p4=plot(x_mi,f_mi,'Color',[0.7 0.7 1],'LineWidth',3);
                                end
                            end
                            
                            
                            %Compute per mouse avearge
                            each_mouse_average_PACpower=[];
                            no_mice_included=0;
                            for mouseNo=1:max(mean_PACpower_mouseNo_per_mouse)
                                if mouse_included(mouseNo)==1
                                    if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
                                        if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo))>0
                                            no_mice_included=no_mice_included+1;
                                            each_mouse_average_PACpower(no_mice_included)=mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo)));
                                        end
                                    end
                                end
                            end
                            
                            %Show average per mouse
                            if no_mice_included>0
                                for jj=1:length(each_mouse_average_PACpower)
                                    this_f_mi=[];
                                    this_f_mi=cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_mi;
                                    
                                    this_x_mi=[];
                                    this_x_mi=cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_mi;
                                    
                                    xii_below=find(this_x_mi<each_mouse_average_PACpower(jj),1,'last');
                                    xii_above=find(this_x_mi>each_mouse_average_PACpower(jj),1,'first');
                                    
                                    slope=(this_f_mi(xii_above)-this_f_mi(xii_below))/(this_x_mi(xii_above)-this_x_mi(xii_below));
                                    intercept=this_f_mi(xii_above)-slope*this_x_mi(xii_above);
                                    
                                    this_f=slope*each_mouse_average_PACpower(jj)+intercept;
                                    
                                    if grNo==1
                                        if per_ii==1
                                            plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
                                        else
                                            plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
                                        end
                                    else
                                        if per_ii==1
                                            plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
                                        else
                                            plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
                                        end
                                    end
                                    
                                end
                            end
                            
                            %Save data for ranksum
                            ii_rank=ii_rank+1;
                            input_data(ii_rank).data=mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo));
                            input_data(ii_rank).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                        end
                    end
                    
                end
                
                title(evTypeLabels{evNo})
                xlabel('Peak PAC power (dB)')
                ylabel('Probability')
                xlim([minpower-0.1*(maxpower-minpower) maxpower+0.1*(maxpower-minpower)])
            end
            
            
            %Now do the ranksums
            fprintf(1, ['Ranksum or t-test p values for peak PAC power for each electrode calculated per mouse for PAC theta' freq_names{pacii+1} '\n'])
            output_data = drgMutiRanksumorTtest(input_data);
            
            fprintf(1, ['\n\n'])
            
            
            %Plot the trough PAC power
            ii_rank=0;
            input_data=[];
            for evNo=1:length(eventType)
                
                subplot(2,2,evNo+2)
                hold on
                
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    
                    
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        
                        if sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
                            
                            [f_mi,x_mi] = drg_ecdf(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)));
                            cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_mi=f_mi;
                            cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_mi=x_mi;
                            if grNo==1
                                if per_ii==1
                                    legends.pacii(pacii).evNo(evNo).p1=plot(x_mi,f_mi,'Color',[1 0 0],'LineWidth',3);
                                else
                                    legends.pacii(pacii).evNo(evNo).p2=plot(x_mi,f_mi,'Color',[0 0 1],'LineWidth',3);
                                end
                            else
                                if per_ii==1
                                    legends.pacii(pacii).evNo(evNo).p3=plot(x_mi,f_mi,'Color',[1 0.7 0.7],'LineWidth',3);
                                else
                                    legends.pacii(pacii).evNo(evNo).p4=plot(x_mi,f_mi,'Color',[0.7 0.7 1],'LineWidth',3);
                                end
                            end
                            
                            %Compute per mouse avearge
                            each_mouse_average_PACpower=[];
                            no_mice_included=0;
                            for mouseNo=1:max(mean_PACpower_mouseNo_per_mouse)
                                if mouse_included(mouseNo)==1
                                    if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
                                        if sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo))>0
                                            no_mice_included=no_mice_included+1;
                                            each_mouse_average_PACpower(no_mice_included)=mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo)));
                                        end
                                    end
                                end
                            end
                            
                            %Show average per mouse
                            if no_mice_included>0
                                for jj=1:length(each_mouse_average_PACpower)
                                    this_f_mi=[];
                                    this_f_mi=cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_mi;
                                    
                                    this_x_mi=[];
                                    this_x_mi=cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_mi;
                                    
                                    xii_below=find(this_x_mi<each_mouse_average_PACpower(jj),1,'last');
                                    xii_above=find(this_x_mi>each_mouse_average_PACpower(jj),1,'first');
                                    
                                    slope=(this_f_mi(xii_above)-this_f_mi(xii_below))/(this_x_mi(xii_above)-this_x_mi(xii_below));
                                    intercept=this_f_mi(xii_above)-slope*this_x_mi(xii_above);
                                    
                                    this_f=slope*each_mouse_average_PACpower(jj)+intercept;
                                    
                                    if grNo==1
                                        if per_ii==1
                                            plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
                                        else
                                            plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
                                        end
                                    else
                                        if per_ii==1
                                            plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
                                        else
                                            plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
                                        end
                                    end
                                    
                                end
                            end
                            
                            %Save data for ranksum
                            ii_rank=ii_rank+1;
                            input_data(ii_rank).data=mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo));
                            input_data(ii_rank).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                        end
                    end
                    
                end
                
                title(evTypeLabels{evNo})
                xlabel('Trough PAC power (dB)')
                ylabel('Probability')
                xlim([minpower-0.1*(maxpower-minpower) maxpower+0.1*(maxpower-minpower)])
            end
            
            suptitle(['Average PAC power for each electrode calculated per  for PAC theta/' freq_names{pacii+1}])
            
            %Now do the ranksums
            fprintf(1, ['Ranksum or t-test p values for trough PAC power for each electrode calculated per mouse for PAC theta' freq_names{pacii+1} '\n'])
            output_data = drgMutiRanksumorTtest(input_data);
            
            fprintf(1, ['\n\n'])
            
        end
        
        
        %Now plot the average PAC power per mouse averaged over electrodes
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            
            %Plot the average PAC power per mouse averaged over electrodes
            try
                close(pacii+6)
            catch
            end
            hFig=figure(pacii+6);
            set(hFig, 'units','normalized','position',[.1 .5 .7 .7])
            
            %Peak PAC power
            subplot(2,1,1)
            hold on
            
            
            %             bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            
            ii_rank=0;
            input_data=[];
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                
                for evNo=1:length(eventType)
                    
                    for per_ii=1:2
                        
                        if sum(eventType==3)>0
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
                        else
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        end
                        
                        these_offsets(per_ii)=bar_offset;
                        
                        if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
                            
                            include_group=1;
                            
                            %Compute per mouse avearge
                            each_mouse_average_PACpower=[];
                            no_mice_included=0;
                            for mouseNo=1:max(mean_PACpower_mouseNo_per_mouse)
                                if mouse_included(mouseNo)==1
                                    if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
                                        if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo))>0
                                            no_mice_included=no_mice_included+1;
                                            each_mouse_average_PACpower(no_mice_included)=mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo)));
                                        end
                                    end
                                end
                            end
                            if no_mice_included>0
                                
                                include_group=1;
                                
                                if per_ii==1
                                    bar(bar_offset,mean(each_mouse_average_PACpower),'r','LineWidth', 3,'EdgeColor','none')
                                else
                                    bar(bar_offset,mean(each_mouse_average_PACpower),'b','LineWidth', 3,'EdgeColor','none')
                                end
                                
                                
                                plot(bar_offset,mean(each_mouse_average_PACpower),'ok','LineWidth', 3)
                                plot((bar_offset)*ones(1,length(each_mouse_average_PACpower)),each_mouse_average_PACpower,'o',...
                                    'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                                
                                if length(each_mouse_average_PACpower)>2
                                    CI = bootci(1000, {@mean, each_mouse_average_PACpower},'type','cper');
                                    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                                end
                                
                                
                                
                                %Save data for ranksum
                                ii_rank=ii_rank+1;
                                input_data(ii_rank).data=each_mouse_average_PACpower;
                                input_data(ii_rank).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                                
                            end
                        end
                    end
                end
                
                
                if include_group==1
                    ii_gr_included=ii_gr_included+1;
                    groups_included(ii_gr_included)=grNo;
                end
            end
            
            title('Peak PAC power')
            
            %Annotations identifying groups
            x_interval=0.8/ii_gr_included;
            for ii=1:ii_gr_included
                annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            
            xticks(sorted_bar_lab_loc)
            xticklabels(sorted_bar_labels)
            
            if sum(eventType==3)==0
                xlabel('Concentration (%)')
            end
            
            ylabel('dB')
            
            %Now do the ranksums
            fprintf(1, ['Ranksum or t-test p values for peak PAC power per mouse averaged over electrodes for PAC theta PAC theta' freq_names{pacii+1} '\n'])
            output_data = drgMutiRanksumorTtest(input_data);
            
            %Trough power
            subplot(2,1,2)
            hold on
            
            
            %             bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            
            ii_rank=0;
            input_data=[];
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                
                for evNo=1:length(eventType)
                    
                    for per_ii=1:2
                        
                        if sum(eventType==3)>0
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
                        else
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        end
                        
                        these_offsets(per_ii)=bar_offset;
                        
                        if sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
                            
                            include_group=1;
                            
                            %Compute per mouse avearge
                            each_mouse_average_PACpower=[];
                            no_mice_included=0;
                            for mouseNo=1:max(mean_PACpower_mouseNo_per_mouse)
                                if mouse_included(mouseNo)==1
                                    if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
                                        if sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo))>0
                                            no_mice_included=no_mice_included+1;
                                            each_mouse_average_PACpower(no_mice_included)=mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo)));
                                        end
                                    end
                                end
                            end
                            if no_mice_included>0
                                
                                include_group=1;
                                
                                if per_ii==1
                                    bar(bar_offset,mean(each_mouse_average_PACpower),'r','LineWidth', 3,'EdgeColor','none')
                                else
                                    bar(bar_offset,mean(each_mouse_average_PACpower),'b','LineWidth', 3,'EdgeColor','none')
                                end
                                
                                
                                plot(bar_offset,mean(each_mouse_average_PACpower),'ok','LineWidth', 3)
                                plot((bar_offset)*ones(1,length(each_mouse_average_PACpower)),each_mouse_average_PACpower,'o',...
                                    'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                                
                                if length(each_mouse_average_PACpower)>2
                                    CI = bootci(1000, {@mean, each_mouse_average_PACpower},'type','cper');
                                    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                                end
                                
                                
                                
                                %Save data for ranksum
                                ii_rank=ii_rank+1;
                                input_data(ii_rank).data=each_mouse_average_PACpower;
                                input_data(ii_rank).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
                                
                            end
                        end
                    end
                end
                
                
                if include_group==1
                    ii_gr_included=ii_gr_included+1;
                    groups_included(ii_gr_included)=grNo;
                end
            end
            
            title('Trough PAC power')
            suptitle(['Average PAC power per mouse averaged over all electrodes for PAC theta/' freq_names{pacii+1}])
            
            %Annotations identifying groups
            x_interval=0.8/ii_gr_included;
            for ii=1:ii_gr_included
                annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %             %x labels
            %             to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
            %             sorted_A=sortrows(to_sort);
            %             sorted_bar_lab_loc=sorted_A(:,1);
            %             for ii=1:length(bar_lab_loc)
            %                 sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
            %             end
            xticks(sorted_bar_lab_loc)
            xticklabels(sorted_bar_labels)
            
            if sum(eventType==3)==0
                xlabel('Concentration (%)')
            end
            
            ylabel('dB')
            
            %Now do the ranksums
            fprintf(1, ['Ranksum or t-test p values for trough PAC power per mouse averaged over electrodes for PAC theta PAC theta' freq_names{pacii+1} '\n'])
            output_data = drgMutiRanksumorTtest(input_data);
            
        end
        
        pfft=1;
        
    case 24
        %24 Oscillatory wavelet power at the peak and trough of the PAC
        
        mean_PACpower_No_per_mouse=0;
        mean_deltaLickFreq_No_per_mouse=0;
        
        mean_PACpower_No=0;
        mean_peakPACpower=[];
        mean_troughPACpower=[];
        mean_PACpower_perii=[];
        mean_PACpower_evNo=[];
        mean_PACpower_pacii=[];
        mean_PACpower_fileNo=[];
        per_session_group_no=[];
        mean_VL=[];
        mean_VA=[];
        mean_PA=[];
        
        handles_out=[];
        handles_out.PRP_ii=0;
        handles_out.LickF_ii=0;
        handles_out.AUC_ii=0;
        
        prof_naive_leg{1}='Proficient';
        prof_naive_leg{2}='Naive';
        
        %Initialize ROC
        no_ROCs=0;
        ROCfileNo=[]
        ROCelec=[];
        ROCgroups=[];
        ROCmouse=[];
        ROCpacii=[];
        ROCper_ii=[];
        ROCEvNo1=[];
        ROCEvNo2=[];
        ROC_between=[];
        auROCpeak=[];
        auROCtrough=[];
        p_valROCpeak=[];
        p_vals_ROCpeak=[];
        p_valROCtrough=[];
        p_vals_ROCtrough=[];
        
        fprintf(1, ['PAC power analysis using wavelet power for Justin''s paper\n\n'])
        p_vals=[];
        no_files=max(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        
        szpc=size(percent_windows);
        for per_ii=1:szpc(1)
            
            
            
            for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
                mouse_has_files=0;
                for elec=1:16
                    if sum(which_electrodes==elec)>0
                        
                        theseEvNos_thisMouse_thisElec=[];
                        for evNo=1:length(eventType)
                            for pacii=1:no_pacii
                                for group_no=1:max(handles_drgb.drgbchoices.group_no)
                                    theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv=0;
                                end
                            end
                        end
                        
                        
                        for fileNo=1:no_files
                            %If this file is in the list of files the user wants to process in drgAnalysisBatchLFP continue
                            if sum(files==fileNo)>0
                                if handles_drgb.drgbchoices.mouse_no(fileNo)==mouseNo
                                    lfpodNo=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                    lfpodRefNo=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                                    
                                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo)))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodRefNo)))
                                        
                                        
                                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo).PAC))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodRefNo).PAC))
                                            
                                            percent_mask=[];
                                            percent_mask=logical((handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).perCorrPAC>=percent_windows(per_ii,1))...
                                                &(handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).perCorrPAC<=percent_windows(per_ii,2)));
                                            
                                            if ~isempty(percent_mask)
                                                
                                                
                                                for evNo=1:length(eventType)
                                                    
                                                    trials_in_event_Ev=[];
                                                    trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).which_eventPAC(eventType(evNo),:)==1)&percent_mask;
                                                    
                                                    if (sum(trials_in_event_Ev)>=1)
                                                        
                                                        %Do per bandwidth analysis
                                                        for pacii=1:no_pacii
                                                            
                                                            group_no=handles_drgb.drgbchoices.group_no(fileNo);
                                                            
                                                            %Enter the PACpower
                                                            this_peakPACpower_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_peakPACpower_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.meanPeakPower(trials_in_event_Ev)-handles_drgb.drgb.lfpevpair(lfpodRefNo).PACwave(pacii).PACwave.meanPeakPower(trials_in_event_Ev);
                                                            this_troughPACpower_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_troughPACpower_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.meanTroughPower(trials_in_event_Ev)-handles_drgb.drgb.lfpevpair(lfpodRefNo).PACwave(pacii).PACwave.meanTroughPower(trials_in_event_Ev);
                                                            %                                                             if (elec==which_electrodes(1))&(pacii==1)
                                                            %                                                                 this_deltaLickFreq_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            %                                                                 this_deltaLickFreq_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PACwave(pacii).PACwave.mean_lick_freq(trials_in_event_Ev)-handles_drgb.drgb.lfpevpair(lfpodRefNo).PACwave(pacii).PACwave.mean_lick_freq(trials_in_event_Ev);
                                                            %                                                             end
                                                            
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_peakPACpower_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1 ...
                                                                :theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=this_peakPACpower_Ev;
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_troughPACpower_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1 ...
                                                                :theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=this_troughPACpower_Ev;
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).whichMouse(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=mouseNo*ones(1,length(this_peakPACpower_Ev));
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev);
                                                            %                                                             if (elec==which_electrodes(1))&(pacii==1)
                                                            %                                                                  theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_deltaLickFreq_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+1 ...
                                                            %                                                                 :theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv+sum(trials_in_event_Ev))=this_deltaLickFreq_Ev;
                                                            %                                                             end
                                                            
                                                            %Save per session value for peak power
                                                            mean_PACpower_No=mean_PACpower_No+1;
                                                            mean_peakPACpower(mean_PACpower_No)=mean(this_peakPACpower_Ev);
                                                            mean_troughPACpower(mean_PACpower_No)=mean(this_troughPACpower_Ev);
                                                            mean_PACpower_perii(mean_PACpower_No)=per_ii;
                                                            mean_PACpower_evNo(mean_PACpower_No)=evNo;
                                                            mean_PACpower_pacii(mean_PACpower_No)=pacii;
                                                            mean_PACpower_fileNo(mean_PACpower_No)=fileNo;
                                                            mean_PACpower_mouse(mean_PACpower_No)=mouseNo;
                                                            per_session_group_no(mean_PACpower_No)=handles_drgb.drgbchoices.group_no(fileNo);
                                                            
                                                            
                                                            mouse_has_files=1;
                                                            
                                                        end
                                                        
                                                        
                                                        fprintf(1, ['%d trials in event No %d succesfully processed for file No %d electrode %d\n'],sum(trials_in_event_Ev), evNo,fileNo,elec);
                                                        
                                                    else
                                                        
                                                        
                                                        fprintf(1, ['%d trials in event No %d fewer than minimum trials per event ' evTypeLabels{evNo} ' for file No %d electrode %d\n'],sum(trials_in_event_Ev), min_trials_per_event,fileNo,elec);
                                                        
                                                        
                                                    end
                                                    
                                                    
                                                end
                                                
                                            else
                                                
                                                fprintf(1, ['Empty percent_mask for file No %d electrode %d\n'],fileNo,elec);
                                                
                                            end
                                            
                                        else
                                            
                                            fprintf(1, ['Empty PAC for file No %d electrode %d\n'],fileNo,elec);
                                            
                                        end
                                        
                                        
                                    else
                                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],fileNo,elec);
                                        
                                        
                                    end
                                end %if mouseNo
                            end
                        end %fileNo
                        
                        if theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).noEv>0
                            
                            %Calculate per mouse PAC power
                            for evNo=1:length(eventType)
                                for pacii=1:no_pacii
                                    for group_no=1:max(handles_drgb.drgbchoices.group_no)
                                        %Calculate per mouse PAC power
                                        mean_PACpower_No_per_mouse=mean_PACpower_No_per_mouse+1;
                                        this_mouse_peakPACpower=[];
                                        this_mouse_peakPACpower=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_peakPACpower_Ev;
                                        this_mouse_troughPACpower=[];
                                        this_mouse_troughPACpower=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_troughPACpower_Ev;
                                        if ~isempty(this_mouse_peakPACpower)
                                            mean_peakPACpower_per_mouse(mean_PACpower_No_per_mouse)=mean(this_mouse_peakPACpower);
                                            mean_troughPACpower_per_mouse(mean_PACpower_No_per_mouse)=mean(this_mouse_troughPACpower);
                                        else
                                            mean_peakPACpower_per_mouse(mean_PACpower_No_per_mouse)=NaN;
                                            mean_troughPACpower_per_mouse(mean_PACpower_No_per_mouse)=NaN;
                                        end
                                        
                                        mean_PACpower_perii_per_mouse(mean_PACpower_No_per_mouse)=per_ii;
                                        mean_PACpower_evNo_per_mouse(mean_PACpower_No_per_mouse)=evNo;
                                        mean_PACpower_pacii_per_mouse(mean_PACpower_No_per_mouse)=pacii;
                                        mean_PACpower_mouseNo_per_mouse(mean_PACpower_No_per_mouse)=mouseNo;
                                        mean_PACpower_electNo_per_mouse(mean_PACpower_No_per_mouse)=elec;
                                        mean_PACpower_group_no_per_mouse(mean_PACpower_No_per_mouse)=group_no;
                                        
                                        %                                         if (elec==which_electrodes(1))&(pacii==1)
                                        %                                             mean_deltaLickFreq_No_per_mouse=mean_deltaLickFreq_No_per_mouse+1;
                                        %                                             this_mouse_deltaLickFreq=[];
                                        %                                             this_mouse_deltaLickFreq=theseEvNos_thisMouse_thisElec(evNo,pacii,group_no).this_deltaLickFreq_Ev;
                                        %                                             if ~isempty(this_mouse_deltaLickFreq)
                                        %                                                 mean_deltaLickFreq_per_mouse(mean_deltaLickFreq_No_per_mouse)=mean(this_mouse_deltaLickFreq);
                                        %                                             else
                                        %                                                 mean_deltaLickFreq_per_mouse(mean_deltaLickFreq_No_per_mouse)=NaN;
                                        %                                             end
                                        %                                             mean_deltaLickFreq_perii_per_mouse(mean_deltaLickFreq_No_per_mouse)=per_ii;
                                        %                                             mean_deltaLickFreq_evNo_per_mouse(mean_deltaLickFreq_No_per_mouse)=evNo;
                                        %                                             mean_deltaLickFreq_pacii_per_mouse(mean_deltaLickFreq_No_per_mouse)=pacii;
                                        %                                             mean_deltaLickFreq_mouseNo_per_mouse(mean_deltaLickFreq_No_per_mouse)=mouseNo;
                                        %                                             mean_deltaLickFreq_electNo_per_mouse(mean_deltaLickFreq_No_per_mouse)=elec;
                                        %                                             mean_deltaLickFreq_group_no_per_mouse(mean_deltaLickFreq_No_per_mouse)=group_no;
                                        %                                         end
                                        
                                    end
                                end
                            end
                            
                            
                            
                            %Calculate per electrode ROC
                            can_calculate_auroc=1;
                            if can_calculate_auroc==1
                                for group_no=1:max(handles_drgb.drgbchoices.group_no)
                                    for evNo1=1:length(eventType)
                                        %                                         if theseEvNos(evNo1).noEv>0
                                        for evNo2=evNo1+1:length(eventType)
                                            %                                                 if theseEvNos(evNo2).noEv>0
                                            if (theseEvNos_thisMouse_thisElec(evNo1,pacii,group_no).noEv>0)&...
                                                    (theseEvNos_thisMouse_thisElec(evNo2,pacii,group_no).noEv>0)
                                                for pacii=1:no_pacii
                                                    
                                                    if (theseEvNos_thisMouse_thisElec(evNo1,pacii,group_no).noEv>=5)&...
                                                            (theseEvNos_thisMouse_thisElec(evNo2,pacii,group_no).noEv>=5)
                                                        
                                                        %Enter Ev1
                                                        this_mouse_peakPACpowerEv1=[];
                                                        this_mouse_peakPACpowerEv1=theseEvNos_thisMouse_thisElec(evNo1,pacii,group_no).this_peakPACpower_Ev;
                                                        this_mouse_troughPACpowerEv1=[];
                                                        this_mouse_troughPACpowerEv1=theseEvNos_thisMouse_thisElec(evNo1,pacii,group_no).this_troughPACpower_Ev;
                                                        trials_in_event_Ev1=theseEvNos_thisMouse_thisElec(evNo1,pacii,group_no).noEv;
                                                        
                                                        roc_data_peak=[];
                                                        roc_data_peak(1:trials_in_event_Ev1,1)=this_mouse_peakPACpowerEv1;
                                                        roc_data_peak(1:trials_in_event_Ev1,2)=zeros(trials_in_event_Ev1,1);
                                                        roc_data_trough=[];
                                                        roc_data_trough(1:trials_in_event_Ev1,1)=this_mouse_troughPACpowerEv1;
                                                        roc_data_trough(1:trials_in_event_Ev1,2)=zeros(trials_in_event_Ev1,1);
                                                        
                                                        %Enter Ev2
                                                        this_mouse_peakPACpowerEv2=[];
                                                        this_mouse_peakPACpowerEv2=theseEvNos_thisMouse_thisElec(evNo2,pacii,group_no).this_peakPACpower_Ev;
                                                        this_mouse_troughPACpowerEv2=[];
                                                        this_mouse_troughPACpowerEv2=theseEvNos_thisMouse_thisElec(evNo2,pacii,group_no).this_troughPACpower_Ev;
                                                        trials_in_event_Ev2=theseEvNos_thisMouse_thisElec(evNo2,pacii,group_no).noEv;
                                                        
                                                        
                                                        roc_data_peak(trials_in_event_Ev1+1:trials_in_event_Ev2+trials_in_event_Ev1,1)=this_mouse_peakPACpowerEv2;
                                                        roc_data_peak(trials_in_event_Ev1+1:trials_in_event_Ev2+trials_in_event_Ev1,2)=ones(trials_in_event_Ev2,1);
                                                        roc_data_trough(trials_in_event_Ev1+1:trials_in_event_Ev2+trials_in_event_Ev1,1)=this_mouse_troughPACpowerEv2;
                                                        roc_data_trough(trials_in_event_Ev1+1:trials_in_event_Ev2+trials_in_event_Ev1,2)=ones(trials_in_event_Ev2,1);
                                                        
                                                        
                                                        %Find  per electrode ROC
                                                        
                                                        no_ROCs=no_ROCs+1;
                                                        
                                                        
                                                        %                                                         ROCfileNo(no_ROCs)=handles_drgb.drgb.lfpevpair(lfpodNo_ref).fileNo;
                                                        ROCelec(no_ROCs)=elec;
                                                        ROCgroups(no_ROCs)=group_no;
                                                        ROCmouse(no_ROCs)=mouseNo;
                                                        ROCpacii(no_ROCs)=pacii;
                                                        ROCper_ii(no_ROCs)=per_ii;
                                                        ROCEvNo1(no_ROCs)=evNo1;
                                                        ROCEvNo2(no_ROCs)=evNo2;
                                                        
                                                        if ((abs(evNo1-2)<=1)&(abs(evNo2-5)<=1))||((abs(evNo1-5)<=1)&(abs(evNo2-2)<=1))
                                                            ROC_between(no_ROCs)=1;
                                                        else
                                                            ROC_between(no_ROCs)=0;
                                                        end
                                                        
                                                        roc=[];
                                                        roc=roc_calc(roc_data_peak,0,0.05,0);
                                                        auROCpeak(no_ROCs)=roc.AUC-0.5;
                                                        p_valROCpeak(no_ROCs)=roc.p;
                                                        p_vals_ROCpeak=[p_vals_ROCpeak roc.p];
                                                        
                                                        roc=[];
                                                        roc=roc_calc(roc_data_trough,0,0.05,0);
                                                        auROCtrough(no_ROCs)=roc.AUC-0.5;
                                                        p_valROCtrough(no_ROCs)=roc.p;
                                                        p_vals_ROCtrough=[p_vals_ROCtrough roc.p];
                                                        
                                                        %I have this code here to plot the ROC
                                                        
                                                        show_roc=0;
                                                        if (show_roc==1)&(mouseNo==4)
                                                            %I have this code here to plot the ROC
                                                            roc=roc_calc(roc_data_peak,0,0.05,1);
                                                            
                                                            %Do the histograms
                                                            try
                                                                close(2)
                                                            catch
                                                            end
                                                            figure(2)
                                                            
                                                            hold on
                                                            
                                                            max_dB=max([max(this_mouse_peakPACpowerEv1) max(this_mouse_peakPACpowerEv2)]);
                                                            min_dB=min([min(this_mouse_peakPACpowerEv1) min(this_mouse_peakPACpowerEv2)]);
                                                            
                                                            edges=[min_dB-0.1*(max_dB-min_dB):(max_dB-min_dB)/20:max_dB+0.1*(max_dB-min_dB)];
                                                            histogram(this_mouse_peakPACpowerEv1,edges,'FaceColor','b','EdgeColor','b')
                                                            histogram(this_mouse_peakPACpowerEv2,edges,'FaceColor','r','EdgeColor','r')
                                                            xlabel('delta power dB')
                                                            title(['Histogram for concentrations ' num2str(concs2(evNo1)) ' and ' num2str(concs2(evNo2))])
                                                            pffft=1;
                                                        end
                                                        
                                                    end
                                                end
                                                %
                                            end
                                        end
                                    end
                                    %                                         end
                                end
                                
                            end
                        end
                        
                    end
                    
                end
                
                if mouse_has_files==1
                    mouse_included(mouseNo)=1;
                else
                    mouse_included(mouseNo)=0;
                end
            end
            
            
        end
        fprintf(1, '\n\n')
        
        fprintf(1, ['The number of mice included in the PAC analysis for this odor pair is %d\n\n\n'], sum(mouse_included))
        
        figureNo = 0;
        %Now plot the average peakPACpower for each electrode calculated per mouse
        %(including all sessions for each mouse)
        edges=[-25:0.5:15];
        
        rand_offset=0.8;
        
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            data_PACpower=[];
            prof_naive=[];
            events=[];
            groups=[];
            mice=[];
            
            %Plot the average
            figureNo = figureNo + 1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            set(hFig, 'units','normalized','position',[.1 .2 .7 .7])
            %             subplot(2,1,1)
            ax=gca;ax.LineWidth=3;
            
            hold on
            
            bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            bar_offset = 0;
            
            %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
            
            include_group=0;
            
            for evNo=1:length(eventType)
                
                for per_ii=2:-1:1
                    
                    for grNo=1:max(handles_drgb.drgbchoices.group_no)
                        bar_offset = bar_offset +1;
                        %                         if sum(eventType==3)>0
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
                        %                         else
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        %                         end
                        
                        %                         these_offsets(per_ii)=bar_offset;
                        bar_offset = bar_offset + 1;
                        
                        if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
                            
                            include_group=1;
                            
                            switch grNo
                                case 1
                                    bar(bar_offset,mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'g','LineWidth', 3,'EdgeColor','none')
                                case 2
                                    bar(bar_offset,mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                                case 3
                                    bar(bar_offset,mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'y','LineWidth', 3,'EdgeColor','none')
                            end
                            
                            %Save data
                            handles_out.PRP_ii=handles_out.PRP_ii+1;
                            handles_out.PRP_values(handles_out.PRP_ii).pacii=pacii;
                            handles_out.PRP_values(handles_out.PRP_ii).evNo=evNo;
                            handles_out.PRP_values(handles_out.PRP_ii).per_ii=per_ii;
                            handles_out.PRP_values(handles_out.PRP_ii).groupNo=grNo;
                            handles_out.PRP_values(handles_out.PRP_ii).peak=1;
                            handles_out.PRP_values(handles_out.PRP_ii).PRP=mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)));
                            
                            %Save the mean per mouse
                            these_mice=mean_PACpower_mouseNo_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo));
                            handles_out.PRP_values(handles_out.PRP_ii).noMice=0;
                            for iiMice=min(these_mice):max(these_mice)
                                if sum(these_mice==iiMice)>0
                                    handles_out.PRP_values(handles_out.PRP_ii).noMice=handles_out.PRP_values(handles_out.PRP_ii).noMice+1;
                                    handles_out.PRP_values(handles_out.PRP_ii).mouseNo(handles_out.PRP_values(handles_out.PRP_ii).noMice)=iiMice;
                                    handles_out.PRP_values(handles_out.PRP_ii).PRP_per_mouse(handles_out.PRP_values(handles_out.PRP_ii).noMice)=mean(mean_peakPACpower_per_mouse((mean_PACpower_mouseNo_per_mouse==iiMice)&(~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)));
                                end
                            end
                            
                            %Violin plot
                            [mean_out, CIout]=drgViolinPoint(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))...
                                ,edges,bar_offset,rand_offset,'k','k',3);
                            
                            %                             plot(bar_offset,mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
                            %                             plot((bar_offset)*ones(1,sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),...
                            %                                 mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)),'o',...
                            %                                 'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            %
                            %                             if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>=2
                            %                                 CI = bootci(1000, {@mean, mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))},'type','cper');
                            %                                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            %                             end
                            
                            % %Save data for anovan
                            % data_PACpower=[data_PACpower mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
                            % prof_naive=[prof_naive per_ii*ones(1,sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            % events=[events evNo*ones(1,sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            % groups=[groups grNo*ones(1,sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            % mice=[mice mean_PACpower_mouseNo_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
                            
                        end
                    end
                    bar_offset = bar_offset + 2;
                    %                     if include_group==1
                    %                         bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                    %                         no_ev_labels=no_ev_labels+1;
                    %                         if sum(eventType==3)>0
                    %                             bar_labels{no_ev_labels}=evTypeLabels{evNo};
                    %                         else
                    %                             bar_labels{no_ev_labels}=num2str(concs(evNo));
                    %                         end
                    %                     end
                end
                bar_offset = bar_offset + 3;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%aqui
                %                 if include_group==1
                %                     ii_gr_included=ii_gr_included+1;
                %                     groups_included(ii_gr_included)=grNo;
                %                 end
                
            end
            
            title(['Peak PAC power theta/' freq_names{pacii+1}])
            
            
            %Annotations identifying groups
            x_interval=0.8/ii_gr_included;
            for ii=1:ii_gr_included
                annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %             %x labels
            %             to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
            %             sorted_A=sortrows(to_sort);
            %             sorted_bar_lab_loc=sorted_A(:,1);
            %             for ii=1:length(bar_lab_loc)
            %                 sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
            %             end
            xticks([2 4 6 10 12 14 21 23 25 29 31 33])
            xticklabels({'nwtS+', 'nHETS+', 'nKOS+', 'pwtS+', 'pHETS+', 'pKOS+', 'nwtS-', 'nHETS-', 'nKOS-', 'pwtS-', 'pHETS-', 'pKOS-'})
            
            if sum(eventType==3)==0
                xlabel('Concentration (%)')
            end
            
            ylabel('dB')
            
            %ylim1(pacii,:)=ylim;
            ylim([-7 5])
            
            
        end
        
        %         %Now plot the average delta lick frequency calculated per mouse
        %         %(including all sessions for each mouse)
        %
        %
        %         data_PACpower=[];
        %         prof_naive=[];
        %         events=[];
        %         groups=[];
        %         mice=[];
        %
        %         %Plot the average
        %         figureNo = figureNo + 1;
        %         try
        %             close(figureNo)
        %         catch
        %         end
        %         hFig=figure(figureNo);
        %
        %         set(hFig, 'units','normalized','position',[.1 .2 .7 .7])
        %         %             subplot(2,1,1)
        %
        %         hold on
        %
        %         bar_lab_loc=[];
        %         no_ev_labels=0;
        %         ii_gr_included=0;
        %         bar_offset = 0;
        %
        %         %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %         include_group=0;
        %
        %         for evNo=1:length(eventType)
        %
        %             for per_ii=2:-1:1
        %
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %                     bar_offset = bar_offset +1;
        %                     %                         if sum(eventType==3)>0
        %                     %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
        %                     %                         else
        %                     %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
        %                     %                         end
        %
        %                     %                         these_offsets(per_ii)=bar_offset;
        %                     bar_offset = bar_offset + 1;
        %
        %                     if sum((~isnan(mean_deltaLickFreq_per_mouse))&(mean_deltaLickFreq_perii_per_mouse==per_ii)&(mean_deltaLickFreq_pacii_per_mouse==pacii)&(mean_deltaLickFreq_evNo_per_mouse==evNo)&(mean_deltaLickFreq_group_no_per_mouse==grNo))>1
        %
        %                         include_group=1;
        %
        %                         switch grNo
        %                             case 1
        %                                 bar(bar_offset,mean(mean_deltaLickFreq_per_mouse((~isnan(mean_deltaLickFreq_per_mouse))&(mean_deltaLickFreq_perii_per_mouse==per_ii)&(mean_deltaLickFreq_pacii_per_mouse==pacii)&(mean_deltaLickFreq_evNo_per_mouse==evNo)&(mean_deltaLickFreq_group_no_per_mouse==grNo))),'g','LineWidth', 3,'EdgeColor','none')
        %                             case 2
        %                                 bar(bar_offset,mean(mean_deltaLickFreq_per_mouse((~isnan(mean_deltaLickFreq_per_mouse))&(mean_deltaLickFreq_perii_per_mouse==per_ii)&(mean_deltaLickFreq_pacii_per_mouse==pacii)&(mean_deltaLickFreq_evNo_per_mouse==evNo)&(mean_deltaLickFreq_group_no_per_mouse==grNo))),'b','LineWidth', 3,'EdgeColor','none')
        %                             case 3
        %                                 bar(bar_offset,mean(mean_deltaLickFreq_per_mouse((~isnan(mean_deltaLickFreq_per_mouse))&(mean_deltaLickFreq_perii_per_mouse==per_ii)&(mean_deltaLickFreq_pacii_per_mouse==pacii)&(mean_deltaLickFreq_evNo_per_mouse==evNo)&(mean_deltaLickFreq_group_no_per_mouse==grNo))),'y','LineWidth', 3,'EdgeColor','none')
        %                         end
        %
        %                         %Save data
        %                         handles_out.LickF_ii=handles_out.LickF_ii+1;
        %                         handles_out.LickF_values(handles_out.LickF_ii).pacii=pacii;
        %                         handles_out.LickF_values(handles_out.LickF_ii).evNo=evNo;
        %                         handles_out.LickF_values(handles_out.LickF_ii).groupNo=grNo;
        %                         handles_out.LickF_values(handles_out.LickF_ii).deltaLickF=mean(mean_deltaLickFreq_per_mouse((~isnan(mean_deltaLickFreq_per_mouse))&(mean_deltaLickFreq_perii_per_mouse==per_ii)&(mean_deltaLickFreq_pacii_per_mouse==pacii)&(mean_deltaLickFreq_evNo_per_mouse==evNo)&(mean_deltaLickFreq_group_no_per_mouse==grNo)));
        %
        %                         %Violin plot
        %                         [mean_out, CIout]=drgViolinPoint(mean_deltaLickFreq_per_mouse((~isnan(mean_deltaLickFreq_per_mouse))&(mean_deltaLickFreq_perii_per_mouse==per_ii)&(mean_deltaLickFreq_pacii_per_mouse==pacii)&(mean_deltaLickFreq_evNo_per_mouse==evNo)&(mean_deltaLickFreq_group_no_per_mouse==grNo))...
        %                             ,edges,bar_offset,rand_offset,'k','k',1);
        %
        %                         %                             plot(bar_offset,mean(mean_deltaLickFreq_per_mouse((~isnan(mean_deltaLickFreq_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
        %                         %                             plot((bar_offset)*ones(1,sum((~isnan(mean_deltaLickFreq_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),...
        %                         %                                 mean_deltaLickFreq_per_mouse((~isnan(mean_deltaLickFreq_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)),'o',...
        %                         %                                 'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
        %                         %
        %                         %                             if sum((~isnan(mean_deltaLickFreq_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>=2
        %                         %                                 CI = bootci(1000, {@mean, mean_deltaLickFreq_per_mouse((~isnan(mean_deltaLickFreq_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))},'type','cper');
        %                         %                                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        %                         %                             end
        %
        %                         % %Save data for anovan
        %                         % data_PACpower=[data_PACpower mean_deltaLickFreq_per_mouse((~isnan(mean_deltaLickFreq_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
        %                         % prof_naive=[prof_naive per_ii*ones(1,sum((~isnan(mean_deltaLickFreq_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
        %                         % events=[events evNo*ones(1,sum((~isnan(mean_deltaLickFreq_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
        %                         % groups=[groups grNo*ones(1,sum((~isnan(mean_deltaLickFreq_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
        %                         % mice=[mice mean_PACpower_mouseNo_per_mouse((~isnan(mean_deltaLickFreq_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
        %
        %                     end
        %                 end
        %                 bar_offset = bar_offset + 2;
        %                 %                     if include_group==1
        %                 %                         bar_lab_loc=[bar_lab_loc mean(these_offsets)];
        %                 %                         no_ev_labels=no_ev_labels+1;
        %                 %                         if sum(eventType==3)>0
        %                 %                             bar_labels{no_ev_labels}=evTypeLabels{evNo};
        %                 %                         else
        %                 %                             bar_labels{no_ev_labels}=num2str(concs(evNo));
        %                 %                         end
        %                 %                     end
        %             end
        %             bar_offset = bar_offset + 3;
        %
        %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%aqui
        %             %                 if include_group==1
        %             %                     ii_gr_included=ii_gr_included+1;
        %             %                     groups_included(ii_gr_included)=grNo;
        %             %                 end
        %
        %         end
        %
        %         title(['Delta lick frequency'])
        %
        %
        %         %Annotations identifying groups
        %         x_interval=0.8/ii_gr_included;
        %         for ii=1:ii_gr_included
        %             annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
        %         end
        %
        %         %Proficient/Naive annotations
        %         annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
        %         annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
        %
        %         %             %x labels
        %         %             to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
        %         %             sorted_A=sortrows(to_sort);
        %         %             sorted_bar_lab_loc=sorted_A(:,1);
        %         %             for ii=1:length(bar_lab_loc)
        %         %                 sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
        %         %             end
        %         xticks([2 4 6 10 12 14 21 23 25 29 31 33])
        %         xticklabels({'nwtS+', 'nHETS+', 'nKOS+', 'pwtS+', 'pHETS+', 'pKOS+', 'nwtS-', 'nHETS-', 'nKOS-', 'pwtS-', 'pHETS-', 'pKOS-'})
        %
        %         if sum(eventType==3)==0
        %             xlabel('Concentration (%)')
        %         end
        %
        %         ylabel('Hz')
        %
        
        
        
        
        
        
        %Now plot the average through PACpower for each electrode calculated per mouse
        %(including all sessions for each mouse)
        edges=[-25:0.5:15];
        
        rand_offset=0.8;
        
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            data_PACpower=[];
            prof_naive=[];
            events=[];
            groups=[];
            mice=[];
            
            %Plot the average
            figureNo = figureNo + 1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            set(hFig, 'units','normalized','position',[.1 .2 .7 .7])
            %             subplot(2,1,1)
            
            hold on
            
            bar_lab_loc=[];
            no_ev_labels=0;
            ii_gr_included=0;
            bar_offset = 0;
            
            %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
            
            include_group=0;
            
            for evNo=1:length(eventType)
                
                for per_ii=2:-1:1
                    
                    for grNo=1:max(handles_drgb.drgbchoices.group_no)
                        bar_offset = bar_offset +1;
                        %                         if sum(eventType==3)>0
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
                        %                         else
                        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        %                         end
                        
                        %                         these_offsets(per_ii)=bar_offset;
                        bar_offset = bar_offset + 1;
                        
                        if sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
                            
                            include_group=1;
                            
                            switch grNo
                                case 1
                                    bar(bar_offset,mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'g','LineWidth', 3,'EdgeColor','none')
                                case 2
                                    bar(bar_offset,mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                                case 3
                                    bar(bar_offset,mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'y','LineWidth', 3,'EdgeColor','none')
                            end
                            
                            %Save data
                            
                            %Save data
                            handles_out.PRP_ii=handles_out.PRP_ii+1;
                            handles_out.PRP_values(handles_out.PRP_ii).pacii=pacii;
                            handles_out.PRP_values(handles_out.PRP_ii).evNo=evNo;
                            handles_out.PRP_values(handles_out.PRP_ii).per_ii=per_ii;
                            handles_out.PRP_values(handles_out.PRP_ii).groupNo=grNo;
                            handles_out.PRP_values(handles_out.PRP_ii).peak=0;
                            handles_out.PRP_values(handles_out.PRP_ii).PRP=mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)));
                            
                            %Save the mean per mouse
                            these_mice=mean_PACpower_mouseNo_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo));
                            handles_out.PRP_values(handles_out.PRP_ii).noMice=0;
                            for iiMice=min(these_mice):max(these_mice)
                                if sum(these_mice==iiMice)>0
                                    handles_out.PRP_values(handles_out.PRP_ii).noMice=handles_out.PRP_values(handles_out.PRP_ii).noMice+1;
                                    handles_out.PRP_values(handles_out.PRP_ii).mouseNo(handles_out.PRP_values(handles_out.PRP_ii).noMice)=iiMice;
                                    handles_out.PRP_values(handles_out.PRP_ii).PRP_per_mouse(handles_out.PRP_values(handles_out.PRP_ii).noMice)=mean(mean_troughPACpower_per_mouse((mean_PACpower_mouseNo_per_mouse==iiMice)&(~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)));
                                end
                            end
                            
                            %Violin plot
                            [mean_out, CIout]=drgViolinPoint(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))...
                                ,edges,bar_offset,rand_offset,'k','k',3);
                            
                            
                            %                             plot(bar_offset,mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
                            %                             plot((bar_offset)*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))),...
                            %                                 mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)),'o',...
                            %                                 'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            %
                            %                             if sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>=2
                            %                                 CI = bootci(1000, {@mean, mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))},'type','cper');
                            %                                 plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            %                             end
                            
                            %                             %Save data for anovan
                            %                             data_PACpower=[data_PACpower mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
                            %                             prof_naive=[prof_naive per_ii*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            %                             events=[events evNo*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            %                             groups=[groups grNo*ones(1,sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))];
                            %                             mice=[mice mean_PACpower_mouseNo_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))];
                            %
                        end
                    end
                    bar_offset = bar_offset + 2;
                    %                     if include_group==1
                    %                         bar_lab_loc=[bar_lab_loc mean(these_offsets)];
                    %                         no_ev_labels=no_ev_labels+1;
                    %                         if sum(eventType==3)>0
                    %                             bar_labels{no_ev_labels}=evTypeLabels{evNo};
                    %                         else
                    %                             bar_labels{no_ev_labels}=num2str(concs(evNo));
                    %                         end
                    %                     end
                end
                bar_offset = bar_offset + 3;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%aqui
                %                 if include_group==1
                %                     ii_gr_included=ii_gr_included+1;
                %                     groups_included(ii_gr_included)=grNo;
                %                 end
                
            end
            
            title(['Trough PAC power theta/' freq_names{pacii+1}])
            
            
            %Annotations identifying groups
            x_interval=0.8/ii_gr_included;
            for ii=1:ii_gr_included
                annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
            end
            
            %Proficient/Naive annotations
            annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
            annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
            
            %x labels
            to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
            sorted_A=sortrows(to_sort);
            sorted_bar_lab_loc=sorted_A(:,1);
            for ii=1:length(bar_lab_loc)
                sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
            end
            xticks(sorted_bar_lab_loc)
            %             xticklabels(sorted_bar_labels)
            
            if sum(eventType==3)==0
                xlabel('Concentration (%)')
            end
            
            ylabel('dB')
            %             ylim1(pacii,:)=ylim;
            ylim([-7 5])
            
        end
        
        %         %Now do the cumulative histograms and ranksums for PAC power for each electrode calculated with all sessons per mouse
        %         pvals=[];
        %         legends=[];
        %
        % %         prof_naive_leg=handles_pars.per_lab;
        %
        %         wave_power=[];
        %
        %         for pacii=1:no_pacii
        %
        %             %Find max and min to use for xlim
        %             maxpower=-200;
        %             minpower=200;
        %             for evNo=1:length(eventType)
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %                         if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
        %                             maxpower=max([maxpower max(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))]);
        %                             minpower=min([minpower min(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)))]);
        %                         end
        %                     end
        %                 end
        %             end
        %              figureNo = figureNo + 1;
        %             try
        %                 close(figureNo)
        %             catch
        %             end
        %             hFig=figure(figureNo);
        %
        %
        %             set(hFig, 'units','normalized','position',[.1 .1 .7 .8])
        %
        %
        %             %Plot the trough PAC power
        %             ii_rank=0;
        %             input_data=[];
        %             ii_rankpt=0;
        %             input_datapt=[];
        %             glm_peakwave=[];
        %             glm_ii=0;
        %             glm_ptwave=[];
        %             glmpt_ii=0;
        %             for evNo=1:length(eventType)
        %
        %                 subplot(2,2,evNo)
        %                 hold on
        %
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                         if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
        %
        %                             [f_mi,x_mi] = drg_ecdf(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)));
        %                             cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_mi=f_mi;
        %                             cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_mi=x_mi;
        %
        %                             switch grNo
        %                             case 1
        %                             if per_ii==1
        %
        %                                legends.pacii(pacii).evNo(evNo).p1=plot(x_mi,f_mi,'Color',[0 0 1],'LineWidth',3);
        %                                 else
        %                                legends.pacii(pacii).evNo(evNo).p2=plot(x_mi,f_mi,'Color',[0.7 0.7 1],'LineWidth',3);
        %                             end
        %
        %                             case 2
        %                             if per_ii==1
        %
        %                                 legends.pacii(pacii).evNo(evNo).p3=plot(x_mi,f_mi,'Color',[1 0 0],'LineWidth',3);
        %                                 else
        %                                legends.pacii(pacii).evNo(evNo).p4=plot(x_mi,f_mi,'Color',[1 0.7 0.7],'LineWidth',3);
        %                             end
        %
        %                              case 3
        %                             if per_ii==1
        %
        %                                legends.pacii(pacii).evNo(evNo).p5=plot(x_mi,f_mi,'Color',[0 1 0],'LineWidth',3);
        %                                 else
        %                                legends.pacii(pacii).evNo(evNo).p6=plot(x_mi,f_mi,'Color',[0.7 1 0.7],'LineWidth',3);
        %                             end
        %
        %
        %
        % %                             if grNo==1
        % %                                 if per_ii==1
        % %                                     legends.pacii(pacii).evNo(evNo).p1=plot(x_mi,f_mi,'Color',[1 0 0],'LineWidth',3);
        % %                                 else
        % %                                     legends.pacii(pacii).evNo(evNo).p2=plot(x_mi,f_mi,'Color',[0 0 1],'LineWidth',3);
        % %                                 end
        % %                             else
        % %                                 if per_ii==1
        % %                                     legends.pacii(pacii).evNo(evNo).p3=plot(x_mi,f_mi,'Color',[1 0.7 0.7],'LineWidth',3);
        % %                                 else
        % %                                     legends.pacii(pacii).evNo(evNo).p4=plot(x_mi,f_mi,'Color',[0.7 0.7 1],'LineWidth',3);
        % %                                 end
        %                             end
        %
        %
        %                             %Compute per mouse avearge
        %                             each_mouse_average_PACpower=[];
        %                             no_mice_included=0;
        %                             for mouseNo=1:max(mean_PACpower_mouseNo_per_mouse)
        %                                 if mouse_included(mouseNo)==1
        %                                     if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
        %                                         if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo))>0
        %                                             no_mice_included=no_mice_included+1;
        %                                             each_mouse_average_PACpower(no_mice_included)=mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo)));
        %                                         end
        %                                     end
        %                                 end
        %                             end
        %
        %                             %Show average per mouse
        %                             if no_mice_included>0
        %                                 for jj=1:length(each_mouse_average_PACpower)
        %                                     this_f_mi=[];
        %                                     this_f_mi=cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_mi;
        %
        %                                     this_x_mi=[];
        %                                     this_x_mi=cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_mi;
        %
        %                                     xii_below=find(this_x_mi<each_mouse_average_PACpower(jj),1,'last');
        %                                     xii_above=find(this_x_mi>each_mouse_average_PACpower(jj),1,'first');
        %
        %                                     slope=(this_f_mi(xii_above)-this_f_mi(xii_below))/(this_x_mi(xii_above)-this_x_mi(xii_below));
        %                                     intercept=this_f_mi(xii_above)-slope*this_x_mi(xii_above);
        %
        %                                     this_f=slope*each_mouse_average_PACpower(jj)+intercept;
        %
        %
        %
        %                                     switch grNo
        %                             case 1
        %                             if per_ii==1
        %
        %                                    plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        %                                      else
        %                                     plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        %                             end
        %
        %                             case 2
        %                             if per_ii==1
        %                                   plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        %                                      else
        %                                     plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
        %                             end
        %
        %                              case 3
        %                             if per_ii==1
        %
        %                                    plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0 1 0],'MarkerEdge',[0 1 0],'MarkerSize',10)
        %                                      else
        %                                    plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0.7 1 0.7],'MarkerEdge',[0.7 1 0.7],'MarkerSize',10)
        %                             end
        %
        %
        % %                                     if grNo==1
        % %                                         if per_ii==1
        % %                                             plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
        % %                                         else
        % %                                             plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
        % %                                         end
        % %                                     else
        % %                                         if per_ii==1
        % %                                             plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
        % %                                         else
        % %                                             plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
        % %                                         end
        %                                     end
        %
        %                                 end
        %                             end
        %
        %                             %Save data for ranksum
        %                             ii_rank=ii_rank+1;
        %                             input_data(ii_rank).data=mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo));
        %                             input_data(ii_rank).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
        %                             input_data(ii_rank).per_ii=per_ii;
        %                             input_data(ii_rank).grNo=grNo;
        %                             input_data(ii_rank).evNo=evNo;
        %
        %                             ii_rankpt=ii_rankpt+1;
        %                             input_datapt(ii_rankpt).data=mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo));
        %                             input_datapt(ii_rankpt).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
        %                             input_datapt(ii_rankpt).per_ii=per_ii;
        %                             input_datapt(ii_rankpt).grNo=grNo;
        %                             input_datapt(ii_rankpt).evNo=evNo;
        %                             input_datapt(ii_rankpt).peak_trough=1;
        %                         end
        %                     end
        %
        %                 end
        %
        %                 title(evTypeLabels{evNo})
        %                 xlabel('Peak PAC power (dB)')
        %                 ylabel('Probability')
        %                 xlim([minpower-0.1*(maxpower-minpower) maxpower+0.1*(maxpower-minpower)])
        %             end
        %
        %
        %             for ii=1:ii_rank
        %                 glm_peakwave.data(glm_ii+1:glm_ii+length(input_data(ii).data))=input_data(ii).data;
        %                 glm_peakwave.group(glm_ii+1:glm_ii+length(input_data(ii).data))=input_data(ii).grNo;
        %                 glm_peakwave.perCorr(glm_ii+1:glm_ii+length(input_data(ii).data))=input_data(ii).per_ii;
        %                 glm_peakwave.event(glm_ii+1:glm_ii+length(input_data(ii).data))=input_data(ii).evNo;
        %                 glm_ii=glm_ii+length(input_data(ii).data);
        %             end
        %
        %             %Perform the glm
        %             fprintf(1, ['\nglm for wavelet power at the peak of theta for each electrode calculated per mouse for PAC theta ' freq_names{pacii+1} '\n'])
        %
        %             if sum(glm_peakwave.group==1)==length(glm_peakwave.group)
        %                 %There is only one group here (e.g. for Justin's paper we only include
        %                 %forward)
        %                 fprintf(1, ['\n\nglm for vector length for Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_peakwave.data',glm_peakwave.perCorr',glm_peakwave.event',...
        %                     'VariableNames',{'Peak_wave','perCorr','event'});
        %                 mdl = fitglm(tbl,'Peak_wave~perCorr+event+perCorr*event'...
        %                     ,'CategoricalVars',[2,3])
        %             else
        %
        %                 fprintf(1, ['\n\nglm for vector length for Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_peakwave.data',glm_peakwave.group',glm_peakwave.perCorr',glm_peakwave.event',...
        %                     'VariableNames',{'Peak_wave','group','perCorr','event'});
        %                 mdl = fitglm(tbl,'Peak_wave~group+perCorr+event+perCorr*group*event'...
        %                     ,'CategoricalVars',[2,3,4])
        %             end
        %
        %             fprintf(1, ['\n\nRanksum or t-test for  wavelet power at the peak of theta for each electrode calculated per mouse for PAC theta ' freq_names{pacii+1} '\n'])
        %             %Now do the ranksums
        %             output_data = drgMutiRanksumorTtest(input_data);
        %
        %             fprintf(1, ['\n\n'])
        %
        %
        %             %Plot the trough PAC power
        %             ii_rank=0;
        %             input_data=[];
        %             glm_troughwave=[];
        %             glm_ii=0;
        %             for evNo=1:length(eventType)
        %
        %                 subplot(2,2,evNo+2)
        %                 hold on
        %
        %                 for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %
        %
        %                     for per_ii=1:2      %performance bins. blue = naive, red = proficient
        %
        %
        %                         if sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
        %
        %                             [f_mi,x_mi] = drg_ecdf(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)));
        %                             cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_mi=f_mi;
        %                             cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_mi=x_mi;
        %
        %                             switch grNo
        %                             case 1
        %                             if per_ii==1
        %
        %                                      legends.pacii(pacii).evNo(evNo).p1=plot(x_mi,f_mi,'Color',[0 0 1],'LineWidth',3);
        %                                 else
        %                                     legends.pacii(pacii).evNo(evNo).p2=plot(x_mi,f_mi,'Color',[0.7 0.7 1],'LineWidth',3);
        %                                 end
        %
        %                             case 2
        %                             if per_ii==1
        %
        %                                   legends.pacii(pacii).evNo(evNo).p3=plot(x_mi,f_mi,'Color',[1 0 0],'LineWidth',3);
        %                                 else
        %                                   legends.pacii(pacii).evNo(evNo).p4=plot(x_mi,f_mi,'Color',[1 0.7 0.7],'LineWidth',3);
        %                                 end
        %
        %                              case 3
        %                             if per_ii==1
        %
        %                                   legends.pacii(pacii).evNo(evNo).p5=plot(x_mi,f_mi,'Color',[0 1 0],'LineWidth',3);
        %                                 else
        %                                   legends.pacii(pacii).evNo(evNo).p6=plot(x_mi,f_mi,'Color',[0.7 1 0.7],'LineWidth',3);
        %                             end
        %
        %
        % %                             if grNo==1
        % %                                 if per_ii==1
        % %                                     legends.pacii(pacii).evNo(evNo).p1=plot(x_mi,f_mi,'Color',[1 0 0],'LineWidth',3);
        % %                                 else
        % %                                     legends.pacii(pacii).evNo(evNo).p2=plot(x_mi,f_mi,'Color',[0 0 1],'LineWidth',3);
        % %                                 end
        % %                             else
        % %                                 if per_ii==1
        % %                                     legends.pacii(pacii).evNo(evNo).p3=plot(x_mi,f_mi,'Color',[1 0.7 0.7],'LineWidth',3);
        % %                                 else
        % %                                     legends.pacii(pacii).evNo(evNo).p4=plot(x_mi,f_mi,'Color',[0.7 0.7 1],'LineWidth',3);
        % %                                 end
        %                             end
        %
        %                             %Compute per mouse avearge
        %                             each_mouse_average_PACpower=[];
        %                             no_mice_included=0;
        %                             for mouseNo=1:max(mean_PACpower_mouseNo_per_mouse)
        %                                 if mouse_included(mouseNo)==1
        %                                     if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
        %                                         if sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo))>0
        %                                             no_mice_included=no_mice_included+1;
        %                                             each_mouse_average_PACpower(no_mice_included)=mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo)));
        %                                         end
        %                                     end
        %                                 end
        %                             end
        %
        %                             %Show average per mouse
        %                             if no_mice_included>0
        %                                 for jj=1:length(each_mouse_average_PACpower)
        %                                     this_f_mi=[];
        %                                     this_f_mi=cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_mi;
        %
        %                                     this_x_mi=[];
        %                                     this_x_mi=cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_mi;
        %
        %                                     xii_below=find(this_x_mi<each_mouse_average_PACpower(jj),1,'last');
        %                                     xii_above=find(this_x_mi>each_mouse_average_PACpower(jj),1,'first');
        %
        %                                     slope=(this_f_mi(xii_above)-this_f_mi(xii_below))/(this_x_mi(xii_above)-this_x_mi(xii_below));
        %                                     intercept=this_f_mi(xii_above)-slope*this_x_mi(xii_above);
        %
        %                                     this_f=slope*each_mouse_average_PACpower(jj)+intercept;
        %                                     p3=[];
        %
        %                                        switch grNo
        %                             case 1
        %                             if per_ii==1
        %
        %                                  p1=plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10);
        %                              else
        %                                  p2=plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10);
        %                             end
        %                             case 2
        %                             if per_ii==1
        %                                  p3=plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10);
        %                              else
        %                                  p4=plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10);
        %                             end
        %                              case 3
        %                             if per_ii==1
        %                                  p5=plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0 1 0 ],'MarkerEdge',[0 1 0],'MarkerSize',10);
        %                              else
        %                                  p6=plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0.7 1 0.7],'MarkerEdge',[0.7 1 0.7],'MarkerSize',10);
        %                             end
        %
        %
        %
        % %                                     if grNo==1
        % %                                         if per_ii==1
        % %                                             p1=plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10);
        % %                                         else
        % %                                             p2=plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10);
        % %                                         end
        % %                                     else
        % %                                         if per_ii==1
        % %                                             p3=plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10);
        % %                                         else
        % %                                             p4=plot(each_mouse_average_PACpower(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10);
        % %                                         end
        %                                     end
        %
        %                                 end
        %                             end
        %
        %                             %Save data for ranksum
        %                             ii_rank=ii_rank+1;
        %                             input_data(ii_rank).data=mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo));
        %                             input_data(ii_rank).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
        %                             input_data(ii_rank).per_ii=per_ii;
        %                             input_data(ii_rank).grNo=grNo;
        %                             input_data(ii_rank).evNo=evNo;
        %
        %                             ii_rankpt=ii_rankpt+1;
        %                             input_datapt(ii_rankpt).data=mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo));
        %                             input_datapt(ii_rankpt).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
        %                             input_datapt(ii_rankpt).per_ii=per_ii;
        %                             input_datapt(ii_rankpt).grNo=grNo;
        %                             input_datapt(ii_rankpt).evNo=evNo;
        %                             input_datapt(ii_rankpt).peak_trough=2;
        %                         end
        %                     end
        %
        %                 end
        %
        %                 title(evTypeLabels{evNo})
        %                 xlabel('Trough PAC power (dB)')
        %                 ylabel('Probability')
        %
        %                 if isempty(p3)
        %                     legend([p1 p2],{'Proficient','Naive'})
        %                 else
        %                     if (~isempty(p3))&(~isempty(p4))
        %
        %                            legend([p1 p2 p3 p4 p5 p6],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
        %                     [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'],...
        %                 [handles_drgb.drgbchoices.group_no_names{3} ' proficient'],[handles_drgb.drgbchoices.group_no_names{3} ' naive'])
        %
        % %                         legend([p1 p2 p3 p4],{['Proficient' handles_drgb.drgbchoices.group_no_names{1}],['Naive' handles_drgb.drgbchoices.group_no_names{1}],...
        % %                             ['Proficient' handles_drgb.drgbchoices.group_no_names{2}],['Naive' handles_drgb.drgbchoices.group_no_names{2}]})
        %                     end
        %                 end
        %                 xlim([minpower-0.1*(maxpower-minpower) maxpower+0.1*(maxpower-minpower)])
        %             end
        %
        %             suptitle(['Average PAC power for each electrode calculated per  for PAC theta/' freq_names{pacii+1}])
        %
        %
        %             %glm for trough
        %             for ii=1:ii_rank
        %                 glm_troughwave.data(glm_ii+1:glm_ii+length(input_data(ii).data))=input_data(ii).data;
        %                 glm_troughwave.group(glm_ii+1:glm_ii+length(input_data(ii).data))=input_data(ii).grNo;
        %                 glm_troughwave.perCorr(glm_ii+1:glm_ii+length(input_data(ii).data))=input_data(ii).per_ii;
        %                 glm_troughwave.event(glm_ii+1:glm_ii+length(input_data(ii).data))=input_data(ii).evNo;
        %                 glm_ii=glm_ii+length(input_data(ii).data);
        %             end
        %
        %             %Perform the glm
        %             fprintf(1, ['\nglm for wavelet power at the trough of theta for each electrode calculated per mouse for PAC theta ' freq_names{pacii+1} '\n'])
        %
        %             if sum(glm_troughwave.group==1)==length(glm_troughwave.group)
        %                 %There is only one group here (e.g. for Justin's paper we only include
        %                 %forward)
        %                 fprintf(1, ['\n\nglm for wavelet power at the trough of Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_troughwave.data',glm_troughwave.perCorr',glm_troughwave.event',...
        %                     'VariableNames',{'Trough_wave','perCorr','event'});
        %                 mdl = fitglm(tbl,'Trough_wave~perCorr+event+perCorr*event'...
        %                     ,'CategoricalVars',[2,3])
        %             else
        %
        %                 fprintf(1, ['\n\nglm for wavelet power at the trough of Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_troughwave.data',glm_troughwave.group',glm_troughwave.perCorr',glm_troughwave.event',...
        %                     'VariableNames',{'Trough_wave','group','perCorr','event'});
        %                 mdl = fitglm(tbl,'Trough_wave~group+perCorr+event+perCorr*group*event'...
        %                     ,'CategoricalVars',[2,3,4])
        %             end
        %
        %             %Now do the ranksums
        %             fprintf(1, ['\n\nRanksum or t-test p values for trough PAC power for each electrode calculated per mouse for PAC theta' freq_names{pacii+1} '\n'])
        %             output_data = drgMutiRanksumorTtest(input_data);
        %
        %             fprintf(1, ['\n\n'])
        %
        %             %glm for peak and trough
        %             %glm for trough
        %             for ii=1:ii_rankpt
        %                 glm_ptwave.data(glmpt_ii+1:glmpt_ii+length(input_datapt(ii).data))=input_datapt(ii).data;
        %                 glm_ptwave.group(glmpt_ii+1:glmpt_ii+length(input_datapt(ii).data))=input_datapt(ii).grNo;
        %                 glm_ptwave.perCorr(glmpt_ii+1:glmpt_ii+length(input_datapt(ii).data))=input_datapt(ii).per_ii;
        %                 glm_ptwave.event(glmpt_ii+1:glmpt_ii+length(input_datapt(ii).data))=input_datapt(ii).evNo;
        %                 glm_ptwave.peak_trough(glmpt_ii+1:glmpt_ii+length(input_datapt(ii).data))=input_datapt(ii).peak_trough;
        %                 glmpt_ii=glmpt_ii+length(input_datapt(ii).data);
        %             end
        %
        %             %Perform the glm
        %             fprintf(1, ['\nglm for wavelet power at the peak or trough of theta for each electrode calculated per mouse for PAC theta ' freq_names{pacii+1} '\n'])
        %
        %             if sum(glm_ptwave.group==1)==length(glm_ptwave.group)
        %                 %There is only one group here (e.g. for Justin's paper we only include
        %                 %forward)
        %                 fprintf(1, ['\n\nglm for wavelet power at the peak or trough of Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_ptwave.data',glm_ptwave.perCorr',glm_ptwave.event',glm_ptwave.peak_trough',...
        %                     'VariableNames',{'wave_power','perCorr','event','peak_vs_trough'});
        %                 mdl = fitglm(tbl,'wave_power~perCorr+event+peak_vs_trough+perCorr*event*peak_vs_trough'...
        %                     ,'CategoricalVars',[2,3])
        %             else
        %
        %                 fprintf(1, ['\n\nglm for wavelet power at the peak or trough of Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_ptwave.data',glm_ptwave.group',glm_ptwave.perCorr',glm_ptwave.event',glm_ptwave.peak_trough',...
        %                     'VariableNames',{'Trough_wave','group','perCorr','event','peak_vs_trough'});
        %                 mdl = fitglm(tbl,'Trough_wave~group+perCorr+event+peak_vs_trough+perCorr*group*event*peak_vs_trough'...
        %                     ,'CategoricalVars',[2,3,4])
        %             end
        %
        %             %Now do the ranksums
        %             fprintf(1, ['\n\nRanksum or t-test p values for peak or trough PAC power for each electrode calculated per mouse for PAC theta' freq_names{pacii+1} '\n'])
        %             output_data = drgMutiRanksumorTtest(input_datapt);
        %
        %             fprintf(1, ['\n\n'])
        %
        %             wave_power(pacii).input_datapt=input_datapt;
        %
        %         end
        
        
        %         %Now plot the average PAC power per mouse averaged over electrodes
        %         for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
        %
        %
        %             %Plot the average PAC power per mouse averaged over electrodes
        %              figureNo = figureNo + 1;
        %             try
        %                 close(figureNo)
        %             catch
        %             end
        %             hFig=figure(figureNo);
        %             set(hFig, 'units','normalized','position',[.1 .2 .7 .7])
        %
        %             %Peak PAC power
        %             subplot(2,1,1)
        %             hold on
        %
        %
        %             %             bar_lab_loc=[];
        %             no_ev_labels=0;
        %             ii_gr_included=0;
        %
        %             ii_rank=0;
        %             input_data=[];
        %
        %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %                 include_group=0;
        %
        %                 for evNo=1:length(eventType)
        %
        %                     for per_ii=1:2
        %
        %                         if sum(eventType==3)>0
        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
        %                         else
        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
        %                         end
        %
        %                         these_offsets(per_ii)=bar_offset;
        %
        %                         if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
        %
        %                             include_group=1;
        %
        %                             %Compute per mouse avearge
        %                             each_mouse_average_PACpower=[];
        %                             no_mice_included=0;
        %                             for mouseNo=1:max(mean_PACpower_mouseNo_per_mouse)
        %                                 if mouse_included(mouseNo)==1
        %                                     if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
        %                                         if sum((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo))>0
        %                                             no_mice_included=no_mice_included+1;
        %                                             each_mouse_average_PACpower(no_mice_included)=mean(mean_peakPACpower_per_mouse((~isnan(mean_peakPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo)));
        %                                         end
        %                                     end
        %                                 end
        %                             end
        %                             if no_mice_included>0
        %
        %                                 include_group=1;
        %
        %                                 if per_ii==1
        %                                     bar(bar_offset,mean(each_mouse_average_PACpower),'r','LineWidth', 3,'EdgeColor','none')
        %                                 else
        %                                     bar(bar_offset,mean(each_mouse_average_PACpower),'b','LineWidth', 3,'EdgeColor','none')
        %                                 end
        %
        %
        %                                 plot(bar_offset,mean(each_mouse_average_PACpower),'ok','LineWidth', 3)
        %                                 plot((bar_offset)*ones(1,length(each_mouse_average_PACpower)),each_mouse_average_PACpower,'o',...
        %                                     'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
        %
        %                                 if length(each_mouse_average_PACpower)>2
        %                                     CI = bootci(1000, {@mean, each_mouse_average_PACpower},'type','cper');
        %                                     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        %                                 end
        %
        %
        %
        %                                 %Save data for ranksum
        %                                 ii_rank=ii_rank+1;
        %                                 input_data(ii_rank).data=each_mouse_average_PACpower;
        %                                 input_data(ii_rank).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
        %
        %                             end
        %                         end
        %                     end
        %                 end
        %
        %
        %                 if include_group==1
        %                     ii_gr_included=ii_gr_included+1;
        %                     groups_included(ii_gr_included)=grNo;
        %                 end
        %             end
        %
        %             title('Peak PAC power')
        %
        %             %Annotations identifying groups
        %             x_interval=0.8/ii_gr_included;
        %             for ii=1:ii_gr_included
        %                 annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
        %             end
        %
        %             %Proficient/Naive annotations
        %             annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
        %             annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
        %
        %
        %             xticks(sorted_bar_lab_loc)
        % %             xticklabels(sorted_bar_labels)
        %
        %             if sum(eventType==3)==0
        %                 xlabel('Concentration (%)')
        %             end
        %
        %             ylabel('dB')
        %
        % %             %Now do the ranksums
        % %             fprintf(1, ['Ranksum or t-test p values for peak PAC power per mouse averaged over electrodes for PAC theta PAC theta' freq_names{pacii+1} '\n'])
        % %             output_data = drgMutiRanksumorTtest(input_data);
        % %
        %             %Trough power
        %             subplot(2,1,2)
        %             hold on
        %
        %
        %             %             bar_lab_loc=[];
        %             no_ev_labels=0;
        %             ii_gr_included=0;
        %
        %             ii_rank=0;
        %             input_data=[];
        %             glm_averagewave=[];
        %             glm_ii=0;
        %
        %             for grNo=1:max(handles_drgb.drgbchoices.group_no)
        %
        %                 include_group=0;
        %
        %                 for evNo=1:length(eventType)
        %
        %                     for per_ii=1:2
        %
        %                         if sum(eventType==3)>0
        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(2-evNo);
        %                         else
        %                             bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
        %                         end
        %
        %                         these_offsets(per_ii)=bar_offset;
        %
        %                         if sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo))>1
        %
        %                             include_group=1;
        %
        %                             %Compute per mouse avearge
        %                             each_mouse_average_PACpower=[];
        %                             no_mice_included=0;
        %                             for mouseNo=1:max(mean_PACpower_mouseNo_per_mouse)
        %                                 if mouse_included(mouseNo)==1
        %                                     if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
        %                                         if sum((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo))>0
        %                                             no_mice_included=no_mice_included+1;
        %                                             each_mouse_average_PACpower(no_mice_included)=mean(mean_troughPACpower_per_mouse((~isnan(mean_troughPACpower_per_mouse))&(mean_PACpower_perii_per_mouse==per_ii)&(mean_PACpower_pacii_per_mouse==pacii)&(mean_PACpower_evNo_per_mouse==evNo)&(mean_PACpower_group_no_per_mouse==grNo)&(mean_PACpower_mouseNo_per_mouse==mouseNo)));
        %                                         end
        %                                     end
        %                                 end
        %                             end
        %                             if no_mice_included>0
        %
        %                                 include_group=1;
        %
        %                                 if per_ii==1
        %                                     bar(bar_offset,mean(each_mouse_average_PACpower),'r','LineWidth', 3,'EdgeColor','none')
        %                                 else
        %                                     bar(bar_offset,mean(each_mouse_average_PACpower),'b','LineWidth', 3,'EdgeColor','none')
        %                                 end
        %
        %
        %                                 plot(bar_offset,mean(each_mouse_average_PACpower),'ok','LineWidth', 3)
        %                                 plot((bar_offset)*ones(1,length(each_mouse_average_PACpower)),each_mouse_average_PACpower,'o',...
        %                                     'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
        %
        %                                 if length(each_mouse_average_PACpower)>2
        %                                     CI = bootci(1000, {@mean, each_mouse_average_PACpower},'type','cper');
        %                                     plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        %                                 end
        %
        %
        %
        %                                 %Save data for ranksum
        %                                 ii_rank=ii_rank+1;
        %                                 input_data(ii_rank).data=each_mouse_average_PACpower;
        %                                 input_data(ii_rank).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evNo} ' ' prof_naive_leg{per_ii}];
        %                                 input_data(ii_rank).per_ii=per_ii;
        %                                 input_data(ii_rank).grNo=grNo;
        %                                 input_data(ii_rank).evNo=evNo;
        %                             end
        %                         end
        %                     end
        %                 end
        %
        %
        %                 if include_group==1
        %                     ii_gr_included=ii_gr_included+1;
        %                     groups_included(ii_gr_included)=grNo;
        %                 end
        %             end
        %
        %             title('Trough PAC power')
        %             suptitle(['Average PAC power per mouse averaged over all electrodes for PAC theta/' freq_names{pacii+1}])
        %
        %             %Annotations identifying groups
        %             x_interval=0.8/ii_gr_included;
        %             for ii=1:ii_gr_included
        %                 annotation('textbox',[0.7*x_interval+x_interval*(ii-1) 0.7 0.3 0.1],'String',handles_drgb.drgbchoices.group_no_names{ groups_included(ii)},'FitBoxToText','on');
        %             end
        %
        %             %Proficient/Naive annotations
        %             annotation('textbox',[0.15 0.8 0.3 0.1],'String','Proficient','FitBoxToText','on','Color','r','LineStyle','none');
        %             annotation('textbox',[0.15 0.75 0.3 0.1],'String','Naive','FitBoxToText','on','Color','b','LineStyle','none');
        %
        %             %             %x labels
        %             %             to_sort=[bar_lab_loc' [1:length(bar_lab_loc)]'];
        %             %             sorted_A=sortrows(to_sort);
        %             %             sorted_bar_lab_loc=sorted_A(:,1);
        %             %             for ii=1:length(bar_lab_loc)
        %             %                 sorted_bar_labels{ii}=bar_labels{sorted_A(ii,2)};
        %             %             end
        %             xticks(sorted_bar_lab_loc)
        % %             xticklabels(sorted_bar_labels)
        %
        %             if sum(eventType==3)==0
        %                 xlabel('Concentration (%)')
        %             end
        %
        %             ylabel('dB')
        %
        %
        %
        %             for ii=1:ii_rank
        %                 glm_averagewave.data(glm_ii+1:glm_ii+length(input_data(ii).data))=input_data(ii).data;
        %                 glm_averagewave.group(glm_ii+1:glm_ii+length(input_data(ii).data))=input_data(ii).grNo;
        %                 glm_averagewave.perCorr(glm_ii+1:glm_ii+length(input_data(ii).data))=input_data(ii).per_ii;
        %                 glm_averagewave.event(glm_ii+1:glm_ii+length(input_data(ii).data))=input_data(ii).evNo;
        %                 glm_ii=glm_ii+length(input_data(ii).data);
        %             end
        %
        %             %Perform the glm
        %             fprintf(1, ['\nglm for average wavelet power for each electrode calculated per mouse ' freq_names{pacii+1} '\n'])
        %
        %             if sum(glm_averagewave.group==1)==length(glm_averagewave.group)
        %                 %There is only one group here (e.g. for Justin's paper we only include
        %                 %forward)
        %                 fprintf(1, ['\n\nglm for average wavelet power for Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_averagewave.data',glm_averagewave.perCorr',glm_averagewave.event',...
        %                     'VariableNames',{'Average_wave','perCorr','event'});
        %                 mdl = fitglm(tbl,'Average_wave~perCorr+event+perCorr*event'...
        %                     ,'CategoricalVars',[2,3])
        %             else
        %
        %                 fprintf(1, ['\n\nglm for average wavelet power for Theta/' freq_names{pacii+1} '\n'])
        %                 tbl = table(glm_averagewave.data',glm_averagewave.group',glm_averagewave.perCorr',glm_averagewave.event',...
        %                     'VariableNames',{'Average_wave','group','perCorr','event'});
        %                 mdl = fitglm(tbl,'Average_wave~group+perCorr+event+perCorr*group*event'...
        %                     ,'CategoricalVars',[2,3,4])
        %             end
        %
        %             %Now do the ranksums
        %             fprintf(1, ['Ranksum or t-test p values for trough PAC power per mouse averaged over electrodes for PAC theta PAC theta' freq_names{pacii+1} '\n'])
        %             output_data = drgMutiRanksumorTtest(input_data);
        %
        %         end
        %
        
        %Display auROC
        edges=[-0.3:0.05:0.5];
        rand_offset=0.8;
        for pacii=1:no_pacii    %for different PACs
            
            ii_roc=0;
            roc_data=[];
            glm_roc=[];
            glm_roc_ii=0;
            
            ii_rocpk=0;
            rocpk_data=[];
            glm_rocpk=[];
            glm_rocpk_ii=0;
            
            ii_roctr=0;
            roctr_data=[];
            glm_roctr=[];
            glm_roctr_ii=0;
            
            %Display the peak auROC
            figureNo=figureNo+1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            
            
            bar_offset=0;
            
            for per_ii=2:-1:1
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    
                    if sum((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))>0
                        
                        bar_offset=bar_offset+1;
                        
                        switch grNo
                            case 1
                                bar(bar_offset,mean(auROCpeak((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))),'g','LineWidth', 3,'EdgeColor','none')
                            case 2
                                bar(bar_offset,mean(auROCpeak((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                            case 3
                                bar(bar_offset,mean(auROCpeak((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))),'y','LineWidth', 3,'EdgeColor','none')
                        end
                        
                        
                        %Save data
                        handles_out.AUC_ii=handles_out.AUC_ii+1;
                        handles_out.AUC_values(handles_out.AUC_ii).pacii=pacii;
                        handles_out.AUC_values(handles_out.AUC_ii).evNo=evNo;
                        handles_out.AUC_values(handles_out.AUC_ii).per_ii=per_ii;
                        handles_out.AUC_values(handles_out.AUC_ii).groupNo=grNo;
                        handles_out.AUC_values(handles_out.AUC_ii).peak=1;
                        handles_out.AUC_values(handles_out.AUC_ii).AUC=mean(auROCpeak((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo)));
                        
                        %Violin plot
                        [mean_out, CIout]=drgViolinPoint(auROCpeak((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))...
                            ,edges,bar_offset,rand_offset,'k','k',3);
                        
                        
                        %Enter data for glm for peak only
                        these_data=auROCpeak((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))';
                        glm_rocpk.data(glm_rocpk_ii+1:glm_rocpk_ii+length(these_data))=these_data;
                        glm_rocpk.group(glm_rocpk_ii+1:glm_rocpk_ii+length(these_data))=grNo;
                        glm_rocpk.perCorr(glm_rocpk_ii+1:glm_rocpk_ii+length(these_data))=per_ii;
                        glm_rocpk_ii=glm_rocpk_ii+length(these_data);
                        
                        %Enter the data for t-test/ranksum
                        ii_rocpk=ii_rocpk+1;
                        rocpk_data(ii_rocpk).data=these_data;
                        rocpk_data(ii_rocpk).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' prof_naive_leg{per_ii}];
                        
                        
                        %Enter data for glm for peak and trough only
                        glm_roc.data(glm_roc_ii+1:glm_roc_ii+length(these_data))=these_data;
                        glm_roc.group(glm_roc_ii+1:glm_roc_ii+length(these_data))=grNo;
                        glm_roc.perCorr(glm_roc_ii+1:glm_roc_ii+length(these_data))=per_ii;
                        glm_roc.peak_trough(glm_roc_ii+1:glm_roc_ii+length(these_data))=1;
                        glm_roc_ii=glm_roc_ii+length(these_data);
                        
                        %Enter the data for t-test/ranksum
                        ii_roc=ii_roc+1;
                        roc_data(ii_roc).data=these_data;
                        roc_data(ii_roc).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' prof_naive_leg{per_ii} ' peak'];
                        
                        
                        
                        
                    end
                    
                    
                end
                
                
                bar_offset=bar_offset+1;
            end
            
            title(['auROC per mouse, per electrode peak-referenced power theta/' freq_names{pacii+1}])
            
            
            
            ylabel('auROC')
            
            pffft=1;
            
            
            
            %Display the trough auROC for all trials per mouse (per electrode) for within vs
            %betweeen
            
            %Plot the average
            figureNo=figureNo+1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);
            
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            
            
            bar_offset=0;
            
            for per_ii=2:-1:1
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    
                    if sum((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))>0
                        
                        bar_offset=bar_offset+1;
                        
                        switch grNo
                            case 1
                                bar(bar_offset,mean(auROCtrough((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))),'g','LineWidth', 3,'EdgeColor','none')
                            case 2
                                bar(bar_offset,mean(auROCtrough((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                            case 3
                                bar(bar_offset,mean(auROCtrough((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))),'y','LineWidth', 3,'EdgeColor','none')
                        end
                        
                        
                        %Save data
                        handles_out.AUC_ii=handles_out.AUC_ii+1;
                        handles_out.AUC_values(handles_out.AUC_ii).pacii=pacii;
                        handles_out.AUC_values(handles_out.AUC_ii).evNo=evNo;
                        handles_out.AUC_values(handles_out.AUC_ii).per_ii=per_ii;
                        handles_out.AUC_values(handles_out.AUC_ii).groupNo=grNo;
                        handles_out.AUC_values(handles_out.AUC_ii).peak=0;
                        handles_out.AUC_values(handles_out.AUC_ii).AUC=mean(auROCtrough((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo)));
                        
                        
                        %Violin plot
                        
                        [mean_out, CIout]=drgViolinPoint(auROCtrough((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))...
                            ,edges,bar_offset,rand_offset,'k','k',3);
                        
                        
                        %Enter data for glm for peak only
                        these_data=auROCtrough((ROCper_ii==per_ii)&(ROCpacii==pacii)&(ROCgroups==grNo))';
                        glm_roctr.data(glm_roctr_ii+1:glm_roctr_ii+length(these_data))=these_data;
                        glm_roctr.group(glm_roctr_ii+1:glm_roctr_ii+length(these_data))=grNo;
                        glm_roctr.perCorr(glm_roctr_ii+1:glm_roctr_ii+length(these_data))=per_ii;
                        glm_roctr_ii=glm_roctr_ii+length(these_data);
                        
                        %Enter the data for t-test/ranksum
                        ii_roctr=ii_roctr+1;
                        roctr_data(ii_roctr).data=these_data;
                        roctr_data(ii_roctr).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' prof_naive_leg{per_ii}];
                        
                        
                        %Enter data for glm for peak and trough only
                        glm_roc.data(glm_roc_ii+1:glm_roc_ii+length(these_data))=these_data;
                        glm_roc.group(glm_roc_ii+1:glm_roc_ii+length(these_data))=grNo;
                        glm_roc.perCorr(glm_roc_ii+1:glm_roc_ii+length(these_data))=per_ii;
                        glm_roc.peak_trough(glm_roc_ii+1:glm_roc_ii+length(these_data))=2;
                        glm_roc_ii=glm_roc_ii+length(these_data);
                        
                        %Enter the data for t-test/ranksum
                        ii_roc=ii_roc+1;
                        roc_data(ii_roc).data=these_data;
                        roc_data(ii_roc).description=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' prof_naive_leg{per_ii} ' peak'];
                        
                        
                        
                    end
                    
                    
                end
                
                
                bar_offset=bar_offset+1;
            end
            
            title(['auROC per mouse, per electrode trough-referenced power theta/' freq_names{pacii+1}])
            
            
            
            ylabel('auROC')
            
            %Do glm for peak/trough
            fprintf(1, ['\n\nglm for auROC peak/trough for Theta/' freq_names{pacii+1} '\n'])
            tbl = table(glm_roc.data',glm_roc.group',glm_roc.perCorr',glm_roc.peak_trough',...
                'VariableNames',{'Peak_wave','group','perCorr','peak_trough'});
            mdl = fitglm(tbl,'Peak_wave~group+perCorr+peak_trough+group*perCorr*peak_trough'...
                ,'CategoricalVars',[2,3,4])
            
            
            fprintf(1, ['\n\nRanksum or t-test for auROC peak for theta ' freq_names{pacii+1} '\n'])
            %Now do the ranksums
            output_data = drgMutiRanksumorTtest(rocpk_data);
            
            fprintf(1, ['\n\nRanksum or t-test for auROC trough for theta ' freq_names{pacii+1} '\n'])
            %Now do the ranksums
            output_data = drgMutiRanksumorTtest(roctr_data);
            
            pffft=1;
            
        end
        
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) handles_pars.output_suffix],'handles_out')
        pfft=1;
        
        
        
        
end


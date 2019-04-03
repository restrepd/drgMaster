function drgAnalysisBatchLFP(handles)

%drgAnalysisBatchLFP displays the LFP power spectrum for drgRunBatch-generated
%data

%Which analysis is performed is determined by the value enterd in the
%variable which_display:
%
%
% 1 ERP analysis compare auROC in the last few trials of pre with first few trials of post. Analyzed per session
%
% 2 Generate delta LFP power and auROC for reversal figure for Daniel's paper. Analyzed per session
%
% 3 For a subset of first files for events 1 and 2 plot the LFP bandwide spectrum,
%   LFP power histograms for different bandwidths for each electrode and LFP auROCs.
%   To generate Fig. 2 for Daniels' LFP power paper enter the proficient files
%   Analyzed per session
%
% 4 Generate LFP power auROC for Fig. 3 for Daniel's paper. first vs last. Analyzed per session
%
%
% 5 Compare auROC in the last few trials of pre with first few trials of post
%    Used for old Fig. 4 of Daniel's paper with acetoethylben_electrode9202017.mat. Analyzed per session
%
% 6 For a subset of first files for events 1 and 2 plot the ERP LFP bandwide spectrum,
%   ERP LFP power histograms for different bandwidths for each electrode and ERP LFP auROCs.
%   Analyzed per session
%
% 7 Generate ERP LFP power auROC. first vs last. . Analyzed per session
%
% 8 Compare auROC for ERP LFP in the last few trials of pre with first few trials of post
%   Used for New Fig. 7 of Daniel's paper. Analyzed per session
%
% 9 Compare auROC for power LFP for two events in two percent windows for
%   all of the files. Analyzed per session
%
% 10 Compare auROC for power LFP for two groups (e.g. NRG1 vs control)
%    within one precent window. Analyzed per session
%
% 11 Compare auROC for ERP LFP powerin between two percent correct windows. Analyzed per session
%
% 12 Justin's analysis of LFP power differences for naive and proficient
% mice. Analyzed per session
%
% 13 Compare auROC for ERP wavelet LFP power between two percent correct
%     windows at different ERP delta t values. Analyzed per session
%
% 14  Justin's multiclass ROC analysis of LFP power differences for naive and proficient
% mice for different epochs (concentrations or S+ vs. S-). Analyzed per mouse
%
% 15  Justin's fitlm analysis of LFP power differences for naive and proficient
% mice. Analyzed per mouse
%
% 16  Justin's fitglm analysis of LFP power differences for naive and proficient
% mice. Analyzed per mouse
%
% 17  Justin's PAC analysis for events (concentrations or S+/S-) for naive and proficient
% Analyzed per mouse (all trials for all sessions for each mouse are used
% in the analysis)
%
% 18 Multiclass ROC analysis for mean angle for PAC for multiple epochs and proficient
% vs naive
%
% 19 PAC MI analysis for events (concentrations or S+/S-) for naive and proficient
% Analyzed per mouse for groups defined by the user
%
% 20  Multiclass ROC analysis of LFP power differences for naive and proficient
% mice for different epochs (concentrations or S+ vs. S-) and different groups. Analyzed per mouse
%
% 21  Multiclass ROC analysis of coherence for naive and proficient
%  mice for different epochs (concentrations or S+ vs. S-) and different groups. Analyzed per mouse
%
% 22 ERWA analysis for LFP power 
%
% 23 Oscillatory power at the peak and trough of the PAC

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
    freq_names={'Theta','Beta','Low gamma','High gamma'};
else
    no_bandwidths=handles_pars.no_bandwidths;
    low_freq=handles_pars.low_freq;
    high_freq=handles_pars.high_freq;
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
load([handles.PathName handles.drgb.outFileName])

if (~isfield(handles_pars,'files'))||(isempty(handles_pars.files))
    files=[1:handles_drgb.drgbchoices.no_files];
end

if isfield(handles_pars,'groupNo')
    handles_drgb.drgbchoices.group_no=handles_pars.groupNo;
end

fprintf(1, ['\ndrgDisplayBatchLFPPowerPairwise run for ' handles.drgb.outFileName '\nwhich_display= = %d\n\n'],which_display);

switch which_display
    
    case {1,6,7,8,11,13}
        
        frequency=handles_drgb.drgb.lfpevpair(1).fERP;
        max_events_per_sec=(handles_drgb.drgbchoices.timeEnd(winNo)-handles_drgb.drgbchoices.timeStart(winNo))*handles_drgb.max_events_per_sec;
    case 21
        frequency=handles_drgb.drgb.lfpevpair(1).f_coh;
    case 22
        frequency=handles_drgb.drgb.lfpevpair(1).wave_fERWA;
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
                                    [p_val(no_dBs,bwii), r_or_t]=drg_ranksum_or_ttest(this_delta_dB_powerpreHit,this_delta_dB_powerpostHit);
                                    p_vals=[p_vals p_val(no_dBs,bwii)];
                                    groupNopre(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                    groupNopost(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                    events(no_dBs)=1;
                                    
                                    
                                    %CR
                                    [p_val(no_dBs+1,bwii),r_or_t]=drg_ranksum_or_ttest(this_delta_dB_powerpreCR,this_delta_dB_powerpostCR);
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
        for bwii=1:no_bandwidths
            
            
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
        %             for bwii=1:no_bandwidths
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
        
        
        for bwii=1:no_bandwidths
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
                [pval_auROCperm, r_or_t]=drg_ranksum_or_ttest(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)), auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                
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
                                    
                                    [p_val(no_dBs,bwii),r_or_t]=drg_ranksum_or_ttest(this_delta_dB_powerfp1Ev1,this_delta_dB_powerfp2Ev1);
                                    p_vals=[p_vals p_val(no_dBs,bwii)];
                                    groupNofp1(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                    groupNofp2(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                    events(no_dBs)=1;
                                    
                                    
                                    %Ev2
                                    [p_val(no_dBs+1,bwii),r_or_t]=drg_ranksum_or_ttest(this_delta_dB_powerfp1Ev2,this_delta_dB_powerfp2Ev2);
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
        for bwii=1:no_bandwidths
            
            
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
        
        for bwii=1:no_bandwidths
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
                    lfppairNo=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfppairNo)))
                        
                        
                        if (~isempty(handles_drgb.drgb.lfpevpair(lfppairNo).allPower))
                            
                            if (length(handles_drgb.drgb.lfpevpair(lfppairNo).which_eventLFPPower(1,:))>=trials_to_process)
                                
                                trials=length(handles_drgb.drgb.lfpevpair(lfppairNo).which_eventLFPPower(1,:));
                                mask=logical([zeros(1,trials-trials_to_process) ones(1,trials_to_process)]);
                                trials_in_eventEv1=(handles_drgb.drgb.lfpevpair(lfppairNo).which_eventLFPPower(event1,:)==1);
                                trials_in_eventEv2=(handles_drgb.drgb.lfpevpair(lfppairNo).which_eventLFPPower(event2,:)==1);
                                
                                if (sum(trials_in_eventEv1)>=min_trials_per_event) & (sum(trials_in_eventEv2)>=min_trials_per_event)
                                    
                                    lfpodNo=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                    
                                    % Ev1
                                    this_dB_powerrefEv1=zeros(sum(trials_in_eventEv1&mask),length(frequency));
                                    this_dB_powerrefEv1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfppairNo).allPower(trials_in_eventEv1&mask,:));
                                    
                                    
                                    this_dB_powerEv1=zeros(sum(trials_in_eventEv1&mask),length(frequency));
                                    this_dB_powerEv1(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_eventEv1&mask,:));
                                    
                                    % Ev2
                                    this_dB_powerrefEv2=zeros(sum(trials_in_eventEv2&mask),length(frequency));
                                    this_dB_powerrefEv2(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfppairNo).allPower(trials_in_eventEv2&mask,:));
                                    
                                    
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
                                        ROCout(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfppairNo).fileNo;
                                        ROCgroupNo(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfppairNo).fileNo);
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
                                
                                fprintf(1, ['%d trials fewer than %d trials to process for file No %d electrode %d\n'],length(handles_drgb.drgb.lfpevpair(lfppairNo).which_eventLFPPower(1,:)),trials_to_process,files(fileNo),elec);
                                
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
        for bwii=1:no_bandwidths
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
        
        for bwii=1:no_bandwidths
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
            legend('auROC not significant','auROC significant')
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
        
        for bwii=1:no_bandwidths
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
                                        
                                        [p_val(no_dBs,bwii),r_or_t]=drg_ranksum_or_ttest(this_delta_dB_powerfp1Ev1,this_delta_dB_powerfp2Ev1);
                                        p_vals=[p_vals p_val(no_dBs,bwii)];
                                        groupNofp1(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                        groupNofp2(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                        events(no_dBs)=1;
                                        
                                        
                                        %Ev2
                                        [p_val(no_dBs+1,bwii),r_or_t]=drg_ranksum_or_ttest(this_delta_dB_powerfp1Ev2,this_delta_dB_powerfp2Ev2);
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
        
        for bwii=1:no_bandwidths
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
        title('Percent significant auROC')
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
                                    
                                    [p_val(no_dBs,bwii),r_or_t]=drg_ranksum_or_ttest(this_delta_dB_powerpreHit,this_delta_dB_powerpostHit);
                                    p_vals=[p_vals p_val(no_dBs,bwii)];
                                    groupNopre(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                    groupNopost(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                    events(no_dBs)=1;
                                    
                                    
                                    %CR
                                    [p_val(no_dBs+1,bwii),r_or_t]=drg_ranksum_or_ttest(this_delta_dB_powerpreCR,this_delta_dB_powerpostCR);
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
        for bwii=1:no_bandwidths
            
            
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
        %             for bwii=1:no_bandwidths
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
        for bwii=1:no_bandwidths
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
                [pval_auROCperm,r_or_t]=drg_ranksum_or_ttest(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)), auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                
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
        for bwii=1:no_bandwidths
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
        
        for bwii=1:no_bandwidths
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
            legend('auROC not significant','auROC significant')
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
        
        for bwii=1:no_bandwidths
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
                                    
                                    [p_val(no_dBs,bwii),r_or_t]=drg_ranksum_or_ttest(this_delta_dB_powerfp1Ev1,this_delta_dB_powerfp2Ev1);
                                    p_vals=[p_vals p_val(no_dBs,bwii)];
                                    groupNofp1(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                    groupNofp2(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                    events(no_dBs)=1;
                                    
                                    
                                    %Ev2
                                    [p_val(no_dBs+1,bwii),r_or_t]=drg_ranksum_or_ttest(this_delta_dB_powerfp1Ev2,this_delta_dB_powerfp2Ev2);
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
        
        for bwii=1:no_bandwidths
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
        title('Percent significant auROC')
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
                                        [p_val(no_dBs,bwii),r_or_t]=drg_ranksum_or_ttest(this_delta_dB_powerpreHit,this_delta_dB_powerpostHit);
                                        p_vals=[p_vals p_val(no_dBs,bwii)];
                                        groupNopre(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                        groupNopost(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                        events(no_dBs)=1;
                                        
                                        
                                        %CR
                                        [p_val(no_dBs+1,bwii), r_or_t]=drg_ranksum_or_ttest(this_delta_dB_powerpreCR,this_delta_dB_powerpostCR);
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
            
            
            for bwii=1:no_bandwidths
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
        for bwii=1:no_bandwidths
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
        %Compare LFP power auROC for two events in two percent windows for all of the files
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
                                percent_mask=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower>=percent_windows(per_ii,1))...
                                    &(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower<=percent_windows(per_ii,2));
                                trials_in_event_Ev1=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(event1,:)==1)&percent_mask;
                                trials_in_event_Ev2=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(event2,:)==1)&percent_mask;
                                
                                
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
        
        
        for bwii=1:no_bandwidths
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
            %             legend(file_label{1},file_label{2})
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
        title('Percent significant auROC')
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
        
        
        for bwii=1:no_bandwidths
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
        title('Percent significant auROC')
        legend(group_names{1},group_names{2})
        ylim([0 100])
        
        pFDRanovanauROC=drsFDRpval(pval_auROCperm);
        fprintf(1, ['\npFDR for premuted anovan p value for difference between\n ' group_names{1} ' and ' group_names{2} ' for auROC = %d\n\n'],pFDRanovanauROC);
        
        
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix],'learn_sig','learn_not_sig','prof_sig','prof_not_sig');
        
        
        pffft=1
        
    case 11
        
        %Compare auROC for ERP LFP powerin between two percent correct windows
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
        
        
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        sz_per=size(percent_windows);
        no_per_win=sz_per(1);
        
        
        
        fprintf(1, ['Pairwise auROC log power LFP ERP analysis for ' evTypeLabels{1} ' and ' evTypeLabels{2} ' LFP power\n\n'])
        p_vals=[];
        
        no_files=length(files);
        
        for fileNo=1:no_files
            
            
            for elec=1:16
                if sum(which_electrodes==elec)>0
                    lfpodNo=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                    
                    if ~isempty(handles_drgb.drgb.lfpevpair(lfpodNo))
                        
                        if ~isempty(handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP)
                            
                            
                            for per_corr_ii=1:no_per_win
                                
                                these_per_corr=(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrERP>=percent_windows(per_corr_ii,1))&...
                                    (handles_drgb.drgb.lfpevpair(lfpodNo).perCorrERP<=percent_windows(per_corr_ii,2));
                                
                                %Which trials have events (i.e. licks)?
                                trials_with_event=(handles_drgb.drgb.lfpevpair(lfpodNo).no_events_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNo).no_events_per_trial<=max_events_per_sec)...
                                    &(handles_drgb.drgb.lfpevpair(lfpodNo).no_ref_evs_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNo).no_ref_evs_per_trial<=max_events_per_sec);
                                
                                
                                trials_in_Ev1=(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventERP(event1,:)==1)&these_per_corr&trials_with_event;
                                trials_in_Ev2=(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventERP(event2,:)==1)&these_per_corr&trials_with_event;
                                
                                if (sum(trials_in_Ev1)>=min_trials_per_event)&(sum(trials_in_Ev2)>=min_trials_per_event)
                                    
                                    %Ev1
                                    this_dB_powerEv1=zeros(sum(trials_in_Ev1),length(frequency));
                                    this_dB_powerEv1(:,:)=handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP(trials_in_Ev1,:,shift_ii);
                                    
                                    these_dB_powerEv1=zeros(sum(trials_in_Ev1),length(frequency),length(handles_drgb.drgb.lfpevpair(1).out_times));
                                    these_dB_powerEv1(:,:,:)=handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP(trials_in_Ev1,:,:);
                                    
                                    %Ev2
                                    this_dB_powerEv2=zeros(sum(trials_in_Ev2),length(frequency));
                                    this_dB_powerEv2(:,:)=handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP(trials_in_Ev2,:,shift_ii);
                                    
                                    
                                    these_dB_powerEv2=zeros(sum(trials_in_Ev2),length(frequency),length(handles_drgb.drgb.lfpevpair(1).out_times));
                                    these_dB_powerEv2(:,:,:)=handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP(trials_in_Ev2,:,:);
                                    
                                    for bwii=1:no_bandwidths
                                        
                                        no_ROCs=no_ROCs+1;
                                        this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                        
                                        %Enter Ev1
                                        this_delta_dB_powerEv1=zeros(sum(trials_in_Ev1),1);
                                        this_delta_dB_powerEv1=mean(this_dB_powerEv1(:,this_band),2);
                                        roc_data=[];
                                        roc_data(1:sum(trials_in_Ev1),1)=this_delta_dB_powerEv1;
                                        roc_data(1:sum(trials_in_Ev1),2)=zeros(sum(trials_in_Ev1),1);
                                        
                                        %Enter Ev2
                                        total_trials=sum(trials_in_Ev1)+sum(trials_in_Ev2);
                                        this_delta_dB_powerEv2=zeros(sum(trials_in_Ev2),1);
                                        this_delta_dB_powerEv2=mean(this_dB_powerEv2(:,this_band),2);
                                        roc_data(sum(trials_in_Ev1)+1:total_trials,1)=this_delta_dB_powerEv2;
                                        roc_data(sum(trials_in_Ev1)+1:total_trials,2)=ones(sum(trials_in_Ev2),1);
                                        
                                        
                                        %Find ROC
                                        ROCout(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                        ROCout(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
                                        ROCper_corr_ii(no_ROCs)=per_corr_ii;
                                        ROCoutwin(no_ROCs)=winNo;
                                        ROCbandwidth(no_ROCs)=bwii;
                                        auROC(no_ROCs)=ROCout(no_ROCs).roc.AUC-0.5;
                                        p_valROC(no_ROCs)=ROCout(no_ROCs).roc.p;
                                        
                                        p_vals_ROC=[p_vals_ROC ROCout(no_ROCs).roc.p];
                                        
                                        dB_power_out_Ev1(no_ROCs,1:length(handles_drgb.drgb.lfpevpair(1).out_times))=mean(mean(these_dB_powerEv1(:,this_band,:),2),1);
                                        dB_power_out_Ev2(no_ROCs,1:length(handles_drgb.drgb.lfpevpair(1).out_times))=mean(mean(these_dB_powerEv2(:,this_band,:),2),1);
                                        
                                        
                                    end
                                    
                                else
                                    
                                    if (sum(trials_in_Ev1)<min_trials_per_event)
                                        fprintf(1, ['%d trials in ' evTypeLabels{1} ' percent window #%d fewer than minimum trials per event= %d for file No %d electrode %d\n'],sum(trials_in_Ev1), per_corr_ii, min_trials_per_event,fileNo,elec);
                                    end
                                    
                                    if (sum(trials_in_Ev2)<min_trials_per_event)
                                        fprintf(1, ['%d trials in ' evTypeLabels{2} ' percent window #%d fewer than minimum trials per event= %d for file No %d electrode %d\n'],sum(trials_in_Ev2), per_corr_ii, min_trials_per_event,fileNo,elec);
                                    end
                                    
                                end
                            end
                        else
                            fprintf(1, ['Empty allPower ERP for file No %d electrode %d\n'],files(fileNo),elec);
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
        %Initialize figure counter
        figNo=0;
        
        
        %Initializethe percent significant auROC graph
        try
            close(5)
        catch
        end
        figure(5)
        hold on
        x=0;
        
        for bwii=1:no_bandwidths
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
            
            h2=histogram(auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==1)),edges);
            h2.FaceColor='r';
            h1=histogram(auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==2)),edges);
            h1.FaceColor='b';
            
            xlabel('auROC')
            ylabel('# of electrodes')
            legend(file_label{1},file_label{2})
            title(['auROC for ' freq_names{bwii}])
            xlim([-0.3 0.6])
            ylim([0 45])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Plot the single electrodes
            pos2=[0.8 0.1 0.1 0.8];
            subplot('Position',pos2)
            hold on
            plot(ones(1,sum((ROCbandwidth==bwii)&(ROCper_corr_ii==1))),auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==1)),'o', 'Color',[0.7 0.7 0.7])
            plot(zeros(1,sum((ROCbandwidth==bwii)&(ROCper_corr_ii==2))),auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==2)),'o', 'Color',[0.7 0.7 0.7])
            
            
            %PLot the mean and 95% CI
            plot([0 1],[mean(auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==2))) mean(auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==1)))],'-k','LineWidth', 3)
            CI = bootci(1000, @mean, auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==2)));
            plot([0 0],CI,'-b','LineWidth',3)
            plot(0,mean(auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==2))),'ob','MarkerSize', 10,'MarkerFace','b')
            CI = bootci(1000, @mean, auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==1)));
            plot([1 1],CI,'-r','LineWidth',3)
            plot(1,mean(auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==1))),'or','MarkerSize', 10,'MarkerFace','r')
            ylabel('auROC')
            ylim([-0.2 0.5])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Do the statistics for auROC differences
            a={auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==1)) auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==2))};
            mode_statcond='perm';
            [F df pval_auROCperm] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for permuted anovan for auROC S+ vs S- ' freq_names{bwii} '= %d\n\n'],  pval_auROCperm);
            pvals_auROCperm=[pvals_auROCperm pval_auROCperm];
            
            
            %Plot the bars in the percent significant auROC graph
            figure(5)
            percent_auROCper1=100*sum(p_valROC((ROCbandwidth==bwii)&(ROCper_corr_ii==1))<=pFDRROC)/sum((ROCbandwidth==bwii)&(ROCper_corr_ii==1));
            bar(x,percent_auROCper1,'r')
            
            learn_sig(bwii)=sum(p_valROC((ROCbandwidth==bwii)&(ROCper_corr_ii==1))<=pFDRROC);
            learn_not_sig(bwii)=sum((ROCbandwidth==bwii)&(ROCper_corr_ii==1)-sum(p_valROC((ROCbandwidth==bwii)&(ROCper_corr_ii==1))<=pFDRROC));
            
            percent_auROCper2=100*sum(p_valROC((ROCbandwidth==bwii)&(ROCper_corr_ii==2))<=pFDRROC)/sum((ROCbandwidth==bwii)&(ROCper_corr_ii==2));
            bar(x+1,percent_auROCper2,'b')
            
            prof_sig(bwii)=sum(p_valROC((ROCbandwidth==bwii)&(ROCper_corr_ii==2))<=pFDRROC);
            prof_not_sig(bwii)=sum((ROCbandwidth==bwii)&(ROCper_corr_ii==2)-sum(p_valROC((ROCbandwidth==bwii)&(ROCper_corr_ii==2))<=pFDRROC));
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            x=x+3;
            
            %Plot the ERP for Ev1
            try
                close(figNo+5)
            catch
            end
            figure(figNo+5)
            
            
            out_times=handles_drgb.drgb.lfpevpair(1).out_times';
            
            hold on
            
            %Event 1
            
            %Naive
            mean_Ev1_naive=mean(dB_power_out_Ev1((ROCbandwidth==bwii)&(ROCper_corr_ii==2),:),1)';
            CIEv1naive = bootci(1000, {@mean, dB_power_out_Ev1((ROCbandwidth==bwii)&(ROCper_corr_ii==2),:)})';
            CIEv1naive(:,1)=mean_Ev1_naive-CIEv1naive(:,1);
            CIEv1naive(:,2)=CIEv1naive(:,2)-mean_Ev1_naive;
            %
            %             [hl1, hp1] = boundedline(out_times,mean_Ev1, CI, '--r','transparency',0.05);
            %             outlinebounds(hl1,hp1)
            
            %Proficient
            mean_Ev1_proficient=mean(dB_power_out_Ev1((ROCbandwidth==bwii)&(ROCper_corr_ii==1),:),1)';
            CIEv1prof = bootci(1000, {@mean, dB_power_out_Ev1((ROCbandwidth==bwii)&(ROCper_corr_ii==1),:)})';
            CIEv1prof(:,1)=mean_Ev1_proficient-CIEv1prof(:,1);
            CIEv1prof(:,2)=CIEv1prof(:,2)-mean_Ev1_proficient;
            
            %             [hl1, hp1] = boundedline(out_times,mean_Ev1, CI, 'r','transparency',0.05);
            %             outlinebounds(hl1,hp1);
            
            %Event 2
            
            %Naive
            mean_Ev2_naive=mean(dB_power_out_Ev2((ROCbandwidth==bwii)&(ROCper_corr_ii==2),:),1)';
            CIEv2naive = bootci(1000, {@mean, dB_power_out_Ev2((ROCbandwidth==bwii)&(ROCper_corr_ii==2),:)})';
            CIEv2naive(:,1)=mean_Ev2_naive-CIEv2naive(:,1);
            CIEv2naive(:,2)=CIEv2naive(:,2)-mean_Ev2_naive;
            
            %             [hl1, hp1] = boundedline(out_times,mean_Ev2, CI, '--b','transparency',0.05);
            %             outlinebounds(hl1,hp1)
            
            %Proficient
            mean_Ev2_proficient=mean(dB_power_out_Ev2((ROCbandwidth==bwii)&(ROCper_corr_ii==1),:),1)';
            CIEv2prof = bootci(1000, {@mean, dB_power_out_Ev2((ROCbandwidth==bwii)&(ROCper_corr_ii==1),:)})';
            CIEv2prof(:,1)=mean_Ev2_proficient-CIEv2prof(:,1);
            CIEv2prof(:,2)=CIEv2prof(:,2)-mean_Ev2_proficient;
            
            %             [hl1, hp1] = boundedline(out_times,mean_Ev2, CI, 'b','transparency',0.05);
            %             outlinebounds(hl1,hp1)
            
            
            [hl1, hp1] = boundedline(out_times,mean_Ev1_naive, CIEv1naive, '--r',out_times,mean_Ev1_proficient, CIEv1prof,'r',...
                out_times,mean_Ev2_naive, CIEv2naive, '--b',out_times,mean_Ev2_proficient, CIEv2prof,'b');
            outlinebounds(hl1,hp1)
            
            title(['ERP power in dB for ' freq_names{bwii}])
            %             legend('',[evTypeLabels{1} ' ' file_label{2}],'',[evTypeLabels{1} ' ' file_label{1}],'',[evTypeLabels{2} ' ' file_label{2}],'',[evTypeLabels{2} ' ' file_label{1}])
            xlabel('delta t (sec')
            %             CI = bootci(1000, @mean, dB_power_out_Ev1((ROCbandwidth==bwii)&(ROCper_corr_ii==2),:));
            %             [hl1, hp1] = boundedline(handles_drgb.drgb.lfpevpair(1).out_times',mean(dB_power_out_Ev1((ROCbandwidth==bwii)&(ROCper_corr_ii==2),:),1)', CI', 'b');
            %
            %
            %             CI = bootci(1000, @mean, dB_power_out_Ev1((ROCbandwidth==bwii)&(ROCper_corr_ii==1),:));
            %             [hl1, hp1] = boundedline(handles_drgb.drgb.lfpevpair(1).out_times',mean(dB_power_out_Ev1((ROCbandwidth==bwii)&(ROCper_corr_ii==1),:),1)', CI', 'r');
            %
        end
        
        figure(5)
        title('Percent significant auROC')
        legend(file_label{1},file_label{2})
        ylim([0 100])
        
        pFDRanovanauROC=drsFDRpval(pval_auROCperm);
        fprintf(1, ['\npFDR for premuted anovan p value for difference between ' file_label{1} ' and ' file_label{2} ' for auROC = %d\n\n'],pFDRanovanauROC);
        
        
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix],'learn_sig','learn_not_sig','prof_sig','prof_not_sig');
        
        pffft=1;
        
    case 12
        %Justin's concentration dependence per session
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
        for evNo=1:length(eventType)
            evNo_out(evNo).noWB=0;
        end
        delta_dB_powerEv1WB=[];
        delta_dB_powerEv2WB=[];
        
        
        fprintf(1, ['Pairwise auROC analysis for Fig 1 of Daniel''s paper\n\n'])
        p_vals=[];
        no_files=length(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        no_ROCs=0;
        
        for fileNo=1:no_files
            if sum(files==fileNo)>0
                for elec=1:16
                    if sum(which_electrodes==elec)>0
                        lfpodNo_ref=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                        
                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref)))
                            
                            
                            if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower))
                                szpc=size(percent_windows);
                                for per_ii=1:szpc(1)
                                    percent_mask=[];
                                    percent_mask=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower>=percent_windows(per_ii,1))...
                                        &(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower<=percent_windows(per_ii,2));
                                    
                                    theseEvNos=[];
                                    for evNo=1:length(eventType)
                                        
                                        noWB_for_evNo(evNo)=-1;
                                        
                                        trials_in_event_Ev=[];
                                        trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(eventType(evNo),:)==1)&percent_mask;
                                        
                                        if (sum(trials_in_event_Ev)>=min_trials_per_event)
                                            
                                            lfpodNo=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                            
                                            % Ev1
                                            this_dB_powerref=zeros(sum(trials_in_event_Ev),length(frequency));
                                            this_dB_powerref(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower(trials_in_event_Ev,:));
                                            
                                            
                                            this_dB_power=zeros(sum(trials_in_event_Ev),length(frequency));
                                            this_dB_power(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_event_Ev,:));
                                            
                                            %Wide band spectrum
                                            evNo_out(evNo).noWB=evNo_out(evNo).noWB+1;
                                            evNo_out(evNo).delta_dB_powerEvWB(evNo_out(evNo).noWB,1:length(frequency))=mean(this_dB_power-this_dB_powerref,1);
                                            evNo_out(evNo).per_ii(evNo_out(evNo).noWB)=per_ii;
                                            
                                            noWB_for_evNo(evNo)=evNo_out(evNo).noWB;
                                            
                                            %Do per bandwidth analysis
                                            for bwii=1:no_bandwidths
                                                
                                                this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                                
                                                %Enter the  Ev1
                                                this_delta_dB_powerEv=zeros(sum(trials_in_event_Ev),1);
                                                this_delta_dB_powerEv=mean(this_dB_power(:,this_band)-this_dB_powerref(:,this_band),2);
                                                theseEvNos(evNo,bwii).this_delta_dB_powerEv=this_delta_dB_powerEv;
                                                evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).noWB,bwii)=mean(this_delta_dB_powerEv);
                                                if evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).noWB,bwii)<-20
                                                    fprintf(1, ['mean dB less than -20 dB in file no %d, electrode %d\n'],files(fileNo),elec);
                                                end
                                                
                                            end
                                            
                                            
                                            fprintf(1, ['%d trials in event No %d succesfully processed for file No %d electrode %d\n'],sum(trials_in_event_Ev), evNo,files(fileNo),elec);
                                            
                                        else
                                            
                                            
                                            fprintf(1, ['%d trials in event No %d fewer than minimum trials per event ' evTypeLabels{evNo} ' for file No %d electrode %d\n'],sum(trials_in_event_Ev), min_trials_per_event,files(fileNo),elec);
                                            
                                            
                                        end
                                        
                                        
                                    end
                                    
                                    for evNo1=1:length(eventType)
                                        for evNo2=evNo1+1:length(eventType)
                                            if (noWB_for_evNo(evNo1)~=-1)&(noWB_for_evNo(evNo2)~=-1)
                                                
                                                for bwii=1:no_bandwidths
                                                    
                                                    %Enter Ev1
                                                    trials_in_event_Ev1=length(theseEvNos(evNo1,bwii).this_delta_dB_powerEv);
                                                    this_delta_dB_powerEv1=zeros(trials_in_event_Ev1,1);
                                                    this_delta_dB_powerEv1=theseEvNos(evNo1,bwii).this_delta_dB_powerEv;
                                                    roc_data=[];
                                                    roc_data(1:sum(trials_in_event_Ev1),1)=this_delta_dB_powerEv1;
                                                    roc_data(1:sum(trials_in_event_Ev1),2)=zeros(sum(trials_in_event_Ev1),1);
                                                    
                                                    %Enter Ev2
                                                    trials_in_event_Ev2=length(theseEvNos(evNo2,bwii).this_delta_dB_powerEv);
                                                    total_trials=trials_in_event_Ev1+trials_in_event_Ev2;
                                                    this_delta_dB_powerEv2=zeros(trials_in_event_Ev2,1);
                                                    this_delta_dB_powerEv2=theseEvNos(evNo2,bwii).this_delta_dB_powerEv;
                                                    roc_data(sum(trials_in_event_Ev1)+1:total_trials,1)=this_delta_dB_powerEv2;
                                                    roc_data(sum(trials_in_event_Ev1)+1:total_trials,2)=ones(sum(trials_in_event_Ev2),1);
                                                    
                                                    
                                                    %Find  ROC
                                                    if (trials_in_event_Ev1>=5)&(trials_in_event_Ev2>=5)
                                                        no_ROCs=no_ROCs+1;
                                                        roc=roc_calc(roc_data,0,0.05,0);
                                                        ROCout(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo_ref).fileNo;
                                                        ROCelec(no_ROCs)=elec;
                                                        ROCbandwidth(no_ROCs)=bwii;
                                                        ROCper_ii(no_ROCs)=per_ii;
                                                        ROCEvNo1(no_ROCs)=evNo1;
                                                        ROCEvNo2(no_ROCs)=evNo2;
                                                        auROC(no_ROCs)=roc.AUC-0.5;
                                                        p_valROC(no_ROCs)=roc.p;
                                                        
                                                        p_vals_ROC=[p_vals_ROC roc.p];
                                                    end
                                                end
                                            end
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
        end
        fprintf(1, '\n\n')
        
        
        %Now plot the bounded line for
        
        %Calculate and plot the mean and 95% CI for each event
        figure(1)
        conc_anno_loc = {[0.15 0.15 0.2 0.2], [0.15 0.15 0.2 0.17], [0.15 0.15 0.2 0.14], [0.15 0.15 0.2 0.11], [0.15 0.15 0.2 0.08], [0.15 0.15 0.2 0.05]};
        for evNo=1:length(eventType)
            dB_Ev_ci=zeros(length(frequency),2);
            dB_Ev_mean=[];
            CI=[];
            dB_Ev_mean=mean(evNo_out(evNo).delta_dB_powerEvWB(evNo_out(evNo).per_ii==1,:));
            CI = bootci(1000, {@mean, evNo_out(evNo).delta_dB_powerEvWB(evNo_out(evNo).per_ii==1,:)},'type','cper');
            [hl1, hp1] = boundedline(frequency,dB_Ev_mean', CI', these_colors{evNo});
            annotation('textbox',conc_anno_loc{evNo},'String',evTypeLabels(evNo),'Color',these_colors{evNo},'EdgeColor','none');
            %             rectangle('Position',[0.15 0.15 0.2 0.2],'FaceColor','w','EdgeColor','k');
        end
        
        
        xlabel('Frequency (Hz)')
        ylabel('delta Power (dB)')
        ylim([-20 20]);
        title('Wideband spectrum proficient mice')
        set(gcf,'OuterPosition',[93 36 576 513]);
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        
        %Calculate and plot the mean and 95% CI for each event
        figure(2)
        for evNo=1:length(eventType)
            dB_Ev_ci=zeros(length(frequency),2);
            dB_Ev_mean=[];
            dB_Ev_mean=mean(evNo_out(evNo).delta_dB_powerEvWB(evNo_out(evNo).per_ii==2,:));
            CI = bootci(1000, {@mean, evNo_out(evNo).delta_dB_powerEvWB(evNo_out(evNo).per_ii==2,:)},'type','cper');
            [hl1, hp1] = boundedline(frequency,dB_Ev_mean', CI', these_colors{evNo});
            annotation('textbox',conc_anno_loc{evNo},'String',evTypeLabels(evNo),'Color',these_colors{evNo},'EdgeColor','none');
            rectangle('Position',[0.15 0.15 0.2 0.2],'FaceColor','w','EdgeColor','k');
        end
        
        xlabel('Frequency (Hz)')
        ylabel('delta Power (dB)')
        ylim([-5 10]);
        title('Wideband spectrum naive mice')
        set(gcf,'OuterPosition',[93 550 576 513]);
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        ylim([-20 20])
        
        fig_pos = {[664 550 576 513],[1233 550 576 513],[664 36 576 513],[1233 36 576 513]};
        
        
        %Now plot the histograms and the average
        for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            %Plot the average
            figure(bwii+2)
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            data_dB=[];
            spm=[];
            conc=[];
            
            for evNo=1:length(eventType)
                
                for per_ii=1:2      %performance bins. blue = naive, red = proficient
                    
                    bar_offset=21-evNo*3+(per_ii-1);
                    if per_ii==1
                        bar(bar_offset,mean(evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).per_ii==per_ii,bwii)),'r','LineWidth', 3)
                    else
                        bar(bar_offset,mean(evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).per_ii==per_ii,bwii)),'b','LineWidth', 3)
                    end
                    plot(bar_offset,mean(evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).per_ii==per_ii,bwii)),'ok','LineWidth', 3)
                    CI = bootci(1000, {@mean, evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).per_ii==per_ii,bwii)},'type','cper');
                    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                    plot((bar_offset)*ones(1,sum(evNo_out(evNo).per_ii==per_ii)),evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).per_ii==per_ii,bwii),'o',...
                        'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                    data_dB=[data_dB evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).per_ii==per_ii,bwii)'];
                    switch evNo
                        case {1,2,3}
                            spm=[spm zeros(1,evNo_out(evNo).noWB)];
                            
                        case {4,5,6}
                            spm=[spm ones(1,evNo_out(evNo).noWB)];
                    end
                    conc=[conc evNo*ones(1,evNo_out(evNo).noWB)];
                    
                    annotation('textbox',conc_anno_loc{per_ii},'String',per_lab(per_ii),'Color',these_colors{per_ii},'EdgeColor','none');
                end
            end
            title(freq_names{bwii})
            set(gcf,'OuterPosition',fig_pos{bwii});
            bar_lab_loc = [3.5 6.5 9.5 12.5 15.5 18.5];
            xticks(bar_lab_loc)
            xticklabels(concs2)
            xlabel('Concentration (%)')
            ylabel('Delta power (dB)')
            %             p=anovan(data_dB,{spm});
        end
        
        
        pFDRauROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for auROC  = %d\n\n'],pFDRauROC);
        %Plot cumulative histos for auROCs
        
        
        
        %Plot percent significant ROC
        
        figNo=0;
        for bwii=1:no_bandwidths
            for pcii=1:szpc(1)
                no_pairs=0;
                EvNo1=[];
                EvNo2=[];
                per_sig=zeros(length(eventType),length(eventType));
                for evNo1=1:length(eventType)
                    for evNo2=evNo1+1:length(eventType)
                        
                        
                        no_pairs=no_pairs+1;
                        these_ROCs=(ROCbandwidth==bwii)&(ROCEvNo1==evNo1)&(ROCEvNo2==evNo2)&(ROCper_ii==pcii);
                        sig(no_pairs)=sum((p_valROC<=pFDRauROC)&these_ROCs);
                        not_sig(no_pairs)=sum(these_ROCs)-sum((p_valROC<=pFDRauROC)&these_ROCs);
                        EvNo1(no_pairs)=evNo1;
                        EvNo2(no_pairs)=evNo2;
                        per_sig(evNo1,evNo2)=100*sig(no_pairs)/(sig(no_pairs)+not_sig(no_pairs));
                        
                    end
                end
                
                for evNo1=1:length(eventType)
                    for evNo2=evNo1+1:length(eventType)
                        
                        if per_sig(evNo1,evNo2)==0
                            per_sig(evNo1,evNo2)=100/64;
                        end
                    end
                end
                
                %Plot the pseudocolor for percent significant auROCs
                figNo = get(gcf,'Number')+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                
                evNos_for_1=1:length(eventType);
                evNos_for_2=[1:length(eventType)]';
                
                drg_pcolor(repmat(evNos_for_1,length(eventType),1),repmat(evNos_for_2,1,length(eventType)),per_sig)
                cmjet=colormap(jet);
                cmjet(1,1)=0.7;
                cmjet(1,2)=0.7;
                cmjet(1,3)=0.7;
                colormap(cmjet)
                
                hold on
                plot([4 4],[1 4],'-w','LineWidth', 5)
                plot([4 7],[4 4],'-w','LineWidth', 5)
                
                ax=gca;
                set(ax,'XTickLabel','')
                ylabel('dB')
                xticks([1.5:1:length(eventType)+1])
                xticklabels(handles_pars.concs2)
                yticks([1.5:1:length(eventType)+1])
                yticklabels(handles_pars.concs2)
                
                title(['Percent auROC significantly different from zero ' freq_names{bwii} ' ' per_lab(pcii)])
                set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                
                pfft=1
            end
        end
        
        figNo = get(gcf,'Number')+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        
        
        set(hFig, 'units','normalized','position',[.83 .1 .05 .3])
        
        prain=[0:100/99:100];
        drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
        colormap jet
        shading interp
        ax=gca;
        set(ax,'XTickLabel','')
        
        ppno=0;
        for pair_no1=1:no_pairs
            for pair_no2=pair_no1+1:no_pairs
                ppno=ppno+1;
                [pChiSq(ppno), Q]= chi2test([sig(pair_no1), not_sig(pair_no1); sig(pair_no2), not_sig(pair_no2)]);
                fprintf(1, ['pchi  = %d\n'],pChiSq(ppno));
            end
        end
        
        pFDRchisq=drsFDRpval(pChiSq);
        fprintf(1, ['pFDR for ChiSquared  = %d\n\n'],pFDRchisq);
        
        pfft=1
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix]);
        
    case 13
        %Compare auROC for ERP wavelet LFP power between two percent correct
        %windows at different ERP delta t values
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
        
        
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        sz_per=size(percent_windows);
        no_per_win=sz_per(1);
        
        
        
        fprintf(1, ['Pairwise auROC log power LFP ERP analysis for ' evTypeLabels{1} ' and ' evTypeLabels{2} ' LFP power\n\n'])
        p_vals=[];
        
        no_files=length(files);
        
        p_vals_ROC=[];
        no_sii=0;
        out_times=handles_drgb.drgb.lfpevpair(1).out_times';
        for shift_ii=1:delta_t_ii:length(handles_drgb.drgb.lfpevpair(1).out_times)
            no_sii=no_sii+1;
            no_pvals(no_sii)=0;
            sii_times(no_sii)=out_times(shift_ii);
        end
        
        for fileNo=1:no_files
            
            
            for elec=1:16
                if sum(which_electrodes==elec)>0
                    lfpodNo=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                    
                    if ~isempty(handles_drgb.drgb.lfpevpair(lfpodNo))
                        
                        if ~isempty(handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP)
                            
                            
                            for per_corr_ii=1:no_per_win
                                
                                these_per_corr=(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrERP>=percent_windows(per_corr_ii,1))&...
                                    (handles_drgb.drgb.lfpevpair(lfpodNo).perCorrERP<=percent_windows(per_corr_ii,2));
                                
                                %Which trials have events (i.e. licks)?
                                trials_with_event=(handles_drgb.drgb.lfpevpair(lfpodNo).no_events_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNo).no_events_per_trial<=max_events_per_sec)...
                                    &(handles_drgb.drgb.lfpevpair(lfpodNo).no_ref_evs_per_trial>0)&(handles_drgb.drgb.lfpevpair(lfpodNo).no_ref_evs_per_trial<=max_events_per_sec);
                                
                                
                                trials_in_Ev1=(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventERP(event1,:)==1)&these_per_corr&trials_with_event;
                                trials_in_Ev2=(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventERP(event2,:)==1)&these_per_corr&trials_with_event;
                                
                                if (sum(trials_in_Ev1)>=min_trials_per_event)&(sum(trials_in_Ev2)>=min_trials_per_event)
                                    
                                    for bwii=1:no_bandwidths
                                        
                                        no_ROCs=no_ROCs+1;
                                        this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                        
                                        these_dB_powerEv1=zeros(sum(trials_in_Ev1),length(frequency),length(handles_drgb.drgb.lfpevpair(1).out_times));
                                        these_dB_powerEv1(:,:,:)=handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP(trials_in_Ev1,:,:);
                                        
                                        these_dB_powerEv2=zeros(sum(trials_in_Ev2),length(frequency),length(handles_drgb.drgb.lfpevpair(1).out_times));
                                        these_dB_powerEv2(:,:,:)=handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP(trials_in_Ev2,:,:);
                                        
                                        dB_power_out_Ev1(no_ROCs,1:length(handles_drgb.drgb.lfpevpair(1).out_times))=mean(mean(these_dB_powerEv1(:,this_band,:),2),1);
                                        dB_power_out_Ev2(no_ROCs,1:length(handles_drgb.drgb.lfpevpair(1).out_times))=mean(mean(these_dB_powerEv2(:,this_band,:),2),1);
                                        
                                        no_dt_ii=0;
                                        for shift_ii=1:delta_t_ii:length(handles_drgb.drgb.lfpevpair(1).out_times)
                                            %Ev1
                                            this_dB_powerEv1=zeros(sum(trials_in_Ev1),length(frequency));
                                            this_dB_powerEv1(:,:)=handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP(trials_in_Ev1,:,shift_ii);
                                            
                                            %Ev2
                                            this_dB_powerEv2=zeros(sum(trials_in_Ev2),length(frequency));
                                            this_dB_powerEv2(:,:)=handles_drgb.drgb.lfpevpair(lfpodNo).log_P_tERP(trials_in_Ev2,:,shift_ii);
                                            
                                            %Enter roc data for Ev1
                                            this_delta_dB_powerEv1=zeros(sum(trials_in_Ev1),1);
                                            this_delta_dB_powerEv1=mean(this_dB_powerEv1(:,this_band),2);
                                            roc_data=[];
                                            roc_data(1:sum(trials_in_Ev1),1)=this_delta_dB_powerEv1;
                                            roc_data(1:sum(trials_in_Ev1),2)=zeros(sum(trials_in_Ev1),1);
                                            
                                            %Enter roc data for Ev2
                                            total_trials=sum(trials_in_Ev1)+sum(trials_in_Ev2);
                                            this_delta_dB_powerEv2=zeros(sum(trials_in_Ev2),1);
                                            this_delta_dB_powerEv2=mean(this_dB_powerEv2(:,this_band),2);
                                            roc_data(sum(trials_in_Ev1)+1:total_trials,1)=this_delta_dB_powerEv2;
                                            roc_data(sum(trials_in_Ev1)+1:total_trials,2)=ones(sum(trials_in_Ev2),1);
                                            
                                            
                                            %Find ROC
                                            roc=roc_calc(roc_data,0,0.05,0);
                                            ROCout(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
                                            ROCper_corr_ii(no_ROCs)=per_corr_ii;
                                            ROCoutwin(no_ROCs)=winNo;
                                            ROCbandwidth(no_ROCs)=bwii;
                                            no_dt_ii=no_dt_ii+1;
                                            auROC(no_ROCs,no_dt_ii)=roc.AUC-0.5;
                                            p_valROC(no_ROCs,no_dt_ii)=roc.p;
                                            
                                            no_pvals(no_dt_ii)=no_pvals(no_dt_ii)+1;
                                            p_vals_ROC(no_dt_ii,no_pvals(no_dt_ii))=roc.p;
                                            
                                        end
                                        
                                    end
                                else
                                    
                                    if (sum(trials_in_Ev1)<min_trials_per_event)
                                        fprintf(1, ['%d trials in ' evTypeLabels{1} ' percent window #%d fewer than minimum trials per event= %d for file No %d electrode %d\n'],sum(trials_in_Ev1), per_corr_ii, min_trials_per_event,fileNo,elec);
                                    end
                                    
                                    if (sum(trials_in_Ev2)<min_trials_per_event)
                                        fprintf(1, ['%d trials in ' evTypeLabels{2} ' percent window #%d fewer than minimum trials per event= %d for file No %d electrode %d\n'],sum(trials_in_Ev2), per_corr_ii, min_trials_per_event,fileNo,elec);
                                    end
                                    
                                end
                            end
                        else
                            fprintf(1, ['Empty allPower ERP for file No %d electrode %d\n'],files(fileNo),elec);
                        end
                        
                    else
                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],files(fileNo),elec);
                    end
                end
            end
        end
        fprintf(1, '\n\n')
        
        
        for sii=1:no_sii
            pFDRROC(sii)=drsFDRpval(p_vals_ROC(sii,1:no_pvals(sii)));
            fprintf(1, ['pFDR for significant difference of auROC p value from 0.5  = %d\n'],pFDRROC(sii));
        end
        
        fprintf(1, '\n\n')
        
        
        %Plot cumulative histos for auROCs
        %Initialize figure counter
        figNo=0;
        
        
        for bwii=1:no_bandwidths
            
            %Plot auROC for proficient and naive
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            
            %Proficient
            mean_auROC_proficient=mean(auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==1),:),1)';
            CIauROCprof = bootci(1000, {@mean, auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==1),:)})';
            CIauROCprof(:,1)=mean_auROC_proficient-CIauROCprof(:,1);
            CIauROCprof(:,2)=CIauROCprof(:,2)-mean_auROC_proficient;
            
            %Naive
            mean_auROC_naive=mean(auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==2),:),1)';
            CIauROCnaive = bootci(1000, {@mean, auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==2),:)})';
            CIauROCnaive(:,1)=mean_auROC_naive-CIauROCnaive(:,1);
            CIauROCnaive(:,2)=CIauROCnaive(:,2)-mean_auROC_naive;
            
            [hl1, hp1] = boundedline(sii_times,mean_auROC_proficient, CIauROCprof, 'r',sii_times,mean_auROC_naive, CIauROCnaive,'b');
            outlinebounds(hl1,hp1)
            
            xlabel('delta t (sec)')
            ylabel('auROC')
            legend(file_label{1}, file_label{2})
            title(['auROC for ' freq_names{bwii}])
            xlim([-0.2 0.2])
            ylim([0 0.5])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            %Plot the percent significant auROCs
            
            try
                close(figNo+4)
            catch
            end
            figure(figNo+4)
            hold on
            
            for sii=1:no_sii
                sig_per_proficient(bwii,sii)=100*sum(p_vals_ROC(sii,(ROCbandwidth==bwii)&(ROCper_corr_ii==1))<=pFDRROC(sii))/sum((ROCbandwidth==bwii)&(ROCper_corr_ii==1));
                sig_per_naive(bwii,sii)=100*sum(p_vals_ROC(sii,(ROCbandwidth==bwii)&(ROCper_corr_ii==2))<=pFDRROC(sii))/sum((ROCbandwidth==bwii)&(ROCper_corr_ii==2));
                
                no_sig_proficient=sum(p_vals_ROC(sii,(ROCbandwidth==bwii)&(ROCper_corr_ii==1))<=pFDRROC(sii));
                no_not_sig_proficient=sum(p_vals_ROC(sii,(ROCbandwidth==bwii)&(ROCper_corr_ii==1))>pFDRROC(sii));
                
                no_sig_naive=sum(p_vals_ROC(sii,(ROCbandwidth==bwii)&(ROCper_corr_ii==2))<=pFDRROC(sii));
                no_not_sig_naive=sum(p_vals_ROC(sii,(ROCbandwidth==bwii)&(ROCper_corr_ii==2))>pFDRROC(sii));
                
                [p(bwii,sii), Q]= chi2test([no_sig_proficient, no_not_sig_proficient; no_sig_naive, no_not_sig_naive]);
            end
            
            plot(sii_times,sig_per_proficient(bwii,:),'-or')
            plot(sii_times,sig_per_naive(bwii,:),'-ob')
            
            xlabel('delta t (sec)')
            ylabel('Percent significant auROC')
            legend(file_label{1}, file_label{2})
            title(['Percent significant auROC for ' freq_names{bwii}])
            xlim([-0.2 0.2])
            ylim([0 100])
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            
            
            
            %             %Do the statistics for auROC differences
            %             a={auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==1)) auROC((ROCbandwidth==bwii)&(ROCper_corr_ii==2))};
            %             mode_statcond='perm';
            %             [F df pval_auROCperm] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            %             fprintf(1, ['p value for permuted anovan for auROC S+ vs S- ' freq_names{bwii} '= %d\n\n'],  pval_auROCperm);
            %             pvals_auROCperm=[pvals_auROCperm pval_auROCperm];
            
            
            
            
            %Plot the ERP for Ev1
            try
                close(figNo+8)
            catch
            end
            figure(figNo+8)
            
            
            out_times=handles_drgb.drgb.lfpevpair(1).out_times';
            
            hold on
            
            %Event 1
            
            %Naive
            mean_Ev1_naive=mean(dB_power_out_Ev1((ROCbandwidth==bwii)&(ROCper_corr_ii==2),:),1)';
            CIEv1naive = bootci(1000, {@mean, dB_power_out_Ev1((ROCbandwidth==bwii)&(ROCper_corr_ii==2),:)})';
            CIEv1naive(:,1)=mean_Ev1_naive-CIEv1naive(:,1);
            CIEv1naive(:,2)=CIEv1naive(:,2)-mean_Ev1_naive;
            %
            %             [hl1, hp1] = boundedline(out_times,mean_Ev1, CI, '--r','transparency',0.05);
            %             outlinebounds(hl1,hp1)
            
            %Proficient
            mean_Ev1_proficient=mean(dB_power_out_Ev1((ROCbandwidth==bwii)&(ROCper_corr_ii==1),:),1)';
            CIEv1prof = bootci(1000, {@mean, dB_power_out_Ev1((ROCbandwidth==bwii)&(ROCper_corr_ii==1),:)})';
            CIEv1prof(:,1)=mean_Ev1_proficient-CIEv1prof(:,1);
            CIEv1prof(:,2)=CIEv1prof(:,2)-mean_Ev1_proficient;
            
            %             [hl1, hp1] = boundedline(out_times,mean_Ev1, CI, 'r','transparency',0.05);
            %             outlinebounds(hl1,hp1);
            
            %Event 2
            
            %Naive
            mean_Ev2_naive=mean(dB_power_out_Ev2((ROCbandwidth==bwii)&(ROCper_corr_ii==2),:),1)';
            CIEv2naive = bootci(1000, {@mean, dB_power_out_Ev2((ROCbandwidth==bwii)&(ROCper_corr_ii==2),:)})';
            CIEv2naive(:,1)=mean_Ev2_naive-CIEv2naive(:,1);
            CIEv2naive(:,2)=CIEv2naive(:,2)-mean_Ev2_naive;
            
            %             [hl1, hp1] = boundedline(out_times,mean_Ev2, CI, '--b','transparency',0.05);
            %             outlinebounds(hl1,hp1)
            
            %Proficient
            mean_Ev2_proficient=mean(dB_power_out_Ev2((ROCbandwidth==bwii)&(ROCper_corr_ii==1),:),1)';
            CIEv2prof = bootci(1000, {@mean, dB_power_out_Ev2((ROCbandwidth==bwii)&(ROCper_corr_ii==1),:)})';
            CIEv2prof(:,1)=mean_Ev2_proficient-CIEv2prof(:,1);
            CIEv2prof(:,2)=CIEv2prof(:,2)-mean_Ev2_proficient;
            
            %             [hl1, hp1] = boundedline(out_times,mean_Ev2, CI, 'b','transparency',0.05);
            %             outlinebounds(hl1,hp1)
            
            
            [hl1, hp1] = boundedline(out_times,mean_Ev1_naive, CIEv1naive, '--r',out_times,mean_Ev1_proficient, CIEv1prof,'r',...
                out_times,mean_Ev2_naive, CIEv2naive, '--b',out_times,mean_Ev2_proficient, CIEv2prof,'b');
            outlinebounds(hl1,hp1)
            
            title(['ERP power in dB for ' freq_names{bwii}])
            %             legend('',[evTypeLabels{1} ' ' file_label{2}],'',[evTypeLabels{1} ' ' file_label{1}],'',[evTypeLabels{2} ' ' file_label{2}],'',[evTypeLabels{2} ' ' file_label{1}])
            xlabel('delta t (sec')
            
        end
        
        for sii=1:no_sii
            pFDRROC_chisq(sii)=drsFDRpval(p(:,sii));
            fprintf(1, ['pFDR for chi squared  = %d\n'],pFDRROC_chisq(sii));
        end
        
        figNo=4;
        for bwii=1:no_bandwidths
            figNo=figNo+1;
            figure(figNo)
            hold on
            
            plot(sii_times(p(bwii,:)<=pFDRROC_chisq),(sig_per_proficient(bwii,p(bwii,:)<=pFDRROC_chisq)+sig_per_naive(bwii,p(bwii,:)<=pFDRROC_chisq))/2,'*k')
        end
        
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix]);
        
        
        
    case 14
        %Justin's per mouse auROC analysis for LFP power
        %For the proficient mice in the first and last sessions
        %plot the LFP spectrum for S+ vs S-, plot LFP power for S+ vs S- for each electrode and plot auROCs
        %NOTE: This does the analysis in all the files and DOES not distinguish between groups!!!
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
                        
                        %                         theseEvNos_thisMouse_thisElec=[];
                        %                         for evNo=1:length(eventType)
                        %                             for bwii=1:no_bandwidths
                        %                                 theseEvNos_thisMouse_thisElec(evNo,bwii).noEv=0;
                        %                             end
                        %                         end
                        
                        theseEvNos=[];
                        for evNo=1:length(eventType)
                            for bwii=1:no_bandwidths
                                theseEvNos(evNo,bwii).noEv=0;
                            end
                        end
                        
                        mouse_has_files=0;
                        for fileNo=1:no_files
                            if sum(files==fileNo)>0
                                if handles_drgb.drgbchoices.mouse_no(fileNo)==mouseNo
                                    
                                    
                                    
                                    lfpodNo_ref=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                                    
                                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref)))
                                        
                                        
                                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower))
                                            
                                            percent_mask=[];
                                            percent_mask=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower>=percent_windows(per_ii,1))...
                                                &(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower<=percent_windows(per_ii,2));
                                            
                                            
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
                                                    evNo_out(evNo).noWB=evNo_out(evNo).noWB+1;
                                                    evNo_out(evNo).delta_dB_powerEvWB(evNo_out(evNo).noWB,1:length(frequency))=mean(this_dB_power-this_dB_powerref,1);
                                                    evNo_out(evNo).per_ii(evNo_out(evNo).noWB)=per_ii;
                                                    
                                                    noWB_for_evNo(evNo)=evNo_out(evNo).noWB;
                                                    
                                                    %Do per bandwidth analysis
                                                    for bwii=1:no_bandwidths
                                                        
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
                                                        theseEvNosPerEl(evNo,bwii,elec).noEv=theseEvNosPerEl(evNo,bwii,elec).noEv+sum(trials_in_event_Ev);
                                                        
                                                        %                                                         theseEvNos_thisMouse_thisElec(evNo,bwii).this_delta_dB_Ev(theseEvNos_thisMouse_thisElec(evNo,bwii).noEv+1:theseEvNos_thisMouse_thisElec(evNo,bwii).noEv+sum(trials_in_event_Ev))=this_delta_dB_powerEv;
                                                        %                                                         theseEvNos_thisMouse_thisElec(evNo,bwii).whichMouse(theseEvNos_thisMouse_thisElec(evNo,bwii).noEv+1:theseEvNos_thisMouse_thisElec(evNo,bwii).noEv+sum(trials_in_event_Ev))=mouseNo*ones(1,length(this_delta_dB_powerEv));
                                                        %
                                                        
                                                        evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).noWB,bwii)=mean(this_delta_dB_powerEv);
                                                        %                                                     if evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).noWB,bwii)<-20
                                                        %                                                         fprintf(1, ['mean dB less than -20 dB in file no %d, electrode %d\n'],files(fileNo),elec);
                                                        %                                                     end
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
                                        
                                        
                                    else
                                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],fileNo,elec);
                                        
                                        
                                    end
                                end %if mouseNo
                            end
                        end %fileNo
                        
                        if mouse_has_files==1
                            mouse_included(mouseNo)=1;
                            %Calculate per mouse per electrode delta_dB
                            for evNo=1:length(eventType)
                                for bwii=1:no_bandwidths
                                    delta_dB_No_per_mouse=delta_dB_No_per_mouse+1;
                                    %                                     this_mouse_delta_dB=[];
                                    %                                     this_mouse_delta_dB=theseEvNos_thisMouse_thisElec(evNo,bwii).this_delta_dB_Ev(theseEvNos_thisMouse_thisElec(evNo,bwii).whichMouse==mouseNo);
                                    %                                     delta_dB_per_mouse(delta_dB_No_per_mouse)=mean(this_mouse_delta_dB);
                                    delta_dB_per_mouse(delta_dB_No_per_mouse)=mean(theseEvNos(evNo,bwii).this_delta_dB_powerEv(1:theseEvNos(evNo,bwii).noEv));
                                    delta_dB_perii_per_mouse(delta_dB_No_per_mouse)=per_ii;
                                    delta_dB_evNo_per_mouse(delta_dB_No_per_mouse)=evNo;
                                    delta_dB_bwii_per_mouse(delta_dB_No_per_mouse)=bwii;
                                    delta_dB_mouseNo_per_mouse(delta_dB_No_per_mouse)=mouseNo;
                                    delta_dB_electrode_per_mouse(delta_dB_No_per_mouse)=elec;
                                end
                            end
                            
                            %Calculate per electrode ROC
                            can_calculate_auroc=1;
                            if can_calculate_auroc==1
                                for evNo1=1:length(eventType)
                                    for evNo2=evNo1+1:length(eventType)
                                        for bwii=1:no_bandwidths
                                            
                                            
                                            %Enter Ev1
                                            trials_in_event_Ev1=length(theseEvNos(evNo1,bwii).this_delta_dB_powerEv);
                                            this_delta_dB_powerEv1=zeros(trials_in_event_Ev1,1);
                                            this_delta_dB_powerEv1=theseEvNos(evNo1,bwii).this_delta_dB_powerEv;
                                            roc_data=[];
                                            roc_data(1:trials_in_event_Ev1,1)=this_delta_dB_powerEv1;
                                            roc_data(1:trials_in_event_Ev1,2)=zeros(trials_in_event_Ev1,1);
                                            
                                            %Enter Ev2
                                            trials_in_event_Ev2=length(theseEvNos(evNo2,bwii).this_delta_dB_powerEv);
                                            total_trials=trials_in_event_Ev1+trials_in_event_Ev2;
                                            this_delta_dB_powerEv2=zeros(trials_in_event_Ev2,1);
                                            this_delta_dB_powerEv2=theseEvNos(evNo2,bwii).this_delta_dB_powerEv;
                                            roc_data(trials_in_event_Ev1+1:total_trials,1)=this_delta_dB_powerEv2;
                                            roc_data(trials_in_event_Ev1+1:total_trials,2)=ones(trials_in_event_Ev2,1);
                                            
                                            
                                            %Find  per electrode ROC
                                            if (trials_in_event_Ev1>=5)&(trials_in_event_Ev2>=5)
                                                no_ROCs=no_ROCs+1;
                                                roc=roc_calc(roc_data,0,0.05,0);
                                                ROCout(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo_ref).fileNo;
                                                ROCelec(no_ROCs)=elec;
                                                ROCbandwidth(no_ROCs)=bwii;
                                                ROCper_ii(no_ROCs)=per_ii;
                                                ROCEvNo1(no_ROCs)=evNo1;
                                                ROCEvNo2(no_ROCs)=evNo2;
                                                if ((abs(evNo1-2)<=1)&(abs(evNo2-5)<=1))||((abs(evNo1-5)<=1)&(abs(evNo2-2)<=1))
                                                    ROC_between(no_ROCs)=1;
                                                else
                                                    ROC_between(no_ROCs)=0;
                                                end
                                                ROC_neighbor(no_ROCs)=abs(evNo1-evNo2);
                                                auROC(no_ROCs)=roc.AUC-0.5;
                                                p_valROC(no_ROCs)=roc.p;
                                                p_vals_ROC=[p_vals_ROC roc.p];
                                                
                                                %                                                 if (per_ii==1)&(bwii==4)&(roc.AUC-0.5>0.3)
                                                %                                                     %This is here to stop and plot the ROC
                                                %                                                     %roc_out=roc_calc(roc_data);
                                                %                                                     pffft=1;
                                                %                                                 end
                                                
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
                            
                        else
                            mouse_included(mouseNo)=0;
                        end
                        
                        
                    end
                    
                end
                
                %Calculate per mouse ROC
                if mouse_has_files==1
                    for evNo1=1:length(eventType)
                        for evNo2=evNo1+1:length(eventType)
                            
                            
                            for bwii=1:no_bandwidths
                                
                                %Enter Ev1
                                trials_in_event_Ev1=length(theseEvNosPerEl(evNo1,bwii,which_electrodes(1)).this_delta_dB_powerEv);
                                this_delta_dB_powerEv1=zeros(trials_in_event_Ev1,1);
                                for elec=which_electrodes
                                    this_delta_dB_powerEv1=this_delta_dB_powerEv1+(theseEvNosPerEl(evNo1,bwii,elec).this_delta_dB_powerEv')/length(which_electrodes);
                                end
                                
                                roc_data=[];
                                roc_data(1:trials_in_event_Ev1,1)=this_delta_dB_powerEv1;
                                roc_data(1:trials_in_event_Ev1,2)=zeros(trials_in_event_Ev1,1);
                                
                                %Enter Ev2
                                trials_in_event_Ev2=length(theseEvNosPerEl(evNo2,bwii,which_electrodes(1)).this_delta_dB_powerEv);
                                total_trials=trials_in_event_Ev1+trials_in_event_Ev2;
                                this_delta_dB_powerEv2=zeros(trials_in_event_Ev2,1);
                                for elec=which_electrodes
                                    this_delta_dB_powerEv2=this_delta_dB_powerEv2+(theseEvNosPerEl(evNo2,bwii,elec).this_delta_dB_powerEv')/length(which_electrodes);
                                end
                                
                                roc_data(trials_in_event_Ev1+1:total_trials,1)=this_delta_dB_powerEv2;
                                roc_data(trials_in_event_Ev1+1:total_trials,2)=ones(trials_in_event_Ev2,1);
                                
                                
                                %Find  per electrode ROC
                                if (trials_in_event_Ev1>=5)&(trials_in_event_Ev2>=5)
                                    per_mouse_no_ROCs=per_mouse_no_ROCs+1;
                                    roc=roc_calc(roc_data,0,0.05,0);
                                    per_mouse_ROCout(per_mouse_no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo_ref).fileNo;
                                    per_mouse_ROCbandwidth(per_mouse_no_ROCs)=bwii;
                                    per_mouse_ROCper_ii(per_mouse_no_ROCs)=per_ii;
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
            
            
            
        end
        fprintf(1, '\n\n')
        
        pFDRauROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for per electrode auROC  = %d\n\n'],pFDRauROC);
        
        per_mouse_pFDRauROC=drsFDRpval(per_mouse_p_vals_ROC);
        fprintf(1, ['pFDR for per mouse auROC  = %d\n\n'],per_mouse_pFDRauROC);
        
        
        
        %Now plot the delta dB power per session per electrode
        conc_anno_loc = {[0.15 0.15 0.2 0.2], [0.15 0.15 0.2 0.17], [0.15 0.15 0.2 0.14], [0.15 0.15 0.2 0.11], [0.15 0.15 0.2 0.08], [0.15 0.15 0.2 0.05]};
        fig_pos = {[664 550 576 513],[1233 550 576 513],[664 36 576 513],[1233 36 576 513]};
        for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            %Plot the average
            
            try
                close(bwii)
            catch
            end
            figure(bwii)
            
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            data_dB=[];
            spm=[];
            conc=[];
            
            
            for evNo=1:length(eventType)
                
                for per_ii=1:2      %performance bins. blue = naive, red = proficient
                    
                    bar_offset=21-evNo*3+(2-per_ii);
                    if per_ii==1
                        bar(bar_offset,mean(evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).per_ii==per_ii,bwii)),'r','LineWidth', 3)
                    else
                        bar(bar_offset,mean(evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).per_ii==per_ii,bwii)),'b','LineWidth', 3)
                    end
                    plot(bar_offset,mean(evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).per_ii==per_ii,bwii)),'ok','LineWidth', 3)
                    CI = bootci(1000, {@mean, evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).per_ii==per_ii,bwii)},'type','cper');
                    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                    plot((bar_offset)*ones(1,sum(evNo_out(evNo).per_ii==per_ii)),evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).per_ii==per_ii,bwii),'o',...
                        'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                    data_dB=[data_dB evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).per_ii==per_ii,bwii)'];
                    switch evNo
                        case {1,2,3}
                            spm=[spm zeros(1,evNo_out(evNo).noWB)];
                            
                        case {4,5,6}
                            spm=[spm ones(1,evNo_out(evNo).noWB)];
                    end
                    conc=[conc evNo*ones(1,evNo_out(evNo).noWB)];
                    
                    annotation('textbox',conc_anno_loc{per_ii},'String',per_lab(per_ii),'Color',these_colors{3-per_ii},'EdgeColor','none');
                end
            end
            title([freq_names{bwii} ' average delta dB power per session, per electrode'])
            set(gcf,'OuterPosition',fig_pos{bwii});
            
            if sum(eventType==3)>0
                xticks([15,16,18,19])
                xticklabels({evTypeLabels{2},evTypeLabels{2},evTypeLabels{1},evTypeLabels{1}})
            else
                bar_lab_loc = [3.5 6.5 9.5 12.5 15.5 18.5];
                xticks(bar_lab_loc)
                xticklabels(concs2)
                xlabel('Concentration (%)')
            end
            
            
            ylabel('Delta power (dB)')
            %             p=anovan(data_dB,{spm});
        end
        
        %Now plot the histograms and the average per mouse LFP power computed per electrode
        conc_anno_loc = {[0.15 0.15 0.2 0.2], [0.15 0.15 0.2 0.17], [0.15 0.15 0.2 0.14], [0.15 0.15 0.2 0.11], [0.15 0.15 0.2 0.08], [0.15 0.15 0.2 0.05]};
        fig_pos = {[664 550 576 513],[1233 550 576 513],[664 36 576 513],[1233 36 576 513]};
        for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            %Plot the average
            
            try
                close(bwii+4)
            catch
            end
            figure(bwii+4)
            
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            data_delta_dB=[];
            prof_naive=[];
            events=[];
            mice=[];
            electrodes=[];
            
            for evNo=1:length(eventType)
                
                for per_ii=1:2      %performance bins. blue = naive, red = proficient
                    
                    bar_offset=21-evNo*3+(2-per_ii);
                    if per_ii==1
                        bar(bar_offset,mean(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii))),'r','LineWidth', 3)
                    else
                        bar(bar_offset,mean(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii))),'b','LineWidth', 3)
                    end
                    
                    %Individual points; in the future add lines linking the
                    %points?
                    plot((bar_offset)*ones(1,sum((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii))),...
                        delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)),'o',...
                        'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                    
                    %Average and CI
                    plot(bar_offset,mean(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii))),'ok','LineWidth', 3)
                    CI = bootci(1000, {@mean, delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii))},'type','cper');
                    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                    
                    annotation('textbox',conc_anno_loc{per_ii},'String',per_lab(per_ii),'Color',these_colors{3-per_ii},'EdgeColor','none');
                    
                    %Save data for anovan
                    data_delta_dB=[data_delta_dB delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii))];
                    prof_naive=[prof_naive per_ii*ones(1,sum((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)))];
                    events=[events evNo*ones(1,sum((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)))];
                    mice=[mice delta_dB_mouseNo_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii))];
                    electrodes=[electrodes delta_dB_electrode_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii))];
                end
            end
            title([freq_names{bwii} ' average delta dB power per mouse, per electrode'])
            set(gcf,'OuterPosition',fig_pos{bwii});
            
            if sum(eventType==3)>0
                xticks([15,16,18,19])
                xticklabels({evTypeLabels{2},evTypeLabels{2},evTypeLabels{1},evTypeLabels{1}})
            else
                bar_lab_loc = [3.5 6.5 9.5 12.5 15.5 18.5];
                xticks(bar_lab_loc)
                xticklabels(concs2)
                xlabel('Concentration (%)')
            end
            
            ylabel('Delta power (dB)')
            
            
            
            
            
            %Calculate anovan for inteaction
            [p,tbl,stats]=anovan(data_delta_dB,{prof_naive events mice electrodes},'varnames',{'proficient_vs_naive','events','mice','electrodes'},'display','off','random',[3 4]);
            fprintf(1, ['p value for anovan delta dB histogram per mouse per electrode for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan delta dB histogram  per mouse per electrode for events ' freq_names{bwii} '= %d \n'],  p(2));
            fprintf(1, ['p value for anovan delta dB histogram  per mouse per electrode for mice for ' freq_names{bwii} '= %d \n'],  p(3));
            fprintf(1, ['p value for anovan delta dB histogram per mouse per electrode for electrodes for ' freq_names{bwii} '= %d \n\n'],  p(4));
            
        end
        
        %Now plot the histograms and the average per mouse LFP power
        %computed per mouse
        conc_anno_loc = {[0.15 0.15 0.2 0.2], [0.15 0.15 0.2 0.17], [0.15 0.15 0.2 0.14], [0.15 0.15 0.2 0.11], [0.15 0.15 0.2 0.08], [0.15 0.15 0.2 0.05]};
        fig_pos = {[664 550 576 513],[1233 550 576 513],[664 36 576 513],[1233 36 576 513]};
        for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            %Plot the average
            
            try
                close(bwii+8)
            catch
            end
            figure(bwii+8)
            
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
            hold on
            
            per_mouse_data_delta_dB=[];
            per_mouse_prof_naive=[];
            per_mouse_events=[];
            
            for evNo=1:length(eventType)
                
                for per_ii=1:2      %performance bins. blue = naive, red = proficient
                    
                    bar_offset=21-evNo*3+(2-per_ii);
                    
                    %Compute per mouse avearge
                    each_mouse_average_delta_dB=[];
                    for mouseNo=1:max(delta_dB_mouseNo_per_mouse)
                        if mouse_included(mouseNo)==1
                            each_mouse_average_delta_dB(mouseNo)=mean(delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_mouseNo_per_mouse==mouseNo)));
                        end
                    end
                    
                    if per_ii==1
                        bar(bar_offset,mean(each_mouse_average_delta_dB(logical(mouse_included))),'r','LineWidth', 3)
                    else
                        bar(bar_offset,mean(each_mouse_average_delta_dB(logical(mouse_included))),'b','LineWidth', 3)
                    end
                    
                    
                    %In the future add lines linking the points
                    plot((bar_offset)*ones(1,sum(mouse_included)),each_mouse_average_delta_dB(logical(mouse_included)),'o',...
                        'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                    
                    %Average and CI
                    plot(bar_offset,mean(each_mouse_average_delta_dB(logical(mouse_included))),'ok','LineWidth', 3)
                    if sum(mouse_included)>2
                        CI = bootci(1000, {@mean, each_mouse_average_delta_dB(logical(mouse_included))},'type','cper');
                        plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                    end
                    
                    annotation('textbox',conc_anno_loc{per_ii},'String',per_lab(per_ii),'Color',these_colors{3-per_ii},'EdgeColor','none');
                    
                    %Save data for anovan
                    per_mouse_data_delta_dB=[per_mouse_data_delta_dB each_mouse_average_delta_dB(logical(mouse_included))];
                    per_mouse_prof_naive=[per_mouse_prof_naive per_ii*ones(1,sum(mouse_included))];
                    per_mouse_events=[per_mouse_events evNo*ones(1,sum(mouse_included))];
                end
            end
            title([freq_names{bwii} ' average delta dB power per mouse, electrode avearage'])
            set(gcf,'OuterPosition',fig_pos{bwii});
            
            if sum(eventType==3)>0
                xticks([15,16,18,19])
                xticklabels({evTypeLabels{2},evTypeLabels{2},evTypeLabels{1},evTypeLabels{1}})
            else
                bar_lab_loc = [3.5 6.5 9.5 12.5 15.5 18.5];
                xticks(bar_lab_loc)
                xticklabels(concs2)
                xlabel('Concentration (%)')
            end
            
            ylabel('Delta power (dB)')
            
            
            
            %Calculate anovan for inteaction
            [p,tbl,stats]=anovan(per_mouse_data_delta_dB,{per_mouse_prof_naive per_mouse_events},'varnames',{'proficient_vs_naive','events'},'display','off');
            fprintf(1, ['p value for anovan delta dB histogram per mouse, electrode avearage for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan delta dB histogram  per mouse, electrode avearage for events ' freq_names{bwii} '= %d \n\n'],  p(2));
            
        end
        
        
        %Do auROC calculations per electrode per mouse
        
        %Plot cumulative histos for auROCs within vs between S+ and S-
        for bwii=1:no_bandwidths
            
            
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            hold on
            
            %Naive between
            if length(auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1)))>3
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1)));
                plot(x_aic,f_aic,'b')
            end
            
            %Proficient between
            if length(auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1)))>3
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1)));
                plot(x_aic,f_aic,'r')
            end
            
            %Naive within
            if length(auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0)))>3
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0)));
                if sum(eventType==3)>0
                    plot(x_aic,f_aic,'b')
                else
                    plot(x_aic,f_aic,'Color',[0.8 0.8 1])
                end
            end
            
            %Proficient within
            if length(auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0)))>3
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0)));
                if sum(eventType==3)>0
                    plot(x_aic,f_aic,'r')
                else
                    plot(x_aic,f_aic,'Color',[1 0.8 0.8])
                end
            end
            
            if sum(eventType==3)>0
                legend('Naive','Proficient')
            else
                legend('Naive between','Proficient between','Naive within','Proficient within')
            end
            xlabel('auROC')
            ylabel('Cumulative probability')
            title(['Per electrode per mouse auROC for ' freq_names{bwii}])
            
            %Save the data for anovan
            data_auROC=[];
            prof_naive=[];
            between=[];
            
            %Naive between
            data_auROC=[data_auROC auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1))];
            prof_naive=[prof_naive zeros(1,sum((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1)))];
            between=[between ones(1,sum((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1)))];
            
            %Proficient between
            data_auROC=[data_auROC auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1))];
            prof_naive=[prof_naive ones(1,sum((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1)))];
            between=[between ones(1,sum((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1)))];
            
            %Naive within
            data_auROC=[data_auROC auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0))];
            prof_naive=[prof_naive zeros(1,sum((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0)))];
            between=[between zeros(1,sum((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0)))];
            
            %Proficient within
            data_auROC=[data_auROC auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0))];
            prof_naive=[prof_naive ones(1,sum((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0)))];
            between=[between zeros(1,sum((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0)))];
            
            %Calculate anovan for inteaction
            [p,tbl,stats]=anovan(data_auROC,{prof_naive between},'model','interaction','varnames',{'proficient_vs_naive','within_vs_between'},'display','off');
            fprintf(1, ['p value for anovan per electrode per mouse auROC for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan per electrode per mouse auROC for within vs between for ' freq_names{bwii} '= %d \n'],  p(2));
            p_anova_np(bwii)=p(1);
            p_anova_wb(bwii)=p(2);
        end
        
        if sum(ROC_between)>0
            
            %Plot cumulative histos for auROCs within vs between S+ and S- using
            %only ROCs for adjacent concentrations
            for bwii=1:no_bandwidths
                
                
                figNo = get(gcf,'Number')+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                
                hold on
                
                %Naive between
                if ~isempty(auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==1)))
                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==1)));
                    plot(x_aic,f_aic,'b')
                end
                
                %Proficient between
                if ~isempty(auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==1)))
                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==1)));
                    plot(x_aic,f_aic,'r')
                end
                
                %Naive within
                if ~isempty(auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==1)))
                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==1)));
                    plot(x_aic,f_aic,'Color',[0.8 0.8 1])
                end
                
                %Proficient within
                if ~isempty(auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==1)))
                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==1)));
                    plot(x_aic,f_aic,'Color',[1 0.8 0.8])
                end
                
                legend('Naive between','Proficient between','Naive within','Proficient within')
                xlabel('auROC')
                ylabel('Cumulative probability')
                title([freq_names{bwii} ' Adjacent concentrations'])
                
                %Save the data for anovan for adjacent ROCs
                data_auROC=[];
                prof_naive=[];
                between=[];
                
                %Naive between
                data_auROC=[data_auROC auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==1))];
                prof_naive=[prof_naive zeros(1,sum((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==1)))];
                between=[between ones(1,sum((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==1)))];
                
                %Proficient between
                data_auROC=[data_auROC auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1&(ROC_neighbor==1)))];
                prof_naive=[prof_naive ones(1,sum((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==1)))];
                between=[between ones(1,sum((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==1)))];
                
                %Naive within
                data_auROC=[data_auROC auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==1))];
                prof_naive=[prof_naive zeros(1,sum((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==1)))];
                between=[between zeros(1,sum((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==1)))];
                
                %Proficient within
                data_auROC=[data_auROC auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==1))];
                prof_naive=[prof_naive ones(1,sum((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==1)))];
                between=[between zeros(1,sum((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==1)))];
                
                %Calculate anovan for inteaction
                [p,tbl,stats]=anovan(data_auROC,{prof_naive between},'model','interaction','varnames',{'proficient_vs_naive','within_vs_between'},'display','off');
                fprintf(1, ['p value for anovan auROC adjacent concentrations for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
                fprintf(1, ['p value for anovan auROC adjacent concentrations  for within vs between for ' freq_names{bwii} '= %d \n'],  p(2));
                p_anova_np_adj(bwii)=p(1);
                p_anova_wb_adj(bwii)=p(2);
            end
            
            %Plot cumulative histos for auROCs within vs between S+ and S- using
            %only ROCs for concentrations separated by two log steps
            for bwii=1:no_bandwidths
                
                
                figNo = get(gcf,'Number')+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                
                hold on
                
                %Naive between
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==2)));
                plot(x_aic,f_aic,'b')
                
                %Proficient between
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==2)));
                plot(x_aic,f_aic,'r')
                
                %Naive within
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==2)));
                plot(x_aic,f_aic,'Color',[0.8 0.8 1])
                
                %Proficient within
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==2)));
                plot(x_aic,f_aic,'Color',[1 0.8 0.8])
                
                legend('Naive between','Proficient between','Naive within','Proficient within')
                xlabel('auROC')
                ylabel('Cumulative probability')
                title([freq_names{bwii} ' concentrations separated by two log steps'])
                
                %Save the data for anovan for adjacent ROCs
                data_auROC=[];
                prof_naive=[];
                between=[];
                
                %Naive between
                data_auROC=[data_auROC auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==2))];
                prof_naive=[prof_naive zeros(1,sum((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==2)))];
                between=[between ones(1,sum((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==2)))];
                
                %Proficient between
                data_auROC=[data_auROC auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1&(ROC_neighbor==2)))];
                prof_naive=[prof_naive ones(1,sum((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==2)))];
                between=[between ones(1,sum((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==2)))];
                
                %Naive within
                data_auROC=[data_auROC auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==2))];
                prof_naive=[prof_naive zeros(1,sum((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==2)))];
                between=[between zeros(1,sum((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==2)))];
                
                %Proficient within
                data_auROC=[data_auROC auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==2))];
                prof_naive=[prof_naive ones(1,sum((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==2)))];
                between=[between zeros(1,sum((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==2)))];
                
                %Calculate anovan for inteaction
                [p,tbl,stats]=anovan(data_auROC,{prof_naive between},'model','interaction','varnames',{'proficient_vs_naive','within_vs_between'},'display','off');
                fprintf(1, ['p value for anovan auROC two log step concentrations for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
                fprintf(1, ['p value for anovan auROC two log step concentrations  for within vs between for ' freq_names{bwii} '= %d \n'],  p(2));
                p_anova_np_adj(bwii)=p(1);
                p_anova_wb_adj(bwii)=p(2);
            end
        end
        
        %Plot percent significant ROC
        no_within=zeros(2,4);
        no_sig_within=zeros(2,4);
        no_between=zeros(2,4);
        no_sig_between=zeros(2,4);
        no_within1=zeros(2,4);
        no_sig_within1=zeros(2,4);
        no_between1=zeros(2,4);
        no_sig_between1=zeros(2,4);
        no_within2=zeros(2,4);
        no_sig_within2=zeros(2,4);
        no_between2=zeros(2,4);
        no_sig_between2=zeros(2,4);
        
        for bwii=1:no_bandwidths
            for pcii=1:szpc(1)
                no_pairs=0;
                EvNo1=[];
                EvNo2=[];
                per_sig=zeros(length(eventType),length(eventType));
                for evNo1=1:length(eventType)
                    for evNo2=evNo1+1:length(eventType)
                        
                        
                        no_pairs=no_pairs+1;
                        these_ROCs=(ROCbandwidth==bwii)&(ROCEvNo1==evNo1)&(ROCEvNo2==evNo2)&(ROCper_ii==pcii);
                        sig(no_pairs)=sum((p_valROC<=pFDRauROC)&these_ROCs);
                        not_sig(no_pairs)=sum(these_ROCs)-sum((p_valROC<=pFDRauROC)&these_ROCs);
                        EvNo1(no_pairs)=evNo1;
                        EvNo2(no_pairs)=evNo2;
                        per_sig(evNo1,evNo2)=100*sig(no_pairs)/(sig(no_pairs)+not_sig(no_pairs));
                        
                        
                        if ((evNo1==1)&(evNo2==2))||((evNo1==2)&(evNo2==1))||((evNo1==1)&(evNo2==3))||((evNo1==3)&(evNo2==1))...
                                ||((evNo1==2)&(evNo2==3))||((evNo1==3)&(evNo2==2))||((evNo1==4)&(evNo2==5))||((evNo1==5)&(evNo2==4))...
                                ||((evNo1==4)&(evNo2==6))||((evNo1==6)&(evNo2==4))||((evNo1==5)&(evNo2==6))||((evNo1==6)&(evNo2==5))
                            %This is within
                            no_sig_within(pcii,bwii)=no_sig_within(pcii,bwii)+sum((p_valROC<=pFDRauROC)&these_ROCs);
                            no_within(pcii,bwii)=no_within(pcii,bwii)+sum(these_ROCs);
                            if abs(evNo1-evNo2)==1
                                no_sig_within1(pcii,bwii)=no_sig_within1(pcii,bwii)+sum((p_valROC<=pFDRauROC)&these_ROCs);
                                no_within1(pcii,bwii)=no_within1(pcii,bwii)+sum(these_ROCs);
                            end
                            if abs(evNo1-evNo2)==2
                                no_sig_within2(pcii,bwii)=no_sig_within2(pcii,bwii)+sum((p_valROC<=pFDRauROC)&these_ROCs);
                                no_within2(pcii,bwii)=no_within2(pcii,bwii)+sum(these_ROCs);
                            end
                        else
                            %This is bewteen
                            no_sig_between(pcii,bwii)=no_sig_between(pcii,bwii)+sum((p_valROC<=pFDRauROC)&these_ROCs);
                            no_between(pcii,bwii)=no_between(pcii,bwii)+sum(these_ROCs);
                            if abs(evNo1-evNo2)==1
                                no_sig_between1(pcii,bwii)=no_sig_between1(pcii,bwii)+sum((p_valROC<=pFDRauROC)&these_ROCs);
                                no_between1(pcii,bwii)=no_between1(pcii,bwii)+sum(these_ROCs);
                            end
                            if abs(evNo1-evNo2)==2
                                no_sig_between2(pcii,bwii)=no_sig_between2(pcii,bwii)+sum((p_valROC<=pFDRauROC)&these_ROCs);
                                no_between2(pcii,bwii)=no_between2(pcii,bwii)+sum(these_ROCs);
                            end
                        end
                        
                    end
                end
                
                for evNo1=1:length(eventType)
                    for evNo2=evNo1+1:length(eventType)
                        
                        if per_sig(evNo1,evNo2)==0
                            per_sig(evNo1,evNo2)=100/64;
                        end
                    end
                end
                
                %Plot the pseudocolor for percent significant auROCs
                figNo = get(gcf,'Number')+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                
                evNos_for_1=1:length(eventType);
                evNos_for_2=[1:length(eventType)]';
                
                drg_pcolor(repmat(evNos_for_1,length(eventType),1),repmat(evNos_for_2,1,length(eventType)),per_sig)
                cmjet=colormap(jet);
                cmjet(1,1)=0.7;
                cmjet(1,2)=0.7;
                cmjet(1,3)=0.7;
                colormap(cmjet)
                
                caxis([0 100]);
                
                hold on
                plot([4 4],[1 4],'-w','LineWidth', 5)
                plot([4 7],[4 4],'-w','LineWidth', 5)
                
                ax=gca;
                set(ax,'XTickLabel','')
                ylabel('dB')
                xticks([1.5:1:length(eventType)+1])
                xticklabels(handles_pars.concs2)
                yticks([1.5:1:length(eventType)+1])
                yticklabels(handles_pars.concs2)
                
                title(['Percent auROC significantly different from zero ' freq_names{bwii} ' ' per_lab(pcii)])
                set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                
                pfft=1;
            end
        end
        
        
        
        if sum(ROC_between)>0
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            
            set(hFig, 'units','normalized','position',[.83 .1 .05 .3])
            
            prain=[0:100/99:100];
            drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
            colormap jet
            shading interp
            ax=gca;
            set(ax,'XTickLabel','')
            
            %Plot percent significant ROCs
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            for bwii=1:no_bandwidths
                subplot(2,2,bwii)
                hold on
                %Plot within naive
                %             bar(1,100*no_sig_within(2,bwii)/no_within(2,bwii),'b')
                %             bar(2,100*no_sig_within1(2,bwii)/no_within1(2,bwii),'b')
                bar(3,100*no_sig_within2(2,bwii)/no_within2(2,bwii),'b')
                
                %Plot within proficient
                %             bar(4,100*no_sig_within(1,bwii)/no_within(1,bwii),'r')
                %             bar(5,100*no_sig_within1(1,bwii)/no_within1(1,bwii),'r')
                bar(6,100*no_sig_within2(1,bwii)/no_within2(1,bwii),'r')
                
                %             %Plot between naive
                %             bar(9,100*no_sig_between(2,bwii)/no_between(2,bwii),'b')
                %             bar(10,100*no_sig_between1(2,bwii)/no_between1(2,bwii),'b')
                bar(11,100*no_sig_between2(2,bwii)/no_between2(2,bwii),'b')
                
                %Plot between proficient
                %             bar(12,100*no_sig_between(1,bwii)/no_between(1,bwii),'r')
                %             bar(13,100*no_sig_between1(1,bwii)/no_between1(1,bwii),'r')
                bar(14,100*no_sig_between2(1,bwii)/no_between2(1,bwii),'r')
                
                title(['Percent significant auROC ' freq_names{bwii} ])
            end
            
        else
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            
            set(hFig, 'units','normalized','position',[.83 .1 .05 .3])
            
            prain=[0:100/99:100];
            drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
            colormap jet
            shading interp
            ax=gca;
            set(ax,'XTickLabel','')
            
            %Plot percent significant ROCs
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            for bwii=1:no_bandwidths
                subplot(2,2,bwii)
                hold on
                %Plot within naive
                bar(1,100*no_sig_within2(2,bwii)/no_within(2,bwii),'b')
                %             bar(2,100*no_sig_within1(2,bwii)/no_within1(2,bwii),'b')
                %           bar(3,100*no_sig_within2(2,bwii)/no_within2(2,bwii),'b')
                
                %Plot within proficient
                bar(2,100*no_sig_within2(1,bwii)/no_within(1,bwii),'r')
                %             bar(5,100*no_sig_within1(1,bwii)/no_within1(1,bwii),'r')
                %           bar(6,100*no_sig_within2(1,bwii)/no_within2(1,bwii),'r')
                
                %             %Plot between naive
                %             bar(9,100*no_sig_between(2,bwii)/no_between(2,bwii),'b')
                %             bar(10,100*no_sig_between1(2,bwii)/no_between1(2,bwii),'b')
                %bar(11,100*no_sig_between2(2,bwii)/no_between2(2,bwii),'b')
                
                %Plot between proficient
                %             bar(12,100*no_sig_between(1,bwii)/no_between(1,bwii),'r')
                %             bar(13,100*no_sig_between1(1,bwii)/no_between1(1,bwii),'r')
                %                 bar(14,100*no_sig_between2(1,bwii)/no_between2(1,bwii),'r')
                
                title(['Percent significant auROC ' freq_names{bwii} ])
                legend('Naive','Proficient')
            end
            
        end
        
        %Do auROC calulations per mouse (avearge of electrodes)
        
        %Plot cumulative histos for auROCs within vs between S+ and S-
        for bwii=1:no_bandwidths
            
            
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            hold on
            
            %Naive between
            if length(per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)))>3
                [f_aic,x_aic] = drg_ecdf(per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)));
                plot(x_aic,f_aic,'b')
            end
            
            %Proficient between
            if length(per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)))>3
                [f_aic,x_aic] = drg_ecdf(per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)));
                plot(x_aic,f_aic,'r')
            end
            
            %Naive within
            if length(per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)))>3
                [f_aic,x_aic] = drg_ecdf(per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)));
                plot(x_aic,f_aic,'Color',[0.8 0.8 1])
            end
            
            %Proficient within
            if length(per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)))>3
                [f_aic,x_aic] = drg_ecdf(per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)));
                plot(x_aic,f_aic,'Color',[1 0.8 0.8])
            end
            
            legend('Naive between','Proficient between','Naive within','Proficient within')
            xlabel('per_mouse_auROC')
            ylabel('Cumulative probability')
            title(['auROC (electrode average) ' freq_names{bwii}])
            
            %Save the data for anovan
            data_per_mouse_auROC=[];
            prof_naive=[];
            between=[];
            
            %Naive between
            data_per_mouse_auROC=[data_per_mouse_auROC per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1))];
            prof_naive=[prof_naive zeros(1,sum((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)))];
            between=[between ones(1,sum((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)))];
            
            %Proficient between
            data_per_mouse_auROC=[data_per_mouse_auROC per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1))];
            prof_naive=[prof_naive ones(1,sum((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)))];
            between=[between ones(1,sum((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)))];
            
            %Naive within
            data_per_mouse_auROC=[data_per_mouse_auROC per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0))];
            prof_naive=[prof_naive zeros(1,sum((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)))];
            between=[between zeros(1,sum((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)))];
            
            %Proficient within
            data_per_mouse_auROC=[data_per_mouse_auROC per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0))];
            prof_naive=[prof_naive ones(1,sum((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)))];
            between=[between zeros(1,sum((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)))];
            
            %Calculate anovan for inteaction
            [p,tbl,stats]=anovan(data_per_mouse_auROC,{prof_naive between},'model','interaction','varnames',{'proficient_vs_naive','within_vs_between'},'display','off');
            fprintf(1, ['p value for anovan auROC electrode average for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan auROC electrode average for within vs between for ' freq_names{bwii} '= %d \n'],  p(2));
            p_anova_np(bwii)=p(1);
            p_anova_wb(bwii)=p(2);
        end
        
        if sum(per_mouse_ROC_between)>0
            
            %Plot cumulative histos for auROCs within vs between S+ and S- using
            %only ROCs for adjacent concentrations
            for bwii=1:no_bandwidths
                
                
                figNo = get(gcf,'Number')+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                
                hold on
                
                %Naive between
                if ~isempty(per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==1)))
                    [f_aic,x_aic] = drg_ecdf(per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==1)));
                    plot(x_aic,f_aic,'b')
                end
                
                %Proficient between
                if ~isempty(per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==1)))
                    [f_aic,x_aic] = drg_ecdf(per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==1)));
                    plot(x_aic,f_aic,'r')
                end
                
                %Naive within
                if ~isempty(per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==1)))
                    [f_aic,x_aic] = drg_ecdf(per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==1)));
                    plot(x_aic,f_aic,'Color',[0.8 0.8 1])
                end
                
                %Proficient within
                if ~isempty(per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==1)))
                    [f_aic,x_aic] = drg_ecdf(per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==1)));
                    plot(x_aic,f_aic,'Color',[1 0.8 0.8])
                end
                
                legend('Naive between','Proficient between','Naive within','Proficient within')
                xlabel('auROC')
                ylabel('Cumulative probability')
                title([freq_names{bwii} ' Adjacent concentrations (electrode average)'])
                
                %Save the data for anovan for adjacent ROCs
                data_per_mouse_auROC=[];
                prof_naive=[];
                between=[];
                
                %Naive between
                data_per_mouse_auROC=[data_per_mouse_auROC per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==1))];
                prof_naive=[prof_naive zeros(1,sum((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==1)))];
                between=[between ones(1,sum((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==1)))];
                
                %Proficient between
                data_per_mouse_auROC=[data_per_mouse_auROC per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1&(per_mouse_ROC_neighbor==1)))];
                prof_naive=[prof_naive ones(1,sum((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==1)))];
                between=[between ones(1,sum((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==1)))];
                
                %Naive within
                data_per_mouse_auROC=[data_per_mouse_auROC per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==1))];
                prof_naive=[prof_naive zeros(1,sum((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==1)))];
                between=[between zeros(1,sum((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==1)))];
                
                %Proficient within
                data_per_mouse_auROC=[data_per_mouse_auROC per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==1))];
                prof_naive=[prof_naive ones(1,sum((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==1)))];
                between=[between zeros(1,sum((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==1)))];
                
                %Calculate anovan for inteaction
                [p,tbl,stats]=anovan(data_per_mouse_auROC,{prof_naive between},'model','interaction','varnames',{'proficient_vs_naive','within_vs_between'},'display','off');
                fprintf(1, ['p value for anovan auROC electrode average adjacent concentrations for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
                fprintf(1, ['p value for anovan auROC electrode average adjacent concentrations  for within vs between for ' freq_names{bwii} '= %d \n'],  p(2));
                p_anova_np_adj(bwii)=p(1);
                p_anova_wb_adj(bwii)=p(2);
            end
            
            %Plot cumulative histos for auROCs within vs between S+ and S- using
            %only ROCs for concentrations separated by two log steps
            for bwii=1:no_bandwidths
                
                
                figNo = get(gcf,'Number')+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                
                hold on
                
                %Naive between
                [f_aic,x_aic] = drg_ecdf(per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==2)));
                plot(x_aic,f_aic,'b')
                
                %Proficient between
                [f_aic,x_aic] = drg_ecdf(per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==2)));
                plot(x_aic,f_aic,'r')
                
                %Naive within
                [f_aic,x_aic] = drg_ecdf(per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==2)));
                plot(x_aic,f_aic,'Color',[0.8 0.8 1])
                
                %Proficient within
                [f_aic,x_aic] = drg_ecdf(per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==2)));
                plot(x_aic,f_aic,'Color',[1 0.8 0.8])
                
                legend('Naive between','Proficient between','Naive within','Proficient within')
                xlabel('auROC')
                ylabel('Cumulative probability')
                title([freq_names{bwii} ' concentrations separated by two log steps (electrode average)'])
                
                %Save the data for anovan for adjacent ROCs
                data_per_mouse_auROC=[];
                prof_naive=[];
                between=[];
                
                %Naive between
                data_per_mouse_auROC=[data_per_mouse_auROC per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==2))];
                prof_naive=[prof_naive zeros(1,sum((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==2)))];
                between=[between ones(1,sum((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==2)))];
                
                %Proficient between
                data_per_mouse_auROC=[data_per_mouse_auROC per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1&(per_mouse_ROC_neighbor==2)))];
                prof_naive=[prof_naive ones(1,sum((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==2)))];
                between=[between ones(1,sum((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==1)&(per_mouse_ROC_neighbor==2)))];
                
                %Naive within
                data_per_mouse_auROC=[data_per_mouse_auROC per_mouse_auROC((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==2))];
                prof_naive=[prof_naive zeros(1,sum((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==2)))];
                between=[between zeros(1,sum((per_mouse_ROCper_ii==2)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==2)))];
                
                %Proficient within
                data_per_mouse_auROC=[data_per_mouse_auROC per_mouse_auROC((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==2))];
                prof_naive=[prof_naive ones(1,sum((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==2)))];
                between=[between zeros(1,sum((per_mouse_ROCper_ii==1)&(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROC_between==0)&(per_mouse_ROC_neighbor==2)))];
                
                %Calculate anovan for inteaction
                [p,tbl,stats]=anovan(data_per_mouse_auROC,{prof_naive between},'model','interaction','varnames',{'proficient_vs_naive','within_vs_between'},'display','off');
                fprintf(1, ['p value for anovan auROC two log step concentrations electrode average for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
                fprintf(1, ['p value for anovan auROC two log step concentrations electrode average for within vs between for ' freq_names{bwii} '= %d \n'],  p(2));
                p_anova_np_adj(bwii)=p(1);
                p_anova_wb_adj(bwii)=p(2);
            end
        end
        %Plot percent significant ROC
        no_within=zeros(2,4);
        no_sig_within=zeros(2,4);
        no_between=zeros(2,4);
        no_sig_between=zeros(2,4);
        no_within1=zeros(2,4);
        no_sig_within1=zeros(2,4);
        no_between1=zeros(2,4);
        no_sig_between1=zeros(2,4);
        no_within2=zeros(2,4);
        no_sig_within2=zeros(2,4);
        no_between2=zeros(2,4);
        no_sig_between2=zeros(2,4);
        
        for bwii=1:no_bandwidths
            for pcii=1:szpc(1)
                no_pairs=0;
                EvNo1=[];
                EvNo2=[];
                per_sig=zeros(length(eventType),length(eventType));
                for evNo1=1:length(eventType)
                    for evNo2=evNo1+1:length(eventType)
                        
                        
                        no_pairs=no_pairs+1;
                        these_ROCs=(per_mouse_ROCbandwidth==bwii)&(per_mouse_ROCEvNo1==evNo1)&(per_mouse_ROCEvNo2==evNo2)&(per_mouse_ROCper_ii==pcii);
                        sig(no_pairs)=sum((per_mouse_p_valROC<=per_mouse_pFDRauROC)&these_ROCs);
                        not_sig(no_pairs)=sum(these_ROCs)-sum((per_mouse_p_valROC<=per_mouse_pFDRauROC)&these_ROCs);
                        EvNo1(no_pairs)=evNo1;
                        EvNo2(no_pairs)=evNo2;
                        per_sig(evNo1,evNo2)=100*sig(no_pairs)/(sig(no_pairs)+not_sig(no_pairs));
                        
                        
                        if ((evNo1==1)&(evNo2==2))||((evNo1==2)&(evNo2==1))||((evNo1==1)&(evNo2==3))||((evNo1==3)&(evNo2==1))...
                                ||((evNo1==2)&(evNo2==3))||((evNo1==3)&(evNo2==2))||((evNo1==4)&(evNo2==5))||((evNo1==5)&(evNo2==4))...
                                ||((evNo1==4)&(evNo2==6))||((evNo1==6)&(evNo2==4))||((evNo1==5)&(evNo2==6))||((evNo1==6)&(evNo2==5))
                            %This is within
                            no_sig_within(pcii,bwii)=no_sig_within(pcii,bwii)+sum((per_mouse_p_valROC<=per_mouse_pFDRauROC)&these_ROCs);
                            no_within(pcii,bwii)=no_within(pcii,bwii)+sum(these_ROCs);
                            if abs(evNo1-evNo2)==1
                                no_sig_within1(pcii,bwii)=no_sig_within1(pcii,bwii)+sum((per_mouse_p_valROC<=per_mouse_pFDRauROC)&these_ROCs);
                                no_within1(pcii,bwii)=no_within1(pcii,bwii)+sum(these_ROCs);
                            end
                            if abs(evNo1-evNo2)==2
                                no_sig_within2(pcii,bwii)=no_sig_within2(pcii,bwii)+sum((per_mouse_p_valROC<=per_mouse_pFDRauROC)&these_ROCs);
                                no_within2(pcii,bwii)=no_within2(pcii,bwii)+sum(these_ROCs);
                            end
                        else
                            %This is bewteen
                            no_sig_between(pcii,bwii)=no_sig_between(pcii,bwii)+sum((per_mouse_p_valROC<=per_mouse_pFDRauROC)&these_ROCs);
                            no_between(pcii,bwii)=no_between(pcii,bwii)+sum(these_ROCs);
                            if abs(evNo1-evNo2)==1
                                no_sig_between1(pcii,bwii)=no_sig_between1(pcii,bwii)+sum((per_mouse_p_valROC<=per_mouse_pFDRauROC)&these_ROCs);
                                no_between1(pcii,bwii)=no_between1(pcii,bwii)+sum(these_ROCs);
                            end
                            if abs(evNo1-evNo2)==2
                                no_sig_between2(pcii,bwii)=no_sig_between2(pcii,bwii)+sum((per_mouse_p_valROC<=per_mouse_pFDRauROC)&these_ROCs);
                                no_between2(pcii,bwii)=no_between2(pcii,bwii)+sum(these_ROCs);
                            end
                        end
                        
                    end
                end
                
                for evNo1=1:length(eventType)
                    for evNo2=evNo1+1:length(eventType)
                        
                        if per_sig(evNo1,evNo2)==0
                            per_sig(evNo1,evNo2)=100/64;
                        end
                    end
                end
                
                %Plot the pseudocolor for percent significant auROCs
                figNo = get(gcf,'Number')+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                
                evNos_for_1=1:length(eventType);
                evNos_for_2=[1:length(eventType)]';
                
                drg_pcolor(repmat(evNos_for_1,length(eventType),1),repmat(evNos_for_2,1,length(eventType)),per_sig)
                cmjet=colormap(jet);
                cmjet(1,1)=0.7;
                cmjet(1,2)=0.7;
                cmjet(1,3)=0.7;
                colormap(cmjet)
                
                caxis([0 100]);
                
                hold on
                plot([4 4],[1 4],'-w','LineWidth', 5)
                plot([4 7],[4 4],'-w','LineWidth', 5)
                
                ax=gca;
                set(ax,'XTickLabel','')
                ylabel('dB')
                xticks([1.5:1:length(eventType)+1])
                xticklabels(handles_pars.concs2)
                yticks([1.5:1:length(eventType)+1])
                yticklabels(handles_pars.concs2)
                
                title(['Percent significant auROC (electrode average) ' freq_names{bwii} ' ' per_lab(pcii)])
                set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                
                pfft=1;
            end
        end
        
        
        
        if sum(per_mouse_ROC_between)>0
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            
            set(hFig, 'units','normalized','position',[.83 .1 .05 .3])
            
            prain=[0:100/99:100];
            drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
            colormap jet
            shading interp
            ax=gca;
            set(ax,'XTickLabel','')
            
            %Plot percent significant ROCs
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            for bwii=1:no_bandwidths
                subplot(2,2,bwii)
                hold on
                %Plot within naive
                %             bar(1,100*no_sig_within(2,bwii)/no_within(2,bwii),'b')
                %             bar(2,100*no_sig_within1(2,bwii)/no_within1(2,bwii),'b')
                bar(3,100*no_sig_within2(2,bwii)/no_within2(2,bwii),'b')
                
                %Plot within proficient
                %             bar(4,100*no_sig_within(1,bwii)/no_within(1,bwii),'r')
                %             bar(5,100*no_sig_within1(1,bwii)/no_within1(1,bwii),'r')
                bar(6,100*no_sig_within2(1,bwii)/no_within2(1,bwii),'r')
                
                %             %Plot between naive
                %             bar(9,100*no_sig_between(2,bwii)/no_between(2,bwii),'b')
                %             bar(10,100*no_sig_between1(2,bwii)/no_between1(2,bwii),'b')
                bar(11,100*no_sig_between2(2,bwii)/no_between2(2,bwii),'b')
                
                %Plot between proficient
                %             bar(12,100*no_sig_between(1,bwii)/no_between(1,bwii),'r')
                %             bar(13,100*no_sig_between1(1,bwii)/no_between1(1,bwii),'r')
                bar(14,100*no_sig_between2(1,bwii)/no_between2(1,bwii),'r')
                
                title(['Percent significant auROC electrode average ' freq_names{bwii} ])
            end
            
        else
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            
            set(hFig, 'units','normalized','position',[.83 .1 .05 .3])
            
            prain=[0:100/99:100];
            drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
            colormap jet
            shading interp
            ax=gca;
            set(ax,'XTickLabel','')
            
            %Plot percent significant ROCs
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            for bwii=1:no_bandwidths
                subplot(2,2,bwii)
                hold on
                %Plot within naive
                bar(1,100*no_sig_within2(2,bwii)/no_within(2,bwii),'b')
                %             bar(2,100*no_sig_within1(2,bwii)/no_within1(2,bwii),'b')
                %           bar(3,100*no_sig_within2(2,bwii)/no_within2(2,bwii),'b')
                
                %Plot within proficient
                bar(2,100*no_sig_within2(1,bwii)/no_within(1,bwii),'r')
                %             bar(5,100*no_sig_within1(1,bwii)/no_within1(1,bwii),'r')
                %           bar(6,100*no_sig_within2(1,bwii)/no_within2(1,bwii),'r')
                
                %             %Plot between naive
                %             bar(9,100*no_sig_between(2,bwii)/no_between(2,bwii),'b')
                %             bar(10,100*no_sig_between1(2,bwii)/no_between1(2,bwii),'b')
                %bar(11,100*no_sig_between2(2,bwii)/no_between2(2,bwii),'b')
                
                %Plot between proficient
                %             bar(12,100*no_sig_between(1,bwii)/no_between(1,bwii),'r')
                %             bar(13,100*no_sig_between1(1,bwii)/no_between1(1,bwii),'r')
                %                 bar(14,100*no_sig_between2(1,bwii)/no_between2(1,bwii),'r')
                
                title(['Percent significant auROC electrode average ' freq_names{bwii} ])
                legend('Naive','Proficient')
            end
            
        end
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix]);
        
    case 15
        %Justin
        %Linear fit of delta power to: spm, concentration, previous reward, percent
        %This does the analysis in all the files and DOES not distinguish between groups!!!
        
        
        
        
        no_fitlms=0;
        fitlm_results=[];
        
        
        
        fprintf(1, ['fitlm analysis for Justin''s paper\n\n'])
        
        no_files=max(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        
        szpc=size(percent_windows);
        for per_ii=1:szpc(1)
            for bwii=1:no_bandwidths
                for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
                    for elec=1:16
                        if sum(which_electrodes==elec)>0
                            
                            no_fitlms=no_fitlms+1;
                            fitlm_results(no_fitlms).per_ii=per_ii;
                            fitlm_results(no_fitlms).mouseNo=mouseNo;
                            fitlm_results(no_fitlms).elec=elec;
                            fitlm_results(no_fitlms).no_trials=0;
                            fitlm_results(no_fitlms).bwii=bwii;
                            fitlm_results(no_fitlms).elec=elec;
                            
                            %For each electrode compute the delta dB and
                            %input parameters for the fitlm
                            for fileNo=1:no_files
                                if sum(files==fileNo)>0
                                    if handles_drgb.drgbchoices.mouse_no(fileNo)==mouseNo
                                        lfpodNo_ref=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                                        
                                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref)))
                                            
                                            
                                            if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower))
                                                
                                                percent_mask=[];
                                                percent_mask=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower>=percent_windows(per_ii,1))...
                                                    &(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower<=percent_windows(per_ii,2));
                                                
                                                for evNo=1:length(eventType)
                                                    
                                                    trials_in_event_Ev=[];
                                                    trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(eventType(evNo),:)==1)&percent_mask;
                                                    
                                                    if (sum(trials_in_event_Ev)>=1)
                                                        
                                                        lfpodNo=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                                        
                                                        % Ev1
                                                        this_dB_powerref=zeros(sum(trials_in_event_Ev),length(frequency));
                                                        this_dB_powerref(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower(trials_in_event_Ev,:));
                                                        
                                                        
                                                        this_dB_power=zeros(sum(trials_in_event_Ev),length(frequency));
                                                        this_dB_power(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_event_Ev,:));
                                                        
                                                        
                                                        this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                                        
                                                        %Enter the  Ev1
                                                        this_delta_dB_powerEv=zeros(sum(trials_in_event_Ev),1);
                                                        this_delta_dB_powerEv=mean(this_dB_power(:,this_band)-this_dB_powerref(:,this_band),2);
                                                        
                                                        fitlm_results(no_fitlms).concs(fitlm_results(no_fitlms).no_trials+1:fitlm_results(no_fitlms).no_trials+sum(trials_in_event_Ev))=...
                                                            handles_pars.concs(evNo);
                                                        fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).no_trials+1:fitlm_results(no_fitlms).no_trials+sum(trials_in_event_Ev))=...
                                                            this_delta_dB_powerEv;
                                                        
                                                        fitlm_results(no_fitlms).percent(fitlm_results(no_fitlms).no_trials+1:fitlm_results(no_fitlms).no_trials+sum(trials_in_event_Ev))=...
                                                            handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower(trials_in_event_Ev);
                                                        
                                                        %Is this S= or S-?
                                                        if evNo<=length(eventType)/2
                                                            %This is S+
                                                            fitlm_results(no_fitlms).spm(fitlm_results(no_fitlms).no_trials+1:fitlm_results(no_fitlms).no_trials+sum(trials_in_event_Ev))=...
                                                                ones(1,sum(trials_in_event_Ev));
                                                        else
                                                            %This is S-
                                                            fitlm_results(no_fitlms).spm(fitlm_results(no_fitlms).no_trials+1:fitlm_results(no_fitlms).no_trials+sum(trials_in_event_Ev))=...
                                                                zeros(1,sum(trials_in_event_Ev));
                                                        end
                                                        
                                                        %Was the mouse rewarded
                                                        %in the last trial?
                                                        Hit_trials=[];
                                                        Hit_trials=handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(2,:)==1;
                                                        last_Hit_trials=zeros(1,length(Hit_trials));
                                                        last_Hit_trials(2:end)=Hit_trials(1:end-1);
                                                        fitlm_results(no_fitlms).last_rewarded(fitlm_results(no_fitlms).no_trials+1:fitlm_results(no_fitlms).no_trials+sum(trials_in_event_Ev))=...
                                                            last_Hit_trials(trials_in_event_Ev);
                                                        fitlm_results(no_fitlms).no_trials=fitlm_results(no_fitlms).no_trials+sum(trials_in_event_Ev);
                                                        
                                                        fprintf(1, ['%d trials in event No %d succesfully processed for file No %d electrode %d\n'],sum(trials_in_event_Ev), evNo,files(fileNo),elec);
                                                    end
                                                    
                                                end
                                            else
                                                
                                                fprintf(1, ['Empty allPower for file No %d electrode %d\n'],fileNo,elec);
                                                
                                                
                                            end
                                            
                                            
                                        else
                                            fprintf(1, ['Empty LFP reference for file No %d electrode %d\n'],fileNo,elec);
                                            
                                            
                                        end
                                        
                                        
                                    end
                                end
                            end %file No
                            
                            %For each electrode compute the fitlm
                            %Do fitlm
                            
                            
                            % Store the variables in a table.
                            tbl = table(fitlm_results(no_fitlms).delta_dB_Power',fitlm_results(no_fitlms).concs',fitlm_results(no_fitlms).spm',fitlm_results(no_fitlms).last_rewarded','VariableNames',{'delta_dB','log10_Odor_Concentration','spm','was_last_rewarded'});
                            
                            % Fit a linear regression model
                            fitlm_results(no_fitlms).lm = fitlm(tbl,'delta_dB~log10_Odor_Concentration+spm+was_last_rewarded');
                            
                            %Save the two p values for concentration and
                            %spm
                            fitlm_results(no_fitlms).p_val_conc=fitlm_results(no_fitlms).lm.Coefficients{2,4};
                            fitlm_results(no_fitlms).p_val_spm=fitlm_results(no_fitlms).lm.Coefficients{3,4};
                            fitlm_results(no_fitlms).p_val_last_reward=fitlm_results(no_fitlms).lm.Coefficients{4,4};
                            fitlm_results(no_fitlms).AIC(7)=fitlm_results(no_fitlms).lm.ModelCriterion.AIC;
                            fitlm_results(no_fitlms).BIC(7)=fitlm_results(no_fitlms).lm.ModelCriterion.BIC;
                            
                            % Fit a linear regression model for delta_dB~log10_Odor_Concentration+was_last_rewarded
                            lm = fitlm(tbl,'delta_dB~log10_Odor_Concentration+was_last_rewarded');
                            fitlm_results(no_fitlms).AIC(6)=lm.ModelCriterion.AIC;
                            fitlm_results(no_fitlms).BIC(6)=lm.ModelCriterion.BIC;
                            
                            % Fit a linear regression model for spm+was_last_rewarded
                            lm = fitlm(tbl,'delta_dB~spm+was_last_rewarded');
                            fitlm_results(no_fitlms).AIC(5)=lm.ModelCriterion.AIC;
                            fitlm_results(no_fitlms).BIC(5)=lm.ModelCriterion.BIC;
                            
                            % Fit a linear regression model for delta_dB~log10_Odor_Concentration+spm
                            lm = fitlm(tbl,'delta_dB~log10_Odor_Concentration+spm');
                            fitlm_results(no_fitlms).AIC(4)=lm.ModelCriterion.AIC;
                            fitlm_results(no_fitlms).BIC(4)=lm.ModelCriterion.BIC;
                            
                            % Fit a linear regression model for last
                            % rewarded
                            lm = fitlm(tbl,'delta_dB~was_last_rewarded');
                            fitlm_results(no_fitlms).AIC(3)=lm.ModelCriterion.AIC;
                            fitlm_results(no_fitlms).BIC(3)=lm.ModelCriterion.BIC;
                            
                            % Fit a linear regression model for spm
                            lm = fitlm(tbl,'delta_dB~spm');
                            fitlm_results(no_fitlms).AIC(2)=lm.ModelCriterion.AIC;
                            fitlm_results(no_fitlms).BIC(2)=lm.ModelCriterion.BIC;
                            
                            % Fit a linear regression model for conc
                            lm = fitlm(tbl,'delta_dB~log10_Odor_Concentration');
                            fitlm_results(no_fitlms).AIC(1)=lm.ModelCriterion.AIC;
                            fitlm_results(no_fitlms).BIC(1)=lm.ModelCriterion.BIC;
                            
                            %Plot dB power for each electrode
                            try
                                close(1)
                            catch
                            end
                            hFig1=figure(1);
                            set(hFig1, 'units','normalized','position',[.07 .1 .75 .3])
                            
                            subplot(1,3,1)
                            plot(fitlm_results(no_fitlms).concs',fitlm_results(no_fitlms).delta_dB_Power','ob')
                            minc=min(fitlm_results(no_fitlms).concs');
                            mindB=lm.Coefficients{1,1}+(lm.Coefficients{2,1})*minc;
                            maxc=max(fitlm_results(no_fitlms).concs');
                            maxdB=lm.Coefficients{1,1}+(lm.Coefficients{2,1})*maxc;
                            hold on
                            plot([minc maxc],[mindB maxdB],'-r')
                            xlabel('log(concentration)')
                            ylabel('dB power')
                            
                            subplot(1,3,2)
                            plot(fitlm_results(no_fitlms).spm',fitlm_results(no_fitlms).delta_dB_Power','ob')
                            hold on
                            %Splus bar
                            bar(1,mean(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).spm==1)))
                            %Confidence intervals
                            SEMsp = std(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).spm==1))...
                                /sqrt(length(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).spm==1)));               % Standard Error
                            tssp = tinv([0.05  0.95],sum(fitlm_results(no_fitlms).spm==1)-1);      % T-Score
                            CIsp = mean(mean(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).spm==1))) + tssp*SEMsp;
                            plot([1 1],CIsp,'-r','Linewidth',3)
                            %Sminus bar
                            bar(0,mean(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).spm==0)))% Confidence Intervals
                            %Confidence intervals
                            SEMsm = std(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).spm==0))...
                                /sqrt(length(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).spm==0)));               % Standard Error
                            tssm = tinv([0.05  0.95],sum(fitlm_results(no_fitlms).spm==0)-1);      % T-Score
                            CIsm = mean(mean(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).spm==0))) + tssp*SEMsp;
                            plot([0 0],CIsm,'-r','Linewidth',3)
                            xlabel('S= or S-')
                            ylabel('db power')
                            title(['no_fitlms ' num2str(no_fitlms)])
                            
                            %last_rewarded
                            subplot(1,3,3)
                            plot(fitlm_results(no_fitlms).last_rewarded',fitlm_results(no_fitlms).delta_dB_Power','ob')
                            hold on
                            %Not rewarded bar
                            bar(1,mean(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).last_rewarded==1)))
                            %Confidence intervals
                            SEMsp = std(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).last_rewarded==1))...
                                /sqrt(length(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).last_rewarded==1)));               % Standard Error
                            tssp = tinv([0.05  0.95],sum(fitlm_results(no_fitlms).last_rewarded==1)-1);      % T-Score
                            CIsp = mean(mean(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).last_rewarded==1))) + tssp*SEMsp;
                            plot([1 1],CIsp,'-r','Linewidth',3)
                            %Rewarded bar
                            bar(0,mean(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).last_rewarded==0)))% Confidence Intervals
                            %Confidence intervals
                            SEMsm = std(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).last_rewarded==0))...
                                /sqrt(length(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).last_rewarded==0)));               % Standard Error
                            tssm = tinv([0.05  0.95],sum(fitlm_results(no_fitlms).last_rewarded==0)-1);      % T-Score
                            CIsm = mean(mean(fitlm_results(no_fitlms).delta_dB_Power(fitlm_results(no_fitlms).last_rewarded==0))) + tssp*SEMsp;
                            plot([0 0],CIsm,'-r','Linewidth',3)
                            
                            
                            xlabel('Last rewarded')
                            ylabel('db power')
                            
                            pfft=1;
                            
                        end
                        
                    end
                end %mouseNo
                
            end %bandwidth
        end %end perii
        
        %Now parse through the p values
        for ii=1:no_fitlms
            periis(ii)=fitlm_results(ii).per_ii;
            bwiis(ii)=fitlm_results(ii).bwii;
            pvc(ii)=fitlm_results(ii).p_val_conc;
            pvspm(ii)=fitlm_results(ii).p_val_spm;
            pv_last_rewards(ii)=fitlm_results(ii).p_val_last_reward;
            for jj=1:7
                aic(ii,jj)=fitlm_results(ii).AIC(jj);
                bic(ii,jj)=fitlm_results(ii).BIC(jj);
            end
        end
        
        %Plot a bar graph showing the percent of significant estimated coefficents for the fitlm regression model
        try
            close(2)
        catch
        end
        hFig2=figure(2);
        hold on
        
        xpos=0;
        for perii=1:2
            for bwii=1:no_bandwidths
                %Concentration dependence
                percent_sig_conc=100*sum(pvc((periis==perii)&(bwiis==bwii))<=0.05)/sum((periis==perii)&(bwiis==bwii));
                fprintf(1, ['percent significant for dependence of delta dB on concentration for ' freq_names{bwii} ' ' per_lab{perii} '= %d, n= %d\n' ]...
                    ,percent_sig_conc,sum((periis==perii)&(bwiis==bwii)));
                xpos=xpos+1;
                bar(xpos,percent_sig_conc,'b')
                
                %spm dependence
                percent_sig_spm=100*sum(pvspm((periis==perii)&(bwiis==bwii))<=0.05)/sum((periis==perii)&(bwiis==bwii));
                fprintf(1, ['percent significant for dependence of delta dB on spm for ' freq_names{bwii} ' ' per_lab{perii} '= %d, n= %d\n' ]...
                    ,percent_sig_spm,sum((periis==perii)&(bwiis==bwii)));
                xpos=xpos+1;
                bar(xpos,percent_sig_spm,'r')
                
                %last reward
                percent_sig_last=100*sum(pv_last_rewards((periis==perii)&(bwiis==bwii))<=0.05)/sum((periis==perii)&(bwiis==bwii));
                fprintf(1, ['percent significant for dependence of delta dB on last_rewards for ' freq_names{bwii} ' ' per_lab{perii} '= %d, n= %d\n\n' ]...
                    ,percent_sig_last,sum((periis==perii)&(bwiis==bwii)));
                xpos=xpos+1;
                bar(xpos,percent_sig_last,'m')
                xpos=xpos+2;
            end
            xpos=xpos+4;
        end
        
        title('blue=concentration, red=spm, magenta=last reward, proficient left, naive right')
        ylabel('Percent significant regression coefficients')
        
        %Plot AICs
        figNo=2;
        for perii=1:2
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2=figure(figNo);
            
            if perii==1
                title('AICs for naive mice')
            else
                title('AICs for proficient mice')
            end
            
            
            
            for bwii=1:no_bandwidths
                subplot(1,4,bwii)
                hold on
                for ii=1:7
                    [f_aic,x_aic] = drg_ecdf(aic((periis==perii)&(bwiis==bwii),ii));
                    plot(x_aic,f_aic,these_colors{ii})
                end
                
            end
            
            
        end
        
        
        %We should use aic() Akaike's Information Criterion for estimated
        %model with either conc or spm
        %https://en.wikipedia.org/wiki/Akaike_information_criterion
        %This is working well with fitlm, but should we try fitglm?
        %https://www.mathworks.com/help/stats/generalized-linear-regression.html
        %Should we make last rewarded 2 if the last two were rewarded?
        
        fprintf(1, '\n\n')
        
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix]);
        
        
    case 16
        %Justin
        %fitglm of delta power to: spm, concentration, previous reward, percent
        %This does the analysis in all the files and DOES not distinguish between groups!!!
        %https://en.wikipedia.org/wiki/Generalized_linear_model
        
        
        
        
        no_fitglms=0;
        fitglm_results=[];
        
        
        
        fprintf(1, ['fitglm analysis for Justin''s paper\n\n'])
        
        no_files=max(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        
        szpc=size(percent_windows);
        for per_ii=1:szpc(1)
            for bwii=1:no_bandwidths
                for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
                    for elec=1:16
                        if sum(which_electrodes==elec)>0
                            
                            no_fitglms=no_fitglms+1;
                            fitglm_results(no_fitglms).per_ii=per_ii;
                            fitglm_results(no_fitglms).mouseNo=mouseNo;
                            fitglm_results(no_fitglms).elec=elec;
                            fitglm_results(no_fitglms).no_trials=0;
                            fitglm_results(no_fitglms).bwii=bwii;
                            fitglm_results(no_fitglms).elec=elec;
                            
                            %For each electrode compute the delta dB and
                            %input parameters for the fitglm
                            for fileNo=1:no_files
                                if sum(files==fileNo)>0
                                    if handles_drgb.drgbchoices.mouse_no(fileNo)==mouseNo
                                        lfpodNo_ref=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                                        
                                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref)))
                                            
                                            
                                            if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower))
                                                
                                                percent_mask=[];
                                                percent_mask=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower>=percent_windows(per_ii,1))...
                                                    &(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower<=percent_windows(per_ii,2));
                                                
                                                for evNo=1:length(eventType)
                                                    
                                                    trials_in_event_Ev=[];
                                                    trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(eventType(evNo),:)==1)&percent_mask;
                                                    
                                                    if (sum(trials_in_event_Ev)>=1)
                                                        
                                                        lfpodNo=find((files_per_lfp==fileNo)&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                                        
                                                        % Ev1
                                                        this_dB_powerref=zeros(sum(trials_in_event_Ev),length(frequency));
                                                        this_dB_powerref(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower(trials_in_event_Ev,:));
                                                        
                                                        
                                                        this_dB_power=zeros(sum(trials_in_event_Ev),length(frequency));
                                                        this_dB_power(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_event_Ev,:));
                                                        
                                                        
                                                        this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                                        
                                                        %Enter the  Ev1
                                                        this_delta_dB_powerEv=zeros(sum(trials_in_event_Ev),1);
                                                        this_delta_dB_powerEv=mean(this_dB_power(:,this_band)-this_dB_powerref(:,this_band),2);
                                                        
                                                        fitglm_results(no_fitglms).concs(fitglm_results(no_fitglms).no_trials+1:fitglm_results(no_fitglms).no_trials+sum(trials_in_event_Ev))=...
                                                            handles_pars.concs(evNo);
                                                        fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).no_trials+1:fitglm_results(no_fitglms).no_trials+sum(trials_in_event_Ev))=...
                                                            this_delta_dB_powerEv;
                                                        
                                                        fitglm_results(no_fitglms).percent(fitglm_results(no_fitglms).no_trials+1:fitglm_results(no_fitglms).no_trials+sum(trials_in_event_Ev))=...
                                                            handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower(trials_in_event_Ev);
                                                        
                                                        %Is this S= or S-?
                                                        if evNo<=length(eventType)/2
                                                            %This is S+
                                                            fitglm_results(no_fitglms).spm(fitglm_results(no_fitglms).no_trials+1:fitglm_results(no_fitglms).no_trials+sum(trials_in_event_Ev))=...
                                                                ones(1,sum(trials_in_event_Ev));
                                                        else
                                                            %This is S-
                                                            fitglm_results(no_fitglms).spm(fitglm_results(no_fitglms).no_trials+1:fitglm_results(no_fitglms).no_trials+sum(trials_in_event_Ev))=...
                                                                zeros(1,sum(trials_in_event_Ev));
                                                        end
                                                        
                                                        %Was the mouse rewarded
                                                        %in the last trial?
                                                        Hit_trials=[];
                                                        Hit_trials=handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(2,:)==1;
                                                        last_Hit_trials=zeros(1,length(Hit_trials));
                                                        last_Hit_trials(2:end)=Hit_trials(1:end-1);
                                                        fitglm_results(no_fitglms).last_rewarded(fitglm_results(no_fitglms).no_trials+1:fitglm_results(no_fitglms).no_trials+sum(trials_in_event_Ev))=...
                                                            last_Hit_trials(trials_in_event_Ev);
                                                        fitglm_results(no_fitglms).no_trials=fitglm_results(no_fitglms).no_trials+sum(trials_in_event_Ev);
                                                        
                                                        fprintf(1, ['%d trials in event No %d succesfully processed for file No %d electrode %d\n'],sum(trials_in_event_Ev), evNo,files(fileNo),elec);
                                                    end
                                                    
                                                end
                                            else
                                                
                                                fprintf(1, ['Empty allPower for file No %d electrode %d\n'],fileNo,elec);
                                                
                                                
                                            end
                                            
                                            
                                        else
                                            fprintf(1, ['Empty LFP reference for file No %d electrode %d\n'],fileNo,elec);
                                            
                                            
                                        end
                                        
                                        
                                    end
                                end
                            end %file No
                            
                            %For each electrode compute the fitglm
                            %Do fitglm
                            
                            
                            % Store the variables in a table.
                            tbl = table(fitglm_results(no_fitglms).delta_dB_Power',fitglm_results(no_fitglms).concs',fitglm_results(no_fitglms).spm',fitglm_results(no_fitglms).last_rewarded','VariableNames',{'delta_dB','log10_Odor_Concentration','spm','was_last_rewarded'});
                            
                            % Fit a linear regression model
                            fitglm_results(no_fitglms).lm = fitglm(tbl,'delta_dB~log10_Odor_Concentration+spm+was_last_rewarded');
                            
                            %Save the two p values for concentration and
                            %spm
                            fitglm_results(no_fitglms).p_val_conc=fitglm_results(no_fitglms).lm.Coefficients{2,4};
                            fitglm_results(no_fitglms).p_val_spm=fitglm_results(no_fitglms).lm.Coefficients{3,4};
                            fitglm_results(no_fitglms).p_val_last_reward=fitglm_results(no_fitglms).lm.Coefficients{4,4};
                            fitglm_results(no_fitglms).AIC(7)=fitglm_results(no_fitglms).lm.ModelCriterion.AIC;
                            fitglm_results(no_fitglms).BIC(7)=fitglm_results(no_fitglms).lm.ModelCriterion.BIC;
                            
                            % Fit a linear regression model for delta_dB~log10_Odor_Concentration+was_last_rewarded
                            lm = fitglm(tbl,'delta_dB~log10_Odor_Concentration+was_last_rewarded');
                            fitglm_results(no_fitglms).AIC(6)=lm.ModelCriterion.AIC;
                            fitglm_results(no_fitglms).BIC(6)=lm.ModelCriterion.BIC;
                            
                            % Fit a linear regression model for spm+was_last_rewarded
                            lm = fitglm(tbl,'delta_dB~spm+was_last_rewarded');
                            fitglm_results(no_fitglms).AIC(5)=lm.ModelCriterion.AIC;
                            fitglm_results(no_fitglms).BIC(5)=lm.ModelCriterion.BIC;
                            
                            % Fit a linear regression model for delta_dB~log10_Odor_Concentration+spm
                            lm = fitglm(tbl,'delta_dB~log10_Odor_Concentration+spm');
                            fitglm_results(no_fitglms).AIC(4)=lm.ModelCriterion.AIC;
                            fitglm_results(no_fitglms).BIC(4)=lm.ModelCriterion.BIC;
                            
                            % Fit a linear regression model for last
                            % rewarded
                            lm = fitglm(tbl,'delta_dB~was_last_rewarded');
                            fitglm_results(no_fitglms).AIC(3)=lm.ModelCriterion.AIC;
                            fitglm_results(no_fitglms).BIC(3)=lm.ModelCriterion.BIC;
                            
                            % Fit a linear regression model for spm
                            lm = fitglm(tbl,'delta_dB~spm');
                            fitglm_results(no_fitglms).AIC(2)=lm.ModelCriterion.AIC;
                            fitglm_results(no_fitglms).BIC(2)=lm.ModelCriterion.BIC;
                            
                            % Fit a linear regression model for conc
                            lm = fitglm(tbl,'delta_dB~log10_Odor_Concentration');
                            fitglm_results(no_fitglms).AIC(1)=lm.ModelCriterion.AIC;
                            fitglm_results(no_fitglms).BIC(1)=lm.ModelCriterion.BIC;
                            
                            %Plot dB power for each electrode
                            try
                                close(1)
                            catch
                            end
                            hFig1=figure(1);
                            set(hFig1, 'units','normalized','position',[.07 .1 .75 .3])
                            
                            subplot(1,3,1)
                            plot(fitglm_results(no_fitglms).concs',fitglm_results(no_fitglms).delta_dB_Power','ob')
                            minc=min(fitglm_results(no_fitglms).concs');
                            mindB=lm.Coefficients{1,1}+(lm.Coefficients{2,1})*minc;
                            maxc=max(fitglm_results(no_fitglms).concs');
                            maxdB=lm.Coefficients{1,1}+(lm.Coefficients{2,1})*maxc;
                            hold on
                            plot([minc maxc],[mindB maxdB],'-r')
                            xlabel('log(concentration)')
                            ylabel('dB power')
                            
                            subplot(1,3,2)
                            plot(fitglm_results(no_fitglms).spm',fitglm_results(no_fitglms).delta_dB_Power','ob')
                            hold on
                            %Splus bar
                            bar(1,mean(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).spm==1)))
                            %Confidence intervals
                            SEMsp = std(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).spm==1))...
                                /sqrt(length(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).spm==1)));               % Standard Error
                            tssp = tinv([0.05  0.95],sum(fitglm_results(no_fitglms).spm==1)-1);      % T-Score
                            CIsp = mean(mean(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).spm==1))) + tssp*SEMsp;
                            plot([1 1],CIsp,'-r','Linewidth',3)
                            %Sminus bar
                            bar(0,mean(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).spm==0)))% Confidence Intervals
                            %Confidence intervals
                            SEMsm = std(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).spm==0))...
                                /sqrt(length(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).spm==0)));               % Standard Error
                            tssm = tinv([0.05  0.95],sum(fitglm_results(no_fitglms).spm==0)-1);      % T-Score
                            CIsm = mean(mean(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).spm==0))) + tssp*SEMsp;
                            plot([0 0],CIsm,'-r','Linewidth',3)
                            xlabel('S= or S-')
                            ylabel('db power')
                            title(['no_fitglms ' num2str(no_fitglms)])
                            
                            %last_rewarded
                            subplot(1,3,3)
                            plot(fitglm_results(no_fitglms).last_rewarded',fitglm_results(no_fitglms).delta_dB_Power','ob')
                            hold on
                            %Not rewarded bar
                            bar(1,mean(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).last_rewarded==1)))
                            %Confidence intervals
                            SEMsp = std(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).last_rewarded==1))...
                                /sqrt(length(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).last_rewarded==1)));               % Standard Error
                            tssp = tinv([0.05  0.95],sum(fitglm_results(no_fitglms).last_rewarded==1)-1);      % T-Score
                            CIsp = mean(mean(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).last_rewarded==1))) + tssp*SEMsp;
                            plot([1 1],CIsp,'-r','Linewidth',3)
                            %Rewarded bar
                            bar(0,mean(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).last_rewarded==0)))% Confidence Intervals
                            %Confidence intervals
                            SEMsm = std(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).last_rewarded==0))...
                                /sqrt(length(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).last_rewarded==0)));               % Standard Error
                            tssm = tinv([0.05  0.95],sum(fitglm_results(no_fitglms).last_rewarded==0)-1);      % T-Score
                            CIsm = mean(mean(fitglm_results(no_fitglms).delta_dB_Power(fitglm_results(no_fitglms).last_rewarded==0))) + tssp*SEMsp;
                            plot([0 0],CIsm,'-r','Linewidth',3)
                            
                            
                            xlabel('Last rewarded')
                            ylabel('db power')
                            
                            pfft=1;
                            
                        end
                        
                    end
                end %mouseNo
                
            end %bandwidth
        end %end perii
        
        %Now parse through the p values
        for ii=1:no_fitglms
            periis(ii)=fitglm_results(ii).per_ii;
            bwiis(ii)=fitglm_results(ii).bwii;
            pvc(ii)=fitglm_results(ii).p_val_conc;
            pvspm(ii)=fitglm_results(ii).p_val_spm;
            pv_last_rewards(ii)=fitglm_results(ii).p_val_last_reward;
            for jj=1:7
                aic(ii,jj)=fitglm_results(ii).AIC(jj);
                bic(ii,jj)=fitglm_results(ii).BIC(jj);
            end
        end
        
        %Plot a bar graph showing the percent of significant estimated coefficents for the fitglm regression model
        try
            close(2)
        catch
        end
        hFig2=figure(2);
        hold on
        
        xpos=0;
        for perii=1:2
            for bwii=1:no_bandwidths
                %Concentration dependence
                percent_sig_conc=100*sum(pvc((periis==perii)&(bwiis==bwii))<=0.05)/sum((periis==perii)&(bwiis==bwii));
                fprintf(1, ['percent significant for dependence of delta dB on concentration for ' freq_names{bwii} ' ' per_lab{perii} '= %d, n= %d\n' ]...
                    ,percent_sig_conc,sum((periis==perii)&(bwiis==bwii)));
                xpos=xpos+1;
                bar(xpos,percent_sig_conc,'b')
                
                %spm dependence
                percent_sig_spm=100*sum(pvspm((periis==perii)&(bwiis==bwii))<=0.05)/sum((periis==perii)&(bwiis==bwii));
                fprintf(1, ['percent significant for dependence of delta dB on spm for ' freq_names{bwii} ' ' per_lab{perii} '= %d, n= %d\n' ]...
                    ,percent_sig_spm,sum((periis==perii)&(bwiis==bwii)));
                xpos=xpos+1;
                bar(xpos,percent_sig_spm,'r')
                
                %last reward
                percent_sig_last=100*sum(pv_last_rewards((periis==perii)&(bwiis==bwii))<=0.05)/sum((periis==perii)&(bwiis==bwii));
                fprintf(1, ['percent significant for dependence of delta dB on last_rewards for ' freq_names{bwii} ' ' per_lab{perii} '= %d, n= %d\n\n' ]...
                    ,percent_sig_last,sum((periis==perii)&(bwiis==bwii)));
                xpos=xpos+1;
                bar(xpos,percent_sig_last,'m')
                xpos=xpos+2;
            end
            xpos=xpos+4;
        end
        
        title('blue=concentration, red=spm, magenta=last reward, proficient left, naive right')
        ylabel('Percent significant regression coefficients')
        
        %Plot AICs
        figNo=2;
        for perii=1:2
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig2=figure(figNo);
            
            if perii==1
                title('AICs for naive mice')
            else
                title('AICs for proficient mice')
            end
            
            
            
            for bwii=1:no_bandwidths
                subplot(1,4,bwii)
                hold on
                for ii=1:7
                    [f_aic,x_aic] = drg_ecdf(aic((periis==perii)&(bwiis==bwii),ii));
                    plot(x_aic,f_aic,these_colors{ii})
                end
                
            end
            
            
        end
        
        
        %We should use aic() Akaike's Information Criterion for estimated
        %model with either conc or spm
        %https://en.wikipedia.org/wiki/Akaike_information_criterion
        %This is working well with fitglm, but should we try fitglm?
        %https://www.mathworks.com/help/stats/generalized-linear-regression.html
        %Should we make last rewarded 2 if the last two were rewarded?
        
        fprintf(1, '\n\n')
        
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix]);
        
        
    case 17
        %Justin's per mouse analysis for PAC generates plots of MI, vector
        %length and mean angle vs events and proficent vs naive
        
        
        mean_MI_No=0;
        mean_MI=[];
        mean_MI_perii=[];
        
        mean_MI_No_per_mouse=0;
        
        
        
        fprintf(1, ['PAC analysis for Justin''s paper\n\n'])
        p_vals=[];
        no_files=max(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        
        szpc=size(percent_windows);
        for per_ii=1:szpc(1)
            
            
            
            for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
                for elec=1:16
                    if sum(which_electrodes==elec)>0
                        
                        theseEvNos_thisMouse_thisElec=[];
                        for evNo=1:length(eventType)
                            for pacii=1:no_pacii
                                theseEvNos_thisMouse_thisElec(evNo,pacii).noEv=0;
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
                                            percent_mask=(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrLFPPower>=percent_windows(per_ii,1))...
                                                &(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrLFPPower<=percent_windows(per_ii,2));
                                            
                                            if ~isempty(percent_mask)
                                                
                                                for evNo=1:length(eventType)
                                                    
                                                    noWB_for_evNo(evNo)=-1;
                                                    
                                                    trials_in_event_Ev=[];
                                                    trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventLFPPower(eventType(evNo),:)==1)&percent_mask;
                                                    
                                                    if (sum(trials_in_event_Ev)>=1)
                                                        
                                                        %Do per bandwidth analysis
                                                        for pacii=1:no_pacii
                                                            
                                                            %Enter the modulation index
                                                            this_MI_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_MI_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).mod_indx(trials_in_event_Ev);
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii).this_MI_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+sum(trials_in_event_Ev))=this_MI_Ev;
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii).whichMouse(theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+sum(trials_in_event_Ev))=mouseNo*ones(1,length(this_MI_Ev));
                                                            
                                                            %Save per session
                                                            %value for MI
                                                            mean_MI_No=mean_MI_No+1;
                                                            mean_MI(mean_MI_No)=mean(this_MI_Ev);
                                                            mean_MI_perii(mean_MI_No)=per_ii;
                                                            mean_MI_evNo(mean_MI_No)=evNo;
                                                            mean_MI_pacii(mean_MI_No)=pacii;
                                                            mean_MI_fileNo(mean_MI_No)=fileNo;
                                                            
                                                            if mean_MI(mean_MI_No)>=0.035
                                                                fprintf(1, ['MI larger than 0.035 for mouse no %d, file no %d, electrode, %d, pac no %d, perii %d, conc, %d\n'],mouseNo, fileNo, elec, pacii, per_ii, evNo);
                                                            end
                                                            
                                                            %Enter the meanVectorLength
                                                            this_meanVectorLength_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_meanVectorLength_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).meanVectorLength(trials_in_event_Ev);
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii).this_meanVectorLength_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+sum(trials_in_event_Ev))=this_meanVectorLength_Ev;
                                                            
                                                            %Enter the meanVectorAngle
                                                            this_meanVectorAngle_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_meanVectorAngle_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).meanVectorAngle(trials_in_event_Ev);
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii).this_meanVectorAngle_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+sum(trials_in_event_Ev))=this_meanVectorAngle_Ev;
                                                            
                                                            %Enter the peakAngle
                                                            this_peakAngle_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_peakAngle_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).peakAngle(trials_in_event_Ev);
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii).this_peakAngle_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+sum(trials_in_event_Ev))=this_peakAngle_Ev;
                                                            
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii).noEv=theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+sum(trials_in_event_Ev);
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
                        
                        
                        %Calculate per mouse PAC measures
                        for evNo=1:length(eventType)
                            for pacii=1:no_pacii
                                
                                if theseEvNos_thisMouse_thisElec(evNo,pacii).noEv>0
                                    %Calculate per mouse MI
                                    mean_MI_No_per_mouse=mean_MI_No_per_mouse+1;
                                    this_mouse_MI=[];
                                    this_mouse_MI=theseEvNos_thisMouse_thisElec(evNo,pacii).this_MI_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).whichMouse==mouseNo);
                                    mean_MI_per_mouse(mean_MI_No_per_mouse)=mean(this_mouse_MI);
                                    
                                    mean_MI_perii_per_mouse(mean_MI_No_per_mouse)=per_ii;
                                    mean_MI_evNo_per_mouse(mean_MI_No_per_mouse)=evNo;
                                    mean_MI_pacii_per_mouse(mean_MI_No_per_mouse)=pacii;
                                    mean_MI_mouseNo_per_mouse(mean_MI_No_per_mouse)=mouseNo;
                                    mean_MI_electNo_per_mouse(mean_MI_No_per_mouse)=elec;
                                    
                                    %Calculate per mouse meanVectorLength
                                    this_mouse_meanVectorLength=[];
                                    this_mouse_meanVectorLength=theseEvNos_thisMouse_thisElec(evNo,pacii).this_meanVectorLength_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).whichMouse==mouseNo);
                                    mean_meanVectorLength_per_mouse(mean_MI_No_per_mouse)=mean(this_mouse_meanVectorLength);
                                    
                                    %Calculate per mouse meanVectorAngle
                                    this_mouse_meanVectorAngle=[];
                                    this_mouse_meanVectorAngle=theseEvNos_thisMouse_thisElec(evNo,pacii).this_meanVectorAngle_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).whichMouse==mouseNo);
                                    mean_meanVectorAngle_per_mouse(mean_MI_No_per_mouse)=(180/pi)*circ_axial(circ_mean(this_mouse_meanVectorAngle'*pi/180)');
                                    
                                    %Calculate per mouse peakAngle
                                    this_mouse_peakAngle=[];
                                    this_mouse_peakAngle=theseEvNos_thisMouse_thisElec(evNo,pacii).this_peakAngle_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).whichMouse==mouseNo);
                                    mean_peakAngle_per_mouse(mean_MI_No_per_mouse)=(180/pi)*circ_axial(circ_mean(this_mouse_peakAngle'*pi/180)');
                                end
                            end
                        end
                        
                    end
                    
                end
            end
            
            
        end
        fprintf(1, '\n\n')
        
        
        %Now plot the average MI per electrode per session (file)
        conc_anno_loc = {[0.15 0.6 0.2 0.2], [0.15 0.6 0.2 0.17], [0.15 0.6 0.2 0.14], [0.15 0.6 0.2 0.11], [0.15 0.6 0.2 0.08], [0.15 0.6 0.2 0.05]};
        fig_pos = {[664 550 576 513],[1233 550 576 513],[664 36 576 513],[1233 36 576 513]};
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            data_MI=[];
            prof_naive=[];
            events=[];
            
            %Plot the average
            try
                close(pacii)
            catch
            end
            figure(pacii)
            hold on
            
            for evNo=1:length(eventType)
                
                for per_ii=1:2      %performance bins. blue = naive, red = proficient
                    
                    bar_offset=21-evNo*3+(2-per_ii);
                    %bar_offset=(2-(per_ii-1))+3*(evNo-1);
                    if per_ii==1
                        bar(bar_offset,mean(mean_MI((mean_MI_perii==per_ii)&(mean_MI_pacii==pacii)&(mean_MI_evNo==evNo))),'r','LineWidth', 3)
                    else
                        bar(bar_offset,mean(mean_MI((mean_MI_perii==per_ii)&(mean_MI_pacii==pacii)&(mean_MI_evNo==evNo))),'b','LineWidth', 3)
                    end
                    plot(bar_offset,mean(mean_MI((mean_MI_perii==per_ii)&(mean_MI_pacii==pacii)&(mean_MI_evNo==evNo))),'ok','LineWidth', 3)
                    plot((bar_offset)*ones(1,sum((mean_MI_perii==per_ii)&(mean_MI_pacii==pacii)&(mean_MI_evNo==evNo))),mean_MI((mean_MI_perii==per_ii)&(mean_MI_pacii==pacii)&(mean_MI_evNo==evNo)),'o',...
                        'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                    CI = bootci(1000, {@mean, mean_MI((mean_MI_perii==per_ii)&(mean_MI_pacii==pacii)&(mean_MI_evNo==evNo))},'type','cper');
                    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                    
                    annotation('textbox',conc_anno_loc{per_ii},'String',per_lab(per_ii),'Color',these_colors{3-per_ii},'EdgeColor','none');
                    
                    
                    %Save data for anovan
                    data_MI=[data_MI mean_MI((mean_MI_perii==per_ii)&(mean_MI_pacii==pacii)&(mean_MI_evNo==evNo))];
                    prof_naive=[prof_naive per_ii*ones(1,sum((mean_MI_perii==per_ii)&(mean_MI_pacii==pacii)&(mean_MI_evNo==evNo)))];
                    events=[events evNo*ones(1,sum((mean_MI_perii==per_ii)&(mean_MI_pacii==pacii)&(mean_MI_evNo==evNo)))];
                end
            end
            title(['Average MI per electrode per session for PAC theta/' freq_names{pacii+1}])
            set(gcf,'OuterPosition',fig_pos{pacii});
            if sum(eventType==3)>0
                xticks([15,16,18,19])
                xticklabels({evTypeLabels{2},evTypeLabels{2},evTypeLabels{1},evTypeLabels{1}})
                
            else
                bar_lab_loc = [3.5 6.5 9.5 12.5 15.5 18.5];
                xticks(bar_lab_loc)
                xticklabels(concs2)
                xlabel('Concentration (%)')
            end
            ylabel('Modulation Index')
            
            
            %Calculate anovan for inteaction
            [p,tbl,stats]=anovan(data_MI,{prof_naive events},'model','interaction','varnames',{'proficient_vs_naive','within_vs_between'},'display','off');
            fprintf(1, ['p value for anovan MI per session for naive vs proficient for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan MI per session for events for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(2));
            fprintf(1, ['p value for anovan MI per session for events x naive-proficient for PAC theta/' freq_names{pacii+1} '= %d \n\n'],  p(3));
            p_anova_np_adj(pacii)=p(1);
            p_anova_wb_adj(pacii)=p(2);
        end
        
        %Now plot the average MI per electrode per mouse
        conc_anno_loc = {[0.15 0.6 0.2 0.2], [0.15 0.6 0.2 0.17], [0.15 0.6 0.2 0.14], [0.15 0.6 0.2 0.11], [0.15 0.6 0.2 0.08], [0.15 0.6 0.2 0.05]};
        fig_pos = {[664 550 576 513],[1233 550 576 513],[664 36 576 513],[1233 36 576 513]};
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            data_MI=[];
            prof_naive=[];
            events=[];
            
            %Plot the average
            try
                close(pacii+3)
            catch
            end
            figure(pacii+3)
            hold on
            
            for evNo=1:length(eventType)
                
                for per_ii=1:2      %performance bins. blue = naive, red = proficient
                    
                    bar_offset=21-evNo*3+(2-per_ii);
                    if per_ii==1
                        bar(bar_offset,mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))),'r','LineWidth', 3)
                    else
                        bar(bar_offset,mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))),'b','LineWidth', 3)
                    end
                    plot(bar_offset,mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))),'ok','LineWidth', 3)
                    plot((bar_offset)*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))),...
                        mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)),'o',...
                        'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                    
                    CI = bootci(1000, {@mean, mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))},'type','cper');
                    
                    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                    
                    annotation('textbox',conc_anno_loc{per_ii},'String',per_lab(per_ii),'Color',these_colors{3-per_ii},'EdgeColor','none');
                    
                    
                    %Save data for anovan
                    data_MI=[data_MI mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))];
                    prof_naive=[prof_naive per_ii*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)))];
                    events=[events evNo*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)))];
                end
            end
            title(['Average MI per electrode per mouse for PAC theta/' freq_names{pacii+1}])
            set(gcf,'OuterPosition',fig_pos{pacii});
            if sum(eventType==3)>0
                xticks([15,16,18,19])
                xticklabels({evTypeLabels{2},evTypeLabels{2},evTypeLabels{1},evTypeLabels{1}})
                
            else
                bar_lab_loc = [3.5 6.5 9.5 12.5 15.5 18.5];
                xticks(bar_lab_loc)
                xticklabels(concs2)
                xlabel('Concentration (%)')
            end
            ylabel('Modulation Index')
            
            
            %Calculate anovan for inteaction
            [p,tbl,stats]=anovan(data_MI,{prof_naive events},'model','interaction','varnames',{'proficient_vs_naive','within_vs_between'},'display','off');
            fprintf(1, ['p value for anovan MI per mouse for naive vs proficient for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan MI per mouse for events for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(2));
            fprintf(1, ['p value for anovan MI per mouse for events x naive-proficient for PAC theta/' freq_names{pacii+1} '= %d \n\n'],  p(3));
            p_anova_np_adj(pacii)=p(1);
            p_anova_wb_adj(pacii)=p(2);
        end
        
        %Now plot the cumulative histogram for MI per electrode per mouse
        
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            data_MI=[];
            prof_naive=[];
            events=[];
            
            %Plot the average
            try
                close(pacii+6)
            catch
            end
            figure(pacii+6)
            hold on
            
            %performance bins. blue = naive, red = proficient
            
            %Proficient S-
            per_ii=1;
            evNo=2;
            if length((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))>3
                [f_aic,x_aic] = drg_ecdf(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)));
                plot(x_aic,f_aic,'Color',[1 0.8 0.8])
            end
            
            %Proficient S+
            per_ii=1;
            evNo=1;
            if length((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))>3
                [f_aic,x_aic] = drg_ecdf(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)));
                plot(x_aic,f_aic,'r')
            end
            
            %Naive S-
            per_ii=2;
            evNo=2;
            if length((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))>3
                [f_aic,x_aic] = drg_ecdf(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)));
                plot(x_aic,f_aic,'Color',[0.8 0.8 1])
            end
            
            %Naive S+
            per_ii=2;
            evNo=1;
            if length((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))>3
                [f_aic,x_aic] = drg_ecdf(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)));
                plot(x_aic,f_aic,'b')
            end
            
            
            legend('S- proficient','S+ proficient','S- naive','S+ naive')
            
            title(['Cumulative histograms for average MI per electrode per mouse for PAC theta/' freq_names{pacii+1}])
            
            xlabel('Average MI')
            ylabel('Modulation Index')
            
            
            
        end
        
        %Plot the meanVectorAngle
        min_MI=0.000;
        
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            data_meanVectorAngle=[];
            prof_naive=[];
            events=[];
            
            try
                close(12+pacii)
            catch
            end
            figure(12+pacii)
            hold on
            
            shift_x=0.03;
            
            if length(eventType)>2
                %These are concentrations
                for per_ii=2:-1:1
                    
                    
                    
                    %Plot lines between individual points
                    
                    for mouseNo=1:max(mean_MI_mouseNo_per_mouse)
                        for elect=1:16
                            if per_ii==1
                                plot(log10(concs2)+shift_x,...
                                    mean_meanVectorAngle_per_mouse((mean_MI_electNo_per_mouse==elect)&(mean_MI_mouseNo_per_mouse==mouseNo)&(mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)),...
                                    '-','lineWidth',0.5,'Color',[1 0.7 0.7])
                            else
                                plot(log10(concs2)-shift_x,...
                                    mean_meanVectorAngle_per_mouse((mean_MI_electNo_per_mouse==elect)&(mean_MI_mouseNo_per_mouse==mouseNo)&(mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)),...
                                    '-','LineWidth',0.5,'Color',[0.7 0.7 1])
                            end
                        end
                    end
                    
                    
                    
                    
                    %Plot individual points
                    
                    %These are concentrations
                    for evNo=1:length(eventType)
                        if per_ii==1
                            plot((log10(concs2(evNo))+shift_x)*ones(1,sum((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))),...
                                mean_meanVectorAngle_per_mouse((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)),...
                                'o','MarkerFaceColor',[1 0.7 0.7],'MarkerEdgeColor',[1 0.7 0.7])
                        else
                            plot((log10(concs2(evNo))-shift_x)*ones(1,sum((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))),...
                                mean_meanVectorAngle_per_mouse((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)),...
                                'o','MarkerFaceColor',[1 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 1])
                        end
                    end
                    
                    
                end
                
                %                 for per_ii=2:-1:1
                %                     if per_ii==1
                %                         plot(log10(concs2)+shift_x,circ_meanVA,'-o','LineWidth',3,'MarkerEdgeColor','r','MarkerFaceColor','r','Color','r');
                %                     else
                %                         plot(log10(concs2)-shift_x,circ_meanVA,'-o','LineWidth',3,'MarkerEdgeColor','b','MarkerFaceColor','b','Color','b');
                %                     end
                %                 end
                
                legend('Proficient','Naive')
                title(['Vector Angle for PAC theta/' freq_names{pacii+1}])
                xlabel('log10(percent dilution)')
                ylabel('Angle')
            else
                
                %This is S+ S-
                for evNo=1:length(eventType)
                    
                    
                    %These are S+/S-
                    for mouseNo=1:max(mean_MI_mouseNo_per_mouse)
                        for elect=1:16
                            if (~isempty(mean_meanVectorAngle_per_mouse((mean_MI_electNo_per_mouse==elect)&(mean_MI_mouseNo_per_mouse==mouseNo)&(mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==1)&(mean_MI_pacii_per_mouse==pacii))))&...
                                    (~isempty(mean_meanVectorAngle_per_mouse((mean_MI_electNo_per_mouse==elect)&(mean_MI_mouseNo_per_mouse==mouseNo)&(mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==2)&(mean_MI_pacii_per_mouse==pacii))))
                                if evNo==2
                                    %S-
                                    plot([4 5],...
                                        [mean_meanVectorAngle_per_mouse((mean_MI_electNo_per_mouse==elect)&(mean_MI_mouseNo_per_mouse==mouseNo)&(mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==2)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))...
                                        mean_meanVectorAngle_per_mouse((mean_MI_electNo_per_mouse==elect)&(mean_MI_mouseNo_per_mouse==mouseNo)&(mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==1)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))],...
                                        '-','lineWidth',0.5,'Color',[0.7 0.7 0.7])
                                else
                                    %S+
                                    plot([1 2],...
                                        [mean_meanVectorAngle_per_mouse((mean_MI_electNo_per_mouse==elect)&(mean_MI_mouseNo_per_mouse==mouseNo)&(mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==2)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))...
                                        mean_meanVectorAngle_per_mouse((mean_MI_electNo_per_mouse==elect)&(mean_MI_mouseNo_per_mouse==mouseNo)&(mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==1)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))],...
                                        '-','lineWidth',0.5,'Color',[0.7 0.7 0.7])
                                end
                            end
                        end
                    end
                    
                end
                
                %Plot individual points
                for per_ii=1:2
                    for evNo=1:length(eventType)
                        if per_ii==1
                            plot((3*evNo-2+1)*ones(1,sum((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))),...
                                mean_meanVectorAngle_per_mouse((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)),...
                                'o','MarkerFaceColor',[1 0.7 0.7],'MarkerEdgeColor',[1 0.7 0.7])
                        else
                            plot((3*evNo-2)*ones(1,sum((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))),...
                                mean_meanVectorAngle_per_mouse((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)),...
                                'o','MarkerFaceColor',[0.7 0.7 1],'MarkerEdgeColor',[0.7 0.7 1])
                        end
                    end
                    
                end
                
                evNo=1;
                per_ii=2;
                circ_meanVA=[];
                CI=zeros(1,2);
                circ_meanVA=180*circ_axial(circ_mean(pi*mean_meanVectorAngle_per_mouse((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))'/180))/pi;
                CI(1,1:2) = 180*bootci(1000, {@circ_mean, pi*mean_meanVectorAngle_per_mouse((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))'/180},'type','cper')/pi;
                plot([1],circ_meanVA,'o','MarkerEdgeColor','b','MarkerFaceColor','b');
                plot([1 1],CI,'-b','LineWidth',3);
                
                evNo=1;
                per_ii=1;
                circ_meanVA=[];
                CI=zeros(1,2);
                circ_meanVA=180*circ_axial(circ_mean(pi*mean_meanVectorAngle_per_mouse((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))'/180))/pi;
                CI(1,1:2) = 180*bootci(1000, {@circ_mean, pi*mean_meanVectorAngle_per_mouse((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))'/180},'type','cper')/pi;
                plot([2],circ_meanVA,'o','MarkerEdgeColor','r','MarkerFaceColor','r');
                plot([2 2],CI,'-r','LineWidth',3);
                
                
                evNo=2;
                per_ii=2;
                circ_meanVA=[];
                CI=zeros(1,2);
                circ_meanVA=180*circ_axial(circ_mean(pi*mean_meanVectorAngle_per_mouse((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))'/180))/pi;
                CI(1,1:2) = 180*bootci(1000, {@circ_mean, pi*mean_meanVectorAngle_per_mouse((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))'/180},'type','cper')/pi;
                plot([4],circ_meanVA,'o','MarkerEdgeColor','b','MarkerFaceColor','b');
                plot([4 4],CI,'-b','LineWidth',3);
                
                evNo=2;
                per_ii=1;
                circ_meanVA=[];
                CI=zeros(1,2);
                circ_meanVA=180*circ_axial(circ_mean(pi*mean_meanVectorAngle_per_mouse((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))'/180))/pi;
                CI(1,1:2) = 180*bootci(1000, {@circ_mean, pi*mean_meanVectorAngle_per_mouse((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))'/180},'type','cper')/pi;
                plot([5],circ_meanVA,'o','MarkerEdgeColor','r','MarkerFaceColor','r');
                plot([5 5],CI,'-r','LineWidth',3);
                
                pffft=1;
                xlim([0 6])
                xticks([1 2 4 5])
                xticklabels({evTypeLabels{1},evTypeLabels{1},evTypeLabels{2},evTypeLabels{2}})
                legend('Proficient','Naive')
                title(['Vector Angle for PAC theta/' freq_names{pacii+1}])
                ylabel('Angle')
            end
            
            for per_ii=2:-1:1
                for evNo=1:length(eventType)
                    circ_meanVA(evNo)=180*circ_axial(circ_mean(pi*mean_meanVectorAngle_per_mouse((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))'/180))/pi;
                    CI(1:2,evNo) = 180*circ_axial(bootci(1000, {@circ_mean, pi*mean_meanVectorAngle_per_mouse((mean_MI_per_mouse>min_MI)&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))'/180},'type','cper'))/pi;
                    
                    %Save data for anovan
                    data_meanVectorAngle=[data_meanVectorAngle mean_meanVectorAngle_per_mouse((mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo))];
                    prof_naive=[prof_naive per_ii*ones(1,sum((mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)))];
                    events=[events evNo*ones(1,sum((mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)))];
                end
            end
            
            
            
            %             %Calculate anovan for inteaction
            %             [p,tbl,stats]=anovan(data_MI,{prof_naive events},'model','interaction','varnames',{'proficient_vs_naive','within_vs_between'},'display','off');
            %             fprintf(1, ['p value for anovan vector angle per mouse for naive vs proficient for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(1));
            %             fprintf(1, ['p value for anovan vector angle per mouse for events for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(2));
            %             fprintf(1, ['p value for anovan vector angle per mouse for events x naive-proficient for PAC theta/' freq_names{pacii+1} '= %d \n\n'],  p(3));
            
        end
        
        
        
        pfft=1
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix]);
        
    case 18
        %Justin's per mouse analysis for PAC
        %NOTE: This does the analysis in all the files and DOES not distinguish between groups!!!
        
        mean_MI_No=0;
        mean_MI=[];
        mean_MI_perii=[];
        mean_MI_No_per_mouse=0;
        
        
        
        fprintf(1, ['Pairwise auROC analysis for Fig 1 of Justin''s paper\n\n'])
        p_vals=[];
        no_files=max(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        %Initialize ROCs
        no_ROCs=0;
        ROCout=[];
        ROCelec=[];
        ROCpacii=[];
        ROCper_ii=[];
        ROCEvNo1=[];
        ROCEvNo2=[];
        ROC_between=[];
        ROC_neighbor=[];;
        auROC=[];
        p_valROC=[];
        p_vals_ROC=[];
        
        szpc=size(percent_windows);
        for per_ii=1:szpc(1)
            
            for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
                for elec=1:16
                    if sum(which_electrodes==elec)>0
                        
                        theseEvNos_thisMouse_thisElec=[];
                        for evNo=1:length(eventType)
                            for pacii=1:no_pacii
                                theseEvNos_thisMouse_thisElec(evNo,pacii).noEv=0;
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
                                            percent_mask=(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrLFPPower>=percent_windows(per_ii,1))...
                                                &(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrLFPPower<=percent_windows(per_ii,2));
                                            
                                            if ~isempty(percent_mask)
                                                for evNo=1:length(eventType)
                                                    
                                                    noWB_for_evNo(evNo)=-1;
                                                    
                                                    trials_in_event_Ev=[];
                                                    trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventLFPPower(eventType(evNo),:)==1)&percent_mask;
                                                    
                                                    if (sum(trials_in_event_Ev)>=1)
                                                        
                                                        %Do per bandwidth analysis
                                                        for pacii=1:no_pacii
                                                            
                                                            %Enter the modulation index
                                                            this_MI_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_MI_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).mod_indx(trials_in_event_Ev);
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii).this_MI_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+sum(trials_in_event_Ev))=this_MI_Ev;
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii).whichMouse(theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+sum(trials_in_event_Ev))=mouseNo*ones(1,length(this_MI_Ev));
                                                            
                                                            %Save per session
                                                            %value for MI
                                                            mean_MI_No=mean_MI_No+1;
                                                            mean_MI(mean_MI_No)=mean(this_MI_Ev);
                                                            mean_MI_perii(mean_MI_No)=per_ii;
                                                            mean_MI_evNo(mean_MI_No)=evNo;
                                                            mean_MI_pacii(mean_MI_No)=pacii;
                                                            mean_MI_fileNo(mean_MI_No)=fileNo;
                                                            
                                                            if mean_MI(mean_MI_No)>=0.035
                                                                fprintf(1, ['MI larger than 0.035 for mouse no %d, file no %d, electrode, %d, pac no %d, perii %d, conc, %d\n'],mouseNo, fileNo, elec, pacii, per_ii, evNo);
                                                            end
                                                            
                                                            %Enter the meanVectorLength
                                                            this_meanVectorLength_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_meanVectorLength_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).meanVectorLength(trials_in_event_Ev);
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii).this_meanVectorLength_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+sum(trials_in_event_Ev))=this_meanVectorLength_Ev;
                                                            
                                                            %Enter the meanVectorAngle
                                                            this_meanVectorAngle_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_meanVectorAngle_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).meanVectorAngle(trials_in_event_Ev);
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii).this_meanVectorAngle_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+sum(trials_in_event_Ev))=this_meanVectorAngle_Ev;
                                                            
                                                            %Enter the peakAngle
                                                            this_peakAngle_Ev=zeros(sum(trials_in_event_Ev),1);
                                                            this_peakAngle_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).peakAngle(trials_in_event_Ev);
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii).this_peakAngle_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+1:theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+sum(trials_in_event_Ev))=this_peakAngle_Ev;
                                                            
                                                            theseEvNos_thisMouse_thisElec(evNo,pacii).noEv=theseEvNos_thisMouse_thisElec(evNo,pacii).noEv+sum(trials_in_event_Ev);
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
                        
                        %Calculate ROC for this electrode and this percent
                        %band
                        for evNo1=1:length(eventType)
                            for evNo2=evNo1+1:length(eventType)
                                
                                
                                for pacii=1:no_pacii
                                    if (theseEvNos_thisMouse_thisElec(evNo1,pacii).noEv>5)&(theseEvNos_thisMouse_thisElec(evNo2,pacii).noEv>5)
                                        %Enter Ev1
                                        trials_in_event_Ev1=theseEvNos_thisMouse_thisElec(evNo1,pacii).noEv;
                                        this_meanAngle_Ev1=zeros(trials_in_event_Ev1,1);
                                        this_meanAngle_Ev1(:,1)=theseEvNos_thisMouse_thisElec(evNo1,pacii).this_meanVectorAngle_Ev;
                                        
                                        %Enter Ev2
                                        trials_in_event_Ev2=theseEvNos_thisMouse_thisElec(evNo2,pacii).noEv;
                                        total_trials=trials_in_event_Ev1+trials_in_event_Ev2;
                                        this_meanAngle_Ev2=zeros(trials_in_event_Ev2,1);
                                        this_meanAngle_Ev2(:,1)=theseEvNos_thisMouse_thisElec(evNo2,pacii).this_meanVectorAngle_Ev;
                                        
                                        
                                        %One has to find the best
                                        %start for the circle, I do
                                        %not know if this is the
                                        %best
                                        %In the future write a ROC that
                                        %goes around the circle using
                                        %parfor
                                        this_mean_angle=(180/pi)*circ_axial(circ_mean([this_meanAngle_Ev1' this_meanAngle_Ev2']'*pi/180)');
                                        
                                        if this_mean_angle<180
                                            this_mean_angle=this_mean_angle+180;
                                        else
                                            this_mean_angle=this_mean_angle-180;
                                        end
                                        
                                        this_meanAngle_Ev1=this_meanAngle_Ev1-this_mean_angle;
                                        this_meanAngle_Ev1(this_meanAngle_Ev1<0)=360+this_meanAngle_Ev1(this_meanAngle_Ev1<0);
                                        this_meanAngle_Ev2=this_meanAngle_Ev2-this_mean_angle;
                                        this_meanAngle_Ev2(this_meanAngle_Ev2<0)=360+this_meanAngle_Ev2(this_meanAngle_Ev2<0);
                                        
                                        roc_data=[];
                                        roc_data(1:sum(trials_in_event_Ev1),1)=this_meanAngle_Ev1;
                                        roc_data(1:sum(trials_in_event_Ev1),2)=zeros(sum(trials_in_event_Ev1),1);
                                        
                                        
                                        roc_data(sum(trials_in_event_Ev1)+1:total_trials,1)=this_meanAngle_Ev2;
                                        roc_data(sum(trials_in_event_Ev1)+1:total_trials,2)=ones(sum(trials_in_event_Ev2),1);
                                        
                                        
                                        %Find  ROC
                                        
                                        no_ROCs=no_ROCs+1;
                                        roc=roc_calc(roc_data,0,0.05,0);
                                        ROCout(no_ROCs).mouseNo=mouseNo;
                                        ROCelec(no_ROCs)=elec;
                                        ROCpacii(no_ROCs)=pacii;
                                        ROCper_ii(no_ROCs)=per_ii;
                                        ROCEvNo1(no_ROCs)=evNo1;
                                        ROCEvNo2(no_ROCs)=evNo2;
                                        if ((abs(evNo1-2)<=1)&(abs(evNo2-5)<=1))||((abs(evNo1-5)<=1)&(abs(evNo2-2)<=1))
                                            ROC_between(no_ROCs)=1;
                                        else
                                            ROC_between(no_ROCs)=0;
                                        end
                                        ROC_neighbor(no_ROCs)=abs(evNo1-evNo2);
                                        auROC(no_ROCs)=roc.AUC-0.5;
                                        p_valROC(no_ROCs)=roc.p;
                                        p_vals_ROC=[p_vals_ROC roc.p];
                                        
                                        
                                    end
                                end
                                
                            end
                        end
                        
                        if theseEvNos_thisMouse_thisElec(evNo,pacii).noEv>0
                            %Calculate per mouse PAC measures
                            for evNo=1:length(eventType)
                                for pacii=1:no_pacii
                                    
                                    %Calculate per mouse MI
                                    mean_MI_No_per_mouse=mean_MI_No_per_mouse+1;
                                    this_mouse_MI=[];
                                    this_mouse_MI=theseEvNos_thisMouse_thisElec(evNo,pacii).this_MI_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).whichMouse==mouseNo);
                                    mean_MI_per_mouse(mean_MI_No_per_mouse)=mean(this_mouse_MI);
                                    
                                    mean_MI_perii_per_mouse(mean_MI_No_per_mouse)=per_ii;
                                    mean_MI_evNo_per_mouse(mean_MI_No_per_mouse)=evNo;
                                    mean_MI_pacii_per_mouse(mean_MI_No_per_mouse)=pacii;
                                    mean_MI_mouseNo_per_mouse(mean_MI_No_per_mouse)=mouseNo;
                                    mean_MI_electNo_per_mouse(mean_MI_No_per_mouse)=elec;
                                    
                                    %Calculate per mouse meanVectorLength
                                    this_mouse_meanVectorLength=[];
                                    this_mouse_meanVectorLength=theseEvNos_thisMouse_thisElec(evNo,pacii).this_meanVectorLength_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).whichMouse==mouseNo);
                                    mean_meanVectorLength_per_mouse(mean_MI_No_per_mouse)=mean(this_mouse_meanVectorLength);
                                    
                                    %Calculate per mouse meanVectorAngle
                                    this_mouse_meanVectorAngle=[];
                                    this_mouse_meanVectorAngle=theseEvNos_thisMouse_thisElec(evNo,pacii).this_meanVectorAngle_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).whichMouse==mouseNo);
                                    if isempty(this_mouse_meanVectorAngle)
                                        mean_meanVectorAngle_per_mouse(mean_MI_No_per_mouse)=NaN;
                                    else
                                        mean_meanVectorAngle_per_mouse(mean_MI_No_per_mouse)=(180/pi)*circ_axial(circ_mean(this_mouse_meanVectorAngle'*pi/180)');
                                    end
                                    
                                    %Calculate per mouse peakAngle
                                    this_mouse_peakAngle=[];
                                    this_mouse_peakAngle=theseEvNos_thisMouse_thisElec(evNo,pacii).this_peakAngle_Ev(theseEvNos_thisMouse_thisElec(evNo,pacii).whichMouse==mouseNo);
                                    if isempty(this_mouse_peakAngle)
                                        mean_peakAngle_per_mouse(mean_MI_No_per_mouse)=NaN;
                                    else
                                        mean_peakAngle_per_mouse(mean_MI_No_per_mouse)=(180/pi)*circ_axial(circ_mean(this_mouse_peakAngle'*pi/180)');
                                    end
                                    
                                end
                            end
                            
                        end
                    end
                    
                end
            end
            
            
        end
        fprintf(1, '\n\n')
        
        pFDRauROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for auROC  = %d\n\n'],pFDRauROC);
        
        
        %Plot cumulative histos for auROCs within vs between S+ and S-
        %Using all epochs
        figNo = get(gcf,'Number');
        for pacii=1:no_pacii
            
            
            
            try
                close(figNo)
            catch
            end
            figure(figNo)
            figNo=figNo+1;
            
            hold on
            
            %Naive between
            if length(auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1)))>3
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1)));
                plot(x_aic,f_aic,'b')
            end
            
            %Proficient between
            if length(auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1)))>3
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1)));
                plot(x_aic,f_aic,'r')
            end
            
            %Naive within
            if length(auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0)))>3
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0)));
                plot(x_aic,f_aic,'Color',[0.8 0.8 1])
            end
            
            %Proficient within
            if length(auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0)))>3
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0)));
                plot(x_aic,f_aic,'Color',[1 0.8 0.8])
            end
            
            legend('Naive between','Proficient between','Naive within','Proficient within')
            xlabel('auROC')
            ylabel('Cumulative probability')
            if length(eventType)>2
                title(['Theta/' freq_names{pacii+1} ' PAC mean angle ROC, all concentrations'])
            else
                title(['Theta/' freq_names{pacii+1} ' PAC mean angle ROC, S+ vs S-'])
            end
            
            %Save the data for anovan
            data_auROC=[];
            prof_naive=[];
            between=[];
            
            %Naive between
            data_auROC=[data_auROC auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1))];
            prof_naive=[prof_naive zeros(1,sum((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1)))];
            between=[between ones(1,sum((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1)))];
            
            %Proficient between
            data_auROC=[data_auROC auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1))];
            prof_naive=[prof_naive ones(1,sum((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1)))];
            between=[between ones(1,sum((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1)))];
            
            %Naive within
            data_auROC=[data_auROC auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0))];
            prof_naive=[prof_naive zeros(1,sum((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0)))];
            between=[between zeros(1,sum((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0)))];
            
            %Proficient within
            data_auROC=[data_auROC auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0))];
            prof_naive=[prof_naive ones(1,sum((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0)))];
            between=[between zeros(1,sum((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0)))];
            
            %Calculate anovan for inteaction
            [p,tbl,stats]=anovan(data_auROC,{prof_naive between},'model','interaction','varnames',{'proficient_vs_naive','within_vs_between'},'display','off');
            fprintf(1, ['p value for anovan mean angle auROC for naive vs proficient for Theta/' freq_names{pacii+1} ' PAC= %d \n'],  p(1));
            fprintf(1, ['p value for anovan mean angle auROC for within vs between for for Theta/' freq_names{pacii+1} ' PAC= %d \n'],  p(2));
        end
        
        if sum(ROC_between)>0
            
            %Plot cumulative histos for auROCs within vs between S+ and S- using
            %only ROCs for adjacent concentrations
            for pacii=1:no_pacii
                
                
                
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                figNo=figNo+1;
                
                hold on
                
                %Naive between
                if ~isempty(auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==1)))
                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==1)));
                    plot(x_aic,f_aic,'b')
                end
                
                %Proficient between
                if ~isempty(auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==1)))
                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==1)));
                    plot(x_aic,f_aic,'r')
                end
                
                %Naive within
                if ~isempty(auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==1)))
                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==1)));
                    plot(x_aic,f_aic,'Color',[0.8 0.8 1])
                end
                
                %Proficient within
                if ~isempty(auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==1)))
                    [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==1)));
                    plot(x_aic,f_aic,'Color',[1 0.8 0.8])
                end
                
                legend('Naive between','Proficient between','Naive within','Proficient within')
                xlabel('auROC')
                ylabel('Cumulative probability')
                title(['Theta/' freq_names{pacii+1} ' PAC mean angle ROC, Adjacent concentrations'])
                
                %Save the data for anovan for adjacent ROCs
                data_auROC=[];
                prof_naive=[];
                between=[];
                
                %Naive between
                data_auROC=[data_auROC auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==1))];
                prof_naive=[prof_naive zeros(1,sum((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==1)))];
                between=[between ones(1,sum((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==1)))];
                
                %Proficient between
                data_auROC=[data_auROC auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1&(ROC_neighbor==1)))];
                prof_naive=[prof_naive ones(1,sum((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==1)))];
                between=[between ones(1,sum((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==1)))];
                
                %Naive within
                data_auROC=[data_auROC auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==1))];
                prof_naive=[prof_naive zeros(1,sum((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==1)))];
                between=[between zeros(1,sum((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==1)))];
                
                %Proficient within
                data_auROC=[data_auROC auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==1))];
                prof_naive=[prof_naive ones(1,sum((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==1)))];
                between=[between zeros(1,sum((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==1)))];
                
                %Calculate anovan for inteaction
                [p,tbl,stats]=anovan(data_auROC,{prof_naive between},'model','interaction','varnames',{'proficient_vs_naive','within_vs_between'},'display','off');
                fprintf(1, ['p value for anovan mean angle auROC adjacent concentrations for naive vs proficient for Theta/' freq_names{pacii+1} ' PAC= %d \n'],  p(1));
                fprintf(1, ['p value for anovan mean angle auROC adjacent concentrations  for within vs between for Theta/' freq_names{pacii+1} ' PAC= %d \n'],  p(2));
                
            end
            
            %Plot cumulative histos for auROCs within vs between S+ and S- using
            %only ROCs for concentrations separated by two log steps
            for pacii=1:no_pacii
                
                
                
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                figNo = figNo+1;
                
                hold on
                
                %Naive between
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==2)));
                plot(x_aic,f_aic,'b')
                
                %Proficient between
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==2)));
                plot(x_aic,f_aic,'r')
                
                %Naive within
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==2)));
                plot(x_aic,f_aic,'Color',[0.8 0.8 1])
                
                %Proficient within
                [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==2)));
                plot(x_aic,f_aic,'Color',[1 0.8 0.8])
                
                legend('Naive between','Proficient between','Naive within','Proficient within')
                xlabel('auROC')
                ylabel('Cumulative probability')
                title(['Theta/' freq_names{pacii+1} ' PAC mean angle ROC, concentrations separated by two log steps'])
                
                %Save the data for anovan for adjacent ROCs
                data_auROC=[];
                prof_naive=[];
                between=[];
                
                %Naive between
                data_auROC=[data_auROC auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==2))];
                prof_naive=[prof_naive zeros(1,sum((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==2)))];
                between=[between ones(1,sum((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==2)))];
                
                %Proficient between
                data_auROC=[data_auROC auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1&(ROC_neighbor==2)))];
                prof_naive=[prof_naive ones(1,sum((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==2)))];
                between=[between ones(1,sum((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==2)))];
                
                %Naive within
                data_auROC=[data_auROC auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==2))];
                prof_naive=[prof_naive zeros(1,sum((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==2)))];
                between=[between zeros(1,sum((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==2)))];
                
                %Proficient within
                data_auROC=[data_auROC auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==2))];
                prof_naive=[prof_naive ones(1,sum((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==2)))];
                between=[between zeros(1,sum((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==2)))];
                
                %Calculate anovan for inteaction
                [p,tbl,stats]=anovan(data_auROC,{prof_naive between},'model','interaction','varnames',{'proficient_vs_naive','within_vs_between'},'display','off');
                fprintf(1, ['p value for anovan mean angle auROC two log step concentrations for naive vs proficient for Theta/' freq_names{pacii+1} ' PAC= %d \n'],  p(1));
                fprintf(1, ['p value for anovan mean angle auROC two log step concentrations  for within vs between for Theta/' freq_names{pacii+1} ' PAC= %d \n'],  p(2));
            end
        end
        
        %Plot percent significant ROC
        no_within=zeros(2,4);
        no_sig_within=zeros(2,4);
        no_between=zeros(2,4);
        no_sig_between=zeros(2,4);
        no_within1=zeros(2,4);
        no_sig_within1=zeros(2,4);
        no_between1=zeros(2,4);
        no_sig_between1=zeros(2,4);
        no_within2=zeros(2,4);
        no_sig_within2=zeros(2,4);
        no_between2=zeros(2,4);
        no_sig_between2=zeros(2,4);
        
        for pacii=1:no_pacii
            for pcii=1:szpc(1)
                no_pairs=0;
                EvNo1=[];
                EvNo2=[];
                per_sig=zeros(length(eventType),length(eventType));
                for evNo1=1:length(eventType)
                    for evNo2=evNo1+1:length(eventType)
                        
                        
                        no_pairs=no_pairs+1;
                        these_ROCs=(ROCpacii==pacii)&(ROCEvNo1==evNo1)&(ROCEvNo2==evNo2)&(ROCper_ii==pcii);
                        sig(no_pairs)=sum((p_valROC<=pFDRauROC)&these_ROCs);
                        not_sig(no_pairs)=sum(these_ROCs)-sum((p_valROC<=pFDRauROC)&these_ROCs);
                        EvNo1(no_pairs)=evNo1;
                        EvNo2(no_pairs)=evNo2;
                        per_sig(evNo1,evNo2)=100*sig(no_pairs)/(sig(no_pairs)+not_sig(no_pairs));
                        
                        
                        if ((evNo1==1)&(evNo2==2))||((evNo1==2)&(evNo2==1))||((evNo1==1)&(evNo2==3))||((evNo1==3)&(evNo2==1))...
                                ||((evNo1==2)&(evNo2==3))||((evNo1==3)&(evNo2==2))||((evNo1==4)&(evNo2==5))||((evNo1==5)&(evNo2==4))...
                                ||((evNo1==4)&(evNo2==6))||((evNo1==6)&(evNo2==4))||((evNo1==5)&(evNo2==6))||((evNo1==6)&(evNo2==5))
                            %This is within
                            no_sig_within(pcii,pacii)=no_sig_within(pcii,pacii)+sum((p_valROC<=pFDRauROC)&these_ROCs);
                            no_within(pcii,pacii)=no_within(pcii,pacii)+sum(these_ROCs);
                            if abs(evNo1-evNo2)==1
                                no_sig_within1(pcii,pacii)=no_sig_within1(pcii,pacii)+sum((p_valROC<=pFDRauROC)&these_ROCs);
                                no_within1(pcii,pacii)=no_within1(pcii,pacii)+sum(these_ROCs);
                            end
                            if abs(evNo1-evNo2)==2
                                no_sig_within2(pcii,pacii)=no_sig_within2(pcii,pacii)+sum((p_valROC<=pFDRauROC)&these_ROCs);
                                no_within2(pcii,pacii)=no_within2(pcii,pacii)+sum(these_ROCs);
                            end
                        else
                            %This is bewteen
                            no_sig_between(pcii,pacii)=no_sig_between(pcii,pacii)+sum((p_valROC<=pFDRauROC)&these_ROCs);
                            no_between(pcii,pacii)=no_between(pcii,pacii)+sum(these_ROCs);
                            if abs(evNo1-evNo2)==1
                                no_sig_between1(pcii,pacii)=no_sig_between1(pcii,pacii)+sum((p_valROC<=pFDRauROC)&these_ROCs);
                                no_between1(pcii,pacii)=no_between1(pcii,pacii)+sum(these_ROCs);
                            end
                            if abs(evNo1-evNo2)==2
                                no_sig_between2(pcii,pacii)=no_sig_between2(pcii,pacii)+sum((p_valROC<=pFDRauROC)&these_ROCs);
                                no_between2(pcii,pacii)=no_between2(pcii,pacii)+sum(these_ROCs);
                            end
                        end
                        
                    end
                end
                
                for evNo1=1:length(eventType)
                    for evNo2=evNo1+1:length(eventType)
                        
                        if per_sig(evNo1,evNo2)==0
                            per_sig(evNo1,evNo2)=100/64;
                        end
                    end
                end
                
                %Plot the pseudocolor for percent significant auROCs
                
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                figNo = figNo+1;
                
                evNos_for_1=1:length(eventType);
                evNos_for_2=[1:length(eventType)]';
                
                drg_pcolor(repmat(evNos_for_1,length(eventType),1),repmat(evNos_for_2,1,length(eventType)),per_sig)
                cmjet=colormap(jet);
                cmjet(1,1)=0.7;
                cmjet(1,2)=0.7;
                cmjet(1,3)=0.7;
                colormap(cmjet)
                
                caxis([0 100]);
                
                hold on
                plot([4 4],[1 4],'-w','LineWidth', 5)
                plot([4 7],[4 4],'-w','LineWidth', 5)
                
                ax=gca;
                set(ax,'XTickLabel','')
                ylabel('dB')
                xticks([1.5:1:length(eventType)+1])
                xticklabels(handles_pars.concs2)
                yticks([1.5:1:length(eventType)+1])
                yticklabels(handles_pars.concs2)
                
                title(['Percent mean angle significant auROC  for Theta/' freq_names{pacii+1} ' PAC, ' per_lab(pcii)])
                set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                
                pfft=1;
            end
        end
        
        
        
        if sum(ROC_between)>0
            
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            
            set(hFig, 'units','normalized','position',[.83 .1 .05 .3])
            figNo = figNo+1;
            
            prain=[0:100/99:100];
            drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
            colormap jet
            shading interp
            caxis([0 100]);
            ax=gca;
            set(ax,'XTickLabel','')
            
            %Plot percent significant ROCs
            
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            figNo = figNo+1;
            
            for pacii=1:no_pacii
                subplot(2,2,pacii)
                hold on
                %Plot within naive
                %             bar(1,100*no_sig_within(2,pacii)/no_within(2,pacii),'b')
                %             bar(2,100*no_sig_within1(2,pacii)/no_within1(2,pacii),'b')
                bar(3,100*no_sig_within2(2,pacii)/no_within2(2,pacii),'b')
                
                %Plot within proficient
                %             bar(4,100*no_sig_within(1,pacii)/no_within(1,pacii),'r')
                %             bar(5,100*no_sig_within1(1,pacii)/no_within1(1,pacii),'r')
                bar(6,100*no_sig_within2(1,pacii)/no_within2(1,pacii),'r')
                
                %             %Plot between naive
                %             bar(9,100*no_sig_between(2,pacii)/no_between(2,pacii),'b')
                %             bar(10,100*no_sig_between1(2,pacii)/no_between1(2,pacii),'b')
                bar(11,100*no_sig_between2(2,pacii)/no_between2(2,pacii),'b')
                
                %Plot between proficient
                %             bar(12,100*no_sig_between(1,pacii)/no_between(1,pacii),'r')
                %             bar(13,100*no_sig_between1(1,pacii)/no_between1(1,pacii),'r')
                bar(14,100*no_sig_between2(1,pacii)/no_between2(1,pacii),'r')
                
                title(['Percent significant auROC ' freq_names{pacii+1} ])
            end
            
        else
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            
            set(hFig, 'units','normalized','position',[.83 .1 .05 .3])
            
            prain=[0:100/99:100];
            drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
            colormap jet
            shading interp
            caxis([0 100]);
            ax=gca;
            set(ax,'XTickLabel','')
            
            %Plot percent significant ROCs
            
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            figNo = figNo+1;
            
            for pacii=1:no_pacii
                subplot(2,2,pacii)
                hold on
                %Plot within naive
                bar(1,100*no_sig_within(2,pacii)/no_within(2,pacii),'b')
                %             bar(2,100*no_sig_within1(2,pacii)/no_within1(2,pacii),'b')
                %           bar(3,100*no_sig_within2(2,pacii)/no_within2(2,pacii),'b')
                
                %Plot within proficient
                bar(2,100*no_sig_within(1,pacii)/no_within(1,pacii),'r')
                %             bar(5,100*no_sig_within1(1,pacii)/no_within1(1,pacii),'r')
                %           bar(6,100*no_sig_within2(1,pacii)/no_within2(1,pacii),'r')
                
                %             %Plot between naive
                %             bar(9,100*no_sig_between(2,pacii)/no_between(2,pacii),'b')
                %             bar(10,100*no_sig_between1(2,pacii)/no_between1(2,pacii),'b')
                %bar(11,100*no_sig_between2(2,pacii)/no_between2(2,pacii),'b')
                
                %Plot between proficient
                %             bar(12,100*no_sig_between(1,pacii)/no_between(1,pacii),'r')
                %             bar(13,100*no_sig_between1(1,pacii)/no_between1(1,pacii),'r')
                %                 bar(14,100*no_sig_between2(1,pacii)/no_between2(1,pacii),'r')
                
                title(['Percent significant auROC ' freq_names{pacii+1} ])
                legend('Naive','Proficient')
            end
            
        end
        
        pffft=1;
        %         save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix]);
        
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
        
        
        %Now plot the average MI for each electrode calculated per mouse
        %(including all sessions for each mouse)
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            data_MI=[];
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
            
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
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
                        
                        if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
                            
                            include_group=1;
                            
                            if per_ii==1
                                bar(bar_offset,mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),'r','LineWidth', 3,'EdgeColor','none')
                            else
                                bar(bar_offset,mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),'b','LineWidth', 3,'EdgeColor','none')
                            end
                            
                            
                            plot(bar_offset,mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
                            plot((bar_offset)*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))),...
                                mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)),'o',...
                                'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            
                            if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>=2
                                CI = bootci(1000, {@mean, mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))},'type','cper');
                                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            end
                            
                            %Save data for anovan
                            data_MI=[data_MI mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))];
                            prof_naive=[prof_naive per_ii*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))];
                            events=[events evNo*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))];
                            groups=[groups grNo*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))];
                            mice=[mice mean_MI_mouseNo_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))];
                            
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
            
            title(['Average MI for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1}])
            
            
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
            
            ylabel('Modulation Index')
            
            
            %Calculate anovan for inteaction
            fprintf(1, ['ANOVAN for average MI for each electrode calculated per mouse with mouse as random factor PAC theta/' freq_names{pacii+1} '\n'])
            [p,tbl,stats]=anovan(data_MI,{prof_naive, events, groups, mice},'varnames',{'proficient_vs_naive','events','groups','mice'},'display','off','random',4);
            fprintf(1, ['p value for anovan MI per mouse per electrode for naive vs proficient for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan MI per mouse per electrode for events for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(2));
            fprintf(1, ['p value for anovan MI per mouse per electrode for groups for PAC theta/' freq_names{pacii+1} '= %d \n\n'],  p(3));
            
            %Calculate anovan for inteaction
            fprintf(1, ['ANOVAN for average MI for each electrode calculated per mouse without mouse as a factor for PAC theta/' freq_names{pacii+1} '\n'])
            [p,tbl,stats]=anovan(data_MI,{prof_naive, events, groups},'varnames',{'proficient_vs_naive','events','groups'},'display','off');
            fprintf(1, ['p value for anovan MI per mouse per electrode for naive vs proficient for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan MI per mouse per electrode for events for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(2));
            fprintf(1, ['p value for anovan MI per mouse per electrode for groups for PAC theta/' freq_names{pacii+1} '= %d \n\n'],  p(3));
            
        end
        
        %Now do the cumulative histograms and ranksums for MI for each electrode calculated with all sessons per mouse
        pvals=[];
        legends=[];
        for pacii=1:no_pacii
            
            ii_rank=0;
            mi_rank=[];
            maxmi=-200;
            minmi=200;
            
            try
                close(pacii+3)
            catch
            end
            hFig=figure(pacii+3);
            
            if length(eventType)>2
                set(hFig, 'units','normalized','position',[.1 .1 .7 .8])
            else
                set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            end
            
            
            for evNo=1:length(eventType)
                
                subplot(ceil(length(eventType)/2),2,evNo)
                hold on
                
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    
                    
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        
                        if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
                            
                            [f_mi,x_mi] = drg_ecdf(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)));
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
                            
                            
                            %Save data for ranksum
                            ii_rank=ii_rank+1;
                            mi_rank(ii_rank).mi=mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
                            mi_rank(ii_rank).per_ii=per_ii;
                            mi_rank(ii_rank).grNo=grNo;
                            mi_rank(ii_rank).evNo=evNo;
                            maxmi=max([maxmi max(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))]);
                            minmi=min([minmi min(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))]);
                        end
                    end
                    
                end
                
                title(evTypeLabels{evNo})
                xlabel('MI')
                ylabel('Probability')
            end
            
            
            
            for evNo=1:length(eventType)
                subplot(ceil(length(eventType)/2),2,evNo)
                xlim([minmi-0.1*(maxmi-minmi) maxmi+0.1*(maxmi-minmi)])
            end
            
            %suptitle(['Average MI for each electrode calculated per  for PAC theta/' freq_names{pacii+1}])
            
            %Now do the ranksums
            fprintf(1, ['Ranksum or t-test p values for average MI for each electrode calculated per mouse for PAC theta' freq_names{pacii+1} '\n'])
            prof_naive_leg{1}='Proficient';
            prof_naive_leg{2}='Naive';
            
            input_data=[];
            for ii=1:ii_rank
                input_data(ii).data=mi_rank(ii).mi;
                input_data(ii).description=[handles_drgb.drgbchoices.group_no_names{mi_rank(ii).grNo} ' ' evTypeLabels{mi_rank(ii).evNo} ' ' prof_naive_leg{mi_rank(ii).per_ii}];
            end
            [output_data] = drgMutiRanksumorTtest(input_data);
            
%             for ii=1:ii_rank
%                 for jj=ii+1:ii_rank
%                     [p, r_or_t]=drg_ranksum_or_ttest(mi_rank(ii).mi,mi_rank(jj).mi);
%                     if r_or_t==0
%                         fprintf(1, ['p value ranksum for ' handles_drgb.drgbchoices.group_no_names{mi_rank(ii).grNo} ' ' evTypeLabels{mi_rank(ii).evNo} ' ' prof_naive_leg{mi_rank(ii).per_ii} ' vs ' ...
%                             handles_drgb.drgbchoices.group_no_names{mi_rank(jj).grNo} ' ' evTypeLabels{mi_rank(jj).evNo} ' ' prof_naive_leg{mi_rank(jj).per_ii} ' =  %d\n'],p)
%                     else
%                         fprintf(1, ['p value t-test for ' handles_drgb.drgbchoices.group_no_names{mi_rank(ii).grNo} ' ' evTypeLabels{mi_rank(ii).evNo} ' ' prof_naive_leg{mi_rank(ii).per_ii} ' vs ' ...
%                             handles_drgb.drgbchoices.group_no_names{mi_rank(jj).grNo} ' ' evTypeLabels{mi_rank(jj).evNo} ' ' prof_naive_leg{mi_rank(jj).per_ii} ' =  %d\n'],p)
%                     end
%                     pvals=[pvals p];
%                 end
%             end
%             fprintf(1, ['\n\n'])
        end
%         
%         pFDR_mi_rank=drsFDRpval( pvals);
%         fprintf(1, ['pFDR for per mi per mouse, per electrode  = %d\n\n'],pFDR_mi_rank);
        
        %Now plot the average MI per mouse averaged over electrodes
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            data_MI_per_mouse=[];
            prof_naive_per_mouse=[];
            events_per_mouse=[];
            groups_per_mouse=[];
            
            %Plot the average
            try
                close(pacii+6)
            catch
            end
            hFig=figure(pacii+6);
            
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            hold on
            
            %             bar_lab_loc=[];
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
                        
                        if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
                            
                            include_group=1;
                            
                            %Compute per mouse avearge
                            each_mouse_average_MI=[];
                            no_mice_included=0;
                            for mouseNo=1:max(mean_MI_mouseNo_per_mouse)
                                if mouse_included(mouseNo)==1
                                    if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
                                        if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo))>0
                                            no_mice_included=no_mice_included+1;
                                            each_mouse_average_MI(no_mice_included)=mean(mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo)));
                                        end
                                    end
                                end
                            end
                            if no_mice_included>0
                                
                                include_group=1;
                                
                                if per_ii==1
                                    bar(bar_offset,mean(each_mouse_average_MI),'r','LineWidth', 3,'EdgeColor','none')
                                else
                                    bar(bar_offset,mean(each_mouse_average_MI),'b','LineWidth', 3,'EdgeColor','none')
                                end
                                
                                
                                plot(bar_offset,mean(each_mouse_average_MI),'ok','LineWidth', 3)
                                plot((bar_offset)*ones(1,length(each_mouse_average_MI)),each_mouse_average_MI,'o',...
                                    'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                                
                                if length(each_mouse_average_MI)>2
                                    CI = bootci(1000, {@mean, each_mouse_average_MI},'type','cper');
                                    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                                end
                                
                                %Show the mean in the cumulative histos
                                figure(pacii+3)
                                subplot(ceil(length(eventType)/2),2,evNo)
                                hold on
                                
                                
                                for jj=1:length(each_mouse_average_MI)
                                    this_f_mi=[];
                                    this_f_mi=cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_mi;
                                    
                                    this_x_mi=[];
                                    this_x_mi=cum_histodB.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_mi;
                                    
                                    xii_below=find(this_x_mi<each_mouse_average_MI(jj),1,'last');
                                    xii_above=find(this_x_mi>each_mouse_average_MI(jj),1,'first');
                                    
                                    slope=(this_f_mi(xii_above)-this_f_mi(xii_below))/(this_x_mi(xii_above)-this_x_mi(xii_below));
                                    intercept=this_f_mi(xii_above)-slope*this_x_mi(xii_above);
                                    
                                    this_f=slope*each_mouse_average_MI(jj)+intercept;
                                    
                                    if grNo==1
                                        if per_ii==1
                                            plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
                                        else
                                            plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
                                        end
                                    else
                                        if per_ii==1
                                            plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
                                        else
                                            plot(each_mouse_average_MI(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
                                        end
                                    end
                                    
                                end
                                
                                
                                figure(pacii+6)
                                hold on
                                
                                %Save data for anovan
                                data_MI_per_mouse=[data_MI_per_mouse mean_MI_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))];
                                prof_naive_per_mouse=[prof_naive_per_mouse per_ii*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))];
                                events_per_mouse=[events_per_mouse evNo*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))];
                                groups_per_mouse=[groups_per_mouse grNo*ones(1,sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)))];
                            end
                        end
                    end
                end
                
                
                if include_group==1
                    ii_gr_included=ii_gr_included+1;
                    groups_included(ii_gr_included)=grNo;
                end
            end
            
            title(['Average MI per mouse averaged over all electrodes for PAC theta/' freq_names{pacii+1}])
            
            
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
            
            ylabel('Modulation Index')
            
            %Calculate anovan for inteaction
            [p,tbl,stats]=anovan(data_MI_per_mouse,{prof_naive_per_mouse events_per_mouse, groups_per_mouse},'model','interaction','varnames',{'proficient_vs_naive','events','groups'},'display','off');
            fprintf(1, ['anovan MI per mouse averaged over electrodes for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(1));
            
            fprintf(1, ['anovan p value for naive vs proficient for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(1));
            fprintf(1, ['anovan p value for events for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(2));
            fprintf(1, ['anovan p value for groups for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(3));
            fprintf(1, ['anovan p value for naive vs proficient * events for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(1));
            fprintf(1, ['anovan p value for naive vs proficient * groups for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(2));
            fprintf(1, ['anovan p value for events * groups for PAC theta/' freq_names{pacii+1} '= %d \n\n'],  p(3));
            
        end
        
        for pacii=1:no_pacii
            for evNo=1:length(eventType)
                figure(pacii+3)
                subplot(ceil(length(eventType)/2),2,evNo)
                hold on
                try
                    legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
                        [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
                catch
                end
            end
            suptitle(['Average MI for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1}])
        end
        
        %Now do the cumulative histograms and ranksums for mean vector angle per electrode per mouse
        pvals=[];
        legends=[];
        cum_histoVA=[];
        mean_all_meansVA=[];
        max_all_shifted_meansVA=[];
        all_means_shifted_meanVAs=[];
        all_CIs_shifted_meanVAs=[];
                            
        for pacii=1:no_pacii
            
            ii_rank=0;
            VA_rank=[];
            maxVA=-2000;
            minVA=2000;
            
            try
                close(pacii+9)
            catch
            end
            hFig=figure(pacii+9);
            
            if length(eventType)>2
                set(hFig, 'units','normalized','position',[.1 .1 .7 .8])
            else
                set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            end
            
            
            for evNo=1:length(eventType)
                subplot(ceil(length(eventType)/2),2,evNo)
                hold on
                
                %Calculate the mean of the mean of each distribution
                all_meansVA=[];
                no_means=0;
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    
                    
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        
                        if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
                            
                            %Note: we have to shift the mean to accomodate
                            %all cumulative histograms in one angle axis
                            if (evNo==2)&(pacii==2)
                                pffft=1;
                            end
                            these_meanVAs=[];
                            these_meanVAs=mean_meanVectorAngle_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
                            
                            sztmVA=size(these_meanVAs);
                            shifted_meanVAs=zeros(sztmVA(1),sztmVA(2));
                            
                            this_meanVA=[];
                            this_meanVA=(180/pi)*circ_axial(circ_mean(these_meanVAs'*pi/180))';
                            
                            no_means=no_means+1;
                            all_meansVA(no_means)=this_meanVA;
                        end
                    end
                end
                
                
                mean_all_meansVA(pacii,evNo)=(180/pi)*circ_axial(circ_mean(all_meansVA'*pi/180))';
                
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    
                    
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        
                        if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
                            
                            %Note: we have to shift the mean to accomodate
                            %all cumulative histograms in one angle axis
                            if (evNo==2)&(pacii==1)
                                pffft=1;
                            end
                            these_meanVAs=[];
                            these_meanVAs=mean_meanVectorAngle_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
                            
                            sztmVA=size(these_meanVAs);
                            shifted_meanVAs=zeros(sztmVA(1),sztmVA(2));
                            
                            this_meanVA=[];
                            this_meanVA=(180/pi)*circ_axial(circ_mean(these_meanVAs'*pi/180))';
                            
                            shifted_meanVAs(these_meanVAs>180+this_meanVA)=-(360-these_meanVAs(these_meanVAs>180+this_meanVA));
                            shifted_meanVAs(these_meanVAs<this_meanVA-180)=360+these_meanVAs(these_meanVAs<this_meanVA-180);
                            shifted_meanVAs((these_meanVAs<=180+this_meanVA)&(these_meanVAs>=this_meanVA-180))=these_meanVAs((these_meanVAs<=180+this_meanVA)&(these_meanVAs>=this_meanVA-180));
                            
                            %Make sure they are all grouped on the same
                            %side of the cumulative histogram
                            if mean_all_meansVA(pacii,evNo)<180
                                if abs(this_meanVA-mean_all_meansVA(pacii,evNo))>abs(this_meanVA-360+mean_all_meansVA(pacii,evNo))
                                    shifted_meanVAs=shifted_meanVAs-360;
                                end
                            else
                                if abs(this_meanVA-mean_all_meansVA(pacii,evNo))<abs(this_meanVA-360+mean_all_meansVA(pacii,evNo))
                                    shifted_meanVAs=shifted_meanVAs-360;
                                end
                            end
                            
                            max_all_shifted_meansVA(pacii,evNo,grNo,per_ii)=max(shifted_meanVAs');
                            if length(eventType)>2
                                all_means_shifted_meanVAs(pacii,evNo,grNo,per_ii)=mean(shifted_meanVAs');
                                all_CIs_shifted_meanVAs(pacii,evNo,grNo,per_ii,1:2) = bootci(1000, {@mean, shifted_meanVAs'},'type','cper');
                            end
                            
                            [f_VA,x_VA] = drg_ecdf(shifted_meanVAs);
                            
                            
                            
                            cum_histoVA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_VA=f_VA;
                            cum_histoVA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_VA=x_VA;
                            if grNo==1
                                if per_ii==1
                                    legends.pacii(pacii).evNo(evNo).p1=plot(x_VA,f_VA,'Color',[1 0 0],'LineWidth',3);
                                else
                                    legends.pacii(pacii).evNo(evNo).p2=plot(x_VA,f_VA,'Color',[0 0 1],'LineWidth',3);
                                end
                            else
                                if per_ii==1
                                    legends.pacii(pacii).evNo(evNo).p3=plot(x_VA,f_VA,'Color',[1 0.7 0.7],'LineWidth',3);
                                else
                                    legends.pacii(pacii).evNo(evNo).p4=plot(x_VA,f_VA,'Color',[0.7 0.7 1],'LineWidth',3);
                                end
                            end
                            
                            
                            %Save data for ranksum
                            ii_rank=ii_rank+1;
                            VA_rank(ii_rank).meanVA=shifted_meanVAs;
                            VA_rank(ii_rank).per_ii=per_ii;
                            VA_rank(ii_rank).grNo=grNo;
                            VA_rank(ii_rank).evNo=evNo;
                            maxVA=max([maxVA max(shifted_meanVAs)]);
                            minVA=min([minVA min(shifted_meanVAs)]);
                        end
                    end
                    
                end
                
                title(evTypeLabels{evNo})
                xlabel('mean vector angle')
                ylabel('Probability')
            end
            
            
            for evNo=1:length(eventType)
                subplot(ceil(length(eventType)/2),2,evNo)
                xlim([minVA-0.1*(maxVA-minVA) maxVA+0.1*(maxVA-minVA)])
            end
            
            %suptitle(['Mean vector angle per mouse for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1} ])
            
            %Now do the Watson-Williams test
            fprintf(1, ['Watson-Williams test p values for mean vector angle for each electrode calculated per mouse for PAC theta ' freq_names{pacii+1} '\n'])
            prof_naive_leg{1}='Proficient';
            prof_naive_leg{2}='Naive';
            
            ww_test_out.ii_out=0;
            pvals=[];
            for ii=1:ii_rank
                for jj=ii+1:ii_rank
                    %p=ranksum(VA_rank(ii).meanVA,VA_rank(jj).meanVA);
                    p=circ_wwtest(pi*VA_rank(ii).meanVA/180,pi*VA_rank(jj).meanVA/180);
%                     fprintf(1, ['p value for ' handles_drgb.drgbchoices.group_no_names{VA_rank(ii).grNo} ' ' evTypeLabels{VA_rank(ii).evNo} ' ' prof_naive_leg{VA_rank(ii).per_ii} ' vs ' ...
%                         handles_drgb.drgbchoices.group_no_names{VA_rank(jj).grNo} ' ' evTypeLabels{VA_rank(jj).evNo} ' ' prof_naive_leg{VA_rank(jj).per_ii} ' =  %d\n'],p)
                    pvals=[pvals p];
                    ww_test_out.ii_out=ww_test_out.ii_out+1;
                    ww_test_out.p(ww_test_out.ii_out)=p;
                    ww_test_out.label{ww_test_out.ii_out}=[handles_drgb.drgbchoices.group_no_names{VA_rank(ii).grNo} ' ' evTypeLabels{VA_rank(ii).evNo} ' ' prof_naive_leg{VA_rank(ii).per_ii} ' vs ' ...
                        handles_drgb.drgbchoices.group_no_names{VA_rank(jj).grNo} ' ' evTypeLabels{VA_rank(jj).evNo} ' ' prof_naive_leg{VA_rank(jj).per_ii}];
                end
            end
            
            pFDR_VA_rank=drsFDRpval(pvals);
            
            %Now sort the data
            these_ii_pairs=[1:ww_test_out.ii_out];
            to_sort=[pvals' these_ii_pairs'];
            sorted_rows=sortrows(to_sort);
            ww_test_out.sorted_ii_pairs=sorted_rows(:,2);
            
            first_above=1;
            for ii_out=1:ww_test_out.ii_out
                if first_above==1
                    if ww_test_out.p(ww_test_out.sorted_ii_pairs(ii_out))>pFDR_VA_rank
                       first_above=0;
                       fprintf(1, ['\n\npFDR for per mean vector angle for each electrode calculated per mouse  = %d\n\n'],pFDR_VA_rank);
                    end
                end
                this_label=ww_test_out.label(ww_test_out.sorted_ii_pairs(ii_out));
                fprintf(1, ['p value for ' this_label{1} ' =  %d\n'],ww_test_out.p(ww_test_out.sorted_ii_pairs(ii_out)));
            end

            fprintf(1, ['\n\n'])
            
            
            
        end
        
        %Now plot the average VA per mouse averaged over electrodes
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            
            ii_gr_included=0;
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                
                for evNo=1:length(eventType)
                    
                    for per_ii=1:2
                        
                        
                        if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
                            
                            include_group=1;
                            
                            %Compute per mouse avearge
                            each_mouse_average_VA=[];
                            no_mice_included=0;
                            for mouseNo=1:max(mean_MI_mouseNo_per_mouse)
                                if mouse_included(mouseNo)==1
                                    if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
                                        if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo))>0
                                            no_mice_included=no_mice_included+1;
                                            these_VA=[];
                                            these_VA=mean_meanVectorAngle_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo));
                                            each_mouse_average_VA(no_mice_included)=(180/pi)*circ_axial(circ_mean(these_VA'*pi/180)');
                                        end
                                    end
                                end
                            end
                            if (pacii==1)&(evNo==2)
                                pffft=1;
                            end
                            if no_mice_included>0
                                
                                include_group=1;
                                
                                for no_mice=1:no_mice_included
                                    if each_mouse_average_VA(no_mice)>max_all_shifted_meansVA(pacii,evNo,grNo,per_ii)
                                        each_mouse_average_VA(no_mice)=each_mouse_average_VA(no_mice)-360;
                                    end
                                end
                                
                                %Show the mean in the cumulative histos
                                figure(pacii+9)
                                subplot(ceil(length(eventType)/2),2,evNo)
                                hold on
                                
                                
                                for jj=1:length(each_mouse_average_VA)
                                    this_f_VA=[];
                                    this_f_VA=cum_histoVA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_VA;
                                    
                                    this_x_VA=[];
                                    this_x_VA=cum_histoVA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_VA;
                                    
                                    xii_below=find(this_x_VA<each_mouse_average_VA(jj),1,'last');
                                    xii_above=find(this_x_VA>each_mouse_average_VA(jj),1,'first');
                                    
                                    slope=(this_f_VA(xii_above)-this_f_VA(xii_below))/(this_x_VA(xii_above)-this_x_VA(xii_below));
                                    intercept=this_f_VA(xii_above)-slope*this_x_VA(xii_above);
                                    
                                    this_f=slope*each_mouse_average_VA(jj)+intercept;
                                    
                                    if each_mouse_average_VA(jj)>max(this_x_VA)
                                        this_f=1;
                                    end
                                    
                                    if each_mouse_average_VA(jj)<min(this_x_VA)
                                        this_f=0;
                                    end
                                    
                                    if grNo==1
                                        if per_ii==1
                                            plot(each_mouse_average_VA(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
                                        else
                                            plot(each_mouse_average_VA(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
                                        end
                                    else
                                        if per_ii==1
                                            plot(each_mouse_average_VA(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
                                        else
                                            plot(each_mouse_average_VA(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
                                        end
                                    end
                                    
                                end
                                
                                
                            end
                        end
                    end
                end
                if include_group==1
                    ii_gr_included=ii_gr_included+1;
                    groups_included(ii_gr_included)=grNo;
                end
            end
        end
        
        
        for pacii=1:no_pacii
            %             figure(pacii+3)
            %             for evNo=1:length(eventType)
            %                 subplot(1,2,evNo)
            %                 hold on
            %                 try
            %                 legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
            %                     [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
            %                 catch
            %                 end
            %
            %             end
            %             suptitle(['Average per mouse MI for each electrode calculated per for PAC theta/' freq_names{pacii+1}])
            
            figure(pacii+9)
            for evNo=1:length(eventType)
                subplot(ceil(length(eventType)/2),2,evNo)
                hold on
                try
                    legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
                        [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
                catch
                end
                
            end
            suptitle(['Mean vector angle per mouse for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1} ])
        end
        
        %Now do the cumulative histograms and ranksums for peak angle per electrode per mouse
        pvals=[];
        legends=[];
        cum_histoPA=[];
        mean_all_meansPA=[];
        max_all_shifted_meansPA=[];
        all_means_shifted_meanPAs=[];
        all_CIs_shifted_meanPAs=[];
                            
        for pacii=1:no_pacii
            
            ii_rank=0;
            PA_rank=[];
            maxPA=-2000;
            minPA=2000;
            
            try
                close(pacii+12)
            catch
            end
            hFig=figure(pacii+12);
            
            if length(eventType)>2
                set(hFig, 'units','normalized','position',[.1 .1 .7 .8])
            else
                set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            end
            
            
            for evNo=1:length(eventType)
                subplot(ceil(length(eventType)/2),2,evNo)
                hold on
                
                %Calculate the mean of the mean of each distribution
                all_meansPA=[];
                no_means=0;
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    
                    
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        
                        if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
                            
                            %Note: we have to shift the mean to accomodate
                            %all cumulative histograms in one angle axis
                            if (evNo==2)&(pacii==2)
                                pffft=1;
                            end
                            these_meanPAs=[];
                            these_meanPAs=mean_peakAngle_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
                            
                            sztmPA=size(these_meanPAs);
                            shifted_meanPAs=zeros(sztmPA(1),sztmPA(2));
                            
                            this_meanPA=[];
                            this_meanPA=(180/pi)*circ_axial(circ_mean(these_meanPAs'*pi/180))';
                            
                            no_means=no_means+1;
                            all_meansPA(no_means)=this_meanPA;
                        end
                    end
                end
                
                
                mean_all_meansPA(pacii,evNo)=(180/pi)*circ_axial(circ_mean(all_meansPA'*pi/180))';
                
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    
                    
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        
                        if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
                            
                            %Note: we have to shift the mean to accomodate
                            %all cumulative histograms in one angle axis
                            if (evNo==2)&(pacii==1)
                                pffft=1;
                            end
                            these_meanPAs=[];
                            these_meanPAs=mean_peakAngle_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo));
                            
                            sztmPA=size(these_meanPAs);
                            shifted_meanPAs=zeros(sztmPA(1),sztmPA(2));
                            
                            this_meanPA=[];
                            this_meanPA=(180/pi)*circ_axial(circ_mean(these_meanPAs'*pi/180))';
                            
                            shifted_meanPAs(these_meanPAs>180+this_meanPA)=-(360-these_meanPAs(these_meanPAs>180+this_meanPA));
                            shifted_meanPAs(these_meanPAs<this_meanPA-180)=360+these_meanPAs(these_meanPAs<this_meanPA-180);
                            shifted_meanPAs((these_meanPAs<=180+this_meanPA)&(these_meanPAs>=this_meanPA-180))=these_meanPAs((these_meanPAs<=180+this_meanPA)&(these_meanPAs>=this_meanPA-180));
                            
                            %Make sure they are all grouped on the same
                            %side of the cumulative histogram
                            if mean_all_meansPA(pacii,evNo)<180
                                if abs(this_meanPA-mean_all_meansPA(pacii,evNo))>abs(this_meanPA-360+mean_all_meansPA(pacii,evNo))
                                    shifted_meanPAs=shifted_meanPAs-360;
                                end
                            else
                                if abs(this_meanPA-mean_all_meansPA(pacii,evNo))<abs(this_meanPA-360+mean_all_meansPA(pacii,evNo))
                                    shifted_meanPAs=shifted_meanPAs-360;
                                end
                            end
                            
                            max_all_shifted_meansPA(pacii,evNo,grNo,per_ii)=max(shifted_meanPAs');
                            if length(eventType)>2
                                all_means_shifted_meanPAs(pacii,evNo,grNo,per_ii)=mean(shifted_meanPAs');
                                all_CIs_shifted_meanPAs(pacii,evNo,grNo,per_ii,1:2) = bootci(1000, {@mean, shifted_meanPAs'},'type','cper');
                            end
                            
                            [f_PA,x_PA] = drg_ecdf(shifted_meanPAs);
                            
                            
                            
                            cum_histoPA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_PA=f_PA;
                            cum_histoPA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_PA=x_PA;
                            if grNo==1
                                if per_ii==1
                                    legends.pacii(pacii).evNo(evNo).p1=plot(x_PA,f_PA,'Color',[1 0 0],'LineWidth',3);
                                else
                                    legends.pacii(pacii).evNo(evNo).p2=plot(x_PA,f_PA,'Color',[0 0 1],'LineWidth',3);
                                end
                            else
                                if per_ii==1
                                    legends.pacii(pacii).evNo(evNo).p3=plot(x_PA,f_PA,'Color',[1 0.7 0.7],'LineWidth',3);
                                else
                                    legends.pacii(pacii).evNo(evNo).p4=plot(x_PA,f_PA,'Color',[0.7 0.7 1],'LineWidth',3);
                                end
                            end
                            
                            
                            %Save data for ranksum
                            ii_rank=ii_rank+1;
                            PA_rank(ii_rank).meanPA=shifted_meanPAs;
                            PA_rank(ii_rank).per_ii=per_ii;
                            PA_rank(ii_rank).grNo=grNo;
                            PA_rank(ii_rank).evNo=evNo;
                            maxPA=max([maxPA max(shifted_meanPAs)]);
                            minPA=min([minPA min(shifted_meanPAs)]);
                        end
                    end
                    
                end
                
                title(evTypeLabels{evNo})
                xlabel('Peak angle')
                ylabel('Probability')
            end
            
            
            for evNo=1:length(eventType)
                subplot(ceil(length(eventType)/2),2,evNo)
                xlim([minPA-0.1*(maxPA-minPA) maxPA+0.1*(maxPA-minPA)])
            end
            
            %suptitle(['Mean vector angle per mouse for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1} ])
            
            %Now do the Watson-Williams test
            fprintf(1, ['Watson-Williams test p values for peak angle for each electrode calculated per mouse for PAC theta ' freq_names{pacii+1} '\n'])
            prof_naive_leg{1}='Proficient';
            prof_naive_leg{2}='Naive';
            
            ww_test_out.ii_out=0;
            pvals=[];
            for ii=1:ii_rank
                for jj=ii+1:ii_rank
                    %p=ranksum(PA_rank(ii).meanPA,PA_rank(jj).meanPA);
                    p=circ_wwtest(pi*PA_rank(ii).meanPA/180,pi*PA_rank(jj).meanPA/180);
%                     fprintf(1, ['p value for ' handles_drgb.drgbchoices.group_no_names{PA_rank(ii).grNo} ' ' evTypeLabels{PA_rank(ii).evNo} ' ' prof_naive_leg{PA_rank(ii).per_ii} ' vs ' ...
%                         handles_drgb.drgbchoices.group_no_names{PA_rank(jj).grNo} ' ' evTypeLabels{PA_rank(jj).evNo} ' ' prof_naive_leg{PA_rank(jj).per_ii} ' =  %d\n'],p)
                    pvals=[pvals p];
                    ww_test_out.ii_out=ww_test_out.ii_out+1;
                    ww_test_out.p(ww_test_out.ii_out)=p;
                    ww_test_out.label{ww_test_out.ii_out}=[handles_drgb.drgbchoices.group_no_names{PA_rank(ii).grNo} ' ' evTypeLabels{PA_rank(ii).evNo} ' ' prof_naive_leg{PA_rank(ii).per_ii} ' vs ' ...
                        handles_drgb.drgbchoices.group_no_names{PA_rank(jj).grNo} ' ' evTypeLabels{PA_rank(jj).evNo} ' ' prof_naive_leg{PA_rank(jj).per_ii}];
                end
            end
            
            pFDR_PA_rank=drsFDRpval(pvals);
            
            %Now sort the data
            these_ii_pairs=[1:ww_test_out.ii_out];
            to_sort=[pvals' these_ii_pairs'];
            sorted_rows=sortrows(to_sort);
            ww_test_out.sorted_ii_pairs=sorted_rows(:,2);
            
            first_above=1;
            for ii_out=1:ww_test_out.ii_out
                if first_above==1
                    if ww_test_out.p(ww_test_out.sorted_ii_pairs(ii_out))>pFDR_PA_rank
                       first_above=0;
                       fprintf(1, ['\n\npFDR for peak angle for each electrode calculated per mouse  = %d\n\n'],pFDR_PA_rank);
                    end
                end
                this_label=ww_test_out.label(ww_test_out.sorted_ii_pairs(ii_out));
                fprintf(1, ['p value for ' this_label{1} ' =  %d\n'],ww_test_out.p(ww_test_out.sorted_ii_pairs(ii_out)));
            end

            fprintf(1, ['\n\n'])
            
            
            
        end
        
        %Now plot the average PA per mouse averaged over electrodes
        for pacii=1:no_pacii    %for amplitude bandwidths (beta, low gamma, high gamma)
            
            
            ii_gr_included=0;
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                
                for evNo=1:length(eventType)
                    
                    for per_ii=1:2
                        
                        
                        if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo))>1
                            
                            include_group=1;
                            
                            %Compute per mouse avearge
                            each_mouse_average_PA=[];
                            no_mice_included=0;
                            for mouseNo=1:max(mean_MI_mouseNo_per_mouse)
                                if mouse_included(mouseNo)==1
                                    if sum(handles_drgb.drgbchoices.group_no(handles_drgb.drgbchoices.mouse_no==mouseNo)==grNo)
                                        if sum((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo))>0
                                            no_mice_included=no_mice_included+1;
                                            these_PA=[];
                                            these_PA=mean_peakAngle_per_mouse((~isnan(mean_MI_per_mouse))&(mean_MI_perii_per_mouse==per_ii)&(mean_MI_pacii_per_mouse==pacii)&(mean_MI_evNo_per_mouse==evNo)&(mean_MI_group_no_per_mouse==grNo)&(mean_MI_mouseNo_per_mouse==mouseNo));
                                            each_mouse_average_PA(no_mice_included)=(180/pi)*circ_axial(circ_mean(these_PA'*pi/180)');
                                        end
                                    end
                                end
                            end
                            if (pacii==1)&(evNo==2)
                                pffft=1;
                            end
                            if no_mice_included>0
                                
                                include_group=1;
                                
                                for no_mice=1:no_mice_included
                                    if each_mouse_average_PA(no_mice)>max_all_shifted_meansPA(pacii,evNo,grNo,per_ii)
                                        each_mouse_average_PA(no_mice)=each_mouse_average_PA(no_mice)-360;
                                    end
                                end
                                
                                %Show the mean in the cumulative histos
                                figure(pacii+12)
                                subplot(ceil(length(eventType)/2),2,evNo)
                                hold on
                                
                                
                                for jj=1:length(each_mouse_average_PA)
                                    this_f_PA=[];
                                    this_f_PA=cum_histoPA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).f_PA;
                                    
                                    this_x_PA=[];
                                    this_x_PA=cum_histoPA.pacii(pacii).grNo(grNo).evNo(evNo).per_ii(per_ii).x_PA;
                                    
                                    xii_below=find(this_x_PA<each_mouse_average_PA(jj),1,'last');
                                    xii_above=find(this_x_PA>each_mouse_average_PA(jj),1,'first');
                                    
                                    slope=(this_f_PA(xii_above)-this_f_PA(xii_below))/(this_x_PA(xii_above)-this_x_PA(xii_below));
                                    intercept=this_f_PA(xii_above)-slope*this_x_PA(xii_above);
                                    
                                    this_f=slope*each_mouse_average_PA(jj)+intercept;
                                    
                                    if each_mouse_average_PA(jj)>max(this_x_PA)
                                        this_f=1;
                                    end
                                    
                                    if each_mouse_average_PA(jj)<min(this_x_PA)
                                        this_f=0;
                                    end
                                    
                                    if grNo==1
                                        if per_ii==1
                                            plot(each_mouse_average_PA(jj),this_f,'o','MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
                                        else
                                            plot(each_mouse_average_PA(jj),this_f,'o','MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
                                        end
                                    else
                                        if per_ii==1
                                            plot(each_mouse_average_PA(jj),this_f,'o','MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
                                        else
                                            plot(each_mouse_average_PA(jj),this_f,'o','MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
                                        end
                                    end
                                    
                                end
                                
                                
                            end
                        end
                    end
                end
                if include_group==1
                    ii_gr_included=ii_gr_included+1;
                    groups_included(ii_gr_included)=grNo;
                end
            end
        end
        
        
        for pacii=1:no_pacii
            %             figure(pacii+3)
            %             for evNo=1:length(eventType)
            %                 subplot(1,2,evNo)
            %                 hold on
            %                 try
            %                 legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
            %                     [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
            %                 catch
            %                 end
            %
            %             end
            %             suptitle(['Average per mouse MI for each electrode calculated per for PAC theta/' freq_names{pacii+1}])
            
            figure(pacii+12)
            for evNo=1:length(eventType)
                subplot(ceil(length(eventType)/2),2,evNo)
                hold on
                try
                    legend([legends.pacii(pacii).evNo(evNo).p1 legends.pacii(pacii).evNo(evNo).p2 legends.pacii(pacii).evNo(evNo).p3 legends.pacii(pacii).evNo(evNo).p4],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
                        [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
                catch
                end
                
            end
            suptitle(['Peak angle per mouse for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1} ])
        end
        
        
        %If this is for concentration plot as a function of concentration
        if length(eventType)>2
            
            
            for pacii=1:no_pacii
                try
                    close(pacii+12)
                catch
                end
                hFig=figure(pacii+12);
                
                set(hFig, 'units','normalized','position',[.1 .1 .4 .4])
                hold on
                
                
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    for per_ii=1:2
                        
                        these_means=zeros(length(eventType),1);
                        these_means(:,1)=all_means_shifted_meanVAs(pacii,:,grNo,per_ii);
                        
                        these_CIs=zeros(length(eventType),2);
                        these_CIs(:,:)=all_CIs_shifted_meanVAs(pacii,:,grNo,per_ii,1:2);
                        
                        
                        if grNo==1
                            if per_ii==1
                                p1=plot(concs,these_means,'-o','Color',[1 0 0], 'MarkerFace',[1 0 0],'MarkerEdge',[1 0 0],'MarkerSize',10)
                            else
                                p2=plot(concs,these_means,'-o','Color',[0 0 1],'MarkerFace',[0 0 1],'MarkerEdge',[0 0 1],'MarkerSize',10)
                            end
                        else
                            if per_ii==1
                                p3=plot(concs,these_means,'-o','Color',[1 0.7 0.7],'MarkerFace',[1 0.7 0.7],'MarkerEdge',[1 0.7 0.7],'MarkerSize',10)
                            else
                                p4=plot(concs,these_means,'-o','Color',[0.7 0.7 1],'MarkerFace',[0.7 0.7 1],'MarkerEdge',[0.7 0.7 1],'MarkerSize',10)
                            end
                        end
                        
                        for evTy=1:length(eventType)
                            plot([concs(evTy) concs(evTy)],these_CIs(evTy,:),'-k')
                        end
                        
                        
                    end
                end
                legend([p1 p2 p3 p4],[handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive'],...
                    [handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
                xlabel('Percent odor in mineral oil')
                ylabel('Phase (deg)')
                title(['Mean vector angle per mouse for each electrode calculated per mouse for PAC theta/' freq_names{pacii+1} ])
                 
            end
        end
        pfft=1;
        
         
    case 20
        % Multiclass ROC analysis of LFP power differences for naive and proficient
        % mice for different epochs (concentrations or S+ vs. S-) and different groups. Analyzed per mouse
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
                                                    evNo_out(evNo).noWB=evNo_out(evNo).noWB+1;
                                                    evNo_out(evNo).delta_dB_powerEvWB(evNo_out(evNo).noWB,1:length(frequency))=mean(this_dB_power-this_dB_powerref,1);
                                                    evNo_out(evNo).per_ii(evNo_out(evNo).noWB)=per_ii;
                                                    evNo_out(evNo).groupNo(evNo_out(evNo).noWB)=handles_drgb.drgbchoices.group_no(fileNo);
                                                    
                                                    noWB_for_evNo(evNo)=evNo_out(evNo).noWB;
                                                    
                                                    
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
                                                        
                                                        evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).noWB,bwii)=mean(this_delta_dB_powerEv);
                                                        
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
        
        pFDRauROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for per electrode auROC  = %d\n\n'],pFDRauROC);
        
        per_mouse_pFDRauROC=drsFDRpval(per_mouse_p_vals_ROC);
        fprintf(1, ['pFDR for per mouse auROC  = %d\n\n'],per_mouse_pFDRauROC);
        
        
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
            
            fprintf(1, ['\n\n'])
            fprintf(1, ['ANOVAN for delta dB power per mouse per electrode, mouse as random factor\n\n'])
            
            
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
            
            
            
            %Calculate anovan for inteaction
            [p,tbl,stats]=anovan(data_delta_dB,{prof_naive events mice electrodes},'varnames',{'proficient_vs_naive','events','groups','mice'},'display','off','random',4);
            fprintf(1, ['p value for anovan delta dB power per mouse per electrode for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan delta dB power  per mouse per electrode for events ' freq_names{bwii} '= %d \n'],  p(2));
            fprintf(1, ['p value for anovan delta dB power  per mouse per electrode for groups for ' freq_names{bwii} '= %d \n\n'],  p(3));
            
            
            %Plot the cumulative histos and do ranksum
            %Plot the average
            try
                close(bwii+4)
            catch
            end
            hFig=figure(bwii+4);
            
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
                            dBpower_rank(ii_rank).delta_dBpower=delta_dB_per_mouse((delta_dB_perii_per_mouse==per_ii)&(delta_dB_evNo_per_mouse==evNo)&(delta_dB_bwii_per_mouse==bwii)&(delta_dB_group_no_per_mouse==grNo));
                            dBpower_rank(ii_rank).per_ii=per_ii;
                            dBpower_rank(ii_rank).grNo=grNo;
                            dBpower_rank(ii_rank).evNo=evNo;
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
            
            %Now do the ranksums
            fprintf(1, ['Ranksum or t-test p values for delta dB power per electrode for ' freq_names{bwii} '\n'])
            prof_naive_leg{1}='Proficient';
            prof_naive_leg{2}='Naive';
            for ii=1:ii_rank
                for jj=ii+1:ii_rank
                    [p, r_or_t]=drg_ranksum_or_ttest(dBpower_rank(ii).delta_dBpower,dBpower_rank(jj).delta_dBpower);
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
            
            fprintf(1, ['\n\n'])
            fprintf(1, ['ANOVAN for delta dB power per mouse, electrode average\n\n'])
            
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
            
            
            
            %Calculate anovan for inteaction
            
            
            [p,tbl,stats]=anovan(per_mouse_data_delta_dB,{per_mouse_prof_naive per_mouse_events per_mouse_groups},'varnames',{'proficient_vs_naive','events','groups'},'display','off');
            fprintf(1, ['p value for anovan delta dB histogram per mouse, electrode avearage for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan delta dB histogram  per mouse, electrode avearage for events ' freq_names{bwii} '= %d \n'],  p(2));
            fprintf(1, ['p value for anovan delta dB histogram  per mouse, electrode avearage for groups ' freq_names{bwii} '= %d \n\n'],  p(2));
            
            
            
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
            
            fprintf(1, ['\n\n'])
            fprintf(1, ['ANOVAN for auROC per mouse per electrode\n\n'])
            
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
            
            if sum(eventType==3)>0
                %S+ S-
                %Calculate anovan for inteaction
                fprintf(1, ['ANOVAN for auROC per mouse per electrode for concentrations separated by two steps, mouse random factor\n'])
                [p,tbl,stats]=anovan(data_auROC,{prof_naive groups mice},'varnames',{'proficient_vs_naive','groups','mice'},'display','off','random',3);
                fprintf(1, ['p value for anovan auROC per mouse per electrode for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
                fprintf(1, ['p value for anovan auROC  per mouse per electrode for groups for ' freq_names{bwii} '= %d \n\n'],  p(2));
                
                %Calculate anovan for inteaction
                fprintf(1, ['ANOVAN for auROC per mouse per electrode for concentrations separated by two steps\n'])
                [p,tbl,stats]=anovan(data_auROC,{prof_naive groups},'varnames',{'proficient_vs_naive','groups'},'display','off');
                fprintf(1, ['p value for anovan auROC per mouse per electrode for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
                fprintf(1, ['p value for anovan auROC  per mouse per electrode for groups for ' freq_names{bwii} '= %d \n\n'],  p(2));
            else
                %Calculate anovan for inteaction
                fprintf(1, ['ANOVAN for auROC per mouse per electrode for concentrations separated by two steps, mouse random factor\n'])
                [p,tbl,stats]=anovan(data_auROC,{prof_naive within_between groups mice},'varnames',{'proficient_vs_naive','within_between','groups','mice'},'display','off','random',4);
                fprintf(1, ['p value for anovan auROC per mouse per electrode for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
                fprintf(1, ['p value for anovan auROC  per mouse per electrode for within vs between ' freq_names{bwii} '= %d \n'],  p(2));
                fprintf(1, ['p value for anovan auROC  per mouse per electrode for groups for ' freq_names{bwii} '= %d \n\n'],  p(3));
                
                
                %Calculate anovan for inteaction
                fprintf(1, ['ANOVAN for auROC per mouse per electrode for concentrations separated by two steps\n'])
                [p,tbl,stats]=anovan(data_auROC,{prof_naive within_between groups},'varnames',{'proficient_vs_naive','within_between','groups'},'display','off');
                fprintf(1, ['p value for anovan auROC per mouse per electrode for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
                fprintf(1, ['p value for anovan auROC  per mouse per electrode for within vs between ' freq_names{bwii} '= %d \n'],  p(2));
                fprintf(1, ['p value for anovan auROC  per mouse per electrode for groups for ' freq_names{bwii} '= %d \n\n'],  p(3));
                
            end
        end
        
        fprintf(1, ['\n\n'])
        
        
        
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
                
                title([freq_names{bwii} ' auROC per mouse, per electrode'])
                
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
        
        
        
        fprintf(1, ['Pairwise auROC analysis for coherence analysis of LFP\n\n'])
        %         p_vals=[];
        no_files=max(files);
        
        
        %         no_ROCs=0;
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
                    for evNo=1:length(eventType)
                        for bwii=1:no_bandwidths
                            theseEvNos(evNo,bwii).noEv=0;
                        end
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
        fprintf(1, '\n\n')
        
        pFDRauROC=drsFDRpval(p_vals_ROC);
        fprintf(1, ['pFDR for per electrode pair auROC  = %d\n\n'],pFDRauROC);
        
        per_mouse_pFDRauROC=drsFDRpval(per_mouse_p_vals_ROC);
        fprintf(1, ['pFDR for per mouse auROC  = %d\n\n'],per_mouse_pFDRauROC);
        
        
        %Now plot the average per electrode, all sessions for each mouse
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
            
            fprintf(1, ['\n\n'])
            fprintf(1, ['ANOVAN for delta dB power per mouse per electrode, mouse as random factor\n\n'])
            
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                
                include_group=0;
                
                for evNo=1:length(eventType)
                    
                    per_included=0;
                    these_dB_per_e=[];
                    there_are_NaNs=0;
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        %bar_offset=21-evNo*3+(2-per_ii);
                        if sum(eventType==3)>0
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(evNo-1);
                        else
                            bar_offset=(grNo-1)*(3.5*length(eventType))+(2-(per_ii-1))+3*(length(eventType)-evNo);
                        end
                        
                        these_offsets(per_ii)=bar_offset;
                        
                        if sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))>0
                            
                            include_group=1;
                            
                            per_included=per_included+1;
                            
                            if per_ii==1
                                bar(bar_offset,mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))),'r','LineWidth', 3)
                            else
                                bar(bar_offset,mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))),'b','LineWidth', 3)
                            end
                            
                            %Individual points; in the future add lines linking the
                            %points?
                            plot((bar_offset)*ones(1,sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))),...
                                deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)),'o',...
                                'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            data_for_lines(per_ii).these_dB_per_e(1:length(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))))=...
                                deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo));
                            data_for_lines(per_ii).these_mice(1:length(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))))=...
                                deltaCxy_mouseNo_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo));
                            data_for_lines(per_ii).these_bar_offsets=bar_offset;
                            
                            
                            %Average and CI
                            plot(bar_offset,mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))),'ok','LineWidth', 3)
                            if sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))>2
                                CI = bootci(1000, {@mean, deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))},'type','cper');
                                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                            end
                            
                            %Save data for anovan
                            data_delta_dB=[data_delta_dB deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))];
                            prof_naive=[prof_naive per_ii*ones(1,sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)))];
                            events=[events evNo*ones(1,sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)))];
                            mice=[mice deltaCxy_mouseNo_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))];
                            electrodes=[electrodes deltaCxy_elec_pair_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))];
                            groups=[groups deltaCxy_group_no_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))];
                            
                        end
                    end
                    
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
            
            title([freq_names{bwii} ' average delta coherence per mouse, per electrode'])
            
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
            
            ylabel('Delta coherence')
            
            
            
            %Calculate anovan for inteaction
            [p,tbl,stats]=anovan(data_delta_dB,{prof_naive events mice electrodes},'varnames',{'proficient_vs_naive','events','groups','mice'},'display','off','random',4);
            fprintf(1, ['p value for anovan delta dB power per mouse per electrode for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan delta dB power  per mouse per electrode for events ' freq_names{bwii} '= %d \n'],  p(2));
            fprintf(1, ['p value for anovan delta dB power  per mouse per electrode for groups for ' freq_names{bwii} '= %d \n\n'],  p(3));
            
            
            
        end
        fprintf(1, ['\n\n'])
        
        %For rm anova use a for loop to make a table with the following columns:
        %proficient_naive
        %group
        %mouse
        %event
        %for each electrode one column with dB power
        
        %I do not think repeated measures anova is granted here. I am
        %commenting this out
        
        
        %
        %         for bwii=1:no_bandwidths
        %
        %             fprintf(1, ['Repeated measures anova for delta coherence during odor' freq_names{bwii} '\n\n'])
        %
        %             %For all the condiitons setup the table
        %             kk=0;
        %             mouse_col=[];
        %             naive_proficient_col=[];
        %             event_col=[];
        %             group_col=[];
        %             deltaCxy_cols=[];
        %             Cxy_cols=[];
        %             for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
        %
        %
        %
        %                 for evNo=1:length(eventType)
        %
        %
        %                     for per_ii=1:2
        %
        %
        %                         this_ii=(deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_mouseNo_per_mouse==mouseNo)&(deltaCxy_elec_pair_per_mouse==1);
        %                         these_elecs=(deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_mouseNo_per_mouse==mouseNo);
        %
        %                         if sum(this_ii)>0
        %                             kk=kk+1;
        %
        %                             mouse_col(kk,1)=deltaCxy_mouseNo_per_mouse(this_ii);
        %                             naive_proficient_col(kk,1)=deltaCxy_perii_per_mouse(this_ii);
        %                             event_col(kk,1)=deltaCxy_evNo_per_mouse(this_ii);
        %                             group_col(kk,1)=deltaCxy_group_no_per_mouse(this_ii);
        %                             deltaCxy_cols(kk,1:max(deltaCxy_elec_pair_per_mouse))=deltaCxy_per_mouse(these_elecs);
        %                             Cxy_cols(kk,1:max(deltaCxy_elec_pair_per_mouse))=Cxy_per_mouse(these_elecs);
        %                         end
        %                     end
        %                 end
        %
        %             end
        %
        %             t = table(group_col,naive_proficient_col,event_col,mouse_col,deltaCxy_cols(:,1),deltaCxy_cols(:,2),...
        %                 deltaCxy_cols(:,3),deltaCxy_cols(:,4),deltaCxy_cols(:,5),deltaCxy_cols(:,6),deltaCxy_cols(:,7),deltaCxy_cols(:,8),...
        %                 deltaCxy_cols(:,9),deltaCxy_cols(:,10),deltaCxy_cols(:,11),deltaCxy_cols(:,12),deltaCxy_cols(:,13),deltaCxy_cols(:,14),...
        %                 deltaCxy_cols(:,15),deltaCxy_cols(:,16),...
        %                 'VariableNames',{'Genotype','Proficiency','Event','Mouse','deltaCxy1','deltaCxy2',...
        %                 'deltaCxy3','deltaCxy4','deltaCxy5','deltaCxy6','deltaCxy7','deltaCxy8','deltaCxy9','deltaCxy10',...
        %                 'deltaCxy11','deltaCxy12','deltaCxy13','deltaCxy14','deltaCxy15','deltaCxy16'});
        %             electrode_pairs=[1:16]'; %these are the electrode pairs
        %             rm = fitrm(t,'deltaCxy1-deltaCxy16 ~ Genotype + Proficiency + Event + Mouse','WithinDesign',electrode_pairs)
        %
        %             % Perform repeated measures analysis of variance.
        % %             ranovatbl = ranova(rm)
        %             anovatbl=anova(rm,'WithinModel',electrode_pairs)
        %             fprintf(1, '\n\n\n\n')
        %         end
        %
        %
        %
        %         for bwii=1:no_bandwidths
        %
        %             fprintf(1, ['Repeated measures anova for coherence before the odor ' freq_names{bwii} '\n\n'])
        %
        %             %For all the condiitons setup the table
        %             kk=0;
        %             mouse_col=[];
        %             naive_proficient_col=[];
        %             event_col=[];
        %             group_col=[];
        %             deltaCxy_cols=[];
        %             Cxy_cols=[];
        %             for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
        %
        %
        %
        %                 for evNo=1:length(eventType)
        %
        %
        %                     for per_ii=1:2
        %
        %
        %                         this_ii=(deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_mouseNo_per_mouse==mouseNo)&(deltaCxy_elec_pair_per_mouse==1);
        %                         these_elecs=(deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_mouseNo_per_mouse==mouseNo);
        %
        %                         if sum(this_ii)>0
        %                             kk=kk+1;
        %
        %                             mouse_col(kk,1)=deltaCxy_mouseNo_per_mouse(this_ii);
        %                             naive_proficient_col(kk,1)=deltaCxy_perii_per_mouse(this_ii);
        %                             event_col(kk,1)=deltaCxy_evNo_per_mouse(this_ii);
        %                             group_col(kk,1)=deltaCxy_group_no_per_mouse(this_ii);
        %                             deltaCxy_cols(kk,1:max(deltaCxy_elec_pair_per_mouse))=deltaCxy_per_mouse(these_elecs);
        %                             Cxy_cols(kk,1:max(deltaCxy_elec_pair_per_mouse))=Cxy_per_mouse(these_elecs);
        %                         end
        %                     end
        %                 end
        %
        %             end
        %
        %             t = table(group_col,naive_proficient_col,event_col,mouse_col,Cxy_cols(:,1),Cxy_cols(:,2),...
        %                 Cxy_cols(:,3),Cxy_cols(:,4),Cxy_cols(:,5),Cxy_cols(:,6),Cxy_cols(:,7),Cxy_cols(:,8),...
        %                 Cxy_cols(:,9),Cxy_cols(:,10),Cxy_cols(:,11),Cxy_cols(:,12),Cxy_cols(:,13),Cxy_cols(:,14),...
        %                 Cxy_cols(:,15),Cxy_cols(:,16),...
        %                 'VariableNames',{'Genotype','Proficiency','Event','Mouse','Cxy1','Cxy2',...
        %                 'Cxy3','Cxy4','Cxy5','Cxy6','Cxy7','Cxy8','Cxy9','Cxy10',...
        %                 'Cxy11','Cxy12','Cxy13','Cxy14','Cxy15','Cxy16'});
        %             electrode_pairs=[1:16]'; %these are the electrode pairs
        %             rm = fitrm(t,'Cxy1-Cxy16 ~ Genotype + Proficiency + Event + Mouse','WithinDesign',electrode_pairs)
        %             % Perform repeated measures analysis of variance.
        % %             ranovatbl = ranova(rm)
        %             anovatbl=anova(rm,'WithinModel',electrode_pairs)
        %             fprintf(1, '\n\n\n\n')
        %         end
        %
        %
        
        
        %Now plot the average per mouse LFP power
        for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
            %Plot the average
            
            try
                close(bwii+4)
            catch
            end
            hFig=figure(bwii+4);
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
            
            fprintf(1, ['\n\n'])
            fprintf(1, ['ANOVAN for delta dB power per mouse, electrode average\n\n'])
            
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
                        for mouseNo=1:max(deltaCxy_mouseNo_per_mouse)
                            if sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_mouseNo_per_mouse==mouseNo)&(deltaCxy_group_no_per_mouse==grNo))>0
                                no_mice_for_this_group=no_mice_for_this_group+1;
                                each_mouse_average_delta_dB(no_mice_for_this_group)=mean(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_mouseNo_per_mouse==mouseNo)&(deltaCxy_group_no_per_mouse==grNo)));
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
                title([freq_names{bwii} ' delta coherence per mouse, electrode average'])
            else
                title([freq_names{bwii} ' delta coherence per mouse, electrode avearage concentrations two steps appart'])
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
            
            
            ylabel('Delta coherence')
            
            
            
            %Calculate anovan for inteaction
            
            
            [p,tbl,stats]=anovan(per_mouse_data_delta_dB,{per_mouse_prof_naive per_mouse_events per_mouse_groups},'varnames',{'proficient_vs_naive','events','groups'},'display','off');
            fprintf(1, ['p value for anovan delta dB histogram per mouse, electrode avearage for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan delta dB histogram  per mouse, electrode avearage for events ' freq_names{bwii} '= %d \n'],  p(2));
            fprintf(1, ['p value for anovan delta dB histogram  per mouse, electrode avearage for groups ' freq_names{bwii} '= %d \n\n'],  p(2));
            
            
            
        end
        
        
        
        %Display cumulative histograms for delta coherence for average per electrode per mouse and do ranksum
        pvals=[];
        ranksum_deltaCxy=[];
        if sum(eventType==3)>0
            for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
                
                %fprintf for ranksums
                fprintf(1, ['Ranksum or t-test p values for delta coherence for ' freq_names{bwii} '\n'])
                
                ii_rank=0;
                for evNo=1:length(eventType)
                    %Plot the average
                    
                    try
                        close(bwii+8+(evNo-1)*4)
                    catch
                    end
                    hFig=figure(bwii+8+(evNo-1)*4);
                    
                    set(hFig, 'units','normalized','position',[.2 .2 .6 .6])
                    
                    set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                    hold on
                    for grNo=1:max(handles_drgb.drgbchoices.group_no)
                        
                        include_group=0;
                        
                        
                        for per_ii=1:2      %performance bins. blue = naive, red = proficient
                            
                            if sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))>0
                                
                                include_group=1;
                                
                                
                                if per_ii==1
                                    if grNo==1
                                        [f_aic,x_aic] = drg_ecdf(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
                                        plot(x_aic,f_aic,'Color',[1 0 0],'LineWidth',3)
                                    else
                                        [f_aic,x_aic] = drg_ecdf(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
                                        plot(x_aic,f_aic,'Color',[0 0 1],'LineWidth',3)
                                    end
                                else
                                    if grNo==1
                                        [f_aic,x_aic] = drg_ecdf(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
                                        plot(x_aic,f_aic,'Color',[1 0 0])
                                    else
                                        [f_aic,x_aic] = drg_ecdf(deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
                                        plot(x_aic,f_aic,'Color',[0 0 1])
                                    end
                                end
                                
                                
                                %Save data for ranksum
                                ii_rank=ii_rank+1;
                                ranksum_deltaCxy(ii_rank).deltaCxy=deltaCxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_evNo_per_mouse==evNo)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo));
                                ranksum_deltaCxy(ii_rank).per_ii=per_ii;
                                ranksum_deltaCxy(ii_rank).grNo=grNo;
                                ranksum_deltaCxy(ii_rank).evNo=evNo;
                            end
                        end
                        
                        if include_group==1
                            ii_gr_included=ii_gr_included+1;
                            groups_included(ii_gr_included)=grNo;
                        end
                    end
                    title([freq_names{bwii} ' delta coherence per electrode (per mouse) for ' evTypeLabels{evNo}])
                    
                    legend([handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive']...
                        ,[handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
                    
                    xlabel('delta coherence')
                    ylabel('Cumulative probability')
                    
                    
                    prof_naive_leg{1}='Proficient';
                    prof_naive_leg{2}='Naive';
                    
                end
                
                
                for ii=1:ii_rank
                    for jj=ii+1:ii_rank
                        [p, r_or_t]=drg_ranksum_or_ttest(ranksum_deltaCxy(ii).deltaCxy,ranksum_deltaCxy(jj).deltaCxy);
                        if r_or_t==0
                            fprintf(1, ['p value ranksum for ' handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(ii).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(ii).per_ii} ' ' evTypeLabels{ranksum_deltaCxy(ii).evNo} ' vs ' ...
                                handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(jj).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(jj).per_ii} ' ' evTypeLabels{ranksum_deltaCxy(jj).evNo} ' =  %d\n'],p)
                        else
                            fprintf(1, ['p value t-test for ' handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(ii).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(ii).per_ii} ' ' evTypeLabels{ranksum_deltaCxy(ii).evNo} ' vs ' ...
                                handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(jj).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(jj).per_ii} ' ' evTypeLabels{ranksum_deltaCxy(jj).evNo} ' =  %d\n'],p)
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
        fprintf(1, ['\n\n'])
        
        
        %Display cumulative histograms for  coherence before odor on for average per electrode per mouse and do ranksum
        pvals=[];
        ranksum_Cxy=[];
        if sum(eventType==3)>0
            for bwii=1:no_bandwidths    %for bandwidths (theta, beta, low gamma, high gamma)
                
                %fprintf for ranksums
                fprintf(1, ['Ranksum or t-test p values for coherence before odor for ' freq_names{bwii} '\n'])
                
                ii_rank=0;
                
                
                try
                    close(bwii+16)
                catch
                end
                hFig=figure(bwii+16);
                
                set(hFig, 'units','normalized','position',[.2 .2 .6 .6])
                
                set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                hold on
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    include_group=0;
                    
                    
                    for per_ii=1:2      %performance bins. blue = naive, red = proficient
                        
                        if sum((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo))>0
                            
                            include_group=1;
                            
                            
                            if per_ii==1
                                if grNo==1
                                    [f_aic,x_aic] = drg_ecdf(Cxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
                                    plot(x_aic,f_aic,'Color',[1 0 0],'LineWidth',3)
                                else
                                    [f_aic,x_aic] = drg_ecdf(Cxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
                                    plot(x_aic,f_aic,'Color',[0 0 1],'LineWidth',3)
                                end
                            else
                                if grNo==1
                                    [f_aic,x_aic] = drg_ecdf(Cxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
                                    plot(x_aic,f_aic,'Color',[1 0 0])
                                else
                                    [f_aic,x_aic] = drg_ecdf(Cxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo)));
                                    plot(x_aic,f_aic,'Color',[0 0 1])
                                end
                            end
                            
                            
                            %Save data for ranksum
                            ii_rank=ii_rank+1;
                            ranksum_Cxy(ii_rank).Cxy=Cxy_per_mouse((deltaCxy_perii_per_mouse==per_ii)&(deltaCxy_bwii_per_mouse==bwii)&(deltaCxy_group_no_per_mouse==grNo));
                            ranksum_Cxy(ii_rank).per_ii=per_ii;
                            ranksum_Cxy(ii_rank).grNo=grNo;
                            
                        end
                    end
                    
                    if include_group==1
                        ii_gr_included=ii_gr_included+1;
                        groups_included(ii_gr_included)=grNo;
                    end
                end
                title([freq_names{bwii} ' coherence before odor per electrode (per mouse)'])
                
                legend([handles_drgb.drgbchoices.group_no_names{1} ' proficient'],[handles_drgb.drgbchoices.group_no_names{1} ' naive']...
                    ,[handles_drgb.drgbchoices.group_no_names{2} ' proficient'],[handles_drgb.drgbchoices.group_no_names{2} ' naive'])
                
                xlabel('coherence')
                ylabel('Cumulative probability')
                
                
                prof_naive_leg{1}='Proficient';
                prof_naive_leg{2}='Naive';
                
                
                
                
                for ii=1:ii_rank
                    for jj=ii+1:ii_rank
                        [p, r_or_t]=drg_ranksum_or_ttest(ranksum_Cxy(ii).Cxy,ranksum_Cxy(jj).Cxy);
                        if r_or_t==0
                            fprintf(1, ['p value ranksum for ' handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(ii).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(ii).per_ii}  ' vs ' ...
                                handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(jj).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(jj).per_ii}  ' =  %d\n'],p)
                        else
                            fprintf(1, ['p value t-test for ' handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(ii).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(ii).per_ii}  ' vs ' ...
                                handles_drgb.drgbchoices.group_no_names{ranksum_deltaCxy(jj).grNo} ' ' prof_naive_leg{ranksum_deltaCxy(jj).per_ii}  ' =  %d\n'],p)
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
        fprintf(1, ['\n\n'])
        
        
        
        %         save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix]);
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
        
        out_times=handles_drgb.drgb.file(1).waveERWA_out_times;
        
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
        
        
        fprintf(1, ['PAC power analysis for Justin''s paper\n\n'])
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
        prof_naive_leg{1}='Naive';
        prof_naive_leg{2}='Proficient';
        
        
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
        
end


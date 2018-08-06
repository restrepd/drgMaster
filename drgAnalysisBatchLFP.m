function drgAnalysisBatchLFP(handles)

%drgAnalysisBatchLFP displays the LFP power spectrum for drgRunBatch-generated
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
% 9 Compare auROC for power LFP for two events in two percent windows for all of the files
%
% 10 Compare auROC for power LFP for two groups (e.g. NRG1 vs control)
% within one precent window
%
% 11 Compare auROC for ERP LFP powerin between two percent correct windows
%
% 12 Justin's analysis of LFP power differences for naive and proficient
% mice
%
% 13 Compare auROC for ERP wavelet LFP power between two percent correct
%     windows at different ERP delta t values
%
% 14  Justin's analysis of LFP power differences for naive and proficient
% mice. Analyzed per mouse
%
% 15  Justin's fitlm analysis of LFP power differences for naive and proficient
% mice. Analyzed per mouse
%
% 16  Justin's fitglm analysis of LFP power differences for naive and proficient
% mice. Analyzed per mouse
%
% 17  Justin's PAC analysis for naive and proficient
% Analyzed per mouse

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
if isfield(handles_pars,'per_lab')
    per_lab=handles_pars.per_lab;
end

if ~isfield(handles_pars,'no_bandwidths')
    no_bandwidths=4;
    low_freq=[6 15 35 65];
    high_freq=[12 30 55 95];
    freq_names={'Theta','Beta','Low gamma','High gamma'};
else
    no_bandwidths=handles_pars.no_bandwidths;
    low_freq=handles_pars.low_freq;
    high_freq=handles_pars.high_freq;
    freq_names=handles_pars.freq_names;
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


fprintf(1, ['\ndrgDisplayBatchLFPPowerPairwise run for ' handles.drgb.outFileName '\nwhich_display= = %d\n\n'],which_display);

switch which_display
    
    case {1,6,7,8,11,13}
        
        frequency=handles_drgb.drgb.lfpevpair(1).fERP;
        max_events_per_sec=(handles_drgb.drgbchoices.timeEnd(winNo)-handles_drgb.drgbchoices.timeStart(winNo))*handles_drgb.max_events_per_sec;
    otherwise
        frequency=handles_drgb.drgb.freq_for_LFPpower;
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
        title('Percent singificant auROC')
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
                                            
                                            
                                            fprintf(1, ['%d trials in event No %d succesfully processed for file No %d electrode %d\n'],sum(trials_in_event_Ev), min_trials_per_event,files(fileNo),elec);
                                            
                                        else
                                            
                                            
                                            fprintf(1, ['%d trials in event No %d fewer than minimum trials per event ' evTypeLabels{evNo} ' for file No %d electrode %d\n'],sum(trials_in_event_Ev), min_trials_per_event,files(fileNo),elec);
                                            
                                            
                                        end
                                        
                                        
                                    end
                                    
                                    for evNo1=1:length(eventType)
                                        for evNo2=evNo1+1:length(eventType)
                                            if (noWB_for_evNo(evNo1)~=-1)&(noWB_for_evNo(evNo2)~=-1)
                                                
                                                for bwii=1:4
                                                    
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
        for bwii=1:4    %for bandwidths (theta, beta, low gamma, high gamma)
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
        for bwii=1:4
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
        hFig=figure(figNo)
        
        
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
        
        
        for bwii=1:4
            
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
        for bwii=1:4
            figNo=figNo+1;
            figure(figNo)
            hold on
            
            plot(sii_times(p(bwii,:)<=pFDRROC_chisq),(sig_per_proficient(bwii,p(bwii,:)<=pFDRROC_chisq)+sig_per_naive(bwii,p(bwii,:)<=pFDRROC_chisq))/2,'*k')
        end
        
        
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix]);
        
        
        
    case 14
        %Justin's per mouse analysis
        %For the proficient mice in the first and last sessions
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
        
        
        
        fprintf(1, ['Pairwise auROC analysis for Fig 1 of Justin''s paper\n\n'])
        p_vals=[];
        no_files=length(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        no_ROCs=0;
        szpc=size(percent_windows);
        for per_ii=1:szpc(1)
            
            
            
            for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
                for elec=1:16
                    if sum(which_electrodes==elec)>0
                        
                        
                        for evNo=1:length(eventType)
                            for bwii=1:4
                                theseEvNos(evNo,bwii).noEv=0;
                            end
                        end
                        
                        for fileNo=1:no_files
                            if sum(files==fileNo)>0
                                if handles_drgb.drgbchoices.mouse_no(fileNo)==mouseNo
                                    lfpodNo_ref=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                                    
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
                                                        
                                                        theseEvNos(evNo,bwii).this_delta_dB_powerEv(theseEvNos(evNo,bwii).noEv+1:theseEvNos(evNo,bwii).noEv+sum(trials_in_event_Ev))=this_delta_dB_powerEv;
                                                        
                                                        theseEvNos(evNo,bwii).noEv=theseEvNos(evNo,bwii).noEv+sum(trials_in_event_Ev);
                                                        
                                                        evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).noWB,bwii)=mean(this_delta_dB_powerEv);
                                                        %                                                     if evNo_out(evNo).mean_delta_dB_powerEvperBW(evNo_out(evNo).noWB,bwii)<-20
                                                        %                                                         fprintf(1, ['mean dB less than -20 dB in file no %d, electrode %d\n'],files(fileNo),elec);
                                                        %                                                     end
                                                        
                                                    end
                                                    
                                                    
                                                    fprintf(1, ['%d trials in event No %d succesfully processed for file No %d electrode %d\n'],sum(trials_in_event_Ev), min_trials_per_event,files(fileNo),elec);
                                                    
                                                else
                                                    
                                                    
                                                    fprintf(1, ['%d trials in event No %d fewer than minimum trials per event ' evTypeLabels{evNo} ' for file No %d electrode %d\n'],sum(trials_in_event_Ev), min_trials_per_event,files(fileNo),elec);
                                                    
                                                    
                                                end
                                                
                                                
                                            end
                                            
                                            
                                        else
                                            
                                            fprintf(1, ['Empty allPower for file No %d electrode %d\n'],files(fileNo),elec);
                                            
                                        end
                                        
                                        
                                    else
                                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],files(fileNo),elec);
                                        
                                        
                                    end
                                end %if mouseNo
                            end
                        end %fileNo
                        
                        
                        
                        can_calculate_auroc=1;
                        if can_calculate_auroc==1
                            for evNo1=1:length(eventType)
                                for evNo2=evNo1+1:length(eventType)
                                    if (noWB_for_evNo(evNo1)~=-1)&(noWB_for_evNo(evNo2)~=-1)
                                        
                                        for bwii=1:4
                                            
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
                                                if ((abs(evNo1-2)<=1)&(abs(evNo2-5)<=1))||((abs(evNo1-5)<=1)&(abs(evNo2-2)<=1))
                                                    ROC_between(no_ROCs)=1;
                                                else
                                                    ROC_between(no_ROCs)=0;
                                                end
                                                ROC_neighbor(no_ROCs)=abs(evNo1-evNo2);
                                                auROC(no_ROCs)=roc.AUC-0.5;
                                                p_valROC(no_ROCs)=roc.p;
                                                p_vals_ROC=[p_vals_ROC roc.p];
                                                
                                                if (per_ii==1)&(bwii==4)&(roc.AUC-0.5>0.3)
                                                    %This is here to stop and plot the ROC
                                                    %roc_out=roc_calc(roc_data);
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
            
            
        end
        fprintf(1, '\n\n')
        
        
        %         %Wide band plots are not useful, and have been commented out
        %
        %         %Calculate and plot the mean and 95% CI for each event
        %         figure(1)
        %         conc_anno_loc = {[0.15 0.15 0.2 0.2], [0.15 0.15 0.2 0.17], [0.15 0.15 0.2 0.14], [0.15 0.15 0.2 0.11], [0.15 0.15 0.2 0.08], [0.15 0.15 0.2 0.05]};
        %         for evNo=1:length(eventType)
        %             dB_Ev_ci=zeros(length(frequency),2);
        %             dB_Ev_mean=[];
        %             CI=[];
        %             dB_Ev_mean=mean(evNo_out(evNo).delta_dB_powerEvWB(evNo_out(evNo).per_ii==1,:));
        %             CI = bootci(1000, {@mean, evNo_out(evNo).delta_dB_powerEvWB(evNo_out(evNo).per_ii==1,:)},'type','cper');
        %             [hl1, hp1] = boundedline(frequency,dB_Ev_mean', CI', these_colors{evNo});
        %             annotation('textbox',conc_anno_loc{evNo},'String',evTypeLabels(evNo),'Color',these_colors{evNo},'EdgeColor','none');
        %             %             rectangle('Position',[0.15 0.15 0.2 0.2],'FaceColor','w','EdgeColor','k');
        %         end
        %
        %
        %         xlabel('Frequency (Hz)')
        %         ylabel('delta Power (dB)')
        %         ylim([-20 20]);
        %         title('Wideband spectrum proficient mice')
        %         set(gcf,'OuterPosition',[93 36 576 513]);
        %         set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        %
        %         %Calculate and plot the mean and 95% CI for each event
        %         figure(2)
        %         for evNo=1:length(eventType)
        %             dB_Ev_ci=zeros(length(frequency),2);
        %             dB_Ev_mean=[];
        %             dB_Ev_mean=mean(evNo_out(evNo).delta_dB_powerEvWB(evNo_out(evNo).per_ii==2,:));
        %             CI = bootci(1000, {@mean, evNo_out(evNo).delta_dB_powerEvWB(evNo_out(evNo).per_ii==2,:)},'type','cper');
        %             [hl1, hp1] = boundedline(frequency,dB_Ev_mean', CI', these_colors{evNo});
        %             annotation('textbox',conc_anno_loc{evNo},'String',evTypeLabels(evNo),'Color',these_colors{evNo},'EdgeColor','none');
        %             rectangle('Position',[0.15 0.15 0.2 0.2],'FaceColor','w','EdgeColor','k');
        %         end
        %
        %         xlabel('Frequency (Hz)')
        %         ylabel('delta Power (dB)')
        %         ylim([-5 10]);
        %         title('Wideband spectrum naive mice')
        %         set(gcf,'OuterPosition',[93 550 576 513]);
        %         set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        %         ylim([-20 20])
        %
        %         fig_pos = {[664 550 576 513],[1233 550 576 513],[664 36 576 513],[1233 36 576 513]};
        %
        
        %Now plot the histograms and the average
        conc_anno_loc = {[0.15 0.15 0.2 0.2], [0.15 0.15 0.2 0.17], [0.15 0.15 0.2 0.14], [0.15 0.15 0.2 0.11], [0.15 0.15 0.2 0.08], [0.15 0.15 0.2 0.05]};
        fig_pos = {[664 550 576 513],[1233 550 576 513],[664 36 576 513],[1233 36 576 513]};
        for bwii=1:4    %for bandwidths (theta, beta, low gamma, high gamma)
            %Plot the average
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
            title([freq_names{bwii} ' average delta dB power per electrode'])
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
        
        %Plot cumulative histos for auROCs within vs between S+ and S-
        for bwii=1:4
            
            
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            hold on
            
            %Naive between
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1)));
            plot(x_aic,f_aic,'b')
            
            %Proficient between
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1)));
            plot(x_aic,f_aic,'r')
            
            %Naive within
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0)));
            plot(x_aic,f_aic,'Color',[0.8 0.8 1])
            
            %Proficient within
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0)));
            plot(x_aic,f_aic,'Color',[1 0.8 0.8])
            
            legend('Naive between','Proficient between','Naive within','Proficient within')
            xlabel('auROC')
            ylabel('Cumulative probability')
            title(freq_names{bwii})
            
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
            fprintf(1, ['p value for anovan auROC for naive vs proficient for ' freq_names{bwii} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan auROC for within vs between for ' freq_names{bwii} '= %d \n'],  p(2));
            p_anova_np(bwii)=p(1);
            p_anova_wb(bwii)=p(2);
        end
        
        %Plot cumulative histos for auROCs within vs between S+ and S- using
        %only ROCs for adjacent concentrations
        for bwii=1:4
            
            
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            hold on
            
            %Naive between
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==1)));
            plot(x_aic,f_aic,'b')
            
            %Proficient between
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==1)&(ROC_neighbor==1)));
            plot(x_aic,f_aic,'r')
            
            %Naive within
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==1)));
            plot(x_aic,f_aic,'Color',[0.8 0.8 1])
            
            %Proficient within
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCbandwidth==bwii)&(ROC_between==0)&(ROC_neighbor==1)));
            plot(x_aic,f_aic,'Color',[1 0.8 0.8])
            
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
        for bwii=1:4
            
            
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
        
        for bwii=1:4
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
        hFig=figure(figNo)
        
        
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
        hFig=figure(figNo)
        
        for bwii=1:4
            subplot(2,2,bwii)
            hold on
            %Plot within naive
            bar(1,100*no_sig_within(2,bwii)/no_within(2,bwii),'b')
            bar(2,100*no_sig_within1(2,bwii)/no_within1(2,bwii),'b')
            bar(3,100*no_sig_within2(2,bwii)/no_within2(2,bwii),'b')
            
            %Plot within proficient
            bar(4,100*no_sig_within(1,bwii)/no_within(1,bwii),'r')
            bar(5,100*no_sig_within1(1,bwii)/no_within1(1,bwii),'r')
            bar(6,100*no_sig_within2(1,bwii)/no_within2(1,bwii),'r')
            
            %Plot between naive
            bar(9,100*no_sig_between(2,bwii)/no_between(2,bwii),'b')
            bar(10,100*no_sig_between1(2,bwii)/no_between1(2,bwii),'b')
            bar(11,100*no_sig_between2(2,bwii)/no_between2(2,bwii),'b')
            
            %Plot between proficient
            bar(12,100*no_sig_between(1,bwii)/no_between(1,bwii),'r')
            bar(13,100*no_sig_between1(1,bwii)/no_between1(1,bwii),'r')
            bar(14,100*no_sig_between2(1,bwii)/no_between2(1,bwii),'r')
            
            title(['Percent significant auROC ' freq_names{bwii} ])
        end
        
        ppno=0;
        for pair_no1=1:no_pairs
            for pair_no2=pair_no1+1:no_pairs
                ppno=ppno+1;
                [pChiSq(ppno), Q]= chi2test([sig(pair_no1), not_sig(pair_no1); sig(pair_no2), not_sig(pair_no2)]);
                fprintf(1, ['pchi  = %d\n'],pChiSq(ppno));
            end
        end
        
        try
            pFDRchisq=drsFDRpval(pChiSq);
            fprintf(1, ['pFDR for ChiSquared  = %d\n\n'],pFDRchisq);
        catch
        end
        
        %We should add a bar graph to compare percent significant ROCs for
        %within vs between (with/wo adjacency requirement)
        
        pfft=1
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix]);
        
    case 15
        %Justin
        %Linear fit of delta power to: spm, concentration, previous reward, percent
        %This does the analysis in all the files and DOES not distinguish between groups!!!
        
        
        
        
        no_fitlms=0;
        fitlm_results=[];
        
        
        
        fprintf(1, ['fitlm analysis for Justin''s paper\n\n'])
        
        no_files=length(files);
        
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
                                        lfpodNo_ref=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                                        
                                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref)))
                                            
                                            
                                            if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower))
                                                
                                                percent_mask=[];
                                                percent_mask=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower>=percent_windows(per_ii,1))...
                                                    &(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower<=percent_windows(per_ii,2));
                                                
                                                for evNo=1:length(eventType)
                                                    
                                                    trials_in_event_Ev=[];
                                                    trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(eventType(evNo),:)==1)&percent_mask;
                                                    
                                                    if (sum(trials_in_event_Ev)>=1)
                                                        
                                                        lfpodNo=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                                        
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
                                                
                                                fprintf(1, ['Empty allPower for file No %d electrode %d\n'],files(fileNo),elec);
                                                
                                                
                                            end
                                            
                                            
                                        else
                                            fprintf(1, ['Empty LFP reference for file No %d electrode %d\n'],files(fileNo),elec);
                                            
                                            
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
            for bwii=1:4
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
            
            
            
            for bwii=1:4
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
        
        no_files=length(files);
        
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
                                        lfpodNo_ref=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                                        
                                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref)))
                                            
                                            
                                            if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo_ref).allPower))
                                                
                                                percent_mask=[];
                                                percent_mask=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower>=percent_windows(per_ii,1))...
                                                    &(handles_drgb.drgb.lfpevpair(lfpodNo_ref).perCorrLFPPower<=percent_windows(per_ii,2));
                                                
                                                for evNo=1:length(eventType)
                                                    
                                                    trials_in_event_Ev=[];
                                                    trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfpodNo_ref).which_eventLFPPower(eventType(evNo),:)==1)&percent_mask;
                                                    
                                                    if (sum(trials_in_event_Ev)>=1)
                                                        
                                                        lfpodNo=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                                        
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
                                                
                                                fprintf(1, ['Empty allPower for file No %d electrode %d\n'],files(fileNo),elec);
                                                
                                                
                                            end
                                            
                                            
                                        else
                                            fprintf(1, ['Empty LFP reference for file No %d electrode %d\n'],files(fileNo),elec);
                                            
                                            
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
            for bwii=1:4
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
            
            
            
            for bwii=1:4
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
        %Justin's per mouse analysis for PAC
        %For the proficient mice in the first and last sessions
        %plot the LFP spectrum for S+ vs S-, plot LFP power for S+ vs S- for each electrode and plot auROCs
        %NOTE: This does the analysis in all the files and DOES not distinguish between groups!!!
        no_ROCs=0;
        ROCout=[];
        p_vals_ROC=[];
        mean_MI_No=0;
        mean_MI=[];
        mean_MI_perii=[];
        
        
        
        fprintf(1, ['Pairwise auROC analysis for Fig 1 of Justin''s paper\n\n'])
        p_vals=[];
        no_files=length(files);
        
        if exist('which_electrodes')==0
            which_electrodes=[1:16];
        end
        
        no_ROCs=0;
        szpc=size(percent_windows);
        for per_ii=1:szpc(1)
            
            
            
            for mouseNo=1:max(handles_drgb.drgbchoices.mouse_no)
                for elec=1:16
                    if sum(which_electrodes==elec)>0
                        
                        
                        for evNo=1:length(eventType)
                            for pacii=1:3
                                theseEvNos(evNo,pacii).noEv=0;
                            end
                        end
                        
                        for fileNo=1:no_files
                            if sum(files==fileNo)>0
                                if handles_drgb.drgbchoices.mouse_no(fileNo)==mouseNo
                                    lfpodNo=find((files_per_lfp==files(fileNo))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                    
                                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo)))
                                        
                                        
                                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNo).PAC))
                                            
                                            percent_mask=[];
                                            percent_mask=(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrLFPPower>=percent_windows(per_ii,1))...
                                                &(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrLFPPower<=percent_windows(per_ii,2));
                                            
                                            
                                            for evNo=1:length(eventType)
                                                
                                                noWB_for_evNo(evNo)=-1;
                                                
                                                trials_in_event_Ev=[];
                                                trials_in_event_Ev=(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventLFPPower(eventType(evNo),:)==1)&percent_mask;
                                                
                                                if (sum(trials_in_event_Ev)>=1)
                                                    
                                                    %Do per bandwidth analysis
                                                    for pacii=1:3
                                                        
                                                        %Enter the modulation index
                                                        this_MI_Ev=zeros(sum(trials_in_event_Ev),1);
                                                        this_MI_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).mod_indx(trials_in_event_Ev);
                                                        theseEvNos(evNo,pacii).this_MI_Ev(theseEvNos(evNo,pacii).noEv+1:theseEvNos(evNo,pacii).noEv+sum(trials_in_event_Ev))=this_MI_Ev;
                                                        mean_MI_No=mean_MI_No+1;
                                                        mean_MI(mean_MI_No)=mean(this_MI_Ev);
                                                        mean_MI_perii(mean_MI_No)=per_ii;
                                                        mean_MI_evNo(mean_MI_No)=evNo;
                                                        mean_MI_pacii(mean_MI_No)=pacii;
                                                        
                                                        if mean_MI(mean_MI_No)>=0.035
                                                            fprintf(1, ['MI larger than 0.035 for mouse no %d, file no %d, electrode, %d, pac no %d, perii %d, conc, %d\n'],mouseNo, fileNo, elec, pacii, per_ii, evNo);
                                                        end
                                                        
                                                        %Enter the meanVectorLength
                                                        this_meanVectorLength_Ev=zeros(sum(trials_in_event_Ev),1);
                                                        this_meanVectorLength_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).meanVectorLength(trials_in_event_Ev);
                                                        theseEvNos(evNo,pacii).this_meanVectorLength_Ev(theseEvNos(evNo,pacii).noEv+1:theseEvNos(evNo,pacii).noEv+sum(trials_in_event_Ev))=this_meanVectorLength_Ev;
                                                        
                                                        %Enter the meanVectorAngle
                                                        this_meanVectorAngle_Ev=zeros(sum(trials_in_event_Ev),1);
                                                        this_meanVectorAngle_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).meanVectorAngle(trials_in_event_Ev);
                                                        theseEvNos(evNo,pacii).this_meanVectorAngle_Ev(theseEvNos(evNo,pacii).noEv+1:theseEvNos(evNo,pacii).noEv+sum(trials_in_event_Ev))=this_meanVectorAngle_Ev;
                                                        
                                                        %Enter the peakAngle
                                                        this_peakAngle_Ev=zeros(sum(trials_in_event_Ev),1);
                                                        this_peakAngle_Ev=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(pacii).peakAngle(trials_in_event_Ev);
                                                        theseEvNos(evNo,pacii).this_peakAngle_Ev(theseEvNos(evNo,pacii).noEv+1:theseEvNos(evNo,pacii).noEv+sum(trials_in_event_Ev))=this_peakAngle_Ev;
                                                        
                                                        theseEvNos(evNo,pacii).noEv=theseEvNos(evNo,pacii).noEv+sum(trials_in_event_Ev);
                                                    end
                                                    
                                                    
                                                    fprintf(1, ['%d trials in event No %d succesfully processed for file No %d electrode %d\n'],sum(trials_in_event_Ev), min_trials_per_event,files(fileNo),elec);
                                                    
                                                else
                                                    
                                                    
                                                    fprintf(1, ['%d trials in event No %d fewer than minimum trials per event ' evTypeLabels{evNo} ' for file No %d electrode %d\n'],sum(trials_in_event_Ev), min_trials_per_event,files(fileNo),elec);
                                                    
                                                    
                                                end
                                                
                                                
                                            end
                                            
                                            
                                        else
                                            
                                            fprintf(1, ['Empty allPower for file No %d electrode %d\n'],files(fileNo),elec);
                                            
                                        end
                                        
                                        
                                    else
                                        fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],files(fileNo),elec);
                                        
                                        
                                    end
                                end %if mouseNo
                            end
                        end %fileNo
                        
                        
                        
                        
                        for evNo1=1:length(eventType)
                            for evNo2=evNo1+1:length(eventType)
                                
                                for pacii=1:3
                                    
                                    %Enter Ev1
                                    trials_in_event_Ev1=length(theseEvNos(evNo1,pacii).this_MI_Ev);
                                    this_MI_Ev1=zeros(trials_in_event_Ev1,1);
                                    this_MI_Ev1=theseEvNos(evNo1,pacii).this_MI_Ev;
                                    roc_data=[];
                                    roc_data(1:sum(trials_in_event_Ev1),1)=this_MI_Ev1;
                                    roc_data(1:sum(trials_in_event_Ev1),2)=zeros(sum(trials_in_event_Ev1),1);
                                    
                                    %Enter Ev2
                                    trials_in_event_Ev2=length(theseEvNos(evNo2,pacii).this_MI_Ev);
                                    total_trials=trials_in_event_Ev1+trials_in_event_Ev2;
                                    this_MI_Ev2=zeros(trials_in_event_Ev2,1);
                                    this_MI_Ev2=theseEvNos(evNo2,pacii).this_MI_Ev;
                                    roc_data(sum(trials_in_event_Ev1)+1:total_trials,1)=this_MI_Ev2;
                                    roc_data(sum(trials_in_event_Ev1)+1:total_trials,2)=ones(sum(trials_in_event_Ev2),1);
                                    
                                    
                                    %Find  ROC
                                    if (trials_in_event_Ev1>=5)&(trials_in_event_Ev2>=5)
                                        no_ROCs=no_ROCs+1;
                                        roc=roc_calc(roc_data,0,0.05,0);
                                        ROCout(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
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
                                        
                                        if (per_ii==1)&(pacii==4)&(roc.AUC-0.5>0.3)
                                            %This is here to stop and plot the ROC
                                            %roc_out=roc_calc(roc_data);
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
        fprintf(1, ['pFDR for auROC  = %d\n\n'],pFDRauROC);
        
        %         %Wide band plots are not useful, and have been commented out
        %
        %         %Calculate and plot the mean and 95% CI for each event
        %         figure(1)
        %         conc_anno_loc = {[0.15 0.15 0.2 0.2], [0.15 0.15 0.2 0.17], [0.15 0.15 0.2 0.14], [0.15 0.15 0.2 0.11], [0.15 0.15 0.2 0.08], [0.15 0.15 0.2 0.05]};
        %         for evNo=1:length(eventType)
        %             dB_Ev_ci=zeros(length(frequency),2);
        %             dB_Ev_mean=[];
        %             CI=[];
        %             dB_Ev_mean=mean(evNo_out(evNo).delta_dB_powerEvWB(evNo_out(evNo).per_ii==1,:));
        %             CI = bootci(1000, {@mean, evNo_out(evNo).delta_dB_powerEvWB(evNo_out(evNo).per_ii==1,:)},'type','cper');
        %             [hl1, hp1] = boundedline(frequency,dB_Ev_mean', CI', these_colors{evNo});
        %             annotation('textbox',conc_anno_loc{evNo},'String',evTypeLabels(evNo),'Color',these_colors{evNo},'EdgeColor','none');
        %             %             rectangle('Position',[0.15 0.15 0.2 0.2],'FaceColor','w','EdgeColor','k');
        %         end
        %
        %
        %         xlabel('Frequency (Hz)')
        %         ylabel('delta Power (dB)')
        %         ylim([-20 20]);
        %         title('Wideband spectrum proficient mice')
        %         set(gcf,'OuterPosition',[93 36 576 513]);
        %         set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        %
        %         %Calculate and plot the mean and 95% CI for each event
        %         figure(2)
        %         for evNo=1:length(eventType)
        %             dB_Ev_ci=zeros(length(frequency),2);
        %             dB_Ev_mean=[];
        %             dB_Ev_mean=mean(evNo_out(evNo).delta_dB_powerEvWB(evNo_out(evNo).per_ii==2,:));
        %             CI = bootci(1000, {@mean, evNo_out(evNo).delta_dB_powerEvWB(evNo_out(evNo).per_ii==2,:)},'type','cper');
        %             [hl1, hp1] = boundedline(frequency,dB_Ev_mean', CI', these_colors{evNo});
        %             annotation('textbox',conc_anno_loc{evNo},'String',evTypeLabels(evNo),'Color',these_colors{evNo},'EdgeColor','none');
        %             rectangle('Position',[0.15 0.15 0.2 0.2],'FaceColor','w','EdgeColor','k');
        %         end
        %
        %         xlabel('Frequency (Hz)')
        %         ylabel('delta Power (dB)')
        %         ylim([-5 10]);
        %         title('Wideband spectrum naive mice')
        %         set(gcf,'OuterPosition',[93 550 576 513]);
        %         set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
        %         ylim([-20 20])
        %
        %         fig_pos = {[664 550 576 513],[1233 550 576 513],[664 36 576 513],[1233 36 576 513]};
        %
        
        %Now plot the average MI
        
        conc_anno_loc = {[0.15 0.15 0.2 0.2], [0.15 0.15 0.2 0.17], [0.15 0.15 0.2 0.14], [0.15 0.15 0.2 0.11], [0.15 0.15 0.2 0.08], [0.15 0.15 0.2 0.05]};
        fig_pos = {[664 550 576 513],[1233 550 576 513],[664 36 576 513],[1233 36 576 513]};
        for pacii=1:3    %for amplitude bandwidths (beta, low gamma, high gamma)
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
                end
            end
            title(['Average MI per electrode for PAC theta/' freq_names{pacii+1}])
            set(gcf,'OuterPosition',fig_pos{pacii});
            bar_lab_loc = [3.5 6.5 9.5 12.5 15.5 18.5];
            xticks(bar_lab_loc)
            xticklabels(concs2)
            xlabel('Concentration (%)')
            ylabel('Modulation Index')
        end
        
        
        
        
        %Plot cumulative histos for auROCs within vs between S+ and S-
        for pacii=1:3
            
            
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            hold on
            
            %Naive between
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1)));
            plot(x_aic,f_aic,'b')
            
            %Proficient between
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1)));
            plot(x_aic,f_aic,'r')
            
            %Naive within
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0)));
            plot(x_aic,f_aic,'Color',[0.8 0.8 1])
            
            %Proficient within
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0)));
            plot(x_aic,f_aic,'Color',[1 0.8 0.8])
            
            legend('Naive between','Proficient between','Naive within','Proficient within')
            xlabel('auROC')
            ylabel('Cumulative probability')
            title(freq_names{pacii})
            
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
            fprintf(1, ['p value for anovan auROC for naive vs proficient for ' freq_names{pacii+1} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan auROC for within vs between for ' freq_names{pacii+1} '= %d \n'],  p(2));
            p_anova_np(pacii)=p(1);
            p_anova_wb(pacii)=p(2);
        end
        
        %Plot cumulative histos for auROCs within vs between S+ and S- using
        %only ROCs for adjacent concentrations
        for pacii=1:3
            
            
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            hold on
            
            %Naive between
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==1)));
            plot(x_aic,f_aic,'b')
            
            %Proficient between
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==1)&(ROC_neighbor==1)));
            plot(x_aic,f_aic,'r')
            
            %Naive within
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==2)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==1)));
            plot(x_aic,f_aic,'Color',[0.8 0.8 1])
            
            %Proficient within
            [f_aic,x_aic] = drg_ecdf(auROC((ROCper_ii==1)&(ROCpacii==pacii)&(ROC_between==0)&(ROC_neighbor==1)));
            plot(x_aic,f_aic,'Color',[1 0.8 0.8])
            
            legend('Naive between','Proficient between','Naive within','Proficient within')
            xlabel('auROC')
            ylabel('Cumulative probability')
            title([freq_names{pacii} ' Adjacent concentrations'])
            
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
            fprintf(1, ['p value for anovan MI auROC adjacent concentrations for naive vs proficient for PAC theta/' freq_names{pacii} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan MI auROC adjacent concentrations  for within vs between for PAC theta/' freq_names{pacii} '= %d \n'],  p(2));
            p_anova_np_adj(pacii)=p(1);
            p_anova_wb_adj(pacii)=p(2);
        end
        
        %Plot cumulative histos for auROCs within vs between S+ and S- using
        %only ROCs for concentrations separated by two log steps
        for pacii=1:3
            
            
            figNo = get(gcf,'Number')+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
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
            title([freq_names{pacii} ' concentrations separated by two log steps'])
            
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
            fprintf(1, ['p value for anovan MI auROC two log step concentrations for naive vs proficient for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(1));
            fprintf(1, ['p value for anovan MI auROC two log step concentrations for within vs between for PAC theta/' freq_names{pacii+1} '= %d \n'],  p(2));
            p_anova_np_adj(pacii)=p(1);
            p_anova_wb_adj(pacii)=p(2);
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
        
        for pacii=1:3
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
                
                title(['% significant MI auROC for PAC theta/' freq_names{pacii+1} ' ' per_lab(pcii)])
                set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2)
                
                pfft=1
            end
        end
        
        
        
        
        figNo = get(gcf,'Number')+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo)
        
        
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
        hFig=figure(figNo)
        
        for pacii=1:3
            subplot(2,2,pacii)
            hold on
            %Plot within naive
            bar(1,100*no_sig_within(2,pacii)/no_within(2,pacii),'b')
            bar(2,100*no_sig_within1(2,pacii)/no_within1(2,pacii),'b')
            bar(3,100*no_sig_within2(2,pacii)/no_within2(2,pacii),'b')
            
            %Plot within proficient
            bar(4,100*no_sig_within(1,pacii)/no_within(1,pacii),'r')
            bar(5,100*no_sig_within1(1,pacii)/no_within1(1,pacii),'r')
            bar(6,100*no_sig_within2(1,pacii)/no_within2(1,pacii),'r')
            
            %Plot between naive
            bar(9,100*no_sig_between(2,pacii)/no_between(2,pacii),'b')
            bar(10,100*no_sig_between1(2,pacii)/no_between1(2,pacii),'b')
            bar(11,100*no_sig_between2(2,pacii)/no_between2(2,pacii),'b')
            
            %Plot between proficient
            bar(12,100*no_sig_between(1,pacii)/no_between(1,pacii),'r')
            bar(13,100*no_sig_between1(1,pacii)/no_between1(1,pacii),'r')
            bar(14,100*no_sig_between2(1,pacii)/no_between2(1,pacii),'r')
            
            title(['% significant MI auROC for theta/' freq_names{pacii+1} ' PAC'])
        end
        
        %         ppno=0;
        %         for pair_no1=1:no_pairs
        %             for pair_no2=pair_no1+1:no_pairs
        %                 ppno=ppno+1;
        %                 [pChiSq(ppno), Q]= chi2test([sig(pair_no1), not_sig(pair_no1); sig(pair_no2), not_sig(pair_no2)]);
        %                 fprintf(1, ['pchi  = %d\n'],pChiSq(ppno));
        %             end
        %         end
        
        %         try
        %             pFDRchisq=drsFDRpval(pChiSq);
        %             fprintf(1, ['pFDR for ChiSquared  = %d\n\n'],pFDRchisq);
        %         catch
        %         end
        
        %We should add a bar graph to compare percent significant ROCs for
        %within vs between (with/wo adjacency requirement)
        
        pfft=1
        save([handles.PathName handles.drgb.outFileName(1:end-4) output_suffix]);
        
  
end


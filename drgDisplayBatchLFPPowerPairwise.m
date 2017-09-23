function drgDisplayBatchLFPPowerPairwise(handles)

%This function displays the LFP power spectrum for drgRunBatch-generated
%data

%Which analysis is performed is determined by the value enterd in the
%variable which_display:
%
%
% 1 Find the electrodes displaying a post vs. pre difference
%
% 2 auROC analysis for post vs. pre differences
%
% 3 Generate figure 1 for Daniels' LFP power paper. For the proficient mice in the firstlast LFP batch analysis
%   Plot the spectrum for S+ vs S-, for each bandwidth plot LFP power for S+ vs S- for each electrode
%   and plot auROCs
%
%
% 4 Generate figure 2 for Daniel's paper. first vs last.




%% First enter the choices for what will be analyzed.
% THESE VALUES ARE IMPORTANT and differ for each user



%
% % For Daniel's APEB tstart
% winNo=2;
% refWin=1;
% which_display=2;
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% eventType=[1 1];
% evTypeLabels={'tstart','tstart'};
%
% %Experiment pairs
% %Important: The first file must be the experiment performed first
% %For example in acetophenone ethyl benzoate no laser is first, laser is
% %second
% file_pairs=[
%     1 7;
%     2 8;
%     3 9
%     4 10;
%     5 11;
%     6 12;
%     13 14;
%     15 16;
%     21 17;
%     23 18;
%     22 19;
%     24 20];
% no_file_pairs=12;
%
% comp_window=10; %works well with 8-12
% comp_window_auROC=20;
%
% grpre=[1 3];
% grpost=[2 4];


% % For Daniel's tabaproprio
% winNo=2;
% refWin=1;
% which_display=2;
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
%
% %Experiment pairs
% %Important: The first file must be the experiment performed first
% %For example in acetophenone ethyl benzoate no laser is first, laser is
% %second
% file_pairs=[
%     1 8;
%     2 9;
%     3 10;
%     4 11;
%     5 12;
%     6 13;
%     7 14];
% no_file_pairs=7;
%
% comp_window=10; %works well with 8-12
% comp_window_auROC=20;
%
% grpre=[1 3];
% grpost=[2 4];


% % For Daniel's isoamyl acetate tstart
winNo=2;
refWin=1;
which_display=5;
% eventType=[2 5];
% evTypeLabels={'Hit','CR'};
eventType=[3 6];
evTypeLabels={'S+','S-'};

%Experiment pairs
%Important: The first file must be the experiment performed first
%For example in acetophenone ethyl benzoate no laser is first, laser is
%second
file_pairs=[
    1 13;
    2 11;
    3 12;
    4 14;
    5 15;
    6 16;
    7 17;
    8 18;
    9 20;
    10 21;
    22 19];
no_file_pairs=11;

comp_window=10; %works well with 8-12
comp_window_auROC=20;

grpre=[1 3];
grpost=[2 4];

% % For Daniel's acetophenone ethyl benzoate
% winNo=2;
% refWin=1;
% which_display=2;
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
%
%
% %Experiment pairs
% %Important: The first file must be the experiment performed first
% %For example in acetophenone ethyl benzoate no laser is first, laser is
% %second
% file_pairs=[
%     1 7;
%     2 8;
%     3 9;
%     4 10;
%     5 11;
%     6 12;
%     13 14;
%     15 16;
%     21 17;
%     23 18;
%     22 19;
%     24 20];
% no_file_pairs=12;
%
% comp_window=15; %works well with 8-12
% comp_window_auROC=30;
%
% grpre=[1 3];
% grpost=[2 4];

% % % For Daniel's APEBEfirstandlast91117
% winNo=2;
% refWin=1;
% which_display=3;
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
%
%
% %Experiment pairs
% %Important: The first file must be the experiment performed first
% %For example in acetophenone ethyl benzoate no laser is first, laser is
% %second
% file_pairs=[
%     1 11;
%     2 7;
%     3 8;
%     4 9;
%     5 10;
%     6 12;
%     13 16;
%     14 17;
%     15 18];
% no_file_pairs=9;
%
% comp_window=15; %works well with 8-12
% comp_window_auROC=30;
%
% grpre=[1 3];
% grpost=[2 4];


% % For Daniel's Fig. 1 run with isomin_firstandlastIAMO91817
% % and which_display=3
% % For Daniel's Fig. 2 run with isomin_firstandlastIAMO91817
% % and which_display=4
% winNo=2;
% refWin=1;
% which_display=3;
% % eventType=[2 5];
% % evTypeLabels={'Hit','CR'};
% eventType=[3 6];
% evTypeLabels={'S+','S-'};
%
%
% %Experiment pairs
% %Important: The first file must be the experiment performed first
% %For example in acetophenone ethyl benzoate no laser is first, laser is
% %second
% file_pairs=[
%     1 5;
%     2 4;
%     3 6;
%     7 11;
%     9 12;
%     10 13];
% no_file_pairs=6;
%
% comp_window=15; %works well with 8-12
% comp_window_auROC=30;
%
% grpre=[1 3];
% grpost=[2 4];


%% The code processing pairwise batch LFP starts here

close all
warning('off')

no_event_types=length(eventType);

%Statistics
%stats_method = 1, FDR, 2, statscond
stat_method=1;

%Mode for percent ANOVA
mode_statcond='perm';
% mode_statcond='bootstrap';

%Which percent correct bins do you want to use?
percent_low=[45 65 80];
percent_high=[65 80 100];
percent_bin_legend={' learning';' between learning and proficient';' proficient'};
no_percent_bins=3;

%Bandwidths
no_bandwidths=4;
low_freq=[6 15 35 65];
high_freq=[12 30 55 95];
freq_names={'Theta','Beta','Low gamma','High gamma'};

evHit=eventType(1);
evCR=eventType(2);


%Ask user for the drgb output .mat file and load those data
[handles.drgb.outFileName,handles.PathName] = uigetfile('*.mat','Select the drgb output file');
load([handles.PathName handles.drgb.outFileName])



frequency=handles_drgb.drgb.freq_for_LFPpower;



figNo=0;

%These are the colors for the different lines

these_colors{1}='b';
these_colors{2}='r';
these_colors{3}='m';
these_colors{4}='g';
these_colors{5}='y';
these_colors{6}='k';
these_colors{7}='c';
these_colors{8}='k';

these_lines{1}='-b';
these_lines{2}='-r';
these_lines{3}='-m';
these_lines{4}='-g';
these_lines{5}='-y';
these_lines{6}='-k';
these_lines{7}='-c';
these_lines{8}='-k';

%Initialize the variables
%Get files and electrode numbers
for lfpodNo=1:handles_drgb.drgb.lfpevpair_no
    files_per_lfp(lfpodNo)=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
    elec_per_lfp(lfpodNo)=handles_drgb.drgb.lfpevpair(lfpodNo).elecNo;
    window_per_lfp(lfpodNo)=handles_drgb.drgb.lfpevpair(lfpodNo).timeWindow;
end

switch which_display
    case 1
        %For each electrode determine whether the LFP power differs between the last few trials of pre with first few trials of post
        
        no_dBs=1;
        delta_dB_power_pre=[];
        no_ROCs=0;
        ROCoutpre=[];
        ROCoutpost=[];
        p_vals_ROC=[];
        delta_dB_powerpreHit=[];
        no_hits=0;
        
        fprintf(1, ['Comparison of LFP power between pre and post\n\n'])
        p_vals=[];
        for fps=1:no_file_pairs
            for elec=1:16
                
                lfpodNopre_ref=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                lfpodNopost_ref=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                
                if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref)))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref)))
                    
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower))
                        
                        
                        trials_in_event_preHit=(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).which_eventLFPPower(evHit,:)==1);
                        trials_in_event_preCR=(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).which_eventLFPPower(evCR,:)==1);
                        trials_in_event_pre=logical(trials_in_event_preHit+trials_in_event_preCR);
                        
                        trials_in_event_postHit=(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).which_eventLFPPower(evHit,:)==1);
                        trials_in_event_postCR=(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).which_eventLFPPower(evCR,:)==1);
                        trials_in_event_post=logical(trials_in_event_postHit+trials_in_event_postCR);
                        
                        if (sum(trials_in_event_preHit)>=comp_window)&(sum(trials_in_event_preCR)>=comp_window)&(sum(trials_in_event_postHit)>=comp_window)&(sum(trials_in_event_postCR)>=comp_window)
                            if (sum(trials_in_event_pre)>=comp_window_auROC)&(sum(trials_in_event_post)>=comp_window_auROC)
                                lfpodNopre=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                lfpodNopost=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                
                                %pre Hits
                                
                                
                                no_trials=0;
                                jj=length(trials_in_event_preHit);
                                while no_trials<comp_window
                                    if trials_in_event_preHit(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    jj=jj-1;
                                end
                                jj=jj+1;
                                these_trials=logical([zeros(1,jj-1) ones(1,length(trials_in_event_preHit)-jj+1)].*trials_in_event_preHit);
                                
                                this_dB_powerprerefHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerprerefHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower(these_trials,:));
                                
                                this_dB_powerpreHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpreHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre).allPower(these_trials,:));
                                
                                
                                %pre CRs
                                no_trials=0;
                                jj=length(trials_in_event_preCR);
                                while no_trials<comp_window
                                    if trials_in_event_preCR(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    jj=jj-1;
                                end
                                jj=jj+1;
                                these_trials=logical([zeros(1,jj-1) ones(1,length(trials_in_event_preCR)-jj+1)].*trials_in_event_preCR);
                                
                                this_dB_powerprerefCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerprerefCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower(these_trials,:));
                                
                                this_dB_powerpreCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpreCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre).allPower(these_trials,:));
                                
                                
                                %post Hits
                                no_trials=0;
                                jj=0;
                                while no_trials<comp_window
                                    jj=jj+1;
                                    if trials_in_event_postHit(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    
                                end
                                
                                these_trials=logical([ones(1,jj) zeros(1,length(trials_in_event_postHit)-jj)].*trials_in_event_postHit);
                                
                                this_dB_powerpostrefHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostrefHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower(these_trials,:));
                                
                                this_dB_powerpostHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost).allPower(these_trials,:));
                                
                                
                                %post CRs
                                no_trials=0;
                                jj=0;
                                while no_trials<comp_window
                                    jj=jj+1;
                                    if trials_in_event_postCR(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    
                                end
                                
                                these_trials=logical([ones(1,jj) zeros(1,length(trials_in_event_postCR)-jj)].*trials_in_event_postCR);
                                
                                this_dB_powerpostrefCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostrefCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower(these_trials,:));
                                
                                this_dB_powerpostCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost).allPower(these_trials,:));
                                
                                for bwii=1:no_bandwidths
                                    
                                    no_ROCs=no_ROCs+1;
                                    this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                    
                                    %Enter the pre Hits
                                    this_delta_dB_powerpreHit=zeros(comp_window,1);
                                    this_delta_dB_powerpreHit=mean(this_dB_powerpreHit(:,this_band)-this_dB_powerprerefHit(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:comp_window,1)=this_delta_dB_powerpreHit;
                                    roc_data(1:comp_window,2)=zeros(comp_window,1);
                                    
                                    %Enter pre CR
                                    this_delta_dB_powerpreCR=zeros(comp_window,1);
                                    this_delta_dB_powerpreCR=mean(this_dB_powerpreCR(:,this_band)-this_dB_powerprerefCR(:,this_band),2);
                                    roc_data(comp_window+1:2*comp_window,1)=this_delta_dB_powerpreCR;
                                    roc_data(comp_window+1:2*comp_window,2)=ones(comp_window,1);
                                    
                                    
                                    %Find pre ROC
                                    ROCoutpre(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCoutpre(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopre_ref).fileNo;
                                    ROCgroupNopre(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).fileNo);
                                    ROCoutpre(no_ROCs).timeWindow=winNo;
                                    ROCbandwidthpre(no_ROCs)=bwii;
                                    auROCpre(no_ROCs)=ROCoutpre(no_ROCs).roc.AUC-0.5;
                                    
                                    p_vals_ROC=[p_vals_ROC ROCoutpre(no_ROCs).roc.p];
                                    
                                    %Enter the post Hits
                                    this_delta_dB_powerpostHit=zeros(comp_window,1);
                                    this_delta_dB_powerpostHit=mean(this_dB_powerpostHit(:,this_band)-this_dB_powerpostrefHit(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:comp_window,1)=this_delta_dB_powerpostHit;
                                    roc_data(1:comp_window,2)=zeros(comp_window,1);
                                    
                                    %Enter post CR
                                    this_delta_dB_powerpostCR=zeros(comp_window,1);
                                    this_delta_dB_powerpostCR=mean(this_dB_powerpostCR(:,this_band)-this_dB_powerpostrefCR(:,this_band),2);
                                    roc_data(comp_window+1:2*comp_window,1)=this_delta_dB_powerpostCR;
                                    roc_data(comp_window+1:2*comp_window,2)=ones(comp_window,1);
                                    
                                    
                                    %Find post ROC
                                    ROCoutpost(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCoutpost(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopost_ref).fileNo;
                                    ROCgroupNopost(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).fileNo);
                                    ROCoutpost(no_ROCs).timeWindow=winNo;
                                    ROCbandwidthpost(no_ROCs)=bwii;
                                    auROCpost(no_ROCs)=ROCoutpost(no_ROCs).roc.AUC-0.5;
                                    
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
                                        %Hit, split
                                        if p_val(no_dBs,bwii)<=0.05
                                            delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                            delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                            figure(bwii)
                                            hold on
                                            plot([0 1],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        else
                                            delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                            delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                            figure(bwii)
                                            hold on
                                            plot([3 4],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            
                                        end
                                        
                                        
                                        %CR, split
                                        if p_val(no_dBs+1,bwii)<=0.05
                                            delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                            delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                            figure(bwii+4)
                                            hold on
                                            plot([0 1],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        else
                                            delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                            delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                            figure(bwii+4)
                                            hold on
                                            plot([3 4],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        end
                                        
                                    else
                                        if groupNopre(no_dBs)==3
                                            %Hit, split
                                            if p_val(no_dBs,bwii)<=0.05
                                                delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                                delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                                figure(bwii)
                                                hold on
                                                plot([6 7],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            else
                                                delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                                delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                                figure(bwii)
                                                hold on
                                                plot([9 10],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            end
                                            
                                            %CR, split
                                            if p_val(no_dBs+1,bwii)<=0.05
                                                delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                                delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                                figure(bwii+4)
                                                hold on
                                                plot([6 7],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            else
                                                delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                                delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                                figure(bwii+4)
                                                hold on
                                                plot([9 10],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            end
                                            
                                        end
                                    end
                                end
                                
                                no_dBs=no_dBs+2;
                                
                            else
                                
                                if (sum(trials_in_event_pre)<comp_window)
                                    fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_pre),comp_window,file_pairs(fps,1),elec);
                                end
                                
                                if (sum(trials_in_event_post)<comp_window)
                                    fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_post),comp_window,file_pairs(fps,2),elec);
                                end
                            end
                            
                        else
                            
                            if (sum(trials_in_event_preHit)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_preHit),comp_window,file_pairs(fps,1),elec);
                            end
                            
                            if (sum(trials_in_event_preCR)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_preCR),comp_window,file_pairs(fps,1),elec);
                            end
                            if (sum(trials_in_event_postHit)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_postHit),comp_window,file_pairs(fps,2),elec);
                            end
                            if (sum(trials_in_event_postCR)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_postCR),comp_window,file_pairs(fps,2),elec);
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
        
        %Now plot the bar graphs and do anovan for LFP power
        p_vals_anovan=[];
        pvals_ancova=[];
        pvals_auROCancova=[];
        for bwii=1:4
            %Hits for NL->L
            figure(bwii)
            
            %Significant changes - Hit, NL->L
            sig_mean_pre=mean(delta_dB_powerpreHit((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
            if sum((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpreHit((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(0,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            plot([0 0],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            
            sig_mean_post=mean(delta_dB_powerpostHit((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
            if sum((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))>1
                pd=fitdist(delta_dB_powerpostHit((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(1,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            plot([1 1],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            
            plot([0 1],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
            
            %Changes not singificant - Hit, NL->L
            sig_mean_pre=mean(delta_dB_powerpreHit((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
            if sum((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpreHit((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(3,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            plot([3 3],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',3,'Color',[0.6 0 0])
            
            sig_mean_post=mean(delta_dB_powerpostHit((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
            if sum((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))>1
                pd=fitdist(delta_dB_powerpostHit((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(4,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            plot([4 4],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            
            plot([3 4],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
            
            %Significant changes - Hit, NLc->Lc
            sig_mean_pre=mean(delta_dB_powerpreHit((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            if sum((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpreHit((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(6,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
            plot([6 6],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            
            sig_mean_post=mean(delta_dB_powerpostHit((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            if sum((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpostHit((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(7,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
            plot([7 7],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            
            plot([6 7],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
            
            %Changes not singificant - Hit, NLc->Lc
            sig_mean_pre=mean(delta_dB_powerpreHit((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            if sum((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpreHit((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            
            plot(9,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
            plot([9 9],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',3,'Color',[0 0 0.6])
            
            sig_mean_post=mean(delta_dB_powerpostHit((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            if sum((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpostHit((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(10,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
            plot([10 10],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            
            plot([9 10],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
            
            title(['Change in ' freq_names{bwii} ' power dB elicited by the laser for Hits'])
            xticks(0:10)
            xticklabels({'','light','','','light','','','light','','','light',''})
            xlim([-0.5 10.5])
            ylm=ylim;
            ylim([ylm(1)-(ylm(2)-ylm(1))/10 ylm(2)+(ylm(2)-ylm(1))/10])
            txt1 = 'DBH-Cre x ROSA26 Halo';
            text(1,ylm(1),txt1)
            txt1 = 'DBH-Cre';
            text(7.5,ylm(1),txt1)
            
            
            %CRs for NL->L
            figure(bwii+4)
            
            %Significant changes - CR, NL->L
            sig_mean_pre=mean(delta_dB_powerpreCR((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
            
            if sum((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpreCR((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(0,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
            plot([0 0],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            
            sig_mean_post=mean(delta_dB_powerpostCR((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
            
            if sum((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))>1
                pd=fitdist(delta_dB_powerpostCR((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(1,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
            plot([1 1],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            
            plot([0 1],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
            
            %Changes not singificant - CR, NL->L
            sig_mean_pre=mean(delta_dB_powerpreCR((dB_power_changeCR==0)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
            
            if sum((dB_power_changeCR==1)&(ROCgroupNopre==0)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpreCR((dB_power_changeCR==0)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(3,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
            plot([3 3],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',3,'Color',[0.6 0 0])
            
            sig_mean_post=mean(delta_dB_powerpostCR((dB_power_changeCR==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
            
            if sum((dB_power_changeCR==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))
                pd=fitdist(delta_dB_powerpostCR((dB_power_changeCR==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(4,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
            plot([4 4],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            
            plot([3 4],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
            
            %Significant changes - CR, NLc->Lc
            sig_mean_pre=mean(delta_dB_powerpreCR((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            
            if sum((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpreCR((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(6,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            plot([6 6],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            
            sig_mean_post=mean(delta_dB_powerpostCR((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpost==bwii)));
            
            if sum((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpost==bwii))>1
                pd=fitdist(delta_dB_powerpostCR((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpost==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(7,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            plot([7 7],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            
            plot([6 7],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
            
            %Changes not singificant - CR, NLc->Lc
            sig_mean_pre=mean(delta_dB_powerpreCR((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            if sum((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpreCR((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(9,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            plot([9 9],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',3,'Color',[0 0 0.6])
            
            sig_mean_post=mean(delta_dB_powerpostCR((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            if sum((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpostCR((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(10,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            plot([10 10],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            
            plot([9 10],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
            
            title(['Change in ' freq_names{bwii} ' power dB elicited by the laser for CRs'])
            xticks(0:10)
            xticklabels({'','light','','','light','','','light','','','light',''})
            xlim([-0.5 10.5])
            ylm=ylim;
            ylim([ylm(1)-(ylm(2)-ylm(1))/10 ylm(2)+(ylm(2)-ylm(1))/10])
            txt1 = 'DBH-Cre x ROSA26 Halo';
            text(1,ylm(1),txt1)
            txt1 = 'DBH-Cre';
            text(7.5,ylm(1),txt1)
            
            
            
            
        end
        
        
        
        p_chi=[];
        for evTN1=1:length(eventType)
            fprintf(1, ['Significant changes in pairwise LFP power analysis for event: ' evTypeLabels{evTN1} '\n\n'])
            for bwii=1:4
                for grs=grpre
                    num_sig(grs)=sum(p_val((events==evTN1)&(groupNopre==grs),bwii)<=0.05);
                    tot_num(grs)=sum((events==evTN1)&(grs==groupNopre));
                    fprintf(1, ['Number significant for ' freq_names{bwii} ' and ' handles_drgb.drgbchoices.group_no_names{grs} ' = %d of %d\n'],num_sig(grs),tot_num(grs));
                end
                [p, Q]= chi2test([num_sig(grpre(1)), tot_num(grpre(1))-num_sig(grpre(1)); num_sig(grpre(2)), tot_num(grpre(2))-num_sig(grpre(2))]);
                fprintf(1, ['Chi squared p value  = %d\n\n'],p);
                p_chi=[p_chi p];
            end
            fprintf(1, '\n\n\n')
        end
        
        pFDRchi=drsFDRpval(p_chi);
        fprintf(1, ['pFDR for Chi squared p value  = %d\n\n'],pFDRchi);
        
        
    case 2
        %Compare auROC in the last few trials of pre with first few trials of post
        no_dBs=1;
        delta_dB_power_pre=[];
        no_ROCs=0;
        ROCoutpre=[];
        ROCoutpost=[];
        p_vals_ROC=[];
        delta_dB_powerpreHit=[];
        no_hits=0;
        
        fprintf(1, ['Pairwise auROC analysis for Hit and CR LFP power\n\n'])
        p_vals=[];
        for fps=1:no_file_pairs
            for elec=1:16
                
                lfpodNopre_ref=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                lfpodNopost_ref=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                
                if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref)))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref)))
                    
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower))
                        
                        
                        trials_in_event_preHit=(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).which_eventLFPPower(evHit,:)==1);
                        trials_in_event_preCR=(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).which_eventLFPPower(evCR,:)==1);
                        trials_in_event_pre=logical(trials_in_event_preHit+trials_in_event_preCR);
                        
                        trials_in_event_postHit=(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).which_eventLFPPower(evHit,:)==1);
                        trials_in_event_postCR=(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).which_eventLFPPower(evCR,:)==1);
                        trials_in_event_post=logical(trials_in_event_postHit+trials_in_event_postCR);
                        
                        if (sum(trials_in_event_preHit)>=comp_window)&(sum(trials_in_event_preCR)>=comp_window)&(sum(trials_in_event_postHit)>=comp_window)&(sum(trials_in_event_postCR)>=comp_window)
                            if (sum(trials_in_event_pre)>=comp_window_auROC)&(sum(trials_in_event_post)>=comp_window_auROC)
                                lfpodNopre=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                lfpodNopost=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                
                                %pre Hits
                                
                                
                                no_trials=0;
                                jj=length(trials_in_event_preHit);
                                while no_trials<comp_window
                                    if trials_in_event_preHit(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    jj=jj-1;
                                end
                                jj=jj+1;
                                these_trials=logical([zeros(1,jj-1) ones(1,length(trials_in_event_preHit)-jj+1)].*trials_in_event_preHit);
                                
                                this_dB_powerprerefHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerprerefHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower(these_trials,:));
                                
                                this_dB_powerpreHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpreHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre).allPower(these_trials,:));
                                
                                
                                %pre CRs
                                no_trials=0;
                                jj=length(trials_in_event_preCR);
                                while no_trials<comp_window
                                    if trials_in_event_preCR(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    jj=jj-1;
                                end
                                jj=jj+1;
                                these_trials=logical([zeros(1,jj-1) ones(1,length(trials_in_event_preCR)-jj+1)].*trials_in_event_preCR);
                                
                                this_dB_powerprerefCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerprerefCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower(these_trials,:));
                                
                                this_dB_powerpreCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpreCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre).allPower(these_trials,:));
                                
                                
                                %post Hits
                                no_trials=0;
                                jj=0;
                                while no_trials<comp_window
                                    jj=jj+1;
                                    if trials_in_event_postHit(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    
                                end
                                
                                these_trials=logical([ones(1,jj) zeros(1,length(trials_in_event_postHit)-jj)].*trials_in_event_postHit);
                                
                                this_dB_powerpostrefHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostrefHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower(these_trials,:));
                                
                                this_dB_powerpostHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost).allPower(these_trials,:));
                                
                                
                                %post CRs
                                no_trials=0;
                                jj=0;
                                while no_trials<comp_window
                                    jj=jj+1;
                                    if trials_in_event_postCR(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    
                                end
                                
                                these_trials=logical([ones(1,jj) zeros(1,length(trials_in_event_postCR)-jj)].*trials_in_event_postCR);
                                
                                this_dB_powerpostrefCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostrefCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower(these_trials,:));
                                
                                this_dB_powerpostCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost).allPower(these_trials,:));
                                
                                for bwii=1:no_bandwidths
                                    
                                    no_ROCs=no_ROCs+1;
                                    this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                    
                                    %Enter the pre Hits
                                    this_delta_dB_powerpreHit=zeros(comp_window,1);
                                    this_delta_dB_powerpreHit=mean(this_dB_powerpreHit(:,this_band)-this_dB_powerprerefHit(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:comp_window,1)=this_delta_dB_powerpreHit;
                                    roc_data(1:comp_window,2)=zeros(comp_window,1);
                                    
                                    %Enter pre CR
                                    this_delta_dB_powerpreCR=zeros(comp_window,1);
                                    this_delta_dB_powerpreCR=mean(this_dB_powerpreCR(:,this_band)-this_dB_powerprerefCR(:,this_band),2);
                                    roc_data(comp_window+1:2*comp_window,1)=this_delta_dB_powerpreCR;
                                    roc_data(comp_window+1:2*comp_window,2)=ones(comp_window,1);
                                    
                                    
                                    %Find pre ROC
                                    ROCoutpre(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCoutpre(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopre_ref).fileNo;
                                    ROCgroupNopre(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).fileNo);
                                    ROCoutpre(no_ROCs).timeWindow=winNo;
                                    ROCbandwidthpre(no_ROCs)=bwii;
                                    auROCpre(no_ROCs)=ROCoutpre(no_ROCs).roc.AUC-0.5;
                                    
                                    p_vals_ROC=[p_vals_ROC ROCoutpre(no_ROCs).roc.p];
                                    
                                    %Enter the post Hits
                                    this_delta_dB_powerpostHit=zeros(comp_window,1);
                                    this_delta_dB_powerpostHit=mean(this_dB_powerpostHit(:,this_band)-this_dB_powerpostrefHit(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:comp_window,1)=this_delta_dB_powerpostHit;
                                    roc_data(1:comp_window,2)=zeros(comp_window,1);
                                    
                                    %Enter post CR
                                    this_delta_dB_powerpostCR=zeros(comp_window,1);
                                    this_delta_dB_powerpostCR=mean(this_dB_powerpostCR(:,this_band)-this_dB_powerpostrefCR(:,this_band),2);
                                    roc_data(comp_window+1:2*comp_window,1)=this_delta_dB_powerpostCR;
                                    roc_data(comp_window+1:2*comp_window,2)=ones(comp_window,1);
                                    
                                    
                                    %Find post ROC
                                    ROCoutpost(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCoutpost(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopost_ref).fileNo;
                                    ROCgroupNopost(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).fileNo);
                                    ROCoutpost(no_ROCs).timeWindow=winNo;
                                    ROCbandwidthpost(no_ROCs)=bwii;
                                    auROCpost(no_ROCs)=ROCoutpost(no_ROCs).roc.AUC-0.5;
                                    
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
                                        %                                         %Hit, split
                                        %                                         if p_val(no_dBs,bwii)<=0.05
                                        %                                             delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                        %                                             delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                        %                                             figure(bwii)
                                        %                                             hold on
                                        %                                             plot([0 1],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        %                                         else
                                        %                                             delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                        %                                             delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                        %                                             figure(bwii)
                                        %                                             hold on
                                        %                                             plot([3 4],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        %
                                        %                                         end
                                        
                                        %Hit, all points
                                        delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                        delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                        figure(bwii+12)
                                        hold on
                                        plot([0 1],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        
                                        %                                         %CR, split
                                        %                                         if p_val(no_dBs+1,bwii)<=0.05
                                        %                                             delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                        %                                             delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                        %                                             figure(bwii+4)
                                        %                                             hold on
                                        %                                             plot([0 1],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        %                                         else
                                        %                                             delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                        %                                             delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                        %                                             figure(bwii+4)
                                        %                                             hold on
                                        %                                             plot([3 4],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        %                                         end
                                        
                                        %CR, all points
                                        delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                        delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                        figure(bwii+4+12)
                                        hold on
                                        plot([0 1],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                    else
                                        if groupNopre(no_dBs)==3
                                            %Hit, split
                                            %                                             if p_val(no_dBs,bwii)<=0.05
                                            %                                                 delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                            %                                                 delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                            %                                                 figure(bwii)
                                            %                                                 hold on
                                            %                                                 plot([6 7],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            %                                             else
                                            %                                                 delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                            %                                                 delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                            %                                                 figure(bwii)
                                            %                                                 hold on
                                            %                                                 plot([9 10],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            %                                             end
                                            
                                            %Hit, all points
                                            delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                            delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                            figure(bwii+12)
                                            hold on
                                            plot([3 4],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            
                                            %                                             %CR, split
                                            %                                             if p_val(no_dBs+1,bwii)<=0.05
                                            %                                                 delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                            %                                                 delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                            %                                                 figure(bwii+4)
                                            %                                                 hold on
                                            %                                                 plot([6 7],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            %                                             else
                                            %                                                 delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                            %                                                 delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                            %                                                 figure(bwii+4)
                                            %                                                 hold on
                                            %                                                 plot([9 10],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            %                                             end
                                            
                                            %CR, all points
                                            delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                            delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                            figure(bwii+4+12)
                                            hold on
                                            plot([3 4],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        end
                                    end
                                end
                                
                                no_dBs=no_dBs+2;
                                
                            else
                                
                                if (sum(trials_in_event_pre)<comp_window)
                                    fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_pre),comp_window,file_pairs(fps,1),elec);
                                end
                                
                                if (sum(trials_in_event_post)<comp_window)
                                    fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_post),comp_window,file_pairs(fps,2),elec);
                                end
                            end
                            
                        else
                            
                            if (sum(trials_in_event_preHit)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_preHit),comp_window,file_pairs(fps,1),elec);
                            end
                            
                            if (sum(trials_in_event_preCR)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_preCR),comp_window,file_pairs(fps,1),elec);
                            end
                            if (sum(trials_in_event_postHit)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_postHit),comp_window,file_pairs(fps,2),elec);
                            end
                            if (sum(trials_in_event_postCR)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_postCR),comp_window,file_pairs(fps,2),elec);
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
        
        %Now plot the bar graphs and do anovan for LFP power
        p_vals_anovan=[];
        pvals_ancova=[];
        pvals_auROCancova=[];
        for bwii=1:4
            %Hits for NL->L
            %             figure(bwii)
            %
            %             %Significant changes - Hit, NL->L
            %             sig_mean_pre=mean(delta_dB_powerpreHit((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreHit((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(0,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            %             plot([0 0],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             sig_mean_post=mean(delta_dB_powerpostHit((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
            %             if sum((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))>1
            %                 pd=fitdist(delta_dB_powerpostHit((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(1,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            %             plot([1 1],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             plot([0 1],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             %Changes not singificant - Hit, NL->L
            %             sig_mean_pre=mean(delta_dB_powerpreHit((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreHit((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(3,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            %             plot([3 3],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',3,'Color',[0.6 0 0])
            %
            %             sig_mean_post=mean(delta_dB_powerpostHit((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
            %             if sum((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))>1
            %                 pd=fitdist(delta_dB_powerpostHit((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(4,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            %             plot([4 4],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             plot([3 4],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             %Significant changes - Hit, NLc->Lc
            %             sig_mean_pre=mean(delta_dB_powerpreHit((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreHit((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(6,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
            %             plot([6 6],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             sig_mean_post=mean(delta_dB_powerpostHit((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpostHit((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(7,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
            %             plot([7 7],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             plot([6 7],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             %Changes not singificant - Hit, NLc->Lc
            %             sig_mean_pre=mean(delta_dB_powerpreHit((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreHit((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %
            %             plot(9,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
            %             plot([9 9],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',3,'Color',[0 0 0.6])
            %
            %             sig_mean_post=mean(delta_dB_powerpostHit((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpostHit((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(10,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
            %             plot([10 10],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             plot([9 10],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             title(['Change in ' freq_names{bwii} ' power dB elicited by the laser for Hits'])
            %             xticks(0:10)
            %             xticklabels({'','light','','','light','','','light','','','light',''})
            %             xlim([-0.5 10.5])
            %             ylm=ylim;
            %             ylim([ylm(1)-(ylm(2)-ylm(1))/10 ylm(2)+(ylm(2)-ylm(1))/10])
            %             txt1 = 'DBH-Cre x ROSA26 Halo';
            %             text(1,ylm(1),txt1)
            %             txt1 = 'DBH-Cre';
            %             text(7.5,ylm(1),txt1)
            %
            %
            %             %CRs for NL->L
            %             figure(bwii+4)
            %
            %             %Significant changes - CR, NL->L
            %             sig_mean_pre=mean(delta_dB_powerpreCR((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
            %
            %             if sum((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreCR((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(0,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
            %             plot([0 0],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             sig_mean_post=mean(delta_dB_powerpostCR((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
            %
            %             if sum((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))>1
            %                 pd=fitdist(delta_dB_powerpostCR((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(1,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
            %             plot([1 1],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             plot([0 1],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             %Changes not singificant - CR, NL->L
            %             sig_mean_pre=mean(delta_dB_powerpreCR((dB_power_changeCR==0)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
            %
            %             if sum((dB_power_changeCR==1)&(ROCgroupNopre==0)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreCR((dB_power_changeCR==0)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(3,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
            %             plot([3 3],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',3,'Color',[0.6 0 0])
            %
            %             sig_mean_post=mean(delta_dB_powerpostCR((dB_power_changeCR==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
            %
            %             if sum((dB_power_changeCR==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))
            %                 pd=fitdist(delta_dB_powerpostCR((dB_power_changeCR==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(4,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
            %             plot([4 4],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             plot([3 4],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             %Significant changes - CR, NLc->Lc
            %             sig_mean_pre=mean(delta_dB_powerpreCR((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            %
            %             if sum((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreCR((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(6,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            %             plot([6 6],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             sig_mean_post=mean(delta_dB_powerpostCR((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpost==bwii)));
            %
            %             if sum((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpost==bwii))>1
            %                 pd=fitdist(delta_dB_powerpostCR((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpost==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(7,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            %             plot([7 7],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             plot([6 7],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             %Changes not singificant - CR, NLc->Lc
            %             sig_mean_pre=mean(delta_dB_powerpreCR((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreCR((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(9,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            %             plot([9 9],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',3,'Color',[0 0 0.6])
            %
            %             sig_mean_post=mean(delta_dB_powerpostCR((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpostCR((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(10,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            %             plot([10 10],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             plot([9 10],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             title(['Change in ' freq_names{bwii} ' power dB elicited by the laser for CRs'])
            %             xticks(0:10)
            %             xticklabels({'','light','','','light','','','light','','','light',''})
            %             xlim([-0.5 10.5])
            %             ylm=ylim;
            %             ylim([ylm(1)-(ylm(2)-ylm(1))/10 ylm(2)+(ylm(2)-ylm(1))/10])
            %             txt1 = 'DBH-Cre x ROSA26 Halo';
            %             text(1,ylm(1),txt1)
            %             txt1 = 'DBH-Cre';
            %             text(7.5,ylm(1),txt1)
            %
            %             pffft=1;
            %ANOVAn
            %Note: This is a first try, we should do a repeated measures...
            
            %For Hits
            all_delta_dB_power=[];
            laser=[];
            halo=[];
            
            %No laser, pre, with halo
            all_delta_dB_power=[all_delta_dB_power delta_dB_powerpreHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))];
            laser=[laser zeros(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
            halo=[halo ones(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
            
            %Laser, post, with halo
            all_delta_dB_power=[all_delta_dB_power delta_dB_powerpostHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))];
            laser=[laser ones(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
            halo=[halo ones(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
            
            %No laser, pre, no halo
            all_delta_dB_power=[all_delta_dB_power delta_dB_powerpreHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))];
            laser=[laser zeros(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
            halo=[halo zeros(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
            
            %Laser, post, with halo
            all_delta_dB_power=[all_delta_dB_power delta_dB_powerpostHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))];
            laser=[laser ones(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
            halo=[halo zeros(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
            
            p=anovan(all_delta_dB_power,{laser,halo},'model','interaction','varnames',{'laser','halo'},'display','off');
            p_vals_anovan=[p_vals_anovan p(3)];
            fprintf(1, ['anovan p value ' freq_names{bwii} ' and ' evTypeLabels{1} ' = %d\n'],p(3));
            
            %For CRs
            all_delta_dB_power=[];
            laser=[];
            halo=[];
            
            %No laser, pre, with halo
            all_delta_dB_power=[all_delta_dB_power delta_dB_powerpreCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))];
            laser=[laser zeros(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
            halo=[halo ones(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
            
            %Laser, post, with halo
            all_delta_dB_power=[all_delta_dB_power delta_dB_powerpostCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))];
            laser=[laser ones(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
            halo=[halo ones(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
            
            %No laser, pre, no halo
            all_delta_dB_power=[all_delta_dB_power delta_dB_powerpreCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))];
            laser=[laser zeros(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
            halo=[halo zeros(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
            
            %Laser, post, with halo
            all_delta_dB_power=[all_delta_dB_power delta_dB_powerpostCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))];
            laser=[laser ones(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
            halo=[halo zeros(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
            
            p=anovan(all_delta_dB_power,{laser,halo},'model','interaction','varnames',{'laser','halo'},'display','off');
            p_vals_anovan=[p_vals_anovan p(3)];
            fprintf(1, ['anovan p value ' freq_names{bwii} ' and ' evTypeLabels{2} ' = %d\n\n'],p(3));
            
            
            %Now plot the LFP power, all points
            %Hits for NL->L
            figure(bwii+12)
            
            %All points - Hit, NL->L
            mean_pre=mean(delta_dB_powerpreHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
            if sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpreHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(0,mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            plot([0 0],[mean_pre-ci95 mean_pre+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            
            mean_post=mean(delta_dB_powerpostHit((ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
            if sum((ROCgroupNopre==1)&(ROCbandwidthpost==bwii))>1
                pd=fitdist(delta_dB_powerpostHit((ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(1,mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            plot([1 1],[mean_post-ci95 mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            
            plot([0 1],[mean_pre mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
            
            
            
            %All points - Hit, NLc->Lc
            mean_pre=mean(delta_dB_powerpreHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            if sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpreHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(3,mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
            plot([3 3],[mean_pre-ci95 mean_pre+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            
            mean_post=mean(delta_dB_powerpostHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            if sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpostHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            ci95=pd.mu-ci_tbl(1,1);
            plot(4,mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
            plot([4 4],[mean_post-ci95 mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            
            plot([3 4],[mean_pre mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
            
            
            title(['Change in ' freq_names{bwii} ' power dB elicited by the laser for ' evTypeLabels{1}])
            xticks(0:10)
            xticklabels({'','light','','','light',''})
            xlim([-0.5 4.5])
            ylm=ylim;
            ylim([ylm(1)-(ylm(2)-ylm(1))/10 ylm(2)+(ylm(2)-ylm(1))/10])
            txt1 = 'DBH-Cre x ROSA26 Halo';
            text(0,ylm(1),txt1)
            txt1 = 'DBH-Cre';
            text(3,ylm(1),txt1)
            
            
            %CRs for NL->L
            figure(bwii+4+12)
            
            %All points - CR, NL->L
            mean_pre=mean(delta_dB_powerpreCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
            
            if sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpreCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(0,mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
            plot([0 0],[mean_pre-ci95 mean_pre+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            
            mean_post=mean(delta_dB_powerpostCR((ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
            
            if sum((ROCgroupNopre==1)&(ROCbandwidthpost==bwii))>1
                pd=fitdist(delta_dB_powerpostCR((ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(1,mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
            plot([1 1],[mean_post-ci95 mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            
            plot([0 1],[mean_pre mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
            
            
            
            %All points - CR, NLc->Lc
            mean_pre=mean(delta_dB_powerpreCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            
            if sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpreCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(3,mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            plot([3 3],[mean_pre-ci95 mean_pre+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            
            mean_post=mean(delta_dB_powerpostCR((ROCgroupNopre==3)&(ROCbandwidthpost==bwii)));
            
            if sum((ROCgroupNopre==3)&(ROCbandwidthpost==bwii))>1
                pd=fitdist(delta_dB_powerpostCR((ROCgroupNopre==3)&(ROCbandwidthpost==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot(4,mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            plot([4 4],[mean_post-ci95 mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            
            plot([3 4],[mean_pre mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
            
            
            
            title(['Change in ' freq_names{bwii} ' power dB elicited by the laser for ' evTypeLabels{2}])
            xticks(0:10)
            xticklabels({'','light','','','light',''})
            xlim([-0.5 4.5])
            ylm=ylim;
            ylim([ylm(1)-(ylm(2)-ylm(1))/10 ylm(2)+(ylm(2)-ylm(1))/10])
            txt1 = 'DBH-Cre x ROSA26 Halo';
            text(0,ylm(1),txt1)
            txt1 = 'DBH-Cre';
            text(3,ylm(1),txt1)
            
            
            %Analysis of covariance
            
            %For Hit
            this_delta_dB_pre_powerHit=[];
            this_delta_dB_pre_powerHit=delta_dB_powerpreHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))';
            this_delta_dB_pre_powerHit=[this_delta_dB_pre_powerHit; delta_dB_powerpreHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))'];
            
            
            this_delta_dB_post_powerHit=[];
            this_delta_dB_post_powerHit=delta_dB_powerpostHit((ROCgroupNopre==1)&(ROCbandwidthpost==bwii))';
            this_delta_dB_post_powerHit=[this_delta_dB_post_powerHit; delta_dB_powerpostHit((ROCgroupNopre==3)&(ROCbandwidthpost==bwii))'];
            
            pre_post=[];
            pre_post=[zeros(sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)),1); ones(sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)),1)];
            
            [h,atab,ctab,stats] = aoctool(this_delta_dB_pre_powerHit,this_delta_dB_post_powerHit,pre_post,0.05,'','','','off');
            
            
            pvals_ancova=[pvals_ancova atab{4,6}];
            fprintf(1, ['ancova delta dB p value ' freq_names{bwii} ' and ' evTypeLabels{1} ' = %d\n'],atab{4,6});
            
            %Do ancova figures for Hits
            figure(20+bwii)
            h1=plot(delta_dB_powerpreHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)),delta_dB_powerpostHit((ROCgroupNopre==1)&(ROCbandwidthpost==bwii)),'or');
            hold on
            h2=plot(delta_dB_powerpreHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)),delta_dB_powerpostHit((ROCgroupNopre==3)&(ROCbandwidthpost==bwii)),'ob');
            
            slope_pre=ctab{5,2}+ctab{6,2};
            int_pre=ctab{2,2}+ctab{3,2};
            min_x=min([min(delta_dB_powerpreHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) min(delta_dB_powerpreHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
            max_x=max([max(delta_dB_powerpreHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) max(delta_dB_powerpreHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
            x=[min_x max_x];
            plot(x,slope_pre*x+int_pre,'-r')
            
            slope_post=ctab{5,2}+ctab{7,2};
            int_post=ctab{2,2}+ctab{4,2};
            x=[min_x max_x];
            plot(x,slope_post*x+int_post,'-b')
            
            title(['post vs pre power dB for ' freq_names{bwii} ' and ' evTypeLabels{1}])
            xlabel('pre dB')
            ylabel('post dB')
            legend([h1 h2],'halo','no halo')
            
            %For CR
            this_delta_dB_pre_powerCR=[];
            this_delta_dB_pre_powerCR=delta_dB_powerpreCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))';
            this_delta_dB_pre_powerCR=[this_delta_dB_pre_powerCR; delta_dB_powerpreCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))'];
            
            
            this_delta_dB_post_powerCR=[];
            this_delta_dB_post_powerCR=delta_dB_powerpostCR((ROCgroupNopre==1)&(ROCbandwidthpost==bwii))';
            this_delta_dB_post_powerCR=[this_delta_dB_post_powerCR; delta_dB_powerpostCR((ROCgroupNopre==3)&(ROCbandwidthpost==bwii))'];
            
            pre_post=[];
            pre_post=[zeros(sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)),1); ones(sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)),1)];
            
            [h,atab,ctab,stats] = aoctool(this_delta_dB_pre_powerCR,this_delta_dB_post_powerCR,pre_post,0.05,'','','','off');
            
            
            pvals_ancova=[pvals_ancova atab{4,6}];
            fprintf(1, ['ancova delta dB p value ' freq_names{bwii} ' and ' evTypeLabels{2} ' = %d\n\n'],atab{4,6});
            
            %Do ancova figures for CR
            figure(24+bwii)
            h1=plot(delta_dB_powerpreCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)),delta_dB_powerpostCR((ROCgroupNopre==1)&(ROCbandwidthpost==bwii)),'or');
            hold on
            h2=plot(delta_dB_powerpreCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)),delta_dB_powerpostCR((ROCgroupNopre==3)&(ROCbandwidthpost==bwii)),'ob');
            
            slope_pre=ctab{5,2}+ctab{6,2};
            int_pre=ctab{2,2}+ctab{3,2};
            min_x=min([min(delta_dB_powerpreCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) min(delta_dB_powerpreCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
            max_x=max([max(delta_dB_powerpreCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) max(delta_dB_powerpreCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
            x=[min_x max_x];
            plot(x,slope_pre*x+int_pre,'-r')
            
            slope_post=ctab{5,2}+ctab{7,2};
            int_post=ctab{2,2}+ctab{4,2};
            x=[min_x max_x];
            plot(x,slope_post*x+int_post,'-b')
            
            title(['post vs pre power dB for ' freq_names{bwii} ' and ' evTypeLabels{2}])
            xlabel('pre dB')
            ylabel('post dB')
            legend([h1 h2],'halo','no halo')
            
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
            figure(28+bwii)
            h1=plot(auROCpre((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)),auROCpost((ROCgroupNopre==1)&(ROCbandwidthpost==bwii)),'or');
            hold on
            h2=plot(auROCpre((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)),auROCpost((ROCgroupNopre==3)&(ROCbandwidthpost==bwii)),'ob');
            
            slope_pre=ctab{5,2}+ctab{6,2};
            int_pre=ctab{2,2}+ctab{3,2};
            min_x=min([min(auROCpre((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) min(auROCpre((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
            max_x=max([max(auROCpre((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) max(auROCpre((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
            x=[min_x max_x];
            plot(x,slope_pre*x+int_pre,'-r')
            
            slope_post=ctab{5,2}+ctab{7,2};
            int_post=ctab{2,2}+ctab{4,2};
            x=[min_x max_x];
            plot(x,slope_post*x+int_post,'-b')
            
            title(['post vs pre auROC for ' freq_names{bwii} ])
            xlabel('pre auROC')
            ylabel('post auROC')
            legend([h1 h2],'halo','no halo')
        end
        
        pFDRanovan=drsFDRpval(p_vals_anovan);
        fprintf(1, ['pFDR for anovan p value  = %d\n\n'],pFDRanovan);
        
        pFDRancova=drsFDRpval(pvals_ancova);
        fprintf(1, ['pFDR for power dB ancova p value  = %d\n\n'], pFDRancova);
        
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
        figNo=8;
        p_val_ROC=[];
        
        for bwii=1:4
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            hold on
            n_cum=0;
            this_legend=[];
            for grs=1:2
                %pre
                [f_auROC,x_auROC] = drg_ecdf(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)));
                n_cum=n_cum+1;
                plot(x_auROC,f_auROC,these_lines{n_cum})
                this_legend=[this_legend '''' handles_drgb.drgbchoices.group_no_names{grpre(grs)}  '''' ','];
                
                %post
                [f_auROC,x_auROC] = drg_ecdf(auROCpost((ROCgroupNopost==grpost(grs))&(ROCbandwidthpost==bwii)));
                n_cum=n_cum+1;
                plot(x_auROC,f_auROC,these_lines{n_cum})
                this_legend=[this_legend '''' handles_drgb.drgbchoices.group_no_names{grpost(grs)}  '''' ','];
                
                pROC=ranksum(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)),auROCpost((ROCgroupNopost==grpost(grs))&(ROCbandwidthpost==bwii)));
                p_val_ROC=[ p_val_ROC pROC];
                fprintf(1, ['p value for auROC for ' freq_names{bwii} ' ' handles_drgb.drgbchoices.group_no_names{grs} ' vs. ' handles_drgb.drgbchoices.group_no_names{grs+1} ' = %d\n'],pROC);
            end
            this_legend=['legend(' this_legend(1:end-1) ')'];
            eval(this_legend)
            title(['Cumulative histogram for ' freq_names{bwii} ' auROC for LFPs'])
        end
        
        pFDRauROC=drsFDRpval(p_val_ROC);
        fprintf(1, ['pFDR for auROC  = %d\n\n'],pFDRauROC);
        
        %Now do anovan for auROC
        p_vals_anovan_auROC=[];
        for bwii=1:4
            
            %ANOVAn
            %Note: This is a first try, we should do a repeated measures...
            
            %For Hits
            all_auROCs=[];
            laser=[];
            halo=[];
            
            %No laser, pre, with halo
            all_auROCs=[all_auROCs auROCpre((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))];
            laser=[laser zeros(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
            halo=[halo ones(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
            
            %Laser, post, with halo
            all_auROCs=[all_auROCs auROCpost((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))];
            laser=[laser ones(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
            halo=[halo ones(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
            
            %No laser, pre, no halo
            all_auROCs=[all_auROCs auROCpre((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))];
            laser=[laser zeros(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
            halo=[halo zeros(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
            
            %Laser, post, with halo
            all_auROCs=[all_auROCs auROCpost((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))];
            laser=[laser ones(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
            halo=[halo zeros(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
            
            p=anovan(all_auROCs,{laser,halo},'model','interaction','varnames',{'laser','halo'},'display','off');
            p_vals_anovan_auROC=[p_vals_anovan_auROC p(3)];
            fprintf(1, ['anovan p value for auROC' freq_names{bwii} ' = %d\n'],p(3));
            
            
            
        end
        
        pFDRanovanauROC=drsFDRpval(p_vals_anovan_auROC);
        fprintf(1, ['pFDR for anovan for auROC p value  = %d\n\n'],pFDRanovanauROC);
        
    case 3
        %Generage figure 1 for Daniels' LFP power paper. For the proficient mice in the firstlast LFP batch analysis
        %Plot the spectrum for S+ vs S-, plot LFP power for S+ vs S- for each electrode and plot auROCs
        no_dBs=1;
        delta_dB_power_pre=[];
        no_ROCs=0;
        ROCoutpre=[];
        ROCoutpost=[];
        p_vals_ROC=[];
        delta_dB_powerpreHit=[];
        no_hits=0;
        noWB=0;
        delta_dB_powerpreHitWB=[];
        delta_dB_powerpreCRWB=[];
        
        fprintf(1, ['Pairwise auROC analysis for Fig 1 of Daniel''s paper\n\n'])
        p_vals=[];
        for fps=1:no_file_pairs
            for elec=1:16
                
                lfpodNopre_ref=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                lfpodNopost_ref=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                
                if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref)))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref)))
                    
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower))
                        
                        
                        trials_in_event_preHit=(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).which_eventLFPPower(evHit,:)==1);
                        trials_in_event_preCR=(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).which_eventLFPPower(evCR,:)==1);
                        trials_in_event_pre=logical(trials_in_event_preHit+trials_in_event_preCR);
                        
                        trials_in_event_postHit=(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).which_eventLFPPower(evHit,:)==1);
                        trials_in_event_postCR=(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).which_eventLFPPower(evCR,:)==1);
                        trials_in_event_post=logical(trials_in_event_postHit+trials_in_event_postCR);
                        
                        if (sum(trials_in_event_preHit)>=comp_window)&(sum(trials_in_event_preCR)>=comp_window)&(sum(trials_in_event_postHit)>=comp_window)&(sum(trials_in_event_postCR)>=comp_window)
                            if (sum(trials_in_event_pre)>=comp_window_auROC)&(sum(trials_in_event_post)>=comp_window_auROC)
                                lfpodNopre=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                lfpodNopost=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                
                                %pre Hits
                                
                                
                                no_trials=0;
                                jj=length(trials_in_event_preHit);
                                while no_trials<comp_window
                                    if trials_in_event_preHit(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    jj=jj-1;
                                end
                                jj=jj+1;
                                these_trials=logical([zeros(1,jj-1) ones(1,length(trials_in_event_preHit)-jj+1)].*trials_in_event_preHit);
                                
                                this_dB_powerprerefHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerprerefHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower(these_trials,:));
                                
                                this_dB_powerpreHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpreHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre).allPower(these_trials,:));
                                
                                
                                %pre CRs
                                no_trials=0;
                                jj=length(trials_in_event_preCR);
                                while no_trials<comp_window
                                    if trials_in_event_preCR(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    jj=jj-1;
                                end
                                jj=jj+1;
                                these_trials=logical([zeros(1,jj-1) ones(1,length(trials_in_event_preCR)-jj+1)].*trials_in_event_preCR);
                                
                                this_dB_powerprerefCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerprerefCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower(these_trials,:));
                                
                                this_dB_powerpreCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpreCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre).allPower(these_trials,:));
                                
                                
                                %post Hits
                                no_trials=0;
                                jj=0;
                                while no_trials<comp_window
                                    jj=jj+1;
                                    if trials_in_event_postHit(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    
                                end
                                
                                these_trials=logical([ones(1,jj) zeros(1,length(trials_in_event_postHit)-jj)].*trials_in_event_postHit);
                                
                                this_dB_powerpostrefHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostrefHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower(these_trials,:));
                                
                                this_dB_powerpostHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost).allPower(these_trials,:));
                                
                                
                                %post CRs
                                no_trials=0;
                                jj=0;
                                while no_trials<comp_window
                                    jj=jj+1;
                                    if trials_in_event_postCR(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    
                                end
                                
                                these_trials=logical([ones(1,jj) zeros(1,length(trials_in_event_postCR)-jj)].*trials_in_event_postCR);
                                
                                this_dB_powerpostrefCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostrefCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower(these_trials,:));
                                
                                this_dB_powerpostCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost).allPower(these_trials,:));
                                
                                %Wide band spectrum
                                noWB=noWB+1;
                                WBgroupNo(noWB)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                delta_dB_powerpreHitWB(noWB,:)=mean(this_dB_powerpreHit-this_dB_powerprerefHit,1);
                                delta_dB_powerpreCRWB(noWB,:)=mean(this_dB_powerpreCR-this_dB_powerprerefCR,1);
                                
                                
                                %Do per badwidth analysis
                                for bwii=1:no_bandwidths
                                    if (handles_drgb.drgbchoices.group_no(file_pairs(fps,1))==1)||(handles_drgb.drgbchoices.group_no(file_pairs(fps,1))==3)
                                        no_ROCs=no_ROCs+1;
                                        this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                        
                                        %Enter the pre Hits
                                        this_delta_dB_powerpreHit=zeros(comp_window,1);
                                        this_delta_dB_powerpreHit=mean(this_dB_powerpreHit(:,this_band)-this_dB_powerprerefHit(:,this_band),2);
                                        roc_data=[];
                                        roc_data(1:comp_window,1)=this_delta_dB_powerpreHit;
                                        roc_data(1:comp_window,2)=zeros(comp_window,1);
                                        
                                        %Enter pre CR
                                        this_delta_dB_powerpreCR=zeros(comp_window,1);
                                        this_delta_dB_powerpreCR=mean(this_dB_powerpreCR(:,this_band)-this_dB_powerprerefCR(:,this_band),2);
                                        roc_data(comp_window+1:2*comp_window,1)=this_delta_dB_powerpreCR;
                                        roc_data(comp_window+1:2*comp_window,2)=ones(comp_window,1);
                                        
                                        
                                        %Find pre ROC
                                        ROCoutpre(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                        ROCoutpre(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopre_ref).fileNo;
                                        ROCgroupNopre(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).fileNo);
                                        ROCoutpre(no_ROCs).timeWindow=winNo;
                                        ROCbandwidthpre(no_ROCs)=bwii;
                                        auROCpre(no_ROCs)=ROCoutpre(no_ROCs).roc.AUC-0.5;
                                        p_valROCpre(no_ROCs)=ROCoutpre(no_ROCs).roc.p;
                                        
                                        p_vals_ROC=[p_vals_ROC ROCoutpre(no_ROCs).roc.p];
                                        
                                        
                                        delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                        delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                        
                                        
                                        %Plot this point
                                        figure(bwii+1)
                                        pos2=[0.8 0.1 0.1 0.8];
                                        subplot('Position',pos2)
                                        hold on
                                        plot([0 1],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpreCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        
                                        pffft=1
                                        
                                    end
                                end
                                
                                
                                
                            else
                                
                                if (sum(trials_in_event_pre)<comp_window)
                                    fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_pre),comp_window,file_pairs(fps,1),elec);
                                end
                                
                                if (sum(trials_in_event_post)<comp_window)
                                    fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_post),comp_window,file_pairs(fps,2),elec);
                                end
                            end
                            
                        else
                            
                            if (sum(trials_in_event_preHit)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_preHit),comp_window,file_pairs(fps,1),elec);
                            end
                            
                            if (sum(trials_in_event_preCR)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_preCR),comp_window,file_pairs(fps,1),elec);
                            end
                            if (sum(trials_in_event_postHit)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_postHit),comp_window,file_pairs(fps,2),elec);
                            end
                            if (sum(trials_in_event_postCR)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_postCR),comp_window,file_pairs(fps,2),elec);
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
        
        
        %Now plot the bounded line for pre
        
        %Calculate the mean and 95% CI for Hit
        for ifreq=1:length(frequency)
            %             pd=fitdist(delta_dB_powerpreHitWB(WBgroupNo==1,ifreq),'Normal');
            %             ci=paramci(pd);
            %             dB_Hit_ci(ifreq)=pd.mu-ci(1,1);
            dB_Hit_mean(ifreq)=mean(delta_dB_powerpreHitWB(WBgroupNo==1,ifreq));
            CI = bootci(1000, @mean, delta_dB_powerpreHitWB(WBgroupNo==1,ifreq));
            dB_Hit_ci(ifreq)=dB_Hit_mean(ifreq)-CI(1);
        end
        
        figure(1)
        [hl1, hp1] = boundedline(frequency,dB_Hit_mean, dB_Hit_ci, 'r');
        
        %Calculate the mean and 95% CI for CR
        for ifreq=1:length(frequency)
            %             pd=fitdist(delta_dB_powerpreCRWB(WBgroupNo==1,ifreq),'Normal');
            %             ci=paramci(pd);
            %             dB_CR_ci(ifreq)=pd.mu-ci(1,1);
            dB_CR_mean(ifreq)=mean(delta_dB_powerpreCRWB(WBgroupNo==1,ifreq));
            CI = bootci(1000, @mean, delta_dB_powerpreCRWB(WBgroupNo==1,ifreq));
            dB_CR_ci(ifreq)=dB_CR_mean(ifreq)-CI(1);
        end
        
        hold on
        [hl2, hp2] = boundedline(frequency,dB_CR_mean, dB_CR_ci, 'b');
        xlabel('Frequency (Hz)')
        ylabel('delta Power (dB)')
        legend([hl1 hl2],'S+','S-')
        
        %Now plot the histograms and the average
        for bwii=1:4
            %Plot the average
            figure(bwii+1)
            pos2=[0.8 0.1 0.1 0.8];
            subplot('Position',pos2)
            hold on
            plot([0 1],[mean(delta_dB_powerpreHit(ROCbandwidthpre==bwii)) mean(delta_dB_powerpreCR(ROCbandwidthpre==bwii))],'-k','LineWidth', 3)
            CI = bootci(1000, @mean, delta_dB_powerpreHit(ROCbandwidthpre==bwii));
            plot([0 0],CI,'-r','LineWidth',3)
            plot(0,mean(delta_dB_powerpreHit(ROCbandwidthpre==bwii)),'or','MarkerSize', 10,'MarkerFace','r')
            CI = bootci(1000, @mean, delta_dB_powerpreCR(ROCbandwidthpre==bwii));
            plot([1 1],CI,'-b','LineWidth',3)
            plot(1,mean(delta_dB_powerpreCR(ROCbandwidthpre==bwii)),'ob','MarkerSize', 10,'MarkerFace','b')
            ylabel('delta Power (dB)')
            ylim([-10 15])
            
            %Plot the histograms
            
            maxdB=max([max(delta_dB_powerpreHit(ROCbandwidthpre==bwii)) max(delta_dB_powerpreCR(ROCbandwidthpre==bwii))]);
            mindB=min([min(delta_dB_powerpreHit(ROCbandwidthpre==bwii)) min(delta_dB_powerpreCR(ROCbandwidthpre==bwii))]);
            edges=[-15:0.5:15];
            pos2=[0.1 0.1 0.6 0.8];
            subplot('Position',pos2)
            hold on
            
            h1=histogram(delta_dB_powerpreCR(ROCbandwidthpre==bwii),edges);
            h1.FaceColor='b';
            h2=histogram(delta_dB_powerpreHit(ROCbandwidthpre==bwii),edges);
            h2.FaceColor='r';
            xlabel('delta Power (dB)')
            ylabel('# of electrodes')
            legend('S-','S+')
            xlim([-15 15])
            ylim([0 40])
            title(freq_names{bwii})
            
            pffft=1
            
            a={ delta_dB_powerpreHit(ROCbandwidthpre==bwii)' delta_dB_powerpreCR(ROCbandwidthpre==bwii)'};
            mode_statcond='perm';
            [F df pvals_perm(bwii)] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for premuted anovan dB delta power S+ vs S- ' freq_names{bwii} '= %d\n'],  pvals_perm(bwii));
            
        end
        
        pFDRanovan=drsFDRpval(pvals_perm);
        fprintf(1, ['pFDR for premuted anovan p value  = %d\n\n'],pFDRanovan);
        
        
        
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
            hold on
            n_cum=0;
            this_legend=[];
            
            histogram(auROCpre(( p_valROCpre>pFDRauROC)&(ROCbandwidthpre==bwii)),edges)
            histogram(auROCpre(( p_valROCpre<=pFDRauROC)&(ROCbandwidthpre==bwii)),edges)
            legend('auROC not singificant','auROC significant')
            title(['Histogram for ' freq_names{bwii} ' auROC for LFPs'])
            xlim([-0.2 0.6])
            ylim([0 35])
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
            bar(bwii,100*sum(( p_valROCpre<=pFDRauROC)&(ROCbandwidthpre==bwii))/sum((ROCbandwidthpre==bwii)))
        end
        title('Percent auROC significantly different from zero')
        
    case 4
        %Compare auROC in the last few trials of pre with first few trials of post
        % Generate figure 2 for Daniel's paper. first vs last.
        no_dBs=1;
        delta_dB_power_pre=[];
        no_ROCs=0;
        ROCoutpre=[];
        ROCoutpost=[];
        p_vals_ROC=[];
        delta_dB_powerpreHit=[];
        no_hits=0;
        pvals_auROCperm=[];
        pvals_dBperm=[];
        
        fprintf(1, ['Pairwise auROC analysis for Hit and CR LFP power\n\n'])
        p_vals=[];
        for fps=1:no_file_pairs
            for elec=1:16
                
                lfpodNopre_ref=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                lfpodNopost_ref=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                
                if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref)))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref)))
                    
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower))
                        
                        
                        trials_in_event_preHit=(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).which_eventLFPPower(evHit,:)==1);
                        trials_in_event_preCR=(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).which_eventLFPPower(evCR,:)==1);
                        trials_in_event_pre=logical(trials_in_event_preHit+trials_in_event_preCR);
                        
                        trials_in_event_postHit=(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).which_eventLFPPower(evHit,:)==1);
                        trials_in_event_postCR=(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).which_eventLFPPower(evCR,:)==1);
                        trials_in_event_post=logical(trials_in_event_postHit+trials_in_event_postCR);
                        
                        if (sum(trials_in_event_preHit)>=comp_window)&(sum(trials_in_event_preCR)>=comp_window)&(sum(trials_in_event_postHit)>=comp_window)&(sum(trials_in_event_postCR)>=comp_window)
                            if (sum(trials_in_event_pre)>=comp_window_auROC)&(sum(trials_in_event_post)>=comp_window_auROC)
                                lfpodNopre=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                lfpodNopost=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                
                                %pre Hits
                                
                                
                                no_trials=0;
                                jj=length(trials_in_event_preHit);
                                while no_trials<comp_window
                                    if trials_in_event_preHit(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    jj=jj-1;
                                end
                                jj=jj+1;
                                these_trials=logical([zeros(1,jj-1) ones(1,length(trials_in_event_preHit)-jj+1)].*trials_in_event_preHit);
                                
                                this_dB_powerprerefHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerprerefHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower(these_trials,:));
                                
                                this_dB_powerpreHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpreHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre).allPower(these_trials,:));
                                
                                
                                %pre CRs
                                no_trials=0;
                                jj=length(trials_in_event_preCR);
                                while no_trials<comp_window
                                    if trials_in_event_preCR(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    jj=jj-1;
                                end
                                jj=jj+1;
                                these_trials=logical([zeros(1,jj-1) ones(1,length(trials_in_event_preCR)-jj+1)].*trials_in_event_preCR);
                                
                                this_dB_powerprerefCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerprerefCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower(these_trials,:));
                                
                                this_dB_powerpreCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpreCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre).allPower(these_trials,:));
                                
                                
                                %post Hits
                                no_trials=0;
                                jj=0;
                                while no_trials<comp_window
                                    jj=jj+1;
                                    if trials_in_event_postHit(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    
                                end
                                
                                these_trials=logical([ones(1,jj) zeros(1,length(trials_in_event_postHit)-jj)].*trials_in_event_postHit);
                                
                                this_dB_powerpostrefHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostrefHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower(these_trials,:));
                                
                                this_dB_powerpostHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost).allPower(these_trials,:));
                                
                                
                                %post CRs
                                no_trials=0;
                                jj=0;
                                while no_trials<comp_window
                                    jj=jj+1;
                                    if trials_in_event_postCR(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    
                                end
                                
                                these_trials=logical([ones(1,jj) zeros(1,length(trials_in_event_postCR)-jj)].*trials_in_event_postCR);
                                
                                this_dB_powerpostrefCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostrefCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower(these_trials,:));
                                
                                this_dB_powerpostCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost).allPower(these_trials,:));
                                
                                for bwii=1:no_bandwidths
                                    
                                    no_ROCs=no_ROCs+1;
                                    this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                    
                                    %Enter the pre Hits
                                    this_delta_dB_powerpreHit=zeros(comp_window,1);
                                    this_delta_dB_powerpreHit=mean(this_dB_powerpreHit(:,this_band)-this_dB_powerprerefHit(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:comp_window,1)=this_delta_dB_powerpreHit;
                                    roc_data(1:comp_window,2)=zeros(comp_window,1);
                                    
                                    %Enter pre CR
                                    this_delta_dB_powerpreCR=zeros(comp_window,1);
                                    this_delta_dB_powerpreCR=mean(this_dB_powerpreCR(:,this_band)-this_dB_powerprerefCR(:,this_band),2);
                                    roc_data(comp_window+1:2*comp_window,1)=this_delta_dB_powerpreCR;
                                    roc_data(comp_window+1:2*comp_window,2)=ones(comp_window,1);
                                    
                                    
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
                                    this_delta_dB_powerpostHit=zeros(comp_window,1);
                                    this_delta_dB_powerpostHit=mean(this_dB_powerpostHit(:,this_band)-this_dB_powerpostrefHit(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:comp_window,1)=this_delta_dB_powerpostHit;
                                    roc_data(1:comp_window,2)=zeros(comp_window,1);
                                    
                                    %Enter post CR
                                    this_delta_dB_powerpostCR=zeros(comp_window,1);
                                    this_delta_dB_powerpostCR=mean(this_dB_powerpostCR(:,this_band)-this_dB_powerpostrefCR(:,this_band),2);
                                    roc_data(comp_window+1:2*comp_window,1)=this_delta_dB_powerpostCR;
                                    roc_data(comp_window+1:2*comp_window,2)=ones(comp_window,1);
                                    
                                    
                                    %Find post ROC
                                    ROCoutpost(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCoutpost(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopost_ref).fileNo;
                                    ROCgroupNopost(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).fileNo);
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
                                    if (groupNopre(no_dBs)==1)||(groupNopre(no_dBs)==3)
                                        %
                                        
                                        %Hit, all points
                                        delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                        delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                        figure(6)
                                        hold on
                                        plot([(bwii-1)*2 (bwii-1)*2+1],[delta_dB_powerpostHit(no_ROCs) delta_dB_powerpreHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        
                                        
                                        %CR, all points
                                        delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                        delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                        figure(7)
                                        hold on
                                        plot([(bwii-1)*2 (bwii-1)*2+1],[delta_dB_powerpostCR(no_ROCs) delta_dB_powerpreCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        
                                    end
                                end
                                
                                no_dBs=no_dBs+2;
                                
                            else
                                
                                if (sum(trials_in_event_pre)<comp_window)
                                    fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_pre),comp_window,file_pairs(fps,1),elec);
                                end
                                
                                if (sum(trials_in_event_post)<comp_window)
                                    fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_post),comp_window,file_pairs(fps,2),elec);
                                end
                            end
                            
                        else
                            
                            if (sum(trials_in_event_preHit)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_preHit),comp_window,file_pairs(fps,1),elec);
                            end
                            
                            if (sum(trials_in_event_preCR)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_preCR),comp_window,file_pairs(fps,1),elec);
                            end
                            if (sum(trials_in_event_postHit)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_postHit),comp_window,file_pairs(fps,2),elec);
                            end
                            if (sum(trials_in_event_postCR)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_postCR),comp_window,file_pairs(fps,2),elec);
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
            
            
            %Now plot the LFP power, all points
            %Hits for NL->L
            
            
            %All points - Hit, NL->L
            figure(6)
            mean_pre=mean(delta_dB_powerpreHit(((ROCgroupNopre==1)|(ROCgroupNopre==3))&logical(ROCbandwidthpre==bwii)));
            if sum(((ROCgroupNopre==1)|(ROCgroupNopre==3))&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpreHit(((ROCgroupNopre==1)|(ROCgroupNopre==3))&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            
            plot((bwii-1)*2+1,mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            plot([(bwii-1)*2+1 (bwii-1)*2+1],[mean_pre-ci95 mean_pre+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            
            mean_post=mean(delta_dB_powerpostHit(((ROCgroupNopre==1)|(ROCgroupNopre==3))&(ROCbandwidthpost==bwii)));
            if sum(((ROCgroupNopre==1)|(ROCgroupNopre==3))&(ROCbandwidthpost==bwii))>1
                pd=fitdist(delta_dB_powerpostHit(((ROCgroupNopre==1)|(ROCgroupNopre==3))&(ROCbandwidthpost==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            
            plot((bwii-1)*2,mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            plot([(bwii-1)*2 (bwii-1)*2],[mean_post-ci95 mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            
            plot([(bwii-1)*2 (bwii-1)*2+1],[mean_post mean_pre],'-','LineWidth',2,'Color',[0.6 0 0])
            
            %Do the statistics for dB differences for S+
            a={delta_dB_powerpreHit(((ROCgroupNopre==1)|(ROCgroupNopre==3))&logical(ROCbandwidthpre==bwii))' delta_dB_powerpostHit(((ROCgroupNopre==1)|(ROCgroupNopre==3))&(ROCbandwidthpost==bwii))'};
            mode_statcond='perm';
            [F df pval_dBperm] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for premuted anovan for dB S+ ' freq_names{bwii} '= %d\n'],  pval_dBperm);
            pvals_dBperm=[pvals_dBperm pval_dBperm];
            
            
            %CRs for NL->L
            figure(7)
            
            %All points - CR, NL->L
            mean_pre=mean(delta_dB_powerpreCR(((ROCgroupNopre==1)|(ROCgroupNopre==3))&(ROCbandwidthpre==bwii)));
            
            if sum(((ROCgroupNopre==1)|(ROCgroupNopre==3))&(ROCbandwidthpre==bwii))>1
                pd=fitdist(delta_dB_powerpreCR(((ROCgroupNopre==1)|(ROCgroupNopre==3))&(ROCbandwidthpre==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot((bwii-1)*2+1,mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            plot([(bwii-1)*2+1 (bwii-1)*2+1],[mean_pre-ci95 mean_pre+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            
            mean_post=mean(delta_dB_powerpostCR(((ROCgroupNopre==1)|(ROCgroupNopre==3))&(ROCbandwidthpost==bwii)));
            
            if sum(((ROCgroupNopre==1)|(ROCgroupNopre==3))&(ROCbandwidthpost==bwii))>1
                pd=fitdist(delta_dB_powerpostCR(((ROCgroupNopre==1)|(ROCgroupNopre==3))&(ROCbandwidthpost==bwii))','Normal');
                ci_tbl=paramci(pd);
                ci95=pd.mu-ci_tbl(1,1);
            else
                ci95=0;
            end
            plot((bwii-1)*2,mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            plot([(bwii-1)*2 (bwii-1)*2],[mean_post-ci95 mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            
            plot([(bwii-1)*2 (bwii-1)*2+1],[mean_post mean_pre],'-','LineWidth',2,'Color',[0 0 0.6])
            
            %Do the statistics for dB differences for S-
            a={delta_dB_powerpreCR(((ROCgroupNopre==1)|(ROCgroupNopre==3))&logical(ROCbandwidthpre==bwii))' delta_dB_powerpostCR(((ROCgroupNopre==1)|(ROCgroupNopre==3))&(ROCbandwidthpost==bwii))'};
            mode_statcond='perm';
            [F df pval_dBperm] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for premuted anovan for dB S- ' freq_names{bwii} '= %d\n'],  pval_dBperm);
            pvals_dBperm=[pvals_dBperm pval_dBperm];
        end
        
        pFDRdBperm=drsFDRpval(pvals_dBperm);
        fprintf(1, ['\npFDR for premuted anovan p value for difference between learning and proficent for dB = %d\n\n'],pFDRdBperm);
        
        
        figure(6)
        title(['Change in power dB from learning to proficient for ' evTypeLabels{1}])
        ylabel('delta Power (dB)')
        ylim([-10 12])
        
        figure(7)
        title(['Change in power dB from learning to proficient for ' evTypeLabels{2}])
        ylabel('delta Power (dB)')
        ylim([-10 12])
        
        
        fprintf(1, '\n\n')
        
        
        %Plot cumulative histos for auROCs
        dB_power_change=logical(dB_power_changeHit+dB_power_changeCR);
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
            maxauROC=max([max(auROCpre(ROCbandwidthpre==bwii)) max(auROCpost(ROCbandwidthpost==bwii))]);
            minauROC=min([min(auROCpre(ROCbandwidthpre==bwii)) min(auROCpost(ROCbandwidthpost==bwii))]);
            edges=[-0.5:0.05:0.5];
            pos2=[0.1 0.1 0.6 0.8];
            subplot('Position',pos2)
            hold on
            
            h2=histogram(auROCpost(ROCbandwidthpre==bwii),edges);
            h2.FaceColor='b';
            h1=histogram(auROCpre(ROCbandwidthpre==bwii),edges);
            h1.FaceColor='r';
            
            xlabel('auROC')
            ylabel('# of electrodes')
            legend('Learning','Proficient')
            title(['auROC for ' freq_names{bwii}])
            xlim([-0.3 0.6])
            ylim([0 30])
            
            %Plot the single electrodes
            pos2=[0.8 0.1 0.1 0.8];
            subplot('Position',pos2)
            hold on
            for ii=1:length(auROCpre)
                if ROCbandwidthpre(ii)==bwii
                    plot([0 1],[auROCpost(ii) auROCpre(ii)],'-o', 'Color',[0.7 0.7 0.7])
                end
            end
            
            
            plot([0 1],[mean(auROCpost(ROCbandwidthpre==bwii)) mean(auROCpre(ROCbandwidthpre==bwii))],'-k','LineWidth', 3)
            CI = bootci(1000, @mean, auROCpost(ROCbandwidthpre==bwii));
            plot([0 0],CI,'-r','LineWidth',3)
            plot(0,mean(auROCpost(ROCbandwidthpre==bwii)),'or','MarkerSize', 10,'MarkerFace','r')
            CI = bootci(1000, @mean, auROCpre(ROCbandwidthpre==bwii));
            plot([1 1],CI,'-b','LineWidth',3)
            plot(1,mean(auROCpre(ROCbandwidthpre==bwii)),'ob','MarkerSize', 10,'MarkerFace','b')
            ylabel('auROC')
            ylim([-0.2 0.5])
            
            %Do the statistics for auROC differences
            a={auROCpost(ROCbandwidthpre==bwii)' auROCpre(ROCbandwidthpre==bwii)'};
            mode_statcond='perm';
            [F df pval_auROCperm] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
            fprintf(1, ['p value for premuted anovan for auROC S+ vs S- ' freq_names{bwii} '= %d\n'],  pval_auROCperm);
            pvals_auROCperm=[pvals_auROCperm pval_auROCperm];
            
            %Figure 5
            figure(5)
            
            percent_auROCpost=100*sum(p_valROCpost(ROCbandwidthpre==bwii)<=pFDRROC)/sum(ROCbandwidthpre==bwii);
            bar(x,percent_auROCpost,'b')
            
            percent_auROCpre=100*sum(p_valROCpre(ROCbandwidthpre==bwii)<=pFDRROC)/sum(ROCbandwidthpre==bwii);
            bar(x+1,percent_auROCpre,'r')
            
            x=x+3;
            
        end
        
        figure(5)
        title('Percent singificant auROC')
        legend('Learning','Proficient')
        
        pFDRanovanauROC=drsFDRpval(pval_auROCperm);
        fprintf(1, ['\npFDR for premuted anovan p value for difference between learning and proficent for auROC = %d\n\n'],pFDRanovanauROC);
        
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
        
        fprintf(1, ['Pairwise auROC analysis for Hit and CR LFP power\n\n'])
        p_vals=[];
        for fps=1:no_file_pairs
            for elec=1:16
                
                lfpodNopre_ref=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                lfpodNopost_ref=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                
                if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref)))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref)))
                    
                    
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower))
                        
                        
                        trials_in_event_preHit=(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).which_eventLFPPower(evHit,:)==1);
                        trials_in_event_preCR=(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).which_eventLFPPower(evCR,:)==1);
                        trials_in_event_pre=logical(trials_in_event_preHit+trials_in_event_preCR);
                        
                        trials_in_event_postHit=(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).which_eventLFPPower(evHit,:)==1);
                        trials_in_event_postCR=(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).which_eventLFPPower(evCR,:)==1);
                        trials_in_event_post=logical(trials_in_event_postHit+trials_in_event_postCR);
                        
                        if (sum(trials_in_event_preHit)>=comp_window)&(sum(trials_in_event_preCR)>=comp_window)&(sum(trials_in_event_postHit)>=comp_window)&(sum(trials_in_event_postCR)>=comp_window)
                            if (sum(trials_in_event_pre)>=comp_window_auROC)&(sum(trials_in_event_post)>=comp_window_auROC)
                                lfpodNopre=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                lfpodNopost=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                
                                %pre Hits
                                
                                
                                no_trials=0;
                                jj=length(trials_in_event_preHit);
                                while no_trials<comp_window
                                    if trials_in_event_preHit(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    jj=jj-1;
                                end
                                jj=jj+1;
                                these_trials=logical([zeros(1,jj-1) ones(1,length(trials_in_event_preHit)-jj+1)].*trials_in_event_preHit);
                                
                                this_dB_powerprerefHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerprerefHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower(these_trials,:));
                                
                                this_dB_powerpreHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpreHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre).allPower(these_trials,:));
                                
                                
                                %pre CRs
                                no_trials=0;
                                jj=length(trials_in_event_preCR);
                                while no_trials<comp_window
                                    if trials_in_event_preCR(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    jj=jj-1;
                                end
                                jj=jj+1;
                                these_trials=logical([zeros(1,jj-1) ones(1,length(trials_in_event_preCR)-jj+1)].*trials_in_event_preCR);
                                
                                this_dB_powerprerefCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerprerefCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower(these_trials,:));
                                
                                this_dB_powerpreCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpreCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre).allPower(these_trials,:));
                                
                                
                                %post Hits
                                no_trials=0;
                                jj=0;
                                while no_trials<comp_window
                                    jj=jj+1;
                                    if trials_in_event_postHit(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    
                                end
                                
                                these_trials=logical([ones(1,jj) zeros(1,length(trials_in_event_postHit)-jj)].*trials_in_event_postHit);
                                
                                this_dB_powerpostrefHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostrefHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower(these_trials,:));
                                
                                this_dB_powerpostHit=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost).allPower(these_trials,:));
                                
                                
                                %post CRs
                                no_trials=0;
                                jj=0;
                                while no_trials<comp_window
                                    jj=jj+1;
                                    if trials_in_event_postCR(jj)==1
                                        no_trials=no_trials+1;
                                    end
                                    
                                end
                                
                                these_trials=logical([ones(1,jj) zeros(1,length(trials_in_event_postCR)-jj)].*trials_in_event_postCR);
                                
                                this_dB_powerpostrefCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostrefCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower(these_trials,:));
                                
                                this_dB_powerpostCR=zeros(comp_window,length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost).allPower(these_trials,:));
                                
                                for bwii=1:no_bandwidths
                                    
                                    no_ROCs=no_ROCs+1;
                                    this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                    
                                    %Enter the pre Hits
                                    this_delta_dB_powerpreHit=zeros(comp_window,1);
                                    this_delta_dB_powerpreHit=mean(this_dB_powerpreHit(:,this_band)-this_dB_powerprerefHit(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:comp_window,1)=this_delta_dB_powerpreHit;
                                    roc_data(1:comp_window,2)=zeros(comp_window,1);
                                    
                                    %Enter pre CR
                                    this_delta_dB_powerpreCR=zeros(comp_window,1);
                                    this_delta_dB_powerpreCR=mean(this_dB_powerpreCR(:,this_band)-this_dB_powerprerefCR(:,this_band),2);
                                    roc_data(comp_window+1:2*comp_window,1)=this_delta_dB_powerpreCR;
                                    roc_data(comp_window+1:2*comp_window,2)=ones(comp_window,1);
                                    
                                    
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
                                    this_delta_dB_powerpostHit=zeros(comp_window,1);
                                    this_delta_dB_powerpostHit=mean(this_dB_powerpostHit(:,this_band)-this_dB_powerpostrefHit(:,this_band),2);
                                    roc_data=[];
                                    roc_data(1:comp_window,1)=this_delta_dB_powerpostHit;
                                    roc_data(1:comp_window,2)=zeros(comp_window,1);
                                    
                                    %Enter post CR
                                    this_delta_dB_powerpostCR=zeros(comp_window,1);
                                    this_delta_dB_powerpostCR=mean(this_dB_powerpostCR(:,this_band)-this_dB_powerpostrefCR(:,this_band),2);
                                    roc_data(comp_window+1:2*comp_window,1)=this_delta_dB_powerpostCR;
                                    roc_data(comp_window+1:2*comp_window,2)=ones(comp_window,1);
                                    
                                    
                                    %Find post ROC
                                    ROCoutpost(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                    ROCoutpost(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopost_ref).fileNo;
                                    ROCgroupNopost(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).fileNo);
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
                                        %                                         %Hit, split
                                        %                                         if p_val(no_dBs,bwii)<=0.05
                                        %                                             delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                        %                                             delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                        %                                             figure(bwii)
                                        %                                             hold on
                                        %                                             plot([0 1],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        %                                         else
                                        %                                             delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                        %                                             delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                        %                                             figure(bwii)
                                        %                                             hold on
                                        %                                             plot([3 4],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        %
                                        %                                         end
                                        
                                        %Hit, all points
                                        delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                        delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
%                                         figure(bwii+12)
%                                         hold on
%                                         plot([0 1],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
%                                         
                                        %                                         %CR, split
                                        %                                         if p_val(no_dBs+1,bwii)<=0.05
                                        %                                             delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                        %                                             delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                        %                                             figure(bwii+4)
                                        %                                             hold on
                                        %                                             plot([0 1],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        %                                         else
                                        %                                             delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                        %                                             delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                        %                                             figure(bwii+4)
                                        %                                             hold on
                                        %                                             plot([3 4],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                        %                                         end
                                        
                                        %CR, all points
                                        delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                        delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
%                                         figure(bwii+4+12)
%                                         hold on
%                                         plot([0 1],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                    else
                                        if groupNopre(no_dBs)==3
                                            %Hit, split
                                            %                                             if p_val(no_dBs,bwii)<=0.05
                                            %                                                 delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                            %                                                 delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                            %                                                 figure(bwii)
                                            %                                                 hold on
                                            %                                                 plot([6 7],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            %                                             else
                                            %                                                 delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                            %                                                 delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
                                            %                                                 figure(bwii)
                                            %                                                 hold on
                                            %                                                 plot([9 10],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            %                                             end
                                            
                                            %Hit, all points
                                            delta_dB_powerpreHit(no_ROCs)=mean(this_delta_dB_powerpreHit);
                                            delta_dB_powerpostHit(no_ROCs)=mean(this_delta_dB_powerpostHit);
%                                             figure(bwii+12)
%                                             hold on
%                                             plot([3 4],[delta_dB_powerpreHit(no_ROCs) delta_dB_powerpostHit(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            
                                            %                                             %CR, split
                                            %                                             if p_val(no_dBs+1,bwii)<=0.05
                                            %                                                 delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                            %                                                 delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                            %                                                 figure(bwii+4)
                                            %                                                 hold on
                                            %                                                 plot([6 7],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            %                                             else
                                            %                                                 delta_dB_powerpreCR(no_ROCs)=mean(this_delta_dB_powerpreCR);
                                            %                                                 delta_dB_powerpostCR(no_ROCs)=mean(this_delta_dB_powerpostCR);
                                            %                                                 figure(bwii+4)
                                            %                                                 hold on
                                            %                                                 plot([9 10],[delta_dB_powerpreCR(no_ROCs) delta_dB_powerpostCR(no_ROCs)],'-o', 'Color',[0.7 0.7 0.7])
                                            %                                             end
                                            
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
                                
                                if (sum(trials_in_event_pre)<comp_window)
                                    fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_pre),comp_window,file_pairs(fps,1),elec);
                                end
                                
                                if (sum(trials_in_event_post)<comp_window)
                                    fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_post),comp_window,file_pairs(fps,2),elec);
                                end
                            end
                            
                        else
                            
                            if (sum(trials_in_event_preHit)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_preHit),comp_window,file_pairs(fps,1),elec);
                            end
                            
                            if (sum(trials_in_event_preCR)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_preCR),comp_window,file_pairs(fps,1),elec);
                            end
                            if (sum(trials_in_event_postHit)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_postHit),comp_window,file_pairs(fps,2),elec);
                            end
                            if (sum(trials_in_event_postCR)<comp_window)
                                fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_postCR),comp_window,file_pairs(fps,2),elec);
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
            %Hits for NL->L
            %             figure(bwii)
            %
            %             %Significant changes - Hit, NL->L
            %             sig_mean_pre=mean(delta_dB_powerpreHit((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreHit((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(0,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            %             plot([0 0],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             sig_mean_post=mean(delta_dB_powerpostHit((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
            %             if sum((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))>1
            %                 pd=fitdist(delta_dB_powerpostHit((dB_power_changeHit==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(1,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            %             plot([1 1],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             plot([0 1],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             %Changes not singificant - Hit, NL->L
            %             sig_mean_pre=mean(delta_dB_powerpreHit((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreHit((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(3,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            %             plot([3 3],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',3,'Color',[0.6 0 0])
            %
            %             sig_mean_post=mean(delta_dB_powerpostHit((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
            %             if sum((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))>1
            %                 pd=fitdist(delta_dB_powerpostHit((dB_power_changeHit==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(4,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
            %             plot([4 4],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             plot([3 4],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             %Significant changes - Hit, NLc->Lc
            %             sig_mean_pre=mean(delta_dB_powerpreHit((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreHit((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(6,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
            %             plot([6 6],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             sig_mean_post=mean(delta_dB_powerpostHit((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpostHit((dB_power_changeHit==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(7,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
            %             plot([7 7],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             plot([6 7],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             %Changes not singificant - Hit, NLc->Lc
            %             sig_mean_pre=mean(delta_dB_powerpreHit((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreHit((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %
            %             plot(9,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
            %             plot([9 9],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',3,'Color',[0 0 0.6])
            %
            %             sig_mean_post=mean(delta_dB_powerpostHit((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpostHit((dB_power_changeHit==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(10,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
            %             plot([10 10],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             plot([9 10],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             title(['Change in ' freq_names{bwii} ' power dB elicited by the laser for Hits'])
            %             xticks(0:10)
            %             xticklabels({'','light','','','light','','','light','','','light',''})
            %             xlim([-0.5 10.5])
            %             ylm=ylim;
            %             ylim([ylm(1)-(ylm(2)-ylm(1))/10 ylm(2)+(ylm(2)-ylm(1))/10])
            %             txt1 = 'DBH-Cre x ROSA26 Halo';
            %             text(1,ylm(1),txt1)
            %             txt1 = 'DBH-Cre';
            %             text(7.5,ylm(1),txt1)
            %
            %
            %             %CRs for NL->L
            %             figure(bwii+4)
            %
            %             %Significant changes - CR, NL->L
            %             sig_mean_pre=mean(delta_dB_powerpreCR((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
            %
            %             if sum((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreCR((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(0,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
            %             plot([0 0],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             sig_mean_post=mean(delta_dB_powerpostCR((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
            %
            %             if sum((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))>1
            %                 pd=fitdist(delta_dB_powerpostCR((dB_power_changeCR==1)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(1,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
            %             plot([1 1],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             plot([0 1],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             %Changes not singificant - CR, NL->L
            %             sig_mean_pre=mean(delta_dB_powerpreCR((dB_power_changeCR==0)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
            %
            %             if sum((dB_power_changeCR==1)&(ROCgroupNopre==0)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreCR((dB_power_changeCR==0)&(ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(3,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
            %             plot([3 3],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',3,'Color',[0.6 0 0])
            %
            %             sig_mean_post=mean(delta_dB_powerpostCR((dB_power_changeCR==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
            %
            %             if sum((dB_power_changeCR==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))
            %                 pd=fitdist(delta_dB_powerpostCR((dB_power_changeCR==0)&(ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(4,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
            %             plot([4 4],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             plot([3 4],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
            %
            %             %Significant changes - CR, NLc->Lc
            %             sig_mean_pre=mean(delta_dB_powerpreCR((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            %
            %             if sum((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreCR((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(6,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            %             plot([6 6],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             sig_mean_post=mean(delta_dB_powerpostCR((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpost==bwii)));
            %
            %             if sum((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpost==bwii))>1
            %                 pd=fitdist(delta_dB_powerpostCR((dB_power_changeCR==1)&(ROCgroupNopre==3)&(ROCbandwidthpost==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(7,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            %             plot([7 7],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             plot([6 7],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             %Changes not singificant - CR, NLc->Lc
            %             sig_mean_pre=mean(delta_dB_powerpreCR((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpreCR((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(9,sig_mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            %             plot([9 9],[sig_mean_pre-ci95 sig_mean_pre+ci95],'-','LineWidth',3,'Color',[0 0 0.6])
            %
            %             sig_mean_post=mean(delta_dB_powerpostCR((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
            %             if sum((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
            %                 pd=fitdist(delta_dB_powerpostCR((dB_power_changeCR==0)&(ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
            %                 ci_tbl=paramci(pd);
            %                 ci95=pd.mu-ci_tbl(1,1);
            %             else
            %                 ci95=0;
            %             end
            %             plot(10,sig_mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
            %             plot([10 10],[sig_mean_post-ci95 sig_mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             plot([9 10],[sig_mean_pre sig_mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
            %
            %             title(['Change in ' freq_names{bwii} ' power dB elicited by the laser for CRs'])
            %             xticks(0:10)
            %             xticklabels({'','light','','','light','','','light','','','light',''})
            %             xlim([-0.5 10.5])
            %             ylm=ylim;
            %             ylim([ylm(1)-(ylm(2)-ylm(1))/10 ylm(2)+(ylm(2)-ylm(1))/10])
            %             txt1 = 'DBH-Cre x ROSA26 Halo';
            %             text(1,ylm(1),txt1)
            %             txt1 = 'DBH-Cre';
            %             text(7.5,ylm(1),txt1)
            %
            %             pffft=1;
            %ANOVAn
            %Note: This is a first try, we should do a repeated measures...
            
%             %For Hits
%             all_delta_dB_power=[];
%             laser=[];
%             halo=[];
%             
%             %No laser, pre, with halo
%             all_delta_dB_power=[all_delta_dB_power delta_dB_powerpreHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))];
%             laser=[laser zeros(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
%             halo=[halo ones(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
%             
%             %Laser, post, with halo
%             all_delta_dB_power=[all_delta_dB_power delta_dB_powerpostHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))];
%             laser=[laser ones(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
%             halo=[halo ones(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
%             
%             %No laser, pre, no halo
%             all_delta_dB_power=[all_delta_dB_power delta_dB_powerpreHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))];
%             laser=[laser zeros(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
%             halo=[halo zeros(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
%             
%             %Laser, post, with halo
%             all_delta_dB_power=[all_delta_dB_power delta_dB_powerpostHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))];
%             laser=[laser ones(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
%             halo=[halo zeros(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
%             
%             p=anovan(all_delta_dB_power,{laser,halo},'model','interaction','varnames',{'laser','halo'},'display','off');
%             p_vals_anovan=[p_vals_anovan p(3)];
%             fprintf(1, ['anovan p value ' freq_names{bwii} ' and ' evTypeLabels{1} ' = %d\n'],p(3));
%             
%             %For CRs
%             all_delta_dB_power=[];
%             laser=[];
%             halo=[];
%             
%             %No laser, pre, with halo
%             all_delta_dB_power=[all_delta_dB_power delta_dB_powerpreCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))];
%             laser=[laser zeros(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
%             halo=[halo ones(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
%             
%             %Laser, post, with halo
%             all_delta_dB_power=[all_delta_dB_power delta_dB_powerpostCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))];
%             laser=[laser ones(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
%             halo=[halo ones(1,sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)))];
%             
%             %No laser, pre, no halo
%             all_delta_dB_power=[all_delta_dB_power delta_dB_powerpreCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))];
%             laser=[laser zeros(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
%             halo=[halo zeros(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
%             
%             %Laser, post, with halo
%             all_delta_dB_power=[all_delta_dB_power delta_dB_powerpostCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))];
%             laser=[laser ones(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
%             halo=[halo zeros(1,sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))];
%             
%             p=anovan(all_delta_dB_power,{laser,halo},'model','interaction','varnames',{'laser','halo'},'display','off');
%             p_vals_anovan=[p_vals_anovan p(3)];
%             fprintf(1, ['anovan p value ' freq_names{bwii} ' and ' evTypeLabels{2} ' = %d\n\n'],p(3));
            
            
%             %Now plot the LFP power, all points
%             %Hits for NL->L
%             figure(bwii+12)
%             
%             %All points - Hit, NL->L
%             mean_pre=mean(delta_dB_powerpreHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
%             if sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))>1
%                 pd=fitdist(delta_dB_powerpreHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
%                 ci_tbl=paramci(pd);
%                 ci95=pd.mu-ci_tbl(1,1);
%             else
%                 ci95=0;
%             end
%             plot(0,mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
%             plot([0 0],[mean_pre-ci95 mean_pre+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
%             
%             mean_post=mean(delta_dB_powerpostHit((ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
%             if sum((ROCgroupNopre==1)&(ROCbandwidthpost==bwii))>1
%                 pd=fitdist(delta_dB_powerpostHit((ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
%                 ci_tbl=paramci(pd);
%                 ci95=pd.mu-ci_tbl(1,1);
%             else
%                 ci95=0;
%             end
%             plot(1,mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0],'MarkerEdgeColor',[0.6 0 0])
%             plot([1 1],[mean_post-ci95 mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
%             
%             plot([0 1],[mean_pre mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
%             
%             
%             
%             %All points - Hit, NLc->Lc
%             mean_pre=mean(delta_dB_powerpreHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
%             if sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
%                 pd=fitdist(delta_dB_powerpreHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
%                 ci95=pd.mu-ci_tbl(1,1);
%             else
%                 ci95=0;
%             end
%             plot(3,mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
%             plot([3 3],[mean_pre-ci95 mean_pre+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
%             
%             mean_post=mean(delta_dB_powerpostHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
%             if sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
%                 pd=fitdist(delta_dB_powerpostHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
%                 ci_tbl=paramci(pd);
%                 ci95=pd.mu-ci_tbl(1,1);
%             else
%                 ci95=0;
%             end
%             ci95=pd.mu-ci_tbl(1,1);
%             plot(4,mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6],'MarkerEdgeColor',[0 0 0.6])
%             plot([4 4],[mean_post-ci95 mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
%             
%             plot([3 4],[mean_pre mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
%             
%             
%             title(['Change in ' freq_names{bwii} ' power dB elicited by the laser for ' evTypeLabels{1}])
%             xticks(0:10)
%             xticklabels({'','light','','','light',''})
%             xlim([-0.5 4.5])
%             ylm=ylim;
%             ylim([ylm(1)-(ylm(2)-ylm(1))/10 ylm(2)+(ylm(2)-ylm(1))/10])
%             txt1 = 'DBH-Cre x ROSA26 Halo';
%             text(0,ylm(1),txt1)
%             txt1 = 'DBH-Cre';
%             text(3,ylm(1),txt1)
%             
%             
%             %CRs for NL->L
%             figure(bwii+4+12)
%             
%             %All points - CR, NL->L
%             mean_pre=mean(delta_dB_powerpreCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)));
%             
%             if sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))>1
%                 pd=fitdist(delta_dB_powerpreCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))','Normal');
%                 ci_tbl=paramci(pd);
%                 ci95=pd.mu-ci_tbl(1,1);
%             else
%                 ci95=0;
%             end
%             plot(0,mean_pre,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
%             plot([0 0],[mean_pre-ci95 mean_pre+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
%             
%             mean_post=mean(delta_dB_powerpostCR((ROCgroupNopre==1)&(ROCbandwidthpost==bwii)));
%             
%             if sum((ROCgroupNopre==1)&(ROCbandwidthpost==bwii))>1
%                 pd=fitdist(delta_dB_powerpostCR((ROCgroupNopre==1)&(ROCbandwidthpost==bwii))','Normal');
%                 ci_tbl=paramci(pd);
%                 ci95=pd.mu-ci_tbl(1,1);
%             else
%                 ci95=0;
%             end
%             plot(1,mean_post,'o','MarkerSize', 10,'MarkerFace',[0.6 0 0])
%             plot([1 1],[mean_post-ci95 mean_post+ci95],'-','LineWidth',2,'Color',[0.6 0 0])
%             
%             plot([0 1],[mean_pre mean_post],'-','LineWidth',2,'Color',[0.6 0 0])
%             
%             
%             
%             %All points - CR, NLc->Lc
%             mean_pre=mean(delta_dB_powerpreCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)));
%             
%             if sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))>1
%                 pd=fitdist(delta_dB_powerpreCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))','Normal');
%                 ci_tbl=paramci(pd);
%                 ci95=pd.mu-ci_tbl(1,1);
%             else
%                 ci95=0;
%             end
%             plot(3,mean_pre,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
%             plot([3 3],[mean_pre-ci95 mean_pre+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
%             
%             mean_post=mean(delta_dB_powerpostCR((ROCgroupNopre==3)&(ROCbandwidthpost==bwii)));
%             
%             if sum((ROCgroupNopre==3)&(ROCbandwidthpost==bwii))>1
%                 pd=fitdist(delta_dB_powerpostCR((ROCgroupNopre==3)&(ROCbandwidthpost==bwii))','Normal');
%                 ci_tbl=paramci(pd);
%                 ci95=pd.mu-ci_tbl(1,1);
%             else
%                 ci95=0;
%             end
%             plot(4,mean_post,'o','MarkerSize', 10,'MarkerFace',[0 0 0.6])
%             plot([4 4],[mean_post-ci95 mean_post+ci95],'-','LineWidth',2,'Color',[0 0 0.6])
%             
%             plot([3 4],[mean_pre mean_post],'-','LineWidth',2,'Color',[0 0 0.6])
%             
%             
%             
%             title(['Change in ' freq_names{bwii} ' power dB elicited by the laser for ' evTypeLabels{2}])
%             xticks(0:10)
%             xticklabels({'','light','','','light',''})
%             xlim([-0.5 4.5])
%             ylm=ylim;
%             ylim([ylm(1)-(ylm(2)-ylm(1))/10 ylm(2)+(ylm(2)-ylm(1))/10])
%             txt1 = 'DBH-Cre x ROSA26 Halo';
%             text(0,ylm(1),txt1)
%             txt1 = 'DBH-Cre';
%             text(3,ylm(1),txt1)
%             
%             
%             %Analysis of covariance
%             
%             %For Hit
%             this_delta_dB_pre_powerHit=[];
%             this_delta_dB_pre_powerHit=delta_dB_powerpreHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))';
%             this_delta_dB_pre_powerHit=[this_delta_dB_pre_powerHit; delta_dB_powerpreHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))'];
%             
%             
%             this_delta_dB_post_powerHit=[];
%             this_delta_dB_post_powerHit=delta_dB_powerpostHit((ROCgroupNopre==1)&(ROCbandwidthpost==bwii))';
%             this_delta_dB_post_powerHit=[this_delta_dB_post_powerHit; delta_dB_powerpostHit((ROCgroupNopre==3)&(ROCbandwidthpost==bwii))'];
%             
%             pre_post=[];
%             pre_post=[zeros(sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)),1); ones(sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)),1)];
%             
%             [h,atab,ctab,stats] = aoctool(this_delta_dB_pre_powerHit,this_delta_dB_post_powerHit,pre_post,0.05,'','','','off');
%             
%             
%             pvals_ancova=[pvals_ancova atab{4,6}];
%             fprintf(1, ['ancova delta dB p value ' freq_names{bwii} ' and ' evTypeLabels{1} ' = %d\n'],atab{4,6});
%             
%             %Do ancova figures for Hits
%             figure(20+bwii)
%             h1=plot(delta_dB_powerpreHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)),delta_dB_powerpostHit((ROCgroupNopre==1)&(ROCbandwidthpost==bwii)),'or');
%             hold on
%             h2=plot(delta_dB_powerpreHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)),delta_dB_powerpostHit((ROCgroupNopre==3)&(ROCbandwidthpost==bwii)),'ob');
%             
%             slope_pre=ctab{5,2}+ctab{6,2};
%             int_pre=ctab{2,2}+ctab{3,2};
%             min_x=min([min(delta_dB_powerpreHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) min(delta_dB_powerpreHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
%             max_x=max([max(delta_dB_powerpreHit((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) max(delta_dB_powerpreHit((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
%             x=[min_x max_x];
%             plot(x,slope_pre*x+int_pre,'-r')
%             
%             slope_post=ctab{5,2}+ctab{7,2};
%             int_post=ctab{2,2}+ctab{4,2};
%             x=[min_x max_x];
%             plot(x,slope_post*x+int_post,'-b')
%             
%             title(['post vs pre power dB for ' freq_names{bwii} ' and ' evTypeLabels{1}])
%             xlabel('pre dB')
%             ylabel('post dB')
%             legend([h1 h2],'halo','no halo')
%             
%             %For CR
%             this_delta_dB_pre_powerCR=[];
%             this_delta_dB_pre_powerCR=delta_dB_powerpreCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))';
%             this_delta_dB_pre_powerCR=[this_delta_dB_pre_powerCR; delta_dB_powerpreCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii))'];
%             
%             
%             this_delta_dB_post_powerCR=[];
%             this_delta_dB_post_powerCR=delta_dB_powerpostCR((ROCgroupNopre==1)&(ROCbandwidthpost==bwii))';
%             this_delta_dB_post_powerCR=[this_delta_dB_post_powerCR; delta_dB_powerpostCR((ROCgroupNopre==3)&(ROCbandwidthpost==bwii))'];
%             
%             pre_post=[];
%             pre_post=[zeros(sum((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)),1); ones(sum((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)),1)];
%             
%             [h,atab,ctab,stats] = aoctool(this_delta_dB_pre_powerCR,this_delta_dB_post_powerCR,pre_post,0.05,'','','','off');
%             
%             
%             pvals_ancova=[pvals_ancova atab{4,6}];
%             fprintf(1, ['ancova delta dB p value ' freq_names{bwii} ' and ' evTypeLabels{2} ' = %d\n\n'],atab{4,6});
%             
%             %Do ancova figures for CR
%             figure(24+bwii)
%             h1=plot(delta_dB_powerpreCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)),delta_dB_powerpostCR((ROCgroupNopre==1)&(ROCbandwidthpost==bwii)),'or');
%             hold on
%             h2=plot(delta_dB_powerpreCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)),delta_dB_powerpostCR((ROCgroupNopre==3)&(ROCbandwidthpost==bwii)),'ob');
%             
%             slope_pre=ctab{5,2}+ctab{6,2};
%             int_pre=ctab{2,2}+ctab{3,2};
%             min_x=min([min(delta_dB_powerpreCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) min(delta_dB_powerpreCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
%             max_x=max([max(delta_dB_powerpreCR((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) max(delta_dB_powerpreCR((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
%             x=[min_x max_x];
%             plot(x,slope_pre*x+int_pre,'-r')
%             
%             slope_post=ctab{5,2}+ctab{7,2};
%             int_post=ctab{2,2}+ctab{4,2};
%             x=[min_x max_x];
%             plot(x,slope_post*x+int_post,'-b')
%             
%             title(['post vs pre power dB for ' freq_names{bwii} ' and ' evTypeLabels{2}])
%             xlabel('pre dB')
%             ylabel('post dB')
%             legend([h1 h2],'halo','no halo')
            
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
            h1=plot(auROCpre((ROCgroupNopre==1)&(ROCbandwidthpre==bwii)),auROCpost((ROCgroupNopre==1)&(ROCbandwidthpost==bwii)),'or');
            hold on
            h2=plot(auROCpre((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)),auROCpost((ROCgroupNopre==3)&(ROCbandwidthpost==bwii)),'ob');
            
            slope_pre=ctab{5,2}+ctab{6,2};
            int_pre=ctab{2,2}+ctab{3,2};
            min_x=min([min(auROCpre((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) min(auROCpre((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
            max_x=max([max(auROCpre((ROCgroupNopre==1)&(ROCbandwidthpre==bwii))) max(auROCpre((ROCgroupNopre==3)&(ROCbandwidthpre==bwii)))]);
            x=[-0.2 0.5];
            plot(x,slope_pre*x+int_pre,'-r','LineWidth',2)
            
            slope_post=ctab{5,2}+ctab{7,2};
            int_post=ctab{2,2}+ctab{4,2};
            x=[min_x max_x];
            plot(x,slope_post*x+int_post,'-b','LineWidth',2)
            
            plot([-0.2 0.5],[-0.2 0.5],'-k','LineWidth',2)
            
            title(['post vs pre auROC for ' freq_names{bwii} ])
            xlabel('pre auROC')
            ylabel('post auROC')
            legend([h1 h2],'halo','no halo')
            xlim([-0.2 0.5])
            ylim([-0.2 0.5])
            
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
        try
            close(9)
        catch
        end
        figure(9)
        try
            close(10)
        catch
        end
        figure(10)
        x=0
        for bwii=1:4
            n_cum=0;
            this_legend=[];
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
                ylim([0 30])
                
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
                
                %Do the statistics for auROC differences
                a={auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))' auROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))'};
                mode_statcond='perm';
                [F df pval_auROCperm] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
                if grs==1
                    fprintf(1, ['p value for premuted anovan for auROC DBH Cre x halo pre vs laser ' freq_names{bwii} '= %d\n'],  pval_auROCperm);
                else
                    fprintf(1, ['p value for premuted anovan for auROC DBH Cre pre vs laser ' freq_names{bwii} '= %d\n'],  pval_auROCperm);
                end
                pvals_auROCperm=[pvals_auROCperm pval_auROCperm];
                
                %Figure 5
                if grs==1
                    figure(9)
                else
                    figure(10)
                end
                hold on
                
                percent_auROCpost=100*sum(p_valROCpost((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))<=pFDRROC)/sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii));
                bar(x,percent_auROCpost,'b')
                
                percent_auROCpre=100*sum(p_valROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii))<=pFDRROC)/sum((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii));
                bar(x+1,percent_auROCpre,'r')
         
            end
            figNo=figNo+2;
             x=x+3;
        end
        
        pFDRauROC=drsFDRpval(pvals_auROCperm);
        fprintf(1, ['pFDR for auROC  = %d\n\n'],pFDRauROC);
        
      pffft=1;
end
pffft=1

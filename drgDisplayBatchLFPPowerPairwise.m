function drgDisplayBatchLFPPowerPairwise(handles)

%This function displays the LFP power spectrum for drgRunBatch-generated
%data



%Do you wnat to subtract the LFP spectrum in a reference window?
subtractRef=1;
refWin=1;

grRef=1;  %p values will be calculated with respect to this group

close all
warning('off')


%Choose the format for the display and the event types
%THESE VALUES ARE IMPORTANT


% which_display
%
% 1 Find the electrodes displaying a post vs. pre difference
%
% 2 auROC analysis for post vs. pre differences
%


% For Daniel's tabaproprio
% winNo=2;
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
% file_pairs=[1 7;
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

% 
% % For Daniel's isoamyl acetate tstart
winNo=2;
which_display=2;
% eventType=[2 5];
% evTypeLabels={'Hit','CR'};
eventType=[3 6];
evTypeLabels={'S+','S-'};

%Experiment pairs
%Important: The first file must be the experiment performed first
%For example in acetophenone ethyl benzoate no laser is first, laser is
%second
file_pairs=[1 7;
    1 14;
    2 12;
    3 13;
    4 15;
    5 16;
    6 17;
    7 18;
    8 19;
    9 20;
    10 21;
    11 22];
no_file_pairs=12;

comp_window=10; %works well with 8-12
comp_window_auROC=20;

grpre=[1 3];
grpost=[2 4];

% For Daniel's acetophenone ethyl benzoate
% winNo=2;
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
% file_pairs=[1 7;
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

% For Daniel's Hit, CR, FA, compare events
% winNo=2;
% which_display=1;
% eventType=[2 5];
% evTypeLabels={'Hit','CR'};

%THESE VALUES ARE IMPORTANT
%VERY IMPORTANT: This is the index for this event in your handles.drgbchoices.evTypeNos
% eventType=[2 5]; %Hit and CR
% evTypeLabels={'Hit';'CR'};


% % For Daniel's delta power Hit - CR, compare groups
% winNo=2;
% which_display=8;
% eventType1=2;
% eventType=2;
% eventTypeRef=5;
% evTypeLabel='Hit';
% evTypeRefLabel='CR';





% eventType=[9 10 11 12 13 14]; %Hit and CR
% no_event_types=6;
% evTypeLabels={'Hi Od1';'Hi Od 2';'Hi Od 3';'Low Od1';'Low Od 2';'Low Od 3'};

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

% evHit=2;
% evCR=5;

evHit=3;
evCR=6;

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
        %Compare last few trials of pre with first few trials of post
        no_dBs=0;
        delta_dB_power_pre=[];
        
        for evTN1=1:length(eventType)
            eventType1=eventType(evTN1);
            fprintf(1, ['Pairwise LFP power analysis for event: ' evTypeLabels{evTN1} '\n\n'])
            p_vals=[];
            for fps=1:no_file_pairs
                for elec=1:16
                    
                    lfpodNopre_ref=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                    lfpodNopost_ref=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==refWin));
                    if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref)))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref)))
                        
                        if (~isempty(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower))&(~isempty(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower))
                            
                            trials_in_event_pre=handles_drgb.drgb.lfpevpair(lfpodNopre_ref).which_eventLFPPower(eventType1,:)==1;
                            trials_in_event_post=handles_drgb.drgb.lfpevpair(lfpodNopost_ref).which_eventLFPPower(eventType1,:)==1;
                            
                            if (sum(trials_in_event_pre)>=comp_window)&(sum(trials_in_event_post)>=comp_window)
                                
                                no_dBs=no_dBs+1;
                                
                                this_dB_powerpreref=zeros(sum(trials_in_event_pre),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpreref(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower(trials_in_event_pre,:));
                                
                                lfpodNopre=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                this_dB_powerpre=zeros(sum(trials_in_event_pre),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpre(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre).allPower(trials_in_event_pre,:));
                                
                                for bwii=1:4
                                    no_dB_power_pre(no_dBs)=sum(trials_in_event_pre);
                                    delta_dB_power_pre(no_dBs,1:sum(trials_in_event_pre),bwii)=mean(this_dB_powerpre(:,(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii)))...
                                        -this_dB_powerpreref(:,(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii))),2);
                                end
                                
                                groupNopre(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,1));
                                
                                %Post file
                                
                                this_dB_powerpostref=zeros(sum(trials_in_event_post),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpostref(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower(trials_in_event_post,:));
                                
                                lfpodNopost=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                                this_dB_powerpost=zeros(sum(trials_in_event_post),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_powerpost(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost).allPower(trials_in_event_post,:));
                                
                                for bwii=1:4
                                    no_dB_power_post(no_dBs)=sum(trials_in_event_post);
                                    delta_dB_power_post(no_dBs,1:sum(trials_in_event_post),bwii)=mean(this_dB_powerpost(:,(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii)))...
                                        -this_dB_powerpostref(:,(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii))),2);
                                end
                                
                                groupNopost(no_dBs)=handles_drgb.drgbchoices.group_no(file_pairs(fps,2));
                                
                                %Are they different?
                                for bwii=1:4
                                    this_dB_power_pre=zeros(1,sum(trials_in_event_pre));
                                    this_dB_power_pre(1,:)=delta_dB_power_pre(no_dBs,1:sum(trials_in_event_pre),bwii);
                                    this_dB_power_post=zeros(1,sum(trials_in_event_post));
                                    this_dB_power_post(1,:)=delta_dB_power_post(no_dBs,1:sum(trials_in_event_post),bwii);
                                    
                                    p_val(no_dBs,bwii)=ranksum(this_dB_power_pre,this_dB_power_post);
                                    p_vals=[p_vals p_val(no_dBs,bwii)];
                                end
                                events(no_dBs)=evTN1;
                            else
                                
                                if (sum(trials_in_event_pre)<comp_window)
                                    fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_pre),comp_window,file_pairs(fps,1),elec);
                                end
                                
                                if (sum(trials_in_event_post)<comp_window)
                                    fprintf(1, ['%d trials in event fewer than comp_window %d for file No %d electrode %d\n'],sum(trials_in_event_post),comp_window,file_pairs(fps,2),elec);
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
                            fprintf(1, ['Empty lfpevpair for file No %d electrode %d\n'],file_pairs(fps,2),elec);
                        end
                        
                    end
                    
                end
                
            end
            fprintf(1, '\n\n')
            for bwii=1:4
                pFDR(evTN1,bwii)=drsFDRpval(p_val(:,bwii));
                fprintf(1, ['pFDR for ' freq_names{bwii} ' in ' evTypeLabels{evTN1} ' = %d\n'],pFDR(evTN1,bwii));
            end
            fprintf(1, '\n\n')
        end
        
        
        fprintf(1, '\n\n')
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
                        
                        if (sum(trials_in_event_pre)>=comp_window)&(sum(trials_in_event_post)>=comp_window_auROC)
                            
                            lfpodNopre=find((files_per_lfp==file_pairs(fps,1))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                            lfpodNopost=find((files_per_lfp==file_pairs(fps,2))&(elec_per_lfp==elec)&(window_per_lfp==winNo));
                            
                            %pre Hits
                            this_dB_powerprerefHit=zeros(sum(trials_in_event_preHit),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_powerprerefHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower(trials_in_event_preHit,:));
                            
                            this_dB_powerpreHit=zeros(sum(trials_in_event_preHit),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_powerpreHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre).allPower(trials_in_event_preHit,:));
                            
                            
                            %pre CRs
                            this_dB_powerprerefCR=zeros(sum(trials_in_event_preCR),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_powerprerefCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).allPower(trials_in_event_preCR,:));
                            
                            this_dB_powerpreCR=zeros(sum(trials_in_event_preCR),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_powerpreCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopre).allPower(trials_in_event_preCR,:));
                            
                            
                            %post Hits
                            this_dB_powerpostrefHit=zeros(sum(trials_in_event_postHit),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_powerpostrefHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower(trials_in_event_postHit,:));
                            
                            this_dB_powerpostHit=zeros(sum(trials_in_event_postHit),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_powerpostHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost).allPower(trials_in_event_postHit,:));
                            
                            
                            %post CRs
                            this_dB_powerpostrefCR=zeros(sum(trials_in_event_postCR),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_powerpostrefCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost_ref).allPower(trials_in_event_postCR,:));
                            
                            this_dB_powerpostCR=zeros(sum(trials_in_event_postCR),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_powerpostCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNopost).allPower(trials_in_event_postCR,:));
                            
                            for bwii=1:no_bandwidths
                                
                                no_ROCs=no_ROCs+1;
                                this_band=(frequency>=low_freq(bwii))&(frequency<=high_freq(bwii));
                                
                                %Enter the pre Hits
                                this_delta_dB_powerpreHit=zeros(sum(trials_in_event_preHit),1);
                                this_delta_dB_powerpreHit=mean(this_dB_powerpreHit(:,this_band)-this_dB_powerprerefHit(:,this_band),2);
                                roc_data=[];
                                roc_data(1:sum(trials_in_event_preHit),1)=this_delta_dB_powerpreHit;
                                roc_data(1:sum(trials_in_event_preHit),2)=zeros(sum(trials_in_event_preHit),1);
                                
                                %Enter pre CR
                                this_delta_dB_powerpreCR=zeros(sum(trials_in_event_preCR),1);
                                this_delta_dB_powerpreCR=mean(this_dB_powerpreCR(:,this_band)-this_dB_powerprerefCR(:,this_band),2);
                                roc_data(sum(trials_in_event_preHit)+1:sum(trials_in_event_preHit)+sum(trials_in_event_preCR),1)=this_delta_dB_powerpreCR;
                                roc_data(sum(trials_in_event_preHit)+1:sum(trials_in_event_preHit)+sum(trials_in_event_preCR),2)=ones(sum(trials_in_event_preCR),1);
                                
                                
                                %Find pre ROC
                                ROCoutpre(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
                                ROCoutpre(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNopre_ref).fileNo;
                                ROCgroupNopre(no_ROCs)=handles_drgb.drgbchoices.group_no(handles_drgb.drgb.lfpevpair(lfpodNopre_ref).fileNo);
                                ROCoutpre(no_ROCs).timeWindow=winNo;
                                ROCbandwidthpre(no_ROCs)=bwii;
                                auROCpre(no_ROCs)=ROCoutpre(no_ROCs).roc.AUC-0.5;
                                
                                p_vals_ROC=[p_vals_ROC ROCoutpre(no_ROCs).roc.p];
                                
                                %Enter the post Hits
                                this_delta_dB_powerpostHit=zeros(sum(trials_in_event_postHit),1);
                                this_delta_dB_powerpostHit=mean(this_dB_powerpostHit(:,this_band)-this_dB_powerpostrefHit(:,this_band),2);
                                roc_data=[];
                                roc_data(1:sum(trials_in_event_postHit),1)=this_delta_dB_powerpostHit;
                                roc_data(1:sum(trials_in_event_postHit),2)=zeros(sum(trials_in_event_postHit),1);
                                
                                %Enter post CR
                                this_delta_dB_powerpostCR=zeros(sum(trials_in_event_postCR),1);
                                this_delta_dB_powerpostCR=mean(this_dB_powerpostCR(:,this_band)-this_dB_powerpostrefCR(:,this_band),2);
                                roc_data(sum(trials_in_event_postHit)+1:sum(trials_in_event_postHit)+sum(trials_in_event_postCR),1)=this_delta_dB_powerpostCR;
                                roc_data(sum(trials_in_event_postHit)+1:sum(trials_in_event_postHit)+sum(trials_in_event_postCR),2)=ones(sum(trials_in_event_postCR),1);
                                
                                
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
                                
                                if (p_val(no_dBs,bwii)<0.05)||(p_val(no_dBs+1,bwii)<0.05)
                                    dB_power_change(no_ROCs)=1;
                                else
                                    dB_power_change(no_ROCs)=0;
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
        
        
        
        
        fprintf(1, '\n\n')
        
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
        
        %Plot cumulative histos for auROCs
        
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
                [f_auROC,x_auROC] = drg_ecdf(auROCpre((ROCgroupNopre==grpre(grs))&(ROCbandwidthpre==bwii)&(dB_power_change==1)));
                n_cum=n_cum+1;
                plot(x_auROC,f_auROC,these_lines{n_cum})
                this_legend=[this_legend '''' handles_drgb.drgbchoices.group_no_names{grpre(grs)}  '''' ','];
                
                %post
                [f_auROC,x_auROC] = drg_ecdf(auROCpost((ROCgroupNopost==grpost(grs))&(ROCbandwidthpost==bwii)&(dB_power_change==1)));
                n_cum=n_cum+1;
                plot(x_auROC,f_auROC,these_lines{n_cum})
                this_legend=[this_legend '''' handles_drgb.drgbchoices.group_no_names{grpost(grs)}  '''' ','];
            end
            this_legend=['legend(' this_legend(1:end-1) ')'];
            eval(this_legend)
            title(['Cumulative histogram for ' freq_names{bwii} ' auROC for LFPs with significant power changes'])
        end
        
        % auROCpre(no_ROCs)=ROCoutpre(no_ROCs).roc.AUC-0.5;
end
pffft=1

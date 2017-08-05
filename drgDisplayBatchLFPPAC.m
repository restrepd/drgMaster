function drgDisplayBatchLFPPAC(handles)

%This function displays the results of PAC analysis for drgRunBatch-generated
%data

close all
clear all
warning('off')

%THESE VALUES ARE IMPORTANT

%Choose the format for the display
which_PA_display=1;

% 1 Show the cumulative histograms for peak angles for the difference between events
%   for each group/ percent correct bin in a separate graph
%
% 2 Show the rose plots for the peak angles
%
% 3 Show the phase angle histograms
%
% 4 show the boxplots for the MI
%
% 5 Show the cumulative histos for the MI
%
% 6 Show the cumulative histograms for mean phase angles for the difference between events
%   for each group/ percent correct bin in a separate graph 

%For the moment let's only do the odor window
winNo=2;

%Analyze hit, miss, CR and FA
% evTypeNos=[2 4 5 7];
% evTypeLabels={'Hit','Miss','CR','FA'};

%Analyze Hit CR
evTypeNos=[2 5];
evTypeLabels={'Hit','CR'};

ev_to_test_against=1;



%Ask user for the drgb output .mat file and load those data
[handles.drgb.outFileName,handles.PathName] = uigetfile('*.mat','Select the drgb output file');
load([handles.PathName handles.drgb.outFileName])


%Initialize the variables

%Percent bins
percent_low=[45 65 80];
percent_high=[65 80 100];
no_percent_bins=3;

% percent_low=[40 51 61 71 81 91];
% percent_high=[50 60 70 80 90 100];
% no_percent_bins=6;

%Determine how many trials for each evTypeNo, time window, group, etc
no_groups=max(handles_drgb.drgbchoices.group_no);
no_events=length(handles_drgb.drgbchoices.evTypeNos);
no_lfpevpairs=handles_drgb.drgb.lfpevpair_no;
no_PACpeaks=handles_drgb.no_PACpeaks;
no_trials=zeros(length(handles_drgb.drgbchoices.evTypeNos),handles_drgb.drgbchoices.noWindows,max(handles_drgb.drgbchoices.group_no),no_percent_bins,handles_drgb.no_PACpeaks);

%In the previous versions of the code these variables were not defined
if isfield(handles_drgb.drgb.lfpevpair,'which_eventPAC')==0
    for lfpodNo=1:no_lfpevpairs
        handles_drgb.drgb.lfpevpair(lfpodNo).which_eventPAC=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).which_eventPAC;
        handles_drgb.drgb.lfpevpair(lfpodNo).perCorrPAC=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).perCorrPAC;
    end
end

for lfpodNo=1:no_lfpevpairs
    fileNo=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
    groupNo=handles_drgb.drgbchoices.group_no(fileNo);
    timeWindow=handles_drgb.drgb.lfpevpair(lfpodNo).timeWindow;
    for evTypeNo=1:no_events
        trials_in_event=handles_drgb.drgb.lfpevpair(lfpodNo).which_eventPAC(evTypeNo,:)==1;
        if sum(trials_in_event)>0
            for percent_bin=1:no_percent_bins
                trials_in_perbin=(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrPAC>percent_low(percent_bin))&(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrPAC<=percent_high(percent_bin));
                if sum(trials_in_event&trials_in_perbin)>0
                    for PACno=1:no_PACpeaks
                        no_trials(evTypeNo,timeWindow,groupNo,percent_bin,PACno)=no_trials(evTypeNo,timeWindow,groupNo,percent_bin,PACno)+sum(trials_in_event&trials_in_perbin);
                    end
                end
            end
        end
    end
end
max_trials=max(no_trials(:));

%Now initialize the variables for MI and phase angle
noWindows=handles_drgb.drgbchoices.noWindows;
mi_values=zeros(no_events,noWindows,no_groups,no_percent_bins,no_PACpeaks,no_lfpevpairs);
pa_values=zeros(no_events,noWindows,no_groups,no_percent_bins,no_PACpeaks,no_lfpevpairs);
mva_values=zeros(no_events,noWindows,no_groups,no_percent_bins,no_PACpeaks,no_lfpevpairs);



no_values=zeros(no_events,noWindows,no_groups,no_percent_bins,no_PACpeaks);

sz_all_phase=size(handles_drgb.drgb.lfpevpair(1).PAC(1).all_phase_histo);
all_phase_histo=zeros(no_events,noWindows,no_groups,no_percent_bins,no_PACpeaks,no_lfpevpairs,sz_all_phase(2));


%Sort out all the MI and phase angle values into the matrix
for lfpodNo=1:no_lfpevpairs
    fileNo=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
    groupNo=handles_drgb.drgbchoices.group_no(fileNo);
    timeWindow=handles_drgb.drgb.lfpevpair(lfpodNo).timeWindow;
    for evTypeNo=1:no_events
        trials_in_event=handles_drgb.drgb.lfpevpair(lfpodNo).which_eventPAC(evTypeNo,:)==1;
        if sum(trials_in_event)>0
            for percent_bin=1:no_percent_bins
                trials_in_perbin=(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrPAC>percent_low(percent_bin))&(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrPAC<=percent_high(percent_bin));
                if sum(trials_in_event&trials_in_perbin)>0
                    for PACno=1:handles_drgb.no_PACpeaks
                        
                        no_trials=sum(trials_in_event&trials_in_perbin);
                        no_values(evTypeNo,timeWindow,groupNo,percent_bin,PACno)=no_values(evTypeNo,timeWindow,groupNo,percent_bin,PACno)+1;
                        
                        these_mi_values=zeros(no_trials,1);
                        these_mi_values=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(PACno).mod_indx(trials_in_event&trials_in_perbin)';
                        mi_values(evTypeNo,timeWindow,groupNo,percent_bin,PACno,no_values(evTypeNo,timeWindow,groupNo,percent_bin,PACno))=...
                            mean(these_mi_values);
                        
               
                        these_pa_values=zeros(no_trials,1);
                        these_pa_values=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(PACno).peakAngle(trials_in_event&trials_in_perbin);
                        pa_values(evTypeNo,timeWindow,groupNo,percent_bin,PACno,no_values(evTypeNo,timeWindow,groupNo,percent_bin,PACno))=...
                            circ_rad2ang(circ_mean(circ_ang2rad(these_pa_values')))+180;
                        
                        these_mva_values=zeros(no_trials,1);
                        these_mva_values=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(PACno).meanVectorAngle(trials_in_event&trials_in_perbin);
                        mva_values(evTypeNo,timeWindow,groupNo,percent_bin,PACno,no_values(evTypeNo,timeWindow,groupNo,percent_bin,PACno))=...
                            circ_rad2ang(circ_mean(circ_ang2rad(these_mva_values')))+180;
                         
                        this_all_phase_histo=zeros(no_trials,sz_all_phase(2));
                        this_all_phase_histo=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(PACno).all_phase_histo(trials_in_event&trials_in_perbin,1:sz_all_phase(2));
                        all_phase_histo(evTypeNo,timeWindow,groupNo,percent_bin,PACno,no_values(evTypeNo,timeWindow,groupNo,percent_bin,PACno),1:sz_all_phase(2))=...
                            mean(this_all_phase_histo,1);
                        
                    end
                end
            end
        end
        
    end
    
end


%Now do the figures and statistical analysis for the MI
%These are the colors for the different lines
these_lines{1}='-b';
these_lines{2}='-r';
these_lines{3}='-m';
these_lines{4}='-g';
these_lines{5}='-y';
these_lines{6}='-k';
these_lines{7}='-c';
these_lines{8}='--k';
these_lines{9}='--b';
these_lines{10}='--r';
these_lines{11}='--m';
these_lines{12}='--g';
these_lines{13}='--y';
these_lines{14}='--c';

these_bars{1}='b';
these_bars{2}='r';
these_bars{3}='m';
these_bars{4}='g';
these_bars{5}='y';
these_bars{6}='k';
these_bars{7}='c';
these_bars{8}='k';
% 
% ii_stat_tests=zeros(1,handles_drgb.no_PACpeaks);
% p_vals=[];
% p_vals_aph=[];
% both_normal=[];
% all_normal=ones(1,handles_drgb.no_PACpeaks);
% 
% 
% %First figure out whether the distributions are normal and calculate the p
% %values for each comparison (yikes!)
% for PACno=1:handles_drgb.no_PACpeaks
%     
%     %pFDR is computed separately for each PAC
%     
%     %ii_stat_tests is an index for the pairwise comparisons
%     fprintf(['\nP values for ' handles_drgb.PACnames{PACno} '\n\n'])
%     
%     for evTN1=1:length(evTypeNos)
%         eventType1=evTypeNos(evTN1);
%         for grNo1=1:no_groups
%             for per_bin1=1:no_percent_bins
%                 if no_values(eventType1,winNo,grNo1,per_bin1,PACno)>1
%                     for evTN2=1:length(evTypeNos)
%                         eventType2=evTypeNos(evTN2);
%                         for grNo2=1:no_groups
%                             for per_bin2=1:no_percent_bins
%                                 if no_values(eventType2,winNo,grNo2,per_bin2,PACno)>1
%                                     if(evTN1==1)&(evTN2==1)&(per_bin1==1)&(per_bin2==1)&(PACno==1)&(grNo1==1)&(grNo2==2)
%                                         pfffft=1
%                                     end
%                                     ii_stat_tests(PACno)=ii_stat_tests(PACno)+1;
%                                     evTypeNo1(ii_stat_tests(PACno),PACno)=evTN1;
%                                     groupNo1(ii_stat_tests(PACno),PACno)=grNo1;
%                                     percent_bin1(ii_stat_tests(PACno),PACno)=per_bin1;
%                                     evTypeNo2(ii_stat_tests(PACno),PACno)=evTN2;
%                                     groupNo2(ii_stat_tests(PACno),PACno)=grNo2;
%                                     percent_bin2(ii_stat_tests(PACno),PACno)=per_bin2;
%                                     
%                                     %MI stats
%                                     mi_v1=zeros(1,no_values(eventType1,winNo,grNo1,per_bin1,PACno));
%                                     mi_v1(1,1:no_values(eventType1,winNo,grNo1,per_bin1,PACno))=mi_values(eventType1,winNo,grNo1,per_bin1,PACno,1:no_values(eventType1,winNo,grNo1,per_bin1,PACno));
%                                     
%                                     mi_v2=zeros(1,no_values(eventType2,winNo,grNo2,per_bin2,PACno));
%                                     mi_v2(1,1:no_values(eventType2,winNo,grNo2,per_bin2,PACno))=mi_values(eventType2,winNo,grNo2,per_bin2,PACno,1:no_values(eventType2,winNo,grNo2,per_bin2,PACno));
%                                     
%                                     if (adtest(mi_v1)==0)&(adtest(mi_v2)==0)
%                                         %Both normal distributions
%                                         both_normal(ii_stat_tests(PACno),PACno)=1;
%                                         [h, p_vals(ii_stat_tests(PACno),PACno)]=ttest2(mi_v1,mi_v2);
%                                         fprintf(['\np value MI t test for ' handles_drgb.drgbchoices.group_no_names{grNo1} ' ' evTypeLabels{evTN1} ' per bin %d vs. '...
%                                             handles_drgb.drgbchoices.group_no_names{grNo2} ' ' evTypeLabels{evTN2} ' per bin %d = %E'],per_bin1,per_bin2,p_vals(ii_stat_tests(PACno),PACno))
%                                         
%                                     else
%                                         %At least one distribution not normal
%                                         both_normal(ii_stat_tests(PACno),PACno)=0;
%                                         all_normal(1,PACno)=0;
%                                         [p_vals(ii_stat_tests(PACno),PACno),h]=ranksum(mi_v1,mi_v2);
%                                         
%                                         fprintf(['\np value MI MWT for ' handles_drgb.drgbchoices.group_no_names{grNo1} ' ' evTypeLabels{evTN1} ' per bin %d vs. '...
%                                             handles_drgb.drgbchoices.group_no_names{grNo2} ' ' evTypeLabels{evTN2} ' per bin %d = %E'],per_bin1,per_bin2,p_vals(ii_stat_tests(PACno),PACno))
%                                         
%                                     end
%                                     
%                                     %all_phase_histo stats
%                                     aph_v1=zeros(no_values(eventType1,winNo,grNo1,per_bin1,PACno),sz_all_phase(2));
%                                     aph_v1(1:no_values(eventType1,winNo,grNo1,per_bin1,PACno),1:sz_all_phase(2))=all_phase_histo(eventType1,winNo,grNo1,per_bin1,PACno,1:no_values(eventType1,winNo,grNo1,per_bin1,PACno),1:sz_all_phase(2));
%                                     
%                                     
%                                     aph_v2=zeros(no_values(eventType2,winNo,grNo2,per_bin2,PACno),sz_all_phase(2));
%                                     aph_v2(1:no_values(eventType2,winNo,grNo2,per_bin2,PACno),1:sz_all_phase(2))=all_phase_histo(eventType2,winNo,grNo2,per_bin2,PACno,1:no_values(eventType2,winNo,grNo2,per_bin2,PACno),:);
%                                     
%                                     %Setup a vector for the ANOVAs
%                                     aph_v=zeros(1,(no_values(eventType1,winNo,grNo1,per_bin1,PACno)+no_values(eventType2,winNo,grNo2,per_bin2,PACno))*sz_all_phase(2));
%                                     group=[];
%                                     angle=[];
%                                     no_entries=0;
%                                     for ii=1:no_values(eventType1,winNo,grNo1,per_bin1,PACno)
%                                         aph_v(no_entries+1:no_entries+sz_all_phase(2))=all_phase_histo(eventType1,winNo,grNo1,per_bin1,PACno,ii,1:sz_all_phase(2));
%                                         group(no_entries+1:no_entries+sz_all_phase(2))=zeros(1,sz_all_phase(2));
%                                         angle(no_entries+1:no_entries+sz_all_phase(2))=[1:sz_all_phase(2)];
%                                         no_entries=no_entries+sz_all_phase(2);
%                                     end
%                                     
%                                     for ii=1:no_values(eventType2,winNo,grNo2,per_bin2,PACno)
%                                         aph_v(no_entries+1:no_entries+sz_all_phase(2))=all_phase_histo(eventType2,winNo,grNo2,per_bin2,PACno,ii,:);
%                                         group(no_entries+1:no_entries+sz_all_phase(2))=ones(1,sz_all_phase(2));
%                                         angle(no_entries+1:no_entries+sz_all_phase(2))=[1:sz_all_phase(2)];
%                                         no_entries=no_entries+sz_all_phase(2);
%                                     end
%                                     
%                                     %I cannot find a nonparametric test. I
%                                     %defult to anovan
%                                     %                                     if (adtest(aph_v1(:))==0)&(adtest(aph_v2(:))==0)
%                                     %Both normal distributions
%                                     both_normal(ii_stat_tests(PACno),PACno)=1;
%                                     
%                                     p = anovan(aph_v,{group angle},'display','off');
%                                     p_vals_aph(ii_stat_tests(PACno),PACno)=p(1);
%                                     
%                                     fprintf(['\np value amplitude phase histogram N-ANOVA for ' handles_drgb.drgbchoices.group_no_names{grNo1} ' ' evTypeLabels{evTN1} ' per bin %d vs. '...
%                                         handles_drgb.drgbchoices.group_no_names{grNo2} ' ' evTypeLabels{evTN2} ' per bin %d = %E'],per_bin1,per_bin2,p_vals(ii_stat_tests(PACno),PACno))
%                                     
%                                     %                                     else
%                                     %                                         %At least one distribution not normal
%                                     %                                         both_normal(ii_stat_tests(PACno),PACno)=0;
%                                     %                                         all_normal(1,PACno)=0;
%                                     %                                         [p_vals(ii_stat_tests(PACno),PACno),h]=ranksum(mi_v1,mi_v2);
%                                     %
%                                     %                                         fprintf(['\np value MI MWT for ' handles_drgb.drgbchoices.group_no_names{grNo1} ' ' evTypeLabels{evTN1} ' per bin %d vs. '...
%                                     %                                             handles_drgb.drgbchoices.group_no_names{grNo2} ' ' evTypeLabels{evTN2} ' per bin %d = %E'],per_bin1,per_bin2,p_vals(ii_stat_tests(PACno),PACno))
%                                     %
%                                     %                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     
%     these_p_vals=zeros(1,ii_stat_tests(PACno));
%     these_p_vals(1,:)=p_vals(1:ii_stat_tests(PACno),PACno);
%     
%     pFDR(PACno)=drsFDRpval(these_p_vals);
%     fprintf('\n pFDR for MI = %E\n',pFDR(PACno))
%     
%     these_p_vals=zeros(1,ii_stat_tests(PACno));
%     these_p_vals(1,:)=p_vals_aph(1:ii_stat_tests(PACno),PACno);
%     pFDR_aph(PACno)=drsFDRpval(these_p_vals);
%     fprintf('\n pFDR for MI = %E\n',pFDR_aph(PACno))
%     
%     fprintf('\n\n\n')
% end
% %mi_values(evTypeNo,timeWindow,groupNo,percent_bin,PACno,:)
% 
% 
% %Now do the analysis for peak phase angle
% 
% %First, determine whetehr the distributions are uniform using a Raileigh
% %test and determine whether the peak angles differ
% ii_stat_ral_tests=zeros(1,handles_drgb.no_PACpeaks);
% p_vals_ral=[];
% 
% for PACno=1:handles_drgb.no_PACpeaks
%     
%     %pFDR is computed separately for each PAC
%     
%     %ii_stat_tests is an index for the pairwise comparisons
%     
%     
%     for evTN1=1:length(evTypeNos)
%         eventType1=evTypeNos(evTN1);
%         for grNo1=1:no_groups
%             for per_bin1=1:no_percent_bins
%                 if no_values(eventType1,winNo,grNo1,per_bin1,PACno)>1
%                     
%                     ii_stat_ral_tests(PACno)=ii_stat_ral_tests(PACno)+1;
%                     evTypeNo_ral(ii_stat_ral_tests(PACno),PACno)=evTN1;
%                     groupNo_ral(ii_stat_tests(PACno),PACno)=grNo1;
%                     percent_bin_ral(ii_stat_tests(PACno),PACno)=per_bin1;
%                     
%                     
%                     pa_v1=zeros(1,no_values(eventType1,winNo,grNo1,per_bin1,PACno));
%                     pa_v1(1,1:no_values(eventType1,winNo,grNo1,per_bin1,PACno))=pa_values(eventType1,winNo,grNo1,per_bin1,PACno,1:no_values(eventType1,winNo,grNo1,per_bin1,PACno));
%                     
%                     % Rayleigh test
%                     p_vals_ral(ii_stat_ral_tests(PACno),PACno) = circ_rtest(circ_ang2rad(pa_v1));
%                     
%                     
%                 end
%             end
%         end
%         
%     end
%     these_p_vals=zeros(1,ii_stat_ral_tests(PACno));
%     these_p_vals(1,:)=p_vals(1:ii_stat_ral_tests(PACno),PACno);
%     pFDR_ral(PACno)=drsFDRpval(these_p_vals);
% end
% 
% 
% %Calculate the p values for mean and median angle for each comparison
% 
% ii_stat_angle_tests=zeros(1,handles_drgb.no_PACpeaks);
% p_vals_mean_ww=[];
% p_vals_median_fisher=[];
% 
% 
% for PACno=1:handles_drgb.no_PACpeaks
%     
%     %pFDR is computed separately for each PAC
%     
%     %ii_stat_tests is an index for the pairwise comparisons
%     
%     fprintf(['\nP values for ' handles_drgb.PACnames{PACno} '\n\n'])
%     
%     for evTN1=1:length(evTypeNos)
%         eventType1=evTypeNos(evTN1);
%         for grNo1=1:no_groups
%             for per_bin1=1:no_percent_bins
%                 if no_values(eventType1,winNo,grNo1,per_bin1,PACno)>1
%                     for evTN2=1:length(evTypeNos)
%                         eventType2=evTypeNos(evTN2);
%                         for grNo2=1:no_groups
%                             for per_bin2=1:no_percent_bins
%                                 if no_values(eventType2,winNo,grNo2,per_bin2,PACno)>1
%                                     ii_stat_angle_tests(PACno)=ii_stat_angle_tests(PACno)+1;
%                                     evTypeNo1_angle(ii_stat_angle_tests(PACno),PACno)=evTN1;
%                                     groupNo1_angle(ii_stat_angle_tests(PACno),PACno)=grNo1;
%                                     percent_bin1_angle(ii_stat_angle_tests(PACno),PACno)=per_bin1;
%                                     evTypeNo2_angle(ii_stat_angle_tests(PACno),PACno)=evTN2;
%                                     groupNo2_angle(ii_stat_angle_tests(PACno),PACno)=grNo2;
%                                     percent_bin2_angle(ii_stat_angle_tests(PACno),PACno)=per_bin2;
%                                     
%                                     pa_v1=zeros(1,no_values(eventType1,winNo,grNo1,per_bin1,PACno));
%                                     pa_v1(1,1:no_values(eventType1,winNo,grNo1,per_bin1,PACno))=pa_values(eventType1,winNo,grNo1,per_bin1,PACno,1:no_values(eventType1,winNo,grNo1,per_bin1,PACno));
%                                     
%                                     pa_v2=zeros(1,no_values(eventType2,winNo,grNo2,per_bin2,PACno));
%                                     pa_v2(1,1:no_values(eventType2,winNo,grNo2,per_bin2,PACno))=pa_values(eventType2,winNo,grNo2,per_bin2,PACno,1:no_values(eventType2,winNo,grNo2,per_bin2,PACno));
%                                     
%                                     try
%                                         p_vals_mean_ww(ii_stat_angle_tests(PACno),PACno)=circ_wwtest(pa_v1, pa_v2);
%                                         p_vals_median_fisher(ii_stat_angle_tests(PACno),PACno)=circ_cmtest(pa_v1, pa_v2);
%                                     catch
%                                         p_vals_mean_ww(ii_stat_angle_tests(PACno),PACno)=1;
%                                         p_vals_median_fisher(ii_stat_angle_tests(PACno),PACno)=1;
%                                     end
%                                     
%                                     
%                                     fprintf(['\np value peak angle welch for ' handles_drgb.drgbchoices.group_no_names{grNo1} ' ' evTypeLabels{evTN1} ' per bin %d vs. '...
%                                         handles_drgb.drgbchoices.group_no_names{grNo2} ' ' evTypeLabels{evTN2} ' per bin %d = %E'],per_bin1,per_bin2,p_vals_mean_ww(ii_stat_angle_tests(PACno),PACno))
%                                     
%                                     fprintf(['\np value peak angle fisher for ' handles_drgb.drgbchoices.group_no_names{grNo1} ' ' evTypeLabels{evTN1} ' per bin %d vs. '...
%                                         handles_drgb.drgbchoices.group_no_names{grNo2} ' ' evTypeLabels{evTN2} ' per bin %d = %E'],per_bin1,per_bin2,p_vals_median_fisher(ii_stat_angle_tests(PACno),PACno))
%                                     
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     
%     these_p_vals=zeros(1,ii_stat_angle_tests(PACno));
%     these_p_vals(1,:)=p_vals_mean_ww(1:ii_stat_angle_tests(PACno),PACno);
%     pFDR_angle_ww(PACno)=drsFDRpval(these_p_vals);
%     fprintf('\n pFDR for peak angle welch = %E\n',pFDR_angle_ww(PACno))
%     
%     these_p_vals=zeros(1,ii_stat_angle_tests(PACno));
%     these_p_vals(1,:)=p_vals_median_fisher(1:ii_stat_angle_tests(PACno),PACno);
%     pFDR_angle_fisher(PACno)=drsFDRpval(these_p_vals);
%     fprintf('\n pFDR for peak angle fisher = %E\n',pFDR_angle_fisher(PACno))
%     
%     fprintf('\n\n\n')
% end
% %mi_values(evTypeNo,timeWindow,groupNo,percent_bin,PACno,:)

%Display the results for each PAC
%Sort by percent bin (rows) and within each row evTypeNo(groupNo)
legend_bar=[];
for grNo=1:no_groups
    legend_bar=[legend_bar handles_drgb.drgbchoices.group_no_names{grNo} ' '];
end


%Plot the phase angle
switch which_PA_display
    
    case 1
        %Cumulative histogram for the peak angles
        figNo=0;
        
        for PACno=1:handles_drgb.no_PACpeaks
            for per_bin=1:no_percent_bins
                for grNo=1:no_groups
                    figNo=figNo+1;
                    try
                        close(figNo)
                    catch
                    end
                    hFig=figure(figNo);
                    set(hFig, 'units','normalized','position',[.05 .05 .7 .7])
                    
                    this_legend=[];
                    hold on
                    n_cum=0;
                    
                    for evTN=1:length(evTypeNos)
                        eventType=evTypeNos(evTN);
                        
                        if no_values(eventType,winNo,grNo,per_bin,PACno)>1
                            pa_v=zeros(1,no_values(eventType,winNo,grNo,per_bin,PACno));
                            pa_v(1,1:no_values(eventType,winNo,grNo,per_bin,PACno))=...
                                pa_values(eventType,winNo,grNo,per_bin,PACno,1:no_values(eventType,winNo,grNo,per_bin,PACno));
                            [f_pa_v,x_pa_v] = drg_ecdf(pa_v);
                            n_cum=n_cum+1;
                            plot(x_pa_v,f_pa_v,these_lines{n_cum})
                            this_legend=[this_legend '''' evTypeLabels{evTN} '''' ','];
                        end
                    end
                    this_legend=['legend(' this_legend(1:end-1) ')'];
                    eval(this_legend)
                    title(['Peak phase angle for ' handles_drgb.PACnames{PACno} '. Percent correct from ' num2str(percent_low(per_bin))...
                        ' to ' num2str(percent_high(per_bin)) ', group No: ' handles_drgb.drgbchoices.group_no_names{grNo}])
                    
                    %Do the KS test
                    evTN=ev_to_test_against;
                    eventType=evTypeNos(evTN);
                    p_vals=[];
                    if no_values(eventType,winNo,grNo,per_bin,PACno)>1
                        pa_ref=zeros(1,no_values(eventType,winNo,grNo,per_bin,PACno));
                        pa_ref(1,1:no_values(eventType,winNo,grNo,per_bin,PACno))=...
                            pa_values(eventType,winNo,grNo,per_bin,PACno,1:no_values(eventType,winNo,grNo,per_bin,PACno));
                        fprintf(1,['KS p values for peak phase angle cumulative probability for ' handles_drgb.PACnames{PACno} '. Percent correct from ' num2str(percent_low(per_bin))...
                            ' to ' num2str(percent_high(per_bin)) ', group No: ' handles_drgb.drgbchoices.group_no_names{grNo} '\n'])
                        
                        for evTN=1:length(evTypeNos)
                            if evTN~=ev_to_test_against
                                eventType=evTypeNos(evTN);
                                if no_values(eventType,winNo,grNo,per_bin,PACno)>1
                                    pa_v=zeros(1,no_values(eventType,winNo,grNo,per_bin,PACno));
                                    pa_v(1,1:no_values(eventType,winNo,grNo,per_bin,PACno))=...
                                        pa_values(eventType,winNo,grNo,per_bin,PACno,1:no_values(eventType,winNo,grNo,per_bin,PACno));
                                    [h,p]=kstest2(pa_v,pa_ref);
                                    p_vals=[p_vals p];
                                    fprintf(1, [evTypeLabels{evTN} ' vs. ' evTypeLabels{ev_to_test_against} '= %d\n'],p);
                                end
                            end
                        end
                        
                        if length(p_vals)>1
                            fprintf(1, 'pFDR= = %d\n',drsFDRpval(p_vals));
                        end
                        fprintf(1, '\n');
                    end
                    
                    if figNo==17
                        pfft=1;
                    end
                    
                end
            end
        end
        
    case 2
        %Now do the rose plots for peak angles
        
        for PACno=1:handles_drgb.no_PACpeaks
            try
                close(PACno+2*handles_drgb.no_PACpeaks)
            catch
            end
            hFig=figure(PACno+2*handles_drgb.no_PACpeaks);
            set(hFig, 'units','normalized','position',[.05 .05 .7 .7])
            
            no_sub=0;
            
            for per_bin=1:no_percent_bins
                this_legend=[];
                no_bars=0;
                first_sub=1;
                
                XTicks=[];
                
                for evTN=1:length(evTypeNos)
                    eventType=evTypeNos(evTN);
                    for grNo=1:no_groups
                        no_sub=no_sub+1;
                        
                        if no_values(eventType,winNo,grNo,per_bin,PACno)>1
                            subplot(no_percent_bins,length(evTypeNos)*no_groups,no_sub)
                            
                            pa_v=zeros(1,no_values(eventType,winNo,grNo,per_bin,PACno));
                            pa_v(1,1:no_values(eventType,winNo,grNo,per_bin,PACno))=pa_values(eventType,winNo,grNo,per_bin,PACno,1:no_values(eventType,winNo,grNo,per_bin,PACno));
                            
                            rose(pi*pa_v/180,12)
                            xlabel([evTypeLabels{evTN} ' ' handles_drgb.drgbchoices.group_no_names{grNo}])
                            
                            if first_sub==1
                                title(['Phase angle for ' handles_drgb.PACnames{PACno} '. Percent correct from ' num2str(percent_low(per_bin)) ' to ' num2str(percent_high(per_bin))])
                                first_sub=0;
                            end
                        end
                    end
                    
                end
                
                
            end
            
            
        end
        
    case 3
        
        %Now plot the PAC probability histogram
        no_bins=sz_all_phase(2)-1;
        phase=[0:360/no_bins:360];
        for PACno=1:handles_drgb.no_PACpeaks
            try
                close(PACno+3*handles_drgb.no_PACpeaks)
            catch
            end
            hFig=figure(PACno+3*handles_drgb.no_PACpeaks);
            set(hFig, 'units','normalized','position',[.05 .05 .7 .7])
            
            no_sub=0;
            
            
            for per_bin=1:no_percent_bins
                this_legend=[];
                no_bars=0;
                first_sub=1;
                no_sub=0;
                %legend_string=['legend('];
                
                legend_string_right=['],{'];
                XTicks=[];
                subplot(no_percent_bins,1,per_bin)
                for evTN=1:length(evTypeNos)
                    eventType=evTypeNos(evTN);
                    %for grNo=1:no_groups
                    for grNo=1:1
                        
                        
                        if no_values(eventType,winNo,grNo,per_bin,PACno)>1
                            
                            no_sub=no_sub+1;
                            
                            phase_histo=zeros(no_values(eventType,winNo,grNo,per_bin,PACno),sz_all_phase(2));
                            phase_histo(:,:)=all_phase_histo(eventType,winNo,grNo,per_bin,PACno,1:no_values(eventType,winNo,grNo,per_bin,PACno),1:sz_all_phase(2));
                            pct5=prctile(phase_histo,5);
                            pct95=prctile(phase_histo,95);
                            pctiles=[pct95;pct5];
                            shadedErrorBar(phase,mean(phase_histo,1),pctiles,these_lines{no_sub})
                            hold on
                            
                            if (PACno==3)&(per_bin==3)&(evTN==1)&(grNo==1)
                                figure(14)
                                for electNo=1:16
                                    plot(phase,phase_histo(electNo,:))
                                    hold on
                                end
                                figure(PACno+3*handles_drgb.no_PACpeaks)
                            end
                            
                            if (PACno==3)&(per_bin==3)&(evTN==2)&(grNo==1)
                                figure(15)
                                for electNo=1:16
                                    plot(phase,phase_histo(electNo,:))
                                    hold on
                                end
                                figure(PACno+3*handles_drgb.no_PACpeaks)
                            end
                            
                            legend_string_right=[legend_string_right [''''] handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evTN} [''''] [',']];
                        end
                    end
                    
                end
                
                %I have no idea why this works backwards, but it does!
                legend_string_left=['legend(['];
                for plotNo=no_sub:-1:1
                    legend_string_left=[legend_string_left ' h(' num2str(plotNo*3-2) ')'];
                end
                h=findobj(gca,'Type','Line');
                legend_string=[legend_string_left legend_string_right(1:end-1) ['})']];
                eval(legend_string)
                
                
                title(['PAC for ' handles_drgb.PACnames{PACno} '. Percent correct from ' num2str(percent_low(per_bin)) ' to ' num2str(percent_high(per_bin))])
                xlim([0 360])
                
                ylabel('Probability')
                xlabel('Phase (deg)')
                
            end
            
            
        end
    case 4
        %Display the results for each PAC
        %Sort by percent bin (rows) and within each row evTypeNo(groupNo)
        legend_bar=[];
        for grNo=1:no_groups
            legend_bar=[legend_bar handles_drgb.drgbchoices.group_no_names{grNo} ' '];
        end
        
        
        for PACno=1:handles_drgb.no_PACpeaks
            %Do the analysis for MI
            try
                close(PACno)
            catch
            end
            hFig=figure(PACno);
            set(hFig, 'units','normalized','position',[.05 .05 .7 .7])
            
            
            %Non-parametric
            %Normal distributions
            for per_bin=1:no_percent_bins
                MIs=[];
                MIgroups=[];
                ii_groups=0;
                
                subplot(no_percent_bins,1,per_bin)
                
                
                XTicks=[];
                
                for evTN=1:length(evTypeNos)
                    eventType=evTypeNos(evTN);
                    for grNo=1:no_groups
                        
                        if no_values(evTN,winNo,grNo,per_bin,PACno)>1
                            
                            mi_v=zeros(1,no_values(eventType,winNo,grNo,per_bin,PACno));
                            mi_v(1,1:no_values(eventType,winNo,grNo,per_bin,PACno))=mi_values(eventType,winNo,grNo,per_bin,PACno,1:no_values(eventType,winNo,grNo,per_bin,PACno));
                            
                            
                            MIs=[MIs mi_v];
                            
                            this_gr=[handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evTN} ' ' num2str(length(mi_v))];
                            for ii=1:no_values(eventType,winNo,grNo,per_bin,PACno)
                                ii_groups=ii_groups+1;
                                MIgroups{ii_groups}=this_gr;
                            end
                            
                        else
                            %                         mi_v=zeros(1,5);
                            %                         MIs=[MIs mi_v];
                            %                         for ii=1:5
                            %                             ii_groups=ii_groups+1;
                            %                             MIgroups{ii_groups}=[' '];
                            %                         end
                            %
                        end
                    end
                    
                end
                
                boxplot(MIs,MIgroups,'symbol','k')
                
                title(['MI for ' handles_drgb.PACnames{PACno} '. Percent correct from ' num2str(percent_low(per_bin)) ' to ' num2str(percent_high(per_bin))])
                
                
            end
            
            
        end
        
        
    case 5
        %Display the cumulative probabilities
        for PACno=1:handles_drgb.no_PACpeaks
            %Do the analysis for MI
            try
                close(PACno+handles_drgb.no_PACpeaks)
            catch
            end
            hFig=figure(PACno+handles_drgb.no_PACpeaks);
            set(hFig, 'units','normalized','position',[.05 .05 .7 .7])
            hold on
            
            per_bin=no_percent_bins;
            n_cum=0;
            legend_string=['legend('];
            for evTN=1:length(evTypeNos)
                eventType=evTypeNos(evTN);
                for grNo=1:no_groups
                    
                    if no_values(eventType,winNo,grNo,per_bin,PACno)>1
                        n_cum=n_cum+1;
                        barNo=(no_groups+1)*(evTN-1)+evTN;
                        mi_v=zeros(1,no_values(eventType,winNo,grNo,per_bin,PACno));
                        mi_v(1,1:no_values(eventType,winNo,grNo,per_bin,PACno))=mi_values(eventType,winNo,grNo,per_bin,PACno,1:no_values(eventType,winNo,grNo,per_bin,PACno));
                        
                        [f_mi_v,x_mi_v] = drg_ecdf(mi_v);
                        plot(x_mi_v,f_mi_v,these_lines{n_cum})
                        handles_drgb.pac_vars(PACno,evTN,per_bin,grNo).mi_v=mi_v;
                        
                    end
                    %             if (grNo==no_groups)&(evTN==length(evTypeNos))
                    %                 legend_string=[legend_string [''''] handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evTN} [''''] [')']];
                    %             else
                    legend_string=[legend_string [''''] handles_drgb.drgbchoices.group_no_names{grNo} ' ' evTypeLabels{evTN} [''''] [',']];
                    
                    %             end
                end
                
            end
            
            legend_string=[legend_string(1:end-1) [')']];
            eval(legend_string)
            title(['MI for ' handles_drgb.PACnames{PACno} '. Percent correct from ' num2str(percent_low(per_bin)) ' to ' num2str(percent_high(per_bin))])
            
        end
    case 6
        %Cumulative histogram for the mean phase angle
        figNo=0;
         
        for PACno=1:handles_drgb.no_PACpeaks
            for per_bin=1:no_percent_bins
                for grNo=1:no_groups
                    figNo=figNo+1;
                    try
                        close(figNo)
                    catch
                    end
                    hFig=figure(figNo);
                    set(hFig, 'units','normalized','position',[.05 .05 .7 .7])
                    
                    this_legend=[];
                    hold on
                    n_cum=0;
                    
                    for evTN=1:length(evTypeNos)
                        eventType=evTypeNos(evTN);
                        
                        if no_values(eventType,winNo,grNo,per_bin,PACno)>1
                            mva_v=zeros(1,no_values(eventType,winNo,grNo,per_bin,PACno));
                            mva_v(1,1:no_values(eventType,winNo,grNo,per_bin,PACno))=...
                                mva_values(eventType,winNo,grNo,per_bin,PACno,1:no_values(eventType,winNo,grNo,per_bin,PACno));
                            [f_mva_v,x_mva_v] = drg_ecdf(mva_v);
                            n_cum=n_cum+1;
                            plot(x_mva_v,f_mva_v,these_lines{n_cum})
                            this_legend=[this_legend '''' evTypeLabels{evTN} '''' ','];
                        end
                    end
                    this_legend=['legend(' this_legend(1:end-1) ')'];
                    eval(this_legend)
                    title(['Mean phase angle for ' handles_drgb.PACnames{PACno} '. Percent correct from ' num2str(percent_low(per_bin))...
                        ' to ' num2str(percent_high(per_bin)) ', group No: ' handles_drgb.drgbchoices.group_no_names{grNo}])
                    
                    %Do the KS test
                    evTN=ev_to_test_against;
                    eventType=evTypeNos(evTN);
                    p_vals=[];
                    if no_values(eventType,winNo,grNo,per_bin,PACno)>1
                        mva_ref=zeros(1,no_values(eventType,winNo,grNo,per_bin,PACno));
                        mva_ref(1,1:no_values(eventType,winNo,grNo,per_bin,PACno))=...
                            mva_values(eventType,winNo,grNo,per_bin,PACno,1:no_values(eventType,winNo,grNo,per_bin,PACno));
                        fprintf(1,['KS p values for mean phase angle cumulative probability for ' handles_drgb.PACnames{PACno} '. Percent correct from ' num2str(percent_low(per_bin))...
                            ' to ' num2str(percent_high(per_bin)) ', group No: ' handles_drgb.drgbchoices.group_no_names{grNo} '\n'])
                        
                        for evTN=1:length(evTypeNos)
                            if evTN~=ev_to_test_against
                                eventType=evTypeNos(evTN);
                                if no_values(eventType,winNo,grNo,per_bin,PACno)>1
                                    mva_v=zeros(1,no_values(eventType,winNo,grNo,per_bin,PACno));
                                    mva_v(1,1:no_values(eventType,winNo,grNo,per_bin,PACno))=...
                                        mva_values(eventType,winNo,grNo,per_bin,PACno,1:no_values(eventType,winNo,grNo,per_bin,PACno));
                                    [h,p]=kstest2(mva_v,mva_ref);
                                    p_vals=[p_vals p];
                                    fprintf(1, [evTypeLabels{evTN} ' vs. ' evTypeLabels{ev_to_test_against} '= %d\n'],p);
                                end
                            end
                        end
                        
                        if length(p_vals)>1
                            fprintf(1, 'pFDR= = %d\n',drsFDRpval(p_vals));
                        end
                        fprintf(1, '\n');
                    end
                    
                end
            end
        end
        
end

save([handles.PathName handles.drgb.outFileName(1:end-4) '_out.mat'])

pffft=1
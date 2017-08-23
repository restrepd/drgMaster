function drgDisplayBatchPerCorrPairwise(handles)

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
% 1 Show the difference of dB vs frequency between events for each group/ percent correct
%  bin in a separate graph
%
% 2 Show the difference of dB vs frequency between groups for each event/percent correct
%  bin in a separate graph
%
% 3 Show the difference of dB vs frequency between learning and proficient
%
% 4 Show boxplots for the difference of dB per bandwidth between learning and proficient
%
% 5 Show boxplots for the difference in dB per badwidht between groups


%For Alexia's tstart learning vs. proficeint
% winNo=3;
% which_display=3;
% eventType=1; %tstart
% evTypeLabels={'tstart'};

% For Alexia's odorOn (CS) learning vs. proficeint
% winNo=2;
% which_display=4;
% eventType=1; %OdorOn
% evTypeLabels={'CS'};


% For Alexia's Hit, CR, FA
% winNo=2;
% which_display=1;
% eventType=[2 5 7];
% evTypeLabels={'Hit','CR','FA'};

% % For Daniel's Hit, CR, compare groups
winNo=2;
% which_display=2;
which_display=5;
eventType=[2 5];
evTypeLabels={'Hit','CR'};

% For Daniel's Hit, CR, FA, compare events
% winNo=2;
% which_display=1;
% eventType=[2 5];
% evTypeLabels={'Hit','CR'};

%THESE VALUES ARE IMPORTANT
%VERY IMPORTANT: This is the index for this event in your handles.drgbchoices.evTypeNos
% eventType=[2 5]; %Hit and CR
% evTypeLabels={'Hit';'CR'};






% eventType=[9 10 11 12 13 14]; %Hit and CR
% no_event_types=6;
% evTypeLabels={'Hi Od1';'Hi Od 2';'Hi Od 3';'Low Od1';'Low Od 2';'Low Od 3'};

no_event_types=length(eventType);

%Statistics
%stats_method = 1, FDR, 2, statscond
stat_method=2;

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

%Pairwise comparisons of percent correct for no laser vs laser
%acetophenone vs ethyl benzoate for Daniel
file_pairs=[1 7;
    2 8;
    3 9
    4 10;
    5 11;
    6 12;
    13 14;
    15 16;
    17 21;
    18 23;
    19 22;
    20 24];
no_file_pairs=12;
    

%Ask user for the drgb output .mat file and load those data
[handles.drgb.outFileName,handles.PathName] = uigetfile('*.mat','Select the drgb output file');
load([handles.PathName handles.drgb.outFileName])


%Initialize the variables

%Determine how many LFPs are included for each evTypeNo, time window, group, etc
no_groups=max(handles_drgb.drgbchoices.group_no);
no_events=length(handles_drgb.drgbchoices.evTypeNos);
no_lfpevpairs=handles_drgb.drgb.lfpevpair_no;
no_lfps=zeros(length(handles_drgb.drgbchoices.evTypeNos),handles_drgb.drgbchoices.noWindows,max(handles_drgb.drgbchoices.group_no),no_percent_bins);

files=[];
noFiles=0;

for lfpodNo=1:no_lfpevpairs
    fileNo=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
    noFiles=noFiles+1;
    files(noFiles)=fileNo;
    groupNo=handles_drgb.drgbchoices.group_no(fileNo);
    group_per_file(fileNo)=groupNo;
    timeWindow=handles_drgb.drgb.lfpevpair(lfpodNo).timeWindow;
    for evTypeNo=1:no_events
        
        try
            trials_in_event=handles_drgb.drgb.lfpevpair(lfpodNo).which_eventLFPPower(evTypeNo,:)==1;
        catch
            trials_in_event=[];
        end
        if sum(trials_in_event)>0
            for percent_bin=1:no_percent_bins
                trials_in_perbin=(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrLFPPower>percent_low(percent_bin))&(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrLFPPower<=percent_high(percent_bin));
                
                if sum(trials_in_event&trials_in_perbin)>0
                    no_lfps(evTypeNo,timeWindow,groupNo,percent_bin)=no_lfps(evTypeNo,timeWindow,groupNo,percent_bin)+1;
                end
            end
        end
    end
end
max_lfps=max(no_lfps(:));


fprintf(1, '\nTotal number of files: = %d\n\n',length(group_per_file));
for grNo=1:no_groups
    fprintf(1, ['Number of files for ' handles_drgb.drgbchoices.group_no_names{grNo} ' = %d\n'],sum(group_per_file==grNo));
end
fprintf(1, '\n\n');

%Now initialize the perCorr variable
perCorr_per_file=[];
no_trials=zeros(1,length(handles_drgb.drgb.file));
file_no=0;

%Sort out all perCorr per file into a vector
for lfpodNo=1:no_lfpevpairs
    fileNo=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
    if fileNo>file_no
        file_no=fileNo;
        perCorr_per_file(fileNo).perCorr=handles_drgb.drgb.lfpevpair(lfpodNo).perCorrLFPPower;
        perCorr_per_file(fileNo).groupNo=handles_drgb.drgbchoices.group_no(fileNo);
    end  
end


%Now do a pairwise comparison of the last and first percent correct values

%Do NL vs L
last_pc_NL=[];
first_pc_L=[];
last_pc_NLc=[];
first_pc_Lc=[];
for fps=1:no_file_pairs
    if (perCorr_per_file(file_pairs(fps,1)).groupNo==1)||(perCorr_per_file(file_pairs(fps,1)).groupNo==2)
        for noFil=1:2
            if perCorr_per_file(file_pairs(fps,noFil)).groupNo==1
                last_pc_NL=[last_pc_NL perCorr_per_file(file_pairs(fps,noFil)).perCorr(end)];
            else
                first_pc_L=[first_pc_L perCorr_per_file(file_pairs(fps,noFil)).perCorr(end)];
            end
            
            
        end
        
    end
    
    if (perCorr_per_file(file_pairs(fps,1)).groupNo==3)||(perCorr_per_file(file_pairs(fps,1)).groupNo==4)
        for noFil=1:2
            
            
            if perCorr_per_file(file_pairs(fps,noFil)).groupNo==3
                last_pc_NLc=[last_pc_NLc perCorr_per_file(file_pairs(fps,noFil)).perCorr(end)];
            else
                first_pc_Lc=[first_pc_Lc perCorr_per_file(file_pairs(fps,noFil)).perCorr(end)];
            end
            
        end
        
    end
end

[h p_NL_L]=ttest(last_pc_NL,first_pc_L)

[h p_NLc_Lc]=ttest(last_pc_NLc,first_pc_Lc)


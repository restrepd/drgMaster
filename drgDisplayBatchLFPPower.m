function drgDisplayBatchLFPPower(handles)

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
% 5 Show boxplots for the difference in dB per badwidth between groups
%
% 6 Analysis of delta power (Hit-CR)
%   Show boxplots for the difference in dB per badwidht between groups
%
% 7 Do between group auROC analysis for Hit vs CR for proficient
%
% 8 d prime analysis for FA
%
% 9 d prime analysis for Fig. 3 of Daniel's pub
% use acetoethylben_electrode9202017.mat


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
% winNo=2;
% % which_display=2;
% which_display=5;
% eventType=[2 5];
% evTypeLabels={'Hit','CR'};
% % eventType=[2 5 7];
% % evTypeLabels={'Hit','CR','FA'};

% For Daniel's Hit, CR, FA, compare events
% winNo=2;
% which_display=1;
% eventType=[2 5];
% evTypeLabels={'Hit','CR'};

% % For Justin's Hit, CR, compare groups
% winNo=2;
% % which_display=2;
% which_display=7;
% eventType=[2 5];
% evTypeLabels={'Hit','CR'};
% % eventType=[2 5 7];
% % evTypeLabels={'Hit','CR','FA'};

%THESE VALUES ARE IMPORTANT
%VERY IMPORTANT: This is the index for this event in your handles.drgbchoices.evTypeNos
% eventType=[2 5]; %Hit and CR
% evTypeLabels={'Hit';'CR'};


% % For Daniel's delta power Hit - CR, compare groups
% which_display = 9 is used for the d prime calculations for Fig. 3
% run with:
% isomin_firstandlastIAMO91817.mat
% acetoethylben_firstandlast91117.mat
% ethylacetatepropylacetatefirstandlast92617.mat
%
winNo=2;
which_display=9;
eventType1=2;
eventType=2;
eventTypeRef=5;
evTypeLabel='Hit';
evTypeRefLabel='CR';

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



%Ask user for the drgb output .mat file and load those data
[handles.drgb.outFileName,handles.PathName] = uigetfile('*.mat','Select the drgb output file');
load([handles.PathName handles.drgb.outFileName])

fprintf(1, ['\ndrgDisplayBatchLFPPower run for ' handles.drgb.outFileName ', which_display= = %d\n\n'],which_display);

evHit=find(handles_drgb.drgbchoices.evTypeNos==3);
%evHit=find(handles_drgb.drgbchoices.evTypeNos==5); %S plus
evCR=find(handles_drgb.drgbchoices.evTypeNos==9);
%evCR=find(handles_drgb.drgbchoices.evTypeNos==11); %S-
evFA=find(handles_drgb.drgbchoices.evTypeNos==13);

frequency=handles_drgb.drgb.freq_for_LFPpower;

%Experiment pairs
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

%Now initialize the variables for db power
noWindows=handles_drgb.drgbchoices.noWindows;
dB_power_values=zeros(no_events,noWindows,no_groups,no_percent_bins,max_lfps,length(handles_drgb.drgb.freq_for_LFPpower));
experimentNo=zeros(no_events,noWindows,no_groups,no_percent_bins,max_lfps);
no_values=zeros(no_events,noWindows,no_groups,no_percent_bins);

dB_power_valuesHitCR=zeros(no_events,noWindows,no_groups,no_percent_bins,max_lfps,length(handles_drgb.drgb.freq_for_LFPpower));
experimentNoHitCR=zeros(no_events,noWindows,no_groups,no_percent_bins,max_lfps);
no_valuesHitCR=zeros(no_events,noWindows,no_groups,no_percent_bins);


%Sort out all power values into a matrix
no_ROCs=0;
p_vals_ROC=[];
for lfpodNo=1:no_lfpevpairs
    fileNo=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
    groupNo=handles_drgb.drgbchoices.group_no(fileNo);
    timeWindow=handles_drgb.drgb.lfpevpair(lfpodNo).timeWindow;
    lfpodorNo=lfpodNo;
    
    try
        trials_in_Hit=handles_drgb.drgb.lfpevpair(lfpodNo).which_eventLFPPower(evHit,:)==1;
    catch
        trials_in_Hit=[];
    end
    
    try
        trials_in_CR=handles_drgb.drgb.lfpevpair(lfpodNo).which_eventLFPPower(evCR,:)==1;
    catch
        trials_in_CR=[];
    end
    
    try
        trials_in_FA=handles_drgb.drgb.lfpevpair(lfpodNo).which_eventLFPPower(evFA,:)==1;
    catch
        trials_in_FA=[];
    end
    
    %Calculate the average dB power
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
                    
                    no_trials=sum(trials_in_event&trials_in_perbin);
                    
                    this_dB_power=zeros(no_trials,length(handles_drgb.drgb.freq_for_LFPpower));
                    this_dB_power(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_event&trials_in_perbin,:));
                    dB_power_values(evTypeNo,timeWindow,groupNo,percent_bin,no_values(evTypeNo,timeWindow,groupNo,percent_bin)+1,:)=...
                        mean(this_dB_power,1);
                    experimentNo(evTypeNo,timeWindow,groupNo,percent_bin,no_values(evTypeNo,timeWindow,groupNo,percent_bin)+1)=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
                    no_values(evTypeNo,timeWindow,groupNo,percent_bin)=no_values(evTypeNo,timeWindow,groupNo,percent_bin)+1;
                    if (sum(trials_in_Hit&trials_in_perbin)>0)&(sum(trials_in_CR&trials_in_perbin)>0)
                        dB_power_valuesHitCR(evTypeNo,timeWindow,groupNo,percent_bin,no_values(evTypeNo,timeWindow,groupNo,percent_bin)+1,:)=...
                            mean(this_dB_power,1);
                        experimentNoHitCR(evTypeNo,timeWindow,groupNo,percent_bin,no_valuesHitCR(evTypeNo,timeWindow,groupNo,percent_bin)+1)=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
                        no_valuesHitCR(evTypeNo,timeWindow,groupNo,percent_bin)=no_valuesHitCR(evTypeNo,timeWindow,groupNo,percent_bin)+1;
                    end
                end
            end
        end
        
    end
    
    %Calculate the ROC for Hit vs CR
    %Do proficient only
    percent_bin=3;
    trials_in_perbin=(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrLFPPower>percent_low(percent_bin))&(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrLFPPower<=percent_high(percent_bin));
    
    if (sum(trials_in_Hit&trials_in_perbin)>5)&(sum(trials_in_CR&trials_in_perbin)>5)&((which_display==7)||(which_display==8)||(which_display==9))
        
        no_trialsHit=sum(trials_in_Hit&trials_in_perbin);
        this_dB_powerHit=zeros(no_trialsHit,length(handles_drgb.drgb.freq_for_LFPpower));
        this_dB_powerHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_Hit&trials_in_perbin,:));
        
        
        no_trialsCR=sum(trials_in_CR&trials_in_perbin);
        this_dB_powerCR=zeros(no_trialsCR,length(handles_drgb.drgb.freq_for_LFPpower));
        this_dB_powerCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_CR&trials_in_perbin,:));
        
        this_dB_powerCRbw=zeros(no_trialsCR,4);
        this_dB_powerHitbw=zeros(no_trialsHit,4);
        roc_data=zeros(no_trialsHit+no_trialsCR,2);
        for ii=1:no_bandwidths
            %Enter the Hits
            this_band=(frequency>=low_freq(ii))&(frequency<=high_freq(ii));
            roc_data(1:no_trialsHit,1)=mean(this_dB_powerHit(:,this_band),2);
            roc_data(1:no_trialsHit,2)=zeros(no_trialsHit,1);
            
            %Enter CR
            roc_data(no_trialsHit+1:no_trialsHit+no_trialsCR,1)=mean(this_dB_powerCR(:,this_band),2);
            roc_data(no_trialsHit+1:no_trialsHit+no_trialsCR,2)=ones(no_trialsCR,1);
            no_ROCs=no_ROCs+1;
            ROCout(no_ROCs).roc=roc_calc(roc_data,0,0.05,0);
            ROCout(no_ROCs).fileNo=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
            ROCout(no_ROCs).groupNo=handles_drgb.drgbchoices.group_no(fileNo);
            ROCout(no_ROCs).timeWindow=handles_drgb.drgb.lfpevpair(lfpodNo).timeWindow;
            ROCout(no_ROCs).bandwidth=ii;
            p_vals_ROC=[p_vals_ROC ROCout(no_ROCs).roc.p];
            if timeWindow==2
                pffft=1;
            end
        end
        
        
    end
    
end


if (which_display==8)||(which_display==9)
    nRocs=0;
    pFDRroc=drsFDRpval(p_vals_ROC);
    d_primeHits=[];
    grHits=[];
    bwHits=[];
    d_primeCRs=[];
    grCRs=[];
    bwCRs=[];
    d_primeFAs=[];
    grFAs=[];
    bwFAs=[];
    for lfpodNo=1:no_lfpevpairs
        fileNo=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
        groupNo=handles_drgb.drgbchoices.group_no(fileNo);
        timeWindow=handles_drgb.drgb.lfpevpair(lfpodNo).timeWindow;
        lfpodorNo=lfpodNo;
        
        try
            trials_in_Hit=handles_drgb.drgb.lfpevpair(lfpodNo).which_eventLFPPower(evHit,:)==1;
        catch
            trials_in_Hit=[];
        end
        
        try
            trials_in_CR=handles_drgb.drgb.lfpevpair(lfpodNo).which_eventLFPPower(evCR,:)==1;
        catch
            trials_in_CR=[];
        end
        
        try
            trials_in_FA=handles_drgb.drgb.lfpevpair(lfpodNo).which_eventLFPPower(evFA,:)==1;
        catch
            trials_in_FA=[];
        end
        
        
        
        %Calculate the ROC for Hit vs CR
        %Do proficient only
        percent_bin=3;
        trials_in_perbin=(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrLFPPower>percent_low(percent_bin))&(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrLFPPower<=percent_high(percent_bin));
        
        if (sum(trials_in_Hit&trials_in_perbin)>5)&(sum(trials_in_CR&trials_in_perbin)>5)
            
            no_trialsHit=sum(trials_in_Hit&trials_in_perbin);
            this_dB_powerHit=zeros(no_trialsHit,length(handles_drgb.drgb.freq_for_LFPpower));
            this_dB_powerHit(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_Hit&trials_in_perbin,:));
            
            
            no_trialsCR=sum(trials_in_CR&trials_in_perbin);
            this_dB_powerCR=zeros(no_trialsCR,length(handles_drgb.drgb.freq_for_LFPpower));
            this_dB_powerCR(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_CR&trials_in_perbin,:));
            
            no_trialsFA=sum(trials_in_FA&trials_in_perbin);
            this_dB_powerFA=zeros(no_trialsFA,length(handles_drgb.drgb.freq_for_LFPpower));
            this_dB_powerFA(:,:)=10*log10(handles_drgb.drgb.lfpevpair(lfpodNo).allPower(trials_in_FA&trials_in_perbin,:));
            
            
            for ii=1:no_bandwidths
                %Enter the Hits
                
                nRocs=nRocs+1;
                
                if timeWindow==2
                    pffft=1;
                    %Calculate d prime if auROC is significant
                    if ROCout(nRocs).roc.p<=pFDRroc
                        %Enter the Hits
                        this_band=(frequency>=low_freq(ii))&(frequency<=high_freq(ii));
                        these_hits=mean(this_dB_powerHit(:,this_band),2);
                        mu_hit=mean(these_hits);
                        SDhit=std(these_hits);
                        these_CRs=mean(this_dB_powerCR(:,this_band),2);
                        mu_CR=mean(these_CRs);
                        SDCR=std(these_CRs);
                        if mu_hit>mu_CR
                            mult_fact=1/sqrt((SDhit^2 +SDCR^2)/2);
                        else
                            mult_fact=-1/sqrt((SDhit^2 +SDCR^2)/2);
                        end
                        
                        d_primeHits=[d_primeHits mult_fact*(these_hits-mu_CR)'];
                        grHits=[grHits ROCout(nRocs).groupNo*ones(1,no_trialsHit)];
                        bwHits=[bwHits ii*ones(1,no_trialsHit)];
                        d_primeCRs=[d_primeCRs mult_fact*(these_CRs-mu_CR)'];
                        grCRs=[grCRs ROCout(nRocs).groupNo*ones(1,no_trialsCR)];
                        bwCRs=[bwCRs ii*ones(1,no_trialsCR)];
                        if no_trialsFA>0
                            these_FAs=mean(this_dB_powerFA(:,this_band),2);
                            d_primeFAs=[d_primeFAs mult_fact*(these_FAs-mu_CR)'];
                            grFAs=[grFAs ROCout(nRocs).groupNo*ones(1,no_trialsFA)];
                            bwFAs=[bwFAs ii*ones(1,no_trialsFA)];
                        end
                    end
                end
            end
            
            
        end
        
    end
end



%Display the results for each badnwidth
%Sort by percent bin (rows) and within each row window No

figNo=1;
%Now do the figures and statistical analysis for the MI
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

these_circles{1}='*b';
these_circles{2}='*r';
these_circles{3}='*m';
these_circles{4}='*g';
these_circles{5}='*y';
these_circles{6}='*k';
these_circles{7}='*c';
these_circles{8}='*k';




switch which_display
    case 1
        %Display the difference between events
        for grNo=1:max(handles_drgb.drgbchoices.group_no)
            dB_p_v1_mean=zeros(no_percent_bins, no_event_types,length(handles_drgb.drgb.freq_for_LFPpower));
            dB_p_v1_SEM=zeros(no_percent_bins, no_event_types,length(handles_drgb.drgb.freq_for_LFPpower));
            
            
            frequency=handles_drgb.drgb.freq_for_LFPpower;
            
            
            kk=0;
            for per_bin=1:no_percent_bins
                
                %Calculate the p value for statistical differences
                evTNref=2; %p values will be calculated with respect to this event
                eventTypeRef=eventType(evTNref);
                these_p_vals=[];
                for evTN1=1:length(eventType)
                    
                    if evTN1~=evTNref
                        
                        eventType1=eventType(evTN1);
                        
                        if (no_values(eventType1,winNo,grNo,per_bin)>1)&(no_values(eventTypeRef,winNo,grNo,per_bin)>1)
                            calc_pval(evTN1)=1;
                            
                            this_dB_p_v1=zeros(no_values(eventType1,winNo,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1(:,:)=dB_power_values(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin),:);
                            if subtractRef==1
                                this_dB_p_v1_ref=zeros(no_values(eventType1,refWin,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_p_v1_ref(:,:)=dB_power_values(eventType1,refWin,grNo,per_bin,1:no_values(eventType1,refWin,grNo,per_bin),:);
                                this_dB_p_v1=this_dB_p_v1-this_dB_p_v1_ref;
                            end
                            
                            this_dB_p_v1ref=zeros(no_values(eventTypeRef,winNo,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1ref(:,:)=dB_power_values(eventTypeRef,winNo,grNo,per_bin,1:no_values(eventTypeRef,winNo,grNo,per_bin),:);
                            if subtractRef==1
                                this_dB_p_v1ref_ref=zeros(no_values(eventTypeRef,refWin,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_p_v1ref_ref(:,:)=dB_power_values(eventTypeRef,refWin,grNo,per_bin,1:no_values(eventTypeRef,refWin,grNo,per_bin),:);
                                this_dB_p_v1ref=this_dB_p_v1ref-this_dB_p_v1ref_ref;
                            end
                            
                            %Calculate whether this this_dB_p_v1 is significantly
                            %different from reference
                            for ifreq=1:length(frequency)
                                p_vals(grNo,per_bin,evTN1,ifreq)=ranksum(this_dB_p_v1(:,ifreq),this_dB_p_v1ref(:,ifreq));
                                these_p_vals=[these_p_vals ranksum(this_dB_p_v1(:,ifreq),this_dB_p_v1ref(:,ifreq))];
                            end
                            
                        else
                            calc_pval(evTN1)=0;
                        end
                        
                    else
                        calc_pval(evTN1)=0;
                    end
                end
                
                pFDR(grNo,per_bin)=drsFDRpval(these_p_vals);
                fprintf(1, 'pFDR for group No: %d, percent bin no: %d = %d\n',grNo,per_bin,pFDR(grNo,per_bin));
                
                %Now calculate the mean et al
                for evTN1=1:length(eventType)
                    eventType1=eventType(evTN1);
                    
                    
                    if no_values(eventType1,winNo,grNo,per_bin)>1
                        
                        this_dB_p_v1=zeros(no_values(eventType1,winNo,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                        this_dB_p_v1(:,:)=dB_power_values(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin),:);
                        if subtractRef==1
                            this_dB_p_v1_ref=zeros(no_values(eventType1,refWin,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1_ref(:,:)=dB_power_values(eventType1,refWin,grNo,per_bin,1:no_values(eventType1,refWin,grNo,per_bin),:);
                            this_dB_p_v1=this_dB_p_v1-this_dB_p_v1_ref;
                        end
                        dB_p_v1_mean(per_bin,evTN1,:)=mean(this_dB_p_v1,1);
                        dB_p_v1_SEM(per_bin,evTN1,:)=std(this_dB_p_v1,0,1)/sqrt(no_values(eventType1,winNo,grNo,per_bin));
                        %Calculate the 95% CI
                        for ifreq=1:length(frequency)
                            pd=fitdist(this_dB_p_v1(:,ifreq),'Normal');
                            ci=paramci(pd);
                            dB_p_v1_ci(per_bin,evTN1,ifreq)=pd.mu-ci(1,1);
                        end
                        these_experiments=experimentNo(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin));
                        
                    end
                    
                end
                
                
                
                
                
            end
            
            
            frequency=handles_drgb.drgb.freq_for_LFPpower;
            
            
            for per_bin=1:no_percent_bins
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                hold on
                
                %First draw lines and add legends
                this_legend=[];
                plot_handles=[];
                for evTN1=1:length(eventType)
                    this_dB_p_v1_mean=zeros(1,length(handles_drgb.drgb.freq_for_LFPpower));
                    this_dB_p_v1_mean(1,:)=dB_p_v1_mean(per_bin,evTN1,:);
                    eval(['p' num2str(evTN1) '=plot(frequency,this_dB_p_v1_mean, these_lines{evTN1});'])
                    %plot(frequency,this_dB_p_v1_mean, these_lines{evTN1});
                    plot_handles=[plot_handles 'p' num2str(evTN1) ' '];
                    this_legend=[this_legend '''' evTypeLabels{evTN1} '''' ','];
                end
                
                this_legend=['legend([' plot_handles(1:end-1) '],' this_legend(1:end-1) ')'];
                
                
                
                %Now do the bounded lines and the significance points
                for evTN1=1:length(eventType)
                    
                    this_dB_p_v1_mean=zeros(1,length(handles_drgb.drgb.freq_for_LFPpower));
                    this_dB_p_v1_mean(1,:)=dB_p_v1_mean(per_bin,evTN1,:);
                    
                    this_dB_p_v1_ci=zeros(1,length(handles_drgb.drgb.freq_for_LFPpower));
                    this_dB_p_v1_ci(1,:)=dB_p_v1_ci(per_bin,evTN1,:);
                    
                    [hl, hp] = boundedline(frequency,this_dB_p_v1_mean, this_dB_p_v1_ci, these_colors{evTN1});
                    
                    if calc_pval(evTN1)==1
                        for ifreq=1:length(frequency)
                            
                            
                            if p_vals(grNo,per_bin,evTN1,ifreq)<=pFDR(grNo,per_bin)
                                plot(frequency(ifreq),this_dB_p_v1_mean(ifreq),these_circles{evTN1})
                            end
                            
                            
                        end
                    end
                    
                end
                
                
                max_y=max(dB_p_v1_mean(:)+abs(dB_p_v1_ci(:)))+0.05*(max(dB_p_v1_mean(:)+abs(dB_p_v1_ci(:)))-min(dB_p_v1_mean(:)-abs(dB_p_v1_ci(:))));
                min_y=min(dB_p_v1_mean(:)-abs(dB_p_v1_ci(:)))-0.05*(max(dB_p_v1_mean(:)+abs(dB_p_v1_ci(:)))-min(dB_p_v1_mean(:)-abs(dB_p_v1_ci(:))));
                ylim([min_y max_y])
                
                
                if subtractRef==0
                    title(['Power (dB) for percent correct' percent_bin_legend{per_bin} ' group: ' handles_drgb.drgbchoices.group_no_names{grNo}])
                else
                    title(['delta Power (dB) for percent correct' percent_bin_legend{per_bin} ' group: ' handles_drgb.drgbchoices.group_no_names{grNo}])
                end
                
                eval(this_legend)
            end
        end
        
    case 2
        %Display the difference between groups
        these_p_vals=[];
        these_p_vals_perm2=[];
        
        dB_p_v1_mean=zeros(length(eventType),no_percent_bins, max(handles_drgb.drgbchoices.group_no),length(handles_drgb.drgb.freq_for_LFPpower));
        dB_p_v1_SEM=zeros(length(eventType),no_percent_bins, max(handles_drgb.drgbchoices.group_no),length(handles_drgb.drgb.freq_for_LFPpower));
        
        
        for evTN1=1:length(eventType)
            eventType1=eventType(evTN1);
            
            
            
            frequency=handles_drgb.drgb.freq_for_LFPpower;
            
            
            kk=0;
            for per_bin=1:no_percent_bins
                
                grRef=1;  %p values will be calculated with respect to this group
                
                
                for grNo=max(handles_drgb.drgbchoices.group_no):-1:1
                    
                    if grNo~=grRef
                        
                        if no_values(eventType1,winNo,grNo,per_bin)>1
                            calc_pval(grNo)=1;
                            
                            this_dB_p_v1=zeros(no_values(eventType1,winNo,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1(:,:)=dB_power_values(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin),:);
                            if subtractRef==1
                                this_dB_p_v1_ref=zeros(no_values(eventType1,refWin,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_p_v1_ref(:,:)=dB_power_values(eventType1,refWin,grNo,per_bin,1:no_values(eventType1,refWin,grNo,per_bin),:);
                                this_dB_p_v1=this_dB_p_v1-this_dB_p_v1_ref;
                            end
                            
                            this_dB_p_v1ref=zeros(no_values(eventType1,winNo,grRef,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1ref(:,:)=dB_power_values(eventType1,winNo,grRef,per_bin,1:no_values(eventType1,winNo,grRef,per_bin),:);
                            if subtractRef==1
                                this_dB_p_v1ref_ref=zeros(no_values(eventType1,refWin,grRef,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_p_v1ref_ref(:,:)=dB_power_values(eventType1,refWin,grRef,per_bin,1:no_values(eventType1,refWin,grRef,per_bin),:);
                                this_dB_p_v1ref=this_dB_p_v1ref-this_dB_p_v1ref_ref;
                            end
                            
                            %Calculate whether this this_dB_p_v1 is significantly
                            %different from reference
                            try
                                if subtractRef==1
                                    for ifreq=1:length(frequency)
                                        p_vals(evTN1,grNo,per_bin,evTN1,ifreq)=ranksum(this_dB_p_v1(:,ifreq),this_dB_p_v1ref(:,ifreq));
                                        these_p_vals=[these_p_vals ranksum(this_dB_p_v1(:,ifreq),this_dB_p_v1ref(:,ifreq))];
                                    end
                                    
                                    a={ this_dB_p_v1' this_dB_p_v1ref'};
                                    [F df pvals_perm2(evTN1,grNo,per_bin,evTN1,:)] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
                                    this_p_perm2=zeros(1,length(frequency));
                                    this_p_perm2(1,:)=pvals_perm2(grNo,per_bin,evTN1,:);
                                    these_p_vals_perm2=[these_p_vals_perm2 this_p_perm2];
                                end
                            catch
                            end
                            
                        else
                            calc_pval(grNo)=0;
                        end
                        
                        
                    else
                        calc_pval(grNo)=0;
                    end
                    
                end
                
                
                
                
                
                %Now calculate the means, CIs, etc
                for grNo=max(handles_drgb.drgbchoices.group_no):-1:1
                    
                    
                    
                    if no_values(eventType1,winNo,grNo,per_bin)>1
                        
                        this_dB_p_v1=zeros(no_values(eventType1,winNo,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                        this_dB_p_v1(:,:)=dB_power_values(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin),:);
                        if subtractRef==1
                            this_dB_p_v1_ref=zeros(no_values(eventType1,refWin,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1_ref(:,:)=dB_power_values(eventType1,refWin,grNo,per_bin,1:no_values(eventType1,refWin,grNo,per_bin),:);
                            this_dB_p_v1=this_dB_p_v1-this_dB_p_v1_ref;
                        end
                        dB_p_v1_mean(evTN1,per_bin,grNo,:)=mean(this_dB_p_v1,1);
                        dB_p_v1_SEM(evTN1,per_bin,grNo,:)=std(this_dB_p_v1,0,1)/sqrt(no_values(eventType1,winNo,grNo,per_bin));
                        %Calculate the 95% CI
                        for ifreq=1:length(frequency)
                            pd=fitdist(this_dB_p_v1(:,ifreq),'Normal');
                            ci=paramci(pd);
                            dB_p_v1_ci(evTN1,per_bin,grNo,ifreq)=pd.mu-ci(1,1);
                        end
                        these_experiments=experimentNo(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin));
                        
                    end
                end
                
                
                
            end
            
        end
        
        pFDR=drsFDRpval(these_p_vals);
        fprintf(1, 'pFDR  = %d\n',pFDR);
        
        
        pFDR_perm2=drsFDRpval(these_p_vals_perm2);
        fprintf(1, 'pFDR perm 2 = %d\n',pFDR_perm2);
        frequency=handles_drgb.drgb.freq_for_LFPpower;
        for evTN1=1:length(eventType)
            eventType1=eventType(evTN1);
            %Now plot the lines
            for per_bin=1:no_percent_bins
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                hold on
                
                %Plot the lines
                this_legend=[];
                plot_handles=[];
                for grNo=max(handles_drgb.drgbchoices.group_no):-1:1
                    
                    
                    
                    
                    %[hl, hp] = boundedline(frequency,dB_p_v1_mean(jj,:), dB_p_v1_SEM(jj,:), 'b');
                    this_dB_p_v1_mean=zeros(1,length(handles_drgb.drgb.freq_for_LFPpower));
                    this_dB_p_v1_mean(1,:)=dB_p_v1_mean(evTN1,per_bin,grNo,:);
                    
                    this_dB_p_v1_ci=zeros(1,length(handles_drgb.drgb.freq_for_LFPpower));
                    this_dB_p_v1_ci(1,:)=dB_p_v1_ci(evTN1,per_bin,grNo,:);
                    
                    %plot(frequency,this_dB_p_v1_mean, these_lines{grNo});
                    eval(['p' num2str(grNo) '=plot(frequency,this_dB_p_v1_mean, these_lines{grNo});'])
                    
                    
                    this_legend=[this_legend '''' handles_drgb.drgbchoices.group_no_names{grNo} '''' ' ,'];
                    plot_handles=[plot_handles 'p' num2str(grNo) ' '];
                    
                end
                
                this_legend=['legend([' plot_handles(1:end-1) '],' this_legend(1:end-1) ')'];
                
                
                
                
                %Now plot the bounded lines
                
                for grNo=max(handles_drgb.drgbchoices.group_no):-1:1
                    
                    
                    this_dB_p_v1_mean=zeros(1,length(handles_drgb.drgb.freq_for_LFPpower));
                    this_dB_p_v1_mean(1,:)=dB_p_v1_mean(evTN1,per_bin,grNo,:);
                    
                    this_dB_p_v1_ci=zeros(1,length(handles_drgb.drgb.freq_for_LFPpower));
                    this_dB_p_v1_ci(1,:)=dB_p_v1_ci(evTN1,per_bin,grNo,:);
                    
                    [hl, hp] = boundedline(frequency,this_dB_p_v1_mean, this_dB_p_v1_ci, these_colors{grNo});
                    
                    if subtractRef==1
                        if calc_pval(grNo)==1
                            for ifreq=1:length(frequency)
                                
                                switch stat_method
                                    case 1
                                        if p_vals(evTN1,grNo,per_bin,evTN1,ifreq)<=pFDR
                                            plot(frequency(ifreq),this_dB_p_v1_mean(ifreq),these_circles{grNo})
                                        end
                                    case 2
                                        if pvals_perm2(evTN1,grNo,per_bin,evTN1,ifreq)<=pFDR_perm2
                                            plot(frequency(ifreq),this_dB_p_v1_mean(ifreq),these_circles{grNo})
                                        end
                                end
                                
                            end
                        end
                    end
                    
                    
                end
                
                eval(this_legend)
                
                max_y=max(dB_p_v1_mean(:)+abs(dB_p_v1_ci(:)))+0.05*(max(dB_p_v1_mean(:)+abs(dB_p_v1_ci(:)))-min(dB_p_v1_mean(:)-abs(dB_p_v1_ci(:))));
                min_y=min(dB_p_v1_mean(:)-abs(dB_p_v1_ci(:)))-0.05*(max(dB_p_v1_mean(:)+abs(dB_p_v1_ci(:)))-min(dB_p_v1_mean(:)-abs(dB_p_v1_ci(:))));
                ylim([min_y max_y])
                
                
                
                if subtractRef==0
                    title(['Power (dB) for percent correct' percent_bin_legend{per_bin} ' event: ' evTypeLabels{evTN1}])
                else
                    title(['delta Power (dB) for percent correct' percent_bin_legend{per_bin} ' event: ' evTypeLabels{evTN1}])
                end
            end
        end
        
    case 3
        %Display the difference in dB power vs frequency between learning and proficient
        for grNo=1:max(handles_drgb.drgbchoices.group_no)
            dB_p_v1_mean=zeros(no_percent_bins, no_event_types,length(handles_drgb.drgb.freq_for_LFPpower));
            dB_p_v1_SEM=zeros(no_percent_bins, no_event_types,length(handles_drgb.drgb.freq_for_LFPpower));
            
            
            frequency=handles_drgb.drgb.freq_for_LFPpower;
            
            for evTN1=1:length(eventType)
                
                %Calculate the p value for statistical differences
                per_bin_ref=1; %p values will be calculated with respect to this percent bin
                eventType1=eventType(evTN1);
                these_p_vals=[];
                
                for per_bin=1:no_percent_bins
                    
                    if per_bin~=per_bin_ref
                        
                        if (no_values(eventType1,winNo,grNo,per_bin)>1)&(no_values(eventType1,winNo,grNo,per_bin_ref)>1)
                            calc_pval(per_bin)=1;
                            
                            this_dB_p_v1=zeros(no_values(eventType1,winNo,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1(:,:)=dB_power_values(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin),:);
                            if subtractRef==1
                                this_dB_p_v1_ref=zeros(no_values(eventType1,refWin,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_p_v1_ref(:,:)=dB_power_values(eventType1,refWin,grNo,per_bin,1:no_values(eventType1,refWin,grNo,per_bin),:);
                                this_dB_p_v1=this_dB_p_v1-this_dB_p_v1_ref;
                            end
                            
                            these_experiments=zeros(1,no_values(eventType1,winNo,grNo,per_bin));
                            these_experiments(:,:)=experimentNo(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin),:);
                            
                            this_dB_p_v1ref=zeros(no_values(eventType1,winNo,grNo,per_bin_ref),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1ref(:,:)=dB_power_values(eventType1,winNo,grNo,per_bin_ref,1:no_values(eventType1,winNo,grNo,per_bin_ref),:);
                            if subtractRef==1
                                this_dB_p_v1ref_ref=zeros(no_values(eventType1,refWin,grNo,per_bin_ref),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_p_v1ref_ref(:,:)=dB_power_values(eventType1,refWin,grNo,per_bin_ref,1:no_values(eventType1,refWin,grNo,per_bin_ref),:);
                                this_dB_p_v1ref=this_dB_p_v1ref-this_dB_p_v1ref_ref;
                            end
                            
                            ref_experiments=zeros(1,no_values(eventType1,winNo,grNo,per_bin_ref));
                            ref_experiments(:,:)=experimentNo(eventType1,winNo,grNo,per_bin_ref,1:no_values(eventType1,winNo,grNo,per_bin_ref),:);
                            
                            %Calculate whether this this_dB_p_v1 is significantly
                            %different from reference
                            for ifreq=1:length(frequency)
                                %Ranksum
                                p_vals(grNo,per_bin,evTN1,ifreq)=ranksum(this_dB_p_v1(:,ifreq),this_dB_p_v1ref(:,ifreq));
                                these_p_vals=[these_p_vals ranksum(this_dB_p_v1(:,ifreq),this_dB_p_v1ref(:,ifreq))];
                            end
                            
                            %Perform permutation test using mult_comp_perm_t1
                            
                            %This is a paired test. Find experiments where both
                            %this LFP and reference LFP were recorded
                            this_paired_dB_p_v1=[];
                            this_paired_dB_p_v1ref=[];
                            for expNo=min([min(ref_experiments) min(these_experiments)]):max([max(ref_experiments) max(these_experiments)])
                                if sum(expNo==these_experiments)>0
                                    if sum(expNo==ref_experiments)>0
                                        this_paired_dB_p_v1=[this_paired_dB_p_v1 this_dB_p_v1(these_experiments==expNo,:)'];
                                        this_paired_dB_p_v1ref=[this_paired_dB_p_v1ref this_dB_p_v1ref(ref_experiments==expNo,:)'];
                                    end
                                end
                            end
                            this_paired_dB_p_v1=this_paired_dB_p_v1';
                            this_paired_dB_p_v1ref=this_paired_dB_p_v1ref';
                            
                            %Perform permutation test
                            dif=this_paired_dB_p_v1-this_paired_dB_p_v1ref;
                            [pvals_perm(grNo,per_bin,evTN1,:), t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t1(dif,50000);
                            
                            %Perform permutation test using statcond()
                            %         a = { rand(2,11) rand(2,10) rand(2,12)+0.5 }; % pseudo 'unpaired'
                            %         [F df pvals] = statcond(a); % perform an unpaired ANOVA
                            %           pvals =
                            %              0.00025 % p-values for difference between columns
                            %              0.00002 % for each data row
                            
                            a={ this_dB_p_v1' this_dB_p_v1ref'};
                            [F df pvals_perm2(grNo,per_bin,evTN1,:)] = statcond(a,'mode',mode_statcond,'naccu', 1000); % perform an unpaired ANOVA
                            
                            %anovan with repeated measures
                            %To use anovan for repeated measures, you have to enter the subject
                            %id's as a new factor and specify that factor as a random variable in
                            %anovan. If you have different subjects in different groups, you must
                            %also specify that the subject index is nested in the group variable
                            %Note: This has not been implemented
                            
                            %                             %Setup anovan variables
                            %                             these_dBs=[];
                            %                             these_bins=[];
                            %                             these_freqs=[];
                            %
                            %                             for ii=1:no_values(eventType1,winNo,grNo,per_bin)
                            %                                 these_dBs=[these_dBs this_dB_p_v1(ii,:)];
                            %                                 these_bins=[these_bins ones(1,length(frequency))];
                            %                                 these_freqs=[these_freqs frequency'];
                            %                             end
                            %
                            %                             for ii=1:no_values(eventType1,winNo,grNo,per_bin_ref)
                            %                                 these_dBs=[these_dBs this_dB_p_v1ref(ii,:)];
                            %                                 these_bins=[these_bins 2*ones(1,length(frequency))];
                            %                                 these_freqs=[these_freqs frequency'];
                            %                             end
                            %
                            %                             [p,tbl,stats]=anovan(these_dBs,{these_bins,these_freqs},'display','off');
                            %                             p_anova(grNo,per_bin,evTN1,ifreq)=p(1);
                            %
                            %                             if per_bin==3
                            %                                 pffft=1;
                            %                             end
                        else
                            calc_pval(per_bin)=0;
                        end
                        
                    else
                        calc_pval(per_bin)=0;
                    end
                end
                
                pFDR(grNo,evTN1)=drsFDRpval(these_p_vals);
                fprintf(1, ['pFDR for group No: %d, event: ' evTypeLabels{evTN1} '= %d\n'], grNo, pFDR(grNo,evTN1));
                
                %Now calculate the mean et al
                for per_bin=1:no_percent_bins
                    
                    if no_values(eventType1,winNo,grNo,per_bin)>1
                        
                        this_dB_p_v1=zeros(no_values(eventType1,winNo,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                        this_dB_p_v1(:,:)=dB_power_values(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin),:);
                        if subtractRef==1
                            this_dB_p_v1_ref=zeros(no_values(eventType1,refWin,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1_ref(:,:)=dB_power_values(eventType1,refWin,grNo,per_bin,1:no_values(eventType1,refWin,grNo,per_bin),:);
                            this_dB_p_v1=this_dB_p_v1-this_dB_p_v1_ref;
                        end
                        dB_p_v1_mean(per_bin,evTN1,:)=mean(this_dB_p_v1,1);
                        dB_p_v1_SEM(per_bin,evTN1,:)=std(this_dB_p_v1,0,1)/sqrt(no_values(eventType1,winNo,grNo,per_bin));
                        %Calculate the 95% CI
                        for ifreq=1:length(frequency)
                            pd=fitdist(this_dB_p_v1(:,ifreq),'Normal');
                            ci=paramci(pd);
                            dB_p_v1_ci(per_bin,evTN1,ifreq)=pd.mu-ci(1,1);
                        end
                        these_experiments=experimentNo(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin));
                        
                    end
                    
                end
                
                
                
                
                
            end
            
            
            frequency=handles_drgb.drgb.freq_for_LFPpower;
            
            
            for evTN1=1:length(eventType)
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                hold on
                
                %First draw lines and add legends
                this_legend=[];
                plot_handles=[];
                for per_bin=[1 3]
                    this_dB_p_v1_mean=zeros(1,length(handles_drgb.drgb.freq_for_LFPpower));
                    this_dB_p_v1_mean(1,:)=dB_p_v1_mean(per_bin,evTN1,:);
                    eval(['p' num2str(per_bin) '=plot(frequency,this_dB_p_v1_mean, these_lines{per_bin});'])
                    %plot(frequency,this_dB_p_v1_mean, these_lines{evTN1});
                    plot_handles=[plot_handles 'p' num2str(per_bin) ' '];
                    this_legend=[this_legend '''' percent_bin_legend{per_bin} '''' ','];
                end
                
                this_legend=['legend([' plot_handles(1:end-1) '],' this_legend(1:end-1) ')'];
                
                
                
                %Now do the bounded lines and the significance points
                for per_bin=[1 3]
                    
                    this_dB_p_v1_mean=zeros(1,length(handles_drgb.drgb.freq_for_LFPpower));
                    this_dB_p_v1_mean(1,:)=dB_p_v1_mean(per_bin,evTN1,:);
                    
                    this_dB_p_v1_ci=zeros(1,length(handles_drgb.drgb.freq_for_LFPpower));
                    this_dB_p_v1_ci(1,:)=dB_p_v1_ci(per_bin,evTN1,:);
                    
                    [hl, hp] = boundedline(frequency,this_dB_p_v1_mean, this_dB_p_v1_ci, these_colors{per_bin});
                    
                    if calc_pval(per_bin)==1
                        for ifreq=1:length(frequency)
                            
                            switch stat_method
                                case 1
                                    %FDR
                                    if p_vals(grNo,per_bin,evTN1,ifreq)<=pFDR(grNo,evTN1)
                                        plot(frequency(ifreq),this_dB_p_v1_mean(ifreq),these_circles{per_bin})
                                    end
                                    
                                case 2
                                    %Permutation test
                                    %                             if pvals_perm(grNo,per_bin,evTN1,ifreq)<=0.05
                                    %                                 plot(frequency(ifreq),this_dB_p_v1_mean(ifreq),these_circles{per_bin})
                                    %                             end
                                    
                                    %Permutation test
                                    if pvals_perm2(grNo,per_bin,evTN1,ifreq)<=0.05
                                        plot(frequency(ifreq),this_dB_p_v1_mean(ifreq),these_circles{per_bin})
                                    end
                            end
                        end
                    end
                    
                end
                
                
                max_y=max(dB_p_v1_mean(:)+abs(dB_p_v1_ci(:)))+0.05*(max(dB_p_v1_mean(:)+abs(dB_p_v1_ci(:)))-min(dB_p_v1_mean(:)-abs(dB_p_v1_ci(:))));
                min_y=min(dB_p_v1_mean(:)-abs(dB_p_v1_ci(:)))-0.05*(max(dB_p_v1_mean(:)+abs(dB_p_v1_ci(:)))-min(dB_p_v1_mean(:)-abs(dB_p_v1_ci(:))));
                ylim([min_y max_y])
                
                
                if subtractRef==0
                    title(['Power (dB) for event ' evTypeLabels{evTN1} ' group: ' handles_drgb.drgbchoices.group_no_names{grNo}])
                else
                    title(['delta Power (dB) for event ' evTypeLabels{evTN1} ' group: ' handles_drgb.drgbchoices.group_no_names{grNo}])
                end
                
                eval(this_legend)
            end
        end
        
    case 4
        
        %Boxplot for dB per badwidth displaying the difference between learning and proficient
        for grNo=1:max(handles_drgb.drgbchoices.group_no)
            dB_p_v1_mean=zeros(no_percent_bins, no_event_types,no_bandwidths);
            dB_p_v1_SEM=zeros(no_percent_bins, no_event_types,no_bandwidths);
            dB_p_v1_all=[];
            
            frequency=handles_drgb.drgb.freq_for_LFPpower;
            
            
            for evTN1=1:length(eventType)
                
                %Calculate the p value for statistical differences
                per_bin_ref=1; %p values will be calculated with respect to this percent bin
                eventType1=eventType(evTN1);
                these_p_vals=[];
                
                for per_bin=1:no_percent_bins
                    
                    if per_bin~=per_bin_ref
                        
                        if (no_values(eventType1,winNo,grNo,per_bin)>1)&(no_values(eventType1,winNo,grNo,per_bin_ref)>1)
                            calc_pval(per_bin)=1;
                            
                            this_dB_p_v1=zeros(no_values(eventType1,winNo,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1(:,:)=dB_power_values(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin),:);
                            if subtractRef==1
                                this_dB_p_v1_ref=zeros(no_values(eventType1,refWin,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_p_v1_ref(:,:)=dB_power_values(eventType1,refWin,grNo,per_bin,1:no_values(eventType1,refWin,grNo,per_bin),:);
                                this_dB_p_v1=this_dB_p_v1-this_dB_p_v1_ref;
                            end
                            
                            this_dB_p_bw=zeros(no_values(eventType1,winNo,grNo,per_bin),no_bandwidths);
                            for ii=1:no_bandwidths
                                this_band=(frequency>=low_freq(ii))&(frequency<=high_freq(ii));
                                this_dB_p_bw(:,ii)=mean(this_dB_p_v1(:,this_band),2);
                            end
                            
                            
                            this_dB_p_v1ref=zeros(no_values(eventType1,winNo,grNo,per_bin_ref),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1ref(:,:)=dB_power_values(eventType1,winNo,grNo,per_bin_ref,1:no_values(eventType1,winNo,grNo,per_bin_ref),:);
                            if subtractRef==1
                                this_dB_p_v1ref_ref=zeros(no_values(eventType1,refWin,grNo,per_bin_ref),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_p_v1ref_ref(:,:)=dB_power_values(eventType1,refWin,grNo,per_bin_ref,1:no_values(eventType1,refWin,grNo,per_bin_ref),:);
                                this_dB_p_v1ref=this_dB_p_v1ref-this_dB_p_v1ref_ref;
                            end
                            
                            this_dB_p_bwref=zeros(no_values(eventType1,refWin,grNo,per_bin_ref),no_bandwidths);
                            for ii=1:no_bandwidths
                                this_band=(frequency>=low_freq(ii))&(frequency<=high_freq(ii));
                                this_dB_p_bwref(:,ii)=mean(this_dB_p_v1ref(:,this_band),2);
                            end
                            
                            %Calculate whether this this_dB_p_v1 is significantly
                            %different from reference
                            for bwii=1:no_bandwidths
                                %Ranksum
                                p_vals(grNo,per_bin,evTN1,bwii)=ranksum(this_dB_p_bw(:,bwii),this_dB_p_bwref(:,bwii));
                                these_p_vals=[these_p_vals ranksum(this_dB_p_v1(:,bwii),this_dB_p_v1ref(:,bwii))];
                                
                                fprintf(1, ['p value for ' freq_names{bwii} ', group No: %d, event: ' evTypeLabels{evTN1} ...
                                    ', ' percent_bin_legend{per_bin} '= %d\n'], grNo, p_vals(grNo,per_bin,evTN1,bwii));
                                
                            end
                            fprintf(1, '\n')
                            
                        else
                            calc_pval(per_bin)=0;
                        end
                        
                    else
                        calc_pval(per_bin)=0;
                    end
                end
                
                pFDR(grNo,evTN1)=drsFDRpval(these_p_vals);
                fprintf(1, ['pFDR for group No: %d, event: ' evTypeLabels{evTN1} '= %d\n\n'], grNo, pFDR(grNo,evTN1));
                
                %Now do the boxplot
                for ii=1:no_bandwidths
                    figNo=figNo+1;
                    try
                        close(figNo)
                    catch
                    end
                    figure(figNo)
                    hold on
                    
                    
                    xtickl=[];
                    kk=0;
                    this_dB_p=[];
                    these_experiments=[];
                    bin_group=[];
                    for per_bin=[1 3]
                        
                        if no_values(eventType1,winNo,grNo,per_bin)>1
                            
                            this_dB_p_v1=zeros(no_values(eventType1,winNo,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1(:,:)=dB_power_values(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin),:);
                            if subtractRef==1
                                this_dB_p_v1_ref=zeros(no_values(eventType1,refWin,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_p_v1_ref(:,:)=dB_power_values(eventType1,refWin,grNo,per_bin,1:no_values(eventType1,refWin,grNo,per_bin),:);
                                this_dB_p_v1=this_dB_p_v1-this_dB_p_v1_ref;
                            end
                            
                            
                            this_band=(frequency>=low_freq(ii))&(frequency<=high_freq(ii));
                            expno=zeros(1,no_values(eventType1,winNo,grNo,per_bin));
                            expno(1,:)=experimentNo(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin));
                            these_experiments=[these_experiments expno];
                            this_dB_p=[this_dB_p mean(this_dB_p_v1(:,this_band),2)'];
                            bin_group=[bin_group per_bin*ones(1,no_values(eventType1,winNo,grNo,per_bin))];
                            
                            kk=kk+1;
                            xtickl{kk}=percent_bin_legend{per_bin};
                            
                            %Now plot the individual points
                            deltax=0.03;
                            x_shift=-deltax*floor((max(these_experiments(bin_group==per_bin))-min(these_experiments(bin_group==per_bin)))/2);
                            for jj=min(these_experiments(bin_group==per_bin)):max(these_experiments(bin_group==per_bin))
                                if x_shift==0
                                    x_shift=x_shift+deltax;
                                end
                                no_lfps=sum((these_experiments==jj)&(bin_group==per_bin));
                                x=(kk+x_shift)*ones(1,no_lfps);
                                plot(x,this_dB_p((these_experiments==jj)&(bin_group==per_bin)),'.k')
                                x_shift=x_shift+deltax;
                            end
                            
                        end
                        
                    end
                    
                    boxplot(this_dB_p,bin_group,'Symbol','')
                    xticklabels(xtickl)
                    if subtractRef==0
                        title([ freq_names{ii} ' power (dB) for event ' evTypeLabels{evTN1} ' group: ' handles_drgb.drgbchoices.group_no_names{grNo}])
                    else
                        title([ freq_names{ii}  ' delta Power (dB) for event ' evTypeLabels{evTN1} ' group: ' handles_drgb.drgbchoices.group_no_names{grNo}])
                    end
                    
                    pffft=1;
                end
                
                
                
            end
            
            
            
            
            
            
        end
        
    case 5
        %Display the difference between groups using boxplots
        %For Daniel's acetophenone vs. ethyl benzoate for per_bin=no_percent_bins:no_percent_bins
        
        %First estimate the p values for differences between groups
        frequency=handles_drgb.drgb.freq_for_LFPpower;
        dB_p_v1_mean=zeros(length(eventType),no_percent_bins, max(handles_drgb.drgbchoices.group_no),length(handles_drgb.drgb.freq_for_LFPpower));
        dB_p_v1_SEM=zeros(length(eventType),no_percent_bins, max(handles_drgb.drgbchoices.group_no),length(handles_drgb.drgb.freq_for_LFPpower));
        these_p_vals=[];
        these_p_vals_perm=[];
        
        
        
        for evTN1=1:length(eventType)
            eventType1=eventType(evTN1);
            
            
            
            frequency=handles_drgb.drgb.freq_for_LFPpower;
            
            
            kk=0;
            for per_bin=no_percent_bins:no_percent_bins
                
                
                these_p_vals=[];
                for grRef=1:max(handles_drgb.drgbchoices.group_no)
                    fprintf(1, ['\n\n\np values for all groups compared to group ' handles_drgb.drgbchoices.group_no_names{grRef} '\n'],  grRef);
                    fprintf(1, ['calculated for mice ' percent_bin_legend{per_bin} ' and for ' evTypeLabels{evTN1}  ' event\n\n']);
                    for grNo=1:max(handles_drgb.drgbchoices.group_no)
                        
                        if grNo~=grRef
                            
                            if subtractRef==1
                                continue_for=(no_values(eventType1,winNo,grNo,per_bin)>1)&(no_values(eventType1,winNo,grRef,per_bin));
                            else
                                continue_for=(no_values(eventType1,winNo,grNo,per_bin)>1);
                            end
                            
                            if continue_for
                                calc_pval(grNo)=1;
                                
                                this_dB_p_v1=zeros(no_values(eventType1,winNo,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_p_v1(:,:)=dB_power_values(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin),:);
                                if subtractRef==1
                                    this_dB_p_v1_ref=zeros(no_values(eventType1,refWin,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                                    this_dB_p_v1_ref(:,:)=dB_power_values(eventType1,refWin,grNo,per_bin,1:no_values(eventType1,refWin,grNo,per_bin),:);
                                    this_dB_p_v1=this_dB_p_v1-this_dB_p_v1_ref;
                                end
                                
                                this_dB_p_bw=zeros(no_values(eventType1,winNo,grNo,per_bin),no_bandwidths);
                                for ii=1:no_bandwidths
                                    this_band=(frequency>=low_freq(ii))&(frequency<=high_freq(ii));
                                    this_dB_p_bw(:,ii)=mean(this_dB_p_v1(:,this_band),2);
                                end
                                
                                this_dB_p_v1ref=zeros(no_values(eventType1,winNo,grRef,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_p_v1ref(:,:)=dB_power_values(eventType1,winNo,grRef,per_bin,1:no_values(eventType1,winNo,grRef,per_bin),:);
                                if subtractRef==1
                                    this_dB_p_v1ref_ref=zeros(no_values(eventType1,refWin,grRef,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                                    this_dB_p_v1ref_ref(:,:)=dB_power_values(eventType1,refWin,grRef,per_bin,1:no_values(eventType1,refWin,grRef,per_bin),:);
                                    this_dB_p_v1ref=this_dB_p_v1ref-this_dB_p_v1ref_ref;
                                end
                                
                                this_dB_p_bwref=zeros(no_values(eventType1,winNo,grRef,per_bin),no_bandwidths);
                                for ii=1:no_bandwidths
                                    this_band=(frequency>=low_freq(ii))&(frequency<=high_freq(ii));
                                    this_dB_p_bwref(:,ii)=mean(this_dB_p_v1ref(:,this_band),2);
                                end
                                
                                %Calculate whether this this_dB_p_v1 is significantly
                                %different from reference
                                
                                switch stat_method
                                    case 1
                                        %Ranksum
                                        for bwii=1:no_bandwidths
                                            
                                            
                                            p_vals(grNo,per_bin,evTN1,bwii)=ranksum(this_dB_p_bw(:,bwii),this_dB_p_bwref(:,bwii));
                                            these_p_vals=[these_p_vals ranksum(this_dB_p_bw(:,bwii),this_dB_p_bwref(:,bwii))];
                                            fprintf(1, ['p value for ' freq_names{bwii} ', ' handles_drgb.drgbchoices.group_no_names{grNo} '= %d\n'],  p_vals(grNo,per_bin,evTN1,bwii));
                                        end
                                    case 2
                                        %statcond
                                        a={ this_dB_p_bw' this_dB_p_bwref'};
                                        [F df pvals_perm2(grNo,per_bin,evTN1,:)] = statcond(a,'mode',mode_statcond); % perform an unpaired ANOVA
                                        for bwii=1:no_bandwidths
                                            fprintf(1, ['p value for ' freq_names{bwii} ', ' handles_drgb.drgbchoices.group_no_names{grNo} '= %d\n'],  pvals_perm2(grNo,per_bin,evTN1,bwii));
                                        end
                                        
                                        this_pv=zeros(1,no_bandwidths);
                                        this_pv(:,:)=pvals_perm2(grNo,per_bin,evTN1,:);
                                        these_p_vals_perm=[these_p_vals_perm this_pv];
                                end
                                
                                fprintf(1, '\n')
                                
                            else
                                calc_pval(grNo)=0;
                            end
                            
                        end
                        
                    end
                end
                
            end
        end
        
        switch stat_method
            case 1
                %Ranksum
                
                pFDR=drsFDRpval(these_p_vals);
                fprintf(1, ['pFDR  = %d\n'],pFDR);
                
            case 2
                %statcond
                pFDR_perm=drsFDRpval(these_p_vals_perm);
                fprintf(1, 'pFDR perm  = %d\n',pFDR_perm);
        end
        
        for evTN1=1:length(eventType)
            eventType1=eventType(evTN1);
            
            %Now do the boxplot
            for per_bin=no_percent_bins:no_percent_bins
                for ii=1:no_bandwidths
                    figNo=figNo+1;
                    try
                        close(figNo)
                    catch
                    end
                    figure(figNo)
                    hold on
                    
                    
                    xtickl=[];
                    kk=0;
                    this_dB_p=[];
                    these_experiments=[];
                    group_no=[];
                    for grNo=1:max(handles_drgb.drgbchoices.group_no)
                        
                        if no_values(eventType1,winNo,grNo,per_bin)>1
                            
                            this_dB_p_v1=zeros(no_values(eventType1,winNo,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1(:,:)=dB_power_values(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin),:);
                            if subtractRef==1
                                this_dB_p_v1_ref=zeros(no_values(eventType1,refWin,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_p_v1_ref(:,:)=dB_power_values(eventType1,refWin,grNo,per_bin,1:no_values(eventType1,refWin,grNo,per_bin),:);
                                this_dB_p_v1=this_dB_p_v1-this_dB_p_v1_ref;
                            end
                            
                            
                            this_band=(frequency>=low_freq(ii))&(frequency<=high_freq(ii));
                            expno=zeros(1,no_values(eventType1,winNo,grNo,per_bin));
                            expno(1,:)=experimentNo(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin));
                            these_experiments=[these_experiments expno];
                            this_dB_p=[this_dB_p mean(this_dB_p_v1(:,this_band),2)'];
                            group_no=[group_no grNo*ones(1,no_values(eventType1,winNo,grNo,per_bin))];
                            
                            kk=kk+1;
                            xtickl{kk}=handles_drgb.drgbchoices.group_no_names{grNo};
                            
                            %Now plot the individual points
                            deltax=0.03;
                            no_these_files=0;
                            for ww=1:length(group_per_file)
                                if sum(ww==these_experiments(group_no==grNo))>0
                                    no_these_files=no_these_files+1;
                                    expNo(no_these_files)=ww;
                                end
                            end
                            x_shift=-deltax*floor(no_these_files/2);
                            for jj=1:no_these_files
                                if x_shift==0
                                    x_shift=x_shift+deltax;
                                end
                                no_lfps=sum((these_experiments==expNo(jj))&(group_no==grNo));
                                x=(kk+x_shift)*ones(1,no_lfps);
                                plot(x,this_dB_p((these_experiments==expNo(jj))&(group_no==grNo)),'.k')
                                x_shift=x_shift+deltax;
                            end
                            
                            
                        end
                        
                    end
                    
                    boxplot(this_dB_p,group_no,'Symbol','')
                    xticklabels(xtickl)
                    if subtractRef==0
                        title([ freq_names{ii} ' power (dB) for event ' evTypeLabels{evTN1} ' ' percent_bin_legend{per_bin}])
                    else
                        title([ freq_names{ii}  ' delta Power (dB) for event ' evTypeLabels{evTN1} ' ' percent_bin_legend{per_bin}])
                    end
                    
                    pffft=1;
                    
                end
            end
            
            
            
            
        end
        
        
        
    case 6
        %Analysis of delta power (Hit-CR)
        %Display the difference between groups using boxplots
        %For Daniel's acetophenone vs. ethyl benzoate for per_bin=no_percent_bins:no_percent_bins
        
        %First estimate the p values for differences between groups
        frequency=handles_drgb.drgb.freq_for_LFPpower;
        dB_p_v1_mean=zeros(length(eventType),no_percent_bins, max(handles_drgb.drgbchoices.group_no),length(handles_drgb.drgb.freq_for_LFPpower));
        dB_p_v1_SEM=zeros(length(eventType),no_percent_bins, max(handles_drgb.drgbchoices.group_no),length(handles_drgb.drgb.freq_for_LFPpower));
        these_p_vals=[];
        these_p_vals_perm=[];
        
        
        
        
        
        
        
        
        frequency=handles_drgb.drgb.freq_for_LFPpower;
        
        
        kk=0;
        for per_bin=no_percent_bins:no_percent_bins
            
            
            these_p_vals=[];
            for grRef=1:max(handles_drgb.drgbchoices.group_no)
                fprintf(1, ['\n\n\np values for all groups compared to group ' handles_drgb.drgbchoices.group_no_names{grRef} '\n'],  grRef);
                fprintf(1, ['calculated for mice ' percent_bin_legend{per_bin} ' and for delta ' evTypeLabel  ' minus ' evTypeRefLabel '\n\n']);
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    if grNo~=grRef
                        
                        if subtractRef==1
                            continue_for=(no_valuesHitCR(eventType1,winNo,grNo,per_bin)>1)&(no_valuesHitCR(eventType1,winNo,grRef,per_bin));
                        else
                            continue_for=(no_valuesHitCR(eventType1,winNo,grNo,per_bin)>1);
                        end
                        
                        if continue_for
                            calc_pval(grNo)=1;
                            
                            this_dB_p_v1=zeros(no_valuesHitCR(eventType1,winNo,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1(:,:)=dB_power_valuesHitCR(eventType1,winNo,grNo,per_bin,1:no_valuesHitCR(eventType1,winNo,grNo,per_bin),:);
                            if subtractRef==1
                                this_dB_p_v1_ref=zeros(no_valuesHitCR(eventTypeRef,winNo,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_p_v1_ref(:,:)=dB_power_valuesHitCR(eventTypeRef,winNo,grNo,per_bin,1:no_valuesHitCR(eventTypeRef,winNo,grNo,per_bin),:);
                                this_dB_p_v1=this_dB_p_v1-this_dB_p_v1_ref;
                            end
                            
                            this_dB_p_bw=zeros(no_valuesHitCR(eventType1,winNo,grNo,per_bin),no_bandwidths);
                            for ii=1:no_bandwidths
                                this_band=(frequency>=low_freq(ii))&(frequency<=high_freq(ii));
                                this_dB_p_bw(:,ii)=mean(this_dB_p_v1(:,this_band),2);
                            end
                            
                            this_dB_p_v1ref=zeros(no_valuesHitCR(eventType1,winNo,grRef,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1ref(:,:)=dB_power_valuesHitCR(eventType1,winNo,grRef,per_bin,1:no_valuesHitCR(eventType1,winNo,grRef,per_bin),:);
                            if subtractRef==1
                                this_dB_p_v1ref_ref=zeros(no_valuesHitCR(eventTypeRef,winNo,grRef,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_p_v1ref_ref(:,:)=dB_power_valuesHitCR(eventTypeRef,winNo,grRef,per_bin,1:no_valuesHitCR(eventTypeRef,winNo,grRef,per_bin),:);
                                this_dB_p_v1ref=this_dB_p_v1ref-this_dB_p_v1ref_ref;
                            end
                            
                            this_dB_p_bwref=zeros(no_valuesHitCR(eventType1,winNo,grRef,per_bin),no_bandwidths);
                            for ii=1:no_bandwidths
                                this_band=(frequency>=low_freq(ii))&(frequency<=high_freq(ii));
                                this_dB_p_bwref(:,ii)=mean(this_dB_p_v1ref(:,this_band),2);
                            end
                            
                            %Calculate whether this this_dB_p_v1 is significantly
                            %different from reference
                            
                            switch stat_method
                                case 1
                                    %Ranksum
                                    for bwii=1:no_bandwidths
                                        
                                        
                                        p_vals(grNo,per_bin,bwii)=ranksum(this_dB_p_bw(:,bwii),this_dB_p_bwref(:,bwii));
                                        these_p_vals=[these_p_vals ranksum(this_dB_p_bw(:,bwii),this_dB_p_bwref(:,bwii))];
                                        fprintf(1, ['p value for ' freq_names{bwii} ', ' handles_drgb.drgbchoices.group_no_names{grNo} '= %d\n'],  p_vals(grNo,per_bin,bwii));
                                    end
                                case 2
                                    %statcond
                                    a={ this_dB_p_bw' this_dB_p_bwref'};
                                    [F df pvals_perm2(grNo,per_bin,:)] = statcond(a,'mode',mode_statcond); % perform an unpaired ANOVA
                                    for bwii=1:no_bandwidths
                                        fprintf(1, ['p value for ' freq_names{bwii} ', ' handles_drgb.drgbchoices.group_no_names{grNo} '= %d\n'],  pvals_perm2(grNo,per_bin,bwii));
                                    end
                                    
                                    this_pv=zeros(1,no_bandwidths);
                                    this_pv(:,:)=pvals_perm2(grNo,per_bin,:);
                                    these_p_vals_perm=[these_p_vals_perm this_pv];
                            end
                            
                            fprintf(1, '\n')
                            
                        else
                            calc_pval(grNo)=0;
                        end
                        
                    end
                    
                end
            end
            
        end
        
        
        switch stat_method
            case 1
                %Ranksum
                for bwii=1:no_bandwidths
                    pFDR=drsFDRpval(these_p_vals);
                    fprintf(1, ['pFDR for ' freq_names{bwii} ' = %d\n'],pFDR);
                end
            case 2
                %statcond
                pFDR_perm=drsFDRpval(these_p_vals_perm);
                fprintf(1, 'pFDR perm  = %d\n',pFDR_perm);
        end
        
        
        
        %Now do the boxplot
        for per_bin=no_percent_bins:no_percent_bins
            for ii=1:no_bandwidths
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                hold on
                
                
                xtickl=[];
                kk=0;
                this_dB_p=[];
                these_experiments=[];
                group_no=[];
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    if no_valuesHitCR(eventType1,winNo,grNo,per_bin)>1
                        
                        this_dB_p_v1=zeros(no_valuesHitCR(eventType1,winNo,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                        this_dB_p_v1(:,:)=dB_power_valuesHitCR(eventType1,winNo,grNo,per_bin,1:no_valuesHitCR(eventType1,winNo,grNo,per_bin),:);
                        if subtractRef==1
                            this_dB_p_v1_ref=zeros(no_valuesHitCR(eventTypeRef,winNo,grNo,per_bin),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1_ref(:,:)=dB_power_valuesHitCR(eventTypeRef,winNo,grNo,per_bin,1:no_valuesHitCR(eventType1,refWin,grNo,per_bin),:);
                            this_dB_p_v1=this_dB_p_v1-this_dB_p_v1_ref;
                        end
                        
                        
                        this_band=(frequency>=low_freq(ii))&(frequency<=high_freq(ii));
                        expno=zeros(1,no_valuesHitCR(eventType1,winNo,grNo,per_bin));
                        expno(1,:)=experimentNo(eventType1,winNo,grNo,per_bin,1:no_valuesHitCR(eventType1,winNo,grNo,per_bin));
                        these_experiments=[these_experiments expno];
                        this_dB_p=[this_dB_p mean(this_dB_p_v1(:,this_band),2)'];
                        group_no=[group_no grNo*ones(1,no_valuesHitCR(eventType1,winNo,grNo,per_bin))];
                        
                        kk=kk+1;
                        xtickl{kk}=handles_drgb.drgbchoices.group_no_names{grNo};
                        
                        %Now plot the individual points
                        deltax=0.03;
                        no_these_files=0;
                        for ww=1:length(group_per_file)
                            if sum(ww==these_experiments(group_no==grNo))>0
                                no_these_files=no_these_files+1;
                                expNo(no_these_files)=ww;
                            end
                        end
                        x_shift=-deltax*floor(no_these_files/2);
                        for jj=1:no_these_files
                            if x_shift==0
                                x_shift=x_shift+deltax;
                            end
                            no_lfps=sum((these_experiments==expNo(jj))&(group_no==grNo));
                            x=(kk+x_shift)*ones(1,no_lfps);
                            plot(x,this_dB_p((these_experiments==expNo(jj))&(group_no==grNo)),'.k')
                            x_shift=x_shift+deltax;
                        end
                        
                        
                    end
                    
                end
                
                boxplot(this_dB_p,group_no,'Symbol','')
                xticklabels(xtickl)
                if subtractRef==0
                    title([ freq_names{ii} ' Delta power (dB) for event '  percent_bin_legend{per_bin}])
                else
                    title([ freq_names{ii}  ' Delta delta Power (dB) for event '  percent_bin_legend{per_bin}])
                end
                
                pffft=1;
                
            end
        end
        
        
    case 7
        %ROC analysis for Hit and CR for proficient mice
        fprintf(1, 'ROC analysis for Hit and CR for proficient mice\n\n');
        pFDRroc=drsFDRpval(p_vals_ROC);
        fprintf(1, 'pFDR for auROC: %d\n\n',pFDRroc);
        
        
        these_p_vals_perm=[];
        these_p_vals=[];
        
        for per_bin=no_percent_bins:no_percent_bins
            
            
            these_p_vals=[];
            for grRef=1:max(handles_drgb.drgbchoices.group_no)
                fprintf(1, ['\n\n\np values for all groups compared to group ' handles_drgb.drgbchoices.group_no_names{grRef} '\n'],  grRef);
                
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    if grNo~=grRef
                        auROC=[];
                        auROCref=[];
                        for bwii=1:no_bandwidths
                            
                            auROCii=0;
                            auROCrefii=0;
                            for roc_no=1:no_ROCs
                                if (ROCout(roc_no).timeWindow==winNo)&(ROCout(roc_no).bandwidth==bwii)
                                    
                                    if ROCout(roc_no).groupNo==grNo
                                        auROCii=auROCii+1;
                                        auROC(auROCii,bwii)=ROCout(roc_no).roc.AUC-0.5;
                                    end
                                    if ROCout(roc_no).groupNo==grRef
                                        auROCrefii=auROCrefii+1;
                                        auROCref(auROCrefii,bwii)=ROCout(roc_no).roc.AUC-0.5;
                                    end
                                    
                                end
                            end
                            
                        end
                        
                        %Calculate whether this this_dB_p_v1 is significantly
                        %different from reference
                        
                        switch stat_method
                            case 1
                                
                                for bwii=1:no_bandwidths
                                    p_vals(grNo,per_bin,bwii)=ranksum(auROC(:,bwii),auROCref(:,bwii));
                                    these_p_vals=[these_p_vals ranksum(auROC(:,bwii),auROCref(:,bwii))];
                                    fprintf(1, ['p value for ' freq_names{bwii} ', ' handles_drgb.drgbchoices.group_no_names{grNo} '= %d\n'],  p_vals(grNo,per_bin,bwii));
                                end
                            case 2
                                %statcond
                                a={ auROC' auROCref'};
                                [F df pvals_perm2(grNo,per_bin,:)] = statcond(a,'mode',mode_statcond); % perform an unpaired ANOVA
                                for bwii=1:no_bandwidths
                                    fprintf(1, ['p value for ' freq_names{bwii} ', ' handles_drgb.drgbchoices.group_no_names{grNo} '= %d\n'],  pvals_perm2(grNo,per_bin,bwii));
                                end
                                
                                this_pv=zeros(1,no_bandwidths);
                                this_pv(:,:)=pvals_perm2(grNo,per_bin,:);
                                these_p_vals_perm=[these_p_vals_perm this_pv];
                        end
                        
                        fprintf(1, '\n')
                        
                        
                    end
                    
                    
                end
            end
            
        end
        
        
        switch stat_method
            case 1
                %Ranksum
                for bwii=1:no_bandwidths
                    pFDR=drsFDRpval(these_p_vals);
                    fprintf(1, ['pFDR for ' freq_names{bwii} ' = %d\n'],pFDR);
                end
            case 2
                %statcond
                pFDR_perm=drsFDRpval(these_p_vals_perm);
                fprintf(1, 'pFDR perm  = %d\n',pFDR_perm);
        end
        
        
        %Now do the boxplot of auROCs
        for per_bin=no_percent_bins:no_percent_bins
            for bwii=1:no_bandwidths
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                hold on
                
                
                xtickl=[];
                kk=0;
                this_auROC=[];
                these_experiments=[];
                group_no=[];
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    
                    auROC=[];
                    expno=[];
                    jj=0;
                    for roc_no=1:no_ROCs
                        if (ROCout(roc_no).timeWindow==winNo)&(ROCout(roc_no).bandwidth==bwii)&(ROCout(roc_no).groupNo==grNo)
                            jj=jj+1;
                            auROC(jj)=ROCout(roc_no).roc.AUC-0.5;
                            expno(jj)=ROCout(roc_no).fileNo;
                            auROCp(jj)=ROCout(roc_no).roc.p;
                        end
                    end
                    
                    these_experiments=[these_experiments expno];
                    this_auROC=[this_auROC auROC];
                    group_no=[group_no grNo*ones(1,length(auROC))];
                    
                    kk=kk+1;
                    xtickl{kk}=handles_drgb.drgbchoices.group_no_names{grNo};
                    
                    %Now plot the individual points
                    deltax=0.03;
                    no_these_files=0;
                    for ww=1:length(group_per_file)
                        if sum(ww==these_experiments(group_no==grNo))>0
                            no_these_files=no_these_files+1;
                            expNo(no_these_files)=ww;
                        end
                    end
                    x_shift=-deltax*floor(no_these_files/2);
                    for jj=1:no_these_files
                        if x_shift==0
                            x_shift=x_shift+deltax;
                        end
                        no_lfps=sum((these_experiments==expNo(jj))&(group_no==grNo));
                        x=(kk+x_shift)*ones(1,no_lfps);
                        plot(x,this_auROC((these_experiments==expNo(jj))&(group_no==grNo)),'.k')
                        x_shift=x_shift+deltax;
                    end
                    
                    
                    
                    
                end
                
                boxplot(this_auROC,group_no,'Symbol','')
                xticklabels(xtickl)
                if subtractRef==0
                    title([ freq_names{bwii} ' Delta power (dB) for event '  percent_bin_legend{per_bin}])
                else
                    title([ freq_names{bwii}  ' Delta delta Power (dB) for event '  percent_bin_legend{per_bin}])
                end
                
                pffft=1;
                
            end
        end
        
        %Now do the histograms for auROCs
        for per_bin=no_percent_bins:no_percent_bins
            for bwii=1:no_bandwidths
                
                
                
                
                
                n_cum=0;
                this_legend=[];
                no_significant=[];
                total_no=[];
                for grNo=1:max(handles_drgb.drgbchoices.group_no)
                    
                    
                    auROC=[];
                    expno=[];
                    auROCp=[];
                    jj=0;
                    for roc_no=1:no_ROCs
                        if (ROCout(roc_no).timeWindow==winNo)&(ROCout(roc_no).bandwidth==bwii)&(ROCout(roc_no).groupNo==grNo)
                            jj=jj+1;
                            auROC(jj)=ROCout(roc_no).roc.AUC-0.5;
                            expno(jj)=ROCout(roc_no).fileNo;
                            auROCp(jj)=ROCout(roc_no).roc.p;
                        end
                    end
                    
                    total_no(grNo)=length(auROC);
                    no_significant(grNo)=sum(auROC<=pFDRroc);
                    
                    %                   [f_auROC,x_auROC] = drg_ecdf(auROC);
                    %                   plot(x_auROC,f_auROC,these_lines{n_cum})
                    figNo=figNo+1;
                    try
                        close(figNo)
                    catch
                    end
                    figure(figNo)
                    n_cum=n_cum+1;
                    edges=[-0.5:0.05:0.5];
                    histogram(auROC(auROCp>pFDRroc),edges,'EdgeColor','k','FaceColor','b','FaceAlpha',0.4);
                    hold on
                    histogram(auROC(auROCp<=pFDRroc),edges,'EdgeColor','k','FaceColor','r','FaceAlpha',0.4);
                    legend('auROC not significant','auROC significant')
                    title([ freq_names{bwii} ' auROC cumulative histograms for '  percent_bin_legend{per_bin}...
                        ' mice ' handles_drgb.drgbchoices.group_no_names{grNo}])
                    
                end
                
                %                 this_legend=['legend(' this_legend(1:end-1) ')'];
                %                 eval(this_legend)
                
                
                
            end
        end
        
    case 8
        fprintf(1, 'd prime analysis for Hit, CR and FA for proficient mice\n\n');
        per_bin=3;
        figNo=0;
        
        
        for bwii=1:no_bandwidths
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                
                max_dprime=prctile(d_primeHits,99);
                min_dprime=prctile(d_primeCRs,1);
                edges=[1.1*min_dprime:1.1*(max_dprime-min_dprime)/100:1.1*max_dprime];
                
                %Plot Hit histogram
                these_dprimeHits=d_primeHits((grHits==grNo)&(bwii==bwHits));
                histogram(these_dprimeHits,edges,'EdgeColor','k','FaceColor','r','FaceAlpha',0.4);
                
                hold on
                
                %Plot CR histogram
                these_dprimeCRs=d_primeCRs((grCRs==grNo)&(bwii==bwCRs));
                histogram(these_dprimeCRs,edges,'EdgeColor','k','FaceColor','b','FaceAlpha',0.4);
                
                %Plot FA histogram
                these_dprimeFAs=d_primeFAs((grFAs==grNo)&(bwii==bwFAs));
                histogram(these_dprimeFAs,edges,'EdgeColor','k','FaceColor','g','FaceAlpha',0.4);
                
                legend('Hit','CR','FA')
                
                title([ freq_names{bwii} ' d pime histograms for '  percent_bin_legend{per_bin}...
                    ' mice ' handles_drgb.drgbchoices.group_no_names{grNo}])
                
                
            end
            
        end
        
        %Now plot the cumulative histograms for each event
        
        
        p_valsKS=[];
        
        %Hit
        for bwii=1:no_bandwidths
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            hold on
            
            this_legend=[];
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                [f_dprime,x_dprime] = drg_ecdf(d_primeHits((grHits==grNo)&(bwii==bwHits)));
                plot(x_dprime,f_dprime,these_colors{grNo})
                this_legend=[this_legend '''' handles_drgb.drgbchoices.group_no_names{grNo} '''' ','];
            end
            
            this_legend=['legend(' this_legend(1:end-1) ')'];
            eval(this_legend)
            title([ freq_names{bwii} ' d pime cumulative histograms for Hit for '  percent_bin_legend{per_bin}...
                ' mice '])
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                for grRef=1:max(handles_drgb.drgbchoices.group_no)
                    if grRef~=grNo
                        %[h,p_val]=kstest2(d_primeHits((grHits==grNo)&(bwii==bwHits)),d_primeHits((grHits==grRef)&(bwii==bwHits)));
                        [p_val, h]=ranksum(d_primeHits((grHits==grNo)&(bwii==bwHits)),d_primeHits((grHits==grRef)&(bwii==bwHits)));
                        p_valsKS=[p_valsKS p_val];
                        
                        fprintf(1, ['p value for '  freq_names{bwii} ' Hit d prime ' ', ' handles_drgb.drgbchoices.group_no_names{grNo}...
                            ' vs. ' handles_drgb.drgbchoices.group_no_names{grRef} ' = %d\n'],  p_val);
                    end
                end
            end
            fprintf(1, '\n')
        end
        
        fprintf(1, '\n\n')
        %CR
        for bwii=1:no_bandwidths
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            hold on
            
            this_legend=[];
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                [f_dprime,x_dprime] = drg_ecdf(d_primeCRs((grCRs==grNo)&(bwii==bwCRs)));
                plot(x_dprime,f_dprime,these_colors{grNo})
                this_legend=[this_legend '''' handles_drgb.drgbchoices.group_no_names{grNo} '''' ','];
            end
            
            this_legend=['legend(' this_legend(1:end-1) ')'];
            eval(this_legend)
            title([ freq_names{bwii} ' d pime cumulative histograms for CR for '  percent_bin_legend{per_bin}...
                ' mice '])
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                for grRef=1:max(handles_drgb.drgbchoices.group_no)
                    if grRef~=grNo
                        %[h,p_val]=kstest2(d_primeCRs((grCRs==grNo)&(bwii==bwCRs)),d_primeCRs((grCRs==grRef)&(bwii==bwCRs)));
                        [p_val, h]=ranksum(d_primeCRs((grCRs==grNo)&(bwii==bwCRs)),d_primeCRs((grCRs==grRef)&(bwii==bwCRs)));
                        p_valsKS=[p_valsKS p_val];
                        
                        fprintf(1, ['p value for '  freq_names{bwii} ' CR d prime ' ', ' handles_drgb.drgbchoices.group_no_names{grNo}...
                            ' vs. ' handles_drgb.drgbchoices.group_no_names{grRef} ' = %d\n'],  p_val);
                    end
                end
            end
            fprintf(1, '\n')
        end
        fprintf(1, '\n\n')
        
        %FA
        for bwii=1:no_bandwidths
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            hold on
            
            this_legend=[];
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                [f_dprime,x_dprime] = drg_ecdf(d_primeFAs((grFAs==grNo)&(bwii==bwFAs)));
                plot(x_dprime,f_dprime,these_colors{grNo})
                this_legend=[this_legend '''' handles_drgb.drgbchoices.group_no_names{grNo} '''' ','];
            end
            
            this_legend=['legend(' this_legend(1:end-1) ')'];
            eval(this_legend)
            title([ freq_names{bwii} ' d pime cumulative histograms for FA for '  percent_bin_legend{per_bin}...
                ' mice '])
            
            for grNo=1:max(handles_drgb.drgbchoices.group_no)
                for grRef=1:max(handles_drgb.drgbchoices.group_no)
                    if grRef~=grNo
                        %[h,p_val]=kstest2(d_primeFAs((grFAs==grNo)&(bwii==bwFAs)),d_primeFAs((grFAs==grRef)&(bwii==bwFAs)));
                        [p_val, h]=ranksum(d_primeFAs((grFAs==grNo)&(bwii==bwFAs)),d_primeFAs((grFAs==grRef)&(bwii==bwFAs)));
                        p_valsKS=[p_valsKS p_val];
                        
                        fprintf(1, ['p value for '  freq_names{bwii} ' FA d prime ' ', ' handles_drgb.drgbchoices.group_no_names{grNo}...
                            ' vs. ' handles_drgb.drgbchoices.group_no_names{grRef} ' = %d\n'],  p_val);
                    end
                end
            end
            fprintf(1, '\n')
        end
        fprintf(1, '\n\n')
        pFDRdprime=drsFDRpval(p_valsKS)
        
    case 9
        %d prime analysis for Fig. 3 of Daniel's pub
        fprintf(1, 'd prime analysis for Hit, CR and FA for proficient mice\n\n');
        per_bin=3;
        figNo=0;
        
        
        for bwii=1:no_bandwidths
            
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            max_dprime=prctile(d_primeHits,99);
            min_dprime=prctile(d_primeCRs,1);
            edges=[1.1*min_dprime:1.1*(max_dprime-min_dprime)/100:1.1*max_dprime];
            
            %Plot Hit histogram
            these_dprimeHits=d_primeHits(((grHits==1)|(grHits==3))&(bwii==bwHits));
            histogram(these_dprimeHits,edges,'EdgeColor','k','FaceColor','r','FaceAlpha',0.4);
            
            hold on
            
            %Plot CR histogram
            these_dprimeCRs=d_primeCRs(((grCRs==1)|(grCRs==3))&(bwii==bwCRs));
            histogram(these_dprimeCRs,edges,'EdgeColor','k','FaceColor','b','FaceAlpha',0.4);
            
            %Plot FA histogram
            these_dprimeFAs=d_primeFAs(((grFAs==1)|(grFAs==3))&(bwii==bwFAs));
            histogram(these_dprimeFAs,edges,'EdgeColor','k','FaceColor','g','FaceAlpha',0.4);
            
            legend('Hit','CR','FA')
            
            title([ freq_names{bwii} ' d pime histograms for '  percent_bin_legend{per_bin}])
            
            %Plot ROC for HIt vs CR
            roc_data=[];
            %Enter d prime for Hits
            roc_data(1:length(these_dprimeHits),1)=these_dprimeHits';
            roc_data(1:length(these_dprimeHits),2)=zeros(length(these_dprimeHits),1);
            
            %Enter dprime for CR
            roc_data(length(these_dprimeHits)+1:length(these_dprimeHits)+length(these_dprimeCRs),1)=these_dprimeCRs';
            roc_data(length(these_dprimeHits)+1:length(these_dprimeHits)+length(these_dprimeCRs),2)=ones(length(these_dprimeCRs),1);
            
            %Find pre ROC
            rocHitCR=roc_calc(roc_data,0,0.05,0);
            
            fprintf(1, ['auROC for '  freq_names{bwii} ' Hit vs. CR = %d\n'],  rocHitCR.AUC-0.5);
            
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            
            hold on
            plot([0 1],[0 1],'k');
            plot(rocHitCR.xr,rocHitCR.yr,'ro','markersize',1,'markeredgecolor','r','markerfacecolor','r');
            
            
            
            %Plot ROC for CR vs FA
            roc_data=[];
            %Enter d prime for CRs
            roc_data(1:length(these_dprimeCRs),1)=these_dprimeCRs';
            roc_data(1:length(these_dprimeCRs),2)=zeros(length(these_dprimeCRs),1);
            
            %Enter dprime for FA
            roc_data(length(these_dprimeCRs)+1:length(these_dprimeCRs)+length(these_dprimeFAs),1)=these_dprimeFAs';
            roc_data(length(these_dprimeCRs)+1:length(these_dprimeCRs)+length(these_dprimeFAs),2)=ones(length(these_dprimeFAs),1);
            
            
            %Find pre ROC
            rocCRFA=roc_calc(roc_data,0,0.05,0);
            fprintf(1, ['auROC for '  freq_names{bwii} ' FA vs. CR = %d\n'],  rocCRFA.AUC-0.5);
            
            plot(rocCRFA.xr,rocCRFA.yr,'bo','markersize',1,'markeredgecolor','b','markerfacecolor','b');
            
            xlabel('False positive rate (1-Specificity)')
            ylabel('True positive rate (Sensitivity)')
            title(['ROC curve for Hit or FA vs CR for ' freq_names{bwii}])
            legend('Hit','FA')
            
            
            
        end
        
        
end
function drgDisplayBatchLFPPowerPerExp(handles)

%This function displays the LFP power spectrum for drgRunBatch-generated
%data


%Which time window do you want displayed
winNo=2;

%Do you wnat to subtract the LFP spectrum in a reference window?
subtractRef=1;
refWin=1;

close all
warning('off')


%Choose the format for the display and the event types 
%THESE VALUES ARE IMPORTANT


% which_display
%
% 1 Show the difference between events for each group/ percent correct
%  bin in a separate graph
%
% 2 Show the difference between groups for each event/percent correct
%  bin in a separate graph
%
% 3 Show the difference between learning and proficient
%


%For Alexia's tstart learning vs. proficeint
% which_display=3;
% eventType=1; %tstart
% evTypeLabels={'tstart'};

% For Alexia's odorOn (CS) learning vs. proficeint
% which_display=3;
% eventType=1; %OdorOn
% evTypeLabels={'CS'};


% For Alexia's Hit, CR, FA
% which_display=1;
% eventType=[2 5 7];
% evTypeLabels={'Hit','CR','FA'};

% For Daniel's Hit, CR
which_display=2;
eventType=[2 5];
evTypeLabels={'Hit','CR'};

%THESE VALUES ARE IMPORTANT
%VERY IMPORTANT: This is the index for this event in your handles.drgbchoices.evTypeNos
% eventType=[2 5]; %Hit and CR
% evTypeLabels={'Hit';'CR'};






% eventType=[9 10 11 12 13 14]; %Hit and CR
% no_event_types=6;
% evTypeLabels={'Hi Od1';'Hi Od 2';'Hi Od 3';'Low Od1';'Low Od 2';'Low Od 3'};

no_event_types=length(eventType);

%Which percent correct bins do you want to use?
percent_low=[45 65 80];
percent_high=[65 80 100];
percent_bin_legend={' 45-65%';' 65-80%';' 80-100%'};
no_percent_bins=3;

%Ask user for the drgb output .mat file and load those data
[handles.drgb.outFileName,handles.PathName] = uigetfile('*.mat','Select the drgb output file');
load([handles.PathName handles.drgb.outFileName])


%Initialize the variables

%Determine how many LFPs are included for each evTypeNo, time window, group, etc
no_groups=max(handles_drgb.drgbchoices.group_no);
no_events=length(handles_drgb.drgbchoices.evTypeNos);
no_lfpevpairs=handles_drgb.drgb.lfp_per_exp_no;
no_lfps=zeros(length(handles_drgb.drgbchoices.evTypeNos),handles_drgb.drgbchoices.noWindows,max(handles_drgb.drgbchoices.group_no),no_percent_bins);

files=[];
noFiles=0;

for lfpodNo=1:no_lfpevpairs
    fileNo=handles_drgb.drgb.lfp_per_exp(lfpodNo).fileNo;
    noFiles=noFiles+1;
    files(noFiles)=fileNo;
    groupNo=handles_drgb.drgbchoices.group_no(fileNo);
    timeWindow=handles_drgb.drgb.lfp_per_exp(lfpodNo).timeWindow;
    for evTypeNo=1:no_events
        
        try
            trials_in_event=handles_drgb.drgb.lfp_per_exp(lfpodNo).which_eventLFPPower(evTypeNo,:)==1;
        catch
            trials_in_event=[];
        end
        if sum(trials_in_event)>0
            for percent_bin=1:no_percent_bins
                trials_in_perbin=(handles_drgb.drgb.lfp_per_exp(lfpodNo).perCorrLFPPower>percent_low(percent_bin))&(handles_drgb.drgb.lfp_per_exp(lfpodNo).perCorrLFPPower<=percent_high(percent_bin));
                
                if sum(trials_in_event&trials_in_perbin)>0
                    no_lfps(evTypeNo,timeWindow,groupNo,percent_bin)=no_lfps(evTypeNo,timeWindow,groupNo,percent_bin)+1;
                end
            end
        end
    end
end
max_lfps=max(no_lfps(:));

%Now initialize the variables for db power
noWindows=handles_drgb.drgbchoices.noWindows;
dB_power_values=zeros(no_events,noWindows,no_groups,no_percent_bins,max_lfps,length(handles_drgb.drgb.freq_for_LFPpower));
experimentNo=zeros(no_events,noWindows,no_groups,no_percent_bins,max_lfps);
no_values=zeros(no_events,noWindows,no_groups,no_percent_bins);


%Sort out all power values into a matrix
for lfpodNo=1:no_lfpevpairs
    fileNo=handles_drgb.drgb.lfp_per_exp(lfpodNo).fileNo;
    groupNo=handles_drgb.drgbchoices.group_no(fileNo);
    timeWindow=handles_drgb.drgb.lfp_per_exp(lfpodNo).timeWindow;
    lfpodorNo=lfpodNo;
    for evTypeNo=1:no_events
        
        try
            trials_in_event=handles_drgb.drgb.lfp_per_exp(lfpodNo).which_eventLFPPower(evTypeNo,:)==1;
        catch
            trials_in_event=[];
        end
        
        if sum(trials_in_event)>0
            for percent_bin=1:no_percent_bins
                trials_in_perbin=(handles_drgb.drgb.lfp_per_exp(lfpodNo).perCorrLFPPower>percent_low(percent_bin))&(handles_drgb.drgb.lfp_per_exp(lfpodNo).perCorrLFPPower<=percent_high(percent_bin));
                
                if sum(trials_in_event&trials_in_perbin)>0
                    
                    no_trials=sum(trials_in_event&trials_in_perbin);
                    
                    this_dB_power=zeros(no_trials,length(handles_drgb.drgb.freq_for_LFPpower));
                    this_dB_power(:,:)=10*log10(handles_drgb.drgb.lfp_per_exp(lfpodNo).allPower(:,trials_in_event&trials_in_perbin))';
                    dB_power_values(evTypeNo,timeWindow,groupNo,percent_bin,no_values(evTypeNo,timeWindow,groupNo,percent_bin)+1,:)=...
                        mean(this_dB_power,1);
                    experimentNo(evTypeNo,timeWindow,groupNo,percent_bin,no_values(evTypeNo,timeWindow,groupNo,percent_bin)+1)=handles_drgb.drgb.lfp_per_exp(lfpodNo).fileNo;
                    no_values(evTypeNo,timeWindow,groupNo,percent_bin)=no_values(evTypeNo,timeWindow,groupNo,percent_bin)+1;

                    
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
                            dB_p_v1_ci(per_bin,evTN1,ifreq)=ci(1,1)-pd.mu;
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
        
        for evTN1=1:length(eventType)
            eventType1=eventType(evTN1);
            
            dB_p_v1_mean=zeros(no_percent_bins, max(handles_drgb.drgbchoices.group_no),length(handles_drgb.drgb.freq_for_LFPpower));
            dB_p_v1_SEM=zeros(no_percent_bins, max(handles_drgb.drgbchoices.group_no),length(handles_drgb.drgb.freq_for_LFPpower));
            
            
            frequency=handles_drgb.drgb.freq_for_LFPpower;
            
            
            kk=0;
            for per_bin=1:no_percent_bins
                
                grRef=1;  %p values will be calculated with respect to this group
                these_p_vals=[];
                
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
                                    p_vals(grNo,per_bin,evTN1,ifreq)=ranksum(this_dB_p_v1(:,ifreq),this_dB_p_v1ref(:,ifreq));
                                    these_p_vals=[these_p_vals ranksum(this_dB_p_v1(:,ifreq),this_dB_p_v1ref(:,ifreq))];
                                end
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
                
                pFDR(evTN1,per_bin)=drsFDRpval(these_p_vals);
                fprintf(1, 'pFDR for event type: %d, percent bin no: %d = %d\n',evTN1,per_bin,pFDR(evTN1,per_bin));
                
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
                        dB_p_v1_mean(per_bin,grNo,:)=mean(this_dB_p_v1,1);
                        dB_p_v1_SEM(per_bin,grNo,:)=std(this_dB_p_v1,0,1)/sqrt(no_values(eventType1,winNo,grNo,per_bin));
                        %Calculate the 95% CI
                        for ifreq=1:length(frequency)
                            pd=fitdist(this_dB_p_v1(:,ifreq),'Normal');
                            ci=paramci(pd);
                            dB_p_v1_ci(per_bin,grNo,ifreq)=ci(1,1)-pd.mu;
                        end
                        these_experiments=experimentNo(eventType1,winNo,grNo,per_bin,1:no_values(eventType1,winNo,grNo,per_bin));
                        
                    end
                end
                
                
                
            end
            
            
            frequency=handles_drgb.drgb.freq_for_LFPpower;
            
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
                    this_dB_p_v1_mean(1,:)=dB_p_v1_mean(per_bin,grNo,:);
                    
                    this_dB_p_v1_ci=zeros(1,length(handles_drgb.drgb.freq_for_LFPpower));
                    this_dB_p_v1_ci(1,:)=dB_p_v1_ci(per_bin,grNo,:);
                    
                    %plot(frequency,this_dB_p_v1_mean, these_lines{grNo});
                    eval(['p' num2str(grNo) '=plot(frequency,this_dB_p_v1_mean, these_lines{grNo});'])
                    
                    
                    this_legend=[this_legend '''' handles_drgb.drgbchoices.group_no_names{grNo} '''' ' ,'];
                    plot_handles=[plot_handles 'p' num2str(grNo) ' '];
                    
                end
                
                this_legend=['legend([' plot_handles(1:end-1) '],' this_legend(1:end-1) ')'];
                
              
                
                
                %Now plot the bounded lines
                
                for grNo=max(handles_drgb.drgbchoices.group_no):-1:1
                    
                    
                    this_dB_p_v1_mean=zeros(1,length(handles_drgb.drgb.freq_for_LFPpower));
                    this_dB_p_v1_mean(1,:)=dB_p_v1_mean(per_bin,grNo,:);
                    
                    this_dB_p_v1_ci=zeros(1,length(handles_drgb.drgb.freq_for_LFPpower));
                    this_dB_p_v1_ci(1,:)=dB_p_v1_ci(per_bin,grNo,:);
                    
                    [hl, hp] = boundedline(frequency,this_dB_p_v1_mean, this_dB_p_v1_ci, these_colors{grNo});
                    
                    if subtractRef==1
                        if calc_pval(grNo)==1
                            for ifreq=1:length(frequency)
                                
                                if p_vals(grNo,per_bin,evTN1,ifreq)<=pFDR(evTN1,per_bin)
                                    plot(frequency(ifreq),this_dB_p_v1_mean(ifreq),these_circles{grNo})
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
        %Display the difference between learning and proficient
        for grNo=1:max(handles_drgb.drgbchoices.group_no)
            dB_p_v1_mean=zeros(no_percent_bins, no_event_types,length(handles_drgb.drgb.freq_for_LFPpower));
            dB_p_v1_SEM=zeros(no_percent_bins, no_event_types,length(handles_drgb.drgb.freq_for_LFPpower));
            
            
            frequency=handles_drgb.drgb.freq_for_LFPpower;
            
            
            kk=0;
            
            
            
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
                            
                            this_dB_p_v1ref=zeros(no_values(eventType1,winNo,grNo,per_bin_ref),length(handles_drgb.drgb.freq_for_LFPpower));
                            this_dB_p_v1ref(:,:)=dB_power_values(eventType1,winNo,grNo,per_bin_ref,1:no_values(eventType1,winNo,grNo,per_bin_ref),:);
                            if subtractRef==1
                                this_dB_p_v1ref_ref=zeros(no_values(eventType1,refWin,grNo,per_bin_ref),length(handles_drgb.drgb.freq_for_LFPpower));
                                this_dB_p_v1ref_ref(:,:)=dB_power_values(eventType1,refWin,grNo,per_bin_ref,1:no_values(eventType1,refWin,grNo,per_bin_ref),:);
                                this_dB_p_v1ref=this_dB_p_v1ref-this_dB_p_v1ref_ref;
                            end
                            
                            %Calculate whether this this_dB_p_v1 is significantly
                            %different from reference
                            for ifreq=1:length(frequency)
                                p_vals(grNo,per_bin,evTN1,ifreq)=ranksum(this_dB_p_v1(:,ifreq),this_dB_p_v1ref(:,ifreq));
                                these_p_vals=[these_p_vals ranksum(this_dB_p_v1(:,ifreq),this_dB_p_v1ref(:,ifreq))];
                            end
                            
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
                            dB_p_v1_ci(per_bin,evTN1,ifreq)=ci(1,1)-pd.mu;
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
                            
                            
                            if p_vals(grNo,per_bin,evTN1,ifreq)<=pFDR(grNo,evTN1)
                                plot(frequency(ifreq),this_dB_p_v1_mean(ifreq),these_circles{per_bin})
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
        
end


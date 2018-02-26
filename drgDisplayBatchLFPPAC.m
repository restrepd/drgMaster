function drgDisplayBatchLFPPAC(handles)

%This function displays the results of PAC analysis for drgRunBatch-generated
%data

close all
clear all
warning('off')

%THESE VALUES ARE IMPORTANT

%Choose the format for the display
which_PA_display=1;

% 1 Show the polar histograms for peak angles for the difference between events
%
% 2 Show the polar histograms for peak angles for the difference between
% groups
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
gr_to_test_against=1;



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
% no_lfpevpairs=682;
no_PACpeaks=handles_drgb.no_PACpeaks;
no_trials=zeros(length(handles_drgb.drgbchoices.evTypeNos),handles_drgb.drgbchoices.noWindows,max(handles_drgb.drgbchoices.group_no),no_percent_bins,handles_drgb.no_PACpeaks);

%In the previous versions of the code these variables were not defined
if isfield(handles_drgb.drgb.lfpevpair,'which_eventPAC')==0
    for lfpodNo=1:no_lfpevpairs
        handles_drgb.drgb.lfpevpair(lfpodNo).which_eventPAC=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).which_eventPAC;
        handles_drgb.drgb.lfpevpair(lfpodNo).perCorrPAC=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(1).perCorrPAC;
    end
end
   
no_empty=0;
for lfpodNo=1:no_lfpevpairs
    fileNo=handles_drgb.drgb.lfpevpair(lfpodNo).fileNo;
    groupNo=handles_drgb.drgbchoices.group_no(fileNo);
    timeWindow=handles_drgb.drgb.lfpevpair(lfpodNo).timeWindow;
    for evTypeNo=1:no_events
        for PACno=1:no_PACpeaks
            if ~isempty(handles_drgb.drgb.lfpevpair(lfpodNo).PAC(PACno).which_eventPAC)
                trials_in_event=handles_drgb.drgb.lfpevpair(lfpodNo).PAC(PACno).which_eventPAC(evTypeNo,:)==1;
                if sum(trials_in_event)>0
                    for percent_bin=1:no_percent_bins
                        trials_in_perbin=(handles_drgb.drgb.lfpevpair(lfpodNo).PAC(PACno).perCorrPAC>percent_low(percent_bin))&(handles_drgb.drgb.lfpevpair(lfpodNo).perCorrPAC<=percent_high(percent_bin));
                        if sum(trials_in_event&trials_in_perbin)>0
                            no_trials(evTypeNo,timeWindow,groupNo,percent_bin,PACno)=no_trials(evTypeNo,timeWindow,groupNo,percent_bin,PACno)+sum(trials_in_event&trials_in_perbin);
                        end
                    end
                end
            else
                no_empty=no_empty+1;
            end
        end
    end
end
max_trials=max(no_trials(:));
fprintf(1,'Of %d LFP event pairs %d were empty\n',no_lfpevpairs,no_empty)
 
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
        if ~isempty(handles_drgb.drgb.lfpevpair(lfpodNo).which_eventPAC)
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
                    
                    
                    this_legend=[];
                    
                    n_cum=0;
                    
                    for evTN=1:length(evTypeNos)
                        eventType=evTypeNos(evTN);
                        
                        if no_values(eventType,winNo,grNo,per_bin,PACno)>1
                            pa_v=zeros(1,no_values(eventType,winNo,grNo,per_bin,PACno));
                            pa_v(1,1:no_values(eventType,winNo,grNo,per_bin,PACno))=...
                                pa_values(eventType,winNo,grNo,per_bin,PACno,1:no_values(eventType,winNo,grNo,per_bin,PACno));
                            n_cum=n_cum+1;
                            polarhistogram(pa_v,36,'FaceColor',these_bars{n_cum},'FaceAlpha',0.3)
                            hold on
                            this_legend=[this_legend '''' evTypeLabels{evTN} '''' ','];
                        end
                    end
                    this_legend=['legend(' this_legend(1:end-1) ')'];
                    eval(this_legend)
                    title(['Peak phase angle for ' handles_drgb.PACnames{PACno} '. Percent correct from ' num2str(percent_low(per_bin))...
                        ' to ' num2str(percent_high(per_bin)) ', group No: ' handles_drgb.drgbchoices.group_no_names{grNo}])
                    
                    %Do the circ_wwtest test
                    evTN=ev_to_test_against;
                    eventType=evTypeNos(evTN);
                    p_vals=[];
                    if no_values(eventType,winNo,grNo,per_bin,PACno)>1
                        pa_ref=zeros(1,no_values(eventType,winNo,grNo,per_bin,PACno));
                        pa_ref(1,1:no_values(eventType,winNo,grNo,per_bin,PACno))=...
                            pa_values(eventType,winNo,grNo,per_bin,PACno,1:no_values(eventType,winNo,grNo,per_bin,PACno));
                        fprintf(1,['circ_wwtest p values for peak phase angle cumulative probability for ' handles_drgb.PACnames{PACno} '. Percent correct from ' num2str(percent_low(per_bin))...
                            ' to ' num2str(percent_high(per_bin)) ', group No: ' handles_drgb.drgbchoices.group_no_names{grNo} '\n'])
                        
                        for evTN=1:length(evTypeNos)
                            if evTN~=ev_to_test_against
                                eventType=evTypeNos(evTN);
                                if no_values(eventType,winNo,grNo,per_bin,PACno)>1
                                    pa_v=zeros(1,no_values(eventType,winNo,grNo,per_bin,PACno));
                                    pa_v(1,1:no_values(eventType,winNo,grNo,per_bin,PACno))=...
                                        pa_values(eventType,winNo,grNo,per_bin,PACno,1:no_values(eventType,winNo,grNo,per_bin,PACno));
                                    p=circ_wwtest(pa_v,pa_ref);
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
        %Now do the polar histograms for the difference between groups
        %Cumulative histogram for the peak angles
        figNo=0;
         p_vals=[];
        for PACno=1:handles_drgb.no_PACpeaks
            for per_bin=1:no_percent_bins
                for evTN=1:length(evTypeNos)
                    eventType=evTypeNos(evTN);
                    figNo=figNo+1;
                    try
                        close(figNo)
                    catch
                    end
                    hFig=figure(figNo);
                    
                    
                    this_legend=[];
                    
                    n_cum=0;
                    
                    for grNo=1:no_groups
                        
                        if no_values(eventType,winNo,grNo,per_bin,PACno)>1
                            pa_v=zeros(1,no_values(eventType,winNo,grNo,per_bin,PACno));
                            pa_v(1,1:no_values(eventType,winNo,grNo,per_bin,PACno))=...
                                pa_values(eventType,winNo,grNo,per_bin,PACno,1:no_values(eventType,winNo,grNo,per_bin,PACno));
                            n_cum=n_cum+1;
                            polarhistogram(pa_v,18,'FaceColor',these_bars{n_cum},'FaceAlpha',0.3)
                            hold on
                            this_legend=[this_legend '''' handles_drgb.drgbchoices.group_no_names{grNo} '''' ','];
                        end
                    end
                    this_legend=['legend(' this_legend(1:end-1) ')'];
                    eval(this_legend)
                    title(['Peak phase angle for ' handles_drgb.PACnames{PACno} '. Percent correct from ' num2str(percent_low(per_bin))...
                        ' to ' num2str(percent_high(per_bin)) ', event: ' evTypeLabels{evTN}])
                    
                    %Do the circ_wwtest test
                    
                   
                    if no_values(eventType,winNo,gr_to_test_against,per_bin,PACno)>1
                        pa_ref=zeros(1,no_values(eventType,winNo,gr_to_test_against,per_bin,PACno));
                        pa_ref(1,1:no_values(eventType,winNo,gr_to_test_against,per_bin,PACno))=...
                            pa_values(eventType,winNo,gr_to_test_against,per_bin,PACno,1:no_values(eventType,winNo,gr_to_test_against,per_bin,PACno));
                        fprintf(1,['circ_wwtest p values for peak phase angle cumulative probability for ' handles_drgb.PACnames{PACno} '\n. Percent correct from ' num2str(percent_low(per_bin))...
                            ' to ' num2str(percent_high(per_bin)) ', event: ' evTypeLabels{evTN} '\n'])
                        
                        for grNo=1:no_groups
                            if grNo~=gr_to_test_against
                                if no_values(eventType,winNo,grNo,per_bin,PACno)>1
                                    pa_v=zeros(1,no_values(eventType,winNo,grNo,per_bin,PACno));
                                    pa_v(1,1:no_values(eventType,winNo,grNo,per_bin,PACno))=...
                                        pa_values(eventType,winNo,grNo,per_bin,PACno,1:no_values(eventType,winNo,grNo,per_bin,PACno));
                                    p=circ_wwtest(pa_v,pa_ref);
                                    p_vals=[p_vals p];
                                    fprintf(1, [handles_drgb.drgbchoices.group_no_names{grNo} ' vs. ' handles_drgb.drgbchoices.group_no_names{gr_to_test_against} '= %d\n'],p);
                                end
                            end
                        end
                       
                        fprintf(1, '\n');
                    end
                    
                end
            end
        end
        
        
        if length(p_vals)>1
            fprintf(1, 'pFDR= = %d\n',drsFDRpval(p_vals));
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
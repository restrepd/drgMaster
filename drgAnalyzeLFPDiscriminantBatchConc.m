function drgAnalyzeLFPDiscriminantBatchConc
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'Discriminant_*.mat' output file from drgLFPDiscriminantBatch
%Performs an analysis of the timecourse for percent correct for LDA and for
%the PCA

% which_display chooses the analysis:
%
%1 Displays average predicton for proficeint vs naive for LDA and PCA for power LFP
%
%2 Displays average predicton for proficeint vs naive for LDA and PCA for angle in PAC
%
%3 Displays average prediction and dimensionality for peak and trough for LDA for wavelet
%power referenced to the phase of PAC and plots PC1 for the PCA.
%These are choices 10 and 11 in drgLFPDiscriminantBatch

warning('off')
close all
clear all



t_odor_arrival=0.1;

which_display=3;
mice_excluded=[];

[fname,pname,nCancel] = uigetfile({'Discriminant_*.mat'},'Select the perceptron LFP batch output file ...');
if nCancel
    inputPath = [pname,fname];
    pnameStart = pname;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end

discriminant_name=[pname fname];
load(discriminant_name)

t_pac=handles_out.t_power';
t_power=handles_out.t_power';
handles.drg=handles_out.drg;

figNo=0;



%Define the windows for analysis
window_start=[-1 1.5];
window_end=[0 2.5];
no_wins=2;

%This is the window for area under the curve case 3
auc_from=0.1;
auc_to=2.5;

switch which_display
    
    case 1
        
        %Plot average percent correct for the LFP power LDA
        t=handles_out.t;
        for groupNo=1:max(handles_out.drgbchoices.group_no)
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            
            set(hFig, 'units','normalized','position',[.2 .2 .7 .7])
        end
        
        ii_plot=0;
        for bwii=1:4
            p_correct_stats=[];
            ii_stats=0;
            for percent_correct_ii=1:2
                ii_plot=bwii+(percent_correct_ii-1)*4;
                for groupNo=1:max(handles_out.drgbchoices.group_no)
                    figure(groupNo)
                    subplot(2,4,ii_plot);
                    hold on
                    no_mice=0;
                    all_discriminant_correct=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    all_discriminant_correct_shuffled=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    switch which_display
                        case 1
                            for mouseNo=1:max(handles_out.drgbchoices.mouse_no)
                                
                                if isfield(handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii),'discriminant_calulated')
                                    %This is here because there was a spelling mistake
                                    %"calulated"
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=...
                                        handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calulated;
                                end
                                if handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                                    per_ii=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).no_trials;
                                    if per_ii>=20
                                        no_mice=no_mice+1;
                                        all_discriminant_correct(no_mice,:)=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).discriminant_correct;
                                        all_discriminant_correct_shuffled(no_mice,:)=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).discriminant_correct_shuffled;
                                    end
                                end
                            end
                        case 7
                            for mouseNo=1:max(handles_out.drgbchoices.mouse_no)
                                
                                if isfield(handles_out.discriminant_PCApower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii),'discriminant_calulated')
                                    %This is here because there was a spelling mistake
                                    %"calulated"
                                    handles_out.discriminant_PCApower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=...
                                        handles_out.discriminant_PCApower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calulated;
                                end
                                if handles_out.discriminant_PCApower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                                    per_ii=handles_out.discriminant_PCApower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).no_trials;
                                    if per_ii>=20
                                        no_mice=no_mice+1;
                                        all_discriminant_correct(no_mice,:)=handles_out.discriminant_PCApower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).discriminant_correct_peak;
                                        all_discriminant_correct_shuffled(no_mice,:)=handles_out.discriminant_PCApower_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).discriminant_correct_shuffled_peak;
                                    end
                                end
                            end
                    end
                    
                    no_mice_per(percent_correct_ii)=no_mice;
                    
                    all_discriminant_correct_shuffled=all_discriminant_correct_shuffled(1:no_mice,:);
                    mean_dcsh=mean(all_discriminant_correct_shuffled,1)';
                    CIdcsh = bootci(1000, {@mean, all_discriminant_correct_shuffled})';
                    CIdcsh(:,1)=mean_dcsh-CIdcsh(:,1);
                    CIdcsh(:,2)=CIdcsh(:,2)-mean_dcsh;
                    [hlCR, hpCR] = boundedline(t,mean_dcsh, CIdcsh, 'b');
                    
                    data=[];
                    for mouseNo=1:no_mice
                        for winNo=1:no_wins
                            data=[data mean(all_discriminant_correct_shuffled(mouseNo,(t>=window_start(winNo))&(t<=window_end(winNo))),2)];
                        end
                    end
                    
                    ii_stats=ii_stats+1;
                    p_correct_stats(ii_stats).data=data;
                    p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' shuffled'];
                    p_correct_stats(ii_stats).winNo=winNo; %Note that I use the shuffled data in both windows
                    p_correct_stats(ii_stats).bwii=bwii;
                    p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                    p_correct_stats(ii_stats).groupNo=groupNo;
                    
                    all_discriminant_correct=all_discriminant_correct(1:no_mice,:);
                    mean_dc=mean(all_discriminant_correct,1)';
                    CIdc = bootci(1000, {@mean, all_discriminant_correct})';
                    CIdc(:,1)=mean_dc-CIdc(:,1);
                    CIdc(:,2)=CIdc(:,2)-mean_dc;
                    [hlCR, hpCR] = boundedline(t,mean_dc, CIdc, 'r');
                    
                    for winNo=1:no_wins
                        data=[];
                        for mouseNo=1:no_mice
                            data=[data mean(all_discriminant_correct(mouseNo,(t>=window_start(winNo))&(t<=window_end(winNo))),2)];
                        end
                        ii_stats=ii_stats+1;
                        p_correct_stats(ii_stats).data=data;
                        p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                            handles_out.drgbchoices.per_lab{percent_correct_ii}...
                            ' from ' num2str(window_start(winNo)) ' to ' num2str(window_end(winNo)) ' sec'];
                        p_correct_stats(ii_stats).winNo=winNo;
                        p_correct_stats(ii_stats).bwii=bwii;
                        p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                        p_correct_stats(ii_stats).groupNo=groupNo;
                    end
                    
                    %Odor on markers
                    plot([0 0],[0 100],'-k')
                    odorhl=plot([0 2.5],[10 10],'-k','LineWidth',5);
                    plot([2.5 2.5],[0 100],'-k')
                    
                    title([handles_out.drgbchoices.bwlabels{bwii} ])
                    
                    xlabel('Time (sec)')
                    ylabel(['% correct '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    
                    
                end
            end
            %Do ranksum/t test
            fprintf(1, ['Ranksum or t-test p values for percent correct for ' handles_out.drgbchoices.bwlabels{bwii} '\n'])
            [output_data] = drgMutiRanksumorTtest(p_correct_stats);
            fprintf(1, '\n\n')
        end
        
        for groupNo=1:max(handles_out.drgbchoices.group_no)
            figure(groupNo)
            suptitle(['All trials. Group: ' handles_out.drgbchoices.group_no_names{groupNo} ' # of mice: ' num2str(no_mice_per(1)) ' ' handles_out.drgbchoices.per_lab{1} ' ' num2str(no_mice_per(2)) ' ' handles_out.drgbchoices.per_lab{2}])
        end
        
        
        %Plot PCA results
        for groupNo=1:max(handles_out.drgbchoices.group_no)
            
            
            
            
            for percent_correct_ii=1:2
                for bwii=1:4
                    
                    %Clear variables for anovan
                    figNo=figNo+1;
                    try
                        close(figNo)
                    catch
                    end
                    hFig=figure(figNo);
                    
                    
                    set(hFig, 'units','normalized','position',[.2 .2 .7 .7])
                    
                    
                    
                    
                    
                    no_mice=0;
                    all_PC1s=zeros(max(handles_out.drgbchoices.mouse_no),length(handles_out.drgbchoices.events_to_discriminate),length(t));
                    
                    for mouseNo=1:max(handles_out.drgbchoices.mouse_no)
                        
                        if isfield(handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii),'discriminant_calulated')
                            %This is here because there was a spelling mistake
                            %"calulated"
                            handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=...
                                handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calulated;
                        end
                        if handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                            per_ii=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).no_trials;
                            if per_ii>=20
                                
                                %Show PCA
                                principal_components=[];
                                principal_components=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).principal_components;
                                N=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).no_trials;
                                these_all_which_events=zeros(length(handles_out.drgbchoices.events_to_discriminate),N);
                                for ii=1:length(handles_out.drgbchoices.events_to_discriminate)
                                    kk=find(handles_out.drgbchoices.evTypeNos==handles_out.drgbchoices.events_to_discriminate(ii));
                                    these_all_which_events(ii,:)= handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(kk,:);
                                end
                                
                                %PCA before odor on
                                these_pcs=zeros(N,length(handles_out.drgbchoices.which_electrodes));
                                these_pcs(:,:)=principal_components(6,:,:);
                                subplot(2,2,3);
                                hold on
                                plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'.r')
                                plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'.b')
                                legend(handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(1)},handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(2)})
                                xlabel('PC1')
                                ylabel('PC2')
                                title('-1 sec')
                                
                                %PCA after odor on
                                subplot(2,2,4)
                                hold on
                                
                                these_pcs=zeros(N,length(handles_out.drgbchoices.which_electrodes));
                                these_pcs(:,:)=principal_components(41,:,:);
                                plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'.r')
                                plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'.b')
                                legend(handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(1)},handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(2)})
                                xlabel('PC1')
                                ylabel('PC2')
                                title('2.5 sec')
                                
                                
                                no_mice=no_mice+1;
                                
                                all_PC1s(no_mice,1,:)=mean(principal_components(:,logical(these_all_which_events(1,:)),1),2);
                                all_PC1s(no_mice,2,:)=mean(principal_components(:,logical(these_all_which_events(2,:)),1),2);
                            end
                        end
                    end
                    
                    no_mice_per(percent_correct_ii)=no_mice;
                    
                    all_PC1s=all_PC1s(1:no_mice,:,:);
                    
                    %Show the timecourse for PC1
                    subplot(2,2,[1,2])
                    hold on
                    
                    %Event 2
                    PC1ev2=zeros(no_mice,length(t));
                    PC1ev2(:,:)=all_PC1s(:,2,:);
                    mean_PC1ev2=mean(PC1ev2,1);
                    CIPC1ev2 = bootci(1000, {@mean, PC1ev2});
                    maxCIPC1ev2=max(CIPC1ev2(:));
                    minCIPC1ev2=min(CIPC1ev2(:));
                    CIPC1ev2(1,:)=mean_PC1ev2-CIPC1ev2(1,:);
                    CIPC1ev2(2,:)=CIPC1ev2(2,:)-mean_PC1ev2;
                    [hlCR, hpCR] = boundedline(t',mean_PC1ev2', CIPC1ev2', 'b');
                    
                    PC1ev1=zeros(no_mice,length(t));
                    PC1ev1(:,:)=all_PC1s(:,1,:);
                    mean_PC1ev1=mean(PC1ev1,1);
                    CIPC1ev1 = bootci(1000, {@mean, PC1ev1});
                    maxCIPC1ev1=max(CIPC1ev1(:));
                    minCIPC1ev1=min(CIPC1ev1(:));
                    CIPC1ev1(1,:)=mean_PC1ev1-CIPC1ev1(1,:);
                    CIPC1ev1(2,:)=CIPC1ev1(2,:)-mean_PC1ev1;
                    [hlCR, hpCR] = boundedline(t',mean_PC1ev1', CIPC1ev1', 'r');
                    
                    maxPC1=max([maxCIPC1ev2 maxCIPC1ev1]);
                    minPC1=min([minCIPC1ev2 minCIPC1ev1]);
                    
                    %Odor on markers
                    plot([0 0],[minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)],'-k')
                    odorhl=plot([0 2.5],[minPC1-0.1*(maxPC1-minPC1) minPC1-0.1*(maxPC1-minPC1)],'-k','LineWidth',5);
                    plot([2.5 2.5],[minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)],'-k')
                    
                    xlim([-2 5])
                    ylim([minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)])
                    text(-1,minPC1+0.9*(maxPC1-minPC1),handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(1)},'Color','r')
                    text(-1,minPC1+0.8*(maxPC1-minPC1),handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(2)},'Color','b')
                    title(['PC1 for ' handles_out.drgbchoices.bwlabels{bwii} ' mouse No ' num2str(mouseNo) ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' ' handles_out.drgbchoices.group_no_names{groupNo}])
                    xlabel('Time (sec)')
                    ylabel('PC1')
                    
                end
            end
            suptitle(['All trials. Group: ' handles_out.drgbchoices.group_no_names{groupNo} ' # of mice: ' num2str(no_mice_per(1)) ' ' handles_out.drgbchoices.per_lab{1} ' ' num2str(no_mice_per(2)) ' ' handles_out.drgbchoices.per_lab{2}])
        end
        
    case 2
        
        %Plot average percent correct for the LFP power LDA
        t=handles_out.t_pac;
        for groupNo=1:max(handles_out.drgbchoices.group_no)
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            
            set(hFig, 'units','normalized','position',[.2 .2 .7 .7])
        end
        
        ii_plot=0;
        for pac_ii=1:3
            p_correct_stats=[];
            ii_stats=0;
            for percent_correct_ii=1:2
                ii_plot=pac_ii+(percent_correct_ii-1)*4;
                for groupNo=1:max(handles_out.drgbchoices.group_no)
                    
                    figure(groupNo)
                    subplot(2,4,ii_plot);
                    hold on
                    no_mice=0;
                    all_discriminant_correct=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    all_discriminant_correct_shuffled=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    
                    for mouseNo=1:length(handles_out.discriminant_per_mousePAC)
                        if groupNo<=length(handles_out.discriminant_per_mousePAC(mouseNo).group)
                            if percent_correct_ii<=length(handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct)
                                if isfield(handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii),'discriminant_calulated')
                                    %This is here because there was a spelling mistake
                                    %"calulated"
                                    handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=...
                                        handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calulated;
                                end
                                if handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                                    per_ii=handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pac_ii).no_trials;
                                    if per_ii>=20
                                        no_mice=no_mice+1;
                                        all_discriminant_correct(no_mice,:)=handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pac_ii).discriminant_correctPAC;
                                        all_discriminant_correct_shuffled(no_mice,:)=handles_out.discriminant_per_mousePAC(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pac_ii).discriminant_correct_shuffledPAC;
                                    end
                                end
                            end
                        end
                    end
                    
                    
                    no_mice_per(percent_correct_ii)=no_mice;
                    
                    all_discriminant_correct_shuffled=all_discriminant_correct_shuffled(1:no_mice,:);
                    mean_dcsh=mean(all_discriminant_correct_shuffled,1)';
                    if size(all_discriminant_correct_shuffled,1)>2
                        CIdcsh = bootci(1000, {@mean, all_discriminant_correct_shuffled})';
                        CIdcsh(:,1)=mean_dcsh-CIdcsh(:,1);
                        CIdcsh(:,2)=CIdcsh(:,2)-mean_dcsh;
                        [hlCR, hpCR] = boundedline(t,mean_dcsh, CIdcsh, 'b');
                    else
                        plot(t,mean_dcsh,'-b')
                    end
                    
                    data=[];
                    for mouseNo=1:no_mice
                        for winNo=1:no_wins
                            data=[data mean(all_discriminant_correct_shuffled(mouseNo,(t>=window_start(winNo))&(t<=window_end(winNo))),2)];
                        end
                    end
                    
                    ii_stats=ii_stats+1;
                    p_correct_stats(ii_stats).data=data;
                    p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' shuffled'];
                    p_correct_stats(ii_stats).winNo=winNo; %Note that I use the shuffled data in both windows
                    p_correct_stats(ii_stats).pac_ii=pac_ii;
                    p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                    p_correct_stats(ii_stats).groupNo=groupNo;
                    
                    all_discriminant_correct=all_discriminant_correct(1:no_mice,:);
                    mean_dc=mean(all_discriminant_correct,1)';
                    if size(all_discriminant_correct,1)>2
                        CIdc = bootci(1000, {@mean, all_discriminant_correct})';
                        CIdc(:,1)=mean_dc-CIdc(:,1);
                        CIdc(:,2)=CIdc(:,2)-mean_dc;
                        [hlCR, hpCR] = boundedline(t,mean_dc, CIdc, 'r');
                    else
                        plot(t,mean_dc,'r')
                    end
                    
                    for winNo=1:no_wins
                        data=[];
                        for mouseNo=1:no_mice
                            data=[data mean(all_discriminant_correct(mouseNo,(t>=window_start(winNo))&(t<=window_end(winNo))),2)];
                        end
                        ii_stats=ii_stats+1;
                        p_correct_stats(ii_stats).data=data;
                        p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                            handles_out.drgbchoices.per_lab{percent_correct_ii}...
                            ' from ' num2str(window_start(winNo)) ' to ' num2str(window_end(winNo)) ' sec'];
                        p_correct_stats(ii_stats).winNo=winNo;
                        p_correct_stats(ii_stats).pac_ii=pac_ii;
                        p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                        p_correct_stats(ii_stats).groupNo=groupNo;
                    end
                    
                    %Odor on markers
                    plot([0 0],[0 100],'-k')
                    odorhl=plot([0 2.5],[10 10],'-k','LineWidth',5);
                    plot([2.5 2.5],[0 100],'-k')
                    
                    title([handles_out.drgbchoices.bwlabels{pac_ii} ])
                    
                    xlabel('Time (sec)')
                    ylabel(['% correct '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    
                    
                end
            end
            %Do ranksum/t test
            fprintf(1, ['Ranksum or t-test p values for percent correct for ' handles_out.drgbchoices.bwlabels{pac_ii} '\n'])
            try
                [output_data] = drgMutiRanksumorTtest(p_correct_stats);
                fprintf(1, '\n\n')
            catch
            end
        end
        
        for groupNo=1:max(handles_out.drgbchoices.group_no)
            figure(groupNo)
            suptitle(['All trials. Group: ' handles_out.drgbchoices.group_no_names{groupNo} ' # of mice: ' num2str(no_mice_per(1)) ' ' handles_out.drgbchoices.per_lab{1} ' ' num2str(no_mice_per(2)) ' ' handles_out.drgbchoices.per_lab{2}])
        end
        
        
        %Plot PCA results
        for groupNo=1:max(handles_out.drgbchoices.group_no)
            
            
            
            
            for percent_correct_ii=1:2
                for pca_ii=1:3
                    
                    %Clear variables for anovan
                    figNo=figNo+1;
                    try
                        close(figNo)
                    catch
                    end
                    hFig=figure(figNo);
                    
                    
                    set(hFig, 'units','normalized','position',[.2 .2 .7 .7])
                    
                    
                    
                    
                    
                    no_mice=0;
                    all_PC1s=zeros(max(handles_out.drgbchoices.mouse_no),length(handles_out.drgbchoices.events_to_discriminate),length(t));
                    
                    for mouseNo=1:length(handles_out.discriminant_per_mouse)
                        
                        
                        if handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                            per_ii=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).no_trials;
                            if per_ii>=20
                                
                                %Show PCA
                                principal_components=[];
                                principal_components=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).principal_componentsPAC;
                                N=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).no_trials;
                                these_all_which_events=zeros(length(handles_out.drgbchoices.events_to_discriminate),N);
                                for ii=1:length(handles_out.drgbchoices.events_to_discriminate)
                                    kk=find(handles_out.drgbchoices.evTypeNos==handles_out.drgbchoices.events_to_discriminate(ii));
                                    these_all_which_events(ii,:)= handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(kk,:);
                                end
                                
                                %PCA before odor on
                                szpcs=size(principal_components);
                                these_pcs=zeros(N,szpcs(3));
                                these_pcs(:,:)=principal_components(6,:,:);
                                subplot(2,2,3);
                                hold on
                                plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'.r')
                                plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'.b')
                                legend(handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(1)},handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(2)})
                                xlabel('PC1')
                                ylabel('PC2')
                                title('-1 sec')
                                
                                %PCA after odor on
                                subplot(2,2,4)
                                hold on
                                
                                szpcs=size(principal_components);
                                these_pcs=zeros(N,szpcs(3));
                                these_pcs(:,:)=principal_components(41,:,:);
                                plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'.r')
                                plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'.b')
                                legend(handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(1)},handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(2)})
                                xlabel('PC1')
                                ylabel('PC2')
                                title('2.5 sec')
                                
                                
                                no_mice=no_mice+1;
                                
                                all_PC1s(no_mice,1,:)=mean(principal_components(:,logical(these_all_which_events(1,:)),1),2);
                                all_PC1s(no_mice,2,:)=mean(principal_components(:,logical(these_all_which_events(2,:)),1),2);
                            end
                        end
                    end
                    
                    no_mice_per(percent_correct_ii)=no_mice;
                    
                    all_PC1s=all_PC1s(1:no_mice,:,:);
                    
                    %Show the timecourse for PC1
                    subplot(2,2,[1,2])
                    hold on
                    
                    
                    %Event 2
                    PC1ev2=zeros(no_mice,length(t));
                    PC1ev2(:,:)=all_PC1s(:,2,:);
                    if size(PC1ev2,1)>2
                        mean_PC1ev2=mean(PC1ev2,1);
                        CIPC1ev2 = bootci(1000, {@mean, PC1ev2});
                        maxCIPC1ev2=max(CIPC1ev2(:));
                        minCIPC1ev2=min(CIPC1ev2(:));
                        CIPC1ev2(1,:)=mean_PC1ev2-CIPC1ev2(1,:);
                        CIPC1ev2(2,:)=CIPC1ev2(2,:)-mean_PC1ev2;
                        [hlCR, hpCR] = boundedline(t',mean_PC1ev2', CIPC1ev2', 'b');
                        
                        PC1ev1=zeros(no_mice,length(t));
                        PC1ev1(:,:)=all_PC1s(:,1,:);
                        mean_PC1ev1=mean(PC1ev1,1);
                        CIPC1ev1 = bootci(1000, {@mean, PC1ev1});
                        maxCIPC1ev1=max(CIPC1ev1(:));
                        minCIPC1ev1=min(CIPC1ev1(:));
                        CIPC1ev1(1,:)=mean_PC1ev1-CIPC1ev1(1,:);
                        CIPC1ev1(2,:)=CIPC1ev1(2,:)-mean_PC1ev1;
                        [hlCR, hpCR] = boundedline(t',mean_PC1ev1', CIPC1ev1', 'r');
                        
                        maxPC1=max([maxCIPC1ev2 maxCIPC1ev1]);
                        minPC1=min([minCIPC1ev2 minCIPC1ev1]);
                        
                        %Odor on markers
                        plot([0 0],[minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)],'-k')
                        odorhl=plot([0 2.5],[minPC1-0.1*(maxPC1-minPC1) minPC1-0.1*(maxPC1-minPC1)],'-k','LineWidth',5);
                        plot([2.5 2.5],[minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)],'-k')
                        
                        xlim([-2 5])
                        ylim([minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)])
                        text(-1,minPC1+0.9*(maxPC1-minPC1),handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(1)},'Color','r')
                        text(-1,minPC1+0.8*(maxPC1-minPC1),handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(2)},'Color','b')
                        title(['PC1 for ' handles_out.drgbchoices.bwlabels{pca_ii} ' mouse No ' num2str(mouseNo) ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' ' handles_out.drgbchoices.group_no_names{groupNo}])
                        xlabel('Time (sec)')
                        ylabel('PC1')
                    end
                    
                end
            end
            suptitle(['All trials. Group: ' handles_out.drgbchoices.group_no_names{groupNo} ' # of mice: ' num2str(no_mice_per(1)) ' ' handles_out.drgbchoices.per_lab{1} ' ' num2str(no_mice_per(2)) ' ' handles_out.drgbchoices.per_lab{2}])
        end
        
    case 3
        
        %Plot average percent correct for the LDA for peak and trough for
        %wavelet power referenced to PAC phase
        t=handles_out.t_power;
        handles_out2=[];
        
        all_concs=[0.032 0.1 0.32 1 3.2 10];
        
        pcorrs_beh=zeros(2,6,20);
        pcorrs_beh_no=zeros(2,6);
        
        for PACii=1:length(handles_out.drgbchoices.PACburstLowF)
            
            glm_auc=[];
            glma_ii=0;
            glm_aucp=[];
            glmp_ii=0;
            glm_auct=[];
            glmt_ii=0;
            
            iip_stats=0;
            aucp_stats=[];
            iit_stats=0;
            auct_stats=[];
            
            p_correct_stats=[];
            ii_stats=0;
            
            glm_pc=[];
            glmpc_ii=0;
            iipc_stats=0;
            pc_stats=[];
            
            glm_pk=[];
            glmpk_ii=0;
            iipk_stats=0;
            pk_stats=[];
            
            glm_th=[];
            glmth_ii=0;
            iith_stats=0;
            th_stats=[];
            
            pcorrs_pk=zeros(2,6,20);
            pcorrs_pk_no=zeros(2,6);
            pcorrs_th=zeros(2,6,20);
            pcorrs_th_no=zeros(2,6);

            
            mice_included=zeros(2,length(handles_out.discriminant_PACwavepower));
            
            
            figNo_peak=figNo+5;
            
            try
                close(figNo_peak)
            catch
            end
            
            try
                close(figNo_peak+1)
            catch
            end
            
            %Percent correct per conc
            figNo_pconc=figNo+7;
            try
                close(figNo_pconc)
            catch
            end
            hFig=figure(figNo_pconc);
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            
            try
                close(figNo_pconc+1)
            catch
            end
            hFig=figure(figNo_pconc+1);
            set(hFig, 'units','normalized','position',[.1 .5 .7 .4])
            
            %Do behavioral percent correct only for PACii==3
            if PACii==1
                figNo_pc=figNo+9;
                try
                    close(figNo_pc)
                catch
                end
                hFig=figure(figNo_pc);
                set(hFig, 'units','normalized','position',[.1 .5 .95 .4])
            end
            
            peak_bar_offfset=0;
            trough_bar_offfset=0;
            
            for groupNo=1:max(handles_out.drgbchoices.group_no)
                
                for percent_correct_ii=2:-1:1
                    
                    %Gather all the data
                    no_mice=0;
                    no_mice_included=0;
                    
                    all_discriminant_correct_peak=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    all_discriminant_correct_trough=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    all_discriminant_correct_shuffled_peak=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    all_discriminant_correct_shuffled_trough=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    all_d_per_trial_correct_peak=zeros(max(handles_out.drgbchoices.mouse_no),1000,length(t));
                    all_d_per_trial_correct_trough=zeros(max(handles_out.drgbchoices.mouse_no),1000,length(t));
                    all_d_per_trial_no_trials=zeros(1,max(handles_out.drgbchoices.mouse_no));
                    all_d_per_trial_conc=zeros(max(handles_out.drgbchoices.mouse_no),1000);
                    all_d_per_trial_correct_behavior=zeros(max(handles_out.drgbchoices.mouse_no),1000);
                    
                    for mouseNo=1:length(handles_out.discriminant_PACwavepower)
                        try
                            if handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                                per_ii=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials;
                                no_mice=no_mice+1; 
                                if (per_ii>=20)&(sum(no_mice==mice_excluded)==0)
                                    no_mice_included=no_mice_included+1;
                                    if percent_correct_ii==1
                                        mice_included(groupNo,no_mice_included)=mouseNo;
                                    end
                                    all_discriminant_correct_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_peak;
                                    all_discriminant_correct_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_trough;
                                    all_discriminant_correct_shuffled_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_peak;
                                    all_discriminant_correct_shuffled_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_trough;
                                    all_dimensionality_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak;
                                    all_dimensionality_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough;
                                    all_d_per_trial_no_trials(1,no_mice_included)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials;
                                    for trNo=1:all_d_per_trial_no_trials(1,no_mice_included)
                                        if (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(1,trNo)==1) ...
                                                || (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(4,trNo)==1)
                                            all_d_per_trial_correct_behavior(no_mice_included,trNo)=1;
                                        else
                                            all_d_per_trial_correct_behavior(no_mice_included,trNo)=0;
                                        end
                                        for ii_t=1:length(t)
                                            if (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(2,trNo)==1) ...
                                                    & (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).test_out_per_timepoint_peak(1,trNo,ii_t)==1)
                                                all_d_per_trial_correct_peak(no_mice_included,trNo,ii_t)=1;
                                            else
                                                if (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(5,trNo)==1) ...
                                                        & (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).test_out_per_timepoint_peak(2,trNo,ii_t)==1)
                                                    all_d_per_trial_correct_peak(no_mice_included,trNo,ii_t)=1;
                                                else
                                                    all_d_per_trial_correct_peak(no_mice_included,trNo,ii_t)=0;
                                                end
                                            end
                                            
                                            if (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(2,trNo)==1) ...
                                                    & (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).test_out_per_timepoint_trough(1,trNo,ii_t)==1)
                                                all_d_per_trial_correct_trough(no_mice_included,trNo,ii_t)=1;
                                            else
                                                if (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(5,trNo)==1) ...
                                                        & (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).test_out_per_timepoint_trough(2,trNo,ii_t)==1)
                                                    all_d_per_trial_correct_trough(no_mice_included,trNo,ii_t)=1;
                                                else
                                                    all_d_per_trial_correct_trough(no_mice_included,trNo,ii_t)=0;
                                                end
                                            end
                                        end
                                        
                                        if (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(7,trNo)==1)
                                            if groupNo==1
                                                all_d_per_trial_conc(no_mice_included,trNo)=10;
                                            else
                                                all_d_per_trial_conc(no_mice_included,trNo)=0.032;
                                            end
                                        else
                                            if (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(8,trNo)==1)
                                                if groupNo==1
                                                    all_d_per_trial_conc(no_mice_included,trNo)=3.2;
                                                else
                                                    all_d_per_trial_conc(no_mice_included,trNo)=0.1;
                                                end
                                            else
                                                if (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(9,trNo)==1)
                                                    if groupNo==1
                                                        all_d_per_trial_conc(no_mice_included,trNo)=1;
                                                    else
                                                        all_d_per_trial_conc(no_mice_included,trNo)=0.32;
                                                    end
                                                else
                                                    if (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(10,trNo)==1)
                                                        if groupNo==1
                                                            all_d_per_trial_conc(no_mice_included,trNo)=0.32;
                                                        else
                                                            all_d_per_trial_conc(no_mice_included,trNo)=1;
                                                        end
                                                    else
                                                        if (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(11,trNo)==1)
                                                            if groupNo==1
                                                                all_d_per_trial_conc(no_mice_included,trNo)=0.1;
                                                            else
                                                                all_d_per_trial_conc(no_mice_included,trNo)=3.2;
                                                            end
                                                        else
                                                            if (handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(12,trNo)==1)
                                                                if groupNo==1
                                                                    all_d_per_trial_conc(no_mice_included,trNo)=0.032;
                                                                else
                                                                    all_d_per_trial_conc(no_mice_included,trNo)=10;
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
                                pffft=1;
                            end
                        catch
                            pffft=1; 
                        end
                    end
                    
                    fprintf(1, ['The number of mice included in the LDA analysis for this odor pair is %d\n\n\n'], no_mice_included)
                    
                    no_mice_per(percent_correct_ii)=no_mice_included;
                    
                    %Plot percent correct for the LDA and save the data for
                    %the ranksum
                    figNo=figNo+1;
                    try
                        close(figNo)
                    catch
                    end
                    hFig=figure(figNo);
                    
                    hold on
                    
                    %Note that I merge the shuffled for both peak and
                    %trough for the average plot
                    all_discriminant_correct_shuffled=zeros(2*no_mice_included,length(t));
                    all_discriminant_correct_shuffled(1:no_mice_included,:)=all_discriminant_correct_shuffled_peak(1:no_mice_included,:);
                    all_discriminant_correct_shuffled(no_mice_included+1:end,:)=all_discriminant_correct_shuffled_trough(1:no_mice_included,:);
                    mean_dcsh=mean(all_discriminant_correct_shuffled,1)';
                    if size(all_discriminant_correct_shuffled,1)>2
                        CIdcsh = bootci(1000, {@mean, all_discriminant_correct_shuffled})';
                        CIdcsh(:,1)=mean_dcsh-CIdcsh(:,1);
                        CIdcsh(:,2)=CIdcsh(:,2)-mean_dcsh;
                        [hlCR, hpCR] = boundedline(t,mean_dcsh, CIdcsh, 'k');
                    else
                        plot(t,mean_dcsh,'-k')
                    end
                    
                    data=[];
                    for mouseNo=1:no_mice_included
                        data=[data (mean(all_discriminant_correct_shuffled_peak(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
                    end
                    
                    ii_stats=ii_stats+1;
                    p_correct_stats(ii_stats).data=data;
                    %                     p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' shuffled peak'];
                    %                     p_correct_stats(ii_stats).data_ii=1;
                    %                     p_correct_stats(ii_stats).PACii=PACii;
                    %                     p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                    %                     p_correct_stats(ii_stats).groupNo=groupNo;
                    %
                    %
                    %                     glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
                    %                     glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
                    %                     glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
                    %                     glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
                    %                     glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=1;
                    %                     glm_correct.peak(glm_ii+1:glm_ii+length(data))=1;
                    %                     glm_ii=glm_ii+length(data);
                    %
                    %                     data=[];
                    %                     for mouseNo=1:no_mice_included
                    %                         data=[data (mean(all_discriminant_correct_shuffled_trough(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
                    %                     end
                    %
                    ii_stats=ii_stats+1;
                    p_correct_stats(ii_stats).data=data;
                    %                     p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' shuffled trough'];
                    %                     p_correct_stats(ii_stats).data_ii=2;
                    %                     p_correct_stats(ii_stats).PACii=PACii;
                    %                     p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                    %                     p_correct_stats(ii_stats).groupNo=groupNo;
                    %
                    %                     glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
                    %                     glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
                    %                     glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
                    %                     glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
                    %                     glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=1;
                    %                     glm_correct.peak(glm_ii+1:glm_ii+length(data))=0;
                    %                     glm_ii=glm_ii+length(data);
                    %
                    
                    
                    %Now plot the percent correct for the trough
                    all_discriminant_correct_trough=all_discriminant_correct_trough(1:no_mice_included,:);
                    mean_dc_trough=mean(all_discriminant_correct_trough,1)';
                    if size(all_discriminant_correct_trough,1)>2
                        CIdc_trough = bootci(1000, {@mean, all_discriminant_correct_trough})';
                        CIdc_trough(:,1)=mean_dc_trough-CIdc_trough(:,1);
                        CIdc_trough(:,2)=CIdc_trough(:,2)-mean_dc_trough;
                        [hlCR, hpCR] = boundedline(t,mean_dc_trough, CIdc_trough, 'b');
                    else
                        plot(t,mean_dc_trough,'r')
                    end
                    
                    
                    data=[];
                    for mouseNo=1:no_mice_included
                        data=[data (mean(all_discriminant_correct_trough(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
                    end
                    ii_stats=ii_stats+1;
                    p_correct_stats(ii_stats).data=data;
                    %                     p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' trough'];
                    %                     p_correct_stats(ii_stats).data_ii=3;
                    %                     p_correct_stats(ii_stats).PACii=PACii;
                    %                     p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                    %                     p_correct_stats(ii_stats).groupNo=groupNo;
                    %
                    %                     glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
                    %                     glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
                    %                     glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
                    %                     glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
                    %                     glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=0;
                    %                     glm_correct.peak(glm_ii+1:glm_ii+length(data))=0;
                    %                     glm_ii=glm_ii+length(data);
                    
                    %Now plot the percent correct for the peak
                    all_discriminant_correct_peak=all_discriminant_correct_peak(1:no_mice_included,:);
                    mean_dc_peak=mean(all_discriminant_correct_peak,1)';
                    if size(all_discriminant_correct_peak,1)>2
                        CIdc_peak = bootci(1000, {@mean, all_discriminant_correct_peak})';
                        CIdc_peak(:,1)=mean_dc_peak-CIdc_peak(:,1);
                        CIdc_peak(:,2)=CIdc_peak(:,2)-mean_dc_peak;
                        [hlCR, hpCR] = boundedline(t,mean_dc_peak, CIdc_peak, 'r');
                    else
                        plot(t,mean_dc_peak,'r')
                    end
                    
                    
                    data=[];
                    for mouseNo=1:no_mice_included
                        data=[data (mean(all_discriminant_correct_peak(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
                    end
                    ii_stats=ii_stats+1;
                    p_correct_stats(ii_stats).data=data;
                    %                     p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                    %                         handles_out.drgbchoices.per_lab{percent_correct_ii}...
                    %                         ' peak'];
                    %                     p_correct_stats(ii_stats).data_ii=4;
                    %                     p_correct_stats(ii_stats).PACii=PACii;
                    %                     p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                    %                     p_correct_stats(ii_stats).groupNo=groupNo;
                    %
                    %                     glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
                    %                     glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
                    %                     glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
                    %                     glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
                    %                     glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=0;
                    %                     glm_correct.peak(glm_ii+1:glm_ii+length(data))=1;
                    %                     glm_ii=glm_ii+length(data);
                    
                    %Odor on markers
                    plot([0 0],[0 100],'-k')
                    odorhl=plot([0 2.5],[10 10],'-k','LineWidth',5);
                    plot([2.5 2.5],[0 100],'-k')
                    
                    title(['LDA for theta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    
                    xlabel('Time (sec)')
                    ylabel(['% correct '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    legend('Shuffled','Trough','Peak')
                    
                    %                     %Plot the peak vs trough graph
                    %                     figNo=figNo+1;
                    %                     try
                    %                         close(figNo)
                    %                     catch
                    %                     end
                    %                     hFig=figure(figNo);
                    %
                    %                     hold on
                    %
                    %                     plot(p_correct_stats(ii_stats-1).data,p_correct_stats(ii_stats).data,'or')
                    %                     plot(p_correct_stats(ii_stats-2).data,p_correct_stats(ii_stats-3).data,'ok')
                    %
                    %
                    %                     plot([-0.1 1],[-0.1 1],'-k')
                    %
                    %                     xlim([-0.1 1])
                    %                     ylim([-0.1 1])
                    %                     ylabel('AUC peak')
                    %                     xlabel('AUC trough')
                    %                     legend('Original','Shuffled')
                    %                     title(['Area under the curve for theta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    %
                    %Plot the bar graphs
                    
                    %Peak bar
                    hFig=figure(figNo_peak);
                    hold on
                    
                    peak_bar_offfset = peak_bar_offfset+1;
                    
                    %Naive or proficient
                    if percent_correct_ii==1
                        bar(peak_bar_offfset,mean(p_correct_stats(ii_stats).data),'LineWidth',3,'FaceColor',[204/255,121/255,167/255],'EdgeColor','none')
                    else
                        bar(peak_bar_offfset,mean(p_correct_stats(ii_stats).data),'LineWidth',3,'FaceColor',[0/255,158/255,115/255],'EdgeColor','none')
                    end
                    
                    CI = bootci(1000, {@mean, p_correct_stats(ii_stats).data},'type','cper');
                    plot([peak_bar_offfset peak_bar_offfset],CI,'-k','LineWidth',3)
                    plot(peak_bar_offfset*ones(1,length(p_correct_stats(ii_stats).data)),p_correct_stats(ii_stats).data,'o','MarkerSize',5,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                    
                    %Save glm
                    these_data=p_correct_stats(ii_stats).data;
                    glm_auc.data(glma_ii+1:glma_ii+length(these_data))=these_data;
                    glm_auc.pcorr(glma_ii+1:glma_ii+length(these_data))=percent_correct_ii*ones(1,length(these_data));
                    glm_auc.rewarded(glma_ii+1:glma_ii+length(these_data))=groupNo*ones(1,length(these_data));
                    glm_auc.peak(glma_ii+1:glma_ii+length(these_data))=ones(1,length(these_data));
                    glma_ii=glma_ii+length(these_data);
                    
                    glm_aucp.data(glmp_ii+1:glmp_ii+length(these_data))=these_data;
                    glm_aucp.pcorr(glmp_ii+1:glmp_ii+length(these_data))=percent_correct_ii*ones(1,length(these_data));
                    glm_aucp.rewarded(glmp_ii+1:glmp_ii+length(these_data))=groupNo*ones(1,length(these_data));
                    glmp_ii=glmp_ii+length(these_data);
                    
                    iip_stats=iip_stats+1;
                    aucp_stats(iip_stats).data=these_data;
                    aucp_stats(iip_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                        handles_out.drgbchoices.per_lab{percent_correct_ii}];
                    
                    %Shuffled
                    if percent_correct_ii==1
                        these_peak_shuffled=[these_peak_shuffled p_correct_stats(ii_stats-3).data];
                        bar(peak_bar_offfset-2,mean(these_peak_shuffled),'LineWidth',3,'FaceColor',[86/255,180/255,233/255],'EdgeColor','none')
                        CI = bootci(1000, {@mean, these_peak_shuffled},'type','cper');
                        plot([peak_bar_offfset-2 peak_bar_offfset-2],CI,'-k','LineWidth',3)
                        plot((peak_bar_offfset-2)*ones(1,length(these_peak_shuffled)),these_peak_shuffled,'o','MarkerSize',5,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                        
                        %Save glm
                        these_data=p_correct_stats(ii_stats-3).data;
                        glm_auc.data(glma_ii+1:glma_ii+length(these_data))=these_data;
                        glm_auc.pcorr(glma_ii+1:glma_ii+length(these_data))=3*ones(1,length(these_data));
                        glm_auc.rewarded(glma_ii+1:glma_ii+length(these_data))=groupNo*ones(1,length(these_data));
                        glm_auc.peak(glma_ii+1:glma_ii+length(these_data))=ones(1,length(these_data));
                        glma_ii=glma_ii+length(these_data);
                        
                        glm_aucp.data(glmp_ii+1:glmp_ii+length(these_data))=these_data;
                        glm_aucp.pcorr(glmp_ii+1:glmp_ii+length(these_data))=3*ones(1,length(these_data));
                        glm_aucp.rewarded(glmp_ii+1:glmp_ii+length(these_data))=groupNo*ones(1,length(these_data));
                        glmp_ii=glmp_ii+length(these_data);
                        
                        iip_stats=iip_stats+1;
                        aucp_stats(iip_stats).data=these_data;
                        aucp_stats(iip_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                            'Shuffled'];
                    else
                        these_peak_shuffled=[];
                        these_peak_shuffled=p_correct_stats(ii_stats-3).data;
                    end
                    
                    
                    
                    %Trough bar
                    hFig=figure(figNo_peak+1);
                    hold on
                    
                    trough_bar_offfset=trough_bar_offfset+1;
                    
                    if percent_correct_ii==1
                        bar(trough_bar_offfset,mean(p_correct_stats(ii_stats-1).data),'LineWidth',3,'FaceColor',[204/255,121/255,167/255],'EdgeColor','none')
                    else
                        bar(trough_bar_offfset,mean(p_correct_stats(ii_stats-1).data),'LineWidth',3,'FaceColor',[0/255,158/255,115/255],'EdgeColor','none')
                    end
                    
                    CI = bootci(1000, {@mean, p_correct_stats(ii_stats-1).data},'type','cper');
                    plot([trough_bar_offfset trough_bar_offfset],CI,'-k','LineWidth',3)
                    plot(trough_bar_offfset*ones(1,length(p_correct_stats(ii_stats-1).data)),p_correct_stats(ii_stats-1).data,'o','MarkerSize',5,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                    
                    %Save glm
                    these_data=p_correct_stats(ii_stats-1).data;
                    glm_auc.data(glma_ii+1:glma_ii+length(these_data))=these_data;
                    glm_auc.pcorr(glma_ii+1:glma_ii+length(these_data))=percent_correct_ii*ones(1,length(these_data));
                    glm_auc.rewarded(glma_ii+1:glma_ii+length(these_data))=groupNo*ones(1,length(these_data));
                    glm_auc.peak(glma_ii+1:glma_ii+length(these_data))=zeros(1,length(these_data));
                    glma_ii=glma_ii+length(these_data);
                    
                    glm_auct.data(glmt_ii+1:glmt_ii+length(these_data))=these_data;
                    glm_auct.pcorr(glmt_ii+1:glmt_ii+length(these_data))=percent_correct_ii*ones(1,length(these_data));
                    glm_auct.rewarded(glmt_ii+1:glmt_ii+length(these_data))=groupNo*ones(1,length(these_data));
                    glmt_ii=glmt_ii+length(these_data);
                    
                    iit_stats=iit_stats+1;
                    auct_stats(iit_stats).data=these_data;
                    auct_stats(iit_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                        handles_out.drgbchoices.per_lab{percent_correct_ii}];
                    
                    %Shuffled
                    if percent_correct_ii==1
                        these_trough_shuffled=[these_trough_shuffled p_correct_stats(ii_stats-2).data];
                        bar(trough_bar_offfset-2,mean(these_trough_shuffled),'LineWidth',3,'FaceColor',[86/255,180/255,233/255],'EdgeColor','none')
                        CI = bootci(1000, {@mean, these_trough_shuffled},'type','cper');
                        plot([trough_bar_offfset-2 trough_bar_offfset-2],CI,'-k','LineWidth',3)
                        plot((trough_bar_offfset-2)*ones(1,length(these_trough_shuffled)),these_trough_shuffled,'o','MarkerSize',5,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                        
                        %Save glm
                        these_data=p_correct_stats(ii_stats-2).data;
                        glm_auc.data(glma_ii+1:glma_ii+length(these_data))=these_data;
                        glm_auc.pcorr(glma_ii+1:glma_ii+length(these_data))=3*ones(1,length(these_data));
                        glm_auc.rewarded(glma_ii+1:glma_ii+length(these_data))=groupNo*ones(1,length(these_data));
                        glm_auc.peak(glma_ii+1:glma_ii+length(these_data))=zeros(1,length(these_data));
                        glma_ii=glma_ii+length(these_data);
                        
                        glm_auct.data(glmt_ii+1:glmt_ii+length(these_data))=these_data;
                        glm_auct.pcorr(glmt_ii+1:glmt_ii+length(these_data))=3*ones(1,length(these_data));
                        glm_auct.rewarded(glmt_ii+1:glmt_ii+length(these_data))=groupNo*ones(1,length(these_data));
                        glmt_ii=glmt_ii+length(these_data);
                        
                        iit_stats=iit_stats+1;
                        auct_stats(iit_stats).data=these_data;
                        auct_stats(iit_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                            'Shuffled'];
                    else
                        these_trough_shuffled=[];
                        these_trough_shuffled=p_correct_stats(ii_stats-2).data;
                    end
                    
                    
                    %                     %Plot dimensionality and save the data for
                    %                     %the ranksum
                    %
                    %                     winNo=1;
                    %                     figNo=figNo+1;
                    %                     try
                    %                         close(figNo)
                    %                     catch
                    %                     end
                    %                     hFig=figure(figNo);
                    %
                    %                     hold on
                    %
                    %
                    %                     %Now plot the dimensionality for the trough
                    %                     all_dimensionality_trough=all_dimensionality_trough(1:no_mice_included,:);
                    %                     mean_dc_trough=mean(all_dimensionality_trough,1)';
                    %                     if size(all_dimensionality_trough,1)>2
                    %                         CIdc_trough = bootci(1000, {@mean, all_dimensionality_trough})';
                    %                         CIdc_trough(:,1)=mean_dc_trough-CIdc_trough(:,1);
                    %                         CIdc_trough(:,2)=CIdc_trough(:,2)-mean_dc_trough;
                    %                         [hltrough, hptrough] = boundedline(t,mean_dc_trough, CIdc_trough, 'b');
                    %                     else
                    %                         plot(t,mean_dc_trough,'b')
                    %                     end
                    %
                    %                     for winNo=1:no_wins
                    %                         data=[];
                    %                         for mouseNo=1:no_mice_included
                    %                             data=[data mean(all_dimensionality_trough(mouseNo,(t>=window_start(winNo))&(t<=window_end(winNo))),2)];
                    %                         end
                    %                         ii_dim_stats=ii_dim_stats+1;
                    %                         dim_stats(ii_dim_stats).data=data;
                    %                         if winNo==1
                    %                             dim_stats(ii_dim_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                    %                                 handles_out.drgbchoices.per_lab{percent_correct_ii} ' pre-odor, trough'];
                    %                         else
                    %                             dim_stats(ii_dim_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                    %                                 handles_out.drgbchoices.per_lab{percent_correct_ii} ' odor, trough'];
                    %                         end
                    %                         dim_stats(ii_dim_stats).winNo=winNo;
                    %                         dim_stats(ii_dim_stats).PACii=PACii;
                    %                         dim_stats(ii_dim_stats).per_corr_ii=percent_correct_ii;
                    %                         dim_stats(ii_dim_stats).groupNo=groupNo;
                    %
                    %                         glm_dim.data(glm_dim_ii+1:glm_dim_ii+length(data))=data;
                    %                         glm_dim.PACii(glm_dim_ii+1:glm_dim_ii+length(data))=PACii;
                    %                         glm_dim.group(glm_dim_ii+1:glm_dim_ii+length(data))=groupNo;
                    %                         glm_dim.odor_window(glm_dim_ii+1:glm_dim_ii+length(data))=winNo;
                    %                         glm_dim.perCorr(glm_dim_ii+1:glm_dim_ii+length(data))=percent_correct_ii;
                    %                         glm_dim.peak(glm_dim_ii+1:glm_dim_ii+length(data))=0;
                    %                         glm_dim_ii=glm_dim_ii+length(data);
                    %                     end
                    %
                    %                     %Now plot the percent correct for the peak
                    %                     all_dimensionality_peak=all_dimensionality_peak(1:no_mice_included,:);
                    %                     mean_dc_peak=mean(all_dimensionality_peak,1)';
                    %                     if size(all_dimensionality_peak,1)>2
                    %                         CIdc_peak = bootci(1000, {@mean, all_dimensionality_peak})';
                    %                         CIdc_peak(:,1)=mean_dc_peak-CIdc_peak(:,1);
                    %                         CIdc_peak(:,2)=CIdc_peak(:,2)-mean_dc_peak;
                    %                         [hlpeak, hppeak] = boundedline(t,mean_dc_peak, CIdc_peak, 'r');
                    %                     else
                    %                         plot(t,mean_dc_peak,'r')
                    %                     end
                    %
                    %                     for winNo=1:no_wins
                    %                         data=[];
                    %                         for mouseNo=1:no_mice_included
                    %                             data=[data mean(all_dimensionality_peak(mouseNo,(t>=window_start(winNo))&(t<=window_end(winNo))),2)];
                    %                         end
                    %                         ii_dim_stats=ii_dim_stats+1;
                    %                         dim_stats(ii_dim_stats).data=data;
                    %                         if winNo==1
                    %                             dim_stats(ii_dim_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                    %                                 handles_out.drgbchoices.per_lab{percent_correct_ii} ' pre-odor, peak'];
                    %                         else
                    %                             dim_stats(ii_dim_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                    %                                 handles_out.drgbchoices.per_lab{percent_correct_ii} ' odor, peak'];
                    %                         end
                    %                         dim_stats(ii_dim_stats).winNo=winNo;
                    %                         dim_stats(ii_dim_stats).PACii=PACii;
                    %                         dim_stats(ii_dim_stats).per_corr_ii=percent_correct_ii;
                    %                         dim_stats(ii_dim_stats).groupNo=groupNo;
                    %
                    %                         glm_dim.data(glm_dim_ii+1:glm_dim_ii+length(data))=data;
                    %                         glm_dim.PACii(glm_dim_ii+1:glm_dim_ii+length(data))=PACii;
                    %                         glm_dim.group(glm_dim_ii+1:glm_dim_ii+length(data))=groupNo;
                    %                         glm_dim.odor_window(glm_dim_ii+1:glm_dim_ii+length(data))=winNo;
                    %                         glm_dim.perCorr(glm_dim_ii+1:glm_dim_ii+length(data))=percent_correct_ii;
                    %                         glm_dim.peak(glm_dim_ii+1:glm_dim_ii+length(data))=1;
                    %                         glm_dim_ii=glm_dim_ii+length(data);
                    %                     end
                    %
                    %                     %Odor on markers
                    %                     this_yl=ylim;
                    %                     plot([0 0],this_yl,'-k')
                    %                     odorhl=plot([0 2.5],[this_yl(1)+0.1*(this_yl(2)-this_yl(1)) this_yl(1)+0.1*(this_yl(2)-this_yl(1))],'-k','LineWidth',5);
                    %                     plot([2.5 2.5],this_yl,'-k')
                    %
                    %                     title(['Dimensionality for theta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    %                     legend([hltrough hlpeak],{'Trough','Peak'})
                    %                     xlabel('Time (sec)')
                    %                     ylabel(['Dimensionality '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    %
                    
                    pffft=1;
                    
                    %Now plot the per concentration percent correct
                    if percent_correct_ii==1
                        min_trials=30;
                        
                        %Peak percent correct
                        edges=[0:5:100];
                        rand_offset=0.8;
                        
                        hFig=figure(figNo_pconc);
                        hold on
                        
                        pcorrs_per_mouse=[];
                        m_num=0;
                        for mouseNo=1:no_mice_included
                            %Include only if there are more than min_trials 
                            if all_d_per_trial_no_trials(mouseNo)>min_trials
                                m_num=m_num+1;
                                pk_mice(groupNo,m_num)=mice_included(groupNo,mouseNo);
                                for concNo=1:6
                                    these_corr_b=ones(1,all_d_per_trial_no_trials(mouseNo));
                                    these_corr_b(1,:)=mean(all_d_per_trial_correct_peak(mouseNo,1:all_d_per_trial_no_trials(mouseNo),(t>=2)&(t<=2.5)),3);
                                    these_concs=[];
                                    these_concs=all_d_per_trial_conc(mouseNo,1:all_d_per_trial_no_trials(mouseNo));
                                    pcorrs_per_mouse(concNo,m_num)=100*sum(these_corr_b(these_concs==all_concs(concNo)))/sum(these_concs==all_concs(concNo));
                                end
                            end
                        end
                        
                         fprintf(1, ['No of mice included for the behavioral percent correct %d\n\n'], m_num)
                   
                        %Add bars
                        for concNo=1:6
                            these_pcorrs=zeros(1,m_num);
                            these_pcorrs(1,:)=pcorrs_per_mouse(concNo,:);
                            pcorrs_pk(groupNo,concNo,1:m_num)=these_pcorrs;
                            pcorrs_pk_no(groupNo,concNo)=m_num;
                            if groupNo==1
                                pc_bar_offset=3*(concNo-1)+1;
                                bar(pc_bar_offset,mean(these_pcorrs),'LineWidth',3,'FaceColor',[213/255,94/255,0/255],'EdgeColor','none')
                            else
                                pc_bar_offset=3*(concNo-1)+2;
                                bar(pc_bar_offset,mean(these_pcorrs),'LineWidth',3,'FaceColor',[86/255,180/255,233/255],'EdgeColor','none')
                            end
                            
                            %                             CI = bootci(1000, {@mean, these_pcorrs},'type','cper');
                            %                             plot([pc_bar_offset pc_bar_offset],CI,'-k','LineWidth',3)
                            %                             plot(pc_bar_offset*ones(1,length(these_pcorrs)),these_pcorrs,'o','MarkerSize',5,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            %
                            
                            [mean_out, CIout]=drgViolinPoint(these_pcorrs...
                                ,edges,pc_bar_offset,rand_offset,'k','k',6);
                            
                            %Save glm
                            glm_pk.data(glmpk_ii+1:glmpk_ii+length(these_pcorrs))=these_pcorrs;
                            glm_pk.group(glmpk_ii+1:glmpk_ii+length(these_pcorrs))=groupNo*ones(1,length(these_pcorrs));
                            glm_pk.rewarded(glmpk_ii+1:glmpk_ii+length(these_pcorrs))=groupNo*ones(1,length(these_pcorrs));
                            glm_pk.concs(glmpk_ii+1:glmpk_ii+length(these_pcorrs))=log10(all_concs(concNo))*ones(1,length(these_pcorrs));
                            glmpk_ii=glmpk_ii+length(these_pcorrs);
                            
                            iipk_stats=iipk_stats+1;
                            pk_stats(iipk_stats).data=these_pcorrs;
                            pk_stats(iipk_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                                num2str(all_concs(concNo))];
                        end
                        
                        
                         min_trials=30;
                        
                        %Trough percent correct
                        edges=[0:5:100];
                        rand_offset=0.8;
                        
                        hFig=figure(figNo_pconc+1);
                        hold on
                        
                        pcorrs_per_mouse=[];
                        m_num=0;
                        for mouseNo=1:no_mice_included
                            %Include only if there are more than min_trials 
                            if all_d_per_trial_no_trials(mouseNo)>min_trials
                                m_num=m_num+1;
                                th_mice(groupNo,m_num)=mice_included(groupNo,mouseNo);
                                for concNo=1:6
                                    these_corr_b=ones(1,all_d_per_trial_no_trials(mouseNo));
                                    these_corr_b(1,:)=mean(all_d_per_trial_correct_trough(mouseNo,1:all_d_per_trial_no_trials(mouseNo),(t>=2)&(t<=2.5)),3);
                                    these_concs=[];
                                    these_concs=all_d_per_trial_conc(mouseNo,1:all_d_per_trial_no_trials(mouseNo));
                                    pcorrs_per_mouse(concNo,m_num)=100*sum(these_corr_b(these_concs==all_concs(concNo)))/sum(these_concs==all_concs(concNo));
                                end
                            end
                        end
                        
                        %Add bars
                        for concNo=1:6
                            these_pcorrs=zeros(1,m_num);
                            these_pcorrs(1,:)=pcorrs_per_mouse(concNo,:);
                            pcorrs_th(groupNo,concNo,1:m_num)=these_pcorrs;
                            pcorrs_th_no(groupNo,concNo)=m_num;
                            if groupNo==1
                                pc_bar_offset=3*(concNo-1)+1;
                                bar(pc_bar_offset,mean(these_pcorrs),'LineWidth',3,'FaceColor',[213/255,94/255,0/255],'EdgeColor','none')
                            else
                                pc_bar_offset=3*(concNo-1)+2;
                                bar(pc_bar_offset,mean(these_pcorrs),'LineWidth',3,'FaceColor',[86/255,180/255,233/255],'EdgeColor','none')
                            end
                            
                            %                             CI = bootci(1000, {@mean, these_pcorrs},'type','cper');
                            %                             plot([pc_bar_offset pc_bar_offset],CI,'-k','LineWidth',3)
                            %                             plot(pc_bar_offset*ones(1,length(these_pcorrs)),these_pcorrs,'o','MarkerSize',5,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                            %
                            
                            [mean_out, CIout]=drgViolinPoint(these_pcorrs...
                                ,edges,pc_bar_offset,rand_offset,'k','k',6);
                            
                            %Save glm
                            glm_th.data(glmth_ii+1:glmth_ii+length(these_pcorrs))=these_pcorrs;
                            glm_th.group(glmth_ii+1:glmth_ii+length(these_pcorrs))=groupNo*ones(1,length(these_pcorrs));
                            glm_th.rewarded(glmth_ii+1:glmth_ii+length(these_pcorrs))=groupNo*ones(1,length(these_pcorrs));
                            glm_th.concs(glmth_ii+1:glmth_ii+length(these_pcorrs))=log10(all_concs(concNo))*ones(1,length(these_pcorrs));
                            glmth_ii=glmth_ii+length(these_pcorrs);
                            
                            iith_stats=iith_stats+1;
                            th_stats(iith_stats).data=these_pcorrs;
                            th_stats(iith_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                                num2str(all_concs(concNo))];
                        end
                        
                        
                        %Behavioral percent correct
                        if PACii==1
                            edges=[0:5:100];
                            rand_offset=0.8;
                            
                            hFig=figure(figNo_pc);
                            hold on
                            
                            pcorrs_per_mouse=[];
                            m_num=0;
                            for mouseNo=1:no_mice_included
                                %Include only if there are mre than min_trials
                                if all_d_per_trial_no_trials(mouseNo)>min_trials
                                    m_num=m_num+1;
                                    pcorrs_mice(groupNo,m_num)=mice_included(groupNo,mouseNo);
                                    for concNo=1:6
                                        these_corr_b=[];
                                        these_corr_b=all_d_per_trial_correct_behavior(mouseNo,1:all_d_per_trial_no_trials(mouseNo));
                                        these_concs=[];
                                        these_concs=all_d_per_trial_conc(mouseNo,1:all_d_per_trial_no_trials(mouseNo));
                                        pcorrs_per_mouse(concNo,m_num)=100*sum(these_corr_b(these_concs==all_concs(concNo)))/sum(these_concs==all_concs(concNo));
                                    end
                                end
                            end
                            
                            %Add bars
                            for concNo=1:6
                                these_pcorrs=zeros(1,m_num);
                                these_pcorrs(1,:)=pcorrs_per_mouse(concNo,:);
                                pcorrs_beh(groupNo,concNo,1:m_num)=these_pcorrs;
                                pcorrs_beh_no(groupNo,concNo)=m_num;
                                if groupNo==1
                                    pc_bar_offset=3*(concNo-1)+1;
                                    bar(pc_bar_offset,mean(these_pcorrs),'LineWidth',3,'FaceColor',[213/255,94/255,0/255],'EdgeColor','none')
                                else
                                    pc_bar_offset=3*(concNo-1)+2;
                                    bar(pc_bar_offset,mean(these_pcorrs),'LineWidth',3,'FaceColor',[86/255,180/255,233/255],'EdgeColor','none')
                                end
                                
                                %                             CI = bootci(1000, {@mean, these_pcorrs},'type','cper');
                                %                             plot([pc_bar_offset pc_bar_offset],CI,'-k','LineWidth',3)
                                %                             plot(pc_bar_offset*ones(1,length(these_pcorrs)),these_pcorrs,'o','MarkerSize',5,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                                %
                                
                                [mean_out, CIout]=drgViolinPoint(these_pcorrs...
                                    ,edges,pc_bar_offset,rand_offset,'k','k',6);
                                
                                %Save glm
                                glm_pc.data(glmpc_ii+1:glmpc_ii+length(these_pcorrs))=these_pcorrs;
                                glm_pc.group(glmpc_ii+1:glmpc_ii+length(these_pcorrs))=groupNo*ones(1,length(these_pcorrs));
                                glm_pc.rewarded(glmpc_ii+1:glmpc_ii+length(these_pcorrs))=groupNo*ones(1,length(these_pcorrs));
                                glm_pc.concs(glmpc_ii+1:glmpc_ii+length(these_pcorrs))=log10(all_concs(concNo))*ones(1,length(these_pcorrs));
                                glmpc_ii=glmpc_ii+length(these_pcorrs);
                                
                                iipc_stats=iipc_stats+1;
                                pc_stats(iipc_stats).data=these_pcorrs;
                                pc_stats(iipc_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                                    num2str(all_concs(concNo))];
                            end
                        end
                    end
                    
                end
                
                peak_bar_offfset=peak_bar_offfset+3;
                trough_bar_offfset=trough_bar_offfset+3;
            end
            
            hFig=figure(figNo_peak);
            title(['AUC peak for theta/' handles_out.drgbchoices.PACnames{PACii}])
            ylabel('AUC')
            xticks([0 1 2 5 6 7])
            xticklabels({'H Sh','H N','H P','L Sh','L N','L P'})
            ylim([-0.1 0.8])
            
            
            hFig=figure(figNo_peak+1);
            title(['AUC trough for theta/' handles_out.drgbchoices.PACnames{PACii}])
            ylabel('AUC')
            xticks([0 1 2 5 6 7])
            xticklabels({'H Sh','H N','H P','L Sh','L N','L P'})
            ylim([-0.1 0.8])
            
            figNo=figNo+2;
            

            %Perform the glm for peak
            fprintf(1, ['\n\nglm for area under the curve for both peak and trough Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            tbl = table(glm_auc.data',glm_auc.pcorr',glm_auc.rewarded',glm_auc.peak',...
                'VariableNames',{'auc','naive_prof_sh','rewarded_stimulus','peak_trough'});
            mdl = fitglm(tbl,'auc~naive_prof_sh+rewarded_stimulus+peak_trough+peak_trough*naive_prof_sh*rewarded_stimulus'...
                ,'CategoricalVars',[2,3,4])
            
            %Perform the glm for  peak
            fprintf(1, ['\n\nglm for area under the curve for peak Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            tbl = table(glm_aucp.data',glm_aucp.pcorr',glm_aucp.rewarded',...
                'VariableNames',{'auc','naive_prof_sh','rewarded_stimulus'});
            mdl = fitglm(tbl,'auc~naive_prof_sh+rewarded_stimulus+naive_prof_sh*rewarded_stimulus'...
                ,'CategoricalVars',[2,3])
            
            %Do ranksum/t test
            fprintf(1, ['\n\nRanksum or t-test p values for auc for peak for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            try
                [output_data] = drgMutiRanksumorTtest(aucp_stats);
                fprintf(1, '\n\n')
            catch
            end
            
            %Perform the glm for trough
            fprintf(1, ['\n\nglm for area under the curve for trough Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            tbl = table(glm_auct.data',glm_auct.pcorr',glm_auct.rewarded',...
                'VariableNames',{'auc','naive_prof_sh','rewarded_stimulus'});
            mdl = fitglm(tbl,'auc~naive_prof_sh+rewarded_stimulus+naive_prof_sh*rewarded_stimulus'...
                ,'CategoricalVars',[2,3])
            
            %Do ranksum/t test
            fprintf(1, ['\n\nRanksum or t-test p values for auc for trough for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            try
                [output_data] = drgMutiRanksumorTtest(auct_stats);
                fprintf(1, '\n\n')
            catch
            end
            
            %Format the pecnt correct peak and trough PRP figures
            hFig=figure(figNo_pconc);
            title(['Performance PRP peak for theta/' handles_out.drgbchoices.PACnames{PACii}])
            ylabel('Performance (%)')
            xticks([1.5 4.5 7.5 10.5 13.5 16.5])
            xticklabels({'0.032','0.1','0.32','1','3.2','10'})
            xlabel('Percent dilution')
            ylim([0 120])
            
            
            hFig=figure(figNo_pconc+1);
            title(['Performance PRP trough for theta/' handles_out.drgbchoices.PACnames{PACii}])
            ylabel('Performence (%)')
            xticks([1.5 4.5 7.5 10.5 13.5 16.5])
            xticklabels({'0.032','0.1','0.32','1','3.2','10'})
            xlabel('Percent dilution')
            ylim([0 120])
            
            figNo=figNo+2;
            
            
            %Perform the glm for performance PRP peak
            fprintf(1, ['\n\nglm for performance PRP peak Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            tbl = table(glm_pk.data',glm_pk.rewarded',glm_pk.concs',...
                'VariableNames',{'percent_correct','rewarded_stimulus','concentration'});
            mdl = fitglm(tbl,'percent_correct~rewarded_stimulus+concentration+rewarded_stimulus*concentration'...
                ,'CategoricalVars',[2])
            
            
            %Do ranksum/t test
            fprintf(1, ['\n\nRanksum or t-test p values for performance PRP peak\n'])
            try
                [output_data] = drgMutiRanksumorTtest(pk_stats);
                fprintf(1, '\n\n')
            catch
            end
            
            
            %Perform the glm for performance PRP trough
            fprintf(1, ['\n\nglm for performance PRP trough Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            tbl = table(glm_th.data',glm_th.rewarded',glm_th.concs',...
                'VariableNames',{'percent_correct','rewarded_stimulus','concentration'});
            mdl = fitglm(tbl,'percent_correct~rewarded_stimulus+concentration+rewarded_stimulus*concentration'...
                ,'CategoricalVars',[2])
            
            %Do ranksum/t test
            fprintf(1, ['\n\nRanksum or t-test p values for performance PRP trough\n'])
            try
                [output_data] = drgMutiRanksumorTtest(th_stats);
                fprintf(1, '\n\n')
            catch
            end
            
            
            if PACii==1
                %Format the perecent correct behavior figure
                hFig=figure(figNo_pc);
                hold on
                
                title(['Behavioral percent correct'])
                ylabel('Percent correct')
                xticks([1.5 4.5 7.5 10.5 13.5 16.5])
                xticklabels({'0.032','0.1','0.32','1','3.2','10'})
                xlabel('Percent dilution')
                ylim([0 120])
                
                figNo=figNo+1;
                
                %Perform the glm for behavioral percent correct
                fprintf(1, ['\n\nglm for behavioral percent correct\n'])
                tbl = table(glm_pc.data',glm_pc.rewarded',glm_pc.concs',...
                    'VariableNames',{'percent_correct','rewarded_stimulus','concentration'});
                mdl = fitglm(tbl,'percent_correct~rewarded_stimulus+concentration+rewarded_stimulus*concentration'...
                    ,'CategoricalVars',[2])
                
                
                %Do ranksum/t test
                fprintf(1, ['\n\nRanksum or t-test p values for for behavioral percent correct\n'])
                try
                    [output_data] = drgMutiRanksumorTtest(pc_stats);
                    fprintf(1, '\n\n')
                catch
                end
            end
            pffft=1;
            
            %Make a list of all per mouse percent correct data
            gr1_mice=zeros(1,size(pk_mice,2));
            gr1_mice(1,:)=pk_mice(1,:);
            
            gr2_mice=zeros(1,size(pk_mice,2));
            gr2_mice(1,:)=pk_mice(2,:);
            
            all_pk=[];
            all_beh=[];
            all_th=[];
            ii_all=0;
            
            no_m_inc=0;
            for mouseNo=1:length(handles_out.discriminant_PACwavepower)
                %Is this mouse found in both groups
                if (~isempty(find(gr1_mice==mouseNo)))&(~isempty(find(gr2_mice==mouseNo)))
                    no_m_inc=no_m_inc+1;
                   ii_gr1=find(gr1_mice==mouseNo);
                   ii_gr2=find(gr2_mice==mouseNo);
                    for concNo=1:6
                        ii_all=ii_all+1;
                        
                        all_pk(1,ii_all)=pcorrs_pk(1,concNo,ii_gr1);
                        all_pk(2,ii_all)=pcorrs_pk(2,concNo,ii_gr2);
                        
                        all_th(1,ii_all)=pcorrs_th(1,concNo,ii_gr1);
                        all_th(2,ii_all)=pcorrs_th(2,concNo,ii_gr2);
                        
                        all_beh(1,ii_all)=pcorrs_beh(1,concNo,ii_gr1);
                        all_beh(2,ii_all)=pcorrs_beh(2,concNo,ii_gr2);
                    end
                end
                
            end
            
            fprintf(1, ['No of mice included for the decpodng performance vs. behavioral percent correct %d\n\n'], no_m_inc)
            
            %Plot the peak data
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            hold on
            plot(all_beh(1,:),all_pk(1,:),'o','MarkerFaceColor',[213/255,94/255,0/255],'MarkerEdgeColor',[213/255,94/255,0/255],'MarkerSize',8)
            plot(all_beh(2,:),all_pk(2,:),'o','MarkerFaceColor',[86/255,180/255,233/255],'MarkerEdgeColor',[86/255,180/255,233/255],'MarkerSize',8)
            X=[all_beh(1,:) all_beh(2,:)]';
            Y=[all_pk(1,:) all_pk(2,:)]';
            fit_all=nlinfit(X,Y,@dr_fitline,[0 1]);
            plot([0 100],dr_fitline(fit_all,[0 100]),'-k','LineWidth',2)

            title(['Performance PRP peak for theta/' handles_out.drgbchoices.PACnames{PACii} ' vs. behavior percent correct'])
            ylabel('Performance (%)')
            xlabel('Behavior percent correect')
            ylim([0 100])
            xlim([0 100])
            
            [rho,pval] = corr(X,Y);
            
            fprintf(1, ['\n\n'])
            fprintf(1, ['Correlation performance PRP peak for theta/' handles_out.drgbchoices.PACnames{PACii} ' vs. behavior percent correct\n'])
            fprintf(1, ['rho = %d\n'],rho)
            fprintf(1, ['pval = %d\n'],pval)
            fprintf(1, ['\n\n'])
            
            %Plot the trough data
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            hold on
            plot(all_beh(1,:),all_th(1,:),'o','MarkerFaceColor',[213/255,94/255,0/255],'MarkerEdgeColor',[213/255,94/255,0/255],'MarkerSize',8)
            plot(all_beh(2,:),all_th(2,:),'o','MarkerFaceColor',[86/255,180/255,233/255],'MarkerEdgeColor',[86/255,180/255,233/255],'MarkerSize',8)
            X=[all_beh(1,:) all_beh(2,:)]';
            Y=[all_th(1,:) all_th(2,:)]';
            fit_all=nlinfit(X,Y,@dr_fitline,[0 1]);
            plot([0 100],dr_fitline(fit_all,[0 100]),'-k','LineWidth',2)

            title(['Performance PRP trough for theta/' handles_out.drgbchoices.PACnames{PACii} ' vs. behavior percent correct'])
            ylabel('Performance (%)')
            xlabel('Behavior percent correect')
            ylim([0 100])
            xlim([0 100])
            
            [rho,pval] = corr(X,Y);
            
            fprintf(1, ['\n\n'])
            fprintf(1, ['Correlation performance PRP trough for theta/' handles_out.drgbchoices.PACnames{PACii} ' vs. behavior percent correct\n'])
            fprintf(1, ['rho = %d\n'],rho)
            fprintf(1, ['pval = %d\n'],pval)
            fprintf(1, ['\n\n'])
                            
                    
        end
        
        %         output_name=[pname 'dim_' fname];
        %         save(output_name,'handles_out2','-v7.3')
        
        
        %         %Plot PCA results
        %         for groupNo=1:max(handles_out.drgbchoices.group_no)
        %
        %             for percent_correct_ii=1:2
        %                 for pca_ii=1:3
        %
        %                     %Clear variables for anovan
        %                     figNo=figNo+1;
        %                     try
        %                         close(figNo)
        %                     catch
        %                     end
        %                     hFig=figure(figNo);
        %
        %
        %                     %                     set(hFig, 'units','normalized','position',[.2 .2 .7 .7])
        %
        %
        %
        %
        %
        %                     no_mice=0;
        %                     all_PC1speak=zeros(max(handles_out.drgbchoices.mouse_no),length(handles_out.drgbchoices.events_to_discriminate),length(t));
        %                     all_PC1strough=zeros(max(handles_out.drgbchoices.mouse_no),length(handles_out.drgbchoices.events_to_discriminate),length(t));
        %                     all_varspeak=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
        %                     all_varsstrough=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
        %
        %                     for mouseNo=1:length(handles_out.discriminant_PACwavepower)
        %
        %                         try
        %                             if handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PCA_calculated==1
        %                                 per_ii=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).no_trials;
        %                                 if per_ii>=20
        %
        %                                     %Show PCA
        %                                     principal_components_peak=[];
        %                                     principal_components_peak=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).principal_components_peak;
        %                                     principal_components_trough=[];
        %                                     principal_components_trough=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).principal_components_trough;
        %                                     N=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).no_trials;
        %                                     these_all_which_events=zeros(length(handles_out.drgbchoices.events_to_discriminate),N);
        %                                     for ii=1:length(handles_out.drgbchoices.events_to_discriminate)
        %                                         kk=find(handles_out.drgbchoices.evTypeNos==handles_out.drgbchoices.events_to_discriminate(ii));
        %                                         these_all_which_events(ii,:)= handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(kk,:);
        %                                     end
        %
        %                                     %                                 %PCA before odor on
        %                                     %                                 szpcs=size(principal_components);
        %                                     %                                 these_pcs=zeros(N,szpcs(3));
        %                                     %                                 these_pcs(:,:)=principal_components(6,:,:);
        %                                     %                                 subplot(2,2,3);
        %                                     %                                 hold on
        %                                     %                                 plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'.r')
        %                                     %                                 plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'.b')
        %                                     %                                 legend(handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(1)},handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(2)})
        %                                     %                                 xlabel('PC1')
        %                                     %                                 ylabel('PC2')
        %                                     %                                 title('-1 sec')
        %                                     %
        %                                     %                                 %PCA after odor on
        %                                     %                                 subplot(2,2,4)
        %                                     %                                 hold on
        %                                     %
        %                                     %                                 szpcs=size(principal_components);
        %                                     %                                 these_pcs=zeros(N,szpcs(3));
        %                                     %                                 these_pcs(:,:)=principal_components(41,:,:);
        %                                     %                                 plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'.r')
        %                                     %                                 plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'.b')
        %                                     %                                 legend(handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(1)},handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(2)})
        %                                     %                                 xlabel('PC1')
        %                                     %                                 ylabel('PC2')
        %                                     %                                 title('2.5 sec')
        %
        %
        %                                     no_mice=no_mice+1;
        %
        %                                     all_PC1speak(no_mice,1,:)=mean(principal_components_peak(:,logical(these_all_which_events(1,:)),1),2);
        %                                     all_PC1speak(no_mice,2,:)=mean(principal_components_peak(:,logical(these_all_which_events(2,:)),1),2);
        %                                     all_PC1strough(no_mice,1,:)=mean(principal_components_trough(:,logical(these_all_which_events(1,:)),1),2);
        %                                     all_PC1strough(no_mice,2,:)=mean(principal_components_trough(:,logical(these_all_which_events(2,:)),1),2);
        %                                     sum_var_peak=sum(handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).PC_variance_peak,2);
        %                                     all_PCvarspeak(no_mice,:)=100*handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).PC_variance_peak(:,1)./sum_var_peak;
        %                                     sum_var_trough=sum(handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).PC_variance_trough,2);
        %                                     all_PCvarstrough(no_mice,:)=100*handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).PC_variance_trough(:,1)./sum_var_trough;
        %                                 end
        %                             end
        %                         catch
        %                         end
        %                     end
        %
        %                     no_mice_per(percent_correct_ii)=no_mice;
        %
        %                     all_PC1speak=all_PC1speak(1:no_mice,:,:);
        %                     all_PC1strough=all_PC1strough(1:no_mice,:,:);
        %
        %                     %Show the timecourse for PC1 peak
        %                     subplot(2,2,1)
        %                     hold on
        %
        %                     %Event 2
        %                     PC1ev2=zeros(no_mice,length(t));
        %                     PC1ev2(:,:)=all_PC1speak(:,2,:);
        %                     if size(PC1ev2,1)>2
        %                         mean_PC1ev2=mean(PC1ev2,1);
        %                         CIPC1ev2 = bootci(1000, {@mean, PC1ev2});
        %                         maxCIPC1ev2=max(CIPC1ev2(:));
        %                         minCIPC1ev2=min(CIPC1ev2(:));
        %                         CIPC1ev2(1,:)=mean_PC1ev2-CIPC1ev2(1,:);
        %                         CIPC1ev2(2,:)=CIPC1ev2(2,:)-mean_PC1ev2;
        %                         [hlCR, hpCR] = boundedline(t',mean_PC1ev2', CIPC1ev2', 'b');
        %
        %                         PC1ev1=zeros(no_mice,length(t));
        %                         PC1ev1(:,:)=all_PC1speak(:,1,:);
        %                         mean_PC1ev1=mean(PC1ev1,1);
        %                         CIPC1ev1 = bootci(1000, {@mean, PC1ev1});
        %                         maxCIPC1ev1=max(CIPC1ev1(:));
        %                         minCIPC1ev1=min(CIPC1ev1(:));
        %                         CIPC1ev1(1,:)=mean_PC1ev1-CIPC1ev1(1,:);
        %                         CIPC1ev1(2,:)=CIPC1ev1(2,:)-mean_PC1ev1;
        %                         [hlCR, hpCR] = boundedline(t',mean_PC1ev1', CIPC1ev1', 'r');
        %
        %                         maxPC1=max([maxCIPC1ev2 maxCIPC1ev1]);
        %                         minPC1=min([minCIPC1ev2 minCIPC1ev1]);
        %
        %                         %Odor on markers
        %                         plot([0 0],[minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)],'-k')
        %                         odorhl=plot([0 2.5],[minPC1-0.1*(maxPC1-minPC1) minPC1-0.1*(maxPC1-minPC1)],'-k','LineWidth',5);
        %                         plot([2.5 2.5],[minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)],'-k')
        %
        %                         xlim([-2 5])
        %                         ylim([minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)])
        %                         text(-1,minPC1+0.9*(maxPC1-minPC1),handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(1)},'Color','r')
        %                         text(-1,minPC1+0.8*(maxPC1-minPC1),handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(2)},'Color','b')
        %                         title(['PC1 peak for Theta/' handles_out.drgbchoices.PACnames{pca_ii} ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' ' handles_out.drgbchoices.group_no_names{groupNo}])
        %                         xlabel('Time (sec)')
        %                         ylabel('PC1 peak')
        %                     end
        %
        %                     %Show the timecourse for PC1 trough
        %                     subplot(2,2,2)
        %                     hold on
        %
        %                     %Event 2
        %                     PC1ev2=zeros(no_mice,length(t));
        %                     PC1ev2(:,:)=all_PC1strough(:,2,:);
        %                     if size(PC1ev2,1)>2
        %                         mean_PC1ev2=mean(PC1ev2,1);
        %                         CIPC1ev2 = bootci(1000, {@mean, PC1ev2});
        %                         maxCIPC1ev2=max(CIPC1ev2(:));
        %                         minCIPC1ev2=min(CIPC1ev2(:));
        %                         CIPC1ev2(1,:)=mean_PC1ev2-CIPC1ev2(1,:);
        %                         CIPC1ev2(2,:)=CIPC1ev2(2,:)-mean_PC1ev2;
        %                         [hlCR, hpCR] = boundedline(t',mean_PC1ev2', CIPC1ev2', 'b');
        %
        %                         PC1ev1=zeros(no_mice,length(t));
        %                         PC1ev1(:,:)=all_PC1strough(:,1,:);
        %                         mean_PC1ev1=mean(PC1ev1,1);
        %                         CIPC1ev1 = bootci(1000, {@mean, PC1ev1});
        %                         maxCIPC1ev1=max(CIPC1ev1(:));
        %                         minCIPC1ev1=min(CIPC1ev1(:));
        %                         CIPC1ev1(1,:)=mean_PC1ev1-CIPC1ev1(1,:);
        %                         CIPC1ev1(2,:)=CIPC1ev1(2,:)-mean_PC1ev1;
        %                         [hlCR, hpCR] = boundedline(t',mean_PC1ev1', CIPC1ev1', 'r');
        %
        %                         %                         maxPC1=max([maxCIPC1ev2 maxCIPC1ev1]);
        %                         %                         minPC1=min([minCIPC1ev2 minCIPC1ev1]);
        %
        %                         %Odor on markers
        %                         plot([0 0],[minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)],'-k')
        %                         odorhl=plot([0 2.5],[minPC1-0.1*(maxPC1-minPC1) minPC1-0.1*(maxPC1-minPC1)],'-k','LineWidth',5);
        %                         plot([2.5 2.5],[minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)],'-k')
        %
        %                         xlim([-2 5])
        %                         ylim([minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)])
        %                         text(-1,minPC1+0.9*(maxPC1-minPC1),handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(1)},'Color','r')
        %                         text(-1,minPC1+0.8*(maxPC1-minPC1),handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(2)},'Color','b')
        %                         title(['PC1 trough for Theta/' handles_out.drgbchoices.PACnames{pca_ii} ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' ' handles_out.drgbchoices.group_no_names{groupNo}])
        %                         xlabel('Time (sec)')
        %                         ylabel('PC1 trough')
        %                     end
        %
        %                     subplot(2,2,3)
        %
        %                     mean_all_PCvarspeak=mean(all_PCvarspeak,1);
        %                     CIall_PCvarspeak = bootci(1000, {@mean, all_PCvarspeak});
        %                     maxCIall_PCvarspeak=max(CIall_PCvarspeak(:));
        %                     minCIall_PCvarspeak=min(CIall_PCvarspeak(:));
        %                     CIall_PCvarspeak(1,:)=mean_all_PCvarspeak-CIall_PCvarspeak(1,:);
        %                     CIall_PCvarspeak(2,:)=CIall_PCvarspeak(2,:)-mean_all_PCvarspeak;
        %                     [hlCR, hpCR] = boundedline(t',mean_all_PCvarspeak', CIall_PCvarspeak', 'b');
        %
        %                     title(['% variance for PC1 peak'])
        %                     xlabel('Time (sec)')
        %                     ylabel('% Variance')
        %                     ylim([0 100])
        %
        %                     subplot(2,2,4)
        %
        %                     mean_all_PCvarstrough=mean(all_PCvarstrough,1);
        %                     CIall_PCvarstrough = bootci(1000, {@mean, all_PCvarstrough});
        %                     maxCIall_PCvarstrough=max(CIall_PCvarstrough(:));
        %                     minCIall_PCvarstrough=min(CIall_PCvarstrough(:));
        %                     CIall_PCvarstrough(1,:)=mean_all_PCvarstrough-CIall_PCvarstrough(1,:);
        %                     CIall_PCvarstrough(2,:)=CIall_PCvarstrough(2,:)-mean_all_PCvarstrough;
        %                     [hlCR, hpCR] = boundedline(t',mean_all_PCvarstrough', CIall_PCvarstrough', 'b');
        %
        %                     title(['% variance for PC1 trough'])
        %                     xlabel('Time (sec)')
        %                     ylabel('% Variance')
        %                     ylim([0 100])
        %
        %                 end
        %             end
        %             suptitle(['All trials. Group: ' handles_out.drgbchoices.group_no_names{groupNo} ' # of mice: ' num2str(no_mice_per(1)) ' ' handles_out.drgbchoices.per_lab{1} ' ' num2str(no_mice_per(2)) ' ' handles_out.drgbchoices.per_lab{2}])
        %         end
        
        %Now plot log(p) and find decision times
        edges=[0:0.1:5];
        rand_offset=0.8;
        
        
        for PACii=1:length(handles_out.drgbchoices.PACburstLowF)
            
            glm_dt=[];
            glm_ii=0;
            glm_dtp=[];
            glmp_ii=0;
            glm_dtt=[];
            glmt_ii=0;
            
            iip_stats=0;
            dtp_stats=[];
            iit_stats=0;
            dtt_stats=[];
            
            p_correct_stats=[];
            ii_stats=0;
            
            figNo_peak=figNo+5;
            
            try
                close(figNo_peak)
            catch
            end
            
            try
                close(figNo_peak+1)
            catch
            end
            
            peak_bar_offfset=0;
            trough_bar_offfset=0;
            
            
            
            for groupNo=1:max(handles_out.drgbchoices.group_no)
                
                for percent_correct_ii=1:2
                    
                    %Gather all the data
                    no_mice=0;
                    no_mice_included=0;
                    all_discriminant_p_val_lick=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    all_discriminant_p_val_peak=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    all_discriminant_p_val_trough=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    
                    for mouseNo=1:length(handles_out.discriminant_PACwavepower)
                        try
                            if handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                                per_ii=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials;
                                no_mice=no_mice+1;
                                if (per_ii>=20)&(sum(no_mice==mice_excluded)==0)
                                    no_mice_included=no_mice_included+1;
                                    all_discriminant_p_val_lick(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_lick;
                                    all_discriminant_p_val_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_peak;
                                    all_discriminant_p_val_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_trough;
                                end
                            end
                        catch
                            pffft=1;
                        end
                    end
                    
                    t_det_stats=[];
                    ii_stats=0;
                    
                    no_mice_per(percent_correct_ii)=no_mice_included;
                    
                    %Plot the p value timecourse and compute the decision time
                    figNo=figNo+1;
                    try
                        close(figNo)
                    catch
                    end
                    hFig=figure(figNo);
                    
                    hold on
                    
                    log10all_discriminant_p_val_lick=log10(all_discriminant_p_val_lick(1:no_mice_included,:));
                    log10all_discriminant_p_val_peak=log10(all_discriminant_p_val_peak(1:no_mice_included,:));
                    log10all_discriminant_p_val_trough=log10(all_discriminant_p_val_trough(1:no_mice_included,:));
                    
                    mean_p_val_licks=nanmean(log10all_discriminant_p_val_lick,1)';
                    if size(log10all_discriminant_p_val_lick,1)>2
                        CIlicks = bootci(1000, {@nanmean, log10all_discriminant_p_val_lick})';
                        CIlicks(:,1)=mean_p_val_licks-CIlicks(:,1);
                        CIlicks(:,2)=CIlicks(:,2)-mean_p_val_licks;
                        [hlick, hpCR] = boundedline(t,mean_p_val_licks, CIlicks, 'k');
                    else
                        plot(t,mean_p_val_licks,'-k')
                    end
                    
                    
                    wing=1;
                    t_detect=zeros(1,no_mice_included);
                    jj_start=find(t>=t_odor_arrival,1,'first');
                    %Find the discrimination time
                    for ii=1:no_mice_included
                        this_p_val_licks=[];
                        this_p_val_licks=all_discriminant_p_val_lick(ii,:);
                        found_disc_t=0;
                        while (found_disc_t==0)&(jj_start<length(t))
                            ii_next=find(this_p_val_licks(1,jj_start:end)<=0.05,1,'first');
                            if isempty(ii_next)
                                jj_start=length(t);
                                ii_next=1;
                            else
                                if jj_start+ii_next+wing>length(t)
                                    found_disc_t=1;
                                else
                                    if sum(this_p_val_licks(1,jj_start+ii_next:jj_start+ii_next+wing)>0.05)>=1
                                        jj_start=jj_start+ii_next;
                                    else
                                        found_disc_t=1;
                                    end
                                end
                            end
                        end
                        t_detect(ii)=t(jj_start+ii_next-1)-t_odor_arrival;
                    end
                    
                    fprintf(1, ['Lick discrimination time (sec) %d\n'],mean(t_detect))
                    
                    lick_t_detect=t_detect;
                    
                    ii_stats=ii_stats+1;
                    t_det_stats(ii_stats).data=t_detect;
                    t_det_stats(ii_stats).description='licks';
                    
                    mean_p_val_troughs=nanmean(log10all_discriminant_p_val_trough,1)';
                    if size(log10all_discriminant_p_val_trough,1)>2
                        CItroughs = bootci(1000, {@nanmean, log10all_discriminant_p_val_trough})';
                        CItroughs(:,1)=mean_p_val_troughs-CItroughs(:,1);
                        CItroughs(:,2)=CItroughs(:,2)-mean_p_val_troughs;
                        [hltrough, hpCR] = boundedline(t,mean_p_val_troughs, CItroughs, 'b');
                    else
                        plot(t,mean_p_val_troughs,'-b')
                    end
                    
                    
                    
                    t_detect=zeros(1,no_mice_included);
                    jj_start=find(t>=t_odor_arrival,1,'first');
                    %Find the discrimination time
                    for ii=1:no_mice_included
                        this_p_val_troughs=[];
                        this_p_val_troughs=all_discriminant_p_val_trough(ii,:);
                        found_disc_t=0;
                        while (found_disc_t==0)&(jj_start<length(t))
                            ii_next=find(this_p_val_troughs(1,jj_start:end)<=0.05,1,'first');
                            if isempty(ii_next)
                                jj_start=length(t);
                                ii_next=1;
                            else
                                if jj_start+ii_next+wing>length(t)
                                    found_disc_t=1;
                                else
                                    if sum(this_p_val_troughs(1,jj_start+ii_next:jj_start+ii_next+wing)>0.05)>=1
                                        jj_start=jj_start+ii_next;
                                    else
                                        found_disc_t=1;
                                    end
                                end
                            end
                        end
                        t_detect(ii)=t(jj_start+ii_next-1)-t_odor_arrival;
                    end
                    
                    trough_t_detect=t_detect;
                    
                    fprintf(1, ['Trough discrimination time (sec) %d\n'],mean(t_detect))
                    
                    ii_stats=ii_stats+1;
                    t_det_stats(ii_stats).data=t_detect;
                    t_det_stats(ii_stats).description='troughs';
                    
                    mean_p_val_peaks=nanmean(log10all_discriminant_p_val_peak,1)';
                    if size(log10all_discriminant_p_val_peak,1)>2
                        CIpeaks = bootci(1000, {@nanmean, log10all_discriminant_p_val_peak})';
                        CIpeaks(:,1)=mean_p_val_peaks-CIpeaks(:,1);
                        CIpeaks(:,2)=CIpeaks(:,2)-mean_p_val_peaks;
                        [hpeak, hpCR] = boundedline(t,mean_p_val_peaks, CIpeaks, 'r');
                    else
                        plot(t,mean_p_val_peaks,'-r')
                    end
                    
                    
                    
                    t_detect=zeros(1,no_mice_included);
                    jj_start=find(t>=t_odor_arrival,1,'first');
                    %Find the discrimination time
                    for ii=1:no_mice_included
                        this_p_val_peaks=[];
                        this_p_val_peaks=all_discriminant_p_val_peak(ii,:);
                        found_disc_t=0;
                        while (found_disc_t==0)&(jj_start<length(t))
                            ii_next=find(this_p_val_peaks(1,jj_start:end)<=0.05,1,'first');
                            if isempty(ii_next)
                                jj_start=length(t);
                                ii_next=1;
                            else
                                if jj_start+ii_next+wing>length(t)
                                    found_disc_t=1;
                                else
                                    if sum(this_p_val_peaks(1,jj_start+ii_next:jj_start+ii_next+wing)>0.05)>=1
                                        jj_start=jj_start+ii_next;
                                    else
                                        found_disc_t=1;
                                    end
                                end
                            end
                        end
                        t_detect(ii)=t(jj_start+ii_next-1)-t_odor_arrival;
                    end
                    
                    fprintf(1, ['Peak discrimination time (sec) %d\n'],mean(t_detect))
                    
                    peak_t_detect=t_detect;
                    
                    ii_stats=ii_stats+1;
                    t_det_stats(ii_stats).data=t_detect;
                    t_det_stats(ii_stats).description='peaks';
                    
                    plot(t,mean_p_val_licks,'-k')
                    plot(t,mean_p_val_troughs,'-b')
                    plot(t,mean_p_val_peaks,'-r')
                    plot([t(1) t(end)],[log10(0.05) log10(0.05)],'-r')
                    
                    
                    %Odor on markers
                    
                    odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                    
                    
                    title(['p value for Theta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    
                    xlabel('Time (sec)')
                    ylabel(['p value '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    
                    legend([hlick hltrough hpeak],{'Licks','Trough','Peak'})
                    
                    
                    
                    %Do ranksum/t test
                    fprintf(1, ['\n\nRanksum or t-test p values for area under the curve for Theta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} '\n'])
                    try
                        [output_data] = drgMutiRanksumorTtest(t_det_stats);
                        fprintf(1, '\n\n')
                    catch
                    end
%                     
%                     %Plot bar graph for decision time reationship
%                     figNo=figNo+1;
%                     try
%                         close(figNo)
%                     catch
%                     end
%                     hFig=figure(figNo);
%                     hold on
%                     
%                     plot(lick_t_detect,peak_t_detect,'or')
%                     plot(lick_t_detect,trough_t_detect,'ob')
%                     
%                     ylim_out=ylim;
%                     xlim_out=xlim;
%                     
%                     ylim([0 max([xlim_out(2) ylim_out(2)])])
%                     xlim([0 max([xlim_out(2) ylim_out(2)])])
%                     
%                     
%                     plot([0 max([xlim_out(2) ylim_out(2)])],[0 max([xlim_out(2) ylim_out(2)])],'-k')
%                     
%                     title(['Decision times for Theta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
%                     xlabel('Decision time licks (sec)')
%                     ylabel('Decision time wavelet power (sec)')
%                     legend('Peak','Trough')
                     
                    %Peak bar
                    hFig=figure(figNo_peak);
                    hold on
                    
                    %lick decision time
                    
                    %Naive or proficient
                    if percent_correct_ii==1
                        bar(peak_bar_offfset+4,mean(lick_t_detect),'LineWidth',3,'FaceColor',[0/255,114/255,168/255],'EdgeColor','none')
                        CI = bootci(1000, {@mean, lick_t_detect},'type','cper');
                        plot([peak_bar_offfset+4 peak_bar_offfset+4],CI,'-k','LineWidth',3)
                        %Violin plot
                        [mean_out, CIout]=drgViolinPoint(lick_t_detect,edges,peak_bar_offfset+4,rand_offset,'k','k',2);
                    else
                        bar(peak_bar_offfset+1,mean(lick_t_detect),'LineWidth',3,'FaceColor',[86/255,130/255,233/255],'EdgeColor','none')
                        CI = bootci(1000, {@mean, lick_t_detect},'type','cper');
                        plot([peak_bar_offfset+1 peak_bar_offfset+1],CI,'-k','LineWidth',3)
                        %Violin plot
                        [mean_out, CIout]=drgViolinPoint(lick_t_detect,edges,peak_bar_offfset+1,rand_offset,'k','k',2);
                    end
                    
                    
                    %PRP peak decision time
                    
                    %Naive or proficient
                    if percent_correct_ii==1
                        bar(peak_bar_offfset+5,mean(peak_t_detect),'LineWidth',3,'FaceColor',[204/255,121/255,167/255],'EdgeColor','none')
                        CI = bootci(1000, {@mean, peak_t_detect},'type','cper');
                        plot([peak_bar_offfset+5 peak_bar_offfset+5],CI,'-k','LineWidth',3)
                        %Violin plot
                        [mean_out, CIout]=drgViolinPoint(peak_t_detect,edges,peak_bar_offfset+5,rand_offset,'k','k',2);
                    else
                        bar(peak_bar_offfset+2,mean(peak_t_detect),'LineWidth',3,'FaceColor',[0/255,158/255,115/255],'EdgeColor','none')
                        CI = bootci(1000, {@mean, peak_t_detect},'type','cper');
                        plot([peak_bar_offfset+2 peak_bar_offfset+2],CI,'-k','LineWidth',3)
                        %Violin plot
                        [mean_out, CIout]=drgViolinPoint(peak_t_detect,edges,peak_bar_offfset+2,rand_offset,'k','k',2);
                    end
                    
                    %Save licks in glm
                    these_data=lick_t_detect;
                    glm_dt.data(glm_ii+1:glm_ii+length(these_data))=these_data;
                    glm_dt.pcorr(glm_ii+1:glm_ii+length(these_data))=percent_correct_ii*ones(1,length(these_data));
                    glm_dt.rewarded(glm_ii+1:glm_ii+length(these_data))=groupNo*ones(1,length(these_data));
                    glm_dt.peak(glm_ii+1:glm_ii+length(these_data))=ones(1,length(these_data));
                    glm_ii=glm_ii+length(these_data);
                    
                    glm_dtp.data(glmp_ii+1:glmp_ii+length(these_data))=these_data;
                    glm_dtp.pcorr(glmp_ii+1:glmp_ii+length(these_data))=percent_correct_ii*ones(1,length(these_data));
                    glm_dtp.rewarded(glmp_ii+1:glmp_ii+length(these_data))=groupNo*ones(1,length(these_data));
                    glm_dtp.lick(glmp_ii+1:glmp_ii+length(these_data))=ones(1,length(these_data));
                    glmp_ii=glmp_ii+length(these_data);
                    
                    iip_stats=iip_stats+1;
                    dtp_stats(iip_stats).data=these_data;
                    dtp_stats(iip_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                        handles_out.drgbchoices.per_lab{percent_correct_ii} ' licks'];
                    
                    %Save prp in glm
                    these_data=peak_t_detect;
                    glm_dt.data(glm_ii+1:glm_ii+length(these_data))=these_data;
                    glm_dt.pcorr(glm_ii+1:glm_ii+length(these_data))=percent_correct_ii*ones(1,length(these_data));
                    glm_dt.rewarded(glm_ii+1:glm_ii+length(these_data))=groupNo*ones(1,length(these_data));
                    glm_dt.peak(glm_ii+1:glm_ii+length(these_data))=zeros(1,length(these_data));
                    glm_ii=glm_ii+length(these_data);
                    
                    glm_dtp.data(glmp_ii+1:glmp_ii+length(these_data))=these_data;
                    glm_dtp.pcorr(glmp_ii+1:glmp_ii+length(these_data))=percent_correct_ii*ones(1,length(these_data));
                    glm_dtp.rewarded(glmp_ii+1:glmp_ii+length(these_data))=groupNo*ones(1,length(these_data));
                    glm_dtp.lick(glmp_ii+1:glmp_ii+length(these_data))=zeros(1,length(these_data));
                    glmp_ii=glmp_ii+length(these_data);
                    
                    iip_stats=iip_stats+1;
                    dtp_stats(iip_stats).data=these_data;
                    dtp_stats(iip_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                        handles_out.drgbchoices.per_lab{percent_correct_ii} ' peak'];
                    
                     
                    
                    %Trough bar
                    hFig=figure(figNo_peak+1);
                    hold on
                    
                    
                    %lick decision time
                    
                    %Naive or proficient
                    if percent_correct_ii==1
                        bar(peak_bar_offfset+4,mean(lick_t_detect),'LineWidth',3,'FaceColor',[0/255,114/255,168/255],'EdgeColor','none')
                        CI = bootci(1000, {@mean, lick_t_detect},'type','cper');
                        plot([peak_bar_offfset+4 peak_bar_offfset+4],CI,'-k','LineWidth',3)
                        %Violin plot
                        [mean_out, CIout]=drgViolinPoint(lick_t_detect,edges,peak_bar_offfset+4,rand_offset,'k','k',2);
                    else
                        bar(peak_bar_offfset+1,mean(lick_t_detect),'LineWidth',3,'FaceColor',[86/255,130/255,233/255],'EdgeColor','none')
                        CI = bootci(1000, {@mean, lick_t_detect},'type','cper');
                        plot([peak_bar_offfset+1 peak_bar_offfset+1],CI,'-k','LineWidth',3)
                        %Violin plot
                        [mean_out, CIout]=drgViolinPoint(lick_t_detect,edges,peak_bar_offfset+1,rand_offset,'k','k',2);
                    end
                    
                    
                    %PRP trough decision time
                    
                    %Naive or proficient
                    if percent_correct_ii==1
                        bar(peak_bar_offfset+5,mean(trough_t_detect),'LineWidth',3,'FaceColor',[204/255,121/255,167/255],'EdgeColor','none')
                        CI = bootci(1000, {@mean, trough_t_detect},'type','cper');
                        plot([peak_bar_offfset+5 peak_bar_offfset+5],CI,'-k','LineWidth',3)
                        %Violin plot
                        [mean_out, CIout]=drgViolinPoint(trough_t_detect,edges,peak_bar_offfset+5,rand_offset,'k','k',2);
                    else
                        bar(peak_bar_offfset+2,mean(trough_t_detect),'LineWidth',3,'FaceColor',[0/255,158/255,115/255],'EdgeColor','none')
                        CI = bootci(1000, {@mean, trough_t_detect},'type','cper');
                        plot([peak_bar_offfset+2 peak_bar_offfset+2],CI,'-k','LineWidth',3)
                        %Violin plot
                        [mean_out, CIout]=drgViolinPoint(trough_t_detect,edges,peak_bar_offfset+2,rand_offset,'k','k',2);
                    end
                    
                    %Save licks in glm
                    these_data=lick_t_detect;
                    
                    glm_dtt.data(glmt_ii+1:glmt_ii+length(these_data))=these_data;
                    glm_dtt.pcorr(glmt_ii+1:glmt_ii+length(these_data))=percent_correct_ii*ones(1,length(these_data));
                    glm_dtt.rewarded(glmt_ii+1:glmt_ii+length(these_data))=groupNo*ones(1,length(these_data));
                    glm_dtt.lick(glmt_ii+1:glmt_ii+length(these_data))=ones(1,length(these_data));
                    glmt_ii=glmt_ii+length(these_data);
                    
                    iit_stats=iit_stats+1;
                    dtt_stats(iit_stats).data=these_data;
                    dtt_stats(iit_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                        handles_out.drgbchoices.per_lab{percent_correct_ii} ' licks'];
                    
                    %Save prp in glm
                    these_data=trough_t_detect;
                    glm_dt.data(glm_ii+1:glm_ii+length(these_data))=these_data;
                    glm_dt.pcorr(glm_ii+1:glm_ii+length(these_data))=percent_correct_ii*ones(1,length(these_data));
                    glm_dt.rewarded(glm_ii+1:glm_ii+length(these_data))=groupNo*ones(1,length(these_data));
                    glm_dt.peak(glm_ii+1:glm_ii+length(these_data))=2*ones(1,length(these_data));
                    glm_ii=glm_ii+length(these_data);
                    
                    glm_dtt.data(glmt_ii+1:glmt_ii+length(these_data))=these_data;
                    glm_dtt.pcorr(glmt_ii+1:glmt_ii+length(these_data))=percent_correct_ii*ones(1,length(these_data));
                    glm_dtt.rewarded(glmt_ii+1:glmt_ii+length(these_data))=groupNo*ones(1,length(these_data));
                    glm_dtt.lick(glmt_ii+1:glmt_ii+length(these_data))=ones(1,length(these_data));
                    glmt_ii=glmt_ii+length(these_data);
                    
                    iit_stats=iit_stats+1;
                    dtt_stats(iit_stats).data=these_data;
                    dtt_stats(iit_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                        handles_out.drgbchoices.per_lab{percent_correct_ii} ' trough'];
                    
                end
                
                peak_bar_offfset=peak_bar_offfset+7;
            end
            
            hFig=figure(figNo_peak);
            title(['Decision time PRP peak for theta/' handles_out.drgbchoices.PACnames{PACii}])
            ylabel('Decision time (sec)')
            xticks([1 2 4 5 8 9 11 12])
            xticklabels({'H L N','H P N','H L P','H P P','L L N','L P N','L L P','L P P'})
            ylim([0 5.2])
            
            
            hFig=figure(figNo_peak+1);
            title(['Decision time PRP trough for theta/' handles_out.drgbchoices.PACnames{PACii}])
            ylabel('Decision time (sec)')
            xticks([1 2 4 5 8 9 11 12])
            xticklabels({'H L N','H P N','H L P','H P P','L L N','L P N','L L P','L P P'})
            ylim([0 5.2])
            
            figNo=figNo+2;
            

            %Perform the glm for both trough and peak
            fprintf(1, ['\n\nglm for decision time for lick, peak and trough Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            tbl = table(glm_dt.data',glm_dt.pcorr',glm_dt.rewarded',glm_dt.peak',...
                'VariableNames',{'decision_time','naive_prof_sh','rewarded_stimulus','peak_trough_lick'});
            mdl = fitglm(tbl,'decision_time~naive_prof_sh+rewarded_stimulus+peak_trough_lick+peak_trough_lick*naive_prof_sh*rewarded_stimulus'...
                ,'CategoricalVars',[2,3,4])
            
            %Perform the glm for  peak
            fprintf(1, ['\n\nglm for decision time for peak and lick Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            tbl = table(glm_dtp.data',glm_dtp.pcorr',glm_dtp.rewarded',glm_dtp.lick',...
                'VariableNames',{'decision_time','naive_prof_sh','rewarded_stimulus','lick_prp'});
            mdl = fitglm(tbl,'decision_time~naive_prof_sh+rewarded_stimulus+lick_prp+naive_prof_sh*rewarded_stimulus*lick_prp'...
                ,'CategoricalVars',[2,3,4])
            
            %Do ranksum/t test
            fprintf(1, ['\n\nRanksum or t-test p values for decision times for peak for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            try
                [output_data] = drgMutiRanksumorTtest(dtp_stats);
                fprintf(1, '\n\n')
            catch
            end
            
             %Perform the glm for trough
            fprintf(1, ['\n\nglm for decision time for trough and lick Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            tbl = table(glm_dtt.data',glm_dtt.pcorr',glm_dtt.rewarded',glm_dtt.lick',...
                'VariableNames',{'decision_time','naive_prof_sh','rewarded_stimulus','lick_prp'});
            mdl = fitglm(tbl,'decision_time~naive_prof_sh+rewarded_stimulus+lick_prp+naive_prof_sh*rewarded_stimulus*lick_prp'...
                ,'CategoricalVars',[2,3,4])
            
            %Do ranksum/t test
            fprintf(1, ['\n\nRanksum or t-test p values for decision times for trough for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            try
                [output_data] = drgMutiRanksumorTtest(dtt_stats);
                fprintf(1, '\n\n')
            catch
            end
            
            pffft=1;
        end
        
        pffft=1;
        

end

fprintf(1, ['Finished processing ' discriminant_name '\n'])
           


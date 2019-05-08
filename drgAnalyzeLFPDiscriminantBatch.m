function drgAnalyzeLFPDiscriminantBatch
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

close all
clear all



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

figNo=0;



%Define the windows for analysis
window_start=[1];
window_end=[2.5];
no_wins=1;

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
        
        for PACii=1:length(handles_out.drgbchoices.PACburstLowF)
            p_correct_stats=[];
            ii_stats=0;
            p_dim_stats=[];
            ii_dim_stats=0;
            glm_ii=0;
            glm_correct=[];
            glm_dim_ii=0;
            glm_dim=[];
            for percent_correct_ii=1:2
                
                for groupNo=1:max(handles_out.drgbchoices.group_no)
                    
                    %Gather all the data
                    no_mice=0;
                    no_mice_included=0;
                    all_discriminant_correct_peak=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    all_discriminant_correct_trough=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    all_discriminant_correct_shuffled_peak=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    all_discriminant_correct_shuffled_trough=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    
                    for mouseNo=1:length(handles_out.discriminant_PACwavepower)
                        try
                            if handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                                per_ii=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials;
                                no_mice=no_mice+1;
                                if (per_ii>=20)&(sum(no_mice==mice_excluded)==0)
                                    no_mice_included=no_mice_included+1;
                                    all_discriminant_correct_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_peak;
                                    all_discriminant_correct_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_trough;
                                    all_discriminant_correct_shuffled_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_peak;
                                    all_discriminant_correct_shuffled_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_trough;
                                    all_dimensionality_peak(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak;
                                    all_dimensionality_trough(no_mice_included,:)=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough;
                                end
                            end
                        catch
                        end
                    end
                    
                    
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
                    p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' shuffled peak'];
                    p_correct_stats(ii_stats).data_ii=1;
                    p_correct_stats(ii_stats).PACii=PACii;
                    p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                    p_correct_stats(ii_stats).groupNo=groupNo;
                    
                    
                    glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
                    glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
                    glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
                    glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
                    glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=1;
                    glm_correct.peak(glm_ii+1:glm_ii+length(data))=1;
                    glm_ii=glm_ii+length(data);
                    
                    data=[];
                    for mouseNo=1:no_mice_included
                        data=[data (mean(all_discriminant_correct_shuffled_trough(mouseNo,(t>=auc_from)&(t<=auc_to)),2)-50)/50];
                    end
                    
                    ii_stats=ii_stats+1;
                    p_correct_stats(ii_stats).data=data;
                    p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' shuffled trough'];
                    p_correct_stats(ii_stats).data_ii=2;
                    p_correct_stats(ii_stats).PACii=PACii;
                    p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                    p_correct_stats(ii_stats).groupNo=groupNo;
                    
                    glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
                    glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
                    glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
                    glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
                    glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=1;
                    glm_correct.peak(glm_ii+1:glm_ii+length(data))=0;
                    glm_ii=glm_ii+length(data);
                    
                    
                    
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
                    p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' trough'];
                    p_correct_stats(ii_stats).data_ii=3;
                    p_correct_stats(ii_stats).PACii=PACii;
                    p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                    p_correct_stats(ii_stats).groupNo=groupNo;
                    
                    glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
                    glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
                    glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
                    glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
                    glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=0;
                    glm_correct.peak(glm_ii+1:glm_ii+length(data))=0;
                    glm_ii=glm_ii+length(data);
                    
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
                    p_correct_stats(ii_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                        handles_out.drgbchoices.per_lab{percent_correct_ii}...
                        ' peak'];
                    p_correct_stats(ii_stats).data_ii=4;
                    p_correct_stats(ii_stats).PACii=PACii;
                    p_correct_stats(ii_stats).per_corr_ii=percent_correct_ii;
                    p_correct_stats(ii_stats).groupNo=groupNo;
                    
                    glm_correct.data(glm_ii+1:glm_ii+length(data))=data;
                    glm_correct.PACii(glm_ii+1:glm_ii+length(data))=PACii;
                    glm_correct.group(glm_ii+1:glm_ii+length(data))=groupNo;
                    glm_correct.perCorr(glm_ii+1:glm_ii+length(data))=percent_correct_ii;
                    glm_correct.shuffled(glm_ii+1:glm_ii+length(data))=0;
                    glm_correct.peak(glm_ii+1:glm_ii+length(data))=1;
                    glm_ii=glm_ii+length(data);
                    
                    %Odor on markers
                    plot([0 0],[0 100],'-k')
                    odorhl=plot([0 2.5],[10 10],'-k','LineWidth',5);
                    plot([2.5 2.5],[0 100],'-k')
                    
                    title(['LDA for teta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    
                    xlabel('Time (sec)')
                    ylabel(['% correct '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    legend('Shuffled','Trough','Peak')
                    
                    %Plot the peak vs trough graph
                    figNo=figNo+1;
                    try
                        close(figNo)
                    catch
                    end
                    hFig=figure(figNo);
                    
                    hold on
                    
                    plot(p_correct_stats(ii_stats-1).data,p_correct_stats(ii_stats).data,'or')
                    plot(p_correct_stats(ii_stats-2).data,p_correct_stats(ii_stats-3).data,'ok')
                    
                    
                    plot([-0.1 1],[-0.1 1],'-k')
                    
                    xlim([-0.1 1])
                    ylim([-0.1 1])
                    ylabel('AUC peak')
                    xlabel('AUC trough')
                    legend('Original','Shuffled')
                    title(['Area under the curve for teta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    
                    
                    
                    %Plot dimensionality and save the data for
                    %the ranksum
                    
                    winNo=1;
                    figNo=figNo+1;
                    try
                        close(figNo)
                    catch
                    end
                    hFig=figure(figNo);
                    
                    hold on
                    
                    
                    %Now plot the dimensionality  for the trough
                    all_dimensionality_trough=all_dimensionality_trough(1:no_mice_included,:);
                    mean_dc_trough=mean(all_dimensionality_trough,1)';
                    if size(all_dimensionality_trough,1)>2
                        CIdc_trough = bootci(1000, {@mean, all_dimensionality_trough})';
                        CIdc_trough(:,1)=mean_dc_trough-CIdc_trough(:,1);
                        CIdc_trough(:,2)=CIdc_trough(:,2)-mean_dc_trough;
                        [hltrough, hptrough] = boundedline(t,mean_dc_trough, CIdc_trough, 'b');
                    else
                        plot(t,mean_dc_trough,'b')
                    end
                    
                    
                    data=[];
                    for mouseNo=1:no_mice_included
                        data=[data mean(all_dimensionality_trough(mouseNo,(t>=window_start(winNo))&(t<=window_end(winNo))),2)];
                    end
                    ii_dim_stats=ii_dim_stats+1;
                    dim_stats(ii_dim_stats).data=data;
                    dim_stats(ii_dim_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                        handles_out.drgbchoices.per_lab{percent_correct_ii}...
                        ' trough'];
                    dim_stats(ii_dim_stats).winNo=winNo;
                    dim_stats(ii_dim_stats).PACii=PACii;
                    dim_stats(ii_dim_stats).per_corr_ii=percent_correct_ii;
                    dim_stats(ii_dim_stats).groupNo=groupNo;
                    
                    glm_dim.data(glm_dim_ii+1:glm_dim_ii+length(data))=data;
                    glm_dim.PACii(glm_dim_ii+1:glm_dim_ii+length(data))=PACii;
                    glm_dim.group(glm_dim_ii+1:glm_dim_ii+length(data))=groupNo;
                    glm_dim.perCorr(glm_dim_ii+1:glm_dim_ii+length(data))=percent_correct_ii;
                    glm_dim.peak(glm_dim_ii+1:glm_dim_ii+length(data))=0;
                    glm_dim_ii=glm_dim_ii+length(data);
                    
                    
                    %Now plot the percent correct for the peak
                    all_dimensionality_peak=all_dimensionality_peak(1:no_mice_included,:);
                    mean_dc_peak=mean(all_dimensionality_peak,1)';
                    if size(all_dimensionality_peak,1)>2
                        CIdc_peak = bootci(1000, {@mean, all_dimensionality_peak})';
                        CIdc_peak(:,1)=mean_dc_peak-CIdc_peak(:,1);
                        CIdc_peak(:,2)=CIdc_peak(:,2)-mean_dc_peak;
                        [hlpeak, hppeak] = boundedline(t,mean_dc_peak, CIdc_peak, 'r');
                    else
                        plot(t,mean_dc_peak,'r')
                    end
                    
                    
                    data=[];
                    for mouseNo=1:no_mice_included
                        data=[data mean(all_dimensionality_peak(mouseNo,(t>=window_start(winNo))&(t<=window_end(winNo))),2)];
                    end
                    ii_dim_stats=ii_dim_stats+1;
                    dim_stats(ii_dim_stats).data=data;
                    dim_stats(ii_dim_stats).description=[handles_out.drgbchoices.group_no_names{groupNo} ' ' ...
                        handles_out.drgbchoices.per_lab{percent_correct_ii}...
                        ' peak'];
                    dim_stats(ii_dim_stats).winNo=winNo;
                    dim_stats(ii_dim_stats).PACii=PACii;
                    dim_stats(ii_dim_stats).per_corr_ii=percent_correct_ii;
                    dim_stats(ii_dim_stats).groupNo=groupNo;
                    
                    glm_dim.data(glm_dim_ii+1:glm_dim_ii+length(data))=data;
                    glm_dim.PACii(glm_dim_ii+1:glm_dim_ii+length(data))=PACii;
                    glm_dim.group(glm_dim_ii+1:glm_dim_ii+length(data))=groupNo;
                    glm_dim.perCorr(glm_dim_ii+1:glm_dim_ii+length(data))=percent_correct_ii;
                    glm_dim.peak(glm_dim_ii+1:glm_dim_ii+length(data))=1;
                    glm_dim_ii=glm_dim_ii+length(data);
                    
                    %Odor on markers
                    this_yl=ylim;
                    plot([0 0],this_yl,'-k')
                    odorhl=plot([0 2.5],[this_yl(1)+0.1*(this_yl(2)-this_yl(1)) this_yl(1)+0.1*(this_yl(2)-this_yl(1))],'-k','LineWidth',5);
                    plot([2.5 2.5],this_yl,'-k')
                    
                    title(['Dimensionality for teta/' handles_out.drgbchoices.PACnames{PACii} ' ' handles_out.drgbchoices.group_no_names{groupNo}  ' ' handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    legend([hltrough hlpeak],{'Trough','Peak'})
                    xlabel('Time (sec)')
                    ylabel(['Dimensionality '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
                    
                    
                    pffft=1;
                end
            end
            
            %Perform the glm
            fprintf(1, ['\n\nglm for area under the curve for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            tbl = table(glm_correct.data',glm_correct.group',glm_correct.perCorr',glm_correct.shuffled',glm_correct.peak',...
                'VariableNames',{'prediction','fwd_rev','perCorr','shuffled','peak_trough'});
            mdl = fitglm(tbl,'prediction~fwd_rev+perCorr+shuffled+peak_trough+peak_trough*shuffled+peak_trough*perCorr+perCorr*shuffled+perCorr*shuffled*peak_trough'...
                ,'CategoricalVars',[2,3,4,5])
            
            %Do ranksum/t test
            fprintf(1, ['\n\nRanksum or t-test p values for area under the curve for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            try
                [output_data] = drgMutiRanksumorTtest(p_correct_stats);
                fprintf(1, '\n\n')
            catch
            end
            
          
            %Perform the glm
            fprintf(1, ['\n\nglm for dimensionality for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            tbl = table(glm_dim.data',glm_dim.group',glm_dim.perCorr',glm_dim.peak',...
                'VariableNames',{'prediction','fwd_rev','perCorr','peak_trough'});
            mdl = fitglm(tbl,'prediction~fwd_rev+perCorr+peak_trough+peak_trough*perCorr'...
                ,'CategoricalVars',[2,3,4])
            
            fprintf(1, ['\n\nRanksum or t-test p values for dimensionality for Theta/' handles_out.drgbchoices.PACnames{PACii} '\n'])
            try
                [output_data] = drgMutiRanksumorTtest(dim_stats);
                fprintf(1, '\n\n')
            catch
            end
            
           
            
            pffft=1;
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
                    
                    
                    %                     set(hFig, 'units','normalized','position',[.2 .2 .7 .7])
                    
                    
                    
                    
                    
                    no_mice=0;
                    all_PC1speak=zeros(max(handles_out.drgbchoices.mouse_no),length(handles_out.drgbchoices.events_to_discriminate),length(t));
                    all_PC1strough=zeros(max(handles_out.drgbchoices.mouse_no),length(handles_out.drgbchoices.events_to_discriminate),length(t));
                    all_varspeak=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    all_varsstrough=zeros(max(handles_out.drgbchoices.mouse_no),length(t));
                    
                    for mouseNo=1:length(handles_out.discriminant_PACwavepower)
                        
                        try
                            if handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PCA_calculated==1
                                per_ii=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).no_trials;
                                if per_ii>=20
                                    
                                    %Show PCA
                                    principal_components_peak=[];
                                    principal_components_peak=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).principal_components_peak;
                                    principal_components_trough=[];
                                    principal_components_trough=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).principal_components_trough;
                                    N=handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).no_trials;
                                    these_all_which_events=zeros(length(handles_out.drgbchoices.events_to_discriminate),N);
                                    for ii=1:length(handles_out.drgbchoices.events_to_discriminate)
                                        kk=find(handles_out.drgbchoices.evTypeNos==handles_out.drgbchoices.events_to_discriminate(ii));
                                        these_all_which_events(ii,:)= handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(kk,:);
                                    end
                                    
                                    %                                 %PCA before odor on
                                    %                                 szpcs=size(principal_components);
                                    %                                 these_pcs=zeros(N,szpcs(3));
                                    %                                 these_pcs(:,:)=principal_components(6,:,:);
                                    %                                 subplot(2,2,3);
                                    %                                 hold on
                                    %                                 plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'.r')
                                    %                                 plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'.b')
                                    %                                 legend(handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(1)},handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(2)})
                                    %                                 xlabel('PC1')
                                    %                                 ylabel('PC2')
                                    %                                 title('-1 sec')
                                    %
                                    %                                 %PCA after odor on
                                    %                                 subplot(2,2,4)
                                    %                                 hold on
                                    %
                                    %                                 szpcs=size(principal_components);
                                    %                                 these_pcs=zeros(N,szpcs(3));
                                    %                                 these_pcs(:,:)=principal_components(41,:,:);
                                    %                                 plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'.r')
                                    %                                 plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'.b')
                                    %                                 legend(handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(1)},handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(2)})
                                    %                                 xlabel('PC1')
                                    %                                 ylabel('PC2')
                                    %                                 title('2.5 sec')
                                    
                                    
                                    no_mice=no_mice+1;
                                    
                                    all_PC1speak(no_mice,1,:)=mean(principal_components_peak(:,logical(these_all_which_events(1,:)),1),2);
                                    all_PC1speak(no_mice,2,:)=mean(principal_components_peak(:,logical(these_all_which_events(2,:)),1),2);
                                    all_PC1strough(no_mice,1,:)=mean(principal_components_trough(:,logical(these_all_which_events(1,:)),1),2);
                                    all_PC1strough(no_mice,2,:)=mean(principal_components_trough(:,logical(these_all_which_events(2,:)),1),2);
                                    sum_var_peak=sum(handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).PC_variance_peak,2);
                                    all_PCvarspeak(no_mice,:)=100*handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).PC_variance_peak(:,1)./sum_var_peak;
                                    sum_var_trough=sum(handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).PC_variance_trough,2);
                                    all_PCvarstrough(no_mice,:)=100*handles_out.discriminant_PACwavepower(mouseNo).group(groupNo).percent_correct(percent_correct_ii).PACii(pca_ii).PC_variance_trough(:,1)./sum_var_trough;
                                end
                            end
                        catch
                        end
                    end
                    
                    no_mice_per(percent_correct_ii)=no_mice;
                    
                    all_PC1speak=all_PC1speak(1:no_mice,:,:);
                    all_PC1strough=all_PC1strough(1:no_mice,:,:);
                    
                    %Show the timecourse for PC1 peak
                    subplot(2,2,1)
                    hold on
                    
                    %Event 2
                    PC1ev2=zeros(no_mice,length(t));
                    PC1ev2(:,:)=all_PC1speak(:,2,:);
                    if size(PC1ev2,1)>2
                        mean_PC1ev2=mean(PC1ev2,1);
                        CIPC1ev2 = bootci(1000, {@mean, PC1ev2});
                        maxCIPC1ev2=max(CIPC1ev2(:));
                        minCIPC1ev2=min(CIPC1ev2(:));
                        CIPC1ev2(1,:)=mean_PC1ev2-CIPC1ev2(1,:);
                        CIPC1ev2(2,:)=CIPC1ev2(2,:)-mean_PC1ev2;
                        [hlCR, hpCR] = boundedline(t',mean_PC1ev2', CIPC1ev2', 'b');
                        
                        PC1ev1=zeros(no_mice,length(t));
                        PC1ev1(:,:)=all_PC1speak(:,1,:);
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
                        title(['PC1 peak for Theta/' handles_out.drgbchoices.PACnames{pca_ii} ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' ' handles_out.drgbchoices.group_no_names{groupNo}])
                        xlabel('Time (sec)')
                        ylabel('PC1 peak')
                    end
                    
                    %Show the timecourse for PC1 trough
                    subplot(2,2,2)
                    hold on
                    
                    %Event 2
                    PC1ev2=zeros(no_mice,length(t));
                    PC1ev2(:,:)=all_PC1strough(:,2,:);
                    if size(PC1ev2,1)>2
                        mean_PC1ev2=mean(PC1ev2,1);
                        CIPC1ev2 = bootci(1000, {@mean, PC1ev2});
                        maxCIPC1ev2=max(CIPC1ev2(:));
                        minCIPC1ev2=min(CIPC1ev2(:));
                        CIPC1ev2(1,:)=mean_PC1ev2-CIPC1ev2(1,:);
                        CIPC1ev2(2,:)=CIPC1ev2(2,:)-mean_PC1ev2;
                        [hlCR, hpCR] = boundedline(t',mean_PC1ev2', CIPC1ev2', 'b');
                        
                        PC1ev1=zeros(no_mice,length(t));
                        PC1ev1(:,:)=all_PC1strough(:,1,:);
                        mean_PC1ev1=mean(PC1ev1,1);
                        CIPC1ev1 = bootci(1000, {@mean, PC1ev1});
                        maxCIPC1ev1=max(CIPC1ev1(:));
                        minCIPC1ev1=min(CIPC1ev1(:));
                        CIPC1ev1(1,:)=mean_PC1ev1-CIPC1ev1(1,:);
                        CIPC1ev1(2,:)=CIPC1ev1(2,:)-mean_PC1ev1;
                        [hlCR, hpCR] = boundedline(t',mean_PC1ev1', CIPC1ev1', 'r');
                        
                        %                         maxPC1=max([maxCIPC1ev2 maxCIPC1ev1]);
                        %                         minPC1=min([minCIPC1ev2 minCIPC1ev1]);
                        
                        %Odor on markers
                        plot([0 0],[minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)],'-k')
                        odorhl=plot([0 2.5],[minPC1-0.1*(maxPC1-minPC1) minPC1-0.1*(maxPC1-minPC1)],'-k','LineWidth',5);
                        plot([2.5 2.5],[minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)],'-k')
                        
                        xlim([-2 5])
                        ylim([minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)])
                        text(-1,minPC1+0.9*(maxPC1-minPC1),handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(1)},'Color','r')
                        text(-1,minPC1+0.8*(maxPC1-minPC1),handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.events_to_discriminate(2)},'Color','b')
                        title(['PC1 trough for Theta/' handles_out.drgbchoices.PACnames{pca_ii} ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' ' handles_out.drgbchoices.group_no_names{groupNo}])
                        xlabel('Time (sec)')
                        ylabel('PC1 trough')
                    end
                    
                    subplot(2,2,3)
                    
                    mean_all_PCvarspeak=mean(all_PCvarspeak,1);
                    CIall_PCvarspeak = bootci(1000, {@mean, all_PCvarspeak});
                    maxCIall_PCvarspeak=max(CIall_PCvarspeak(:));
                    minCIall_PCvarspeak=min(CIall_PCvarspeak(:));
                    CIall_PCvarspeak(1,:)=mean_all_PCvarspeak-CIall_PCvarspeak(1,:);
                    CIall_PCvarspeak(2,:)=CIall_PCvarspeak(2,:)-mean_all_PCvarspeak;
                    [hlCR, hpCR] = boundedline(t',mean_all_PCvarspeak', CIall_PCvarspeak', 'b');
                    
                    title(['% variance for PC1 peak'])
                    xlabel('Time (sec)')
                    ylabel('% Variance')
                    ylim([0 100])
                    
                    subplot(2,2,4)
                    
                    mean_all_PCvarstrough=mean(all_PCvarstrough,1);
                    CIall_PCvarstrough = bootci(1000, {@mean, all_PCvarstrough});
                    maxCIall_PCvarstrough=max(CIall_PCvarstrough(:));
                    minCIall_PCvarstrough=min(CIall_PCvarstrough(:));
                    CIall_PCvarstrough(1,:)=mean_all_PCvarstrough-CIall_PCvarstrough(1,:);
                    CIall_PCvarstrough(2,:)=CIall_PCvarstrough(2,:)-mean_all_PCvarstrough;
                    [hlCR, hpCR] = boundedline(t',mean_all_PCvarstrough', CIall_PCvarstrough', 'b');
                    
                    title(['% variance for PC1 trough'])
                    xlabel('Time (sec)')
                    ylabel('% Variance')
                    ylim([0 100])
                    
                end
            end
            suptitle(['All trials. Group: ' handles_out.drgbchoices.group_no_names{groupNo} ' # of mice: ' num2str(no_mice_per(1)) ' ' handles_out.drgbchoices.per_lab{1} ' ' num2str(no_mice_per(2)) ' ' handles_out.drgbchoices.per_lab{2}])
        end
        

end




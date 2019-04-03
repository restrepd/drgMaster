function drgAnalyzeLFPDiscriminantBatch
%Analyzes the linear discriminant analysis performed by drgLFPDiscriminantBatch
%Takes as a in input the 'Discriminant_*.mat' output file from drgLFPDiscriminantBatch
%Performs an analysis of the timecourse for percent correct for LDA and for
%the PCA

close all
clear all

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
 
t=handles_out.t;

%Define the windows for analysis
window_start=[1];
window_end=[2.5];
no_wins=1;



%Plot average percent correct for the LDA
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

pffft=1;


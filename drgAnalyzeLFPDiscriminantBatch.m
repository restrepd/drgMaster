function drgAnalyzeLFPDiscriminantBatch

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

%Plot average percent correct for the LDA
for groupNo=1:max(handles_out.drgbchoices.group_no)
    figNo=figNo+1
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    
    set(hFig, 'units','normalized','position',[.2 .2 .7 .7])
    
     
    ii_plot=0;
    for percent_correct_ii=1:2
        for bwii=1:4
            
            %Clear variables for anovan
            
            
            
            ii_plot=ii_plot+1;
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
            
            all_discriminant_correct=all_discriminant_correct(1:no_mice,:);
            mean_dc=mean(all_discriminant_correct,1)';
            CIdc = bootci(1000, {@mean, all_discriminant_correct})';
            CIdc(:,1)=mean_dc-CIdc(:,1);
            CIdc(:,2)=CIdc(:,2)-mean_dc;
            [hlCR, hpCR] = boundedline(t,mean_dc, CIdc, 'r');
            
            
            %Odor on markers
            plot([0 0],[0 100],'-k')
            odorhl=plot([0 2.5],[20 20],'-k','LineWidth',5);
            plot([2.5 2.5],[0 100],'-k')
            
            title([handles_out.drgbchoices.bwlabels{bwii} ])
            
            xlabel('Time (sec)')
            ylabel(['% correct '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
            
            
        end
    end
    suptitle(['All trials. Group: ' handles_out.drgbchoices.group_no_names{groupNo} ' # of mice: ' num2str(no_mice_per(1)) ' ' handles_out.drgbchoices.per_lab{1} ' ' num2str(no_mice_per(2)) ' ' handles_out.drgbchoices.per_lab{2}])
end

%Do anovan
percent_correct_ii=1;
t_from=t(1);
t_to=t(end);
for bwii=1:4
    %Clear variables
    a_time=[];
    a_group=[];
    a_dcapc=[];
    
    for groupNo=1:max(handles_out.drgbchoices.group_no)
        
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
                    a_time=[a_time t((t>=t_from)&(t<=t_to))'];
                    a_group=[a_group groupNo*ones(1,length(t((t>=t_from)&(t<=t_to))'))];
                    a_dcapc=[a_dcapc handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).discriminant_correct((t>=t_from)&(t<=t_to))];
                end
            end
            
            
        end
    end
    %Calculate anovan for inteaction
    [p,tbl,stats]=anovan(a_dcapc,{a_group a_time},'varnames',{'groups','times'},'display','off');
    fprintf(1, ['p value for anovan LDA percent correct for groups for ' handles_out.drgbchoices.bwlabels{bwii} '= %d \n'],  p(1));
    fprintf(1, ['p value for anovan LDA percent correct for time for ' handles_out.drgbchoices.bwlabels{bwii} '= %d \n'],  p(2));
    
    
end

%Now do the analysis for each event

%Which tested event do all events belong to? For example Hit and Miss
%belong to S+
compare_to_event=[5 5 5 11 11 11];

if ~isfield(handles_out,'drg')
    handles_out.drg.draq_d.eventlabels{3}='Hit';
    handles_out.drg.draq_d.eventlabels{5}='S+';
    handles_out.drg.draq_d.eventlabels{7}='Miss';
    handles_out.drg.draq_d.eventlabels{9}='CR';
    handles_out.drg.draq_d.eventlabels{11}='S-';
    handles_out.drg.draq_d.eventlabels{13}='FA';
end



for evNo=1:length(handles_out.drgbchoices.evTypeNos)
    compare_event_ii=find(handles_out.drgbchoices.events_to_discriminate==compare_to_event(evNo));
    
    for groupNo=1:max(handles_out.drgbchoices.group_no)
        figNo=figNo+1
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        
        
        set(hFig, 'units','normalized','position',[.2 .2 .7 .7])
        
        
        ii_plot=0;
        for percent_correct_ii=1:2
            for bwii=1:4
                ii_plot=ii_plot+1;
                subplot(2,4,ii_plot);
                hold on
                
                decoding_targets=zeros(10000,1);
                test_out=zeros(10000,length(t));
                no_trials=0;
                
                for mouseNo=1:max(handles_out.drgbchoices.mouse_no)
                    
                    if handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                        per_ii=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).no_trials;
                        if per_ii>=20
                            these_events=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(evNo,:);
                            if bwii==1
                                fprintf(1, ['For ' handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.evTypeNos(evNo)} ', mouse no %d ' handles_out.drgbchoices.group_no_names{groupNo} ' ' handles_out.drgbchoices.per_lab{percent_correct_ii} ' there are %d trials\n'], mouseNo, sum(these_events));
                            end
                            decoding_targets(no_trials+1:no_trials+sum(these_events),1)=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).per_targets(compare_event_ii,logical(these_events));
                            for tt=1:length(t)
                                test_out(no_trials+1:no_trials+sum(these_events),tt)=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).test_out_per_timepoint(compare_event_ii,logical(these_events),tt);
                            end
                            no_trials=no_trials+sum(these_events);
                        end
                    end
                end
                
                percent_correct=zeros(1,length(t));
                for tt=1:length(t)
                    percent_correct(tt)=100*sum(test_out(1:no_trials,tt).*decoding_targets(1:no_trials,1))/sum(decoding_targets);
                end
                plot(t,percent_correct,'-r')
                
                %Try to bootstrap a CI here
                
                
                %Odor on markers
                plot([0 0],[0 100],'-k')
                odorhl=plot([0 2.5],[20 20],'-k','LineWidth',5);
                plot([2.5 2.5],[0 100],'-k')
                
                title([handles_out.drgbchoices.bwlabels{bwii} ])
                
                xlabel('Time (sec)')
                ylabel(['% correct '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
            end
        end
        
        suptitle([handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.evTypeNos(evNo)} ' Group: ' handles_out.drgbchoices.group_no_names{groupNo} ' # of mice: ' num2str(no_mice_per(1)) ' ' handles_out.drgbchoices.per_lab{1} ' ' num2str(no_mice_per(2)) ' ' handles_out.drgbchoices.per_lab{2}])
        
    end
    
    
end

for evNo=1:length(handles_out.drgbchoices.evTypeNos)
    compare_event_ii=find(handles_out.drgbchoices.events_to_discriminate==compare_to_event(evNo));
    
    for groupNo=1:max(handles_out.drgbchoices.group_no)
        figNo=figNo+1
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        
        
        set(hFig, 'units','normalized','position',[.2 .2 .7 .7])
        
        
        ii_plot=0;
        for percent_correct_ii=1:2
            bwii=1;
            ii_plot=ii_plot+1;
            subplot(2,1,ii_plot);
            hold on
            
            lick_rate_per_trial=zeros(10000,length(t));
            no_trials=0;
            
            for mouseNo=1:max(handles_out.drgbchoices.mouse_no)
                
                if handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated==1
                    per_ii=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).no_trials;
                    if per_ii>=20
                        these_events=handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events(evNo,:);
                        for eventNo=1:length(these_events)
                            if these_events(eventNo)==1
                                no_trials=no_trials+1;
                                for lick_ii=1:handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii(eventNo)
                                    this_lick_t= handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times(eventNo,lick_ii);
                                    dt=t(2)-t(1);
                                    if (this_lick_t>=t(1)-dt/2)&(this_lick_t<=t(end)+dt/2)
                                        tt=ceil((this_lick_t-t(1))/dt);
                                        if (tt>=1)&(tt<=length(t))
                                            lick_rate_per_trial(no_trials,tt)=lick_rate_per_trial(no_trials,tt)+(1/dt);
                                        end
                                    end
                                end
                            end
                        end
                        
                    end
                end
            end
            
            try
                mean_lick_rate=mean(lick_rate_per_trial(1:no_trials,:),1)';
                CImlr = bootci(1000, {@mean, lick_rate_per_trial(1:no_trials,:)})';
                CImlr(:,1)=mean_lick_rate-CImlr(:,1);
                CImlr(:,2)=CImlr(:,2)-mean_lick_rate;
                [hlCR, hpCR] = boundedline(t,mean_lick_rate, CImlr, 'r');
                
            catch
            end
            %Odor on markers
            plot([0 0],[0 10],'-k')
            odorhl=plot([0 2.5],[1 1],'-k','LineWidth',5);
            plot([2.5 2.5],[0 10],'-k')
            
            
            xlabel('Time (sec)')
            ylabel(['Lick rate (Hz) '  handles_out.drgbchoices.per_lab{percent_correct_ii}])
            
            
        end
        
        suptitle([handles_out.drg.draq_d.eventlabels{handles_out.drgbchoices.evTypeNos(evNo)} ' Group: ' handles_out.drgbchoices.group_no_names{groupNo} ' # of mice: ' num2str(no_mice_per(1)) ' ' handles_out.drgbchoices.per_lab{1} ' ' num2str(no_mice_per(2)) ' ' handles_out.drgbchoices.per_lab{2}])
        
    end
    
    
end


%Plot PCA results
for groupNo=1:max(handles_out.drgbchoices.group_no)
    
    
    
    
    for percent_correct_ii=1:2
        for bwii=1:4
            
            %Clear variables for anovan
            figNo=figNo+1
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


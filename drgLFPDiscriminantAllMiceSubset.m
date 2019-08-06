function drgLFPDiscriminantAllMiceSubset

%Performs post hoc LDA analysis using data from all mice


close all
clear all



[fname,pname,nCancel] = uigetfile({'Discriminant_*.mat'},'Select the drgLFPDiscriminantBatchAllMice output file ...');
if nCancel
    inputPath = [pname,fname];
    pnameStart = pname;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end

output_name=[pname 'Discriminant_2' fname(14:end)];

discriminant_name=[pname fname];
load(discriminant_name)

figNo=0;
handles=[];
handles.drgbchoices=handles_out.drgbchoices;
t_pac=handles_out.t_power';
t_power=handles_out.t_power';
handles.drg=handles_out.drg;
overwrite_out_file=0;

%Variables for all mouse processing
min_trials=20;

%You can choose to process only a subset of mice
mice_included=[1:max(handles.drgbchoices.mouse_no)];

%Used to troubleshoot the IAMO JL
% mice_included=[1 3 4];


%Calculate discriminant analysis for wavelet power calculated at
%different PAC phases for data for all mice
t=t_pac;
for groupNo=1:max(handles.drgbchoices.group_no)
    for PACii=1:length(handles.drgbchoices.PACburstLowF)
        
        for percent_correct_ii=1:length(handles.drgbchoices.per_lab)
            
            %Choose teh analysis type
            %1=use the minimum number of trials
            %2=use all available trials by processing LDA for multile batches of 30 trials
            which_allm_LDA_analysis=2;
            
            switch which_allm_LDA_analysis
                case 1
                    %First determine how many trials per mouse
                    Nmax=0;
                    Nmin=20000;
                    mouse_included=[];
                    max_mouse=[];
                    min_mouse=[];
                    these_all_licks_per_tPACwave=[];
                    these_all_stamped_lick_ii=[];
                    these_all_which_events_licks=[];
                    N_licks=0;
                    
                    for mouseNo=1:max(handles.drgbchoices.mouse_no)
                        if sum(mice_included==mouseNo)>0
                            mouse_included(mouseNo)=0;
                            try
                                if handles_out.all_mouse_wav(mouseNo).group(groupNo).mouse_has_data==1
                                    group_has_data=1;
                                else
                                    group_has_data=0;
                                end
                            catch
                                group_has_data=0;
                            end
                            if group_has_data==1
                                these_per_corr=(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave>=handles.drgbchoices.percent_windows(percent_correct_ii,1))...
                                    &(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave<=handles.drgbchoices.percent_windows(percent_correct_ii,2));
                                thisN=sum(these_per_corr)
                                this_mouseNo=mouseNo
                                
                                if sum(these_per_corr)>=min_trials
                                    mouse_included(mouseNo)=1;
                                    if Nmax<sum(these_per_corr)
                                        max_mouse=mouseNo;
                                        Nmax=sum(these_per_corr);
                                    end
                                    
                                    if Nmin>sum(these_per_corr)
                                        min_mouse=mouseNo;
                                        Nmin=sum(these_per_corr);
                                    end
                                else
                                    mouse_included(mouseNo)=0;
                                end
                                
                                
                                N_licks=N_licks+sum(these_per_corr);
                                
                            end
                        end
                    end
                    
                    %             mouseNo=max_mouse;
                    %             reference_mouse=max_mouse;
                    %             N=Nmax;
                    
                    if (groupNo==1)&(PACii==3)&(percent_correct_ii==1)
                        pffft=1;
                    end
                    
                    mouseNo=min_mouse;
                    reference_mouse=min_mouse;
                    N=Nmin;
                    
                    these_all_licks_per_tPACwave=zeros(N_licks,length(t));
                    these_all_stamped_lick_ii=zeros(1,N_licks);
                    these_all_which_events_licks=zeros(1,N_licks);
                    ii_licks=1;
                    
                    %Now find these_all_log_P_timecoursePACwave and licks
                    these_all_log_P_timecoursePACwave_peak=zeros(length(handles.drgbchoices.which_electrodes)*sum(mouse_included),N,length(t));
                    these_all_log_P_timecoursePACwave_trough=zeros(length(handles.drgbchoices.which_electrodes)*sum(mouse_included),N,length(t));
                    
                    
                    no_mice_included=0;
                    handles_out.group(groupNo).percent_correct(percent_correct_ii).max_mouse=max_mouse;
                    fprintf(1, ['\n\nNumber of trials for each mouse for ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} '\n']);
                    
                    
                    %Compute these_variables for the max_mouse with the largest number of trials
                    
                    
                    no_mice_included=no_mice_included+1;
                    first_electrode=1;
                    
                    these_per_corr=(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave>=handles.drgbchoices.percent_windows(percent_correct_ii,1))...
                        &(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave<=handles.drgbchoices.percent_windows(percent_correct_ii,2));
                    
                    these_all_log_P_timecoursePACwave_peak(first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,:,:)=...
                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_log_P_timecoursePACwavepeak(PACii,:,these_per_corr,:);
                    
                    these_all_log_P_timecoursePACwave_trough(first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,:,:)=...
                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_log_P_timecoursePACwavetrough(PACii,:,these_per_corr,:);
                    
                    %Save these data and notify user of the number
                    %of trials
                    handles_out.all_mouse_wav(mouseNo).group(groupNo).percent_correct(percent_correct_ii).no_trials=sum(these_per_corr);
                    fprintf(1, ['Number of trials for mouse No %d = %d\n'],mouseNo,sum(these_per_corr));
                    
                    
                    these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                    for ii=1:length(handles.drgbchoices.events_to_discriminate)
                        kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                        these_all_which_events(ii,:)= handles_out.all_mouse_wav(mouseNo).group(groupNo).all_which_eventsPACwave(kk,these_per_corr);
                    end
                    
                    
                    all_which_eventsPACwave=handles_out.all_mouse_wav(mouseNo).group(groupNo).all_which_eventsPACwave;
                    
                    
                    these_ii_splus=[];
                    ii_s=0;
                    these_ii_sminus=[];
                    ii_m=0;
                    for ii=1:length(these_all_which_events(1,:))
                        if these_all_which_events(1,ii)==1
                            ii_s=ii_s+1;
                            these_ii_splus(ii_s)=ii;
                        else
                            ii_m=ii_m+1;
                            these_ii_sminus(ii_m)=ii;
                        end
                    end
                    
                    these_all_licks_per_tPACwave(ii_licks:ii_licks+ size(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_licks_per_tPACwave(these_per_corr,:),1)-1,:)=...
                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_licks_per_tPACwave(these_per_corr,:);
                    
                    these_all_stamped_lick_ii(1,ii_licks:ii_licks+ length(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,these_per_corr))-1)=...
                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,these_per_corr);
                    
                    these_all_which_events_licks(1,ii_licks:ii_licks+ length(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,these_per_corr))-1)=...
                        these_all_which_events(1,:);
                    
                    ii_licks=ii_licks+length(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,these_per_corr));
                    
                    first_electrode=first_electrode+length(handles.drgbchoices.which_electrodes);
                    
                    %Now enter the data for the rest of the mice
                    for mouseNo=1:max(handles.drgbchoices.mouse_no)
                        if sum(mice_included==mouseNo)>0
                            if reference_mouse~=mouseNo
                                if mouse_included(mouseNo)==1
                                    %Enter the existing trials
                                    no_mice_included=no_mice_included+1;
                                    
                                    those_per_corr=(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave>=handles.drgbchoices.percent_windows(percent_correct_ii,1))...
                                        &(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave<=handles.drgbchoices.percent_windows(percent_correct_ii,2));
                                    
                                    those_all_log_P_timecoursePACwave_peak=zeros(length(handles.drgbchoices.which_electrodes),sum(those_per_corr),length(t));
                                    those_all_log_P_timecoursePACwave_peak(:,:,:)=...
                                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_log_P_timecoursePACwavepeak(PACii,:,those_per_corr,:);
                                    
                                    those_all_log_P_timecoursePACwave_trough=zeros(length(handles.drgbchoices.which_electrodes),sum(those_per_corr),length(t));
                                    those_all_log_P_timecoursePACwave_trough(:,:,:)=...
                                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_log_P_timecoursePACwavetrough(PACii,:,those_per_corr,:);
                                    
                                    handles_out.all_mouse_wav(mouseNo).group(groupNo).percent_correct(percent_correct_ii).no_trials=sum(those_per_corr);
                                    fprintf(1, ['Number of trials for mouse No %d = %d\n'],mouseNo,sum(those_per_corr));
                                    
                                    those_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(those_per_corr));
                                    
                                    for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                        kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                        those_all_which_events(ii,:)= handles_out.all_mouse_wav(mouseNo).group(groupNo).all_which_eventsPACwave(kk,those_per_corr);
                                    end
                                    
                                    those_ii_splus=[];
                                    ii_s=0;
                                    those_ii_sminus=[];
                                    ii_m=0;
                                    for ii=1:length(those_all_which_events(1,:))
                                        if those_all_which_events(1,ii)==1
                                            ii_s=ii_s+1;
                                            those_ii_splus(ii_s)=ii;
                                        else
                                            ii_m=ii_m+1;
                                            those_ii_sminus(ii_m)=ii;
                                        end
                                    end
                                    
                                    %First transfer S+
                                    jj=1;
                                    for ii=1:length(these_ii_splus)
                                        these_all_log_P_timecoursePACwave_peak(first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,these_ii_splus(ii),:)=...
                                            those_all_log_P_timecoursePACwave_peak(:,those_ii_splus(jj),:);
                                        these_all_log_P_timecoursePACwave_trough(first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,these_ii_splus(ii),:)=...
                                            those_all_log_P_timecoursePACwave_trough(:,those_ii_splus(jj),:);
                                        jj=jj+1;
                                        if jj>length(those_ii_splus)
                                            jj=1;
                                        end
                                    end
                                    
                                    %Then transfer S-
                                    jj=1;
                                    for ii=1:length(these_ii_sminus)
                                        these_all_log_P_timecoursePACwave_peak(first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,these_ii_sminus(ii),:)=...
                                            those_all_log_P_timecoursePACwave_peak(:,those_ii_sminus(jj),:);
                                        these_all_log_P_timecoursePACwave_trough(first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,these_ii_sminus(ii),:)=...
                                            those_all_log_P_timecoursePACwave_trough(:,those_ii_sminus(jj),:);
                                        jj=jj+1;
                                        if jj>length(those_ii_sminus)
                                            jj=1;
                                        end
                                    end
                                    
                                    these_all_licks_per_tPACwave(ii_licks:ii_licks+ size(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_licks_per_tPACwave(those_per_corr,:),1)-1,:)=...
                                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_licks_per_tPACwave(those_per_corr,:);
                                    
                                    these_all_stamped_lick_ii(1,ii_licks:ii_licks+ length(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,those_per_corr))-1)=...
                                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,those_per_corr);
                                    
                                    these_all_which_events_licks(1,ii_licks:ii_licks+ length(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,those_per_corr))-1)=...
                                        those_all_which_events(1,:);
                                    
                                    ii_licks=ii_licks+length(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,those_per_corr));
                                    
                                    first_electrode=first_electrode+length(handles.drgbchoices.which_electrodes);
                                end
                            end
                        end
                    end
                    
                    No_electrodes=first_electrode-1;
                    
                    
                    discriminant_correct=zeros(1,length(t));
                    discriminant_correct_shuffled=zeros(1,length(t));
                    auROC=zeros(1,length(t));
                    
                    
                    %Do the analysis only if there are more than 20 trials
                    if N>=20
                        
                        if (sum(handles.drgbchoices.which_discriminant==15)>0)
                            %Linear discriminant analysis for
                            %peak PAC power al mice
                            
                            
                            fprintf(1, ['\nLDA processed for peak PAC wavelet power for all mice for ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],N);
                            fprintf(1,'For these events: ')
                            for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                            end
                            fprintf(1,'\n')
                            
                            
                            par_out=[];
                            test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                            shuffled_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                            discriminant_correct=zeros(1,length(t));
                            discriminant_correct_shuffled=zeros(1,length(t));
                            auROC=zeros(1,length(t));
                            dimensionality=zeros(1,length(t));
                            per_targets=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                            per_targets=these_all_which_events;
                            
                            parfor time_point=1:length(t)
                                %                             for time_point=1:length(t)
                                %
                                %LFP power per trial per electrode
                                measurements=zeros(N,No_electrodes);
                                measurements(:,:)=these_all_log_P_timecoursePACwave_peak(:,:,time_point)';
                                %                         measurements=(measurements-repmat(mean(measurements),N,1))./repmat(std(measurements,0,1),N,1);
                                
                                %Dimensionality
                                %Rows: trials, Columns: electrodes
                                Signal=measurements;
                                dimensionality(time_point) = nansum(eig(cov(Signal)))^2/nansum(eig(cov(Signal)).^2);
                                
                                %Enter strings labeling each event (one event for
                                %each trial)
                                events=[];
                                
                                for ii=1:N
                                    this_event=zeros(length(handles.drgbchoices.events_to_discriminate),1);
                                    this_event=these_all_which_events(:,ii);
                                    
                                    events{ii,1}=handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(logical(this_event))};
                                    
                                end
                                
                                tested_events=[];
                                shuffled_tested_events=[];
                                scores=[];
                                for ii=1:N
                                    %Partition the data into training and test sets.
                                    
                                    %Create input and target vectors leaving one trial out
                                    %For per_input each column has the dF/F for one trial
                                    %each row is a single time point for dF/F for one of the cells
                                    %For per_target the top row is 1 if the odor is S+ and 0 if it is
                                    %S-, and row 2 has 1 for S-
                                    idxTrn=ones(N,1);
                                    idxTrn(ii)=0;
                                    idxTest=zeros(N,1);
                                    idxTest(ii)=1;
                                    
                                    %Store the training data in a table.
                                    tblTrn=[];
                                    tblTrn = array2table(measurements(logical(idxTrn),:));
                                    tblTrn.Y = events(logical(idxTrn));
                                    
                                    %Train a discriminant analysis model using the training set and default options.
                                    %By default this is a regularized linear discriminant analysis (LDA)
                                    Mdl = fitcdiscr(tblTrn,'Y');
                                    
                                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                                    [label,score] = predict(Mdl,measurements(logical(idxTest),:));
                                    
                                    tested_events{ii,1}=label{1};
                                    scores(ii)=score(2);
                                    
                                    %Do LDA with shuffled trials
                                    shuffled_measurements=zeros(N,No_electrodes);
                                    shuffled_measurements(:,:)=measurements(randperm(N),:);
                                    
                                    %Store the training data in a table.
                                    sh_tblTrn=[];
                                    sh_tblTrn = array2table(shuffled_measurements(logical(idxTrn),:));
                                    sh_tblTrn.Y = events(logical(idxTrn));
                                    
                                    %Train a discriminant analysis model using the training set and default options.
                                    %By default this is a regularized linear discriminant analysis (LDA)
                                    sh_Mdl = fitcdiscr(sh_tblTrn,'Y');
                                    
                                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                                    sh_label = predict(Mdl,shuffled_measurements(logical(idxTest),:));
                                    
                                    shuffled_tested_events{ii,1}=sh_label{1};
                                    
                                end
                                
                                %Calculate auROC
                                [X,Y,T,AUC] = perfcurve(events,scores',handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)});
                                auROC(1,time_point)=AUC-0.5;
                                %test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                                
                                test_out=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                                shuffled_out=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                                for ii=1:N
                                    for jj=1:length(handles.drgbchoices.events_to_discriminate)
                                        if strcmp(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(jj)},tested_events{ii})
                                            test_out(jj,ii)=1;
                                        else
                                            test_out(jj,ii)=0;
                                        end
                                        if strcmp(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(jj)},shuffled_tested_events{ii})
                                            shuffled_out(jj,ii)=1;
                                        else
                                            shuffled_out(jj,ii)=0;
                                        end
                                    end
                                end
                                
                                test_out_per_timepoint(:,:,time_point)=test_out;
                                shuffled_out_per_timepoint(:,:,time_point)=shuffled_out;
                                discriminant_correct(1,time_point)=100*sum(per_targets(1,:)==test_out(1,:))/N;
                                discriminant_correct_shuffled(1,time_point)=100*sum(per_targets(1,:)==shuffled_out(1,:))/N;
                                fprintf(1, 'LDA for peak PAC wavelet power percent correct classification %d (for timepoint %d out of %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,length(t));
                                
                            end
                            
                            
                            
                            figNo=figNo+1;
                            try
                                close(figNo)
                            catch
                            end
                            
                            hFig=figure(figNo)
                            set(hFig, 'units','normalized','position',[.1 .4 .75 .47])
                            
                            subplot(2,5,1)
                            hold on
                            
                            per95=prctile(discriminant_correct_shuffled(1,:),95);
                            per5=prctile(discriminant_correct_shuffled(1,:),5);
                            CIsh=[mean(discriminant_correct_shuffled(1,:))-per5 per95-mean(discriminant_correct_shuffled(1,:))]';
                            [hlCR, hpCR] = boundedline([t(1) t(end)],[mean(discriminant_correct_shuffled(1,:)) mean(discriminant_correct_shuffled(1,:))], CIsh', 'r');
                            
                            plot(t',discriminant_correct(1,:),'-k')
                            
                            %Odor on markers
                            plot([0 0],[0 100],'-k')
                            odorhl=plot([0 2.5],[20 20],'-k','LineWidth',5);
                            plot([2.5 2.5],[0 100],'-k')
                            
                            %title(['LDA % correct for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                            xlabel('Time (sec)')
                            ylabel('% correct peak')
                            
                            subplot(2,5,2)
                            hold on
                            
                            plot(t,auROC,'-b')
                            
                            %Odor on markers
                            plot([0 0],[-0.3 0.5],'-k')
                            odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                            plot([2.5 2.5],[0 0.5],'-k')
                            
                            %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                            xlabel('Time (sec)')
                            ylabel('auROC peak')
                            ylim([-0.3 0.6])
                            
                            
                            subplot(2,5,3)
                            hold on
                            
                            plot(t,dimensionality,'-b')
                            
                            maxdim=max(dimensionality);
                            mindim=min(dimensionality);
                            
                            %Odor on markers
                            plot([0 0],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                            odorhl=plot([0 2.5],[mindim mindim],'-k','LineWidth',5);
                            plot([2.5 2.5],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                            
                            %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                            xlabel('Time (sec)')
                            ylabel('dimensionality peak')
                            ylim([mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)])
                            
                            
                            %Calculate p value for the wavelet power
                            p_val_peak=zeros(1,length(t));
                            for ii_t=1:length(t)
                                splus_out=zeros(1,sum(these_all_which_events(1,:)==1));
                                splus_out(1,:)=test_out_per_timepoint(1,these_all_which_events(1,:)==1,ii_t);
                                sminus_out=zeros(1,sum(these_all_which_events(1,:)==0));
                                sminus_out(1,:)=test_out_per_timepoint(1,these_all_which_events(1,:)==0,ii_t);
                                p_val_peak(ii_t)=ranksum(splus_out,sminus_out);
                            end
                            
                            
                            %Calculate p value for the licks
                            p_val_lick=zeros(1,length(t));
                            for ii_t=1:length(t)
                                splus_licks=zeros(1,sum((these_all_which_events_licks==1)&(these_all_stamped_lick_ii>0)));
                                splus_licks(1,:)=these_all_licks_per_tPACwave((these_all_which_events_licks==1)&(these_all_stamped_lick_ii>0),ii_t);
                                sminus_licks=zeros(1,sum((these_all_which_events_licks==0)&(these_all_stamped_lick_ii>0)));
                                sminus_licks(1,:)=these_all_licks_per_tPACwave((these_all_which_events_licks==0)&(these_all_stamped_lick_ii>0),ii_t);
                                p_val_lick(ii_t)=ranksum(splus_licks,sminus_licks);
                            end
                            
                            
                            %suptitle(['PAC power LDA analysis for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC, mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).all_stamped_lick_ii=these_all_stamped_lick_ii;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).all_licks_per_tPACwave=these_all_licks_per_tPACwave;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_peak=p_val_peak;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_lick=p_val_lick;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).these_all_which_events_licks=these_all_which_events_licks;
                            
                            handles_out.discriminant_PACpower_per_mouse_all_mice.group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_peak=zeros(1,length(t));
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_peak(1,:)=discriminant_correct(1,:);
                            
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak=zeros(1,length(t));
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak(1,:)=dimensionality(1,:);
                            
                            
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_peak=zeros(1,length(t));
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_peak=auROC;
                            
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_peak=zeros(1,length(t));
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_peak(1,:)=discriminant_correct_shuffled(1,:);
                            
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).test_out_per_timepoint_peak=test_out_per_timepoint;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).shuffled_out_per_timepoint_peak=shuffled_out_per_timepoint;
                            
                            handles_out.t_power=t';
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).per_targets=per_targets;
                            
                            these_all_which_eventsPACwave=[];
                            these_all_which_eventsPACwave=all_which_eventsPACwave(:,these_per_corr);
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_eventsPACwave;
                            
                            %                                         these_all_stamped_lick_times=[];
                            %                                         these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                            %
                            %                                         these_all_stamped_lick_ii=[];
                            %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                            
                            %Linear discriminant analysis for
                            %trough PAC power
                            
                            
                            
                            fprintf(1, ['LDA processed for trough PAC wavelet power for all mice for ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],N);
                            fprintf(1,'For these events: ')
                            for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                            end
                            fprintf(1,'\n')
                            
                            
                            par_out=[];
                            test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                            shuffled_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                            discriminant_correct=zeros(1,length(t));
                            discriminant_correct_shuffled=zeros(1,length(t));
                            auROC=zeros(1,length(t));
                            dimensionality=zeros(1,length(t));
                            per_targets=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                            
                            parfor time_point=1:length(t)
                                %                                 for time_point=1:length(t)
                                
                                %LFP power per trial per electrode
                                measurements=zeros(N,No_electrodes);
                                measurements(:,:)=these_all_log_P_timecoursePACwave_trough(:,:,time_point)';
                                %                         measurements=(measurements-repmat(mean(measurements),N,1))./repmat(std(measurements,0,1),N,1);
                                
                                %Dimensionality
                                %Rows: trials, Columns: electrodes
                                Signal=measurements;
                                dimensionality(time_point) = nansum(eig(cov(Signal)))^2/nansum(eig(cov(Signal)).^2);
                                
                                %Enter strings labeling each event (one event for
                                %each trial)
                                events=[];
                                
                                for ii=1:N
                                    this_event=zeros(length(handles.drgbchoices.events_to_discriminate),1);
                                    this_event=these_all_which_events(:,ii);
                                    
                                    events{ii,1}=handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(logical(this_event))};
                                    
                                end
                                
                                tested_events=[];
                                shuffled_tested_events=[];
                                scores=[];
                                for ii=1:N
                                    %Partition the data into training and test sets.
                                    
                                    %Create input and target vectors leaving one trial out
                                    %For per_input each column has the dF/F for one trial
                                    %each row is a single time point for dF/F for one of the cells
                                    %For per_target the top row is 1 if the odor is S+ and 0 if it is
                                    %S-, and row 2 has 1 for S-
                                    idxTrn=ones(N,1);
                                    idxTrn(ii)=0;
                                    idxTest=zeros(N,1);
                                    idxTest(ii)=1;
                                    
                                    %Store the training data in a table.
                                    tblTrn=[];
                                    tblTrn = array2table(measurements(logical(idxTrn),:));
                                    tblTrn.Y = events(logical(idxTrn));
                                    
                                    %Train a discriminant analysis model using the training set and default options.
                                    %By default this is a regularized linear discriminant analysis (LDA)
                                    Mdl = fitcdiscr(tblTrn,'Y');
                                    
                                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                                    [label,score] = predict(Mdl,measurements(logical(idxTest),:));
                                    
                                    tested_events{ii,1}=label{1};
                                    scores(ii)=score(2);
                                    
                                    %Do LDA with shuffled trials
                                    shuffled_measurements=zeros(N,No_electrodes);
                                    shuffled_measurements(:,:)=measurements(randperm(N),:);
                                    
                                    %Store the training data in a table.
                                    sh_tblTrn=[];
                                    sh_tblTrn = array2table(shuffled_measurements(logical(idxTrn),:));
                                    sh_tblTrn.Y = events(logical(idxTrn));
                                    
                                    %Train a discriminant analysis model using the training set and default options.
                                    %By default this is a regularized linear discriminant analysis (LDA)
                                    sh_Mdl = fitcdiscr(sh_tblTrn,'Y');
                                    
                                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                                    sh_label = predict(Mdl,shuffled_measurements(logical(idxTest),:));
                                    
                                    shuffled_tested_events{ii,1}=sh_label{1};
                                    
                                end
                                
                                %Calculate auROC
                                [X,Y,T,AUC] = perfcurve(events,scores',handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)});
                                auROC(1,time_point)=AUC-0.5;
                                %test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                                
                                
                                per_targets=these_all_which_events;
                                
                                test_out=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                                shuffled_out=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                                for ii=1:N
                                    for jj=1:length(handles.drgbchoices.events_to_discriminate)
                                        if strcmp(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(jj)},tested_events{ii})
                                            test_out(jj,ii)=1;
                                        else
                                            test_out(jj,ii)=0;
                                        end
                                        if strcmp(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(jj)},shuffled_tested_events{ii})
                                            shuffled_out(jj,ii)=1;
                                        else
                                            shuffled_out(jj,ii)=0;
                                        end
                                    end
                                end
                                
                                test_out_per_timepoint(:,:,time_point)=test_out;
                                shuffled_out_per_timepoint(:,:,time_point)=shuffled_out;
                                discriminant_correct(1,time_point)=100*sum(sum(test_out.*per_targets))/N;
                                discriminant_correct_shuffled(1,time_point)=100*sum(sum(shuffled_out.*per_targets))/N;
                                fprintf(1, 'LDA for trough PAC wavelet power percent correct classification %d (for timepoint %d out of %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,length(t));
                                
                            end
                            
                            p_val_trough=zeros(1,length(t));
                            for ii_t=1:length(t)
                                splus_out=zeros(1,sum(these_all_which_events(1,:)==1));
                                splus_out(1,:)=test_out_per_timepoint(1,these_all_which_events(1,:)==1,ii_t);
                                sminus_out=zeros(1,sum(these_all_which_events(1,:)==0));
                                sminus_out(1,:)=test_out_per_timepoint(1,these_all_which_events(1,:)==0,ii_t);
                                p_val_trough(ii_t)=ranksum(splus_out,sminus_out);
                            end
                            
                            
                            subplot(2,5,6)
                            hold on
                            
                            per95=prctile(discriminant_correct_shuffled(1,:),95);
                            per5=prctile(discriminant_correct_shuffled(1,:),5);
                            CIsh=[mean(discriminant_correct_shuffled(1,:))-per5 per95-mean(discriminant_correct_shuffled(1,:))]';
                            [hlCR, hpCR] = boundedline([t(1) t(end)],[mean(discriminant_correct_shuffled(1,:)) mean(discriminant_correct_shuffled(1,:))], CIsh', 'r');
                            
                            plot(t',discriminant_correct(1,:),'-k')
                            
                            %Odor on markers
                            plot([0 0],[0 100],'-k')
                            odorhl=plot([0 2.5],[20 20],'-k','LineWidth',5);
                            plot([2.5 2.5],[0 100],'-k')
                            
                            %title(['LDA % correct for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str_all_mice ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                            xlabel('Time (sec)')
                            ylabel('% correct trough')
                            
                            subplot(2,5,7)
                            hold on
                            
                            plot(t,auROC,'-b')
                            
                            %Odor on markers
                            plot([0 0],[0 0.5],'-k')
                            odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                            plot([2.5 2.5],[0 0.5],'-k')
                            
                            %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str_all_mice ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                            xlabel('Time (sec)')
                            ylabel('auROC trough')
                            ylim([-0.3 0.6])
                            
                            %Plot dimensionality
                            subplot(2,5,8)
                            hold on
                            
                            plot(t,dimensionality,'-b')
                            
                            mindim=min(dimensionality);
                            maxdim=max(dimensionality);
                            
                            %Odor on markers
                            plot([0 0],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                            odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                            plot([2.5 2.5],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                            
                            %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str_all_mice ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                            xlabel('Time (sec)')
                            ylabel('Dimensionality trough')
                            xlim([-2 5])
                            
                            
                            %Plot p value
                            subplot(2,5,[4,5,9,10])
                            hold on
                            
                            p1=plot(t,log10(p_val_lick),'-k');
                            p2=plot(t,log10(p_val_trough),'-b');
                            p3=plot(t,log10(p_val_peak),'-m');
                            plot([t(1) t(end)],[log10(0.05) log10(0.05)],'-r')
                            legend([p1 p2 p3],{'Licks','Trough','Peak'})
                            ylabel('log(p)')
                            xlabel('Time (sec)')
                            xlim([-2 5])
                            ylim([-10 0])
                            
                            suptitle(['PAC wavelet power LDA analysis for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                            
                            
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_trough=p_val_trough;
                            
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_trough=zeros(1,length(t));
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_trough(1,:)=discriminant_correct(1,:);
                            
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough=zeros(1,length(t));
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough(1,:)=dimensionality(1,:);
                            
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_trough=zeros(1,length(t));
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_trough=auROC;
                            
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_trough=zeros(1,length(t));
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_trough(1,:)=discriminant_correct_shuffled(1,:);
                            
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).test_out_per_timepoint_trough=test_out_per_timepoint;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).shuffled_out_per_timepoint_trough=shuffled_out_per_timepoint;
                            
                            
                            %                                         handles_out.t_power=t';
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).per_targets=per_targets;
                            %
                            %                                         these_all_which_events=[];
                            %                                         these_all_which_events=all_which_eventsPACwave(:,these_per_corr);
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                            %
                            %                                         these_all_stamped_lick_times=[];
                            %                                         these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                            %
                            %                                         these_all_stamped_lick_ii=[];
                            %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                            
                            
                        end
                        
                        if sum(handles.drgbchoices.which_discriminant==16)>0
                            
                            %PCA for peak PAC wavelet power for all mice
                            
                            
                            
                            fprintf(1, ['PCA processed for peak PAC wavelet power for mall mice for ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],N);
                            fprintf(1,'For these events: ')
                            for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                            end
                            fprintf(1,'\n')
                            
                            for time_point=1:length(t)
                                par_t_out(time_point).principal_components=zeros(N,No_electrodes);
                            end
                            
                            
                            for time_point=1:length(t)
                                
                                %LFP power per trial per electrode
                                measurements=zeros(N,No_electrodes);
                                measurements(:,:)=these_all_log_P_timecoursePACwave_peak(:,:,time_point)';
                                
                                %Do the PCA
                                [coeff,par_t_out(time_point).principal_components,par_t_out(time_point).PC_variance]=pca(measurements);
                                
                            end
                            
                            %NOTE: Useful info from MATLAB answers. MATLAB PCA normalizes the input raw
                            %data so that the normalized data has zero mean (does not scale it for standard deviation).
                            %Because of this the following code holds true
                            % mydata = 10 + randn(20,5); %Random data 20 observations, 5 variables
                            % [coeff,scores_a] = pca(mydata); %Do PCA
                            % mydata_mean = mean(mydata); %Find mean of data (columns)
                            % mydata_mean = repmat(mydata_mean,20,1); %Replicate mean vector to matrix for subtraction
                            % my_data_norm = mydata - mydata_mean; % Normalize data to zero mean y subtraction
                            % scores_b = my_data_norm*coeff; %Manually calculate scores using PCA coeff and normalized data
                            % err = max(max((abs(scores_a - scores_b)))) %Calculate error as the max of absolute difference in 2 methods
                            % For my random data runs, err was of the order of 1e-15
                            
                            %Show a figure of the PCA and record the output
                            
                            principal_components=zeros(length(t),N,size(par_t_out(1).principal_components,2));
                            PC_variance=zeros(length(t),size(par_t_out(1).principal_components,2));
                            for time_point=1:length(t)
                                principal_components(time_point,:,:)=par_t_out(time_point).principal_components;
                                PC_variance(time_point,:)=par_t_out(time_point).PC_variance;
                            end
                            
                            %Show the result of the PCA
                            
                            figNo=figNo+1;
                            try
                                close(figNo)
                            catch
                            end
                            
                            figure(figNo)
                            
                            %Show PCA before odor on
                            subplot(2,4,5)
                            hold on
                            
                            these_pcs=zeros(N,size(par_t_out(1).principal_components,2));
                            these_pcs(:,:)=principal_components(6,:,:);
                            plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'or')
                            plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'ob')
                            legend(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(1)},handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)})
                            xlabel('PC1')
                            ylabel('PC2')
                            title('-1 sec peak')
                            
                            %Show PCA after odor
                            subplot(2,4,6)
                            hold on
                            
                            these_pcs=zeros(N,size(par_t_out(1).principal_components,2));
                            these_pcs(:,:)=principal_components(41,:,:);
                            plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'or')
                            plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'ob')
                            legend(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(1)},handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)})
                            xlabel('PC1')
                            ylabel('PC2')
                            title('2.5 sec peak')
                            
                            %Show the timecourse for PC1
                            subplot(2,4,[1,2])
                            hold on
                            
                            PC1ev1=zeros(length(t),sum(these_all_which_events(1,:)));
                            PC1ev1(:,:)=principal_components(:,logical(these_all_which_events(1,:)),1);
                            
                            PC1ev2=zeros(length(t),sum(these_all_which_events(2,:)));
                            PC1ev2(:,:)=principal_components(:,logical(these_all_which_events(2,:)),1);
                            
                            mean_PC1ev2=mean(PC1ev2,2)';
                            CIPC1ev2 = bootci(1000, {@mean, PC1ev2'});
                            maxCIPC1ev2=max(CIPC1ev2(:));
                            minCIPC1ev2=min(CIPC1ev2(:));
                            CIPC1ev2(1,:)=mean_PC1ev2-CIPC1ev2(1,:);
                            CIPC1ev2(2,:)=CIPC1ev2(2,:)-mean_PC1ev2;
                            [hlCR, hpCR] = boundedline(t',mean_PC1ev2', CIPC1ev2', 'b');
                            
                            mean_PC1ev1=mean(PC1ev1,2)';
                            CIPC1ev1 = bootci(1000, {@mean, PC1ev1'});
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
                            text(-1,minPC1+0.9*(maxPC1-minPC1),handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(1)},'Color','r')
                            text(-1,minPC1+0.8*(maxPC1-minPC1),handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)},'Color','b')
                            title('PC1 for peak')
                            xlabel('Time (sec)')
                            ylabel('PC1')
                            
                            
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PCA_calculated=1;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_peak=zeros(length(t),N,size(par_t_out(1).principal_components,2));
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_peak=principal_components;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).PC_variance_peak=zeros(length(t),size(par_t_out(1).principal_components,2));
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).PC_variance_peak=PC_variance;
                            
                            
                            handles_out.t_power=t_power';
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                            
                            these_all_which_eventsPACwave=[];
                            these_all_which_eventsPACwave=all_which_eventsPACwave(:,these_per_corr);
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_eventsPACwave;
                            
                            %                                         these_all_stamped_lick_times=[];
                            %                                         these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                            %
                            %                                         these_all_stamped_lick_ii=[];
                            %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                            
                            
                            %PCA for trough PAC wavelet power
                            
                            
                            
                            fprintf(1, ['PCA processed for trough PAC wavelet power trough for all mice for ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],N);
                            fprintf(1,'For these events: ')
                            for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                            end
                            fprintf(1,'\n')
                            
                            for time_point=1:length(t)
                                par_t_out(time_point).principal_components=zeros(N,size(par_t_out(1).principal_components,2));
                            end
                            
                            for time_point=1:length(t)
                                par_t_out(time_point).principal_components=[];
                            end
                            
                            for time_point=1:length(t)
                                
                                %LFP power per trial per electrode
                                measurements=zeros(N,No_electrodes);
                                measurements(:,:)=these_all_log_P_timecoursePACwave_trough(:,:,time_point)';
                                
                                %Do the PCA
                                [coeff,par_t_out(time_point).principal_components,par_t_out(time_point).PC_variance]=pca(measurements);
                                
                            end
                            
                            %NOTE: Useful info from MATLAB answers. MATLAB PCA normalizes the input raw
                            %data so that the normalized data has zero mean (does not scale it for standard deviation).
                            %Because of this the following code holds true
                            % mydata = 10 + randn(20,5); %Random data 20 observations, 5 variables
                            % [coeff,scores_a] = pca(mydata); %Do PCA
                            % mydata_mean = mean(mydata); %Find mean of data (columns)
                            % mydata_mean = repmat(mydata_mean,20,1); %Replicate mean vector to matrix for subtraction
                            % my_data_norm = mydata - mydata_mean; % Normalize data to zero mean y subtraction
                            % scores_b = my_data_norm*coeff; %Manually calculate scores using PCA coeff and normalized data
                            % err = max(max((abs(scores_a - scores_b)))) %Calculate error as the max of absolute difference in 2 methods
                            % For my random data runs, err was of the order of 1e-15
                            
                            %Show a figure of the PCA and record the output
                            
                            principal_components=zeros(length(t),N,size(par_t_out(1).principal_components,2));
                            PC_variance=zeros(length(t),size(par_t_out(1).principal_components,2));
                            for time_point=1:length(t)
                                principal_components(time_point,:,:)=par_t_out(time_point).principal_components;
                                PC_variance(time_point,:)=par_t_out(time_point).PC_variance;
                            end
                            
                            
                            
                            %Show PCA before odor on
                            subplot(2,4,7)
                            hold on
                            
                            these_pcs=zeros(N,size(par_t_out(1).principal_components,2));
                            these_pcs(:,:)=principal_components(6,:,:);
                            plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'or')
                            plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'ob')
                            legend(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(1)},handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)})
                            xlabel('PC1')
                            ylabel('PC2')
                            title('-1 sec trough')
                            
                            %Show PCA after odor
                            subplot(2,4,8)
                            hold on
                            
                            these_pcs=zeros(N,size(par_t_out(1).principal_components,2));
                            these_pcs(:,:)=principal_components(41,:,:);
                            plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'or')
                            plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'ob')
                            legend(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(1)},handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)})
                            xlabel('PC1')
                            ylabel('PC2')
                            title('2.5 sec trough')
                            
                            %Show the timecourse for PC1
                            subplot(2,4,[3,4])
                            hold on
                            
                            PC1ev1=zeros(length(t),sum(these_all_which_events(1,:)));
                            PC1ev1(:,:)=principal_components(:,logical(these_all_which_events(1,:)),1);
                            
                            PC1ev2=zeros(length(t),sum(these_all_which_events(2,:)));
                            PC1ev2(:,:)=principal_components(:,logical(these_all_which_events(2,:)),1);
                            
                            mean_PC1ev2=mean(PC1ev2,2)';
                            CIPC1ev2 = bootci(1000, {@mean, PC1ev2'});
                            maxCIPC1ev2=max(CIPC1ev2(:));
                            minCIPC1ev2=min(CIPC1ev2(:));
                            CIPC1ev2(1,:)=mean_PC1ev2-CIPC1ev2(1,:);
                            CIPC1ev2(2,:)=CIPC1ev2(2,:)-mean_PC1ev2;
                            [hlCR, hpCR] = boundedline(t',mean_PC1ev2', CIPC1ev2', 'b');
                            
                            mean_PC1ev1=mean(PC1ev1,2)';
                            CIPC1ev1 = bootci(1000, {@mean, PC1ev1'});
                            maxCIPC1ev1=max(CIPC1ev1(:));
                            minCIPC1ev1=min(CIPC1ev1(:));
                            CIPC1ev1(1,:)=mean_PC1ev1-CIPC1ev1(1,:);
                            CIPC1ev1(2,:)=CIPC1ev1(2,:)-mean_PC1ev1;
                            [hlCR, hpCR] = boundedline(t',mean_PC1ev1', CIPC1ev1', 'r');
                            
                            %                                         maxPC1=max([maxCIPC1ev2 maxCIPC1ev1]);
                            %                                         minPC1=min([minCIPC1ev2 minCIPC1ev1]);
                            
                            %Odor on markers
                            plot([0 0],[minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)],'-k')
                            odorhl=plot([0 2.5],[minPC1-0.1*(maxPC1-minPC1) minPC1-0.1*(maxPC1-minPC1)],'-k','LineWidth',5);
                            plot([2.5 2.5],[minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)],'-k')
                            
                            xlim([-2 5])
                            ylim([minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)])
                            text(-1,minPC1+0.9*(maxPC1-minPC1),handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(1)},'Color','r')
                            text(-1,minPC1+0.8*(maxPC1-minPC1),handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)},'Color','b')
                            suptitle(['PAC wavelet power PC1 for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC, for all mice  ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                            title('PC1 for trough')
                            xlabel('Time (sec)')
                            ylabel('PC1')
                            
                            
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_trough=zeros(length(t),N,size(par_t_out(1).principal_components,2));
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).principal_components_trough=principal_components;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).PC_variance_trough=zeros(length(t),size(par_t_out(1).principal_components,2));
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).PC_variance_trough=PC_variance;
                            
                            %                                         handles_out.t_power=t';
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                            %
                            %                                         these_all_which_events=[];
                            %                                         these_all_which_events=all_which_eventsPACwave(:,these_per_corr);
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                            %
                            %                                         these_all_stamped_lick_times=[];
                            %                                         these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                            %
                            %                                         these_all_stamped_lick_ii=[];
                            %                                         these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                            %                                         handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                            
                        end
                        
                        if sum(handles.drgbchoices.which_discriminant==17)>0
                            
                            %First do peak PAC wavelet power
                            %LDA for subsets of electrodes for PAC wavelet power
                            
                            t_from=2.15;
                            t_to=2.5;
                            subt=t((t>=t_from)&(t<=t_to));
                            
                            
                            
                            fprintf(1, ['LDA processed for PAC wavelet power for all mice ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],N);
                            fprintf(1,'For these events: ')
                            for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                            end
                            fprintf(1,'\n')
                            
                            
                            par_out=[];
                            max_combs=25;
                            for no_elect=1:No_electrodes
                                par_out(no_elect).no_elect=0;
                                par_out(no_elect).no_timepoints=0;
                                par_out(no_elect).is_tetrode(1:max_combs+4)=0;
                                par_out(no_elect).discriminant_correct(1:250)=0;
                                par_out(no_elect).no_samples=0;
                                par_out(no_elect).discriminant_correct_shuffled(1:250)=0;
                                par_out(no_elect).is_tetrode_per_sample(1:250)=0;
                                par_out(no_elect).auROC(1:250)=0;
                                par_out(no_elect).no_elect_combs=0;
                            end
                            
                            %parfor no_elect=1:length(handles.drgbchoices.which_electrodes)
                            parfor no_elect=1:No_electrodes
                                
                                par_out(no_elect).no_elect=No_electrodes;
                                par_out(no_elect).no_timepoints=length(subt);
                                %Choose electrode combinations
                                %(max =25)
                                
                                elect_combs=nchoosek(1:No_electrodes,no_elect);
                                if size(elect_combs,1)>max_combs
                                    no_chosen=0;
                                    chosen_elecs=[];
                                    while no_chosen<max_combs
                                        add_these=randi(size(elect_combs,1),1,25-no_chosen);
                                        chosen_elecs=unique([chosen_elecs add_these]);
                                        no_chosen=length(chosen_elecs);
                                    end
                                    elect_combs=elect_combs(chosen_elecs,:);
                                end
                                
                                no_el_combs=size(elect_combs,1);
                                par_out(no_elect).no_elect_combs=no_el_combs;
                                par_out(no_elect).is_tetrode(1:no_el_combs)=0;
                                
                                %if no_elect=4 also enter each
                                %tetrode
                                if no_elect==4
                                    for tetNo=1:length(handles.drgbchoices.which_electrodes)/4
                                        %Find if tetrode is included
                                        tetrode_found=0;
                                        for ii=1:no_el_combs
                                            if sum(elect_combs(ii,:)==[(tetNo-1)*4+1:(tetNo-1)*4+4])==no_elect
                                                tetrode_found=1;
                                                par_out(no_elect).is_tetrode(ii)=1;
                                            end
                                        end
                                        if tetrode_found==0
                                            %Add tetrode
                                            no_el_combs=no_el_combs+1;
                                            elect_combs(no_el_combs,:)=[(tetNo-1)*4+1:(tetNo-1)*4+4];
                                            par_out(no_elect).is_tetrode(no_el_combs)=1;
                                        end
                                        
                                    end
                                    
                                end
                                
                                par_out(no_elect).no_elect_combs=no_el_combs;
                                
                                no_samples=0;
                                for noelc=1:no_el_combs
                                    for time_point=1:length(subt)
                                        
                                        %LFP power per trial per electrode
                                        measurements=zeros(N,length(elect_combs(noelc,:)));
                                        measurements(:,:)=these_all_log_P_timecoursePACwave(elect_combs(noelc,:),:,time_point)';
                                        
                                        %Enter strings labeling each event (one event for
                                        %each trial)
                                        events=[];
                                        
                                        for ii=1:N
                                            this_event=zeros(length(handles.drgbchoices.events_to_discriminate),1);
                                            this_event=these_all_which_events(:,ii);
                                            events{ii,1}=handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(logical(this_event))};
                                        end
                                        
                                        tested_events=[];
                                        shuffled_tested_events=[];
                                        scores=[];
                                        for ii=1:N
                                            %Partition the data into training and test sets.
                                            
                                            %Create input and target vectors leaving one trial out
                                            %For per_input each column has the dF/F for one trial
                                            %each row is a single time point for dF/F for one of the cells
                                            %For per_target the top row is 1 if the odor is S+ and 0 if it is
                                            %S-, and row 2 has 1 for S-
                                            idxTrn=ones(N,1);
                                            idxTrn(ii)=0;
                                            idxTest=zeros(N,1);
                                            idxTest(ii)=1;
                                            
                                            %Store the training data in a table.
                                            tblTrn=[];
                                            tblTrn = array2table(measurements(logical(idxTrn),:));
                                            tblTrn.Y = events(logical(idxTrn));
                                            
                                            %Train a discriminant analysis model using the training set and default options.
                                            %By default this is a regularized linear discriminant analysis (LDA)
                                            Mdl = fitcdiscr(tblTrn,'Y');
                                            
                                            %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                                            [label,score] = predict(Mdl,measurements(logical(idxTest),:));
                                            
                                            tested_events{ii,1}=label{1};
                                            scores(ii)=score(2);
                                            
                                            %Do LDA with shuffled trials
                                            shuffled_measurements=zeros(N,size(measurements,2));
                                            shuffled_measurements(:,:)=measurements(randperm(N),:);
                                            
                                            %Store the training data in a table.
                                            sh_tblTrn=[];
                                            sh_tblTrn = array2table(shuffled_measurements(logical(idxTrn),:));
                                            sh_tblTrn.Y = events(logical(idxTrn));
                                            
                                            %Train a discriminant analysis model using the training set and default options.
                                            %By default this is a regularized linear discriminant analysis (LDA)
                                            sh_Mdl = fitcdiscr(sh_tblTrn,'Y');
                                            
                                            %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                                            sh_label = predict(Mdl,shuffled_measurements(logical(idxTest),:));
                                            
                                            shuffled_tested_events{ii,1}=sh_label{1};
                                            
                                        end
                                        
                                        %Calculate auROC
                                        [X,Y,T,AUC] = perfcurve(events,scores',handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)});
                                        %                                                     auROC(1,time_point)=AUC-0.5;
                                        %test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                                        
                                        
                                        per_targets=these_all_which_events;
                                        
                                        test_out=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                                        shuffled_out=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                                        for ii=1:N
                                            for jj=1:length(handles.drgbchoices.events_to_discriminate)
                                                if strcmp(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(jj)},tested_events{ii})
                                                    test_out(jj,ii)=1;
                                                else
                                                    test_out(jj,ii)=0;
                                                end
                                                if strcmp(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(jj)},shuffled_tested_events{ii})
                                                    shuffled_out(jj,ii)=1;
                                                else
                                                    shuffled_out(jj,ii)=0;
                                                end
                                            end
                                        end
                                        
                                        %                                                     test_out_per_timepoint(:,:,time_point)=test_out;
                                        %                                                     shuffled_out_per_timepoint(:,:,time_point)=shuffled_out;
                                        no_samples=no_samples+1;
                                        par_out(no_elect).no_samples=no_samples;
                                        if par_out(no_elect).is_tetrode(noelc)==1
                                            par_out(no_elect).is_tetrode_per_sample(no_samples)=1;
                                        else
                                            par_out(no_elect).is_tetrode_per_sample(no_samples)=0;
                                        end
                                        par_out(no_elect).discriminant_correct(no_samples)=100*sum(sum(test_out.*per_targets))/N;
                                        par_out(no_elect).discriminant_correct_shuffled(no_samples)=100*sum(sum(shuffled_out.*per_targets))/N;
                                        par_out(no_elect).auROC(no_samples)=AUC-0.5;
                                        fprintf(1, 'LDA PAC wavelet power percent correct classification %d (for timepoint %d and number of electrodes %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,no_elect);
                                        
                                    end
                                end
                            end
                            
                            if PACii==3
                                figNo=figNo+1;
                                try
                                    close(figNo)
                                catch
                                end
                                
                                figure(figNo)
                                
                                subplot(1,2,1)
                                hold on
                                
                                for elNo=1:par_out(1).no_elect
                                    
                                    mean_dcsh(elNo)=mean(par_out(elNo).discriminant_correct_shuffled(1:par_out(elNo).no_samples));
                                    tempCIdcsh = bootci(1000, {@mean, par_out(elNo).discriminant_correct_shuffled(1:par_out(elNo).no_samples)})';
                                    CIdcsh(elNo,1)=mean_dcsh(elNo)-tempCIdcsh(1);
                                    CIdcsh(elNo,2)=tempCIdcsh(2)-mean_dcsh(elNo);
                                    
                                    
                                    mean_dc(elNo)=mean(par_out(elNo).discriminant_correct(1:par_out(elNo).no_samples));
                                    tempCIdc = bootci(1000, {@mean, par_out(elNo).discriminant_correct(1:par_out(elNo).no_samples)})';
                                    CIdc(elNo,1)=mean_dc(elNo)-tempCIdc(1);
                                    CIdc(elNo,2)=tempCIdc(2)-mean_dc(elNo);
                                    
                                end
                                
                                
                                [hlCR, hpCR] = boundedline(1:par_out(1).no_elect,mean_dcsh, CIdcsh, 'b');
                                [hlCR, hpCR] = boundedline(1:par_out(1).no_elect,mean_dc, CIdc, 'r');
                                
                                %Plot tetrodes here
                                elNo=4;
                                tetNo=0;
                                for sampNo=1:par_out(elNo).no_timepoints:par_out(elNo).no_samples
                                    if par_out(elNo).is_tetrode_per_sample(sampNo)==1
                                        tetNo=tetNo+1;
                                        mean_pcorr(tetNo)=mean(par_out(elNo).discriminant_correct(sampNo:sampNo+par_out(elNo).no_timepoints-1));
                                    end
                                end
                                
                                plot(elNo*ones(1,tetNo),mean_pcorr,'or')
                                
                                xlim([1 par_out(1).no_elect])
                                ylim([40 110])
                                
                                
                                xlabel('Number of electrodes used')
                                ylabel('Percent correct')
                                
                                subplot(1,2,2)
                                hold on
                                
                                
                                for elNo=1:par_out(1).no_elect
                                    mean_auROC(elNo)=mean(par_out(elNo).auROC(1:par_out(elNo).no_samples));
                                    tempCIauROC = bootci(1000, {@mean, par_out(elNo).auROC(1:par_out(elNo).no_samples)})';
                                    CIauROC(elNo,1)=mean_auROC(elNo)-tempCIauROC(1);
                                    CIauROC(elNo,2)=tempCIauROC(2)-mean_auROC(elNo);
                                end
                                
                                [hlCR, hpCR] = boundedline(1:par_out(1).no_elect,mean_auROC, CIauROC, 'r');
                                
                                xlim([1 par_out(1).no_elect])
                                ylim([0 0.5])
                                
                                xlabel('Number of electrodes used')
                                ylabel('auROC')
                                ylim([-0.3 0.6])
                                
                                suptitle(['PAC peak wavelet power LDA ' handles.drgbchoices.PACnames{PACii} ' mouse No '  ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                            end
                            
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).ecomb_discriminant_calculated_peak=1;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).ecomb_par_out_peak=par_out;
                            
                            %First do trough PAC wavelet power
                            %LDA for subsets of electrodes for PAC wavelet power
                            
                            t_from=2.15;
                            t_to=2.5;
                            subt=t((t>=t_from)&(t<=t_to));
                            
                            these_all_log_P_timecoursePACwave=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(subt));
                            these_all_log_P_timecoursePACwave(:,:,:)=all_log_P_timecoursePACwavetrough(PACii,:,these_per_corr,(t>=t_from)&(t<=t_to));
                            
                            these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                            for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                these_all_which_events(ii,:)= all_which_eventsPACwave(kk,these_per_corr);
                            end
                            
                            fprintf(1, ['LDA processed for PAC wavelet power for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
                            fprintf(1,'For these events: ')
                            for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                            end
                            fprintf(1,'\n')
                            
                            
                            par_out=[];
                            max_combs=25;
                            for no_elect=1:length(handles.drgbchoices.which_electrodes)
                                par_out(no_elect).no_elect=0;
                                par_out(no_elect).no_timepoints=0;
                                par_out(no_elect).is_tetrode(1:max_combs+4)=0;
                                par_out(no_elect).discriminant_correct(1:250)=0;
                                par_out(no_elect).no_samples=0;
                                par_out(no_elect).discriminant_correct_shuffled(1:250)=0;
                                par_out(no_elect).is_tetrode_per_sample(1:250)=0;
                                par_out(no_elect).auROC(1:250)=0;
                                par_out(no_elect).no_elect_combs=0;
                            end
                            
                            %parfor no_elect=1:length(handles.drgbchoices.which_electrodes)
                            parfor no_elect=1:length(handles.drgbchoices.which_electrodes)
                                
                                par_out(no_elect).no_elect=length(handles.drgbchoices.which_electrodes);
                                par_out(no_elect).no_timepoints=length(subt);
                                %Choose electrode combinations
                                %(max =25)
                                
                                elect_combs=nchoosek(1:length(handles.drgbchoices.which_electrodes),no_elect);
                                if size(elect_combs,1)>max_combs
                                    no_chosen=0;
                                    chosen_elecs=[];
                                    while no_chosen<max_combs
                                        add_these=randi(size(elect_combs,1),1,25-no_chosen);
                                        chosen_elecs=unique([chosen_elecs add_these]);
                                        no_chosen=length(chosen_elecs);
                                    end
                                    elect_combs=elect_combs(chosen_elecs,:);
                                end
                                
                                no_el_combs=size(elect_combs,1);
                                par_out(no_elect).no_elect_combs=no_el_combs;
                                par_out(no_elect).is_tetrode(1:no_el_combs)=0;
                                
                                %if no_elect=4 also enter each
                                %tetrode
                                if no_elect==4
                                    for tetNo=1:length(handles.drgbchoices.which_electrodes)/4
                                        %Find if tetrode is included
                                        tetrode_found=0;
                                        for ii=1:no_el_combs
                                            if sum(elect_combs(ii,:)==[(tetNo-1)*4+1:(tetNo-1)*4+4])==no_elect
                                                tetrode_found=1;
                                                par_out(no_elect).is_tetrode(ii)=1;
                                            end
                                        end
                                        if tetrode_found==0
                                            %Add tetrode
                                            no_el_combs=no_el_combs+1;
                                            elect_combs(no_el_combs,:)=[(tetNo-1)*4+1:(tetNo-1)*4+4];
                                            par_out(no_elect).is_tetrode(no_el_combs)=1;
                                        end
                                        
                                    end
                                    
                                end
                                
                                par_out(no_elect).no_elect_combs=no_el_combs;
                                
                                no_samples=0;
                                for noelc=1:no_el_combs
                                    for time_point=1:length(subt)
                                        
                                        %LFP power per trial per electrode
                                        measurements=zeros(N,length(elect_combs(noelc,:)));
                                        measurements(:,:)=these_all_log_P_timecoursePACwave(elect_combs(noelc,:),:,time_point)';
                                        
                                        %Enter strings labeling each event (one event for
                                        %each trial)
                                        events=[];
                                        
                                        for ii=1:N
                                            this_event=zeros(length(handles.drgbchoices.events_to_discriminate),1);
                                            this_event=these_all_which_events(:,ii);
                                            events{ii,1}=handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(logical(this_event))};
                                        end
                                        
                                        tested_events=[];
                                        shuffled_tested_events=[];
                                        scores=[];
                                        for ii=1:N
                                            %Partition the data into training and test sets.
                                            
                                            %Create input and target vectors leaving one trial out
                                            %For per_input each column has the dF/F for one trial
                                            %each row is a single time point for dF/F for one of the cells
                                            %For per_target the top row is 1 if the odor is S+ and 0 if it is
                                            %S-, and row 2 has 1 for S-
                                            idxTrn=ones(N,1);
                                            idxTrn(ii)=0;
                                            idxTest=zeros(N,1);
                                            idxTest(ii)=1;
                                            
                                            %Store the training data in a table.
                                            tblTrn=[];
                                            tblTrn = array2table(measurements(logical(idxTrn),:));
                                            tblTrn.Y = events(logical(idxTrn));
                                            
                                            %Train a discriminant analysis model using the training set and default options.
                                            %By default this is a regularized linear discriminant analysis (LDA)
                                            Mdl = fitcdiscr(tblTrn,'Y');
                                            
                                            %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                                            [label,score] = predict(Mdl,measurements(logical(idxTest),:));
                                            
                                            tested_events{ii,1}=label{1};
                                            scores(ii)=score(2);
                                            
                                            %Do LDA with shuffled trials
                                            shuffled_measurements=zeros(N,size(measurements,2));
                                            shuffled_measurements(:,:)=measurements(randperm(N),:);
                                            
                                            %Store the training data in a table.
                                            sh_tblTrn=[];
                                            sh_tblTrn = array2table(shuffled_measurements(logical(idxTrn),:));
                                            sh_tblTrn.Y = events(logical(idxTrn));
                                            
                                            %Train a discriminant analysis model using the training set and default options.
                                            %By default this is a regularized linear discriminant analysis (LDA)
                                            sh_Mdl = fitcdiscr(sh_tblTrn,'Y');
                                            
                                            %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                                            sh_label = predict(Mdl,shuffled_measurements(logical(idxTest),:));
                                            
                                            shuffled_tested_events{ii,1}=sh_label{1};
                                            
                                        end
                                        
                                        %Calculate auROC
                                        [X,Y,T,AUC] = perfcurve(events,scores',handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)});
                                        %                                                     auROC(1,time_point)=AUC-0.5;
                                        %test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                                        
                                        
                                        per_targets=these_all_which_events;
                                        
                                        test_out=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                                        shuffled_out=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                                        for ii=1:N
                                            for jj=1:length(handles.drgbchoices.events_to_discriminate)
                                                if strcmp(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(jj)},tested_events{ii})
                                                    test_out(jj,ii)=1;
                                                else
                                                    test_out(jj,ii)=0;
                                                end
                                                if strcmp(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(jj)},shuffled_tested_events{ii})
                                                    shuffled_out(jj,ii)=1;
                                                else
                                                    shuffled_out(jj,ii)=0;
                                                end
                                            end
                                        end
                                        
                                        %                                                     test_out_per_timepoint(:,:,time_point)=test_out;
                                        %                                                     shuffled_out_per_timepoint(:,:,time_point)=shuffled_out;
                                        no_samples=no_samples+1;
                                        par_out(no_elect).no_samples=no_samples;
                                        if par_out(no_elect).is_tetrode(noelc)==1
                                            par_out(no_elect).is_tetrode_per_sample(no_samples)=1;
                                        else
                                            par_out(no_elect).is_tetrode_per_sample(no_samples)=0;
                                        end
                                        par_out(no_elect).discriminant_correct(no_samples)=100*sum(sum(test_out.*per_targets))/N;
                                        par_out(no_elect).discriminant_correct_shuffled(no_samples)=100*sum(sum(shuffled_out.*per_targets))/N;
                                        par_out(no_elect).auROC(no_samples)=AUC-0.5;
                                        fprintf(1, 'LDA PAC wavelet power percent correct classification %d (for timepoint %d and number of electrodes %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,no_elect);
                                        
                                    end
                                end
                            end
                            
                            if PACii==3
                                figNo=figNo+1;
                                try
                                    close(figNo)
                                catch
                                end
                                
                                figure(figNo)
                                
                                subplot(1,2,1)
                                hold on
                                
                                for elNo=1:par_out(1).no_elect
                                    
                                    mean_dcsh(elNo)=mean(par_out(elNo).discriminant_correct_shuffled(1:par_out(elNo).no_samples));
                                    tempCIdcsh = bootci(1000, {@mean, par_out(elNo).discriminant_correct_shuffled(1:par_out(elNo).no_samples)})';
                                    CIdcsh(elNo,1)=mean_dcsh(elNo)-tempCIdcsh(1);
                                    CIdcsh(elNo,2)=tempCIdcsh(2)-mean_dcsh(elNo);
                                    
                                    
                                    mean_dc(elNo)=mean(par_out(elNo).discriminant_correct(1:par_out(elNo).no_samples));
                                    tempCIdc = bootci(1000, {@mean, par_out(elNo).discriminant_correct(1:par_out(elNo).no_samples)})';
                                    CIdc(elNo,1)=mean_dc(elNo)-tempCIdc(1);
                                    CIdc(elNo,2)=tempCIdc(2)-mean_dc(elNo);
                                    
                                end
                                
                                
                                [hlCR, hpCR] = boundedline(1:par_out(1).no_elect,mean_dcsh, CIdcsh, 'b');
                                [hlCR, hpCR] = boundedline(1:par_out(1).no_elect,mean_dc, CIdc, 'r');
                                
                                %Plot tetrodes here
                                elNo=4;
                                tetNo=0;
                                for sampNo=1:par_out(elNo).no_timepoints:par_out(elNo).no_samples
                                    if par_out(elNo).is_tetrode_per_sample(sampNo)==1
                                        tetNo=tetNo+1;
                                        mean_pcorr(tetNo)=mean(par_out(elNo).discriminant_correct(sampNo:sampNo+par_out(elNo).no_timepoints-1));
                                    end
                                end
                                
                                plot(elNo*ones(1,tetNo),mean_pcorr,'or')
                                
                                xlim([1 par_out(1).no_elect])
                                ylim([40 110])
                                
                                
                                xlabel('Number of electrodes used')
                                ylabel('Percent correct')
                                
                                subplot(1,2,2)
                                hold on
                                
                                
                                for elNo=1:par_out(1).no_elect
                                    mean_auROC(elNo)=mean(par_out(elNo).auROC(1:par_out(elNo).no_samples));
                                    tempCIauROC = bootci(1000, {@mean, par_out(elNo).auROC(1:par_out(elNo).no_samples)})';
                                    CIauROC(elNo,1)=mean_auROC(elNo)-tempCIauROC(1);
                                    CIauROC(elNo,2)=tempCIauROC(2)-mean_auROC(elNo);
                                end
                                
                                [hlCR, hpCR] = boundedline(1:par_out(1).no_elect,mean_auROC, CIauROC, 'r');
                                
                                xlim([1 par_out(1).no_elect])
                                ylim([0 0.5])
                                
                                xlabel('Number of electrodes used')
                                ylabel('auROC')
                                ylim([-0.3 0.6])
                                
                                suptitle(['PAC wavelet power LDA ' handles.drgbchoices.PACnames{PACii} ' mouse No ' num2str_all_mice ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                            end
                            
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).ecomb_discriminant_calculated_trough=1;
                            handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).ecomb_par_out_trough=par_out;
                        end
                        save(output_name,'handles_out','-v7.3')
                    else
                        handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=0;
                        fprintf(1, ['LDA/PCA not processed for PAC wavelet power for all mice ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' because there were only %d trials (fewer than 20 trials)\n'],N);
                    end
                    
                case 2
                    %First determine how many trials per mouse
                    Nmax=0;
                    Nmin=20000;
                    mouse_included=[];
                    max_mouse=[];
                    min_mouse=[];
                    these_all_licks_per_tPACwave=[];
                    these_all_stamped_lick_ii=[];
                    these_all_which_events_licks=[];
                    N_licks=0;
                    n_trials_per_mouse=[];
                    

                    for mouseNo=1:max(handles.drgbchoices.mouse_no)
                        if sum(mice_included==mouseNo)>0
                            mouse_included(mouseNo)=0;
                            try
                                if handles_out.all_mouse_wav(mouseNo).group(groupNo).mouse_has_data==1
                                    group_has_data=1;
                                else
                                    group_has_data=0;
                                end
                            catch
                                group_has_data=0;
                            end
                            if group_has_data==1
                                these_per_corr=(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave>=handles.drgbchoices.percent_windows(percent_correct_ii,1))...
                                    &(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave<=handles.drgbchoices.percent_windows(percent_correct_ii,2));
                                thisN=sum(these_per_corr);
                                this_mouseNo=mouseNo;
                                n_trials_per_mouse(mouseNo)=thisN;
                                
                                if sum(these_per_corr)>=min_trials
                                    mouse_included(mouseNo)=1;
                                    if Nmax<sum(these_per_corr)
                                        max_mouse=mouseNo;
                                        Nmax=sum(these_per_corr);
                                    end
                                    
                                    if Nmin>sum(these_per_corr)
                                        min_mouse=mouseNo;
                                        Nmin=sum(these_per_corr);
                                    end
                                else
                                    mouse_included(mouseNo)=0;
                                end
                                
                                
                                N_licks=N_licks+sum(these_per_corr);
                                
                            end
                        end
                    end
                    
                    reference_mouse=max_mouse;
                    mouseNo=max_mouse;
                    
                    
                    %Now process LDA for these_all_log_P_timecoursePACwave
                    %in batches of 30 trials each until you get to the end of the
                    %number of trials available for max_mouse
                    

                    
                    %Find all trials for this batch for max_mouse
                    

                    %Now find these_all_log_P_timecoursePACwave
                    all_these_all_log_P_timecoursePACwave_peak=zeros(length(handles.drgbchoices.which_electrodes),Nmax,length(t));
                    all_these_all_log_P_timecoursePACwave_trough=zeros(length(handles.drgbchoices.which_electrodes),Nmax,length(t));
                    
                    
                    no_mice_included=0;
                    handles_out.group(groupNo).percent_correct(percent_correct_ii).max_mouse=max_mouse;
                    fprintf(1, ['\n\nNumber of trials for each mouse for ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} '\n']);
                    
                    
                    %Compute these_variables for the max_mouse with the largest number of trials
                    
                    
                    no_mice_included=no_mice_included+1;
                    first_electrode=1;
                    
                    these_per_corr=(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave>=handles.drgbchoices.percent_windows(percent_correct_ii,1))...
                        &(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave<=handles.drgbchoices.percent_windows(percent_correct_ii,2));
                    
                    all_these_all_log_P_timecoursePACwave_peak(1:length(handles.drgbchoices.which_electrodes),:,:)=...
                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_log_P_timecoursePACwavepeak(PACii,:,these_per_corr,:);
                    
                    all_these_all_log_P_timecoursePACwave_trough(1:length(handles.drgbchoices.which_electrodes),:,:)=...
                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_log_P_timecoursePACwavetrough(PACii,:,these_per_corr,:);
                    
                    %Save these data and notify user of the number
                    %of trials
                    handles_out.all_mouse_wav(mouseNo).group(groupNo).percent_correct(percent_correct_ii).no_trials=sum(these_per_corr);
                    fprintf(1, ['Number of trials for mouse No %d = %d\n'],mouseNo,sum(these_per_corr));
                    
                    
                    all_these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                    for ii=1:length(handles.drgbchoices.events_to_discriminate)
                        kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                        all_these_all_which_events(ii,:)= handles_out.all_mouse_wav(mouseNo).group(groupNo).all_which_eventsPACwave(kk,these_per_corr);
                    end
                    
                    
                    %                         all_which_eventsPACwave=handles_out.all_mouse_wav(mouseNo).group(groupNo).all_which_eventsPACwave;
                    
                    
                    these_ii_splus=[];
                    ii_s=0;
                    these_ii_sminus=[];
                    ii_m=0;
                    for ii=1:length(all_these_all_which_events(1,:))
                        if all_these_all_which_events(1,ii)==1
                            ii_s=ii_s+1;
                            these_ii_splus(ii_s)=ii;
                        else
                            ii_m=ii_m+1;
                            these_ii_sminus(ii_m)=ii;
                        end
                    end
                    
                    %Do lick accounting
                    these_all_licks_per_tPACwave=zeros(N_licks,length(t));
                    these_all_stamped_lick_ii=zeros(1,N_licks);
                    these_all_which_events_licks=zeros(1,N_licks);
                    ii_licks=1;
                    
                    these_all_licks_per_tPACwave(ii_licks:ii_licks+ size(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_licks_per_tPACwave(these_per_corr,:),1)-1,:)=...
                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_licks_per_tPACwave(these_per_corr,:);
                    
                    these_all_stamped_lick_ii(1,ii_licks:ii_licks+ length(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,these_per_corr))-1)=...
                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,these_per_corr);
                    
                    these_all_which_events_licks(1,ii_licks:ii_licks+ length(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,these_per_corr))-1)=...
                        all_these_all_which_events(1,:);
                    
                    ii_licks=ii_licks+length(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,these_per_corr));
                    
                    N=30;
                    no_batches=min([floor(length(these_ii_splus)/(N/2)) floor(length(these_ii_sminus)/(N/2))])
                    
                    these_all_which_events=zeros(2,N);
                    these_all_which_events(1,1:N/2)=1;
                    these_all_which_events(2,N/2+1:end)=1;
                    
                    %Now find these_all_log_P_timecoursePACwave
                    these_all_log_P_timecoursePACwave_peak=zeros(no_batches,length(handles.drgbchoices.which_electrodes)*sum(mouse_included),N,length(t));
                    these_all_log_P_timecoursePACwave_trough=zeros(no_batches,length(handles.drgbchoices.which_electrodes)*sum(mouse_included),N,length(t));
                    
                    ii_splus=0;
                    ii_sminus=0;
                    first_electrode=1;
                    for ii_batch=1:no_batches
                        %Load the S+
                        %First transfer S+
                        for ii=1:N/2
                            ii_splus=ii_splus+1;
                            these_all_log_P_timecoursePACwave_peak(ii_batch,first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,ii,:)=...
                                all_these_all_log_P_timecoursePACwave_peak(:,these_ii_splus(ii_splus),:);
                            these_all_log_P_timecoursePACwave_trough(ii_batch,first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,ii,:)=...
                                all_these_all_log_P_timecoursePACwave_trough(:,these_ii_splus(ii_splus),:);
                        end
                        
                        %Then transfer S-
                        for ii=1:N/2
                            ii_sminus=ii_sminus+1;
                            these_all_log_P_timecoursePACwave_peak(ii_batch,first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,ii+(N/2),:)=...
                                all_these_all_log_P_timecoursePACwave_peak(:,these_ii_sminus(ii_sminus),:);
                            these_all_log_P_timecoursePACwave_trough(ii_batch,first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,ii+(N/2),:)=...
                                all_these_all_log_P_timecoursePACwave_trough(:,these_ii_sminus(ii_sminus),:);
                        end
                        
                    end
                    
                    first_electrode=first_electrode+length(handles.drgbchoices.which_electrodes);
                    
                    %Now enter the data for the rest of the mice
                    for mouseNo=1:max(handles.drgbchoices.mouse_no)
                        if sum(mice_included==mouseNo)>0
                            if reference_mouse~=mouseNo
                                if (mouse_included(mouseNo)==1)
                                    
                                    %Enter the existing trials
                                    no_mice_included=no_mice_included+1;
                                    
                                    those_per_corr=(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave>=handles.drgbchoices.percent_windows(percent_correct_ii,1))...
                                        &(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_perCorr_pertrPACwave<=handles.drgbchoices.percent_windows(percent_correct_ii,2));
                                    
                                    those_all_log_P_timecoursePACwave_peak=zeros(length(handles.drgbchoices.which_electrodes),sum(those_per_corr),length(t));
                                    those_all_log_P_timecoursePACwave_peak(:,:,:)=...
                                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_log_P_timecoursePACwavepeak(PACii,:,those_per_corr,:);
                                    
                                    those_all_log_P_timecoursePACwave_trough=zeros(length(handles.drgbchoices.which_electrodes),sum(those_per_corr),length(t));
                                    those_all_log_P_timecoursePACwave_trough(:,:,:)=...
                                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_log_P_timecoursePACwavetrough(PACii,:,those_per_corr,:);
                                    
                                    handles_out.all_mouse_wav(mouseNo).group(groupNo).percent_correct(percent_correct_ii).no_trials=sum(those_per_corr);
                                    fprintf(1, ['Number of trials for mouse No %d = %d\n'],mouseNo,sum(those_per_corr));
                                    
                                    those_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(those_per_corr));
                                    
                                    for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                        kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                        those_all_which_events(ii,:)= handles_out.all_mouse_wav(mouseNo).group(groupNo).all_which_eventsPACwave(kk,those_per_corr);
                                    end
                                    
                                    those_ii_splus=[];
                                    ii_s=0;
                                    those_ii_sminus=[];
                                    ii_m=0;
                                    for ii=1:length(those_all_which_events(1,:))
                                        if those_all_which_events(1,ii)==1
                                            ii_s=ii_s+1;
                                            those_ii_splus(ii_s)=ii;
                                        else
                                            ii_m=ii_m+1;
                                            those_ii_sminus(ii_m)=ii;
                                        end
                                    end
                                    
                                    %First transfer S+
                                    jj=1;
                                    for ii_batch=1:no_batches
                                        for ii=1:(N/2)
                                            these_all_log_P_timecoursePACwave_peak(ii_batch,first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,ii,:)=...
                                                those_all_log_P_timecoursePACwave_peak(:,those_ii_splus(jj),:);
                                            these_all_log_P_timecoursePACwave_trough(ii_batch,first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,ii,:)=...
                                                those_all_log_P_timecoursePACwave_trough(:,those_ii_splus(jj),:);
                                            jj=jj+1;
                                            if jj>length(those_ii_splus)
                                                jj=1;
                                            end
                                        end
                                    end
                                    
                                    %Then transfer S-
                                    jj=1;
                                    for ii_batch=1:no_batches
                                        for ii=1:(N/2)
                                            these_all_log_P_timecoursePACwave_peak(ii_batch,first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,ii+(N/2),:)=...
                                                those_all_log_P_timecoursePACwave_peak(:,those_ii_sminus(jj),:);
                                            these_all_log_P_timecoursePACwave_trough(ii_batch,first_electrode:first_electrode+length(handles.drgbchoices.which_electrodes)-1,ii+(N/2),:)=...
                                                those_all_log_P_timecoursePACwave_trough(:,those_ii_sminus(jj),:);
                                            jj=jj+1;
                                            if jj>length(those_ii_sminus)
                                                jj=1;
                                            end
                                        end
                                    end
                                    
                                    
                                    
                                    these_all_licks_per_tPACwave(ii_licks:ii_licks+ size(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_licks_per_tPACwave(those_per_corr,:),1)-1,:)=...
                                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_licks_per_tPACwave(those_per_corr,:);
                                    
                                    these_all_stamped_lick_ii(1,ii_licks:ii_licks+ length(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,those_per_corr))-1)=...
                                        handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,those_per_corr);
                                    
                                    these_all_which_events_licks(1,ii_licks:ii_licks+ length(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,those_per_corr))-1)=...
                                        those_all_which_events(1,:);
                                    
                                    ii_licks=ii_licks+length(handles_out.all_mouse_wav(mouseNo).group(groupNo).all_stamped_lick_iiPACwave(1,those_per_corr));
                                    
                                    first_electrode=first_electrode+length(handles.drgbchoices.which_electrodes);
                                end
                            end
                        end
                    end
                    
                    No_electrodes=first_electrode-1;
                    
                    %Calculate p value for licks and save
                    
                    %Calculate p value for the licks
                    p_val_lick=zeros(1,length(t));
                    for ii_t=1:length(t)
                        splus_licks=zeros(1,sum((these_all_which_events_licks==1)&(these_all_stamped_lick_ii>0)));
                        splus_licks(1,:)=these_all_licks_per_tPACwave((these_all_which_events_licks==1)&(these_all_stamped_lick_ii>0),ii_t);
                        sminus_licks=zeros(1,sum((these_all_which_events_licks==0)&(these_all_stamped_lick_ii>0)));
                        sminus_licks(1,:)=these_all_licks_per_tPACwave((these_all_which_events_licks==0)&(these_all_stamped_lick_ii>0),ii_t);
                        p_val_lick(ii_t)=ranksum(splus_licks,sminus_licks);
                    end
                    
                    handles_out.discriminant_PACwavepower_all_mice.no_mice_included=no_mice_included;
                    
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).all_stamped_lick_ii=these_all_stamped_lick_ii;
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).all_licks_per_tPACwave=these_all_licks_per_tPACwave;
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_lick=p_val_lick;
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).these_all_which_events_licks=these_all_which_events_licks;
                    
                    
                    
                    batch_discriminant_correct_peak=zeros(no_batches,length(t));
                    batch_dimensionality_peak=zeros(no_batches,length(t));
                    batch_auROC_peak=zeros(no_batches,length(t));
                    batch_discriminant_correct_shuffled_peak=zeros(1,length(t));
                    batch_test_out_per_timepoint_peak=zeros(no_batches,2,N,length(t));
                    batch_shuffled_out_per_timepoint_peak=zeros(no_batches,2,N,length(t));
                    
                    
                    
                    batch_discriminant_correct_trough=zeros(no_batches,length(t));
                    batch_dimensionality_trough=zeros(no_batches,length(t));
                    batch_auROC_trough=zeros(no_batches,length(t));
                    batch_discriminant_correct_shuffled_trough=zeros(1,length(t));
                    batch_test_out_per_timepoint_trough=zeros(no_batches,2,N,length(t));
                    batch_shuffled_out_per_timepoint_trough=zeros(no_batches,2,N,length(t));
                    
                    handles_out.t_power=t';
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).no_trials=N;
                    
                    
                    for ii_batch=1:no_batches
                        discriminant_correct=zeros(1,length(t));
                        discriminant_correct_shuffled=zeros(1,length(t));
                        auROC=zeros(1,length(t));
                        
                        
                        %Do the analysis only if there are more than 20 trials
                        
                        
                        
                        %Linear discriminant analysis for
                        %peak PAC power al mice
                        
                        
                        fprintf(1, ['\nLDA processed for peak PAC wavelet power for all mice for ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],N);
                        fprintf(1,'For these events: ')
                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                            fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                        end
                        fprintf(1,'\n')
                        
                        
                        par_out=[];
                        test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                        shuffled_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                        discriminant_correct=zeros(1,length(t));
                        discriminant_correct_shuffled=zeros(1,length(t));
                        auROC=zeros(1,length(t));
                        dimensionality=zeros(1,length(t));
                        per_targets=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                        per_targets=these_all_which_events;
                        
                        parfor time_point=1:length(t)
                            %                                                         for time_point=1:length(t)
                            
                            %LFP power per trial per electrode
                            measurements=zeros(No_electrodes,N);
                            measurements(:,:)=these_all_log_P_timecoursePACwave_peak(ii_batch,:,:,time_point);
                            measurements=measurements';
                            
                            %Dimensionality
                            %Rows: trials, Columns: electrodes
                            Signal=measurements;
                            dimensionality(time_point) = nansum(eig(cov(Signal)))^2/nansum(eig(cov(Signal)).^2);
                            
                            %Enter strings labeling each event (one event for
                            %each trial)
                            events=[];
                            
                            for ii=1:N
                                this_event=zeros(length(handles.drgbchoices.events_to_discriminate),1);
                                this_event=these_all_which_events(:,ii);
                                
                                events{ii,1}=handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(logical(this_event))};
                                
                            end
                            
                            tested_events=[];
                            shuffled_tested_events=[];
                            scores=[];
                            for ii=1:N
                                %Partition the data into training and test sets.
                                
                                %Create input and target vectors leaving one trial out
                                %For per_input each column has the dF/F for one trial
                                %each row is a single time point for dF/F for one of the cells
                                %For per_target the top row is 1 if the odor is S+ and 0 if it is
                                %S-, and row 2 has 1 for S-
                                idxTrn=ones(N,1);
                                idxTrn(ii)=0;
                                idxTest=zeros(N,1);
                                idxTest(ii)=1;
                                
                                %Store the training data in a table.
                                tblTrn=[];
                                tblTrn = array2table(measurements(logical(idxTrn),:));
                                tblTrn.Y = events(logical(idxTrn));
                                
                                %Train a discriminant analysis model using the training set and default options.
                                %By default this is a regularized linear discriminant analysis (LDA)
                                Mdl = fitcdiscr(tblTrn,'Y');
                                
                                %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                                [label,score] = predict(Mdl,measurements(logical(idxTest),:));
                                
                                tested_events{ii,1}=label{1};
                                scores(ii)=score(2);
                                
                                %Do LDA with shuffled trials
                                shuffled_measurements=zeros(N,No_electrodes);
                                shuffled_measurements(:,:)=measurements(randperm(N),:);
                                
                                %Store the training data in a table.
                                sh_tblTrn=[];
                                sh_tblTrn = array2table(shuffled_measurements(logical(idxTrn),:));
                                sh_tblTrn.Y = events(logical(idxTrn));
                                
                                %Train a discriminant analysis model using the training set and default options.
                                %By default this is a regularized linear discriminant analysis (LDA)
                                sh_Mdl = fitcdiscr(sh_tblTrn,'Y');
                                
                                %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                                sh_label = predict(Mdl,shuffled_measurements(logical(idxTest),:));
                                
                                shuffled_tested_events{ii,1}=sh_label{1};
                                
                            end
                            
                            %Calculate auROC
                            [X,Y,T,AUC] = perfcurve(events,scores',handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)});
                            auROC(1,time_point)=AUC-0.5;
                            %test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                            
                            test_out=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                            shuffled_out=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                            for ii=1:N
                                for jj=1:length(handles.drgbchoices.events_to_discriminate)
                                    if strcmp(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(jj)},tested_events{ii})
                                        test_out(jj,ii)=1;
                                    else
                                        test_out(jj,ii)=0;
                                    end
                                    if strcmp(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(jj)},shuffled_tested_events{ii})
                                        shuffled_out(jj,ii)=1;
                                    else
                                        shuffled_out(jj,ii)=0;
                                    end
                                end
                            end
                            
                            test_out_per_timepoint(:,:,time_point)=test_out;
                            shuffled_out_per_timepoint(:,:,time_point)=shuffled_out;
                            discriminant_correct(1,time_point)=100*sum(per_targets(1,:)==test_out(1,:))/N;
                            discriminant_correct_shuffled(1,time_point)=100*sum(per_targets(1,:)==shuffled_out(1,:))/N;
                            fprintf(1, 'LDA for peak PAC wavelet power percent correct classification %d (for timepoint %d out of %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,length(t));
                            
                        end
                        
                        batch_discriminant_correct_peak(ii_batch,:)=discriminant_correct(1,:);
                        batch_dimensionality_peak(ii_batch,:)=dimensionality(1,:);
                        batch_auROC_peak(ii_batch,:)=auROC;
                        batch_discriminant_correct_shuffled_peak(ii_batch,:)=discriminant_correct_shuffled(1,:);
                        batch_test_out_per_timepoint_peak(ii_batch,:,:,:)=test_out_per_timepoint;
                        batch_shuffled_out_per_timepoint_peak(ii_batch,:,:,:)=shuffled_out_per_timepoint;
                        
                        fprintf(1, ['LDA processed for trough PAC wavelet power for all mice for ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],N);
                        fprintf(1,'For these events: ')
                        for ii=1:length(handles.drgbchoices.events_to_discriminate)
                            fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                        end
                        fprintf(1,'\n')
                        
                        
                        par_out=[];
                        test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                        shuffled_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                        discriminant_correct=zeros(1,length(t));
                        discriminant_correct_shuffled=zeros(1,length(t));
                        auROC=zeros(1,length(t));
                        dimensionality=zeros(1,length(t));
                        per_targets=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                        
                        parfor time_point=1:length(t)
                            %                                 for time_point=1:length(t)
                            
                            %LFP power per trial per electrode
                            measurements=zeros(No_electrodes,N);
                            measurements(:,:)=these_all_log_P_timecoursePACwave_trough(ii_batch,:,:,time_point);
                            measurements=measurements';
                            
                            %Dimensionality
                            %Rows: trials, Columns: electrodes
                            Signal=measurements;
                            dimensionality(time_point) = nansum(eig(cov(Signal)))^2/nansum(eig(cov(Signal)).^2);
                            
                            %Enter strings labeling each event (one event for
                            %each trial)
                            events=[];
                            
                            for ii=1:N
                                this_event=zeros(length(handles.drgbchoices.events_to_discriminate),1);
                                this_event=these_all_which_events(:,ii);
                                
                                events{ii,1}=handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(logical(this_event))};
                                
                            end
                            
                            tested_events=[];
                            shuffled_tested_events=[];
                            scores=[];
                            for ii=1:N
                                %Partition the data into training and test sets.
                                
                                %Create input and target vectors leaving one trial out
                                %For per_input each column has the dF/F for one trial
                                %each row is a single time point for dF/F for one of the cells
                                %For per_target the top row is 1 if the odor is S+ and 0 if it is
                                %S-, and row 2 has 1 for S-
                                idxTrn=ones(N,1);
                                idxTrn(ii)=0;
                                idxTest=zeros(N,1);
                                idxTest(ii)=1;
                                
                                %Store the training data in a table.
                                tblTrn=[];
                                tblTrn = array2table(measurements(logical(idxTrn),:));
                                tblTrn.Y = events(logical(idxTrn));
                                
                                %Train a discriminant analysis model using the training set and default options.
                                %By default this is a regularized linear discriminant analysis (LDA)
                                Mdl = fitcdiscr(tblTrn,'Y');
                                
                                %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                                [label,score] = predict(Mdl,measurements(logical(idxTest),:));
                                
                                tested_events{ii,1}=label{1};
                                scores(ii)=score(2);
                                
                                %Do LDA with shuffled trials
                                shuffled_measurements=zeros(N,No_electrodes);
                                shuffled_measurements(:,:)=measurements(randperm(N),:);
                                
                                %Store the training data in a table.
                                sh_tblTrn=[];
                                sh_tblTrn = array2table(shuffled_measurements(logical(idxTrn),:));
                                sh_tblTrn.Y = events(logical(idxTrn));
                                
                                %Train a discriminant analysis model using the training set and default options.
                                %By default this is a regularized linear discriminant analysis (LDA)
                                sh_Mdl = fitcdiscr(sh_tblTrn,'Y');
                                
                                %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                                sh_label = predict(Mdl,shuffled_measurements(logical(idxTest),:));
                                
                                shuffled_tested_events{ii,1}=sh_label{1};
                                
                            end
                            
                            %Calculate auROC
                            [X,Y,T,AUC] = perfcurve(events,scores',handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)});
                            auROC(1,time_point)=AUC-0.5;
                            %test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                            
                            
                            per_targets=these_all_which_events;
                            
                            test_out=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                            shuffled_out=zeros(length(handles.drgbchoices.events_to_discriminate),N);
                            for ii=1:N
                                for jj=1:length(handles.drgbchoices.events_to_discriminate)
                                    if strcmp(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(jj)},tested_events{ii})
                                        test_out(jj,ii)=1;
                                    else
                                        test_out(jj,ii)=0;
                                    end
                                    if strcmp(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(jj)},shuffled_tested_events{ii})
                                        shuffled_out(jj,ii)=1;
                                    else
                                        shuffled_out(jj,ii)=0;
                                    end
                                end
                            end
                            
                            test_out_per_timepoint(:,:,time_point)=test_out;
                            shuffled_out_per_timepoint(:,:,time_point)=shuffled_out;
                            discriminant_correct(1,time_point)=100*sum(sum(test_out.*per_targets))/N;
                            discriminant_correct_shuffled(1,time_point)=100*sum(sum(shuffled_out.*per_targets))/N;
                            fprintf(1, 'LDA for trough PAC wavelet power percent correct classification %d (for timepoint %d out of %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,length(t));
                            
                        end
                        
                        
                        batch_discriminant_correct_trough(ii_batch,:)=discriminant_correct(1,:);
                        batch_dimensionality_trough(ii_batch,:)=dimensionality(1,:);
                        batch_auROC_trough(ii_batch,:)=auROC;
                        batch_discriminant_correct_shuffled_trough(ii_batch,:)=discriminant_correct_shuffled(1,:);
                        batch_test_out_per_timepoint_trough(ii_batch,:,:,:)=test_out_per_timepoint;
                        batch_shuffled_out_per_timepoint_trough(ii_batch,:,:,:)=shuffled_out_per_timepoint;
                        
                        

                    end
                    
                    

                    
                    %Now save the data
                      %Calculate p value for the wavelet power
                    p_val_peak=zeros(1,length(t));
                    
                    for ii_t=1:length(t)
                        splus_out=zeros(1,no_batches*(N/2));
                        sminus_out=zeros(1,no_batches*(N/2));
                        
                        for ii_batch=1:no_batches
                            splus_out(1,(ii_batch-1)*(N/2)+1:ii_batch*(N/2))=batch_test_out_per_timepoint_peak(ii_batch,1,1:(N/2),ii_t);
                            sminus_out(1,(ii_batch-1)*(N/2)+1:ii_batch*(N/2))=batch_test_out_per_timepoint_peak(ii_batch,1,(N/2)+1:N,ii_t);
                        end
                        
                        p_val_peak(ii_t)=ranksum(splus_out,sminus_out);
                    end
                    
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_peak=p_val_peak;
  
                    handles_out.discriminant_PACpower_per_mouse_all_mice.group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                    
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_peak=zeros(1,length(t));
                    discriminant_correct=mean(batch_discriminant_correct_peak,1);
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_peak(1,:)=discriminant_correct(1,:);
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_discriminant_correct_peak=zeros(no_batches,length(t));
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_discriminant_correct_peak(:,:)=batch_discriminant_correct_peak;
                    
                    
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak=zeros(1,length(t));
                    dimensionality=mean(batch_dimensionality_peak,1);
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_peak(1,:)=dimensionality(1,:);
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_dimensionality_peak=zeros(no_batches,length(t));
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_dimensionality_peak(:,:)=batch_dimensionality_peak;
                    
                    
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_peak=zeros(1,length(t));
                    auROC=mean(batch_auROC_peak,1);
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_peak=auROC;
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_auROC_peak=zeros(no_batches,length(t));
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_auROC_peak(:,:)=batch_auROC_peak;
                    
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_peak=zeros(1,length(t));
                    discriminant_correct_shuffled=mean(batch_discriminant_correct_shuffled_peak,1);
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_peak(1,:)=discriminant_correct_shuffled(1,:);
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_peak=zeros(no_batches,length(t));
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_peak(:,:)=batch_discriminant_correct_shuffled_peak;
                    
                    
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_test_out_per_timepoint_peak=batch_test_out_per_timepoint_peak;
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_shuffled_out_per_timepoint_peak=batch_shuffled_out_per_timepoint_peak;
                    
                    
                    
                    figNo=figNo+1;
                    try
                        close(figNo)
                    catch
                    end
                    
                    hFig=figure(figNo)
                    set(hFig, 'units','normalized','position',[.1 .4 .75 .47])
                    
                    subplot(2,5,1)
                    hold on
                    
                    per95=prctile(discriminant_correct_shuffled(1,:),95);
                    per5=prctile(discriminant_correct_shuffled(1,:),5);
                    CIsh=[mean(discriminant_correct_shuffled(1,:))-per5 per95-mean(discriminant_correct_shuffled(1,:))]';
                    [hlCR, hpCR] = boundedline([t(1) t(end)],[mean(discriminant_correct_shuffled(1,:)) mean(discriminant_correct_shuffled(1,:))], CIsh', 'r');
                    
                    plot(t',discriminant_correct(1,:),'-k')
                    
                    %Odor on markers
                    plot([0 0],[0 100],'-k')
                    odorhl=plot([0 2.5],[20 20],'-k','LineWidth',5);
                    plot([2.5 2.5],[0 100],'-k')
                    
                    %title(['LDA % correct for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                    xlabel('Time (sec)')
                    ylabel('% correct peak')
                    
                    subplot(2,5,2)
                    hold on
                    
                    plot(t,auROC,'-b')
                    
                    %Odor on markers
                    plot([0 0],[-0.3 0.5],'-k')
                    odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                    plot([2.5 2.5],[0 0.5],'-k')
                    
                    %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                    xlabel('Time (sec)')
                    ylabel('auROC peak')
                    ylim([-0.3 0.6])
                    
                    
                    subplot(2,5,3)
                    hold on
                    
                    plot(t,dimensionality,'-b')
                    
                    maxdim=max(dimensionality);
                    mindim=min(dimensionality);
                    
                    %Odor on markers
                    plot([0 0],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                    odorhl=plot([0 2.5],[mindim mindim],'-k','LineWidth',5);
                    plot([2.5 2.5],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                    
                    %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                    xlabel('Time (sec)')
                    ylabel('dimensionality peak')
                    ylim([mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)])
                    

                    %Now do troughs
                      p_val_trough=zeros(1,length(t));
                    for ii_t=1:length(t)
                        splus_out=zeros(1,sum(these_all_which_events(1,:)==1));
                        splus_out(1,:)=test_out_per_timepoint(1,these_all_which_events(1,:)==1,ii_t);
                        sminus_out=zeros(1,sum(these_all_which_events(1,:)==0));
                        sminus_out(1,:)=test_out_per_timepoint(1,these_all_which_events(1,:)==0,ii_t);
                        p_val_trough(ii_t)=ranksum(splus_out,sminus_out);
                    end
                    
                    
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).p_val_trough=p_val_trough;
  
                    handles_out.discriminant_PACpower_per_mouse_all_mice.group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                    
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_trough=zeros(1,length(t));
                    discriminant_correct=mean(batch_discriminant_correct_trough,1);
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_trough(1,:)=discriminant_correct(1,:);
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_discriminant_correct_trough=zeros(no_batches,length(t));
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_discriminant_correct_trough(:,:)=batch_discriminant_correct_trough;
                    
                    
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough=zeros(1,length(t));
                    dimensionality=mean(batch_dimensionality_trough,1);
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).dimensionality_trough(1,:)=dimensionality(1,:);
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_dimensionality_trough=zeros(no_batches,length(t));
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_dimensionality_trough(:,:)=batch_dimensionality_trough;
                    
                    
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_trough=zeros(1,length(t));
                    auROC=mean(batch_auROC_trough,1);
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).auROC_trough=auROC;
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_auROC_trough=zeros(no_batches,length(t));
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_auROC_trough(:,:)=batch_auROC_trough;
                    
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_trough=zeros(1,length(t));
                    discriminant_correct_shuffled=mean(batch_discriminant_correct_shuffled_trough,1);
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_trough(1,:)=discriminant_correct_shuffled(1,:);
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_trough=zeros(no_batches,length(t));
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).discriminant_correct_shuffled_trough(:,:)=batch_discriminant_correct_shuffled_trough;
                    
                    
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_test_out_per_timepoint_trough=batch_test_out_per_timepoint_trough;
                    handles_out.discriminant_PACwavepower_all_mice.group(groupNo).percent_correct(percent_correct_ii).PACii(PACii).batch_shuffled_out_per_timepoint_trough=batch_shuffled_out_per_timepoint_trough;
          
                    subplot(2,5,6)
                    hold on
                    
                    per95=prctile(discriminant_correct_shuffled(1,:),95);
                    per5=prctile(discriminant_correct_shuffled(1,:),5);
                    CIsh=[mean(discriminant_correct_shuffled(1,:))-per5 per95-mean(discriminant_correct_shuffled(1,:))]';
                    [hlCR, hpCR] = boundedline([t(1) t(end)],[mean(discriminant_correct_shuffled(1,:)) mean(discriminant_correct_shuffled(1,:))], CIsh', 'r');
                    
                    plot(t',discriminant_correct(1,:),'-k')
                    
                    %Odor on markers
                    plot([0 0],[0 100],'-k')
                    odorhl=plot([0 2.5],[20 20],'-k','LineWidth',5);
                    plot([2.5 2.5],[0 100],'-k')
                    
                    %title(['LDA % correct for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str_all_mice ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                    xlabel('Time (sec)')
                    ylabel('% correct trough')
                    
                    subplot(2,5,7)
                    hold on
                    
                    plot(t,auROC,'-b')
                    
                    %Odor on markers
                    plot([0 0],[0 0.5],'-k')
                    odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                    plot([2.5 2.5],[0 0.5],'-k')
                    
                    %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str_all_mice ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                    xlabel('Time (sec)')
                    ylabel('auROC trough')
                    ylim([-0.3 0.6])
                    
                    %Plot dimensionality
                    subplot(2,5,8)
                    hold on
                    
                    plot(t,dimensionality,'-b')
                    
                    mindim=min(dimensionality);
                    maxdim=max(dimensionality);
                    
                    %Odor on markers
                    plot([0 0],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                    odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                    plot([2.5 2.5],[mindim-0.1*(maxdim-mindim) maxdim+0.1*(maxdim-mindim)],'-k')
                    
                    %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{PACii} ' mouse No ' num2str_all_mice ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                    xlabel('Time (sec)')
                    ylabel('Dimensionality trough')
                    xlim([-2 5])
                    
                    
                    %Plot p value
                    subplot(2,5,[4,5,9,10])
                    hold on
                    
                    p1=plot(t,log10(p_val_lick),'-k');
                    p2=plot(t,log10(p_val_trough),'-b');
                    p3=plot(t,log10(p_val_peak),'-m');
                    plot([t(1) t(end)],[log10(0.05) log10(0.05)],'-r')
                    legend([p1 p2 p3],{'Licks','Trough','Peak'})
                    ylabel('log(p)')
                    xlabel('Time (sec)')
                    xlim([-2 5])
                    ylim([-10 0])
                    
                    suptitle(['PAC wavelet power LDA analysis for Theta/' handles.drgbchoices.PACnames{PACii} ' PAC ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                    
                     if (groupNo==1)&(PACii==3)&(percent_correct_ii==1)
                        pffft=1;
                    end
            end
        end
    end
    save(output_name,'handles_out','-v7.3')
fprintf(1, ['LDA for phase-referenced wavelet power processed for' output_name 'groupNo %d'],groupNo);
end

save(output_name,'handles_out','-v7.3')
fprintf(1, ['LDA for phase-referenced wavelet power saved for' output_name]);
                            




pffft1=1;


function drgLFPDiscriminantBatch
%Perceptron analysis

close all
clear all

first_file=1;
figNo=0;

[choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoicesDiscriminant*.m'},'Select the .m file with all the choices for discriminant analysis');
fprintf(1, ['\ndrgLFPDiscriminantBatch run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

new_no_files=handles.drgbchoices.no_files;
choicePathName=handles.drgbchoices.PathName;
choiceFileName=handles.drgbchoices.FileName;
%Very, very important!
handles.evTypeNo=handles.drgbchoices.referenceEvent;

%Find out if there is already an output file
overwrite_out_file=1;
if exist([handles.drgb.outPathName 'Discriminant_' handles.drgb.outFileName])==2
    answer = questdlg('An output file exists, do you want to overwrite?', ...
        'Output file exists', ...
        'Yes','No','No');
    if strcmp(answer,'Yes')
        overwrite_out_file=1;
    else
        overwrite_out_file=0;
        load([handles.drgb.outPathName 'Discriminant_' handles.drgb.outFileName])
    end
end

%Parallel batch processing for each file
all_files_present=1;
for filNum=first_file:handles.drgbchoices.no_files
    
    %Make sure that all the files exist
    jtFileName=handles.drgbchoices.FileName{filNum};
    if iscell(handles.drgbchoices.PathName)
        jtPathName=handles.drgbchoices.PathName{filNum};
    else
        jtPathName=handles.drgbchoices.PathName;
    end
    if exist([jtPathName jtFileName])==0
        fprintf(1, ['Program will be terminated because file No %d, ' jtPathName jtFileName ' does not exist\n'],filNum);
        all_files_present=0;
    end
    
    if (exist( [jtPathName jtFileName(10:end-4) '.dg'])==0)&(exist( [jtPathName jtFileName(10:end-4) '.rhd'])==0)
        fprintf(1, ['Program will be terminated because neither dg or rhd files for file No %d, ' [jtPathName jtFileName(10:end-4)] ' does not exist\n'],filNum);
        all_files_present=0;
    end
    
end

if all_files_present==1
    
    
    handles_out.drgbchoices=handles.drgbchoices;
    
    handles.displayData=0;
    
    handles.subtractRef=handles.drgbchoices.subtractRef;
    handles.evTypeNo=handles.drgbchoices.referenceEvent; %Process all OdorOn trials
    
    no_files=handles.drgbchoices.no_files;
    
    for mouseNo=1:max(handles.drgbchoices.mouse_no)
        for groupNo=1:max(handles.drgbchoices.group_no)
            
            fprintf(1, ['Mouse no %d ' handles.drgbchoices.group_no_names{groupNo} '\n\n'],mouseNo);
            
            if overwrite_out_file==1
                process_data_for_mouse=1;
            else
                %Find out whether this mouse/group was processed
                try
                    dcalc=handles_out.discriminant_per_mouse(mouseNo).group(groupNo);
                    process_data_for_mouse=0;
                catch
                    process_data_for_mouse=1;
                end
            end
            
            if process_data_for_mouse==1
                no_trials=0;
                mouse_has_data=0;
                
                for filNum=first_file:no_files
                    if (handles.drgbchoices.mouse_no(filNum)==mouseNo)&(handles.drgbchoices.group_no(filNum)==groupNo)
                        %     for filNum=first_file:handles.drgbchoices.no_files
                        if mouse_has_data==0
                            first_file_for_this_mouse=filNum;
                        end
                        mouse_has_data=1;
                        file_no=filNum
                        
                        this_jt=handles.drgbchoices.FileName{filNum};
                        
                        
                        %Othrwise read the jt_times and do processing
                        %read the jt_times file
                        jtFileName=handles.drgbchoices.FileName{filNum};
                        if iscell(handles.drgbchoices.PathName)
                            jtPathName=handles.drgbchoices.PathName{filNum};
                        else
                            jtPathName=handles.drgbchoices.PathName;
                        end
                        
                        drgRead_jt_times(jtPathName,jtFileName);
                        FileName=[jtFileName(10:end-4) '_drg.mat'];
                        fullName=[jtPathName,FileName];
                        my_drg={'drg'};
                        S=load(fullName,my_drg{:});
                        handles.drg=S.drg;
                        handles_out.drg=handles.drg;
                        
                        switch handles.drg.session(handles.sessionNo).draq_p.dgordra
                            case 1
                            case 2
                                handles.drg.drta_p.fullName=[jtPathName jtFileName(10:end-4) '.dg'];
                            case 3
                                handles.drg.drta_p.fullName=[jtPathName jtFileName(10:end-4) '.rhd'];
                        end
                        
                        %Set the last trial to the last trial in the session
                        
                        handles.time_start=min(handles.drgbchoices.timeStart-handles.time_pad);
                        handles.time_end=max(handles.drgbchoices.timeEnd+handles.time_pad);
                        
                        handles.burstLowF=handles.LFPPowerSpectrumLowF;
                        handles.burstHighF=handles.LFPPowerSpectrumHighF;
                        handles.lastTrialNo=handles.drg.session(handles.sessionNo).events(2).noTimes;
                        handles.trialNo=1;
                        
                        
                        no_elect=length(handles.drgbchoices.which_electrodes);
                        for LFPNo=1:no_elect
                            for bwii=1:length(handles.drgbchoices.lowF)
                                par_out(LFPNo).bwii(bwii).log_P_timecourse=[];
                            end
                            par_out(LFPNo).length_trial_no=0;
                            par_out(LFPNo).t=[];
                            par_out(LFPNo).length_trial_no=0;
                            par_out(LFPNo).this_trialNo=[];
                            par_out(LFPNo).which_event=[];
                            par_out(LFPNo).perCorr_pertr=[];
                            par_out(LFPNo).valid_trials=[];
                        end
                        
                        gcp
                        
                        parfor LFPNo=1:no_elect
                            handlespf=struct();
                            handlespf=handles;
                            this_peakLFPNo=handlespf.drgbchoices.which_electrodes(LFPNo);
                            handlespf.peakLFPNo=this_peakLFPNo;
                            all_Power=[];
                            all_Power_ref=[];
                            all_Power_timecourse=[];
                            perCorr_pertr=[];
                            which_event=[];
                            [t,f,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event]=drgGetLFPPowerForThisEvTypeNo(handlespf);
                            
                            %Note: This is here because for a very small number
                            %of trials a tiral is included that is neigher S+
                            %nor S-
                            these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),length(this_trialNo));
                            for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                these_all_which_events(ii,:)= which_event(kk,:);
                            end
                            
                            valid_trials=logical(sum(these_all_which_events,1));
                            
                            old_all_Power=all_Power;
                            szap=size(old_all_Power);
                            all_Power=zeros(sum(valid_trials),szap(2));
                            all_Power(:,:)=old_all_Power(valid_trials,:);
                            
                            old_all_Power_ref=all_Power_ref;
                            szapr=size(old_all_Power_ref);
                            all_Power_ref=zeros(sum(valid_trials),szapr(2));
                            all_Power_ref(:,:)=old_all_Power_ref(valid_trials,:);
                            
                            old_all_Power_timecourse=all_Power_timecourse;
                            szapt=size(old_all_Power_timecourse);
                            all_Power_timecourse=zeros(sum(valid_trials),szapt(2),szapt(3));
                            all_Power_timecourse(:,:,:)= old_all_Power_timecourse(valid_trials,:,:);
                            
                            old_this_trialNo=this_trialNo;
                            this_trialNo=zeros(1,sum(valid_trials));
                            this_trialNo(:,:)=old_this_trialNo(1,valid_trials);
                            
                            old_perCorr_pertr=perCorr_pertr;
                            perCorr_pertr=zeros(1,sum(valid_trials));
                            perCorr_pertr(:,:)=old_perCorr_pertr(1,valid_trials);
                            
                            old_which_event=which_event;
                            szwe=size(old_which_event);
                            which_event=zeros(szwe(1),sum(valid_trials));
                            which_event(:,:)=old_which_event(:,valid_trials);
                            
                            
                            
                            fprintf(1, 'For file no %d, LFP no %d the number of trials included is %d (out of %d)\n',filNum,this_peakLFPNo,length(this_trialNo),handlespf.lastTrialNo);
                            
                            freq=f';
                            
                            par_out(LFPNo).length_trial_no=length(this_trialNo);
                            par_out(LFPNo).this_trialNo=this_trialNo;
                            par_out(LFPNo).t=t;
                            par_out(LFPNo).which_event=which_event;
                            par_out(LFPNo).perCorr_pertr=perCorr_pertr;
                            par_out(LFPNo).valid_trials=valid_trials;
                            
                            %Save timecourse
                            for bwii=1:length(handles.drgbchoices.lowF)
                                if handlespf.subtractRef==0
                                    par_out(LFPNo).bwii(bwii).log_P_timecourse=zeros(length(this_trialNo),length(t));
                                    par_out(LFPNo).bwii(bwii).log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse(:,(freq>handlespf.drgbchoices.lowF(bwii))&(freq<=handlespf.drgbchoices.highF(bwii)),:)),2);
                                else
                                    par_out(LFPNo).bwii(bwii).log_P_timecourse=zeros(length(this_trialNo),length(t));
                                    par_out(LFPNo).bwii(bwii).log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse(:,(freq>handlespf.drgbchoices.lowF(bwii))&(freq<=handlespf.drgbchoices.highF(bwii)),:)),2);
                                    log_P_timecourse_ref=zeros(length(this_trialNo),length(t));
                                    log_P_timecourse_ref(:,:)=repmat(mean(10*log10(all_Power_ref(:,(freq>handlespf.drgbchoices.lowF(bwii))&(freq<=handlespf.drgbchoices.highF(bwii)))),2),1,length(t));
                                    par_out(LFPNo).bwii(bwii).log_P_timecourse=par_out(LFPNo).bwii(bwii).log_P_timecourse-log_P_timecourse_ref;
                                end
                            end
                            
                        end
                        
                        t=par_out(1).t;
                        
                        for LFPNo=1:length(handles.drgbchoices.which_electrodes)
                            
%                             this_peakLFPNo=handles.drgbchoices.which_electrodes(LFPNo);
                            if LFPNo>1
                                if last_No_trials~=par_out(LFPNo).length_trial_no
                                    fprintf(1, 'WARNING. For file no %d the number of trials differ between different LFPs\n ',filNum);
                                end
                            end
                            last_No_trials=par_out(LFPNo).length_trial_no;
                            
                            if LFPNo==1
                                if no_trials==0
                                    all_log_P_timecourse=zeros(length(handles.drgbchoices.lowF),length(handles.drgbchoices.which_electrodes),par_out(LFPNo).length_trial_no,length(t));
                                else
                                    this_all_log_P_timecourse=[];
                                    this_all_log_P_timecourse=all_log_P_timecourse;
                                    all_log_P_timecourse=zeros(length(handles.drgbchoices.lowF),length(handles.drgbchoices.which_electrodes),par_out(LFPNo).length_trial_no+no_trials,length(t));
                                    all_log_P_timecourse(:,:,1:no_trials,:)=this_all_log_P_timecourse;
                                end
                            end
                            for bwii=1:length(handles.drgbchoices.lowF)
                                all_log_P_timecourse(bwii,LFPNo,no_trials+1:no_trials+par_out(LFPNo).length_trial_no,:)=par_out(LFPNo).bwii(bwii).log_P_timecourse;
                            end
                        end
                        
                        %Calculate licks
                        stamped_lick_ii=[];
                        these_stamped_lick_times=[];
                        no_trials_l=[];
                        trials_included_l=[];
                        
                        [lick_freq,times_lick_freq,lick_traces,CIlickf,lick_trace_times,stamped_lick_ii,these_stamped_lick_times,no_trials_l,trials_included_l]=drgGetLicks(handles);
                        
                        
                        
                        %Save all_which_events and all_perCorr_pertr
                        if filNum==first_file_for_this_mouse
                            all_which_events=zeros(length(handles.drgbchoices.evTypeNos),par_out(1).length_trial_no);
                            all_which_events(:,:)=par_out(1).which_event;
                            
                            all_perCorr_pertr=zeros(1,par_out(1).length_trial_no);
                            all_perCorr_pertr(:,:)=par_out(1).perCorr_pertr;
                            
                            all_stamped_lick_times=zeros(par_out(1).length_trial_no,250);
                            sztslt=size(these_stamped_lick_times);
                            all_stamped_lick_ii=zeros(1,par_out(1).length_trial_no);
                            for ii=1:par_out(1).length_trial_no
                                kk=find(trials_included_l==par_out(1).this_trialNo(ii));
                                all_stamped_lick_times(ii,1:sztslt(2))=these_stamped_lick_times(kk,:);
                                all_stamped_lick_ii(1,ii)=stamped_lick_ii(kk);
                            end
                            
                        else
                            this_all_which_events=[];
                            this_all_which_events=all_which_events;
                            all_which_events=zeros(length(handles.drgbchoices.evTypeNos),par_out(1).length_trial_no+no_trials);
                            all_which_events(:,1:no_trials,:)=this_all_which_events;
                            all_which_events(:,no_trials+1:no_trials+par_out(1).length_trial_no,:)=par_out(1).which_event;
                            
                            this_all_perCorr_pertr=[];
                            this_all_perCorr_pertr=all_perCorr_pertr;
                            all_perCorr_pertr=zeros(1,par_out(1).length_trial_no+no_trials);
                            all_perCorr_pertr(:,1:no_trials)=this_all_perCorr_pertr;
                            all_perCorr_pertr(:,no_trials+1:no_trials+par_out(1).length_trial_no)=par_out(1).perCorr_pertr;
                            
                            this_all_stamped_lick_times=[];
                            this_all_stamped_lick_times=all_stamped_lick_times;
                            all_stamped_lick_times=zeros(par_out(1).length_trial_no+no_trials,250);
                            all_stamped_lick_times(1:no_trials,:)=this_all_stamped_lick_times;
                            sztslt=size(these_stamped_lick_times);
                            
                            this_all_stamped_lick_ii=[];
                            this_all_stamped_lick_ii=all_stamped_lick_ii;
                            all_stamped_lick_ii=zeros(1,par_out(1).length_trial_no+no_trials);
                            all_stamped_lick_ii(1,1:no_trials)=this_all_stamped_lick_ii;
                            
                            for ii=1:par_out(1).length_trial_no
                                kk=find(trials_included_l==par_out(1).this_trialNo(ii));
                                all_stamped_lick_times(ii+no_trials,1:sztslt(2))=these_stamped_lick_times(kk,:);
                                all_stamped_lick_ii(1,no_trials+ii)=stamped_lick_ii(1,kk);
                            end
                            
                        end
                        no_trials=no_trials+par_out(1).length_trial_no;
                        
                    end
                end
                
                
                %Do discriminant analysis for the percent correct windows
                if mouse_has_data==1
                    
                    %Calculate percent correct for each bandwidth
                    for bwii=1:length(handles.drgbchoices.lowF)
                        
                        for percent_correct_ii=1:length(handles.drgbchoices.per_lab)
                            
                            discriminant_correct=zeros(1,length(t));
                            discriminant_correct_shuffled=zeros(1,length(t));
                            auROC=zeros(1,length(t));
                            
                            these_per_corr=(all_perCorr_pertr>=handles.drgbchoices.percent_windows(percent_correct_ii,1))...
                                &(all_perCorr_pertr<=handles.drgbchoices.percent_windows(percent_correct_ii,2));
                            N=sum(these_per_corr);
                            
                            %Do the analysis only if there are more than 20 trials
                            if N>=20
                                
                                if sum(handles.drgbchoices.which_discriminant==1)>0
                                    
                                    
                                    %Do perceptron prediction analysis for every point in the timecourse
                                    num_traces=length(handles.drgbchoices.which_electrodes); %Number of electrodes
                                    num_trials=sum(these_per_corr);
                                    per_targets=zeros(2,sum(these_per_corr));
                                    %S+
                                    per_targets(1,:)=all_which_events(2,these_per_corr);
                                    %S-
                                    per_targets(2,:)=all_which_events(5,these_per_corr);
                                    
                                    
                                    these_all_log_P_timecourse=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(t));
                                    these_all_log_P_timecourse(:,:,:)=all_log_P_timecourse(bwii,:,these_per_corr,:);
                                    
                                    test_out_per_timepoint=zeros(2,N,length(t));
                                    shuffled_out_per_timepoint=zeros(2,N,length(t));
                                    
                                    
                                    gcp
                                    
                                    fprintf(1, ['Perceptron processed for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
                                    
                                    parfor time_point=1:length(t)
                                        %         for time_point=1:length(t)
                                        
                                        per_input=zeros(length(handles.drgbchoices.which_electrodes),N);
                                        per_input(:,:)=these_all_log_P_timecourse(:,:,time_point);
                                        fprintf(1, '\nTime point %d: ',time_point);
                                        
                                        %Perceptron
                                        %leave one out
                                        test_out=zeros(2,N);
                                        shuffled_out=zeros(2,N);
                                        
                                        %Do perceptron analysis for each trial
                                        for ii=1:N
                                            fprintf(1, '%d ',ii);
                                            %Create input and target vectors leaving one trial out
                                            %For per_input each column has the dF/F for one trial
                                            %each row is a single time point for dF/F for one of the cells
                                            %For per_target the top row is 1 if the odor is S+ and 0 if it is
                                            %S-, and row 2 has 1 for S-
                                            this_per_input=[];
                                            this_per_targets=[];
                                            if ii==1
                                                this_per_input=per_input(:,2:end);
                                                this_per_targets=per_targets(:,2:end);
                                            else
                                                if ii==N
                                                    this_per_input=per_input(:,1:end-1);
                                                    this_per_targets=per_targets(:,1:end-1);
                                                else
                                                    this_per_input=[per_input(:,1:ii-1) per_input(:,ii+1:end)];
                                                    this_per_targets=[per_targets(:,1:ii-1) per_targets(:,ii+1:end)];
                                                end
                                            end
                                            
                                            
                                            
                                            %Create a net with the default perceptron
                                            net=perceptron;
                                            
                                            % Set up Division of Data for Training, Validation, Testing
                                            net.divideParam.trainRatio = 1;
                                            net.divideParam.valRatio = 0;
                                            net.divideParam.testRatio = 0;
                                            net.trainParam.showWindow = 0;
                                            
                                            % Train the Network
                                            [net,tr] = train(net,this_per_input,this_per_targets);
                                            
                                            %Calculate the trial that was left out
                                            one_out = per_input(:,ii);
                                            test_out(:,ii) = net(one_out);
                                            
                                            %Calculate a shuffled trial
                                            
                                            
                                            shuffled_trials=ceil(N*rand(1,num_traces));
                                            
                                            one_shuffled=zeros(num_traces,1);
                                            for jj=1:num_traces
                                                one_shuffled(jj,1)=per_input(jj,shuffled_trials(jj));
                                            end
                                            
                                            shuffled_out(:,ii) = net(one_shuffled);
                                            
                                        end
                                        test_out_per_timepoint(:,:,time_point)=test_out;
                                        shuffled_out_per_timepoint(:,:,time_point)=shuffled_out;
                                        discriminant_correct(1,time_point)=100*sum(per_targets(1,:)==test_out(1,:))/N;
                                        discriminant_correct_shuffled(1,time_point)=100*sum(per_targets(1,:)==shuffled_out(1,:))/N;
                                        fprintf(1, '\nPerceptron percent correct classification %d (for timepoint %d out of %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,length(t));
                                    end
                                end
                                
                                if sum(handles.drgbchoices.which_discriminant==2)>0
                                    %Linear discriminant analysis
                                    
                                    %Number of trials
                                    %                                         N=sum(these_per_corr);
                                    
                                    these_all_log_P_timecourse=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(t));
                                    these_all_log_P_timecourse(:,:,:)=all_log_P_timecourse(bwii,:,these_per_corr,:);
                                    
                                    these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                                    for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                        kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                        these_all_which_events(ii,:)= all_which_events(kk,these_per_corr);
                                    end
                                    
                                    test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                                    shuffled_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                                    
                                    fprintf(1, ['LDA processed for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
                                    fprintf(1,'For these events: ')
                                    for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                        fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                                    end
                                    fprintf(1,'\n')
                                    
                                    gcp
                                    
                                    for time_point=1:length(t)
                                        
                                        %LFP power per trial per electrode
                                        measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                        measurements(:,:)=these_all_log_P_timecourse(:,:,time_point)';
                                        
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
                                            shuffled_measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
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
                                        
                                        per_targets=zeros(length(handles.drgbchoices.events_to_discriminate),N);
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
                                        fprintf(1, 'LDA percent correct classification %d (for timepoint %d out of %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/N,time_point,length(t));
                                        
                                    end
                                end
                                
                                
                                
                                if (sum(handles.drgbchoices.which_discriminant==1)>0)||(sum(handles.drgbchoices.which_discriminant==2)>0)
                                    
                                    
                                    figNo=figNo+1
                                    try
                                        close(figNo)
                                    catch
                                    end
                                    
                                    figure(figNo)
                                    
                                    subplot(1,2,1)
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
                                    
                                    %title(['LDA % correct for ' handles.drgbchoices.bwlabels{bwii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                    xlabel('Time (sec)')
                                    ylabel('Percent correct')
                                    
                                    subplot(1,2,2)
                                    hold on
                                    
                                    plot(t,auROC)
                                    
                                    %Odor on markers
                                    plot([0 0],[0 0.5],'-k')
                                    odorhl=plot([0 2.5],[0 0],'-k','LineWidth',5);
                                    plot([2.5 2.5],[0 0.5],'-k')
                                    
                                    %title(['auROC for LDA for ' handles.drgbchoices.bwlabels{bwii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                    xlabel('Time (sec)')
                                    ylabel('auROC')
                                    
                                    suptitle(['LDA analysis for ' handles.drgbchoices.bwlabels{bwii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                    
                                    
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).discriminant_correct=zeros(1,length(t));
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).discriminant_correct(1,:)=discriminant_correct(1,:);
                                    
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).auROC=zeros(1,length(t));
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).auROC=auROC;
                                    
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).discriminant_correct_shuffled=zeros(1,length(t));
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).discriminant_correct_shuffled(1,:)=discriminant_correct_shuffled(1,:);
                                    
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).test_out_per_timepoint=test_out_per_timepoint;
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).shuffled_out_per_timepoint=shuffled_out_per_timepoint;
                                    
                                    handles_out.t=t';
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).no_trials=N;
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).per_targets=per_targets;
                                    
                                    these_all_which_events=[];
                                    these_all_which_events=all_which_events(:,these_per_corr);
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                                    
                                    these_all_stamped_lick_times=[];
                                    these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                                    
                                    these_all_stamped_lick_ii=[];
                                    these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                end
                                
                                if sum(handles.drgbchoices.which_discriminant==3)>0
                                    
                                    %PCA
                                    
                                    these_all_log_P_timecourse=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(t));
                                    these_all_log_P_timecourse(:,:,:)=all_log_P_timecourse(bwii,:,these_per_corr,:);
                                    
                                    these_all_which_events=zeros(length(handles.drgbchoices.events_to_discriminate),sum(these_per_corr));
                                    for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                        kk=find(handles.drgbchoices.evTypeNos==handles.drgbchoices.events_to_discriminate(ii));
                                        these_all_which_events(ii,:)= all_which_events(kk,these_per_corr);
                                    end
                                    
                                    test_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                                    shuffled_out_per_timepoint=zeros(length(handles.drgbchoices.events_to_discriminate),N,length(t));
                                    
                                    fprintf(1, ['PCA processed for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d trials\n'],mouseNo,N);
                                    fprintf(1,'For these events: ')
                                    for ii=1:length(handles.drgbchoices.events_to_discriminate)
                                        fprintf(1,[handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(ii)} ' '])
                                    end
                                    fprintf(1,'\n')
                                    
                                    for time_point=1:length(t)
                                        par_t_out(time_point).principal_components=zeros(N,length(handles.drgbchoices.which_electrodes));
                                    end
                                    
                                    gcp
                                    
                                    for time_point=1:length(t)
                                        
                                        %LFP power per trial per electrode
                                        measurements=zeros(N,length(handles.drgbchoices.which_electrodes));
                                        measurements(:,:)=these_all_log_P_timecourse(:,:,time_point)';
                                        
                                        %Do the PCA
                                        [coeff,par_t_out(time_point).principal_components,latent]=pca(measurements);
                                        
                                    end
                                    
                                    
                                    %Show a figure of the PCA and record the output
                                    
                                    principal_components=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                    for time_point=1:length(t)
                                        principal_components(time_point,:,:)=par_t_out(time_point).principal_components;
                                    end
                                    
                                    %Show the result of the PCA
                                    figNo=figNo+1
                                    try
                                        close(figNo)
                                    catch
                                    end
                                    
                                    figure(figNo)
                                    
                                    %Show PCA before odor on
                                    subplot(2,2,3)
                                    hold on
                                    
                                    these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
                                    these_pcs(:,:)=principal_components(6,:,:);
                                    plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'or')
                                    plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'ob')
                                    legend(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(1)},handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)})
                                    xlabel('PC1')
                                    ylabel('PC2')
                                    title('-1 sec')
                                    
                                    %Show PCA after odor
                                    subplot(2,2,4)
                                    hold on
                                    
                                    these_pcs=zeros(N,length(handles.drgbchoices.which_electrodes));
                                    these_pcs(:,:)=principal_components(41,:,:);
                                    plot(these_pcs(logical(these_all_which_events(1,:)),1),these_pcs(logical(these_all_which_events(1,:)),2),'or')
                                    plot(these_pcs(logical(these_all_which_events(2,:)),1),these_pcs(logical(these_all_which_events(2,:)),2),'ob')
                                    legend(handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(1)},handles.drg.draq_d.eventlabels{handles.drgbchoices.events_to_discriminate(2)})
                                    xlabel('PC1')
                                    ylabel('PC2')
                                    title('2.5 sec')
                                    
                                    %Show the timecourse for PC1
                                    subplot(2,2,[1,2])
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
                                    title(['PC1 for ' handles.drgbchoices.bwlabels{bwii} ' mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' ' handles.drgbchoices.group_no_names{groupNo}])
                                    xlabel('Time (sec)')
                                    ylabel('PC1')
                                    
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=1;
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).principal_components=zeros(length(t),N,length(handles.drgbchoices.which_electrodes));
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).principal_components=principal_components;
                                    
                                    handles_out.t=t';
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).bwii(bwii).no_trials=N;
                                    
                                    these_all_which_events=[];
                                    these_all_which_events=all_which_events(:,these_per_corr);
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                                    
                                    these_all_stamped_lick_times=[];
                                    these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                                    
                                    these_all_stamped_lick_ii=[];
                                    these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                                    
                                end
                            else
                                handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=0;
                                fprintf(1, ['LDA/PCA not processed for mouse No %d ' handles.drgbchoices.group_no_names{groupNo} ' ' handles.drgbchoices.per_lab{percent_correct_ii} ' because there were only %d trials (fewer than 20 trials)\n'],mouseNo,N);
                            end
                        end
                    end
                else
                    handles_out.discriminant_per_mouse(mouseNo).group(groupNo).percent_correct(percent_correct_ii).discriminant_calculated=0;
                end
                save([handles.drgb.outPathName 'Discriminant_' handles.drgb.outFileName],'handles_out','-v7.3')
            end
        end
    end
end



pffft1=1;


function drgLFPspectPerceptronBatch
%Perceptron analysis

first_file=1;
figNo=0;

[choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoicesPerceptron*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgRunBatchLFP run for ' choiceFileName '\n\n']);

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
    handles_out.no_mice_with_data=0;
    handles_out.drgbchoices=handles.drgbchoices;
    handles.displayData=0;
    
    handles.subtractRef=handles.drgbchoices.subtractRef;
    handles.evTypeNo=handles.drgbchoices.referenceEvent; %Process all OdorOn trials
    
    no_files=handles.drgbchoices.no_files;
    
    for mouseNo=1:max(handles.drgbchoices.mouse_no)
        fprintf(1, 'Mouse no %d\n\n',mouseNo);
                   
        no_trials=0;
        mouse_has_data=0;
        for filNum=first_file:no_files
            if handles.drgbchoices.mouse_no(filNum)==mouseNo
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
                
                for LFPNo=1:length(handles.drgbchoices.which_electrodes)
                    this_peakLFPNo=handles.drgbchoices.which_electrodes(LFPNo);
                    handles.peakLFPNo=this_peakLFPNo;
                    %[t,freq,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event]=drgGetLFPwavePowerForThisEvTypeNo(handles);
                    all_Power=[];
                    all_Power_ref=[];
                    all_Power_timecourse=[];
                    perCorr_pertr=[];
                    which_event=[];
                    [t,f,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event]=drgGetLFPPowerForThisEvTypeNo(handles);
                    
                    fprintf(1, 'For file no %d, LFP no %d the number of trials included is %d (out of %d)\n',filNum,this_peakLFPNo,length(this_trialNo),handles.lastTrialNo);
                    if this_peakLFPNo>1
                        if last_No_trials~=length(this_trialNo)
                            fprintf(1, 'WARNING. For file no %d the number of trials differ between different LFPs\n ',filNum,length(this_trialNo));
                        end
                    end
                    last_No_trials=length(this_trialNo);
                    freq=f';
                    
                    %Save timecourse
                    if handles.subtractRef==0
                        log_P_timecourse=zeros(length(this_trialNo),length(t));
                        log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),2);
                    else
                        log_P_timecourse=zeros(length(this_trialNo),length(t));
                        log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),2);
                        log_P_timecourse_ref=zeros(length(this_trialNo),length(t));
                        log_P_timecourse_ref(:,:)=repmat(mean(10*log10(all_Power_ref),2),1,length(t));
                        log_P_timecourse=log_P_timecourse-log_P_timecourse_ref;
                    end
                    if this_peakLFPNo==1
                        if no_trials==0
                            all_log_P_timecourse=zeros(length(handles.drgbchoices.which_electrodes),length(this_trialNo),length(t));
                        else
                            this_all_log_P_timecourse=[];
                            this_all_log_P_timecourse=all_log_P_timecourse;
                            all_log_P_timecourse=zeros(length(handles.drgbchoices.which_electrodes),length(this_trialNo)+no_trials,length(t));
                            all_log_P_timecourse(:,1:no_trials,:)=this_all_log_P_timecourse;
                        end
                    end
                    all_log_P_timecourse(this_peakLFPNo,no_trials+1:no_trials+length(this_trialNo),:)=log_P_timecourse;  
                end
                
                %Calculate licks
                [lick_freq,times_lick_freq,lick_traces,CIlickf,lick_trace_times,stamped_lick_ii,these_stamped_lick_times,no_trials_l,trials_included_l]=drgGetLicks(handles);
                
                %Save all_which_events and all_perCorr_pertr
                
                if filNum==first_file_for_this_mouse
                    all_which_events=zeros(length(handles.drgbchoices.evTypeNos),length(this_trialNo));
                    all_which_events(:,:)=which_event;
                    
                    all_perCorr_pertr=zeros(1,length(this_trialNo));
                    all_perCorr_pertr(:,:)=perCorr_pertr;
                    
                    all_stamped_lick_times=zeros(length(this_trialNo),250);
                    sztslt=size(these_stamped_lick_times);
                    all_stamped_lick_ii=zeros(1,length(this_trialNo));
                    for ii=1:length(this_trialNo)
                        kk=find(trials_included_l==this_trialNo(ii));
                        all_stamped_lick_times(ii,1:sztslt(2))=these_stamped_lick_times(kk,:);
                        all_stamped_lick_ii(1,ii)=stamped_lick_ii(kk);
                    end
                    
                else
                    this_all_which_events=[];
                    this_all_which_events=all_which_events;
                    all_which_events=zeros(length(handles.drgbchoices.evTypeNos),length(this_trialNo)+no_trials);
                    all_which_events(:,1:no_trials,:)=this_all_which_events;
                    all_which_events(:,no_trials+1:no_trials+length(this_trialNo),:)=which_event;
                    
                    this_all_perCorr_pertr=[];
                    this_all_perCorr_pertr=all_perCorr_pertr;
                    all_perCorr_pertr=zeros(1,length(this_trialNo)+no_trials);
                    all_perCorr_pertr(:,1:no_trials)=this_all_perCorr_pertr;
                    all_perCorr_pertr(:,no_trials+1:no_trials+length(this_trialNo))=perCorr_pertr;
                    
                    this_all_stamped_lick_times=[];
                    this_all_stamped_lick_times=all_stamped_lick_times;
                    all_stamped_lick_times=zeros(length(this_trialNo)+no_trials,250);
                    all_stamped_lick_times(1:no_trials,:)=this_all_stamped_lick_times;
                    sztslt=size(these_stamped_lick_times);
                    
                    this_all_stamped_lick_ii=[];
                    this_all_stamped_lick_ii=all_stamped_lick_ii;
                    all_stamped_lick_ii=zeros(1,length(this_trialNo)+no_trials);
                    all_stamped_lick_ii(1,1:no_trials)=this_all_stamped_lick_ii;
                    
                    for ii=1:length(this_trialNo)
                        kk=find(trials_included_l==this_trialNo(ii));
                        all_stamped_lick_times(ii+no_trials,1:sztslt(2))=these_stamped_lick_times(kk,:);
                        all_stamped_lick_ii(1,no_trials+ii)=stamped_lick_ii(1,kk);
                    end
                    
                end
                no_trials=no_trials+length(this_trialNo);
                
                
                
            end
        end
        
        %Now do perceptron analysis for the percent correct windows
        
        
        if mouse_has_data==1
            
            handles_out.no_mice_with_data=handles_out.no_mice_with_data+1;
            handles_out.mice_with_data(handles_out.no_mice_with_data)=mouseNo;
            perceptron_correct=zeros(length(handles.drgbchoices.per_lab),length(t));
            perceptron_correct_shuffled=zeros(length(handles.drgbchoices.per_lab),length(t));
            for percent_correct_ii=1:length(handles.drgbchoices.per_lab)
                
                these_per_corr=(all_perCorr_pertr>=handles.drgbchoices.percent_windows(percent_correct_ii,1))...
                    &(all_perCorr_pertr<=handles.drgbchoices.percent_windows(percent_correct_ii,2));
                
                %Do perceptron prediction analysis for every point in the timecourse
                num_traces=length(handles.drgbchoices.which_electrodes); %Number of electrodes
                num_trials=sum(these_per_corr);
                per_targets=zeros(2,sum(these_per_corr));
                %S+
                per_targets(1,:)=all_which_events(2,these_per_corr);
                %S-
                per_targets(2,:)=all_which_events(5,these_per_corr);
                per_ii=sum(these_per_corr);
                
                these_all_log_P_timecourse=zeros(length(handles.drgbchoices.which_electrodes),sum(these_per_corr),length(t));
                these_all_log_P_timecourse=all_log_P_timecourse(:,these_per_corr,:);
                
                test_out_per_timepoint=zeros(2,per_ii,length(t));
                shuffled_out_per_timepoint=zeros(2,per_ii,length(t));
                
                if per_ii>=20
                    gcp
                    
                    fprintf(1, ['Perceptron processed for mouse No %d ' handles.drgbchoices.per_lab{percent_correct_ii} ' with %d tirals\n'],mouseNo,per_ii);
                    
                    parfor time_point=1:length(t)
                        %         for time_point=1:length(t)
                        
                        per_input=zeros(length(handles.drgbchoices.which_electrodes),per_ii);
                        per_input(:,:)=these_all_log_P_timecourse(:,:,time_point);
                        fprintf(1, '\nTime point %d: ',time_point);
                        
                        %Perceptron
                        %leave one out
                        test_out=zeros(2,per_ii);
                        shuffled_out=zeros(2,per_ii);
                        
                        %Do perceptron analysis for each trial
                        for ii=1:per_ii
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
                                if ii==per_ii
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
                            
                            
                            shuffled_trials=ceil(per_ii*rand(1,num_traces));
                            
                            one_shuffled=zeros(num_traces,1);
                            for jj=1:num_traces
                                one_shuffled(jj,1)=per_input(jj,shuffled_trials(jj));
                            end
                            
                            shuffled_out(:,ii) = net(one_shuffled);
                            
                        end
                        test_out_per_timepoint(:,:,time_point)=test_out;
                        shuffled_out_per_timepoint(:,:,time_point)=shuffled_out;
                        perceptron_correct(percent_correct_ii,time_point)=100*sum(per_targets(1,:)==test_out(1,:))/per_ii;
                        perceptron_correct_shuffled(percent_correct_ii,time_point)=100*sum(per_targets(1,:)==shuffled_out(1,:))/per_ii;
                        fprintf(1, '\nPerceptron percent correct classification %d (for timepoint %d out of %d)\n',100*sum(per_targets(1,:)==test_out(1,:))/per_ii,time_point,length(t));
                    end
                    
                    
                    figNo=figNo+1
                    try
                        close(figNo)
                    catch
                    end
                    
                    figure(figNo)
                    
                    hold on
                    
                    per95=prctile(perceptron_correct_shuffled(percent_correct_ii,:),95);
                    per5=prctile(perceptron_correct_shuffled(percent_correct_ii,:),5);
                    CIsh=[mean(perceptron_correct_shuffled(percent_correct_ii,:))-per5 per95-mean(perceptron_correct_shuffled(percent_correct_ii,:))]';
                    [hlCR, hpCR] = boundedline([t(1) t(end)],[mean(perceptron_correct_shuffled(percent_correct_ii,:)) mean(perceptron_correct_shuffled(percent_correct_ii,:))], CIsh', 'r');
                    
                    plot(t',perceptron_correct(percent_correct_ii,:),'-k')
                    
                    %Odor on markers
                    plot([0 0],[0 100],'-k')
                    odorhl=plot([0 2.5],[20 20],'-k','LineWidth',5);
                    plot([2.5 2.5],[0 100],'-k')
                    
                    title(['Percent correct prediction by perceptron for mouse No ' num2str(mouseNo) ' ' handles.drgbchoices.per_lab{percent_correct_ii}])
                    xlabel('Time (sec)')
                    ylabel('Percent correct')
                    
                    handles_out.perceptron_per_mouse(handles_out.no_mice_with_data).percent_correct(percent_correct_ii).perceptron_correct=zeros(1,length(t));
                    handles_out.perceptron_per_mouse(handles_out.no_mice_with_data).percent_correct(percent_correct_ii).perceptron_correct(1,:)=perceptron_correct(percent_correct_ii,:);
                    
                    handles_out.perceptron_per_mouse(handles_out.no_mice_with_data).percent_correct(percent_correct_ii).perceptron_correct_shuffled=zeros(1,length(t));
                    handles_out.perceptron_per_mouse(handles_out.no_mice_with_data).percent_correct(percent_correct_ii).perceptron_correct_shuffled(1,:)=perceptron_correct_shuffled(percent_correct_ii,:);
                    handles_out.perceptron_per_mouse(handles_out.no_mice_with_data).percent_correct(percent_correct_ii).test_out_per_timepoint=test_out_per_timepoint;
                    handles_out.perceptron_per_mouse(handles_out.no_mice_with_data).percent_correct(percent_correct_ii).shuffled_out_per_timepoint=shuffled_out_per_timepoint;
                    
                    handles_out.t=t';
                    handles_out.perceptron_per_mouse(handles_out.no_mice_with_data).percent_correct(percent_correct_ii).no_trials=per_ii;
                    handles_out.perceptron_per_mouse(handles_out.no_mice_with_data).percent_correct(percent_correct_ii).per_targets=per_targets;
                    
                    these_all_which_events=[];
                    these_all_which_events=all_which_events(:,these_per_corr);
                    handles_out.perceptron_per_mouse(handles_out.no_mice_with_data).percent_correct(percent_correct_ii).which_events=these_all_which_events;
                    
                    these_all_stamped_lick_times=[];
                    these_all_stamped_lick_times=all_stamped_lick_times(these_per_corr,:);
                    handles_out.perceptron_per_mouse(handles_out.no_mice_with_data).percent_correct(percent_correct_ii).stamped_lick_times=these_all_stamped_lick_times;
                    
                    these_all_stamped_lick_ii=[];
                    these_all_stamped_lick_ii=all_stamped_lick_ii(1,these_per_corr);
                    handles_out.perceptron_per_mouse(handles_out.no_mice_with_data).percent_correct(percent_correct_ii).stamped_lick_ii=these_all_stamped_lick_ii;
                    
                else
                    fprintf(1, ['Perceptron not processed for mouse No %d ' handles.drgbchoices.per_lab{percent_correct_ii} ' because there were only %d trials (fewer than 20 trials)\n'],mouseNo,per_ii);
                end
                
                
                pffft=1;
            end
        end
        save([handles.drgb.outPathName 'LFP_per_' handles.drgb.outFileName],'handles_out','-v7.3')
    end
end


 
pffft1=1;


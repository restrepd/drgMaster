function drgLFPspectPerceptron(handles)
%Perceptron analysis

% for ii=1:5
%     try
%         close(ii)
%     catch
%     end
% end

handles.displayData=0;
%Int he future have a dialogue so that the user can choose which events to
%process
evTypeNo=handles.evTypeNo;
% burstLowF=handles.burstLowF;
% burstHighF=handles.burstHighF;
peakLFPNo=handles.peakLFPNo;


handles.drgbchoices.evTypeNos=[5 11];
%S+ and S-
% handles.burstLowF=1;
% handles.burstHighF=100;

handles.drgbchoices.referenceEvent=2;
handles.drgbchoices.timeStart=0.5;
handles.drgbchoices.timeEnd=2.5;


handles.evTypeNo=handles.drgbchoices.referenceEvent; %Process all OdorOn trials

for this_peakLFPNo=1:16
    handles.peakLFPNo=this_peakLFPNo;
    %[t,freq,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event]=drgGetLFPwavePowerForThisEvTypeNo(handles);
    [t,f,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event]=drgGetLFPPowerForThisEvTypeNo(handles);
    freq=f';
    
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
       all_log_P_timecourse=zeros(16,length(this_trialNo),length(t)); 
    end
    all_log_P_timecourse(this_peakLFPNo,:,:)=log_P_timecourse;
end

 
%Do perceptron prediction analysis for every point in the timecourse

num_traces=16;
num_trials=length(this_trialNo);
per_targets=zeros(2,length(this_trialNo));
per_targets(:,:)=which_event;
per_ii=length(this_trialNo);
no_traces=16;

gcp

parfor time_point=1:length(t)
% for time_point=1:length(t)
    per_input=zeros(16,num_trials);
    per_input(:,:)=all_log_P_timecourse(:,:,time_point);
    fprintf(1, '\nTime point %d: ',time_point);
    
    %Perceptron
    %leave one out
    test_out=[];
    shuffled_out=[];
    
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

        
        shuffled_trials=ceil(per_ii*rand(1,no_traces));
        
        one_shuffled=zeros(num_traces,1);
        for jj=1:num_traces
           one_shuffled(jj,1)=per_input(jj,shuffled_trials(jj)); 
        end
        
        shuffled_out(:,ii) = net(one_shuffled);
        
    end
    perceptron_correct(time_point)=100*sum(per_targets(1,:)==test_out(1,:))/per_ii;
    perceptron_correct_shuffled(time_point)=100*sum(per_targets(1,:)==shuffled_out(1,:))/per_ii;
    fprintf(1, '\nPercent correct classification by perceptron = %d, samples processed %d of %d\n',100*sum(per_targets(1,:)==test_out(1,:))/per_ii,time_point,length(t));
end

%save([CaimAn_name(1:end-12) '_pre.mat'])
% save([CaimAn_name(1:end-12) '_perceptron.mat'])


figure(5)

hold on
%plot(time_to_event',perceptron_correct_shuffled,'-b')

per95=prctile(perceptron_correct_shuffled,95);
per5=prctile(perceptron_correct_shuffled,5);
CIsh=[mean(perceptron_correct_shuffled)-per5 per95-mean(perceptron_correct_shuffled)]';
[hlCR, hpCR] = boundedline([time_to_event(1) time_to_event(end)],[mean(perceptron_correct_shuffled) mean(perceptron_correct_shuffled)], CIsh', 'r');

plot(time_to_event',perceptron_correct,'-k')

%Odor on markers
plot([0 0],[0 100],'-k')
odorhl=plot([0 mean(delta_odor)],[20 20],'-k','LineWidth',5);
plot([mean(delta_odor) mean(delta_odor)],[0 100],'-k')

%Reinforcement markers
plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 100],'-r')
reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[20 20],'-r','LineWidth',5);
plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 100],'-r')

title("Percent correct prediction by perceptron")
xlabel('Time (sec)')
ylabel('Percent correct')


pffft1=1;


%Timecourse doing average after log
%Get max and min

if handles.subtractRef==0
    log_P_timecourse=zeros(length(freq),length(t));
    log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),1);
    
    %Per trial power plot
    log_P_per_trial_timecourse=zeros(length(freq)*length(this_trialNo),length(t));
    y_shift=0;
    for trialNo=1:length(this_trialNo)
        this_log_P_timecourse=zeros(length(freq),length(t));
        this_log_P_timecourse(:,:)=10*log10(all_Power_timecourse(trialNo,:,:));
        log_P_per_trial_timecourse(y_shift+1:y_shift+length(freq),:)=this_log_P_timecourse;
        this_power=zeros(1,1,length(freq));
        this_power(1,1,:)=mean(this_log_P_timecourse(:,(t>=handles.drgbchoices.timeStart)&(t>=handles.drgbchoices.timeEnd)),2);
        avg_power_in_window(this_peakLFPNo,trialNo,1:length(freq))=this_power;
        y_shift=y_shift+length(freq);
    end
    
    if handles.autoscale==1
        maxLogPper=prctile(log_P_timecourse(:),99);
        minLogPper=prctile(log_P_timecourse(:),1);
        %Note: Diego added this on purpose to limit the range to 10 dB
        %This results in emphasizing changes in the top 10 dB
        if maxLogP-minLogP>12
            minLogPper=maxLogP-12;
        end
    else
        maxLogPper=handles.maxLogP;
        minLogPper=handles.minLogP;
    end
else
    log_P_timecourse=zeros(length(freq),length(t));
    log_P_timecourse(:,:)=mean(10*log10(all_Power_timecourse),1);
    log_P_timecourse_ref=zeros(length(freq),length(t));
    log_P_timecourse_ref(:,:)=repmat(mean(10*log10(all_Power_ref),1)',1,length(t));
    
    %Per trial power plot and get the power for the time chosen by the user
    log_P_per_trial_timecourse_sub=zeros(length(freq)*length(this_trialNo),length(t));
    y_shift=0;
    sy_shift=0;
    shifted_freq=[];
    for trialNo=1:length(this_trialNo)
        this_log_P_timecourse=zeros(length(freq),length(t));
        this_log_P_timecourse(:,:)=10*log10(all_Power_timecourse(trialNo,:,:));
        this_log_P_timecourse_ref=zeros(length(freq),length(t));
        this_log_P_timecourse_ref(:,:)=repmat(mean(10*log10(all_Power_ref(trialNo,:)),1)',1,length(t));
        log_P_per_trial_timecourse_sub(y_shift+1:y_shift+length(freq),:)=this_log_P_timecourse-this_log_P_timecourse_ref;
        this_power=zeros(1,1,length(freq));
        this_power(1,1,:)=mean(this_log_P_timecourse(:,(t>=handles.drgbchoices.timeStart)&(t>=handles.drgbchoices.timeEnd))...
            -this_log_P_timecourse_ref(:,(t>=handles.drgbchoices.timeStart)&(t>=handles.drgbchoices.timeEnd)),2);
        avg_power_in_window(this_peakLFPNo,trialNo,1:length(freq))=this_power;
        shifted_freq(1,y_shift+1:y_shift+length(freq))=freq+(trialNo-1)*freq(end);
        y_shift=y_shift+length(freq);
    end
    
    max_delta=16;
    if handles.autoscale==1
        maxLogP=prctile(log_P_timecourse(:)-log_P_timecourse_ref(:),99);
        minLogP=prctile(log_P_timecourse(:)-log_P_timecourse_ref(:),1);
        %Note: Diego added this on purpose to limit the range to 10 dB
        %This results in emphasizing changes in the top 10 dB
        if maxLogP-minLogP>max_delta
            minLogP=maxLogP-max_delta;
        end
        
        maxLogPper=prctile(log_P_per_trial_timecourse_sub(:),99);
        minLogPper=prctile(log_P_per_trial_timecourse_sub(:),1);
        %Note: Diego added this on purpose to limit the range to 10 dB
        %This results in emphasizing changes in the top 10 dB
        if maxLogPper-minLogPper>max_delta
            minLogPper=maxLogPper-max_delta;
        end
        
    else
        maxLogPper=handles.maxLogP;
        minLogPper=handles.minLogP;
    end
end

%Calculate the licks
%[lick_freq,times_lick_freq,lick_traces,CIlickf,lick_trace_times,stamped_lick_ii,these_stamped_lick_times,no_trials,trials_included]=drgGetLicks(handles);


if ~isempty(this_trialNo)
    
    
    for bwii=1:4
        figure(bwii+1)
        subplot(4,4,this_peakLFPNo)
        hold on
        
        for trialNo=1:length(this_trialNo)
            %Hit
            if which_event(1,trialNo)==1
                plot(trialNo,mean(avg_power_in_window(this_peakLFPNo,trialNo,(freq>=low_freq(bwii))&(freq<=high_freq(bwii))),3),'or')
            end
            %Miss
            if which_event(2,trialNo)==1
                plot(trialNo,mean(avg_power_in_window(this_peakLFPNo,trialNo,(freq>=low_freq(bwii))&(freq<=high_freq(bwii))),3),'oc')
            end
            %CR
            if which_event(3,trialNo)==1
                plot(trialNo,mean(avg_power_in_window(this_peakLFPNo,trialNo,(freq>=low_freq(bwii))&(freq<=high_freq(bwii))),3),'ob')
            end
            %FA
            if which_event(4,trialNo)==1
                plot(trialNo,mean(avg_power_in_window(this_peakLFPNo,trialNo,(freq>=low_freq(bwii))&(freq<=high_freq(bwii))),3),'om')
            end
        end
        
        %             xlabel('trial No')
        ylabel(num2str(this_peakLFPNo))
        
    end
    
    
    
end


for bwii=1:4
    
    
    figure(bwii+1);
    
    suptitle([freq_names{bwii} ' LFP delta power in dB for ' handles.jtFileName(10:end-4)])
end

%Plot the percent correct
try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.07 0.1 .7 .2])



hold on
trials=1:length(perCorr_pertr);

%Plot in different colors
plot(trials((perCorr_pertr>65)&(perCorr_pertr<80)),perCorr_pertr((perCorr_pertr>65)&(perCorr_pertr<80)),'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
hold on
plot(trials(perCorr_pertr>=80),perCorr_pertr(perCorr_pertr>=80),'or')
plot(trials(perCorr_pertr<=65),perCorr_pertr(perCorr_pertr<=65),'ob')
title('Percent correct per trial')
xlabel('trial No')
ylabel('Percent correct')


handles.evTypeNo=evTypeNo;
handles.burstLowF=burstLowF;
handles.burstHighF=burstHighF;
handles.peakLFPNo=peakLFPNo;
handles.displayData=1;

pffft=1




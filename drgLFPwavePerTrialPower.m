function drgLFPwavePerTrialPower(handles)
%Generates a timecourse of the LFP power in decibels 10*log10(Power)

for ii=1:5
    try
        close(ii)
    catch
    end
end

handles.displayData=0;
%Int he future have a dialogue so that the user can choose which events to
%process
evTypeNo=handles.evTypeNo;
burstLowF=handles.burstLowF;
burstHighF=handles.burstHighF;
peakLFPNo=handles.peakLFPNo;

handles.drgbchoices.evTypeNos=[3 7 9 13];
%Hit, Miss, CR, FA
handles.burstLowF=1;
handles.burstHighF=100;

handles.drgbchoices.referenceEvent=2;
handles.drgbchoices.timeStart=0.5;
handles.drgbchoices.timeEnd=2.5;

low_freq=[6 15 35 65];
high_freq=[12 30 55 95];
freq_names={'Theta','Beta','Low gamma','High gamma'};

handles.evTypeNo=handles.drgbchoices.referenceEvent; %Process all OdorOn trials

for bwii=1:4
    try
        close bwii+1
    catch
    end
    
    hFig = figure(bwii+1);
    set(hFig, 'units','normalized','position',[.07 0.1 .7 .7])
    
end

for this_peakLFPNo=1:16
    handles.peakLFPNo=this_peakLFPNo;
    [t,freq,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo, perCorr_pertr, which_event]=drgGetLFPwavePowerForThisEvTypeNo(handles);
    
    if this_peakLFPNo==1
        avg_power_in_window=zeros(16,length(this_trialNo),length(freq));
    end
    
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
    pffft=1
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




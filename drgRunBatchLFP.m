function drgRunBatchLFP

%Ask user for the m file that contains information on what the user wants the analysis to be
%This file has all the information on what the user wants done, which files
%to process, what groups they fall into, etc
%
% An example of this file: drgbChoicesDanielPrelim
%
%


[choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoices*.m'},'Select the .m file with all the choices for analysis');
addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

new_no_files=handles.drgbchoices.no_files;
choicePathName=handles.drgbchoices.PathName;
choiceFileName=handles.drgbchoices.FileName;

%Very, very important!
handles.evTypeNo=handles.drgbchoices.referenceEvent;

%If you want to skip files that have already been processed enter the number of the first file
first_file=handles.drgb.first_file;

if first_file==1
    handles.drgb.lfpevpair_no=0;
    handles.drgb.lfp_per_exp_no=0;
else
    load([handles.drgb.outPathName handles.drgb.outFileName])
    handles.drgb=handles_drgb.drgb;
    %The user may add new files
    handles.drgbchoices.no_files=new_no_files;
    handles.drgbchoices.PathName=choicePathName;
    handles.drgbchoices.FileName=choiceFileName;
end

test_batch=handles.drgbchoices.test_batch;

first_out=1;

%Do batch processing for each file
for filNum=first_file:handles.drgbchoices.no_files
    
    file_no=filNum
    
    if test_batch==1
        if handles.drgbchoices.group_no(filNum)==1
            handles.data_vs_simulate=5;
        else
            handles.data_vs_simulate=6;
        end
    end
    
    %read the jt_times files
    handles.jtfullName=[handles.drgbchoices.PathName{filNum},handles.drgbchoices.FileName{filNum}];
    handles.jtFileName=handles.drgbchoices.FileName{filNum};
    handles.jtPathName=handles.drgbchoices.PathName{filNum};
    
    drgRead_jt_times(handles.jtPathName,handles.jtFileName);
    FileName=[handles.jtFileName(10:end-4) '_drg.mat'];
    handles.fullName=[handles.jtPathName,FileName];
    handles.FileName=FileName;
    handles.PathName=handles.jtPathName;
    
    load(handles.fullName);
    handles.drg=drg;
    
    if handles.read_entire_file==1
        handles=drgReadAllDraOrDg(handles);
    end
    
    switch handles.drg.session(handles.sessionNo).draq_p.dgordra
        case 1
        case 2
            handles.drg.drta_p.fullName=[handles.jtPathName handles.jtFileName(10:end-4) '.dg'];
        case 3
            handles.drg.drta_p.fullName=[handles.jtPathName handles.jtFileName(10:end-4) '.rhd'];
    end
    
    
    
    %Set the last trial to the last trial in the session
    handles.lastTrialNo=handles.drg.session(handles.sessionNo).events(2).noTimes;
    
    %Save information for this file
    handles.drgb.filNum=filNum;
    handles.drgb.file(filNum).FileName=handles.FileName;
    handles.drgb.file(filNum).PathName=handles.PathName;
    
    
    handles.drgb.file(filNum).drg=handles.drg;
    
    %Run the analysis for each window
    for winNo=1:handles.drgbchoices.noWindows
        
        
        handles.drgb.lfp_per_exp_no=handles.drgb.lfp_per_exp_no+1;
        
        handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).fileNo=filNum;
        handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).referenceEvent=handles.drgbchoices.referenceEvent;
        handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).timeWindow=winNo;
        handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).timeStart=handles.drgbchoices.timeStart(winNo);
        handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).timeEnd=handles.drgbchoices.timeEnd(winNo);
        
        
        if sum(handles.drgbchoices.analyses==2)>0
            handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).allPower=[];
            handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).all_Power_ref=[];
            handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).which_eventLFPPower=[];
            handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).perCorrLFPPower=[];
        end
        
        if sum(handles.drgbchoices.analyses==1)>0
            for ii=1:handles.no_PACpeaks
                handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).no_trials=0;
                handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).meanVectorLength=[];
                handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).meanVectorAngle=[];
                handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).peakAngle=[];
                handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).mod_indx=[];
                handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).all_phase_histo=[];
                handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).perCorrPAC=[];
                handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).which_eventPAC=[];
            end
        end
        
        %Now run the analysis for each lfp
        for lfpNo=1:handles.drgbchoices.no_LFP_elect
            
            
            handles.peakLFPNo=lfpNo;
            handles.burstLFPNo=lfpNo;
            handles.referenceEvent=handles.drgbchoices.referenceEvent;
            
            %Run the analysis for each time window
            
            %First prepare the spectrogram
            %Get LFP power per trial
            if sum(handles.drgbchoices.analyses==2)>0
                
                %Subtracting and adding handles.window/2 gives the correct
                %times
                handles.time_start=min(handles.drgbchoices.timeStart-handles.time_pad-handles.window/2);
                handles.time_end=max(handles.drgbchoices.timeEnd+handles.time_pad+handles.window/2); %2.2 for Shane
                
                handles.burstLowF=handles.LFPPowerSpectrumLowF;
                handles.burstHighF=handles.LFPPowerSpectrumHighF;
                handles.lastTrialNo=handles.drg.session(handles.sessionNo).events(handles.referenceEvent).noTimes;
                handles.trialNo=1;
                all_Power=[];
                all_Power_timecourse=[];
                all_Power_ref=[];
                perCorr=[];
                which_event=[];
                
                %Please note this is the same function called in drgMaster
                %when you choose LFP Power Timecourse Trial Range, which calls drgLFPspectTimecourse
                [t,f,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo, perCorr,which_event]=drgGetLFPPowerForThisEvTypeNo(handles);
                
                if (filNum==1)&(lfpNo==1)
                    handles.drgb.freq_for_LFPpower=f;
                end
            end
            
            fprintf(1, 'File number: %d, window number: %d, lfp number: %d\n',filNum,winNo,lfpNo);
            
            handles.time_start=handles.drgbchoices.timeStart(winNo)-handles.time_pad;
            handles.time_end=handles.drgbchoices.timeEnd(winNo)+handles.time_pad; %2.2 for Shane
            
            
            handles.drgb.lfpevpair_no=handles.drgb.lfpevpair_no+1;
            
            if filNum==3
               pffft=1; 
            end
            
            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).fileNo=filNum;
            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).elecNo=lfpNo;
            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).referenceEvent=handles.drgbchoices.referenceEvent;
            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).timeWindow=winNo;
            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).timeStart=handles.drgbchoices.timeStart(winNo);
            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).timeEnd=handles.drgbchoices.timeEnd(winNo);
            
            
            
            %Get PAC, percent correct and percent lick per trial
            if sum(handles.drgbchoices.analyses==1)>0
                for ii=1:handles.no_PACpeaks
                    handles.n_peak=ii;
                    handles.peakLowF=handles.PACpeakLowF;
                    handles.peakHighF=handles.PACpeakHighF;
                    handles.burstLowF=handles.PACburstLowF(ii);
                    handles.burstHighF=handles.PACburstHighF(ii);
                    
                    
                    %Please note this is the same function called by
                    %drgMaster when the user chooses Phase Amplitude
                    %Coupling
                    handles=drgThetaAmpPhaseTrialRange(handles);

                    
                    %Enter the per LFP values
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).PAC(ii).no_trials=handles.drgb.PAC.no_trials;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).PAC(ii).meanVectorLength=handles.drgb.PAC.meanVectorLength;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).PAC(ii).meanVectorAngle=handles.drgb.PAC.meanVectorAngle;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).PAC(ii).peakAngle=handles.drgb.PAC.peakAngle;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).PAC(ii).mod_indx=handles.drgb.PAC.mod_indx;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).PAC(ii).all_phase_histo=handles.drgb.PAC.all_phase_histo;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).PAC(ii).perCorrPAC=handles.drgb.PAC.perCorr;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).PAC(ii).which_eventPAC=handles.drgb.PAC.which_event;
                    
                   
                    if isempty(handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).PAC(ii).meanVectorLength)
                       pffft=1; 
                    end
                    
                    %Enter the per experiment values
                    handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).no_trials=handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).no_trials+handles.drgb.PAC.no_trials;
                    handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).meanVectorLength=[handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).meanVectorLength handles.drgb.PAC.meanVectorLength];
                    handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).meanVectorAngle=[handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).meanVectorAngle handles.drgb.PAC.meanVectorAngle];
                    handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).peakAngle=[handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).peakAngle handles.drgb.PAC.peakAngle];
                    handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).mod_indx=[handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).mod_indx handles.drgb.PAC.mod_indx];
                    handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).all_phase_histo=[handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).all_phase_histo handles.drgb.PAC.all_phase_histo'];
                    handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).perCorrPAC=[handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).perCorrPAC handles.drgb.PAC.perCorr];
                    handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).which_eventPAC=[handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).PAC(ii).which_eventPAC handles.drgb.PAC.which_event];
                    
                    
                    %handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).percent_lick=handles.drgb.PAC.percent_lick;
                end
            end
            
            %Get LFP power per trial
            if sum(handles.drgbchoices.analyses==2)>0
                sz_all_power=size(all_Power);
                this_all_power=zeros(sz_all_power(1),sz_all_power(2));
                this_all_power(:,:)=mean(all_Power_timecourse(:,:,(t>=handles.drgbchoices.timeStart(winNo))&(t<=handles.drgbchoices.timeEnd(winNo))),3);
                
                %Enter the per LFP values
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).allPower=this_all_power;
                %handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).all_Power_ref=all_Power_ref;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).which_eventLFPPower=which_event;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).perCorrLFPPower=perCorr;
                
                %Enter the per experiment values
                handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).allPower=[handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).allPower this_all_power'];
                %handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).all_Power_ref=[handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).all_Power_ref all_Power_ref'];
                handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).which_eventLFPPower=[handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).which_eventLFPPower which_event];
                handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).perCorrLFPPower=[handles.drgb.lfp_per_exp(handles.drgb.lfp_per_exp_no).perCorrLFPPower perCorr];
                
 
            end
            
            %Do the event-related LFP analysis
            if sum(handles.drgbchoices.analyses==3)>0
                %This was written to answer a reviewer's question on
                %lick-related theta LFP. Because of this I am defaulting
                %the bandwidth of the phase to theta
                handles.peakLowF=6;
                handles.peakHighF=12;
                handles.burstLowF=6;
                handles.burstHighF=12;
                
                handles.peakLFPNo=19; %These are licks
                [log_P_t,no_trials_w_event,which_event,f,out_times,times,phase_per_trial,no_trials,no_events_per_trial,t_per_event_per_trial]=drgEventRelatedAnalysis(handles);
                %[log_P_t,no_trials_w_event,which_event,f,out_times,times,ERLFP_per_trial,phase_per_trial]=drgEventRelatedAnalysis(handles);
                
                if first_out==1
                    handles.drgb.lfpevpair.out_times=out_times;
                    first_out=0;
                end
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).no_trials_w_eventERP=no_trials_w_event;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).which_eventERP=which_event;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).fERP=f;
%                 lpt=zeros(no_trials_w_event,length(f));
%                 lpt(:,:)=log_P_t(:,:,floor(length(out_times)/2)+1+2);
%                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).log_P_tERP=lpt;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).log_P_tERP=log_P_t;
                %handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).timesERP=times; 
                %handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).ERLFP_per_trialERP=ERLFP_per_trial; 
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).phase_per_trialERP=phase_per_trial;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).no_trials=no_trials;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).no_events_per_trial=no_events_per_trial;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no).t_per_event_per_trial=t_per_event_per_trial;
                
            end
        end
        
    end
    
    %Save output file
    handles_drgb=handles;
    if isfield(handles,'data_dg')
        handles_drgb=rmfield(handles_drgb,'data_dg');
    end
    save([handles.drgb.outPathName handles.drgb.outFileName],'handles_drgb','-v7.3')
    
end

pffft=1;





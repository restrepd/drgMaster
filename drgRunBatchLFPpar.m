function drgRunBatchLFPpar

%drgRunBatchLFPpar performs batch analysis of LFP power through fourier
%transform, phase reference power using wavelet analyisis (PRP) and
%phase amplitude coupling (PAC). drgRunBatchLFPpar asks the user for a drgbChoices .m file 
%that has all the information on what the user wants done, which files
%to process, what groups they fall into, etc. The processed data are saved in a .mat file and the analysis can be
%displayed with drgAnalysisBatchLFP
%
%
% handles.drgbchoices.analyses chooses the analysis performed
%
% 1 PAC for case 19 (MI, angle, etc) in drgAnalysisBatchLFP, vetted
% This was used for Figure 2 of Lossaco, Ramirez-Gordillo et al eLife 2020
%
% 2 LFP power for case 20 in drgAnalysisBatchLFP, vetted
%
% 3 Lick-related LFP analysis (fourier transform LFP power)
%
% 4 power for LFP wavelet analysis
%
% 5 event-related LFP wavelet analysis for case 22 of drgAnalysisBatchLFP,
% vetted
%
% 6 phase comparison analysis
%
% 7 phase comparison analysis
%
% 8 Phase-referenced power: LFP wavelet power computed at the peak and through of the low frequency wave
% phase is calculated with a Hilbert transform (PAC), vetted
% This was used for Figure 3 of Lossaco, Ramirez-Gordillo et al eLife 2020

tic

first_file=1;

[choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoices*.m'},'Select the .m file with all the choices for analysis');
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
lfp_per_file=[];
all_files_present=1;
for filNum=first_file:handles.drgbchoices.no_files
    lfp_per_file(filNum).lfp_per_exp_no=0;
    lfp_per_file(filNum).lfpevpair_no=0;
    lfp_per_file(filNum).lfp_per_exp=[];
    lfp_per_file(filNum).lfpevpair=[];
    lfp_per_file(filNum).f=[];
    lfp_per_file(filNum).out_times=[];
    lfp_per_file(filNum).wave_f=[];
    lfp_per_file(filNum).wave_out_times=[];
    lfp_per_file(filNum).eventlabels=[];
     
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
    this_jt=handles.drgbchoices.FileName{filNum};
    handles.temp_exist(filNum)=exist([handles.drgb.outPathName tempDirName '/temp_' this_jt(10:end)]);
    
    if handles.temp_exist(filNum)==2
        %If it was processed load the temp result
        load([handles.drgb.outPathName tempDirName '/temp_' this_jt(10:end)])
        lfp_per_file(filNum)=this_lfp_per_file; 

    end
end


if all_files_present==1
    
    gcp;
    
    no_files=handles.drgbchoices.no_files;
     
    parfor filNum=first_file:no_files
%             for filNum=first_file:no_files
        %         try
        
        file_no=filNum
        handlespf=struct();
            handlespf=handles;
            
            this_jt=handlespf.drgbchoices.FileName{filNum};
            
            %Was the file already processed?
            if handlespf.temp_exist(filNum)==2
                %If it was processed load the temp result
                %If it was processed load the temp result
                %             load([handlespf.drgb.outPathName 'temp/temp_' this_jt(10:end)])
                %             lfp_per_file(filNum)=this_lfp_per_file;
                
                fprintf(1, 'Loaded data from temp file for file number: %d\n',filNum);
            else
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
                handlespf.drg=S.drg;
                lfp_per_file(filNum).eventlabels=handlespf.drg.draq_d.eventlabels;
                
                %         if handles.read_entire_file==1
                %             handlespf=drgReadAllDraOrDg(handlespf);
                %         end
                
                switch handlespf.drg.session(handlespf.sessionNo).draq_p.dgordra
                    case 1
                    case 2
                        handlespf.drg.drta_p.fullName=[jtPathName jtFileName(10:end-4) '.dg'];
                    case 3
                        handlespf.drg.drta_p.fullName=[jtPathName jtFileName(10:end-4) '.rhd'];
                end
                
                %Set the last trial to the last trial in the session
                %handlespf.lastTrialNo=handlespf.drg.session(handlespf.sessionNo).events(2).noTimes;
                handlespf.lastTrialNo=handlespf.drg.session(1).noTrials;
                
                %Run the analysis for each window
                for winNo=1:handles.drgbchoices.noWindows
                    
                    
                    lfp_per_file(filNum).lfp_per_exp_no=lfp_per_file(filNum).lfp_per_exp_no+1;
                    
                    lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).fileNo=filNum;
                    lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).referenceEvent=handles.drgbchoices.referenceEvent;
                    lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).timeWindow=winNo;
                    lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).timeStart=handles.drgbchoices.timeStart(winNo);
                    lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).timeEnd=handles.drgbchoices.timeEnd(winNo);
                    
                    
                    if (sum(handles.drgbchoices.analyses==2)>0)
                        lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).allPower=[];
                        lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).all_Power_ref=[];
                        lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).which_eventLFPPower=[];
                        lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).perCorrLFPPower=[];
                    end
                    
                    if (sum(handles.drgbchoices.analyses==4)>0)
                        lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).wave_allPower=[];
                        lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).wave_all_Power_ref=[];
                        lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).wave_which_eventLFPPower=[];
                        lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).wave_perCorrLFPPower=[];
                    end
                    
                    
                    
                    %Set the lfp_per_file for the PAC analysis
                    if sum(handles.drgbchoices.analyses==1)>0
                        for ii=1:handles.no_PACpeaks
                            lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).no_trials=0;
                            lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).meanVectorLength=[];
                            lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).meanVectorAngle=[];
                            lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).peakAngle=[];
                            lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).mod_indx=[];
                            lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).all_phase_histo=[];
                            lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).perCorrPAC=[];
                            lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).which_eventPAC=[];
                        end
                    end
                    
                    %Set the lfp_per_file for the analysis of wavelet LFP power at peak and
                    %trough of low frequency wave
                    if sum(handles.drgbchoices.analyses==8)>0
                        for ii=1:handles.no_PACpeaks
                            lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC=[];
                            lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PACwave=[]; 
                        end
                    end
                    
                    %Now run the analysis for each lfp
                    for lfpNo=1:handles.drgbchoices.no_LFP_elect
                        
                        
                        handlespf.peakLFPNo=lfpNo;
                        handlespf.burstLFPNo=lfpNo;
                        handlespf.referenceEvent=handles.drgbchoices.referenceEvent;
                        handlespf.subtractRef=handles.drgbchoices.subtractRef;
                        
                        %Run the analysis for each time window
                        
                        %First prepare the spectrogram
                        %Get LFP power per trial
                        if sum(handles.drgbchoices.analyses==2)>0
                            
                            %Subtracting and adding handles.window/2 gives the correct
                            %times
                            handlespf.time_start=min(handles.drgbchoices.timeStart-handles.time_pad-handles.window/2);
                            handlespf.time_end=max(handles.drgbchoices.timeEnd+handles.time_pad+handles.window/2); %2.2 for Shane
                            
                            handlespf.burstLowF=handlespf.LFPPowerSpectrumLowF;
                            handlespf.burstHighF=handlespf.LFPPowerSpectrumHighF;
                            handlespf.lastTrialNo=handlespf.drg.session(1).noTrials;
                            %handlespf.lastTrialNo=handlespf.drg.session(handlespf.sessionNo).events(handlespf.referenceEvent).noTimes;
                            handlespf.trialNo=1;
                            all_Power=[];
                            all_Power_timecourse=[];
                            all_Power_ref=[];
                            perCorr=[];
                            which_event=[];
                            
                            %Please note this is the same function called in drgMaster
                            %when you choose LFP Power Timecourse Trial Range, which calls drgLFPspectTimecourse
                            [t,f,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo, perCorr,which_event]=drgGetLFPPowerForThisEvTypeNo(handlespf);
                            
                            
                            lfp_per_file(filNum).f=f;
                            
                        end
                        
                        
                        
                        handlespf.time_start=handles.drgbchoices.timeStart(winNo)-handles.time_pad;
                        handlespf.time_end=handles.drgbchoices.timeEnd(winNo)+handles.time_pad; %2.2 for Shane
                        
                        
                        lfp_per_file(filNum).lfpevpair_no=lfp_per_file(filNum).lfpevpair_no+1;
                        lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).fileNo=filNum;
                        lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).elecNo=lfpNo;
                        lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).referenceEvent=handles.drgbchoices.referenceEvent;
                        lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).timeWindow=winNo;
                        lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).timeStart=handles.drgbchoices.timeStart(winNo);
                        lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).timeEnd=handles.drgbchoices.timeEnd(winNo);
                        
                        
                        
                        %Get PAC, percent correct and percent lick per trial
                        if sum(handles.drgbchoices.analyses==1)>0
                            for ii=1:handles.no_PACpeaks
                                handlespf.n_peak=ii;
                                handlespf.peakLowF=handles.PACpeakLowF;
                                handlespf.peakHighF=handles.PACpeakHighF;
                                handlespf.burstLowF=handles.PACburstLowF(ii);
                                handlespf.burstHighF=handles.PACburstHighF(ii);
                                
                                
                                %Please note this is the same function called by
                                %drgMaster when the user chooses Phase Amplitude
                                %Coupling
                                handlespf=drgThetaAmpPhaseTrialRange(handlespf);
                                
                                
                                %Enter the per LFP values
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).no_trials=handlespf.drgb.PAC.no_trials;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).meanVectorLength=handlespf.drgb.PAC.meanVectorLength;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).meanVectorAngle=handlespf.drgb.PAC.meanVectorAngle;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).peakAngle=handlespf.drgb.PAC.peakAngle;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).mod_indx=handlespf.drgb.PAC.mod_indx;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).all_phase_histo=handlespf.drgb.PAC.all_phase_histo;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).perCorrPAC=handlespf.drgb.PAC.perCorr;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).which_eventPAC=handlespf.drgb.PAC.which_event;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).peakAngleForPower=handlespf.drgb.PAC.peakAngleForPower;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).troughAngleForPower=handlespf.drgb.PAC.troughAngleForPower;
                                
                                %Enter the per experiment values
                                lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).no_trials=lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).no_trials+handlespf.drgb.PAC.no_trials;
                                lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).meanVectorLength=[lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).meanVectorLength handlespf.drgb.PAC.meanVectorLength];
                                lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).meanVectorAngle=[lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).meanVectorAngle handlespf.drgb.PAC.meanVectorAngle];
                                lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).peakAngle=[lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).peakAngle handlespf.drgb.PAC.peakAngle];
                                lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).mod_indx=[lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).mod_indx handlespf.drgb.PAC.mod_indx];
                                lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).all_phase_histo=[lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).all_phase_histo handlespf.drgb.PAC.all_phase_histo'];
                                lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).perCorrPAC=[lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).perCorrPAC handlespf.drgb.PAC.perCorr];
                                lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).which_eventPAC=[lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).PAC(ii).which_eventPAC handlespf.drgb.PAC.which_event];
                                
                            end
                        end
                        
                        %Get LFP wavelet power at the trough and peak of low frequency wave
                        if sum(handles.drgbchoices.analyses==8)>0
                            for ii=1:handles.no_PACpeaks
                                handlespf.n_peak=ii;
                                handlespf.peakLowF=handles.PACpeakLowF;
                                handlespf.peakHighF=handles.PACpeakHighF;
                                handlespf.burstLowF=handles.PACburstLowF(ii);
                                handlespf.burstHighF=handles.PACburstHighF(ii);
                                this_subtractRef=handlespf.subtractRef;
                                handlespf.subtractRef=0;
                                
                                %Please note this is the same function called by
                                %drgMaster when the user chooses Phase Amplitude
                                %Coupling
                                
                                %Calculate licks only for the first
                                %electrode
                                if (handlespf.peakLFPNo==1)&(ii==1)
                                    handlespf.calculate_lick=1;
                                else
                                    handlespf.calculate_lick=0;
                                end
                                 
                                handlespf=drgLFPwaveTimecourse(handlespf);
                                
                                handlespf.subtractRef=this_subtractRef;
                                
                                %Enter the per LFP values
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).no_trials=handlespf.drgb.PAC.no_trials;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).meanVectorLength=handlespf.drgb.PAC.meanVectorLength;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).meanVectorAngle=handlespf.drgb.PAC.meanVectorAngle;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).peakAngle=handlespf.drgb.PAC.peakAngle;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).mod_indx=handlespf.drgb.PAC.mod_indx;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).all_phase_histo=handlespf.drgb.PAC.all_phase_histo;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).perCorrPAC=handlespf.drgb.PAC.perCorr;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).which_eventPAC=handlespf.drgb.PAC.which_event;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).peakAngleForPower=handlespf.drgb.PAC.peakAngleForPower;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PAC(ii).troughAngleForPower=handlespf.drgb.PAC.troughAngleForPower;
                                
                                %Now save PACwave
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).t_pac=handlespf.drgb.PACwave.t_pac;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).PACtimecourse=handlespf.drgb.PACwave.PACtimecourse;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).meanPeakPower=handlespf.drgb.PACwave.meanPeakPower;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).meanTroughPower=handlespf.drgb.PACwave.meanTroughPower;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).meanPower=handlespf.drgb.PACwave.meanPower;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).mean_lick_freq=handlespf.drgb.PACwave.mean_lick_freq;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).meanPeakPower_per_lick_trial=handlespf.drgb.PACwave.meanPeakPower_per_lick_trial;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).meanTroughPower_per_lick_trial=handlespf.drgb.PACwave.meanTroughPower_per_lick_trial;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).trNo_lick_and_PRP=handlespf.drgb.PACwave.trNo_lick_and_PRP;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).trialNos_PRP=handlespf.drgb.PACwave.trialNos_PRP;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).trialNos_PAC=handlespf.drgb.PACwave.trialNos_PAC;
                                
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).lick_trials_included=handlespf.drgb.PACwave.lick_trials_included;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).lick_freq_per_trial=handlespf.drgb.PACwave.lick_freq_per_trial;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).lick_freq_per_trial=handlespf.drgb.PACwave.lick_freq_per_trial;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).stamped_lick_ii=handlespf.drgb.PACwave.stamped_lick_ii;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).stamped_lick_ii=handlespf.drgb.PACwave.stamped_lick_ii;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).stamped_lick_times=handlespf.drgb.PACwave.these_stamped_lick_times;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).meanLickPower=handlespf.drgb.PACwave.meanLickPower;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).lick_triggered_wpower=handlespf.drgb.PACwave.lick_triggered_wpower;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).lick_triggered_wpower_t=handlespf.drgb.PACwave.lick_triggered_wpower_t;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).lick_triggered_wpower_no=handlespf.drgb.PACwave.lick_triggered_wpower_no;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).lick_timecourse=handlespf.drgb.PACwave.lick_timecourse;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).no_lickPower_trials=handlespf.drgb.PACwave.no_lickPower_trials;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).lickPower_trials=handlespf.drgb.PACwave.lickPower_trials;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).peakPowerSpectrum=handlespf.drgb.PACwave.peakPowerSpectrum;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).troughPowerSpectrum=handlespf.drgb.PACwave.troughPowerSpectrum;
%                                 lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(ii).which_event_licks=handlespf.drgb.PACwave.which_event_licks;
                                
                                %Save the lick data
                                if (handlespf.calculate_lick==1)&(ii==1)
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).times_lick_freq=handlespf.drgb.PACwave.times_lick_freq;
%                                     lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(1).t_pac=handlespf.drgb.PACwave.t_pac;
%                                     lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).PACwave(1).PACtimecourse=handlespf.drgb.PACwave.PACtimecourse;
%                                     lfp_per_file(filNum).lfpevpair(1).PACwave(1).meanPeakPower=handlespf.drgb.PACwave.meanPeakPower;
%                                     lfp_per_file(filNum).lfpevpair(1).PACwave(1).meanTroughPower=handlespf.drgb.PACwave.meanTroughPower;
                                    %                                 lfp_per_file(filNum).lfpevpair(1).PACwave(1).meanPower=handlespf.drgb.PACwave.meanPower;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).mean_lick_freq=handlespf.drgb.PACwave.mean_lick_freq;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).meanPeakPower_per_lick_trial=handlespf.drgb.PACwave.meanPeakPower_per_lick_trial;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).meanTroughPower_per_lick_trial=handlespf.drgb.PACwave.meanTroughPower_per_lick_trial;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).trNo_lick_and_PRP=handlespf.drgb.PACwave.trNo_lick_and_PRP;
%                                     lfp_per_file(filNum).lfpevpair(1).PACwave(1).trialNos_PRP=handlespf.drgb.PACwave.trialNos_PRP;
%                                     lfp_per_file(filNum).lfpevpair(1).PACwave(1).trialNos_PAC=handlespf.drgb.PACwave.trialNos_PAC;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).times_lick_freq=handlespf.drgb.PACwave.times_lick_freq;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).lick_trials_included=handlespf.drgb.PACwave.lick_trials_included;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).lick_freq_per_trial=handlespf.drgb.PACwave.lick_freq_per_trial;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).lick_freq_per_trial=handlespf.drgb.PACwave.lick_freq_per_trial;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).stamped_lick_ii=handlespf.drgb.PACwave.stamped_lick_ii;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).stamped_lick_ii=handlespf.drgb.PACwave.stamped_lick_ii;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).stamped_lick_times=handlespf.drgb.PACwave.these_stamped_lick_times;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).meanLickPower=handlespf.drgb.PACwave.meanLickPower;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).lick_triggered_wpower=handlespf.drgb.PACwave.lick_triggered_wpower;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).lick_triggered_wpower_t=handlespf.drgb.PACwave.lick_triggered_wpower_t;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).lick_triggered_wpower_no=handlespf.drgb.PACwave.lick_triggered_wpower_no;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).lick_timecourse=handlespf.drgb.PACwave.lick_timecourse;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).no_lickPower_trials=handlespf.drgb.PACwave.no_lickPower_trials;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).lickPower_trials=handlespf.drgb.PACwave.lickPower_trials;
%                                     lfp_per_file(filNum).lfpevpair(1).PACwave(1).peakPowerSpectrum=handlespf.drgb.PACwave.peakPowerSpectrum;
%                                     lfp_per_file(filNum).lfpevpair(1).PACwave(1).troughPowerSpectrum=handlespf.drgb.PACwave.troughPowerSpectrum;
                                    lfp_per_file(filNum).lfpevpair(1).PACwave(1).which_event_licks=handlespf.drgb.PACwave.which_event_licks;
                                    
                                end
                            end
                        end
                        
                        %Get LFP power per trial
                        if sum(handles.drgbchoices.analyses==2)>0
                            sz_all_power=size(all_Power);
                            this_all_power=zeros(sz_all_power(1),sz_all_power(2));
                            this_all_power(:,:)=mean(all_Power_timecourse(:,:,(t>=handles.drgbchoices.timeStart(winNo))&(t<=handles.drgbchoices.timeEnd(winNo))),3);
                            
                            %Enter the per LFP values
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).allPower=this_all_power;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).which_eventLFPPower=which_event;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).perCorrLFPPower=perCorr;
                            
                            %Enter the per experiment values
                            lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).allPower=[lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).allPower this_all_power'];
                            lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).which_eventLFPPower=[lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).which_eventLFPPower which_event];
                            lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).perCorrLFPPower=[lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).perCorrLFPPower perCorr];
                            
                            
                        end
                        
                        %Do the event-related LFP analysis
                        if sum(handles.drgbchoices.analyses==3)>0
                            %This was written to answer a reviewer's question on
                            %lick-related theta LFP.
                            handlespf.peakLowF=6;
                            handlespf.peakHighF=12;
                            handlespf.burstLowF=handlespf.LFPPowerSpectrumLowF;
                            handlespf.burstHighF=handlespf.LFPPowerSpectrumHighF;
                            
                            handlespf.peakLFPNo=19; %These are licks
                            [log_P_t,no_trials_w_event,which_event,f,out_times,times,phase_per_trial,no_trials,no_events_per_trial,t_per_event_per_trial,trial_map,perCorr,no_ref_evs_per_trial]=drgEventRelatedAnalysis(handlespf);
                            
                            lfp_per_file(filNum).out_times=out_times;
                            
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).no_trials_w_eventERP=no_trials_w_event;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).which_eventERP=which_event;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).fERP=f;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).log_P_tERP=log_P_t;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).phase_per_trialERP=phase_per_trial;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).no_trials=no_trials;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).no_events_per_trial=no_events_per_trial;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).t_per_event_per_trial=t_per_event_per_trial;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).trial_map=trial_map;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).perCorrERP=perCorr;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).no_ref_evs_per_trial=no_ref_evs_per_trial;
                            
                        end
                        
                        %Do the wavelet LFP power
                        if sum(handles.drgbchoices.analyses==4)>0
                            %Subtracting and adding handles.window/2 gives the correct
                            %times
                            handlespf.time_start=min(handles.drgbchoices.timeStart-handles.time_pad-handles.window/2);
                            handlespf.time_end=max(handles.drgbchoices.timeEnd+handles.time_pad+handles.window/2); %2.2 for Shane
                            
                            handlespf.burstLowF=handlespf.LFPPowerSpectrumLowF;
                            handlespf.burstHighF=handlespf.LFPPowerSpectrumHighF;
                            handlespf.lastTrialNo=handlespf.drg.session(1).noTrials;
                            %handlespf.lastTrialNo=handlespf.drg.session(handlespf.sessionNo).events(handlespf.referenceEvent).noTimes;
                            handlespf.trialNo=1;
                            all_Power=[];
                            all_Power_timecourse=[];
                            all_Power_ref=[];
                            perCorr=[];
                            which_event=[];
                            out_times=[];
                            f=[];
                            all_Power_timecourse=[];
                            
                            %Please note this is the same function called in drgMaster
                            [t,f,all_Power,all_Power_ref, all_Power_timecourse, this_trialNo, perCorr,which_event]=drgGetLFPwavePowerForThisEvTypeNo(handlespf);
                            
                            lfp_per_file(filNum).wave_f=f;
                            
                            sz_all_power=size(all_Power);
                            this_all_power=zeros(sz_all_power(1),sz_all_power(2));
                            this_all_power(:,:)=mean(all_Power_timecourse(:,:,(t>=handles.drgbchoices.timeStart(winNo))&(t<=handles.drgbchoices.timeEnd(winNo))),3);
                            
                            %Enter the per LFP values
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_allPower=this_all_power;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_refPower=all_Power_ref;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_trialNo=this_trialNo;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_which_eventLFPPower=which_event;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_perCorrLFPPower=perCorr;
                            
                            %                         %Enter the per experiment values
                            %                         lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).wave_allPower=[lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).wave_allPower this_all_power];
                            %                         lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).wave_which_eventLFPPower=[lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).which_eventLFPPower which_event];
                            %                         lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).wave_perCorrLFPPower=[lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).perCorrLFPPower perCorr];
                            
                        end
                        
                        %Do the event-related LFP wavelet analysis
                        if sum(handles.drgbchoices.analyses==5)>0
                            handlespf.peakLowF=6;
                            handlespf.peakHighF=12;
                            handlespf.burstLowF=handlespf.LFPPowerSpectrumLowF;
                            handlespf.burstHighF=handlespf.LFPPowerSpectrumHighF;
                            
                            handlespf.peakLFPNo=19; %These are licks
                            
                            log_P_t=[];
                            which_event=[];
                            freq=[];
                            out_times=[];
                            perCorrERP=[];
                            no_events_per_trial=[];
                            trials_with_wavelick=[];
                            f_lick=[];
                            mean_lick_phase=[];
                            no_trials=[];
                            
                            [log_P_t,which_event,freq,out_times,perCorrERP,no_events_per_trial,trials_with_wavelick,f_lick,mean_lick_phase,no_trials]=drgLFP_ERWA(handlespf);
                            
                            
                            %                             lfp_per_file(filNum).wave_out_times=out_times;
                            
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_out_times=out_times;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_no_trialsERWA=no_trials;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_which_eventERWA=which_event;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_fERWA=freq;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_trialNoERWA=trials_with_wavelick;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_perCorrERWA=perCorrERP;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_f_lick=f_lick;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_mean_lick_phase=mean_lick_phase;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_trials_with_wavelick=trials_with_wavelick;
                            
                            %Lick-referenced data
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_log_P_t_lick_referenced=log_P_t;
                            lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).wave_no_licks_per_trial=no_events_per_trial;
                            
                        end
                        
                        %Do the phase comparison analysis
                        if sum(handles.drgbchoices.analyses==6)>0
                            low_freq=handles.drgbchoices.phaseBandLow;
                            high_freq=handles.drgbchoices.phaseBandHigh;
                            
                            for freq_ii=1:length(low_freq)
                                handlespf.peakLowF=low_freq(freq_ii);
                                handlespf.peakHighF=high_freq(freq_ii);
                                handlespf.burstLowF=low_freq(freq_ii);
                                handlespf.burstHighF=high_freq(freq_ii);
                                
                                elec2=lfpNo;
                                for elec1=1:handles.drgbchoices.no_LFP_elect
                                    
                                    handlespf.peakLFPNo=elec1;
                                    handlespf.burstLFPNo=elec2;
                                    
                                    [rho,which_event,no_trials,angleLFPexp,filtLFPexp,filtLFPref]=drgComparePhases(handlespf);
                                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).no_trials=no_trials;
                                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).which_event=which_event;
                                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).rho=rho;
                                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).angleLFPexp=angleLFPexp;
                                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).filtLFPexp=filtLFPexp;
                                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).filtLFPref=filtLFPref;
                                    
                                end
                            end
                            
                        end
                        
                        %Do the phase comparison analysis
                        if sum(handles.drgbchoices.analyses==7)>0
                            low_freq=handles.drgbchoices.phaseBandLow;
                            high_freq=handles.drgbchoices.phaseBandHigh;
                            
                            for freq_ii=1:length(low_freq)
                                handlespf.peakLowF=low_freq(freq_ii);
                                handlespf.peakHighF=high_freq(freq_ii);
                                handlespf.burstLowF=low_freq(freq_ii);
                                handlespf.burstHighF=high_freq(freq_ii);
                                
                                elect1=handles.drgbchoices.refPhase;
                                elect2=lfpNo;
                                
                                handlespf.peakLFPNo=elect1;
                                handlespf.burstLFPNo=elect2;
                                
                                [rho,which_event,no_trials,angleLFPexp,filtLFPexp,filtLFPref]=drgComparePhases(handlespf);
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elect1,elect2).no_trials=no_trials;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elect1,elect2).which_event=which_event;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elect1,elect2).rho=rho;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elect1,elect2).angleLFPexp=angleLFPexp;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elect1,elect2).filtLFPexp=filtLFPexp;
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elect1,elect2).filtLFPref=filtLFPref;
                                
                            end
                            
                        end
                        
                        
                        
                        
                        fprintf(1, 'File number: %d, window number: %d, lfp number: %d\n',filNum,winNo,lfpNo);
                        
                    end
                    
                end
                
                %Save this temp file
                drgSavePar([handlespf.drgb.outPathName tempDirName '/'],this_jt(10:end),lfp_per_file(filNum),filNum)
                
                
            end
%         catch
%             fprintf(1, '\n\nProcessing failed for file number: %d\n\n',filNum);
%             
%             %Save this failed file
%             drgSaveParFail([handlespf.drgb.outPathName tempDirName '/'],this_jt(10:end),filNum,handlespf)
%             
%         end
    end
    
    
    %Now enter the results in the output structures
    %Initialize lfpevpair
    
    handles.drgb.lfpevpair=[];
    handles.drgb.lfpevpair_no=0;
    
    if (sum(handles.drgbchoices.analyses==3)>0)||(sum(handles.drgbchoices.analyses==5)>0)
        handles.drgb.lfpevpair.out_times=lfp_per_file(1).out_times;
    end
    
    for filNum=first_file:handles.drgbchoices.no_files
        
        found_results=1;
        
        %If there are no results skip this file
        if (sum(handles.drgbchoices.analyses==8)>0)
            if ~isfield(lfp_per_file(filNum).lfpevpair(1),'PAC')
                found_results=0;
            end
        end
        
        if found_results==1
            jtFileName=handles.drgbchoices.FileName{filNum};
            
            if iscell(handles.drgbchoices.PathName)
                jtPathName=handles.drgbchoices.PathName{filNum};
            else
                jtPathName=handles.drgbchoices.PathName;
            end
            
            
            %Save information for this file
            handles.drgb.filNum=filNum;
            handles.drgb.file(filNum).FileName=[jtFileName(10:end-4) '_drg.mat'];
            handles.drgb.file(filNum).PathName=jtPathName;
            
            %Save LFP power structures
            if (sum(handles.drgbchoices.analyses==2)>0)||(sum(handles.drgbchoices.analyses==1)>0)||(sum(handles.drgbchoices.analyses==4)>0)
                handles.drgb.file(filNum).freq_for_LFPpower=lfp_per_file(filNum).f;
            end
            
            handles.drgb.file(filNum).eventlabels=lfp_per_file(filNum).eventlabels;
            
            
            
            for ii=1:lfp_per_file(filNum).lfpevpair_no
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).fileNo=...
                    lfp_per_file(filNum).lfpevpair(ii).fileNo;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).elecNo=...
                    lfp_per_file(filNum).lfpevpair(ii).elecNo;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).referenceEvent=...
                    lfp_per_file(filNum).lfpevpair(ii).referenceEvent;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).timeWindow=...
                    lfp_per_file(filNum).lfpevpair(ii).timeWindow;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).timeStart=...
                    lfp_per_file(filNum).lfpevpair(ii).timeStart;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).timeEnd=...
                    lfp_per_file(filNum).lfpevpair(ii).timeEnd;
                
                %LFP power
                if (sum(handles.drgbchoices.analyses==2)>0)
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).allPower=...
                        lfp_per_file(filNum).lfpevpair(ii).allPower;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).which_eventLFPPower=...
                        lfp_per_file(filNum).lfpevpair(ii).which_eventLFPPower;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).perCorrLFPPower=...
                        lfp_per_file(filNum).lfpevpair(ii).perCorrLFPPower;
                end
                
                %LFP PAC
                if (sum(handles.drgbchoices.analyses==1)>0)
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).PAC=...
                        lfp_per_file(filNum).lfpevpair(ii).PAC;
                end
                
                %LFP wave power at peak and trough of PAC
                if (sum(handles.drgbchoices.analyses==8)>0)
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).PAC=...
                        lfp_per_file(filNum).lfpevpair(ii).PAC;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).PACwave=...
                        lfp_per_file(filNum).lfpevpair(ii).PACwave;
                end
                
                if (sum(handles.drgbchoices.analyses==4)>0)
                    handles.drgb.file(filNum).freq_for_LFPpower=lfp_per_file(filNum).wave_f;
                end
                
                %LFP wavelet power
                if (sum(handles.drgbchoices.analyses==4)>0)
                    %Enter the per LFP values
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_allPower=...
                        lfp_per_file(filNum).lfpevpair(ii).wave_allPower;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_trialNo=...
                        lfp_per_file(filNum).lfpevpair(ii).wave_trialNo;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_which_eventLFPPower=...
                        lfp_per_file(filNum).lfpevpair(ii).wave_which_eventLFPPower;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_perCorrLFPPower=...
                        lfp_per_file(filNum).lfpevpair(ii).wave_perCorrLFPPower;
                end
                
                %              %LFP ERWA
                %              if (sum(handles.drgbchoices.analyses==5)>0)
                %                  handles.drgb.file(filNum).waveERWA_out_times=lfp_per_file(filNum).wave_out_times;
                %              end
                
                if (sum(handles.drgbchoices.analyses==5)>0)
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_out_times=...
                        lfp_per_file(filNum).lfpevpair(ii).wave_out_times;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_no_trialsERWA=...
                        lfp_per_file(filNum).lfpevpair(ii).wave_no_trialsERWA;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_which_eventERWA=...
                        lfp_per_file(filNum).lfpevpair(ii).wave_which_eventERWA;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_fERWA=...
                        lfp_per_file(filNum).lfpevpair(ii).wave_fERWA;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_trialNoERWA=...
                        lfp_per_file(filNum).lfpevpair(ii).wave_trialNoERWA;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_perCorrERWA=...
                        lfp_per_file(filNum).lfpevpair(ii).wave_perCorrERWA;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_f_lick=...
                        lfp_per_file(filNum).lfpevpair(ii).wave_f_lick;
                    
                    
                    %Lick-referenced data
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_log_P_t_lick_referenced=...
                        lfp_per_file(filNum).lfpevpair(ii).wave_log_P_t_lick_referenced;
                    %                  handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_lick_phase_per_trialERP=...
                    %                      lfp_per_file(filNum).lfpevpair(ii).wave_lick_phase_per_trialERP;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_no_licks_per_trial=...
                        lfp_per_file(filNum).lfpevpair(ii).wave_no_licks_per_trial;
                    %                  handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_t_per_lick_per_trial=...
                    %                      lfp_per_file(filNum).lfpevpair(ii).wave_t_per_lick_per_trial;
                    %                  handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_no_ref_licks_per_trial=...
                    %                      lfp_per_file(filNum).lfpevpair(ii).wave_no_ref_licks_per_trial;
                    %                  handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).delta_lick_times=...
                    %                      lfp_per_file(filNum).lfpevpair(ii).delta_lick_times;
                    %                  handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).lick_referenced_LFP_per_trial=...
                    %                      lfp_per_file(filNum).lfpevpair(ii).lick_referenced_LFP_per_trial;
                    
                    %Theta-referenced data
                    %                  handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_log_P_t_ltheta_related=...
                    %                      lfp_per_file(filNum).lfpevpair(ii).wave_log_P_t_ltheta_related;
                    %                  handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_no_thetas_per_trial=...
                    %                      lfp_per_file(filNum).lfpevpair(ii).wave_no_thetas_per_trial;
                    %                  handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_t_per_theta_per_trial=...
                    %                      lfp_per_file(filNum).lfpevpair(ii).wave_t_per_theta_per_trial;
                    % %                  handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_theta_referenced_LFP_per_trial=...
                    % %                      lfp_per_file(filNum).lfpevpair(ii).wave_theta_referenced_LFP_per_trial;
                    %                  handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_no_ref_evs_per_trial_th=...
                    %                      lfp_per_file(filNum).lfpevpair(ii).wave_no_ref_evs_per_trial_th;
                    %                  handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).delta_lick_times_theta_ref=...
                    %                      lfp_per_file(filNum).lfpevpair(ii).delta_lick_times_theta_ref;
                    
                end
                
                %LFP ERP
                if (sum(handles.drgbchoices.analyses==3)>0)
                    
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).no_trials_w_eventERP=...
                        lfp_per_file(filNum).lfpevpair(ii).no_trials_w_eventERP;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).which_eventERP=...
                        lfp_per_file(filNum).lfpevpair(ii).which_eventERP;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).fERP=...
                        lfp_per_file(filNum).lfpevpair(ii).fERP;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).log_P_tERP=...
                        lfp_per_file(filNum).lfpevpair(ii).log_P_tERP;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).phase_per_trialERP=...
                        lfp_per_file(filNum).lfpevpair(ii).phase_per_trialERP;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).no_trials=...
                        lfp_per_file(filNum).lfpevpair(ii).no_trials;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).no_events_per_trial=...
                        lfp_per_file(filNum).lfpevpair(ii).no_events_per_trial;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).t_per_event_per_trial=...
                        lfp_per_file(filNum).lfpevpair(ii).t_per_event_per_trial;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).trial_map=...
                        lfp_per_file(filNum).lfpevpair(ii).trial_map;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).perCorrERP=...
                        lfp_per_file(filNum).lfpevpair(ii).perCorrERP;
                    handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).no_ref_evs_per_trial=...
                        lfp_per_file(filNum).lfpevpair(ii).no_ref_evs_per_trial;
                end
                
                %LFP phase
                if (sum(handles.drgbchoices.analyses==6)>0)
                    
                    for freq_ii=1:length(low_freq)
                        
                        for elec1=1:16
                            for elec2=elec1+1:16
                                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii,freq_ii,elec1,elec2).no_trials=...
                                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).no_trials;
                                
                                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii,freq_ii,elec1,elec2).which_event=...
                                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).which_event;
                                
                                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii,freq_ii,elec1,elec2).rho=...
                                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).rho;
                                
                                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii,freq_ii,elec1,elec2).angleLFPexp=...
                                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).angleLFPexp;
                                
                                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii,freq_ii,elec1,elec2).filtLFPexp=...
                                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).filtLFPexp;
                                
                                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii,freq_ii,elec1,elec2).filtLFPref=...
                                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).filtLFPref;
                            end
                        end
                    end
                    
                    
                end
                
                if (sum(handles.drgbchoices.analyses==7)>0)
                    
                    low_freq=handles.drgbchoices.phaseBandLow;
                    
                    for freq_ii=1:length(low_freq)
                        
                        elec1=handles.drgbchoices.refPhase;
                        for elec2=1:16
                            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii,freq_ii,elec1,elec2).no_trials=...
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).no_trials;
                            
                            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii,freq_ii,elec1,elec2).which_event=...
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).which_event;
                            
                            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii,freq_ii,elec1,elec2).rho=...
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).rho;
                            
                            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii,freq_ii,elec1,elec2).angleLFPexp=...
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).angleLFPexp;
                            
                            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii,freq_ii,elec1,elec2).filtLFPexp=...
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).filtLFPexp;
                            
                            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii,freq_ii,elec1,elec2).filtLFPref=...
                                lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no,freq_ii,elec1,elec2).filtLFPref;
                        end
                    end
                    
                    
                    
                end
                
                
            end
            
            handles.drgb.lfpevpair_no=lfp_per_file(filNum).lfpevpair_no+handles.drgb.lfpevpair_no;
        end
        
    end
    
    
    %Save output file
    handles_drgb=handles;
    if isfield(handles,'data_dg')
        handles_drgb=rmfield(handles_drgb,'data_dg');
    end
    save([handles.drgb.outPathName handles.drgb.outFileName],'handles_drgb','-v7.3')
    
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
    
end







function drgRunBatchLFPxLFPpar

%Ask user for the m file that contains information on what the user wants the analysis to be
%This file has all the information on what the user wants done, which files
%to process, what groups they fall into, etc
%
% An example of this file: drgbChoicesDanielPrelim
%
%

tic

first_file=1;

[choiceFileName,choiceBatchPathName] = uigetfile({'drgbChoices*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgRunBatchLFP run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

new_no_files=handles.drgbchoices.no_files;
choicePathName=handles.drgbchoices.PathName;
choiceFileName=handles.drgbchoices.FileName;
%Very, very important!
handles.evTypeNo=handles.drgbchoices.referenceEvent;

test_batch=handles.drgbchoices.test_batch;

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
    
%     gcp;
    
    no_files=handles.drgbchoices.no_files;
%     parfor filNum=first_file:no_files
        
            for filNum=first_file:handles.drgbchoices.no_files
        
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
            lfp_per_file(filNum).drg=handlespf.drg;
            
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
            handlespf.lastTrialNo=handlespf.drg.session(handlespf.sessionNo).events(2).noTimes;
            
            %Run the analysis for each window
            for winNo=1:handles.drgbchoices.noWindows
                
                
                lfp_per_file(filNum).lfp_per_exp_no=lfp_per_file(filNum).lfp_per_exp_no+1;
                lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).fileNo=filNum;
                lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).referenceEvent=handles.drgbchoices.referenceEvent;
                lfp_per_file(filNum).lfp_per_exp(lfp_per_file(filNum).lfp_per_exp_no).timeWindow=winNo;
                
                handlespf.time_start=handles.drgbchoices.timeStart(winNo)-handles.time_pad;
                handlespf.time_end=handles.drgbchoices.timeEnd(winNo)+handles.time_pad; %2.2 for Shane
                
                
                %Now run the analysis for each lfp
                shuffled_elec=randperm(handles.drgbchoices.no_LFP_elect_pairs);
                
                for lfppairNo=1:handles.drgbchoices.no_LFP_elect_pairs
                    
                    
                    handlespf.peakLFPNo=handles.drgbchoices.which_electrodes1(lfppairNo);
                    handlespf.burstLFPNo=handles.drgbchoices.which_electrodes2(shuffled_elec(lfppairNo));
                    handlespf.referenceEvent=handles.drgbchoices.referenceEvent;
                 
                    
                    lfp_per_file(filNum).lfpevpair_no=lfp_per_file(filNum).lfpevpair_no+1;
                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).fileNo=filNum;
                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).elec_pair_No=lfppairNo;
                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).elec1=handlespf.peakLFPNo;
                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).elec2= handlespf.burstLFPNo;
                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).referenceEvent=handles.drgbchoices.referenceEvent;
                    lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).timeWindow=winNo;
                     
                    %Do the LFP coherence analysis
                    if sum(handles.drgbchoices.analyses==8)>0
                        %This is coherence analysis
                        handlespf.peakLowF=6;
                        handlespf.peakHighF=12;
                        handlespf.burstLowF=handlespf.LFPPowerSpectrumLowF;
                        handlespf.burstHighF=handlespf.LFPPowerSpectrumHighF;
                         
                        [out_times,f, all_Cxy_timecourse, trial_numbers, perCorr_pertr, which_event]=drgGetLFPCoherenceForThisEvTypeNo(handlespf);
                        
                        lfp_per_file(filNum).out_times=out_times;
                        
                        lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).out_times_coh=out_times;
                        lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).no_trials_coh=length(trial_numbers);
                        lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).which_event_coh=which_event;
                        lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).f_coh=f;
                        lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).all_Cxy_timecourse=all_Cxy_timecourse;
                        lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).trial_numbers=trial_numbers;
                        lfp_per_file(filNum).lfpevpair(lfp_per_file(filNum).lfpevpair_no).perCorrCoh=perCorr_pertr;
                    end
                    
                    fprintf(1, 'File number: %d, window number: %d, lfp number: %d\n',filNum,winNo,lfppairNo);
                    
                    
                end
                
            end
            
            %Save this temp file
            drgSavePar([handlespf.drgb.outPathName tempDirName '/'],this_jt(10:end),lfp_per_file(filNum),filNum)
            
            
        end
    end
    
    
    %Now enter the results in the output structures
    %Initialize lfpevpair
    
    handles.drgb.lfpevpair=[];
    handles.drgb.lfpevpair_no=0;
    
    if (sum(handles.drgbchoices.analyses==3)>0)||(sum(handles.drgbchoices.analyses==5)>0)
        handles.drgb.lfpevpair.out_times=lfp_per_file(1).out_times;
    end
    
    for filNum=first_file:handles.drgbchoices.no_files
        
        
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
            handles.drgb.freq_for_LFPpower=lfp_per_file(filNum).f;
        end
        
        handles.drgb.file(filNum).drg=lfp_per_file(filNum).drg;
        
        for ii=1:lfp_per_file(filNum).lfpevpair_no
            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).fileNo=...
                lfp_per_file(filNum).lfpevpair(ii).fileNo;
            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).elec_pair_No=...
                lfp_per_file(filNum).lfpevpair(ii).elec_pair_No;
              handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).elec1=...
                lfp_per_file(filNum).lfpevpair(ii).elec1;
              handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).elec2=...
                lfp_per_file(filNum).lfpevpair(ii).elec2;
            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).referenceEvent=...
                lfp_per_file(filNum).lfpevpair(ii).referenceEvent;
            handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).timeWindow=...
                lfp_per_file(filNum).lfpevpair(ii).timeWindow;
            
            %Coherence
            if (sum(handles.drgbchoices.analyses==2)>0)||(sum(handles.drgbchoices.analyses==8)>0)
                
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).out_times_coh=...
                    lfp_per_file(filNum).lfpevpair(ii).out_times_coh;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).no_trials_coh=...
                    lfp_per_file(filNum).lfpevpair(ii).no_trials_coh;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).which_event_coh=...
                    lfp_per_file(filNum).lfpevpair(ii).which_event_coh;
                
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).f_coh=...
                    lfp_per_file(filNum).lfpevpair(ii).f_coh;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).all_Cxy_timecourse=...
                    lfp_per_file(filNum).lfpevpair(ii).all_Cxy_timecourse;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).trial_numbers=...
                    lfp_per_file(filNum).lfpevpair(ii).trial_numbers;
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).perCorrCoh=...
                    lfp_per_file(filNum).lfpevpair(ii).perCorrCoh;
           
                
            end
            
            %LFP PAC
            if (sum(handles.drgbchoices.analyses==1)>0)
                handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).PAC=...
                    lfp_per_file(filNum).lfpevpair(ii).PAC;
            end
            
            %LFP ERP
            if (sum(handles.drgbchoices.analyses==3)>0)||(sum(handles.drgbchoices.analyses==5)>0)
                
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
    
    
    %Save output file
    handles_drgb=handles;
    if isfield(handles,'data_dg')
        handles_drgb=rmfield(handles_drgb,'data_dg');
    end
    save([handles.drgb.outPathName handles.drgb.outFileName],'handles_drgb','-v7.3')
    
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
    
end







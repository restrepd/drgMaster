function drgRunBatchLFP_end_partial_run

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

%If you want to skip files that have already been processed enter the number of the first file
% first_file=handles.drgb.first_file;

% %NOTE: For the moment because of parallel processing I am defaulting to
% %start with file 1
% first_file=1;
% 
% if first_file==1
%     handles.drgb.lfpevpair_no=0;
%     handles.drgb.lfp_per_exp_no=0;
%     first_out=1;
% else
%     load([handles.drgb.outPathName handles.drgb.outFileName])
%     handles.drgb=handles_drgb.drgb;
%     %The user may add new files
%     handles.drgbchoices.no_files=new_no_files;
%     handles.drgbchoices.PathName=choicePathName;
%     handles.drgbchoices.FileName=choiceFileName;
%     first_out=0;
% end

% test_batch=handles.drgbchoices.test_batch;

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
    
    files_with_temp=[];
    for filNum=first_file:no_files
%     for filNum=first_file:no_files
        
%     for filNum=first_file:handles.drgbchoices.no_files
        
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
            files_with_temp=filNum;
                
        end
    end
    
    
    %Now enter the results in the output structures
    %Initialize lfpevpair
    
    handles.drgb.lfpevpair=[];
    handles.drgb.lfpevpair_no=0;
    
    if (sum(handles.drgbchoices.analyses==3)>0)||(sum(handles.drgbchoices.analyses==5)>0)
        handles.drgb.lfpevpair.out_times=lfp_per_file(1).out_times;
    end
    
    for filNum=files_with_temp
        
        
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
             
             %LFP ERWA
             if (sum(handles.drgbchoices.analyses==5)>0)
                 
                 
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_no_trialsERWA=...
                     lfp_per_file(filNum).lfpevpair(ii).wave_no_trialsERWA;
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_which_eventERWA=...
                     lfp_per_file(filNum).lfpevpair(ii).wave_which_eventERWA;
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_fERWA=...
                     lfp_per_file(filNum).lfpevpair(ii).wave_fERWA;
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_trialNoERWA=...
                     lfp_per_file(filNum).lfpevpair(ii).wave_trialNoERWA;
                 
                 
                 %Lick-referenced data
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_log_P_t_lick_referenced=...
                     lfp_per_file(filNum).lfpevpair(ii).wave_log_P_t_lick_referenced;
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_lick_phase_per_trialERP=...
                     lfp_per_file(filNum).lfpevpair(ii).wave_lick_phase_per_trialERP;
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_no_licks_per_trial=...
                     lfp_per_file(filNum).lfpevpair(ii).wave_no_licks_per_trial;
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_t_per_lick_per_trial=...
                     lfp_per_file(filNum).lfpevpair(ii).wave_t_per_lick_per_trial;
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_no_ref_licks_per_trial=...
                     lfp_per_file(filNum).lfpevpair(ii).wave_no_ref_licks_per_trial;
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).delta_lick_times=...
                     lfp_per_file(filNum).lfpevpair(ii).delta_lick_times;
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).lick_referenced_LFP_per_trial=...
                     lfp_per_file(filNum).lfpevpair(ii).lick_referenced_LFP_per_trial;
                 
                 %Theta-referenced data
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_log_P_t_ltheta_related=...
                     lfp_per_file(filNum).lfpevpair(ii).wave_log_P_t_ltheta_related;
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_no_thetas_per_trial=...
                     lfp_per_file(filNum).lfpevpair(ii).wave_no_thetas_per_trial;
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_t_per_theta_per_trial=...
                     lfp_per_file(filNum).lfpevpair(ii).wave_t_per_theta_per_trial;
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_theta_referenced_LFP_per_trial=...
                     lfp_per_file(filNum).lfpevpair(ii).wave_theta_referenced_LFP_per_trial;
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).wave_no_ref_evs_per_trial_th=...
                     lfp_per_file(filNum).lfpevpair(ii).wave_no_ref_evs_per_trial_th;
                 handles.drgb.lfpevpair(handles.drgb.lfpevpair_no+ii).delta_lick_times_theta_ref=...
                     lfp_per_file(filNum).lfpevpair(ii).delta_lick_times_theta_ref;
                 
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
    
    
    %Save output file
    handles_drgb=handles;
    if isfield(handles,'data_dg')
        handles_drgb=rmfield(handles_drgb,'data_dg');
    end
    save([handles.drgb.outPathName handles.drgb.outFileName],'handles_drgb','-v7.3')
    
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
    
end






